"""WP3-01: CPTAC-PDAC'ta omik başına tümör/normal sınıflandırması — iç içe, hasta gruplu CV.

Senaryolar (config):
  free   -> tüm özellikler, iki yön (dokudaki üst sınır performans)
  blood  -> yalnızca tümörde artan özellikler + asiner program dışlaması (K-012);
            panel adayları buradan türetilir
Modeller: L1-lojistik (referans), XGBoost, LightGBM; iç 3-katlı grup-CV'de Optuna.

Çıktılar (results/wp3/):
  oof_<view>.parquet        dış katman tahminleri (tekrar x örnek)
  features_<view>.parquet   her dış katmanda seçilen özellikler ve önemleri
  results/tables/wp3_tabular_summary.tsv

Çalıştırma:  .venv\\Scripts\\python.exe scripts/wp3_01_tabular_cv.py [--repeats 5] [--trials 15] [--views protein mirna ...]
"""

from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from fisilti.config import PROJECT_ROOT, path, set_seed
from fisilti.cv import acinar_program, auc_with_ci, grouped_folds
from fisilti.models import MODELS, feature_importance, tune_and_fit

CONFIGS = {
    "free": {"direction": "both"},
    "blood": {"direction": "up"},
}


def load_view(view: str):
    proc = path("processed") / "cptac_pdac"
    lab = pd.read_csv(proc / "labels.tsv", sep="\t", index_col=0)
    m = pd.read_parquet(proc / f"{view}.parquet")
    if view == "methylation":
        m = m.clip(1e-3, 1 - 1e-3)
        m = np.log2(m / (1 - m))  # M değeri
    X = m.T
    X = X.loc[:, X.notna().mean() >= 0.7]
    lab = lab.loc[X.index]
    y = (lab.group == "PDAC").astype(int).to_numpy()
    return X, y, lab


def run_fold(view, cfg, model, r, k, tr, te, X, y, groups, acinar, k_feats, trials):
    filt = {"k": k_feats, "direction": CONFIGS[cfg]["direction"],
            "exclude": acinar if cfg == "blood" else frozenset()}
    fitted, params, inner_auc = tune_and_fit(model, X.iloc[tr], y[tr], groups[tr], filt,
                                             n_trials=trials, seed=42 + r)
    score = fitted.predict_proba(X.iloc[te])[:, 1]
    imp = feature_importance(fitted)
    oof = pd.DataFrame({"view": view, "config": cfg, "model": model, "repeat": r, "fold": k,
                        "sample": X.index[te], "y": y[te], "score": score})
    feats = pd.DataFrame({"view": view, "config": cfg, "model": model, "repeat": r, "fold": k,
                          "feature": imp.index, "importance": imp.values,
                          "train_auc_univariate": fitted.named_steps["filter"].auc_.reindex(imp.index).values})
    return oof, feats, {"view": view, "config": cfg, "model": model, "repeat": r, "fold": k,
                        "inner_auc": inner_auc, **{f"p_{a}": b for a, b in params.items()}}


def summarize(oof: pd.DataFrame, lab: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (view, cfg, model), d in oof.groupby(["view", "config", "model"]):
        per_rep = []
        for _, dr in d.groupby("repeat"):
            from sklearn.metrics import roc_auc_score
            per_rep.append(roc_auc_score(dr.y, dr.score))
        # Tekrarlar üzerinden ortalama skor -> bootstrap CI
        avg = d.groupby("sample").agg(y=("y", "first"), score=("score", "mean"))
        auc, lo, hi = auc_with_ci(avg.y.to_numpy(), avg.score.to_numpy())
        early = lab.loc[avg.index, "early"]
        sub = avg[(avg.y == 0) | (early == 1)]
        e_auc, e_lo, e_hi = auc_with_ci(sub.y.to_numpy(), sub.score.to_numpy())
        late = avg[(avg.y == 0) | (early == 0)]
        l_auc = auc_with_ci(late.y.to_numpy(), late.score.to_numpy(), n_boot=200)[0]
        rows.append({"view": view, "config": cfg, "model": model,
                     "n_tumor": int((avg.y == 1).sum()), "n_normal": int((avg.y == 0).sum()),
                     "auc_mean_over_repeats": np.mean(per_rep), "auc_sd_over_repeats": np.std(per_rep),
                     "auc": auc, "auc_ci_lo": lo, "auc_ci_hi": hi,
                     "auc_early_vs_normal": e_auc, "early_ci_lo": e_lo, "early_ci_hi": e_hi,
                     "auc_late_vs_normal": l_auc})
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--views", nargs="+", default=["protein", "mirna", "mrna", "methylation"])
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--trials", type=int, default=15)
    ap.add_argument("--k", type=int, default=200, help="ön filtrede tutulacak özellik sayısı")
    ap.add_argument("--jobs", type=int, default=-1)
    a = ap.parse_args()
    set_seed()
    out = PROJECT_ROOT / "results" / "wp3"
    out.mkdir(parents=True, exist_ok=True)
    summaries, logs = [], []
    for view in a.views:
        t0 = time.time()
        X, y, lab = load_view(view)
        groups = lab.patient.to_numpy()
        tumors = lab.group == "PDAC"
        acinar = acinar_program(X[tumors.values], lab.loc[tumors, "Acinar_fraction"])
        folds = grouped_folds(y, groups, n_splits=a.folds, repeats=a.repeats)
        print(f"[..] {view}: {X.shape[1]} özellik, {int(y.sum())} tümör / {int((1 - y).sum())} normal, "
              f"{groups.size} örnek / {len(set(groups))} hasta | asiner program: {len(acinar)} özellik | "
              f"{len(folds)} dış katman x {len(CONFIGS)} senaryo x {len(MODELS)} model")
        jobs = [delayed(run_fold)(view, cfg, model, r, k, tr, te, X, y, groups, acinar, a.k, a.trials)
                for cfg in CONFIGS for model in MODELS for r, k, tr, te in folds]
        res = Parallel(n_jobs=a.jobs, verbose=0)(jobs)
        oof = pd.concat([r[0] for r in res], ignore_index=True)
        feats = pd.concat([r[1] for r in res], ignore_index=True)
        logs += [r[2] for r in res]
        oof.to_parquet(out / f"oof_{view}.parquet")
        feats.to_parquet(out / f"features_{view}.parquet")
        pd.Series(sorted(acinar), name="feature").to_csv(out / f"acinar_program_{view}.tsv", sep="\t", index=False)
        s = summarize(oof, lab)
        summaries.append(s)
        print(s[["config", "model", "auc", "auc_ci_lo", "auc_ci_hi", "auc_sd_over_repeats",
                 "auc_early_vs_normal"]].round(3).to_string(index=False))
        print(f"[ok] {view}: {time.time() - t0:.0f} sn")
    pd.concat(summaries).round(4).to_csv(path("tables") / "wp3_tabular_summary.tsv", sep="\t", index=False)
    pd.DataFrame(logs).to_csv(out / "tuning_log.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
