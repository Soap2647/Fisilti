"""WP3-02: Attention-MIL ile çoklu-omik füzyon — tablo modelleriyle AYNI dış katmanlar.

Her dış katmanda:
  1. Her görünüm için 'blood' ön filtresi (tümörde artan + asiner dışlama) yalnızca
     o görünümü ölçülmüş EĞİTİM örneklerinde öğrenilir; ölçekleme de eğitimden.
  2. Eğitim hastalarının %20'si (gruplu) erken durdurma doğrulaması için ayrılır.
  3. 3 tohumlu topluluk (ensemble); test skoru ortalamadır.
Çıktılar: results/wp3/oof_mil.parquet, results/wp3/attention_mil.parquet,
          results/tables/wp3_mil_summary.tsv
Çalıştırma:  .venv\\Scripts\\python.exe scripts/wp3_02_attention_mil.py [--repeats 5]
"""

from __future__ import annotations

import argparse
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupShuffleSplit

from fisilti.config import PROJECT_ROOT, path, set_seed
from fisilti.cv import UnivariateFilter, acinar_program, auc_with_ci, grouped_folds
from fisilti.data.cptac import OMICS
from fisilti.mil import predict_mil, train_mil

K = {"protein": 50, "mirna": 30, "mrna": 50, "methylation": 50}


def load_all(complete_only: bool = False):
    proc = path("processed") / "cptac_pdac"
    lab = pd.read_csv(proc / "labels.tsv", sep="\t", index_col=0)
    if complete_only:
        # K-014: ölçüm deseni etiketle örtüşür (1-2 görünümlü örneklerin tamamı normal);
        # kısayol öğrenmeyi önlemek için yalnızca 4 görünümü tam örnekler
        pres = pd.read_csv(proc / "omic_presence.tsv", sep="\t", index_col=0)
        lab = lab.loc[pres.index[pres.all(axis=1)]]
    views = {}
    for v in OMICS:
        m = pd.read_parquet(proc / f"{v}.parquet")
        if v == "methylation":
            m = m.clip(1e-3, 1 - 1e-3)
            m = np.log2(m / (1 - m))
        X = m.T.reindex(lab.index).dropna(how="all")
        views[v] = X.loc[:, X.notna().mean() >= 0.7]
    samples = lab.index
    y = (lab.group == "PDAC").astype(int).to_numpy()
    tumors = lab.group == "PDAC"
    acinar = {v: acinar_program(X.reindex(lab.index[tumors]).dropna(how="all"),
                                lab.loc[tumors, "Acinar_fraction"]) for v, X in views.items()}
    return views, lab, samples, y, acinar


def build_fold_inputs(views, samples, y, tr, acinar):
    """Görünüm başına eğitimde öğrenilen filtre + z-ölçekleme; ölçülmemiş örnek = 0 ve maskeli."""
    xs, mask, selected = {}, [], {}
    for v, X in views.items():
        present = samples.isin(X.index)
        tr_v = samples[tr][present[tr]]
        f = UnivariateFilter(k=K[v], direction="up", exclude=acinar[v]).fit(X.loc[tr_v], y[samples.get_indexer(tr_v)])
        Z = f.transform(X.reindex(samples[present]))
        mu, sd = Z.loc[tr_v].mean(), Z.loc[tr_v].std().replace(0, 1)
        full = np.zeros((len(samples), len(f.selected_)), dtype=np.float32)
        full[present] = ((Z - mu) / sd).to_numpy()
        xs[v] = full
        mask.append(present.astype(np.float32))
        selected[v] = f.selected_
    return xs, np.stack(mask, axis=1), selected


def run_fold(r, k, tr, te, views, lab, samples, y, acinar):
    warnings.filterwarnings("ignore")
    import torch
    torch.set_num_threads(1)
    xs, mask, selected = build_fold_inputs(views, samples, y, tr, acinar)
    groups = lab.patient.to_numpy()
    gss = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=100 + r * 10 + k)
    itr, ival = next(gss.split(tr, y[tr], groups[tr]))
    scores, atts = [], []
    for seed in (0, 1, 2):
        model = train_mil(xs, mask, y, tr[itr], tr[ival], seed=seed + 7 * r + 31 * k)
        s, a = predict_mil(model, xs, mask, te)
        scores.append(s)
        atts.append(a)
    score, att = np.mean(scores, 0), np.mean(atts, 0)
    oof = pd.DataFrame({"view": "fusion", "config": "blood", "model": "attention_mil", "repeat": r,
                        "fold": k, "sample": samples[te], "y": y[te], "score": score})
    attd = pd.DataFrame(att, columns=list(views), index=samples[te]).assign(repeat=r, fold=k)
    return oof, attd


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--jobs", type=int, default=-1)
    ap.add_argument("--complete-only", action="store_true",
                    help="yalnızca 4 omiği tam örnekler (K-014)")
    a = ap.parse_args()
    set_seed()
    views, lab, samples, y, acinar = load_all(a.complete_only)
    tag = "_complete" if a.complete_only else ""
    groups = lab.patient.to_numpy()
    folds = grouped_folds(y, groups, n_splits=5, repeats=a.repeats)
    print(f"[..] {len(samples)} örnek ({int(y.sum())} T / {int((1 - y).sum())} N), {len(folds)} dış katman")
    res = Parallel(n_jobs=a.jobs)(delayed(run_fold)(r, k, tr, te, views, lab, samples, y, acinar)
                                  for r, k, tr, te in folds)
    oof = pd.concat([x[0] for x in res], ignore_index=True)
    att = pd.concat([x[1] for x in res])
    out = PROJECT_ROOT / "results" / "wp3"
    oof = oof.assign(model="attention_mil" + tag)
    oof.to_parquet(out / f"oof_mil{tag}.parquet")
    att.to_parquet(out / f"attention_mil{tag}.parquet")

    per_rep = [roc_auc_score(d.y, d.score) for _, d in oof.groupby("repeat")]
    avg = oof.groupby("sample").agg(y=("y", "first"), score=("score", "mean"))
    auc, lo, hi = auc_with_ci(avg.y.to_numpy(), avg.score.to_numpy())
    early = lab.loc[avg.index, "early"]
    sub = avg[(avg.y == 0) | (early == 1)]
    e = auc_with_ci(sub.y.to_numpy(), sub.score.to_numpy())
    # Yalnızca tek bir görünümü ölçülmüş örneklerde performans (eksik veri dayanıklılığı)
    pres = pd.DataFrame({v: samples.isin(X.index) for v, X in views.items()}, index=samples)
    nview = pres.loc[avg.index].sum(axis=1)
    rows = [{"view": "fusion", "config": "blood", "model": "attention_mil" + tag,
             "n_tumor": int(avg.y.sum()), "n_normal": int((1 - avg.y).sum()),
             "auc_mean_over_repeats": np.mean(per_rep), "auc_sd_over_repeats": np.std(per_rep),
             "auc": auc, "auc_ci_lo": lo, "auc_ci_hi": hi,
             "auc_early_vs_normal": e[0], "early_ci_lo": e[1], "early_ci_hi": e[2]}]
    s = pd.DataFrame(rows)
    s.round(4).to_csv(path("tables") / f"wp3_mil_summary{tag}.tsv", sep="\t", index=False)
    print(s.round(3).to_string(index=False))
    for n in sorted(nview.unique()):
        d = avg[nview == n]
        if d.y.nunique() == 2:
            print(f"  {n} görünümü ölçülmüş örnekler (n={len(d)}): AUC {roc_auc_score(d.y, d.score):.3f}")
    ma = att.groupby(level=0).mean(numeric_only=True)[list(views)]
    ma["group"] = lab.loc[ma.index, "group"]
    full = pres.loc[ma.index].all(axis=1)
    print("[ok] 4 görünümü tam örneklerde ortalama dikkat ağırlıkları:")
    print(ma[full].groupby("group").mean().round(3).to_string())


if __name__ == "__main__":
    main()
