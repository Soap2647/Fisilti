"""WP4-01: Dış katman (out-of-sample) SHAP analizi.

Her omik için 'blood' senaryosunda XGBoost, WP3 ile aynı dış katmanlarda (1. tekrar,
5 katman) iç CV + Optuna ile yeniden eğitilir; SHAP değerleri yalnızca modelin
GÖRMEDİĞİ test örneklerinde hesaplanır (XGBoost yerleşik kesin TreeSHAP,
`pred_contribs=True`). Böylece açıklamalar da ezberden arındırılmış olur.

Bir özellik bir katmanda modele girmemişse o katmandaki test örnekleri için
katkısı tanım gereği 0'dır.

Çıktılar:
  results/tables/wp4_shap_<view>.tsv   özellik başına ortalama |SHAP|, yön, kararlılık
  results/figures/wp4_shap_<view>.png  en etkili 20 özellik (beeswarm)
Çalıştırma:  .venv\\Scripts\\python.exe scripts/wp4_01_shap.py
"""

from __future__ import annotations

import sys
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap
import xgboost as xgb
from scipy.stats import spearmanr

from fisilti.config import PROJECT_ROOT, path, set_seed
from fisilti.cv import acinar_program, grouped_folds
from fisilti.models import tune_and_fit

sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from wp3_01_tabular_cv import load_view  # noqa: E402

VIEW_LABEL = {"protein": "Protein", "mirna": "miRNA", "mrna": "mRNA", "methylation": "Metilasyon"}


def main() -> None:
    warnings.filterwarnings("ignore")
    set_seed()
    names = pd.read_csv(path("processed") / "geo" / "mirna_names.tsv", sep="\t", index_col=0)["name"]
    for view in ["protein", "mirna", "mrna", "methylation"]:
        X, y, lab = load_view(view)
        groups = lab.patient.to_numpy()
        tumors = lab.group == "PDAC"
        acinar = acinar_program(X[tumors.values], lab.loc[tumors, "Acinar_fraction"])
        filt = {"k": 200, "direction": "up", "exclude": acinar}
        folds = [f for f in grouped_folds(y, groups, n_splits=5, repeats=1)]
        shap_parts, n_in_model = [], pd.Series(0, index=X.columns, dtype=int)
        for r, k, tr, te in folds:
            model, _, _ = tune_and_fit("xgboost", X.iloc[tr], y[tr], groups[tr], filt, n_trials=15)
            feats = model.named_steps["filter"].selected_
            Xte = model.named_steps["filter"].transform(X.iloc[te])
            contrib = model.named_steps["clf"].get_booster().predict(
                xgb.DMatrix(Xte, feature_names=list(feats)), pred_contribs=True)[:, :-1]  # son sütun: taban değer
            shap_parts.append(pd.DataFrame(contrib, index=X.index[te], columns=feats))
            used = np.abs(contrib).sum(axis=0) > 0
            n_in_model[np.array(feats)[used]] += 1
        S = pd.concat(shap_parts).fillna(0.0)                        # örnek x özellik (dış katman)
        V = X.loc[S.index, S.columns]
        V = V.fillna(V.median())

        tab = pd.DataFrame({
            "mean_abs_shap": S.abs().mean(),
            "folds_in_model": n_in_model.reindex(S.columns).fillna(0).astype(int),
        })
        # Yön: özellik değeri arttıkça PDAC olasılığına katkı artıyor mu? (Spearman değer ~ SHAP)
        with np.errstate(invalid="ignore"):
            tab["direction_rho"] = [spearmanr(V[c], S[c])[0] if S[c].std() > 0 else np.nan for c in S.columns]
        if view == "mirna":
            tab.insert(0, "name", [names.get(i, i) for i in tab.index])
        else:
            tab.insert(0, "name", tab.index)
        tab = tab.sort_values("mean_abs_shap", ascending=False)
        tab.round(5).to_csv(path("tables") / f"wp4_shap_{view}.tsv", sep="\t", index_label="feature")

        top = tab.index[:20]
        plt.figure()
        shap.summary_plot(S[top].to_numpy(), V[top].to_numpy(), feature_names=list(tab.loc[top, "name"]),
                          max_display=20, show=False, plot_size=(6.5, 6.5))
        plt.title(f"{VIEW_LABEL[view]} · dış katman SHAP (n={len(S)} örnek, XGBoost)", fontsize=10)
        plt.tight_layout()
        plt.savefig(path("figures") / f"wp4_shap_{view}.png", dpi=160)
        plt.close("all")

        pos = (tab.direction_rho.head(20) > 0).mean()
        print(f"[ok] {view}: {len(S)} dış katman örneği, {S.shape[1]} özellik; "
              f"en etkili 20'nin %{pos * 100:.0f}'i 'yüksek değer -> PDAC' yönünde")
        print(tab.head(10)[["name", "mean_abs_shap", "folds_in_model", "direction_rho"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
