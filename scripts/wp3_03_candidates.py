"""WP3-03: Kararlılık seçimiyle biyobelirteç aday listesi (yalnızca 'blood' senaryosu).

Her görünüm için, tüm dış katmanlarda (tekrar x katman) ve üç modelde eğitilen
modeller üzerinden:
  sel_freq      : özelliğin modelde sıfırdan büyük önemle yer aldığı model oranı
  mean_imp      : seçildiği modellerde ortalama normalize önem
  auc_all       : tüm veride tümör vs normal tek değişkenli AUC (betimleyici)
  auc_early     : erken evre tümör vs normal AUC
  log2fc        : medyan(tümör) - medyan(normal)  (log ölçekte)
  rho_acinar    : tümörlerde asiner fraksiyonla Spearman ρ
Sıralama: sel_freq (azalan), eşitlikte mean_imp.

Not: auc_all/auc_early/log2fc TÜM veride hesaplanır ve yalnızca betimleyicidir;
performans tahmini değildir (performans: dış katman tahminleri, wp3_01/02).

Çıktı: results/tables/wp3_candidates_<view>.tsv
Çalıştırma:  .venv\\Scripts\\python.exe scripts/wp3_03_candidates.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from fisilti.config import PROJECT_ROOT, path
from fisilti.cv import fast_auc


def main() -> None:
    w3 = PROJECT_ROOT / "results" / "wp3"
    lab = pd.read_csv(path("processed") / "cptac_pdac" / "labels.tsv", sep="\t", index_col=0)
    names = pd.read_csv(path("processed") / "geo" / "mirna_names.tsv", sep="\t", index_col=0)["name"]
    for f in sorted(w3.glob("features_*.parquet")):
        view = f.stem.split("_", 1)[1]
        feats = pd.read_parquet(f)
        feats = feats[feats.config == "blood"]
        n_models = feats.groupby(["model", "repeat", "fold"]).ngroups
        used = feats[feats.importance > 0]
        agg = used.groupby("feature").agg(n_sel=("importance", "size"), mean_imp=("importance", "mean"))
        agg["sel_freq"] = agg.n_sel / n_models
        per_model = used.groupby(["feature", "model"]).size().unstack(fill_value=0)
        per_model = per_model.div(feats.groupby("model")[["repeat", "fold"]].apply(
            lambda d: len(d.drop_duplicates())), axis=1)
        agg = agg.join(per_model.add_prefix("freq_"))

        # Betimleyici istatistikler (tüm veri)
        m = pd.read_parquet(path("processed") / "cptac_pdac" / f"{view}.parquet")
        if view == "methylation":
            m = np.log2(m.clip(1e-3, 1 - 1e-3) / (1 - m.clip(1e-3, 1 - 1e-3)))
        X = m.T.loc[:, agg.index]
        l = lab.loc[X.index]
        y = (l.group == "PDAC").astype(int).to_numpy()
        Xf = X.fillna(X.median())
        agg["auc_all"] = fast_auc(Xf.to_numpy(), y)
        e = ((l.group == "normal") | (l.early == 1)).to_numpy()
        agg["auc_early"] = fast_auc(Xf.to_numpy()[e], y[e])
        agg["log2fc"] = X[y == 1].median() - X[y == 0].median()
        t = (l.group == "PDAC") & l.Acinar_fraction.notna()
        agg["rho_acinar"] = [spearmanr(X.loc[t, c], l.loc[t, "Acinar_fraction"], nan_policy="omit")[0]
                             for c in agg.index]
        agg.insert(0, "name", [names.get(i, i) for i in agg.index] if view == "mirna" else agg.index)
        agg = agg.sort_values(["sel_freq", "mean_imp"], ascending=False)
        agg.round(4).to_csv(path("tables") / f"wp3_candidates_{view}.tsv", sep="\t", index_label="feature")
        top = agg.head(15)
        print(f"== {view}: {len(agg)} özellik en az bir kez seçildi; {n_models} model; "
              f"sel_freq ≥ 0.8: {(agg.sel_freq >= 0.8).sum()}")
        print(top[["name", "sel_freq", "mean_imp", "auc_all", "auc_early", "log2fc", "rho_acinar"]]
              .round(3).to_string(index=False))


if __name__ == "__main__":
    main()
