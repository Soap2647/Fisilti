"""WP6-02 (KEŞİF AMAÇLI; ön-kayıtlı değildir, K-019 sonrası): başarısızlığın nedenleri.

E1  Doku-serum uyumu: her miRNA için doku log2FC (CPTAC tümör-normal) ile serum
    etki büyüklüğü (PDAC-sağlıklı, Cohen d) arasındaki Spearman korelasyonu.
E2  Serumda tekrarlanabilir sinyal var mı? Kohortlar arası çapraz doğrulama:
    bir serum kohortunda (her iki platformda ortak + ölçülebilir miRNA'lar)
    L1-lojistik model eğitilir (iç CV ile C seçimi), DİĞER kohortta test edilir;
    sonra tersi. Özellikler kohort-içi z-skorlanır.
E3  Aynı kohortlar arası tasarım, özellikler doku adaylarıyla (WP3, sel_freq ≥ 0.5)
    sınırlanarak: doku bilgisi serum modeline katkı sağlıyor mu?
Çıktı: results/tables/wp6_exploratory.tsv, results/figures/wp6_tissue_vs_serum.png
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import roc_auc_score

from fisilti.config import path, set_seed
from fisilti.cv import auc_with_ci

FLOOR = 2.612


def load(gse):
    e = pd.read_parquet(path("processed") / "geo" / f"{gse}_expr.parquet")
    l = pd.read_csv(path("processed") / "geo" / f"{gse}_labels.tsv", sep="\t", index_col=0).loc[e.columns]
    return e, l


def cohen_d(a, b):
    return (np.nanmean(a) - np.nanmean(b)) / np.sqrt((np.nanvar(a) + np.nanvar(b)) / 2)


def main() -> None:
    set_seed()
    proc = path("processed") / "cptac_pdac"
    lab = pd.read_csv(proc / "labels.tsv", sep="\t", index_col=0)
    tis = pd.read_parquet(proc / "mirna.parquet")
    g = lab.loc[tis.columns, "group"]
    tissue_fc = tis.loc[:, g == "PDAC"].median(axis=1) - tis.loc[:, g == "normal"].median(axis=1)

    e59, l59 = load("GSE59856")
    e85, l85 = load("GSE85589")
    det = (e59 > FLOOR + 1e-6).mean(axis=1) >= 0.5
    common = sorted(set(det[det].index) & set(e85.index))
    rows = []

    # E1
    fig, axes = plt.subplots(1, 2, figsize=(9, 4), dpi=160)
    for ax, (gse, e, l) in zip(axes, [("GSE85589", e85, l85), ("GSE59856", e59, l59)]):
        pos, neg = l.group == "PDAC", l.group == "healthy"
        feats = [m for m in tissue_fc.index if m in e.index and m in common]
        sd = pd.Series({m: cohen_d(e.loc[m, pos], e.loc[m, neg]) for m in feats})
        rho, p = spearmanr(tissue_fc[feats], sd)
        rows.append({"analysis": "E1_doku_serum_uyumu", "cohort": gse, "n_mirna": len(feats),
                     "value": rho, "p": p, "note": "Spearman(doku log2FC, serum Cohen d)"})
        ax.scatter(tissue_fc[feats], sd, s=8, alpha=0.5, c="#555")
        ax.axhline(0, c="#bbb", lw=0.8)
        ax.axvline(0, c="#bbb", lw=0.8)
        ax.set_xlabel("Doku: log2 kat değişimi (tümör − normal)")
        ax.set_ylabel("Serum: Cohen d (PDAC − sağlıklı)")
        ax.set_title(f"{gse} · ρ = {rho:.2f} (p = {p:.2g}, n = {len(feats)})", fontsize=9)
    fig.suptitle("Dokuda artan miRNA'lar serumda artıyor mu?", fontsize=10)
    fig.tight_layout()
    fig.savefig(path("figures") / "wp6_tissue_vs_serum.png")
    plt.close(fig)

    # E2 / E3
    cand = pd.read_csv(path("tables") / "wp3_candidates_mirna.tsv", sep="\t", index_col=0)
    tissue_set = [m for m in cand.index[cand.sel_freq >= 0.5] if m in common]
    cohorts = {"GSE59856": (e59, l59), "GSE85589": (e85, l85)}

    def xy(gse, feats):
        e, l = cohorts[gse]
        m = l.group.isin(["PDAC", "healthy"])
        z = e.loc[feats]
        z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1), axis=0).fillna(0)
        return z.loc[:, m].T.to_numpy(), (l.group[m] == "PDAC").astype(int).to_numpy()

    for tag, feats in [("E2_serum_tum_ortak", common), ("E3_doku_adaylari", tissue_set)]:
        for tr, te in [("GSE59856", "GSE85589"), ("GSE85589", "GSE59856")]:
            Xtr, ytr = xy(tr, feats)
            Xte, yte = xy(te, feats)
            clf = LogisticRegressionCV(Cs=np.logspace(-3, 1, 12), penalty="l1", solver="liblinear",
                                       cv=5, scoring="roc_auc", class_weight="balanced", max_iter=5000,
                                       random_state=42).fit(Xtr, ytr)
            s = clf.decision_function(Xte)
            auc, lo, hi = auc_with_ci(yte, s)
            nz = int((clf.coef_ != 0).sum())
            rows.append({"analysis": tag, "cohort": f"{tr} -> {te}", "n_mirna": len(feats),
                         "value": auc, "p": np.nan, "note": f"dış AUC {lo:.3f}-{hi:.3f}; sıfırdan farklı katsayı {nz}; "
                                                               f"eğitim içi CV AUC {clf.scores_[1].mean(axis=0).max():.3f}"})
    t = pd.DataFrame(rows)
    t.round(4).to_csv(path("tables") / "wp6_exploratory.tsv", sep="\t", index=False)
    print(t.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
