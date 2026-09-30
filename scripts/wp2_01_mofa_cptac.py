"""WP2-01: CPTAC-PDAC üzerinde MOFA+ çoklu-omik faktör analizi (keşif, denetimsiz).

Adımlar
  1. 4 görünüm (mRNA, protein, miRNA, metilasyon) varyansa göre ön seçilir.
  2. MOFA+ (15 faktör üst sınırı, ARD + spike-slab) tüm CPTAC örneklerinde eğitilir.
  3. Faktörler model kurulduktan SONRA klinik değişkenlerle ilişkilendirilir:
       - tümör vs normal (AUC, Mann-Whitney)
       - erken (I-II) vs geç (III-IV) evre, yalnızca tümörler
       - doku bileşimi (neoplastik / asiner / stromal fraksiyon) — K-007 karıştırıcı denetimi
  4. Her faktörün görünüm başına en yüksek ağırlıklı özellikleri kaydedilir.

Çıktılar: results/tables/wp2_*.tsv, results/figures/wp2_*.png, results/models/mofa_cptac.hdf5
Çalıştırma:  .venv\\Scripts\\python.exe scripts/wp2_01_mofa_cptac.py
"""

from __future__ import annotations

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.metrics import roc_auc_score
from statsmodels.stats.multitest import multipletests

from fisilti.config import PROJECT_ROOT, path, set_seed
from fisilti.data.cptac import OMICS
from fisilti.mofa import fit_mofa, load_mofa, prepare_views

VIEW_LABEL = {"mrna": "mRNA", "protein": "Protein", "mirna": "miRNA", "methylation": "Metilasyon"}


def associations(factors: pd.DataFrame, lab: pd.DataFrame) -> pd.DataFrame:
    rows = []
    tumor = lab.group == "PDAC"
    staged = tumor & lab.early.notna()
    for f in factors.columns:
        z = factors[f]
        r = {"factor": f}
        # Tümör vs normal
        y = tumor.astype(int)
        r["auc_tumor_vs_normal"] = roc_auc_score(y, z)
        r["p_tumor_vs_normal"] = mannwhitneyu(z[tumor], z[~tumor]).pvalue
        # Erken vs geç (yalnızca tümör)
        ye = lab.loc[staged, "early"].astype(int)
        r["auc_early_vs_late"] = roc_auc_score(ye, z[staged])
        r["p_early_vs_late"] = mannwhitneyu(z[staged][ye == 1], z[staged][ye == 0]).pvalue
        # Doku bileşimi (yalnızca tümör)
        for c, short in [("Neoplastic_cellularity", "neoplastic"), ("Acinar_fraction", "acinar"),
                         ("Stromal_fraction", "stromal")]:
            ok = tumor & lab[c].notna()
            rho, p = spearmanr(z[ok], lab.loc[ok, c])
            r[f"rho_{short}"], r[f"p_{short}"] = rho, p
        rows.append(r)
    t = pd.DataFrame(rows).set_index("factor")
    pcols = [c for c in t.columns if c.startswith("p_")]
    # Tüm testler üzerinden Benjamini-Hochberg FDR
    q = multipletests(t[pcols].to_numpy().ravel(), method="fdr_bh")[1].reshape(t[pcols].shape)
    for i, c in enumerate(pcols):
        t["q_" + c[2:]] = q[:, i]
    return t


def plot_r2(r2: pd.DataFrame, file) -> None:
    d = r2.drop(index="toplam").rename(columns=VIEW_LABEL)
    fig, ax = plt.subplots(figsize=(4.2, 0.32 * len(d) + 1.4), dpi=160)
    im = ax.imshow(d.to_numpy(), cmap="Purples", aspect="auto", vmin=0)
    ax.set_xticks(range(d.shape[1]), d.columns, rotation=30, ha="right", fontsize=8)
    ax.set_yticks(range(d.shape[0]), d.index, fontsize=8)
    for i in range(d.shape[0]):
        for j in range(d.shape[1]):
            v = d.iat[i, j]
            ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=6.5,
                    color="white" if v > d.to_numpy().max() * 0.6 else "black")
    tot = r2.loc["toplam"].rename(VIEW_LABEL)
    ax.set_title("MOFA+ açıklanan varyans (%)\nToplam: " +
                 ", ".join(f"{k} {v:.0f}" for k, v in tot.items()), fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.05)
    fig.tight_layout()
    fig.savefig(file)
    plt.close(fig)


def plot_scatter(factors, lab, fx, fy, file) -> None:
    fig, ax = plt.subplots(figsize=(5, 4.2), dpi=160)
    groups = [("normal", "#2e86c1", "Normal"), ("late", "#7b241c", "PDAC geç (III-IV)"),
              ("early", "#e74c3c", "PDAC erken (I-II)"), ("nostage", "#bbbbbb", "PDAC evre ?")]
    cat = np.select([lab.group == "normal", lab.early == 1, lab.early == 0],
                    ["normal", "early", "late"], "nostage")
    for key, col, name in groups:
        m = cat == key
        ax.scatter(factors.loc[m, fx], factors.loc[m, fy], s=14, c=col, alpha=0.8,
                   edgecolors="none", label=f"{name} (n={m.sum()})")
    ax.set_xlabel(fx)
    ax.set_ylabel(fy)
    ax.legend(fontsize=7, frameon=False)
    ax.set_title("CPTAC-PDAC · MOFA+ faktör uzayı", fontsize=9)
    fig.tight_layout()
    fig.savefig(file)
    plt.close(fig)


def main() -> None:
    set_seed()
    proc = path("processed") / "cptac_pdac"
    lab = pd.read_csv(proc / "labels.tsv", sep="\t", index_col=0)
    raw = {o: pd.read_parquet(proc / f"{o}.parquet") for o in OMICS}
    views = prepare_views(raw, samples=list(lab.index))
    for v, x in views.items():
        print(f"[..] {v:<11} {x.shape[1]:>5} özellik | gözlenen örnek {int(x.notna().any(axis=1).sum())}")

    model = PROJECT_ROOT / "results" / "models" / "mofa_cptac.hdf5"
    if model.exists() and "--refit" not in sys.argv:
        print(f"[..] kayıtlı model yükleniyor: {model.name} (yeniden eğitmek için --refit)")
        res, r2_mofa = load_mofa(model, views)
        # Bağımsız doğrulama: kendi R² hesabımız MOFA+'ın iç hesabıyla örtüşmeli
        diff = (res.r2.drop(index="toplam") - r2_mofa).abs().to_numpy().max()
        print(f"[ok] R² çapraz kontrol: MOFA+ iç hesabıyla en büyük fark = {diff:.3f} puan")
    else:
        res = fit_mofa(views, n_factors=15, seed=42, outfile=model, verbose=True)
    print(f"[ok] MOFA+: {res.factors.shape[1]} aktif faktör")

    tables, figs = path("tables"), path("figures")
    res.factors.to_csv(tables / "wp2_mofa_factors.tsv", sep="\t", index_label="sample")
    res.r2.round(2).to_csv(tables / "wp2_mofa_r2.tsv", sep="\t", index_label="factor")
    print(res.r2.round(1).to_string())

    assoc = associations(res.factors, lab.loc[res.factors.index])
    assoc.round(4).to_csv(tables / "wp2_mofa_factor_associations.tsv", sep="\t")
    show = ["auc_tumor_vs_normal", "q_tumor_vs_normal", "auc_early_vs_late", "q_early_vs_late",
            "rho_neoplastic", "rho_acinar", "rho_stromal"]
    print(assoc[show].round(3).to_string())

    # Görünüm başına en yüksek |ağırlık| taşıyan 25 özellik
    rows = []
    for v, w in res.weights.items():
        for f in w.columns:
            top = w[f].abs().sort_values(ascending=False).index[:25]
            for rank, feat in enumerate(top, 1):
                rows.append({"factor": f, "view": v, "rank": rank,
                             "feature": feat.split(":", 1)[1], "weight": w.at[feat, f]})
    pd.DataFrame(rows).to_csv(tables / "wp2_mofa_top_weights.tsv", sep="\t", index=False)

    plot_r2(res.r2, figs / "wp2_mofa_variance_explained.png")
    # Tümör/normali en iyi ayıran iki faktör
    best = (assoc.auc_tumor_vs_normal - 0.5).abs().sort_values(ascending=False).index[:2]
    plot_scatter(res.factors, lab.loc[res.factors.index], best[0], best[1],
                 figs / "wp2_mofa_factor_scatter.png")
    print(f"[ok] çıktılar yazıldı; tümör/normal için en ayırt edici faktörler: {list(best)}")


if __name__ == "__main__":
    main()
