"""WP1-08: Kalite kontrol — kohort özet tablosu ve PCA grafikleri.

Amaç: (1) etiketlerin doğruluğunu görsel olarak sınamak (tümör/normal ayrışmalı),
(2) aykırı örnekleri ve olası batch etkilerini erkenden yakalamak.

Çıktılar:
  results/tables/wp1_cohort_summary.tsv
  results/figures/wp1_pca_<kohort>_<omik>.png
Çalıştırma:  .venv\\Scripts\\python.exe scripts/wp1_08_qc.py
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from fisilti.config import path

PALETTE = {"PDAC": "#c0392b", "normal": "#2e86c1", "healthy": "#27ae60",
           "benign": "#f39c12", "other_cancer": "#95a5a6"}


def cohorts():
    proc = path("processed")
    for gse in ["GSE15471", "GSE62452", "GSE16515", "GSE85589", "GSE59856"]:
        e = pd.read_parquet(proc / "geo" / f"{gse}_expr.parquet")
        l = pd.read_csv(proc / "geo" / f"{gse}_labels.tsv", sep="\t", index_col=0)
        yield gse, "mirna" if gse in ("GSE85589", "GSE59856") else "mrna", e, l
    for name in ["cptac_pdac", "tcga_paad"]:
        l = pd.read_csv(proc / name / "labels.tsv", sep="\t", index_col=0)
        for f in sorted((proc / name).glob("*.parquet")):
            e = pd.read_parquet(f)
            yield name, f.stem, e, l.loc[e.columns]


def pca_plot(name: str, omic: str, e: pd.DataFrame, l: pd.DataFrame) -> dict:
    # numpy ile: eksik oranı <%20 özellikler -> en değişken 2000 -> eksikler özellik
    # ortalamasıyla doldurulur (yalnızca görselleştirme; 486k metilasyon probunda da hızlı)
    v = e.to_numpy(dtype=np.float32)
    ok = np.isnan(v).mean(axis=1) < 0.2
    v = v[ok]
    var = np.nanvar(v, axis=1)
    idx = np.argsort(var)[::-1][:2000]
    idx = idx[var[idx] > 0]
    x = v[idx]
    mu = np.nanmean(x, axis=1, keepdims=True)
    x = np.where(np.isnan(x), mu, x)
    z = ((x - mu) / x.std(axis=1, keepdims=True)).T
    top = idx
    pca = PCA(n_components=2, random_state=0)
    pc = pca.fit_transform(z)

    fig, ax = plt.subplots(figsize=(5.5, 4.5), dpi=150)
    for g in [g for g in PALETTE if g in set(l.group)]:
        m = (l.group == g).values
        ax.scatter(pc[m, 0], pc[m, 1], s=14, alpha=0.75, c=PALETTE[g], label=f"{g} (n={m.sum()})",
                   edgecolors="none")
    ev = pca.explained_variance_ratio_ * 100
    ax.set_xlabel(f"PC1 (%{ev[0]:.1f})")
    ax.set_ylabel(f"PC2 (%{ev[1]:.1f})")
    ax.set_title(f"{name} · {omic} · {len(top)} en değişken özellik", fontsize=9)
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(path("figures") / f"wp1_pca_{name}_{omic}.png")
    plt.close(fig)

    # Aykırı örnek: PC1-PC2 uzayında merkeze uzaklığı > 4 SD
    d = np.sqrt(((pc - pc.mean(0)) / pc.std(0)) ** 2).sum(1)
    outliers = list(e.columns[d > 4 * np.sqrt(2)])
    # Tümör/normal ayrışması: PC1 üzerinde gruplar arası etki büyüklüğü (Cohen's d)
    sep = np.nan
    if {"PDAC"} < set(l.group) and len(set(l.group)) > 1:
        a, b = pc[(l.group == "PDAC").values, 0], pc[(l.group != "PDAC").values, 0]
        sep = abs(a.mean() - b.mean()) / np.sqrt((a.var() + b.var()) / 2)
    return {"pc1_var_%": round(float(ev[0]), 1), "pc1_cohen_d_pdac_vs_rest": round(float(sep), 2),
            "outliers": ",".join(outliers)}


def main() -> None:
    rows = []
    for name, omic, e, l in cohorts():
        g = l.group.value_counts()
        pdac = l[l.group == "PDAC"].drop_duplicates("patient")
        row = {
            "cohort": name, "omic": omic, "features": e.shape[0], "samples": e.shape[1],
            "patients": l.patient.nunique(),
            **{f"n_{k}": int(g.get(k, 0)) for k in PALETTE},
            "pdac_early": int((pdac.early == 1).sum()), "pdac_late": int((pdac.early == 0).sum()),
            "missing_%": round(e.isna().mean().mean() * 100, 1),
        }
        row.update(pca_plot(name, omic, e, l))
        rows.append(row)
        print(f"[ok] {name:<10} {omic:<11} PC1 %{row['pc1_var_%']:<5} d={row['pc1_cohen_d_pdac_vs_rest']:<5} "
              f"aykırı={row['outliers'] or '-'}")
    t = pd.DataFrame(rows)
    t.to_csv(path("tables") / "wp1_cohort_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
