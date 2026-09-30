"""WP6-03: Kilitli protein aday paneli v2'nin bağımsız doku kohortlarında dış doğrulaması.

ÖN-KAYITLI ANALİZ PLANI (bu dosya çalıştırılmadan ÖNCE commit edilmiştir; K-021)
------------------------------------------------------------------------------
Panel: results/panel/protein_panel_v2.json — DEĞİŞTİRİLMEZ (8 gen, tümü 'up').
Ölçüm: panel genlerinin mRNA ifadesi (farklı platform ve kohortlar; CPTAC'tan bağımsız).
Skor: panel genlerinin kohort-içi z-skor ortalaması (etiketsiz, tüm örnekler).
Kohortlar: GSE15471 (U133 Plus 2, eşleşmiş T/N), GSE62452 (HuGene 1.0 ST),
           GSE16515 (U133 Plus 2).

Birincil sonlanım:
  T1  PDAC vs normal doku AUC (%95 bootstrap GA), her kohort ayrı.
İkincil:
  T2  Erken evre (I-II) PDAC vs normal [GSE62452]
  T3  Geç evre (III-IV) PDAC vs normal [GSE62452]
  T4  Tek gen AUC'leri (yön tutarlılığı: kaç gen AUC > 0.5)
  T5  Rastgele 8'li gen paneli kontrolü (platformdaki tüm genler, 1000 panel;
      ampirik p). Not: dokuda birçok gen tümör/normali ayırdığından bu kontrol
      katıdır ve panelin rastgele genlerden daha iyi olup olmadığını sınar.
Yorum sınırı: Bu doğrulama DOKU düzeyindedir; kan testi performansını göstermez.

Çıktılar: results/tables/wp6_tissue_external.tsv, results/figures/wp6_tissue_external_roc.png
Çalıştırma:  .venv\\Scripts\\python.exe scripts/wp6_03_tissue_external.py
"""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve

from fisilti.config import PROJECT_ROOT, path, set_seed
from fisilti.cv import auc_with_ci

COHORTS = ["GSE15471", "GSE62452", "GSE16515"]
N_RANDOM = 1000


def main() -> None:
    set_seed()
    rng = np.random.default_rng(42)
    lock = json.loads((PROJECT_ROOT / "results" / "panel" / "protein_panel_v2.json").read_text(encoding="utf-8"))
    genes = [m["gene"] for m in lock["members"]]
    print(f"[..] kilitli protein paneli v2: {genes}")
    rows, singles = [], []
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), dpi=160)
    for ax, gse in zip(axes, COHORTS):
        e = pd.read_parquet(path("processed") / "geo" / f"{gse}_expr.parquet")
        l = pd.read_csv(path("processed") / "geo" / f"{gse}_labels.tsv", sep="\t", index_col=0).loc[e.columns]
        present = [g for g in genes if g in e.index]
        z = e.sub(e.mean(axis=1), axis=0).div(e.std(axis=1).replace(0, np.nan), axis=0)
        score = z.loc[present].mean(axis=0)
        y = (l.group == "PDAC").astype(int)
        comps = {"T1_PDAC_vs_normal": (l.group == "PDAC", l.group == "normal")}
        if l.early.notna().any():
            comps["T2_erken_PDAC_vs_normal"] = ((l.group == "PDAC") & (l.early == 1), l.group == "normal")
            comps["T3_gec_PDAC_vs_normal"] = ((l.group == "PDAC") & (l.early == 0), l.group == "normal")
        for name, (pos, neg) in comps.items():
            m = pos | neg
            yy = pos[m].astype(int).to_numpy()
            auc, lo, hi = auc_with_ci(yy, score[m].to_numpy())
            r = {"cohort": gse, "comparison": name, "n_pos": int(yy.sum()), "n_neg": int((1 - yy).sum()),
                 "genes_on_platform": f"{len(present)}/{len(genes)}", "auc": auc, "ci_lo": lo, "ci_hi": hi}
            if name.startswith("T1"):
                zz = z.dropna().to_numpy()
                ra = np.array([roc_auc_score(y, zz[rng.choice(len(zz), len(present), replace=False)].mean(0))
                               for _ in range(N_RANDOM)])
                r["random_median"] = float(np.median(ra))
                r["random_p95"] = float(np.percentile(ra, 95))
                r["empirical_p"] = float((np.sum(ra >= auc) + 1) / (N_RANDOM + 1))
                fpr, tpr, _ = roc_curve(yy, score[m])
                ax.plot(fpr, tpr, c="#6c3483", lw=2, label=f"Tüm PDAC (AUC {auc:.2f}; n={yy.sum()}/{(1 - yy).sum()})")
            elif name.startswith("T2"):
                fpr, tpr, _ = roc_curve(yy, score[m])
                ax.plot(fpr, tpr, c="#e67e22", lw=1.6, label=f"Erken evre (AUC {auc:.2f}; n={yy.sum()})")
            rows.append(r)
        for g in present:
            singles.append({"cohort": gse, "gene": g, "auc": roc_auc_score(y, z.loc[g].fillna(0))})
        ax.plot([0, 1], [0, 1], c="#bbb", lw=0.8, ls=":")
        ax.set_title(f"{gse} · doku mRNA", fontsize=10)
        ax.set_xlabel("1 − özgüllük")
        ax.set_ylabel("Duyarlılık")
        ax.legend(fontsize=7, frameon=False, loc="lower right")
    fig.suptitle("Kilitli protein aday paneli v2 · bağımsız doku kohortlarında dış doğrulama", fontsize=10)
    fig.tight_layout()
    fig.savefig(path("figures") / "wp6_tissue_external_roc.png")
    plt.close(fig)
    t = pd.DataFrame(rows)
    t.round(4).to_csv(path("tables") / "wp6_tissue_external.tsv", sep="\t", index=False)
    s = pd.DataFrame(singles).pivot(index="gene", columns="cohort", values="auc")
    s.round(4).to_csv(path("tables") / "wp6_tissue_external_single_genes.tsv", sep="\t")
    print(t.round(3).to_string(index=False))
    print("== Tek gen AUC (PDAC vs normal)")
    print(s.round(3).to_string())


if __name__ == "__main__":
    main()
