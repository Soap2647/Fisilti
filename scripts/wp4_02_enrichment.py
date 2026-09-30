"""WP4-02: Yolak (pathway) zenginleştirme analizi.

1. Ön-sıralı GSEA (CPTAC protein ve mRNA; tümör vs normal Welch t istatistiği ile
   tüm genler sıralanır) — MSigDB Hallmark 2020, KEGG 2021, Reactome 2022.
2. Aday genlerde aşırı temsil (ORA; hipergeometrik) — arka plan = o omikte test
   edilen TÜM genler (yanlış arka plan kullanımı yaygın bir hatadır).
   Aday = WP3 kararlılık seçiminde sel_freq ≥ 0.5 (blood senaryosu).
3. Aday miRNA'ların hedefleri (miRTarBase 2017): en az 2 aday miRNA tarafından
   hedeflenen genlerde ORA (keşif amaçlı).

Gen seti kütüphaneleri ilk çalıştırmada data/raw/genesets/ altına JSON olarak
kaydedilir; sonraki çalıştırmalar yerel kopyayı kullanır (tekrarlanabilirlik).

Çıktılar: results/tables/wp4_gsea_<omik>.tsv, wp4_ora_<küme>.tsv
          results/figures/wp4_gsea_hallmark.png
Çalıştırma:  .venv\\Scripts\\python.exe scripts/wp4_02_enrichment.py
"""

from __future__ import annotations

import json
import warnings
from datetime import date

import matplotlib

matplotlib.use("Agg")
import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

from fisilti.config import path, set_seed

LIBS = ["MSigDB_Hallmark_2020", "KEGG_2021_Human", "Reactome_2022", "GO_Biological_Process_2023"]
GSEA_LIBS = LIBS[:3]


def library(name: str) -> dict[str, list[str]]:
    d = path("raw") / "genesets"
    d.mkdir(parents=True, exist_ok=True)
    f = d / f"{name}.json"
    if not f.exists():
        lib = gp.get_library(name, organism="Human")
        f.write_text(json.dumps({"_retrieved": str(date.today()), "sets": lib}), encoding="utf-8")
    return json.loads(f.read_text(encoding="utf-8"))["sets"]


def ranking(view: str) -> pd.Series:
    proc = path("processed") / "cptac_pdac"
    lab = pd.read_csv(proc / "labels.tsv", sep="\t", index_col=0)
    m = pd.read_parquet(proc / f"{view}.parquet")
    m = m.loc[m.notna().mean(axis=1) >= 0.7]
    g = lab.loc[m.columns, "group"]
    t = ttest_ind(m.loc[:, g == "PDAC"], m.loc[:, g == "normal"], axis=1, equal_var=False, nan_policy="omit").statistic
    r = pd.Series(np.asarray(t), index=m.index).dropna()
    return r[~r.index.duplicated()].sort_values(ascending=False)


def run_gsea(view: str) -> pd.DataFrame:
    rnk = ranking(view)
    out = []
    for lib in GSEA_LIBS:
        res = gp.prerank(rnk=rnk, gene_sets=library(lib), permutation_num=1000, min_size=15, max_size=500,
                         seed=42, threads=8, outdir=None, verbose=False).res2d
        res["library"] = lib
        out.append(res)
    t = pd.concat(out)[["library", "Term", "NES", "NOM p-val", "FDR q-val", "Lead_genes"]]
    t = t.astype({"NES": float, "NOM p-val": float, "FDR q-val": float}).sort_values("FDR q-val")
    t.to_csv(path("tables") / f"wp4_gsea_{view}.tsv", sep="\t", index=False)
    return t


def run_ora(genes: list[str], background: list[str], tag: str) -> pd.DataFrame:
    out = []
    for lib in LIBS:
        res = gp.enrich(gene_list=genes, gene_sets=library(lib), background=background, outdir=None,
                        cutoff=1.0, verbose=False).results
        res["library"] = lib
        out.append(res)
    t = pd.concat(out)[["library", "Term", "Overlap", "P-value", "Adjusted P-value", "Genes"]]
    t = t.sort_values("Adjusted P-value")
    t.to_csv(path("tables") / f"wp4_ora_{tag}.tsv", sep="\t", index=False)
    return t


def plot_hallmark(gsea: dict[str, pd.DataFrame], file) -> None:
    rows = []
    for view, t in gsea.items():
        h = t[(t.library == "MSigDB_Hallmark_2020") & (t["FDR q-val"] < 0.05)]
        rows.append(h.assign(view=view))
    d = pd.concat(rows)
    terms = (d.groupby("Term").NES.apply(lambda s: s.abs().max()).sort_values(ascending=False).index[:18])
    d = d[d.Term.isin(terms)]
    order = d.groupby("Term").NES.mean().sort_values().index
    fig, ax = plt.subplots(figsize=(6.4, 0.3 * len(order) + 1.3), dpi=160)
    for i, (view, col) in enumerate([("protein", "#6c3483"), ("mrna", "#1f618d")]):
        s = d[d.view == view].set_index("Term").reindex(order)
        ax.scatter(s.NES, np.arange(len(order)) + (i - 0.5) * 0.25, c=col, s=28,
                   label="Protein" if view == "protein" else "mRNA")
    ax.axvline(0, c="#888", lw=0.8)
    ax.set_yticks(range(len(order)), [o.replace("_", " ") for o in order], fontsize=7)
    ax.set_xlabel("NES (pozitif = tümörde aktif)")
    ax.set_title("GSEA · MSigDB Hallmark · CPTAC tümör vs normal (FDR < 0,05)", fontsize=9)
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(file)
    plt.close(fig)


def main() -> None:
    warnings.filterwarnings("ignore")
    set_seed()
    gsea = {}
    for view in ["protein", "mrna"]:
        t = run_gsea(view)
        gsea[view] = t
        h = t[t.library == "MSigDB_Hallmark_2020"].head(8)
        print(f"== GSEA {view}: FDR<0,05 terim sayısı {(t['FDR q-val'] < 0.05).sum()}")
        print(h[["Term", "NES", "FDR q-val"]].round(4).to_string(index=False))
    plot_hallmark(gsea, path("figures") / "wp4_gsea_hallmark.png")

    # Aday gen ORA (protein + mRNA ayrı)
    for view in ["protein", "mrna"]:
        cand = pd.read_csv(path("tables") / f"wp3_candidates_{view}.tsv", sep="\t", index_col=0)
        genes = list(cand.index[cand.sel_freq >= 0.5])
        m = pd.read_parquet(path("processed") / "cptac_pdac" / f"{view}.parquet")
        bg = list(m.index[m.notna().mean(axis=1) >= 0.7])
        t = run_ora(genes, bg, f"{view}_candidates")
        print(f"== ORA {view} adayları (n={len(genes)}, arka plan {len(bg)}): FDR<0,05 terim {(t['Adjusted P-value'] < 0.05).sum()}")
        print(t.head(6)[["library", "Term", "Overlap", "Adjusted P-value"]].round(4).to_string(index=False))

    # miRNA hedefleri
    cand = pd.read_csv(path("tables") / "wp3_candidates_mirna.tsv", sep="\t", index_col=0)
    mirs = list(cand.loc[cand.sel_freq >= 0.5, "name"])
    mt = library("miRTarBase_2017")
    hits = pd.Series([g for mir in mirs for g in mt.get(mir, [])]).value_counts()
    targets = list(hits.index[hits >= 2])
    bg = sorted({g for gs in mt.values() for g in gs})
    covered = [m for m in mirs if m in mt]
    print(f"== miRNA adayları {len(mirs)} (miRTarBase'de bulunan {len(covered)}); ≥2 aday tarafından hedeflenen gen {len(targets)}")
    t = run_ora(targets, bg, "mirna_targets")
    print(t.head(8)[["library", "Term", "Overlap", "Adjusted P-value"]].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
