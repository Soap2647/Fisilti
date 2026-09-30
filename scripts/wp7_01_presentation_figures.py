"""WP7-01: Sunum/poster kalitesinde grafikler (Türkçe, büyük yazı, ortak palet).

Tüm sayılar results/ altındaki tablolardan ve işlenmiş veriden okunur; burada
hiçbir sonuç yeniden hesaplanıp değiştirilmez, yalnızca görselleştirilir.
Çıktı: results/figures/presentation/*.png
"""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve

from fisilti.config import PROJECT_ROOT, path

BG = "#F7F5F0"
INK = "#14213D"
MUTED = "#5B6475"
BLUE = "#3B6EA5"
ORANGE = "#D9693F"
GREY = "#A7AEB9"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 18, "axes.titlesize": 20, "axes.labelsize": 18,
    "xtick.labelsize": 16, "ytick.labelsize": 16, "legend.fontsize": 15,
    "axes.facecolor": BG, "figure.facecolor": BG, "savefig.facecolor": BG,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
    "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
})
OUT = path("figures") / "presentation"
OUT.mkdir(parents=True, exist_ok=True)


def zscore(e):
    return e.sub(e.mean(axis=1), axis=0).div(e.std(axis=1).replace(0, np.nan), axis=0)


def fig_external_roc():
    lock = json.loads((PROJECT_ROOT / "results/panel/protein_panel_v2.json").read_text(encoding="utf-8"))
    genes = [m["gene"] for m in lock["members"]]
    res = pd.read_csv(path("tables") / "wp6_tissue_external.tsv", sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(18, 6.2))
    for ax, gse in zip(axes, ["GSE15471", "GSE62452", "GSE16515"]):
        e = pd.read_parquet(path("processed") / "geo" / f"{gse}_expr.parquet")
        l = pd.read_csv(path("processed") / "geo" / f"{gse}_labels.tsv", sep="\t", index_col=0).loc[e.columns]
        s = zscore(e).loc[genes].mean(axis=0)
        m = l.group.isin(["PDAC", "normal"])
        fpr, tpr, _ = roc_curve(l.group[m] == "PDAC", s[m])
        r = res[(res.cohort == gse) & (res.comparison == "T1_PDAC_vs_normal")].iloc[0]
        ax.plot(fpr, tpr, c=ORANGE, lw=4, label=f"Tüm PDAC  AUC {r.auc:.2f}")
        if l.early.notna().any():
            me = (l.group == "normal") | ((l.group == "PDAC") & (l.early == 1))
            fpr, tpr, _ = roc_curve(l.group[me] == "PDAC", s[me])
            re = res[(res.cohort == gse) & (res.comparison == "T2_erken_PDAC_vs_normal")].iloc[0]
            ax.plot(fpr, tpr, c=BLUE, lw=3, ls="--", label=f"Erken evre  AUC {re.auc:.2f}")
        ax.plot([0, 1], [0, 1], c=GREY, lw=1.5, ls=":")
        ax.set_title(f"{gse}\n{int(r.n_pos)} tümör · {int(r.n_neg)} normal", color=INK)
        ax.set_xlabel("1 − özgüllük")
        ax.set_ylabel("Duyarlılık")
        ax.legend(loc="lower right", frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "external_roc.png", dpi=150)
    plt.close(fig)


def fig_mofa():
    z = pd.read_csv(path("tables") / "wp2_mofa_factors.tsv", sep="\t", index_col=0)
    l = pd.read_csv(path("processed") / "cptac_pdac" / "labels.tsv", sep="\t", index_col=0).loc[z.index]
    fig, ax = plt.subplots(figsize=(10, 7))
    for mask, col, lab in [(l.group == "normal", BLUE, "Normal doku"),
                           ((l.group == "PDAC") & (l.early == 1), ORANGE, "PDAC erken evre (I–II)"),
                           ((l.group == "PDAC") & (l.early == 0), "#7A2E14", "PDAC geç evre (III–IV)")]:
        ax.scatter(z.loc[mask, "F1"], z.loc[mask, "F2"], s=60, c=col, alpha=0.8, edgecolors="none",
                   label=f"{lab} (n={mask.sum()})")
    ax.set_xlabel("F1 · asiner sindirim enzimleri (protein, metilasyon)")
    ax.set_ylabel("F2 · asiner program (mRNA, miRNA)")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT / "mofa_factors.png", dpi=150)
    plt.close(fig)


def fig_tissue_auc():
    s = pd.read_csv(path("tables") / "wp3_tabular_summary.tsv", sep="\t")
    s = s[s.config == "blood"].sort_values("auc", ascending=False).drop_duplicates("view")
    label = {"protein": "Protein", "mrna": "mRNA", "methylation": "Metilasyon", "mirna": "miRNA"}
    s = s.set_index("view").loc[["protein", "mrna", "methylation", "mirna"]]
    fig, ax = plt.subplots(figsize=(10, 6))
    y = np.arange(len(s))[::-1]
    ax.errorbar(s.auc, y, xerr=[s.auc - s.auc_ci_lo, s.auc_ci_hi - s.auc], fmt="o", c=ORANGE, ms=14,
                lw=3, capsize=8, label="Tüm PDAC vs normal")
    ax.scatter(s.auc_early_vs_normal, y - 0.22, marker="D", s=110, c=BLUE, label="Erken evre vs normal", zorder=3)
    ax.set_yticks(y, [label[v] for v in s.index])
    ax.set_xlim(0.85, 1.0)
    ax.set_xlabel("AUC (5×5 iç içe, hasta gruplu çapraz doğrulama)")
    ax.legend(frameon=False, loc="lower left")
    fig.tight_layout()
    fig.savefig(OUT / "tissue_auc.png", dpi=150)
    plt.close(fig)


def fig_serum():
    ex = pd.read_csv(path("tables") / "wp6_exploratory.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(18, 6.8), gridspec_kw={"width_ratios": [1, 1.15]})
    ax = axes[0]
    notes = ex[ex.analysis == "E2_serum_tum_ortak"]
    # Kohort içi CV AUC, wp6_02'nin not sütununda "eğitim içi CV AUC x.xxx" olarak kayıtlı
    within = notes.note.str.extract(r"eğitim içi CV AUC ([0-9.]+)", expand=False).astype(float).tolist()
    cross = notes.value.tolist()
    names = ["GSE59856 → GSE85589", "GSE85589 → GSE59856"]
    x = np.arange(2)
    ax.bar(x - 0.2, within, 0.38, color=GREY, label="Kendi kohortu içinde")
    ax.bar(x + 0.2, cross, 0.38, color=ORANGE, label="Diğer kohortta")
    for i in range(2):
        ax.text(x[i] - 0.2, within[i] + 0.02, f"{within[i]:.2f}", ha="center", fontsize=16)
        ax.text(x[i] + 0.2, cross[i] + 0.02, f"{cross[i]:.2f}", ha="center", fontsize=16, color=ORANGE)
    ax.axhline(0.5, c=MUTED, ls=":", lw=1.5)
    ax.text(0.5, 0.52, "şans düzeyi", color=MUTED, fontsize=14, ha="center")
    ax.set_xticks(x, names)
    ax.set_ylim(0, 1.32)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_ylabel("AUC (PDAC vs sağlıklı)")
    ax.set_title("Serum miRNA modeli başka kohorta taşınıyor mu?")
    ax.legend(frameon=False, loc="upper center", ncol=2)

    ax = axes[1]
    proc = path("processed") / "cptac_pdac"
    lab = pd.read_csv(proc / "labels.tsv", sep="\t", index_col=0)
    tis = pd.read_parquet(proc / "mirna.parquet")
    g = lab.loc[tis.columns, "group"]
    fc = tis.loc[:, g == "PDAC"].median(axis=1) - tis.loc[:, g == "normal"].median(axis=1)
    e = pd.read_parquet(path("processed") / "geo" / "GSE59856_expr.parquet")
    l = pd.read_csv(path("processed") / "geo" / "GSE59856_labels.tsv", sep="\t", index_col=0).loc[e.columns]
    det = (e > 2.612 + 1e-6).mean(axis=1) >= 0.5
    feats = [m for m in fc.index if m in e.index and det.get(m, False)]
    a, b = e.loc[feats, l.group == "PDAC"], e.loc[feats, l.group == "healthy"]
    d = (a.mean(axis=1) - b.mean(axis=1)) / np.sqrt((a.var(axis=1) + b.var(axis=1)) / 2)
    ax.scatter(fc[feats], d, s=40, c=BLUE, alpha=0.6, edgecolors="none")
    panel = [m["mimat"] for m in json.loads((PROJECT_ROOT / "results/panel/mirna_panel_v2.json").read_text(encoding="utf-8"))["members"]]
    pm = [p for p in panel if p in feats]
    ax.scatter(fc[pm], d[pm], s=160, c=ORANGE, edgecolors=INK, lw=1, label="Kilitli panel", zorder=3)
    ax.axhline(0, c=GREY, lw=1)
    ax.axvline(0, c=GREY, lw=1)
    rho = ex[(ex.analysis == "E1_doku_serum_uyumu") & (ex.cohort == "GSE59856")].value.iloc[0]
    ax.set_title(f"Dokuda artan, kanda artıyor mu?  ρ = {rho:.2f}")
    ax.set_xlabel("Doku: log2 kat değişimi (tümör − normal)")
    ax.set_ylabel("Serum: etki büyüklüğü (PDAC − sağlıklı)")
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT / "serum_failure.png", dpi=150)
    plt.close(fig)


def fig_gsea():
    rows = []
    for view in ["protein", "mrna"]:
        t = pd.read_csv(path("tables") / f"wp4_gsea_{view}.tsv", sep="\t")
        t = t[(t.library == "MSigDB_Hallmark_2020") & (t["FDR q-val"] < 0.05)]
        rows.append(t.assign(view=view))
    d = pd.concat(rows)
    keep = ["KRAS Signaling Up", "Epithelial Mesenchymal Transition", "Coagulation", "E2F Targets",
            "G2-M Checkpoint", "Interferon Gamma Response", "Hypoxia", "p53 Pathway",
            "Oxidative Phosphorylation", "Pancreas Beta Cells"]
    tr = {"KRAS Signaling Up": "KRAS sinyali", "Epithelial Mesenchymal Transition": "Epitel-mezenkimal geçiş",
          "Coagulation": "Pıhtılaşma", "E2F Targets": "Hücre çoğalması (E2F)", "G2-M Checkpoint": "Hücre döngüsü G2-M",
          "Interferon Gamma Response": "İnterferon-γ yanıtı", "Hypoxia": "Hipoksi", "p53 Pathway": "p53 yolağı",
          "Oxidative Phosphorylation": "Oksidatif fosforilasyon", "Pancreas Beta Cells": "Pankreas beta hücreleri"}
    d = d[d.Term.isin(keep)]
    order = [k for k in keep if k in set(d.Term)][::-1]
    fig, ax = plt.subplots(figsize=(11, 7))
    for view, col, lab, off in [("protein", ORANGE, "Protein", 0.15), ("mrna", BLUE, "mRNA", -0.15)]:
        s = d[d.view == view].set_index("Term").reindex(order)
        ax.scatter(s.NES, np.arange(len(order)) + off, s=170, c=col, label=lab, zorder=3)
    ax.axvline(0, c=GREY, lw=1.2)
    ax.set_yticks(range(len(order)), [tr[o] for o in order])
    ax.set_xlabel("NES  ·  sağ = tümörde aktif, sol = baskılanmış")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT / "gsea.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    for f in [fig_external_roc, fig_mofa, fig_tissue_auc, fig_serum, fig_gsea]:
        f()
        print("[ok]", f.__name__)
