"""WP6-01: Kilitli miRNA panelinin (v2) bağımsız serum kohortlarında doğrulanması.

ÖN-KAYITLI ANALİZ PLANI (bu dosya çalıştırılmadan ÖNCE commit edilmiştir; K-018)
------------------------------------------------------------------------------
Panel ve skor: results/panel/mirna_panel_v2.json — DEĞİŞTİRİLMEZ.
  skor = panel miRNA'larının kohort-içi z-skorlarının ortalaması;
  z-skor parametreleri kohorttaki TÜM örneklerden, etiketsiz hesaplanır.
Kohortlar: GSE85589 (Affymetrix miRNA 4.0), GSE59856 (3D-Gene).

Birincil sonlanım:
  P1  PDAC vs sağlıklı — AUC (%95 bootstrap GA), her kohort ayrı.
İkincil sonlanımlar:
  S1  Erken evre (I-II) PDAC vs sağlıklı            [GSE59856]
  S2  Ameliyat edilebilir PDAC vs sağlıklı          [GSE59856]
  S3  Geç evre (III-IV) PDAC vs sağlıklı            [GSE59856]
  S4  PDAC vs benign pankreas/safra hastalıkları     [her iki kohort]
  S5  PDAC vs diğer kanserler                        [her iki kohort]
  S6  %90 ve %95 özgüllükte duyarlılık (PDAC vs sağlıklı)
  S7  CA19-9 kıyası [GSE85589]: CA19-9 (>37 U/mL) duyarlılık/özgüllüğü; panelin
      aynı özgüllükteki duyarlılığı; CA19-9 negatif PDAC'ta panel pozitifliği
      (%95 özgüllük eşiği); "CA19-9 VEYA panel" birleşik kuralı.
  S8  Rastgele panel kontrolü: serumda ölçülebilir miRNA havuzundan 1000 rastgele
      5'li panel (aynı skorlama); ampirik p = P(rastgele AUC ≥ panel AUC).
  S9  Karar eğrisi analizi (DCA) [GSE59856, PDAC vs sağlıklı + benign]; skor,
      kohort içi 5 katlı CV ile lojistik kalibre edilir (sıralama değişmez).
Tüm eşikler kohort içinde belirlenir (özgüllük hedefli); bu, eşik bağımsızlığı
açısından iyimserdir ve raporda belirtilir.

Çıktılar: results/tables/wp6_serum_*.tsv, results/figures/wp6_*.png
Çalıştırma:  .venv\\Scripts\\python.exe scripts/wp6_01_serum_validation.py
"""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedKFold, cross_val_predict

from fisilti.config import PROJECT_ROOT, path, set_seed
from fisilti.cv import auc_with_ci

FLOOR_59856 = 2.612
COHORTS = ["GSE85589", "GSE59856"]
N_RANDOM = 1000


def load(gse: str):
    e = pd.read_parquet(path("processed") / "geo" / f"{gse}_expr.parquet")
    l = pd.read_csv(path("processed") / "geo" / f"{gse}_labels.tsv", sep="\t", index_col=0)
    return e, l.loc[e.columns]


def zscore_rows(e: pd.DataFrame) -> pd.DataFrame:
    sd = e.std(axis=1).replace(0, np.nan)
    return e.sub(e.mean(axis=1), axis=0).div(sd, axis=0)


def panel_score(z: pd.DataFrame, members: list[str]) -> pd.Series:
    return z.loc[members].mean(axis=0)


def sens_at_spec(y: np.ndarray, s: np.ndarray, spec: float) -> tuple[float, float]:
    """Hedef özgüllüğü sağlayan en düşük eşikte duyarlılık ve eşik."""
    thr = np.quantile(s[y == 0], spec)
    return float((s[y == 1] > thr).mean()), float(thr)


def comparisons(l: pd.DataFrame) -> dict[str, tuple[pd.Series, pd.Series]]:
    pdac = l.group == "PDAC"
    out = {
        "P1_PDAC_vs_saglikli": (pdac, l.group == "healthy"),
        "S4_PDAC_vs_benign": (pdac, l.group == "benign"),
        "S5_PDAC_vs_diger_kanser": (pdac, l.group == "other_cancer"),
    }
    if l.early.notna().any():
        out["S1_erken_PDAC_vs_saglikli"] = (pdac & (l.early == 1), l.group == "healthy")
        out["S3_gec_PDAC_vs_saglikli"] = (pdac & (l.early == 0), l.group == "healthy")
    if "operable" in l and l.operable.notna().any():
        out["S2_ameliyat_edilebilir_PDAC_vs_saglikli"] = (pdac & (l.operable == 1), l.group == "healthy")
    return out


def evaluate(score: pd.Series, l: pd.DataFrame, gse: str) -> list[dict]:
    rows = []
    for name, (pos, neg) in comparisons(l).items():
        m = pos | neg
        if pos.sum() < 5 or neg.sum() < 5:
            continue
        y = pos[m].astype(int).to_numpy()
        s = score[m].to_numpy()
        auc, lo, hi = auc_with_ci(y, s)
        r = {"cohort": gse, "comparison": name, "n_pos": int(y.sum()), "n_neg": int((1 - y).sum()),
             "auc": auc, "ci_lo": lo, "ci_hi": hi}
        if name.startswith("P1"):
            r["sens_at_spec90"] = sens_at_spec(y, s, 0.90)[0]
            r["sens_at_spec95"] = sens_at_spec(y, s, 0.95)[0]
        rows.append(r)
    return rows


def random_panels(z: pd.DataFrame, l: pd.DataFrame, pool: list[str], k: int, rng) -> np.ndarray:
    m = l.group.isin(["PDAC", "healthy"])
    y = (l.group[m] == "PDAC").astype(int).to_numpy()
    zz = z.loc[pool, m[m].index].to_numpy()
    aucs = np.empty(N_RANDOM)
    for i in range(N_RANDOM):
        idx = rng.choice(len(pool), k, replace=False)
        aucs[i] = roc_auc_score(y, np.nanmean(zz[idx], axis=0))  # panel_score ile aynı kural (NaN atlanır)
    return aucs


def ca19_9(score: pd.Series, l: pd.DataFrame) -> dict:
    m = l.group.isin(["PDAC", "healthy"]) & l.ca19_9_pos.notna()
    y = (l.group[m] == "PDAC").astype(int).to_numpy()
    ca = l.ca19_9_pos[m].to_numpy()
    s = score[m].to_numpy()
    ca_sens, ca_spec = ca[y == 1].mean(), 1 - ca[y == 0].mean()
    # Panelin CA19-9 ile aynı özgüllükteki duyarlılığı
    p_sens_same_spec = sens_at_spec(y, s, min(ca_spec, 0.999))[0]
    sens95, thr95 = sens_at_spec(y, s, 0.95)
    neg_pdac = (y == 1) & (ca == 0)
    comb = (ca == 1) | (s > thr95)
    return {
        "n_pdac": int(y.sum()), "n_healthy": int((1 - y).sum()),
        "ca19_9_sens": ca_sens, "ca19_9_spec": ca_spec,
        "panel_sens_at_ca19_9_spec": p_sens_same_spec,
        "panel_auc": roc_auc_score(y, s),
        "n_ca19_9_negative_pdac": int(neg_pdac.sum()),
        "panel_pos_in_ca19_9_negative_pdac_at_spec95": float((s[neg_pdac] > thr95).mean()) if neg_pdac.any() else np.nan,
        "combined_sens": float(comb[y == 1].mean()), "combined_spec": float(1 - comb[y == 0].mean()),
    }


def dca(score: pd.Series, l: pd.DataFrame, file) -> pd.DataFrame:
    m = l.group.isin(["PDAC", "healthy", "benign"])
    y = (l.group[m] == "PDAC").astype(int).to_numpy()
    x = score[m].to_numpy().reshape(-1, 1)
    p = cross_val_predict(LogisticRegression(), x, y, cv=StratifiedKFold(5, shuffle=True, random_state=42),
                          method="predict_proba")[:, 1]
    ths = np.linspace(0.01, 0.6, 60)
    prev = y.mean()
    nb_model = [((p >= t) & (y == 1)).mean() - ((p >= t) & (y == 0)).mean() * t / (1 - t) for t in ths]
    nb_all = [prev - (1 - prev) * t / (1 - t) for t in ths]
    d = pd.DataFrame({"threshold": ths, "net_benefit_panel": nb_model, "net_benefit_treat_all": nb_all,
                      "net_benefit_treat_none": 0.0})
    fig, ax = plt.subplots(figsize=(5, 3.8), dpi=160)
    ax.plot(ths, nb_model, c="#c0392b", lw=2, label="Fısıltı miRNA paneli")
    ax.plot(ths, nb_all, c="#7f8c8d", lw=1, ls="--", label="Herkesi ileri incelemeye yönlendir")
    ax.axhline(0, c="black", lw=0.8, label="Kimseyi yönlendirme")
    ax.set_ylim(-0.05, max(prev, max(nb_model)) + 0.05)
    ax.set_xlabel("Eşik olasılığı")
    ax.set_ylabel("Net fayda")
    ax.set_title(f"Karar eğrisi · GSE59856 (prevalans %{prev * 100:.0f}, vaka-kontrol)", fontsize=9)
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(file)
    plt.close(fig)
    return d


def plot_roc(results: dict, file) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.2), dpi=160)
    for ax, gse in zip(axes, COHORTS):
        score, l = results[gse]
        for name, col in [("P1_PDAC_vs_saglikli", "#c0392b"), ("S1_erken_PDAC_vs_saglikli", "#e67e22"),
                          ("S4_PDAC_vs_benign", "#8e44ad"), ("S5_PDAC_vs_diger_kanser", "#7f8c8d")]:
            comp = comparisons(l).get(name)
            if comp is None:
                continue
            pos, neg = comp
            m = pos | neg
            if pos.sum() < 5 or neg.sum() < 5:
                continue
            y = pos[m].astype(int)
            fpr, tpr, _ = roc_curve(y, score[m])
            label = name.split("_", 1)[1].replace("_", " ")
            ax.plot(fpr, tpr, c=col, lw=1.8, label=f"{label} (AUC {roc_auc_score(y, score[m]):.2f}, n={pos.sum()}/{neg.sum()})")
        ax.plot([0, 1], [0, 1], c="#bbb", lw=0.8, ls=":")
        ax.set_title(f"{gse} · serum", fontsize=10)
        ax.set_xlabel("1 − özgüllük")
        ax.set_ylabel("Duyarlılık")
        ax.legend(fontsize=6.5, frameon=False, loc="lower right")
    fig.suptitle("Kilitli Fısıltı miRNA paneli v2 · bağımsız serum doğrulaması", fontsize=10)
    fig.tight_layout()
    fig.savefig(file)
    plt.close(fig)


def main() -> None:
    set_seed()
    rng = np.random.default_rng(42)
    lock = json.loads((PROJECT_ROOT / "results" / "panel" / "mirna_panel_v2.json").read_text(encoding="utf-8"))
    members = [m["mimat"] for m in lock["members"]]
    names = {m["mimat"]: m["name"] for m in lock["members"]}
    print(f"[..] kilitli panel: {[names[m] for m in members]} (kilit {lock['locked_utc']})")

    e59, _ = load("GSE59856")
    e85, _ = load("GSE85589")
    detectable = (e59 > FLOOR_59856 + 1e-6).mean(axis=1) >= 0.5
    pool = sorted(set(detectable[detectable].index) & set(e85.index))

    rows, singles, rand_rows, results = [], [], [], {}
    for gse in COHORTS:
        e, l = load(gse)
        missing = [m for m in members if m not in e.index]
        assert not missing, f"{gse}: panel üyesi platformda yok: {missing}"
        z = zscore_rows(e)
        score = panel_score(z, members)
        results[gse] = (score, l)
        rows += evaluate(score, l, gse)
        m = l.group.isin(["PDAC", "healthy"])
        y = (l.group[m] == "PDAC").astype(int)
        for mm in members:
            auc, lo, hi = auc_with_ci(y.to_numpy(), z.loc[mm, m].to_numpy(), n_boot=1000)
            singles.append({"cohort": gse, "mirna": names[mm], "auc_pdac_vs_healthy": auc, "ci_lo": lo, "ci_hi": hi})
        ra = random_panels(z, l, [p for p in pool if p in z.index], len(members), rng)
        panel_auc = roc_auc_score(y, score[m])
        rand_rows.append({"cohort": gse, "panel_auc": panel_auc, "random_median": float(np.median(ra)),
                          "random_p95": float(np.percentile(ra, 95)), "pool_size": len(pool),
                          "empirical_p": float((np.sum(ra >= panel_auc) + 1) / (N_RANDOM + 1))})

    res = pd.DataFrame(rows)
    tab = path("tables")
    res.round(4).to_csv(tab / "wp6_serum_primary_secondary.tsv", sep="\t", index=False)
    pd.DataFrame(singles).round(4).to_csv(tab / "wp6_serum_single_mirna.tsv", sep="\t", index=False)
    pd.DataFrame(rand_rows).round(4).to_csv(tab / "wp6_serum_random_panels.tsv", sep="\t", index=False)
    ca = ca19_9(*results["GSE85589"])
    pd.Series(ca).to_frame("value").round(4).to_csv(tab / "wp6_serum_ca19_9.tsv", sep="\t")
    d = dca(*results["GSE59856"], path("figures") / "wp6_dca_GSE59856.png")
    d.round(4).to_csv(tab / "wp6_serum_dca.tsv", sep="\t", index=False)
    plot_roc(results, path("figures") / "wp6_serum_roc.png")

    print("== Birincil + ikincil sonlanımlar")
    print(res.round(3).to_string(index=False))
    print("== Tek miRNA AUC (PDAC vs sağlıklı)")
    print(pd.DataFrame(singles).round(3).to_string(index=False))
    print("== Rastgele panel kontrolü (S8)")
    print(pd.DataFrame(rand_rows).round(3).to_string(index=False))
    print("== CA19-9 kıyası (S7, GSE85589)")
    print(pd.Series(ca).round(3).to_string())


if __name__ == "__main__":
    main()
