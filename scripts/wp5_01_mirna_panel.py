"""WP5-01: Serum RT-qPCR'a yönelik miRNA panelinin tasarımı ve KİLİTLENMESİ.

İlke (K-016): Panel üyeleri ve skorlama kuralı serum kohortlarının ETİKETLERİNE
bakılmadan belirlenir. Serum verisi yalnızca etiketsiz bir ölçülebilirlik
kontrolü için kullanılır. Kilitlenen panel git'e işlenir; WP6 doğrulaması bu
kilitli dosyayı değiştirmeden kullanır.

Adımlar
  1. Doku adayları: WP3 kararlılık seçimi (blood), sel_freq ≥ 0.5
  2. Asiner program güvenlik filtresi: tümörlerde ρ(asiner) < 0.3
  3. Kan hücresi / hemoliz kökenli miRNA'ların dışlanması (literatür listesi)
  4. Serumda ölçülebilirlik (etiketsiz): GSE59856 (3D-Gene) örneklerinin ≥ %50'sinde
     tespit sınırı üzerinde; GSE85589 platformunda bulunma
  5. Fazlalık giderme: dokuda |Spearman r| > 0.7 olan daha düşük sıralı üye atılır
  6. Panel boyutu (5–10): dokuda hasta gruplu CV ile L2-lojistik AUC eğrisi;
     en yüksek AUC'nin 0.005 altına ilk ulaşılan boyut (en az 5)
Skorlama kuralı (kilitli): panel miRNA'larının kohort-içi z-skorlarının ortalaması
(tüm üyeler tümörde artan yönde). Serumda hiçbir ağırlık öğrenilmez.

Çıktılar: results/panel/mirna_panel_v1.json, results/tables/wp5_mirna_panel_steps.tsv
Çalıştırma:  .venv\\Scripts\\python.exe scripts/wp5_01_mirna_panel.py
"""

from __future__ import annotations

import json
import subprocess
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from fisilti.config import PROJECT_ROOT, path, set_seed
from fisilti.cv import grouped_folds

# Kan hücresi kökenli / hemolize duyarlı miRNA'lar
# Alyuvar: Kirschner ve ark. 2011 PLoS One; Blondal ve ark. 2013 Methods
# Lökosit/trombosit: Pritchard ve ark. 2012 Cancer Prev Res
BLOOD_CELL_MIRNAS = {
    "hsa-miR-451a", "hsa-miR-16-5p", "hsa-miR-486-5p", "hsa-miR-92a-3p", "hsa-miR-144-3p",
    "hsa-miR-144-5p", "hsa-miR-15b-5p", "hsa-miR-150-5p", "hsa-miR-142-3p", "hsa-miR-142-5p",
    "hsa-miR-223-3p", "hsa-miR-126-3p",
}
FLOOR_59856 = 2.612  # 3D-Gene tespit sınırı değeri (tüm değerlerin ~%64'ü)
SEL_MIN, RHO_ACINAR_MAX, DETECT_MIN, REDUND_R = 0.5, 0.3, 0.5, 0.7
SIZE_MIN, SIZE_MAX = 5, 10


def git_commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, text=True).strip()
    except Exception:
        return "bilinmiyor"


def tissue_cv_auc(X: pd.DataFrame, y: np.ndarray, groups: np.ndarray, members: list[str]) -> float:
    aucs = []
    for _, _, tr, te in grouped_folds(y, groups, n_splits=5, repeats=5):
        m = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, class_weight="balanced", max_iter=5000))
        m.fit(X.iloc[tr][members], y[tr])
        aucs.append(roc_auc_score(y[te], m.predict_proba(X.iloc[te][members])[:, 1]))
    return float(np.mean(aucs))


def main() -> None:
    set_seed()
    tab = path("tables")
    cand = pd.read_csv(tab / "wp3_candidates_mirna.tsv", sep="\t", index_col=0)
    shap = pd.read_csv(tab / "wp4_shap_mirna.tsv", sep="\t", index_col=0)
    cand["mean_abs_shap"] = shap.mean_abs_shap.reindex(cand.index).fillna(0)
    steps = cand[["name", "sel_freq", "mean_abs_shap", "auc_all", "auc_early", "log2fc", "rho_acinar"]].copy()

    # 4. Serumda ölçülebilirlik (etiketsiz)
    e59 = pd.read_parquet(path("processed") / "geo" / "GSE59856_expr.parquet")
    e85 = pd.read_parquet(path("processed") / "geo" / "GSE85589_expr.parquet")
    steps["detect_frac_GSE59856"] = (e59.reindex(steps.index) > FLOOR_59856 + 1e-6).mean(axis=1)
    steps["on_GSE85589"] = steps.index.isin(e85.index)
    steps["median_pct_GSE85589"] = e85.median(axis=1).rank(pct=True).reindex(steps.index)

    steps["pass_sel"] = steps.sel_freq >= SEL_MIN
    steps["pass_acinar"] = steps.rho_acinar < RHO_ACINAR_MAX
    steps["pass_bloodcell"] = ~steps.name.isin(BLOOD_CELL_MIRNAS)
    steps["pass_serum"] = (steps.detect_frac_GSE59856 >= DETECT_MIN) & steps.on_GSE85589
    ok = steps.pass_sel & steps.pass_acinar & steps.pass_bloodcell & steps.pass_serum
    ranked = steps[ok].sort_values(["sel_freq", "mean_abs_shap"], ascending=False)

    # 5. Fazlalık giderme (doku korelasyonu)
    proc = path("processed") / "cptac_pdac"
    lab = pd.read_csv(proc / "labels.tsv", sep="\t", index_col=0)
    X = pd.read_parquet(proc / "mirna.parquet").T
    lab = lab.loc[X.index]
    y = (lab.group == "PDAC").astype(int).to_numpy()
    groups = lab.patient.to_numpy()
    corr = X[ranked.index].corr(method="spearman")
    kept, dropped_for = [], {}
    for m in ranked.index:
        hit = [k for k in kept if abs(corr.at[m, k]) > REDUND_R]
        if hit:
            dropped_for[m] = hit[0]
        else:
            kept.append(m)
    steps["redundant_with"] = pd.Series({m: steps.at[k, "name"] for m, k in dropped_for.items()})

    # 6. Panel boyutu: doku CV AUC eğrisi
    curve = {n: tissue_cv_auc(X, y, groups, kept[:n]) for n in range(1, min(SIZE_MAX, len(kept)) + 1)}
    best = max(curve.values())
    size = next(n for n in sorted(curve) if n >= SIZE_MIN and curve[n] >= best - 0.005) \
        if len(kept) >= SIZE_MIN else len(kept)
    panel = kept[:size]
    steps["in_panel"] = steps.index.isin(panel)
    steps.round(4).to_csv(tab / "wp5_mirna_panel_steps.tsv", sep="\t", index_label="mimat")

    print("== Filtre hunisi:")
    print(f"   aday (sel_freq ≥ {SEL_MIN}): {steps.pass_sel.sum()}")
    f = steps.pass_sel
    for col, label in [("pass_acinar", "asiner güvenliği"), ("pass_bloodcell", "kan hücresi dışlama"),
                       ("pass_serum", "serumda ölçülebilir")]:
        f = f & steps[col]
        print(f"   + {label}: {f.sum()}")
    print(f"   + fazlalık giderme: {len(kept)}  (atılan: "
          f"{ {steps.at[m, 'name']: steps.at[k, 'name'] for m, k in dropped_for.items()} })")
    print("== Doku CV AUC eğrisi (panel boyutu -> AUC): " + ", ".join(f"{n}:{a:.3f}" for n, a in curve.items()))
    print(f"== KİLİTLİ PANEL ({size} miRNA):")
    print(steps.loc[panel, ["name", "sel_freq", "mean_abs_shap", "auc_early", "log2fc",
                            "detect_frac_GSE59856", "median_pct_GSE85589"]].round(3).to_string())

    lock = {
        "panel": "Fısıltı miRNA paneli v1",
        "locked_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "code_commit_at_lock": git_commit(),
        "members": [{"mimat": m, "name": steps.at[m, "name"], "direction": "up"} for m in panel],
        "scoring_rule": "kohort-içi z-skor (tüm örnekler, etiketsiz) ortalaması; yüksek skor = PDAC",
        "selection_used_serum_labels": False,
        "criteria": {"sel_freq_min": SEL_MIN, "rho_acinar_max": RHO_ACINAR_MAX,
                     "serum_detect_frac_min_GSE59856": DETECT_MIN, "floor_GSE59856": FLOOR_59856,
                     "redundancy_abs_spearman": REDUND_R, "size_range": [SIZE_MIN, SIZE_MAX],
                     "blood_cell_exclusion": sorted(BLOOD_CELL_MIRNAS)},
        "tissue_cv_auc_curve": {str(k): round(v, 4) for k, v in curve.items()},
    }
    out = PROJECT_ROOT / "results" / "panel"
    out.mkdir(parents=True, exist_ok=True)
    (out / "mirna_panel_v1.json").write_text(json.dumps(lock, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[ok] kilit dosyası: results/panel/mirna_panel_v1.json")


if __name__ == "__main__":
    main()
