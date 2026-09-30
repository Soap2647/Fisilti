"""WP5-02: Plazma ELISA'ya yönelik protein aday panelinin tasarımı.

Kanda doğrulama verisi bulunmadığından bu panel "aday panel" olarak raporlanır
(K-017). Salgılanma/plazma kanıtı iki bağımsız kaynaktan toplanır:
  - Human Protein Atlas (HPA): salgılanma konumu (secretome location), insan
    kanında ölçülmüş konsantrasyon (kütle spektrometrisi = MS, immünoanaliz = IM)
  - UniProt (Swiss-Prot): sinyal peptidi, hücre-altı konum ('Secreted')
Başvurudaki SignalP 6.0 / SecretomeP yerine bu kaynaklar kullanıldı: UniProt
sinyal peptidi anotasyonları deneysel + öngörü kanıtını, HPA salgılanma sınıfı ise
klasik dışı salgılanmayı da kapsar; ikisi de programatik ve tekrarlanabilirdir.

Düzey 1: (HPA 'Secreted to blood' VEYA UniProt Secreted + sinyal peptidi) VE kanda ölçülmüş
Düzey 2: kanda ölçülmüş (MS veya IM) ancak klasik salgılanma kanıtı yok
Dışlanan: kanda ölçüm kanıtı yok
Panel: Düzey 1 > Düzey 2; IM (immünoanaliz) ölçümü olanlar öne (ELISA uygulanabilirliği);
       fazlalık giderme ve doku CV ile boyut (5–10), miRNA paneliyle aynı kurallar.

Çıktılar: results/panel/protein_panel_v2.json, results/tables/wp5_protein_annotation.tsv
          data/raw/annotations/{hpa,uniprot}/<gen>.json (ham yanıtlar)
Çalıştırma:  .venv\\Scripts\\python.exe scripts/wp5_02_protein_panel.py
"""

from __future__ import annotations

import json
import sys
import time
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import requests

from fisilti.config import PROJECT_ROOT, path, set_seed

sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from wp5_01_mirna_panel import REDUND_R, RHO_ACINAR_MAX, SEL_MIN, SIZE_MAX, SIZE_MIN, git_commit, tissue_cv_auc  # noqa: E402

HPA_COLS = "g,eg,up,scl,secl,blconcim,blconcms,pe"


def cached_json(kind: str, gene: str, fetch) -> object:
    d = path("raw") / "annotations" / kind
    d.mkdir(parents=True, exist_ok=True)
    f = d / f"{gene}.json"
    if f.exists():
        return json.loads(f.read_text(encoding="utf-8"))
    data = fetch()
    f.write_text(json.dumps({"_retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                             "data": data}, ensure_ascii=False), encoding="utf-8")
    time.sleep(0.2)
    return {"data": data}


def hpa(gene: str) -> dict:
    def fetch():
        r = requests.get("https://www.proteinatlas.org/api/search_download.php",
                         params={"search": gene, "format": "json", "columns": HPA_COLS, "compress": "no"}, timeout=60)
        r.raise_for_status()
        return [x for x in r.json() if x.get("Gene") == gene]
    d = cached_json("hpa", gene, fetch)["data"]
    return d[0] if d else {}


def uniprot(gene: str) -> dict:
    def fetch():
        r = requests.get("https://rest.uniprot.org/uniprotkb/search",
                         params={"query": f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true",
                                 "fields": "accession,gene_primary,ft_signal,cc_subcellular_location",
                                 "format": "json"}, timeout=60)
        r.raise_for_status()
        return r.json().get("results", [])
    d = cached_json("uniprot", gene, fetch)["data"]
    return d[0] if d else {}


def annotate(gene: str) -> dict:
    h, u = hpa(gene), uniprot(gene)
    locs = [l["location"]["value"] for c in u.get("comments", []) if c.get("commentType") == "SUBCELLULAR LOCATION"
            for l in c.get("subcellularLocations", [])]
    signal = any(f.get("type") == "Signal" for f in u.get("features", []))
    secl = h.get("Secretome location") or ""
    ms = h.get("Blood concentration - Conc. blood MS [pg/L]")
    im = h.get("Blood concentration - Conc. blood IM [pg/L]")
    return {
        "uniprot": u.get("primaryAccession"),
        "hpa_secretome": secl,
        "uniprot_secreted": any(l.startswith("Secreted") for l in locs),
        "signal_peptide": signal,
        "blood_ms_pg_per_L": float(ms) if ms else np.nan,
        "blood_im_pg_per_L": float(im) if im else np.nan,
    }


def main() -> None:
    set_seed()
    tab = path("tables")
    cand = pd.read_csv(tab / "wp3_candidates_protein.tsv", sep="\t", index_col=0)
    shap = pd.read_csv(tab / "wp4_shap_protein.tsv", sep="\t", index_col=0)
    cand["mean_abs_shap"] = shap.mean_abs_shap.reindex(cand.index).fillna(0)
    c = cand[(cand.sel_freq >= SEL_MIN) & (cand.rho_acinar < RHO_ACINAR_MAX)].copy()
    ann = pd.DataFrame({g: annotate(g) for g in c.index}).T
    c = c.join(ann)
    in_blood = c.blood_ms_pg_per_L.notna() | c.blood_im_pg_per_L.notna()
    secreted = c.hpa_secretome.str.contains("blood", case=False) | (c.uniprot_secreted & c.signal_peptide)
    c["tier"] = np.select([secreted & in_blood, in_blood], ["1", "2"], "dışlandı")
    c["elisa_evidence_IM"] = c.blood_im_pg_per_L.notna()
    c.round(4).to_csv(tab / "wp5_protein_annotation.tsv", sep="\t", index_label="gene")

    ranked = c[c.tier != "dışlandı"].sort_values(
        ["tier", "elisa_evidence_IM", "sel_freq", "mean_abs_shap"], ascending=[True, False, False, False])
    proc = path("processed") / "cptac_pdac"
    lab = pd.read_csv(proc / "labels.tsv", sep="\t", index_col=0)
    X = pd.read_parquet(proc / "protein.parquet").T
    lab = lab.loc[X.index]
    X = X.fillna(X.median())
    y = (lab.group == "PDAC").astype(int).to_numpy()
    groups = lab.patient.to_numpy()
    corr = X.loc[y == 1, ranked.index].corr(method="spearman")  # yalnızca tümörler (K-017)
    kept = []
    for g in ranked.index:
        if not any(abs(corr.at[g, k]) > REDUND_R for k in kept):
            kept.append(g)
    curve = {n: tissue_cv_auc(X, y, groups, kept[:n]) for n in range(1, min(SIZE_MAX, len(kept)) + 1)}
    best = max(curve.values())
    size = next((n for n in sorted(curve) if n >= SIZE_MIN and curve[n] >= best - 0.005), len(kept))
    panel = kept[:size]

    print(f"== {len(c)} protein adayı: Düzey 1 {int((c.tier == '1').sum())}, Düzey 2 {int((c.tier == '2').sum())}, "
          f"dışlanan {int((c.tier == 'dışlandı').sum())}")
    print(c.sort_values(["tier", "sel_freq"], ascending=[True, False])[
        ["tier", "sel_freq", "hpa_secretome", "signal_peptide", "blood_ms_pg_per_L", "blood_im_pg_per_L"]].to_string())
    print("== Doku CV AUC eğrisi: " + ", ".join(f"{n}:{a:.3f}" for n, a in curve.items()))
    print(f"== PROTEIN ADAY PANELİ ({size}): {panel}")

    lock = {
        "panel": "Fısıltı protein aday paneli v2 (plazma ELISA; kanda doğrulanmamış)",
        "locked_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "code_commit_at_lock": git_commit(),
        "members": [{"gene": g, "uniprot": c.at[g, "uniprot"], "tier": c.at[g, "tier"], "direction": "up",
                     "hpa_secretome": c.at[g, "hpa_secretome"],
                     "blood_ms_pg_per_L": None if pd.isna(c.at[g, "blood_ms_pg_per_L"]) else c.at[g, "blood_ms_pg_per_L"],
                     "blood_im_pg_per_L": None if pd.isna(c.at[g, "blood_im_pg_per_L"]) else c.at[g, "blood_im_pg_per_L"]}
                    for g in panel],
        "blood_validation": "yok — elimizde plazma proteomik kohortu bulunmuyor (K-017)",
        "tissue_cv_auc_curve": {str(k): round(v, 4) for k, v in curve.items()},
    }
    (PROJECT_ROOT / "results" / "panel" / "protein_panel_v2.json").write_text(
        json.dumps(lock, ensure_ascii=False, indent=2), encoding="utf-8")
    print("[ok] results/panel/protein_panel_v2.json")


if __name__ == "__main__":
    main()
