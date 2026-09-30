"""CPTAC-PDAC çoklu-omik verisinin standart biçime getirilmesi.

Kaynaklar: mRNA, protein, metilasyon -> LinkedOmics; miRNA -> GDC isoform dosyaları (K-010).
Örnek kimliği: '<hasta>_T' (tümör) veya '<hasta>_N' (tümör-komşu normal).
Böylece aynı hastanın örnekleri `patient` sütunu üzerinden gruplanabilir (K-004).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from fisilti.config import load_config, path
from fisilti.data.geo_preprocess import _stage_norm

OMICS = ["mrna", "protein", "mirna", "methylation"]
EARLY = set(load_config()["analysis"]["early_stages"])


def _read(key: str) -> pd.DataFrame:
    fname = load_config()["sources"]["cptac_pdac"]["files"][key]
    return pd.read_csv(path("raw") / "cptac_pdac" / fname, sep="\t", index_col=0)


def mirna_log_filter(rpm: pd.DataFrame, min_log: float = 1.0, min_frac: float = 0.1) -> pd.DataFrame:
    """RPM -> log2(RPM+1); örneklerin en az %10'unda log2 ≥ 1 olmayan miRNA'lar atılır.

    Filtre etiket kullanmaz (denetimsiz), bu nedenle modelleme öncesi uygulanması
    sızıntı yaratmaz.
    """
    m = np.log2(rpm + 1)
    return m[(m >= min_log).mean(axis=1) >= min_frac]


def omic_matrix(omic: str, patients: pd.Index | None = None) -> pd.DataFrame:
    """Tümör ve normal örnekleri tek matriste birleştirir (özellik x örnek).

    `patients` verilirse yalnızca bu hastaların örnekleri tutulur.
    """
    if omic == "mirna":
        m = pd.read_parquet(path("interim") / "cptac_pdac" / "mirna_mature.parquet")
    else:
        t = _read(f"{omic}_tumor")
        n = _read(f"{omic}_normal")
        t.columns = [f"{c}_T" for c in t.columns]
        n.columns = [f"{c}_N" for c in n.columns]
        m = pd.concat([t, n], axis=1)
    if patients is not None:
        # Tümörler çalışma hastalarıyla sınırlanır; tümörü çalışılmamış bireylerden
        # alınan bağımsız normal örnekler (protein/metilasyonda 8 kişi) korunur.
        m = m.loc[:, m.columns.str.rsplit("_", n=1).str[0].isin(patients) | m.columns.str.endswith("_N")]
    if omic == "mirna":
        m = mirna_log_filter(m)
        m.index.name = "mimat"
    else:
        m.index.name = "gene"
    return m


def clinical() -> pd.DataFrame:
    c = _read("clinical")
    c.index.name = "patient"
    c["stage"] = c["tumor_stage_pathological"].map(_stage_norm)
    c["early"] = c["stage"].map(lambda s: np.nan if not isinstance(s, str) else float(s in EARLY))
    return c


def sample_labels(samples: list[str], clin: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "stage", "early", "age", "sex", "histology_diagnosis", "tumor_site",
        "Neoplastic_cellularity", "Acinar_fraction", "Islet_fraction", "Stromal_fraction",
        "Non_neoplastic_duct", "Inflammation_fraction", "follow_up_days", "vital_status",
    ]
    lab = pd.DataFrame(index=pd.Index(samples, name="sample"))
    lab["patient"] = [s.rsplit("_", 1)[0] for s in samples]
    lab["group"] = ["PDAC" if s.endswith("_T") else "normal" for s in samples]
    lab = lab.join(clin[keep], on="patient")
    # Evre ve doku fraksiyonları tümöre aittir; normal örneklerde anlamsızdır.
    tumor_only = ["stage", "early", "Neoplastic_cellularity", "Acinar_fraction", "Islet_fraction",
                  "Stromal_fraction", "Non_neoplastic_duct", "Inflammation_fraction"]
    lab.loc[lab.group == "normal", tumor_only] = np.nan
    return lab
