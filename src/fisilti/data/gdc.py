"""NCI Genomic Data Commons (GDC) API üzerinden TCGA verisi indirme ve matris oluşturma.

Yalnızca açık erişimli (open) dosyalar kullanılır; her dosya GDC'nin verdiği
MD5 ile doğrulanır.
"""

from __future__ import annotations

import hashlib
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

from fisilti.config import path
from fisilti.data.download import download

API = "https://api.gdc.cancer.gov"

# data_type -> (GDC filtre alanları, çıktı adı)
OMICS = {
    "mrna": {"data_type": "Gene Expression Quantification", "experimental_strategy": "RNA-Seq"},
    "mirna": {"data_type": "miRNA Expression Quantification", "experimental_strategy": "miRNA-Seq"},
    # Olgun miRNA düzeyi (MIMAT); serum kohortlarıyla ortak anahtar (K-010)
    "mirna_mature": {"data_type": "Isoform Expression Quantification", "experimental_strategy": "miRNA-Seq"},
    "methylation": {"data_type": "Methylation Beta Value"},
    "protein": {"data_type": "Protein Expression Quantification"},
}


def _filters(project: str, extra: dict[str, str]) -> dict:
    content = [
        {"op": "in", "content": {"field": "cases.project.project_id", "value": [project]}},
        {"op": "in", "content": {"field": "access", "value": ["open"]}},
    ]
    for field, value in extra.items():
        content.append({"op": "in", "content": {"field": field, "value": [value]}})
    return {"op": "and", "content": content}


def query_files(project: str, omic: str, primary_site: str | None = None) -> pd.DataFrame:
    """Bir projedeki belirli omik türünün açık dosyalarını örnek bilgisiyle listeler.

    `primary_site` çok organlı projelerde (ör. CPTAC-3 -> 'Pancreas') filtre uygular.
    """
    extra = dict(OMICS[omic])
    if primary_site:
        extra["cases.primary_site"] = primary_site
    fields = [
        "file_id", "file_name", "md5sum", "file_size", "platform",
        "cases.submitter_id", "cases.samples.submitter_id", "cases.samples.sample_type",
    ]
    params = {
        "filters": json.dumps(_filters(project, extra)),
        "fields": ",".join(fields),
        "format": "JSON",
        "size": 5000,
    }
    hits = requests.get(f"{API}/files", params=params, timeout=120).json()["data"]["hits"]
    rows = []
    for h in hits:
        case = h["cases"][0]
        sample = case["samples"][0]
        rows.append(
            {
                "file_id": h["file_id"],
                "file_name": h["file_name"],
                "md5sum": h["md5sum"],
                "file_size": h["file_size"],
                "platform": h.get("platform"),
                "case": case["submitter_id"],
                "sample": sample["submitter_id"],
                "sample_type": sample["sample_type"],
            }
        )
    return pd.DataFrame(rows).sort_values("sample").reset_index(drop=True)


def _md5(file: Path) -> str:
    h = hashlib.md5()
    with open(file, "rb") as f:
        while block := f.read(1 << 20):
            h.update(block)
    return h.hexdigest()


def download_files(files: pd.DataFrame, out_dir: Path, workers: int = 8) -> None:
    """Dosyaları paralel indirir ve MD5 ile doğrular."""

    def one(row) -> str:
        dest = out_dir / row.file_id / row.file_name
        download(f"{API}/data/{row.file_id}", dest)
        if _md5(dest) != row.md5sum:
            dest.unlink()
            raise IOError(f"MD5 uyuşmazlığı: {row.file_name}")
        return row.file_name

    with ThreadPoolExecutor(workers) as ex:
        futures = [ex.submit(one, r) for r in files.itertuples()]
        for f in tqdm(as_completed(futures), total=len(futures), desc=out_dir.name):
            f.result()


# ---------------------------------------------------------------------------
# Dosya okuyucular: her biri tek örnek için pd.Series döndürür
# ---------------------------------------------------------------------------

def _read_mrna(f: Path) -> pd.Series:
    # STAR-Counts: ilk satır yorum, ardından N_unmapped vb. 4 özet satırı
    df = pd.read_csv(f, sep="\t", comment="#")
    df = df[df["gene_id"].str.startswith("ENSG")]
    return df.set_index("gene_id")["tpm_unstranded"]


def _read_mirna(f: Path) -> pd.Series:
    df = pd.read_csv(f, sep="\t")
    return df.set_index("miRNA_ID")["reads_per_million_miRNA_mapped"]


def _read_mirna_mature(f: Path) -> pd.Series:
    """isoforms.quantification: 'mature,MIMATxxxx' bölgeli okumaların RPM'i MIMAT başına toplanır."""
    df = pd.read_csv(f, sep="\t")
    df = df[df["miRNA_region"].str.startswith("mature,")]
    mimat = df["miRNA_region"].str.split(",", n=1).str[1]
    return df.groupby(mimat.values)["reads_per_million_miRNA_mapped"].sum()


def _read_methylation(f: Path) -> pd.Series:
    return pd.read_csv(f, sep="\t", header=None, index_col=0).iloc[:, 0]


def _read_protein(f: Path) -> pd.Series:
    df = pd.read_csv(f, sep="\t")
    return df.set_index("peptide_target")["protein_expression"]


READERS = {
    "mrna": _read_mrna,
    "mirna": _read_mirna,
    "mirna_mature": _read_mirna_mature,
    "methylation": _read_methylation,
    "protein": _read_protein,
}


def build_matrix(files: pd.DataFrame, out_dir: Path, omic: str) -> pd.DataFrame:
    """İndirilen dosyaları özellik x örnek matrisine birleştirir (sütun = örnek barkodu)."""
    reader = READERS[omic]
    dup = files["sample"].duplicated(keep=False)
    if dup.any():
        # Aynı örneğe ait tekrar ölçümler: dosya adına göre ilki tutulur (deterministik).
        print(f"[uyarı] {omic}: {files.loc[dup, 'sample'].nunique()} örnekte tekrar dosya; ilki tutuldu")
        files = files.sort_values(["sample", "file_name"]).drop_duplicates("sample")
    cols = {}
    for r in files.itertuples():
        cols[r.sample] = reader(out_dir / r.file_id / r.file_name)
    mat = pd.DataFrame(cols)
    mat.columns.name = "sample"
    return mat


def gene_annotation(files: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    """STAR-Counts dosyasından Ensembl ID -> gen adı / tipi eşlemesi."""
    r = files.iloc[0]
    df = pd.read_csv(out_dir / r.file_id / r.file_name, sep="\t", comment="#")
    df = df[df["gene_id"].str.startswith("ENSG")]
    return df.set_index("gene_id")[["gene_name", "gene_type"]]


def _primary_diagnosis(diagnoses: list[dict]) -> dict:
    """Bir vakanın birincil tümör tanısını seçer.

    GDC bazı vakalarda birden fazla tanı kaydı tutar (ör. sonradan gelişen
    metastaz). Liste sırası güvenilir değildir; evre bilgisi birincil kayıttadır.
    """
    for d in diagnoses:
        if d.get("diagnosis_is_primary_disease") is True or d.get("classification_of_tumor") == "primary":
            return d
    return diagnoses[0] if diagnoses else {}


def clinical(project: str) -> pd.DataFrame:
    """Vaka düzeyinde klinik veri (evre, yaş, cinsiyet, sağkalım)."""
    fields = [
        "submitter_id",
        "demographic.gender", "demographic.vital_status", "demographic.days_to_death",
        "diagnoses.age_at_diagnosis", "diagnoses.ajcc_pathologic_stage",
        "diagnoses.ajcc_pathologic_t", "diagnoses.ajcc_pathologic_n", "diagnoses.ajcc_pathologic_m",
        "diagnoses.primary_diagnosis", "diagnoses.tumor_grade",
        "diagnoses.days_to_last_follow_up",
        "diagnoses.diagnosis_is_primary_disease", "diagnoses.classification_of_tumor",
    ]
    params = {
        "filters": json.dumps(
            {"op": "in", "content": {"field": "project.project_id", "value": [project]}}
        ),
        "fields": ",".join(fields),
        "format": "JSON",
        "size": 5000,
    }
    hits = requests.get(f"{API}/cases", params=params, timeout=120).json()["data"]["hits"]
    rows = []
    for h in hits:
        d = _primary_diagnosis(h.get("diagnoses") or [])
        demo = h.get("demographic") or {}
        rows.append(
            {
                "case": h["submitter_id"],
                "gender": demo.get("gender"),
                "vital_status": demo.get("vital_status"),
                "days_to_death": demo.get("days_to_death"),
                "days_to_last_follow_up": d.get("days_to_last_follow_up"),
                "age_years": (d["age_at_diagnosis"] / 365.25) if d.get("age_at_diagnosis") else None,
                "stage": d.get("ajcc_pathologic_stage"),
                "t": d.get("ajcc_pathologic_t"),
                "n": d.get("ajcc_pathologic_n"),
                "m": d.get("ajcc_pathologic_m"),
                "grade": d.get("tumor_grade"),
                "primary_diagnosis": d.get("primary_diagnosis"),
            }
        )
    return pd.DataFrame(rows).set_index("case").sort_index()
