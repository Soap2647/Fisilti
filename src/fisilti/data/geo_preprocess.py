"""GEO kohortlarının ön işlemesi: prob -> gen/miRNA eşlemesi, log dönüşümü,
teknik tekrarların ortalanması ve standart etiket tablosu.

Her kohort için çıktı:
  expr   : özellik x örnek (log2 ölçek)
  labels : örnek başına standart sütunlar
           group      -> 'PDAC' | 'normal' | 'healthy' | 'benign' | 'other_cancer'
           patient    -> hasta kimliği (gruplu CV için)
           stage      -> AJCC evresi (varsa)
           early      -> True/False/NaN  (Evre I-II)
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd

from fisilti.config import load_config, path
from fisilti.data.geo import parse_characteristics, read_series_matrix

EARLY = set(load_config()["analysis"]["early_stages"])


def _series(gse: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    f = next((path("raw") / "geo" / gse).glob("*_series_matrix.txt.gz"))
    expr, samples = read_series_matrix(f)
    return expr, samples, parse_characteristics(samples)


def _platform_table(file: Path) -> pd.DataFrame:
    with gzip.open(file, "rt", encoding="utf-8", errors="replace") as f:
        lines = []
        inside = False
        for line in f:
            if line.startswith("!platform_table_begin"):
                inside = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if inside:
                lines.append(line.rstrip("\n").split("\t"))
    return pd.DataFrame(lines[1:], columns=lines[0])


def _stage_norm(s: str | float) -> str | float:
    """'pStage IIB' / 'cStage IV' / 'Stage IIA' -> 'IIB' / 'IV' / 'IIA'."""
    if not isinstance(s, str) or ">" in s:  # '>IIB' gibi belirsiz evreler
        return np.nan
    m = re.search(r"\b(I{1,3}V?|IV)([ABC]?)\b", s.replace("Stage", " "))
    return m.group(0) if m else np.nan


def _early(stage) -> float:
    return np.nan if not isinstance(stage, str) else float(stage in EARLY)


# ---------------------------------------------------------------------------
# mRNA mikroarray -> gen sembolü
# ---------------------------------------------------------------------------

def probes_to_genes(expr: pd.DataFrame, gpl: str) -> pd.DataFrame:
    """Probları gen sembolüne indirger.

    Birden fazla gene eşlenen (/// içeren) veya eşlenmeyen problar atılır;
    aynı gene ait birden fazla probdan örnekler arası ortalama ifadesi en yüksek
    olan tutulur (yaygın 'maxMean' yaklaşımı).
    """
    ann = _platform_table(path("raw") / "geo" / "platforms" / f"{gpl}.annot.gz")
    sym = ann.set_index("ID")["Gene symbol"]
    sym = sym[(sym != "") & ~sym.str.contains("///", regex=False)]
    common = expr.index.intersection(sym.index)
    e = expr.loc[common]
    best = e.mean(axis=1).groupby(sym.loc[common]).idxmax()
    out = e.loc[best.values]
    out.index = best.index
    out.index.name = "gene"
    return out


def _avg_replicates(expr: pd.DataFrame, labels: pd.DataFrame, key: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aynı anahtara (hasta + doku tipi) sahip teknik tekrarları ortalar."""
    groups = labels.groupby(key, sort=False).groups
    cols, rows = {}, []
    for _, idx in groups.items():
        first = idx[0]
        cols[first] = expr[list(idx)].mean(axis=1)
        row = labels.loc[first].copy()
        row["n_replicates"] = len(idx)
        rows.append(row)
    return pd.DataFrame(cols), pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Doku kohortları
# ---------------------------------------------------------------------------

def gse15471() -> tuple[pd.DataFrame, pd.DataFrame]:
    expr, samples, ch = _series("GSE15471")
    labels = pd.DataFrame(index=expr.columns)
    labels["group"] = ch["sample"].map({"tumor": "PDAC", "normal": "normal"})
    labels["patient"] = "GSE15471_" + ch["patient"].astype(str)
    labels["stage"] = np.nan
    expr, labels = _avg_replicates(expr, labels, ["patient", "group"])
    return probes_to_genes(expr, "GPL570"), labels


def gse62452() -> tuple[pd.DataFrame, pd.DataFrame]:
    expr, samples, ch = _series("GSE62452")
    labels = pd.DataFrame(index=expr.columns)
    labels["group"] = ch["tissue"].map({"Pancreatic tumor": "PDAC", "adjacent pancreatic non-tumor": "normal"})
    # Başlık 'E005-Tc8' / 'E005-Tn' biçiminde; hasta kodu 'E' + sayı (baştaki sıfırlar normalize).
    pid = samples["title"].str.extract(r"^E0*(\d+)", expand=False)
    labels["patient"] = "GSE62452_" + pid.fillna(pd.Series(samples.index, index=samples.index))
    labels["stage"] = ch["stage"].map(_stage_norm)
    labels.loc[labels.group != "PDAC", "stage"] = np.nan
    return probes_to_genes(expr, "GPL6244"), labels


def gse16515() -> tuple[pd.DataFrame, pd.DataFrame]:
    expr, samples, ch = _series("GSE16515")
    expr = np.log2(expr.clip(lower=1.0))  # ham MAS5 yoğunlukları -> log2
    labels = pd.DataFrame(index=expr.columns)
    labels["group"] = np.where(ch["tissue"].str.startswith("Tumor"), "PDAC", "normal")
    pid = samples["title"].str.extract(r"Sample (\d+)-", expand=False)
    # Bir örneğin başlığında hasta no. boş ('Sample -Tumor'): GSM kimliği kullanılır.
    labels["patient"] = "GSE16515_" + pid.fillna(pd.Series(samples.index, index=samples.index))
    labels["stage"] = np.nan
    labels["sex"] = ch["sex"]
    labels["age"] = pd.to_numeric(ch["age"], errors="coerce")
    return probes_to_genes(expr, "GPL570"), labels


# ---------------------------------------------------------------------------
# Serum miRNA kohortları -> indeks = MIMAT erişim numarası (miRBase sürümünden bağımsız)
# ---------------------------------------------------------------------------

def _mimat_names(gpl: str) -> pd.Series:
    """MIMAT -> olgun miRNA adı (yalnızca insan)."""
    t = _platform_table(path("raw") / "geo" / "platforms" / f"{gpl}_family.soft.gz")
    if gpl == "GPL19117":
        t = t[t["Species Scientific Name"] == "Homo sapiens"]
        return pd.Series(t["Transcript ID(Array Design)"].values, index=t["Accession"].values)
    return pd.Series(t["miRNA_ID_LIST"].values, index=t["miRNA_ID"].values)


def gse85589() -> tuple[pd.DataFrame, pd.DataFrame]:
    expr, samples, ch = _series("GSE85589")
    expr.index = expr.index.str.replace("_st", "", regex=False)
    names = _mimat_names("GPL19117")
    expr = expr.loc[expr.index.intersection(names.index)]
    expr.index.name = "mimat"
    src = samples["source_name_ch1"]
    group = np.select(
        [src.str.contains("pancreatic cancer"), src.str.contains("healthy"), src.str.contains("cholelithiasis")],
        ["PDAC", "healthy", "benign"],
        default="other_cancer",
    )
    labels = pd.DataFrame(index=expr.columns)
    labels["group"] = group
    labels["patient"] = "GSE85589_" + labels.index
    labels["stage"] = np.nan
    labels["ca19_9_pos"] = ch["ca19-9"].map({">37": 1.0, "<37": 0.0})
    labels["subtype"] = src.str.extract(r"\((.*)\)", expand=False)
    labels["sex"] = ch["gender"]
    labels["age"] = pd.to_numeric(ch["age"], errors="coerce")
    return expr, labels


def gse59856() -> tuple[pd.DataFrame, pd.DataFrame]:
    expr, samples, ch = _series("GSE59856")
    names = _mimat_names("GPL18941")
    expr = expr.loc[expr.index.intersection(names.index)]
    expr.index.name = "mimat"
    ds = ch["disease_state"]
    group = ds.map(
        {
            "pancreatic cancer": "PDAC",
            "healthy control": "healthy",
            "benign pancreatic or biliary tract diseases": "benign",
        }
    ).fillna("other_cancer")
    labels = pd.DataFrame(index=expr.columns)
    labels["group"] = group
    labels["patient"] = "GSE59856_" + labels.index
    labels["stage"] = ch["tumor_stage"].map(_stage_norm)
    labels.loc[labels.group != "PDAC", "stage"] = np.nan
    labels["operable"] = ch["operability"].map({"operable": 1.0, "inoperable": 0.0})
    labels["subtype"] = ds
    labels["sex"] = ch["gender"]
    labels["age"] = pd.to_numeric(ch["age"], errors="coerce")
    return expr, labels


BUILDERS = {
    "GSE15471": gse15471,
    "GSE62452": gse62452,
    "GSE16515": gse16515,
    "GSE85589": gse85589,
    "GSE59856": gse59856,
}


def mirna_name_map() -> pd.Series:
    """İki serum platformunun birleşik MIMAT -> ad eşlemesi."""
    return pd.concat([_mimat_names("GPL19117"), _mimat_names("GPL18941")]).groupby(level=0).first()


def finalize_labels(labels: pd.DataFrame) -> pd.DataFrame:
    labels = labels.copy()
    labels["early"] = labels["stage"].map(_early)
    labels.index.name = "sample"
    return labels
