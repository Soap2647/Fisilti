"""GEO seri matrislerini indirme ve ayrıştırma."""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import pandas as pd
import requests

from fisilti.config import path
from fisilti.data.download import download

GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo"


def _stub(acc: str) -> str:
    """GSE15471 -> GSE15nnn (NCBI FTP klasör düzeni)."""
    return acc[:-3] + "nnn" if len(acc) > 6 else acc[:3] + "nnn"


def matrix_urls(gse: str) -> list[str]:
    """Bir serinin tüm series_matrix dosyalarının URL'lerini döndürür
    (çok platformlu serilerde birden fazla dosya olabilir)."""
    base = f"{GEO_FTP}/series/{_stub(gse)}/{gse}/matrix/"
    html = requests.get(base, timeout=60).text
    names = sorted(set(re.findall(rf'href="({gse}[^"]*_series_matrix\.txt\.gz)"', html)))
    if not names:
        raise RuntimeError(f"{gse} için series_matrix bulunamadı: {base}")
    return [base + n for n in names]


def download_series(gse: str) -> list[Path]:
    out = path("raw") / "geo" / gse
    return [download(u, out / u.rsplit("/", 1)[1]) for u in matrix_urls(gse)]


def download_platform(gpl: str) -> Path:
    """Platform anotasyon tablosunu (prob -> gen) indirir."""
    url = f"{GEO_FTP}/platforms/{_stub(gpl)}/{gpl}/annot/{gpl}.annot.gz"
    dest = path("raw") / "geo" / "platforms" / f"{gpl}.annot.gz"
    try:
        return download(url, dest)
    except RuntimeError:
        # .annot her platformda yok; tam SOFT dosyasına düş.
        url = f"{GEO_FTP}/platforms/{_stub(gpl)}/{gpl}/soft/{gpl}_family.soft.gz"
        return download(url, dest.with_name(f"{gpl}_family.soft.gz"))


def read_series_matrix(file: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """series_matrix.txt.gz dosyasını (ifade, örnek-metaverisi) olarak okur.

    ifade: satır = prob, sütun = GSM
    metaveri: satır = GSM, sütun = !Sample_* alanları (tekrarlı alanlar _1, _2 ...)
    """
    meta: dict[str, list[list[str]]] = {}
    platform = None
    rows: list[list[str]] = []
    in_table = False
    with gzip.open(file, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                rows.append([c.strip('"') for c in line.split("\t")])
            elif line.startswith("!Sample_"):
                key, *vals = line.split("\t")
                meta.setdefault(key[len("!Sample_"):], []).append([v.strip('"') for v in vals])
            elif line.startswith("!Series_platform_id"):
                platform = line.split("\t")[1].strip('"')

    header, body = rows[0], rows[1:]
    expr = pd.DataFrame(body, columns=header).set_index(header[0])
    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr.index.name = "probe"

    cols = {}
    for key, occurrences in meta.items():
        for i, vals in enumerate(occurrences):
            cols[key if len(occurrences) == 1 else f"{key}_{i + 1}"] = vals
    samples = pd.DataFrame(cols)
    samples.index = samples["geo_accession"]
    samples.index.name = "gsm"
    samples.attrs["platform"] = platform
    return expr, samples


def parse_characteristics(samples: pd.DataFrame) -> pd.DataFrame:
    """'characteristics_ch1_*' sütunlarındaki 'anahtar: değer' çiftlerini ayrı sütunlara açar."""
    out = pd.DataFrame(index=samples.index)
    for c in [c for c in samples.columns if c.startswith("characteristics_ch1")]:
        for gsm, v in samples[c].items():
            if ":" in v:
                k, val = v.split(":", 1)
                out.loc[gsm, k.strip().lower().replace(" ", "_")] = val.strip()
    return out
