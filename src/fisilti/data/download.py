"""Güvenli dosya indirme: yeniden denemeli, atomik yazımlı, SHA-256 kayıtlı.

Her indirilen dosya `data/raw/MANIFEST.tsv` dosyasına kaynak URL, boyut, SHA-256
ve indirme zamanıyla kaydedilir. Böylece hangi verinin ne zaman nereden
alındığı her an kanıtlanabilir (tekrarlanabilirlik).
"""

from __future__ import annotations

import csv
import hashlib
import threading
import time
from datetime import datetime, timezone
from pathlib import Path

import requests
from tqdm import tqdm

from fisilti.config import path

MANIFEST_FIELDS = ["file", "url", "bytes", "sha256", "downloaded_utc"]


def sha256sum(file: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(file, "rb") as f:
        while block := f.read(chunk):
            h.update(block)
    return h.hexdigest()


_MANIFEST_LOCK = threading.Lock()


def _record(file: Path, url: str) -> None:
    with _MANIFEST_LOCK:
        _record_unlocked(file, url)


def _record_unlocked(file: Path, url: str) -> None:
    raw = path("raw")
    manifest = raw / "MANIFEST.tsv"
    rel = file.relative_to(raw).as_posix()
    rows = []
    if manifest.exists():
        with open(manifest, encoding="utf-8", newline="") as f:
            rows = [r for r in csv.DictReader(f, delimiter="\t") if r["file"] != rel]
    rows.append(
        {
            "file": rel,
            "url": url,
            "bytes": file.stat().st_size,
            "sha256": sha256sum(file),
            "downloaded_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        }
    )
    rows.sort(key=lambda r: r["file"])
    with open(manifest, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=MANIFEST_FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def download(url: str, dest: Path, retries: int = 4, overwrite: bool = False) -> Path:
    """`url` adresini `dest` yoluna indirir. Dosya varsa atlar (overwrite=False)."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0 and not overwrite:
        return dest

    tmp = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(1, retries + 1):
        try:
            with requests.get(url, stream=True, timeout=60) as r:
                r.raise_for_status()
                # Sunucu gzip ile gönderirse Content-Length sıkıştırılmış boyuttur;
                # requests içeriği açtığı için o durumda boyut karşılaştırılmaz.
                encoded = "Content-Encoding" in r.headers
                total = None if encoded else (int(r.headers.get("Content-Length", 0)) or None)
                with open(tmp, "wb") as f, tqdm(
                    total=total, unit="B", unit_scale=True, desc=dest.name, leave=False
                ) as bar:
                    for block in r.iter_content(1 << 16):
                        f.write(block)
                        bar.update(len(block))
            if total is not None and tmp.stat().st_size != total:
                raise IOError(f"Eksik indirme: {tmp.stat().st_size}/{total} bayt")
            tmp.replace(dest)
            _record(dest, url)
            return dest
        except (requests.RequestException, IOError) as e:
            client_error = (
                isinstance(e, requests.HTTPError)
                and e.response is not None
                and 400 <= e.response.status_code < 500
            )
            if attempt == retries or client_error:
                tmp.unlink(missing_ok=True)
                raise RuntimeError(f"{url} indirilemedi: {e}") from e
            time.sleep(2**attempt)
    return dest
