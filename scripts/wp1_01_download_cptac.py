"""WP1-01: CPTAC-PDAC işlenmiş çoklu-omik tablolarını LinkedOmics'ten indirir.

Çalıştırma:  .venv\\Scripts\\python.exe scripts/wp1_01_download_cptac.py
"""

from fisilti.config import load_config, path
from fisilti.data.download import download


def main() -> None:
    src = load_config()["sources"]["cptac_pdac"]
    out_dir = path("raw") / "cptac_pdac"
    for key, fname in src["files"].items():
        dest = download(src["base_url"] + fname, out_dir / fname)
        print(f"[ok] {key:<18} {dest.stat().st_size / 1e6:8.2f} MB  {fname}")


if __name__ == "__main__":
    main()
