"""WP1-02: GEO doku ve serum kohortlarını ve platform anotasyonlarını indirir.

Çalıştırma:  .venv\\Scripts\\python.exe scripts/wp1_02_download_geo.py
"""

from fisilti.config import load_config
from fisilti.data.geo import download_platform, download_series, read_series_matrix


def main() -> None:
    sources = load_config()["sources"]
    geo = [k for k in sources if k.startswith("GSE")]
    platforms: set[str] = set()
    for gse in geo:
        for f in download_series(gse):
            expr, samples = read_series_matrix(f)
            platforms.add(samples.attrs["platform"])
            print(
                f"[ok] {gse:<9} {f.name:<40} platform={samples.attrs['platform']:<9} "
                f"prob={expr.shape[0]:>6}  örnek={expr.shape[1]:>4}  "
                f"NaN oranı={expr.isna().mean().mean():.3f}"
            )
    for gpl in sorted(platforms):
        p = download_platform(gpl)
        print(f"[ok] platform {gpl:<9} {p.name}")


if __name__ == "__main__":
    main()
