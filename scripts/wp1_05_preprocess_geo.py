"""WP1-05: GEO kohortlarını ön işler ve data/processed/geo/ altına yazar.

Çıktılar: <GSE>_expr.parquet, <GSE>_labels.tsv, mirna_names.tsv
Çalıştırma:  .venv\\Scripts\\python.exe scripts/wp1_05_preprocess_geo.py
"""

from fisilti.config import path
from fisilti.data import geo_preprocess as gp


def main() -> None:
    out = path("processed") / "geo"
    out.mkdir(parents=True, exist_ok=True)
    for gse, build in gp.BUILDERS.items():
        expr, labels = build()
        labels = gp.finalize_labels(labels)
        assert expr.columns.equals(labels.index), f"{gse}: ifade/etiket sırası uyuşmuyor"
        assert labels["group"].notna().all(), f"{gse}: etiketlenmemiş örnek var"
        expr.to_parquet(out / f"{gse}_expr.parquet")
        labels.to_csv(out / f"{gse}_labels.tsv", sep="\t")
        early = labels.loc[labels.group == "PDAC", "early"].value_counts().to_dict()
        print(
            f"[ok] {gse}: {expr.shape[0]:>5} özellik x {expr.shape[1]:>3} örnek | "
            f"gruplar {labels.group.value_counts().to_dict()} | "
            f"hasta {labels.patient.nunique()} | PDAC erken(1)/geç(0) {early}"
        )
    gp.mirna_name_map().rename("name").to_csv(out / "mirna_names.tsv", sep="\t", index_label="mimat")


if __name__ == "__main__":
    main()
