"""WP1-03: TCGA-PAAD açık erişimli çoklu-omik veriyi GDC'den indirir ve matrislere dönüştürür.

Çıktılar (data/interim/tcga_paad/):
  <omik>.parquet            özellik x örnek matrisi
  <omik>_files.tsv          dosya -> örnek barkodu / örnek tipi eşlemesi
  clinical.tsv              vaka düzeyinde klinik veri
  gene_annotation.tsv       Ensembl ID -> gen adı / tipi

Çalıştırma:  .venv\\Scripts\\python.exe scripts/wp1_03_download_tcga.py [omik ...]
"""

import sys

from fisilti.config import load_config, path
from fisilti.data import gdc


def main(omics: list[str]) -> None:
    project = load_config()["sources"]["tcga_paad"]["gdc_project"]
    raw = path("raw") / "tcga_paad"
    out = path("interim") / "tcga_paad"
    out.mkdir(parents=True, exist_ok=True)

    clin = gdc.clinical(project)
    clin.to_csv(out / "clinical.tsv", sep="\t")
    print(f"[ok] klinik: {len(clin)} vaka")

    for omic in omics:
        files = gdc.query_files(project, omic)
        files.to_csv(out / f"{omic}_files.tsv", sep="\t", index=False)
        print(
            f"[..] {omic}: {len(files)} dosya, {files.file_size.sum() / 1e9:.2f} GB; "
            f"örnek tipleri: {files.sample_type.value_counts().to_dict()}"
        )
        gdc.download_files(files, raw / omic)
        mat = gdc.build_matrix(files, raw / omic, omic)
        mat.to_parquet(out / f"{omic}.parquet")
        print(f"[ok] {omic}: matris {mat.shape[0]} özellik x {mat.shape[1]} örnek")
        if omic == "mrna":
            gdc.gene_annotation(files, raw / omic).to_csv(out / "gene_annotation.tsv", sep="\t")


if __name__ == "__main__":
    main(sys.argv[1:] or list(gdc.OMICS))
