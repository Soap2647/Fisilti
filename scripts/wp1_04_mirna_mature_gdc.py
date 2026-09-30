"""WP1-04: CPTAC-PDAC ve TCGA-PAAD için olgun miRNA (MIMAT) matrislerini GDC'den oluşturur (K-010).

LinkedOmics CPTAC-PDAC normal miRNA tablosu tümör verisinin kopyası olduğundan
CPTAC miRNA'sı GDC'nin harmonize isoform dosyalarından yeniden türetilir. Aynı
işlem TCGA için de yapılır; böylece doku ve serum kohortları MIMAT anahtarıyla
ortaklaşır.

Çıktı: data/interim/<kohort>/mirna_mature.parquet  (MIMAT x örnek, RPM; filtrelenmemiş)
       data/interim/<kohort>/mirna_mature_files.tsv
Örnek adı: CPTAC -> '<hasta>_T|N', TCGA -> örnek barkodu.
Çalıştırma:  .venv\\Scripts\\python.exe scripts/wp1_04_mirna_mature_gdc.py
"""

from fisilti.config import path
from fisilti.data import gdc

JOBS = [
    # (kohort klasörü, GDC projesi, primary_site)
    ("cptac_pdac", "CPTAC-3", "Pancreas"),
    ("tcga_paad", "TCGA-PAAD", None),
]


def main() -> None:
    for cohort, project, site in JOBS:
        files = gdc.query_files(project, "mirna_mature", primary_site=site)
        files = files[files.sample_type.isin(["Primary Tumor", "Solid Tissue Normal"])]
        raw = path("raw") / cohort / "mirna_mature"
        interim = path("interim") / cohort
        interim.mkdir(parents=True, exist_ok=True)
        print(f"[..] {cohort}: {len(files)} dosya | {files.sample_type.value_counts().to_dict()}")
        gdc.download_files(files, raw)

        if cohort == "cptac_pdac":
            tag = files.sample_type.map({"Primary Tumor": "T", "Solid Tissue Normal": "N"})
            files = files.assign(sample=files["case"] + "_" + tag)
        files.to_csv(interim / "mirna_mature_files.tsv", sep="\t", index=False)

        m = gdc.build_matrix(files, raw, "mirna_mature").fillna(0.0)
        m.index.name = "mimat"
        m.to_parquet(interim / "mirna_mature.parquet")
        print(f"[ok] {cohort}: {m.shape[0]} olgun miRNA x {m.shape[1]} örnek")


if __name__ == "__main__":
    main()
