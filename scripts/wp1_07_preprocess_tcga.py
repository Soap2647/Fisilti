"""WP1-07: TCGA-PAAD matrislerini standart biçime getirir (data/processed/tcga_paad/).

- Yalnızca duktal adenokarsinom histolojisi tutulur (nöroendokrin vb. dışlanır, K-008).
- mRNA: TPM -> log2(TPM+1), Ensembl -> gen sembolü (protein kodlayan + lncRNA; aynı sembolde en yüksek ortalama).
- miRNA: olgun düzey (MIMAT, wp1_04), RPM -> log2(RPM+1) + düşük ifade filtresi.
- Örnek adı: TCGA örnek barkodu (ör. TCGA-2J-AAB1-01A); patient = ilk 12 karakter.

Çalıştırma:  .venv\\Scripts\\python.exe scripts/wp1_07_preprocess_tcga.py
"""

import numpy as np
import pandas as pd

from fisilti.config import load_config, path
from fisilti.data.cptac import mirna_log_filter

EARLY = set(load_config()["analysis"]["early_stages"])
# TCGA primary_diagnosis değerlerinden duktal adenokarsinom kabul edilenler
DUCTAL = {
    "Infiltrating duct carcinoma, NOS",
    "Adenocarcinoma, NOS",
    "Mucinous adenocarcinoma",  # kolloid; PDAC alt tipi olarak TCGA'da tutulur
}


def main() -> None:
    src = path("interim") / "tcga_paad"
    out = path("processed") / "tcga_paad"
    out.mkdir(parents=True, exist_ok=True)

    clin = pd.read_csv(src / "clinical.tsv", sep="\t", index_col=0)
    print(f"[..] primary_diagnosis: {clin.primary_diagnosis.value_counts().to_dict()}")
    clin["stage"] = clin["stage"].str.replace("Stage ", "", regex=False)
    clin["early"] = clin["stage"].map(lambda s: np.nan if not isinstance(s, str) else float(s in EARLY))
    ductal = set(clin.index[clin.primary_diagnosis.isin(DUCTAL)])

    frames = {}
    # omik -> (ara dosya adı); miRNA olgun düzeyde (wp1_04) okunur
    sources = {"mrna": "mrna", "mirna": "mirna_mature", "methylation": "methylation", "protein": "protein"}
    for omic, stem in sources.items():
        f = src / f"{stem}.parquet"
        if not f.exists():
            print(f"[atla] {omic}: henüz indirilmedi")
            continue
        m = pd.read_parquet(f)
        files = pd.read_csv(src / f"{stem}_files.tsv", sep="\t").drop_duplicates("sample").set_index("sample")
        stype = files.loc[m.columns, "sample_type"]
        patient = m.columns.str[:12]
        keep = (stype.isin(["Primary Tumor", "Solid Tissue Normal"]).values) & patient.isin(ductal)
        m = m.loc[:, keep]

        if omic == "mrna":
            ann = pd.read_csv(src / "gene_annotation.tsv", sep="\t", index_col=0)
            ann = ann[ann.gene_type.isin(["protein_coding", "lncRNA"])]
            m = np.log2(m.loc[m.index.intersection(ann.index)] + 1)
            best = m.mean(axis=1).groupby(ann.loc[m.index, "gene_name"]).idxmax()
            m = m.loc[best.values]
            m.index = best.index
            m.index.name = "gene"
        elif omic == "mirna":
            m = mirna_log_filter(m)
            m.index.name = "mimat"
        m.to_parquet(out / f"{omic}.parquet")
        frames[omic] = m
        print(f"[ok] {omic:<11} {m.shape[0]:>6} özellik x {m.shape[1]:>3} örnek | "
              f"{stype[keep].value_counts().to_dict()}")

    samples = sorted(set().union(*[set(m.columns) for m in frames.values()]))
    lab = pd.DataFrame(index=pd.Index(samples, name="sample"))
    lab["patient"] = lab.index.str[:12]
    lab["group"] = np.where(lab.index.str[13:15] == "11", "normal", "PDAC")  # 01=tümör, 11=normal
    lab = lab.join(clin[["stage", "early", "age_years", "gender", "grade", "vital_status",
                         "days_to_death", "days_to_last_follow_up", "primary_diagnosis"]], on="patient")
    lab.loc[lab.group == "normal", ["stage", "early"]] = np.nan
    lab.to_csv(out / "labels.tsv", sep="\t")
    p = lab[lab.group == "PDAC"].drop_duplicates("patient")
    print(f"[ok] etiketler: {lab.group.value_counts().to_dict()} | duktal dışı dışlanan vaka: "
          f"{len(clin) - len(ductal)} | PDAC evre: {p.stage.value_counts(dropna=False).to_dict()}")


if __name__ == "__main__":
    main()
