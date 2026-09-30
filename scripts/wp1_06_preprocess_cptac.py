"""WP1-06: CPTAC-PDAC omik matrislerini ve örnek etiketlerini data/processed/cptac_pdac/ altına yazar.

Önkoşul: wp1_01 (LinkedOmics) ve wp1_04 (GDC miRNA).
Tümör örnekleri çalışmanın 140 hastasıyla sınırlandırılır.
Çalıştırma:  .venv\\Scripts\\python.exe scripts/wp1_06_preprocess_cptac.py
"""

import pandas as pd

from fisilti.config import path
from fisilti.data import cptac


def main() -> None:
    out = path("processed") / "cptac_pdac"
    out.mkdir(parents=True, exist_ok=True)
    clin = cptac.clinical()
    all_samples: set[str] = set()
    for omic in cptac.OMICS:
        m = cptac.omic_matrix(omic, patients=clin.index)
        m.to_parquet(out / f"{omic}.parquet")
        all_samples |= set(m.columns)
        lab = cptac.sample_labels(list(m.columns), clin)
        print(
            f"[ok] {omic:<11} {m.shape[0]:>6} özellik x {m.shape[1]:>3} örnek | "
            f"{lab.group.value_counts().to_dict()} | eksik oranı {m.isna().mean().mean():.3f}"
        )
    labels = cptac.sample_labels(sorted(all_samples), clin)
    labels.to_csv(out / "labels.tsv", sep="\t")

    # Omik katmanları arası örnek kesişimi (MOFA+ eksik görünümlere izin verir, yine de raporlanır)
    pres = pd.DataFrame(
        {o: labels.index.isin(pd.read_parquet(out / f"{o}.parquet").columns) for o in cptac.OMICS},
        index=labels.index,
    )
    pres.to_csv(out / "omic_presence.tsv", sep="\t")
    for g in ["PDAC", "normal"]:
        p = pres[labels.group == g]
        print(f"[ok] {g}: toplam {len(p)} örnek; 4 omiğin tamamı: {p.all(axis=1).sum()}; "
              f"omik başına {p.sum().to_dict()}")
    pdac = labels[labels.group == "PDAC"]
    print(f"[ok] PDAC evre: {pdac.stage.value_counts(dropna=False).to_dict()}; erken={int(pdac.early.sum())}")


if __name__ == "__main__":
    main()
