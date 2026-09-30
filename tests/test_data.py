import gzip

import numpy as np
import pandas as pd
import pytest

from fisilti.data.cptac import sample_labels
from fisilti.data.gdc import _primary_diagnosis
from fisilti.data.geo import _stub, parse_characteristics, read_series_matrix
from fisilti.data.geo_preprocess import _avg_replicates, _early, _stage_norm


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("IIB", "IIB"),
        ("Stage IIA", "IIA"),
        ("pStage IIB", "IIB"),
        ("ypStage IB", "IB"),
        ("cStage IV", "IV"),
        ("IVA", "IVA"),
        ("Stage III", "III"),
        ("uncertain", np.nan),
        (">IIB", np.nan),
        (np.nan, np.nan),
    ],
)
def test_stage_norm(raw, expected):
    out = _stage_norm(raw)
    if isinstance(expected, str):
        assert out == expected
    else:
        assert pd.isna(out)


def test_early_flags():
    assert _early("IA") == 1.0
    assert _early("IIB") == 1.0
    assert _early("III") == 0.0
    assert _early("IV") == 0.0
    assert np.isnan(_early(np.nan))


@pytest.mark.parametrize(
    "acc, stub",
    [("GSE15471", "GSE15nnn"), ("GSE1234", "GSE1nnn"), ("GSE123", "GSEnnn"),
     ("GPL570", "GPLnnn"), ("GPL19117", "GPL19nnn")],
)
def test_geo_ftp_stub(acc, stub):
    assert _stub(acc) == stub


def test_read_series_matrix(tmp_path):
    content = "\n".join(
        [
            '!Series_platform_id\t"GPL1"',
            '!Sample_title\t"A"\t"B"',
            '!Sample_geo_accession\t"GSM1"\t"GSM2"',
            '!Sample_characteristics_ch1\t"tissue: tumor"\t"tissue: normal"',
            '!Sample_characteristics_ch1\t"stage: IIB"\t"stage: NA"',
            "!series_matrix_table_begin",
            '"ID_REF"\t"GSM1"\t"GSM2"',
            '"p1"\t1.5\t2.5',
            '"p2"\t3\tnull',
            "!series_matrix_table_end",
        ]
    )
    f = tmp_path / "x_series_matrix.txt.gz"
    with gzip.open(f, "wt") as g:
        g.write(content)
    expr, samples = read_series_matrix(f)
    assert samples.attrs["platform"] == "GPL1"
    assert list(expr.columns) == ["GSM1", "GSM2"]
    assert expr.loc["p1", "GSM2"] == 2.5
    assert np.isnan(expr.loc["p2", "GSM2"])
    ch = parse_characteristics(samples)
    assert ch.loc["GSM1", "tissue"] == "tumor"
    assert ch.loc["GSM2", "stage"] == "NA"


def test_avg_replicates_merges_only_same_patient_and_group():
    expr = pd.DataFrame({"a": [1.0, 2.0], "a_rep": [3.0, 4.0], "b": [5.0, 6.0]}, index=["g1", "g2"])
    labels = pd.DataFrame(
        {"patient": ["p1", "p1", "p1"], "group": ["PDAC", "PDAC", "normal"]}, index=expr.columns
    )
    e, l = _avg_replicates(expr, labels, ["patient", "group"])
    assert e.shape == (2, 2)
    assert e["a"].tolist() == [2.0, 3.0]
    assert l.loc["a", "n_replicates"] == 2
    assert l.loc["b", "n_replicates"] == 1


def test_primary_diagnosis_skips_metastasis_record():
    # Gerçek vaka TCGA-2J-AAB9: ilk kayıt metastaz, evre ikinci kayıtta.
    dx = [
        {"classification_of_tumor": "metastasis", "diagnosis_is_primary_disease": False},
        {"classification_of_tumor": "primary", "diagnosis_is_primary_disease": True,
         "ajcc_pathologic_stage": "Stage IIB"},
    ]
    assert _primary_diagnosis(dx)["ajcc_pathologic_stage"] == "Stage IIB"
    assert _primary_diagnosis([]) == {}


def test_mean_of_pieces():
    from fisilti.data.cptac import _mean_of_pieces

    assert _mean_of_pieces("15;5;20") == pytest.approx(40 / 3)
    assert _mean_of_pieces("40") == 40.0
    assert np.isnan(_mean_of_pieces(np.nan))
    assert np.isnan(_mean_of_pieces("NA"))


def test_cptac_labels_mask_tumor_only_fields_for_normals():
    clin = pd.DataFrame(
        {"stage": ["IIB"], "early": [1.0], "age": [60], "sex": ["Male"], "histology_diagnosis": ["PDAC"],
         "tumor_site": ["head"], "Neoplastic_cellularity": [40.0], "Acinar_fraction": [5.0],
         "Islet_fraction": [0.0], "Stromal_fraction": [50.0], "Non_neoplastic_duct": [0.0],
         "Inflammation_fraction": [5.0], "follow_up_days": [300], "vital_status": ["Living"]},
        index=pd.Index(["C3L-1"], name="patient"),
    )
    lab = sample_labels(["C3L-1_T", "C3L-1_N"], clin)
    assert lab.loc["C3L-1_T", "group"] == "PDAC"
    assert lab.loc["C3L-1_N", "group"] == "normal"
    assert (lab["patient"] == "C3L-1").all()
    assert lab.loc["C3L-1_T", "stage"] == "IIB"
    assert pd.isna(lab.loc["C3L-1_N", "stage"])
    assert lab.loc["C3L-1_N", "age"] == 60  # hasta düzeyi bilgi korunur
