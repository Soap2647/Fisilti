import numpy as np
import pandas as pd

from fisilti.mofa import beta_to_m, prepare_views, select_top_variance, variance_explained


def test_beta_to_m_symmetry_and_clipping():
    b = pd.DataFrame({"s": [0.5, 0.0, 1.0, 0.8]})
    m = beta_to_m(b)["s"]
    assert m.iloc[0] == 0.0
    assert np.isfinite(m).all()                   # 0 ve 1 kırpılır
    assert np.isclose(m.iloc[3], np.log2(4))      # 0.8/0.2 = 4


def test_select_top_variance_respects_missingness():
    # satır = özellik, sütun = örnek
    m = pd.DataFrame(
        [[0, 1, 0, 1], [0, np.nan, np.nan, np.nan], [10, 20, 30, 40]],
        index=["low", "mostly_nan", "high"],
        columns=["a", "b", "c", "d"],
    )
    out = select_top_variance(m, n=5, max_missing=0.5)
    assert "mostly_nan" not in out.index          # %75 eksik -> atılır
    assert list(out.index) == ["high", "low"]      # varyansa göre sıralı


def test_prepare_views_aligns_samples_with_nan_rows():
    v1 = pd.DataFrame(np.random.default_rng(0).normal(size=(4, 3)), index=list("wxyz"), columns=["s1", "s2", "s3"])
    v2 = pd.DataFrame(np.random.default_rng(1).normal(size=(4, 2)), index=list("wxyz"), columns=["s2", "s3"])
    out = prepare_views({"mrna": v1, "protein": v2}, top={"mrna": 4, "protein": 4})
    assert list(out["mrna"].index) == ["s1", "s2", "s3"]
    assert list(out["protein"].index) == ["s1", "s2", "s3"]
    assert out["protein"].loc["s1"].isna().all()   # ölçülmemiş görünüm -> NaN satır
    assert all(c.startswith("mrna:") for c in out["mrna"].columns)


def test_variance_explained_recovers_known_structure():
    rng = np.random.default_rng(42)
    n, d = 300, 50
    z = rng.normal(size=(n, 2))
    w = rng.normal(size=(d, 2))
    w[:, 1] *= 0.3                                   # F2 daha zayıf
    y = z @ w.T + 0.1 * rng.normal(size=(n, d))
    y -= y.mean(axis=0)
    y /= y.std()
    # Model tahminleri ölçeklenmiş veriye göre: en küçük kareler ile W
    w_hat = np.linalg.lstsq(z, y, rcond=None)[0].T
    views = {"v": pd.DataFrame(y)}
    factors = pd.DataFrame(z, columns=["F1", "F2"])
    weights = {"v": pd.DataFrame(w_hat, columns=["F1", "F2"])}
    r2 = variance_explained(views, factors, weights)
    assert r2.loc["toplam", "v"] > 95
    assert r2.loc["F1", "v"] > r2.loc["F2", "v"] > 0
