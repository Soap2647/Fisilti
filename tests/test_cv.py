import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from fisilti.cv import (UnivariateFilter, _STATS_CACHE, acinar_program, auc_with_ci, fast_auc,
                        grouped_folds)


def _toy(n=60, seed=0):
    rng = np.random.default_rng(seed)
    y = np.array([1] * (n // 2) + [0] * (n // 2))
    X = pd.DataFrame(rng.normal(size=(n, 6)), columns=list("abcdef"),
                     index=[f"s{i}" for i in range(n)])
    X["up"] = y * 2.0 + rng.normal(size=n) * 0.5        # tümörde artan
    X["down"] = -y * 2.0 + rng.normal(size=n) * 0.5     # tümörde azalan
    return X, y


def test_grouped_folds_never_split_a_patient():
    y = np.array([1, 0] * 20)
    groups = np.repeat(np.arange(20), 2)                 # her hastada T + N
    for _, _, tr, te in grouped_folds(y, groups, n_splits=5, repeats=3):
        assert not set(groups[tr]) & set(groups[te])


def test_fast_auc_matches_sklearn():
    X, y = _toy()
    auc = fast_auc(X.to_numpy(), y)
    for j, c in enumerate(X.columns):
        assert np.isclose(auc[j], roc_auc_score(y, X[c]))


def test_filter_direction_up_excludes_down_features():
    X, y = _toy()
    f = UnivariateFilter(k=3, direction="up").fit(X, y)
    assert f.selected_[0] == "up"
    assert "down" not in f.selected_
    f2 = UnivariateFilter(k=2, direction="both").fit(X, y)
    assert set(f2.selected_) == {"up", "down"}


def test_filter_exclude_and_train_only_medians():
    X, y = _toy()
    X.iloc[0, 0] = np.nan
    f = UnivariateFilter(k=8, direction="both", exclude=frozenset({"up"})).fit(X, y)
    assert "up" not in f.selected_
    test = X.iloc[:2].copy()
    test.iloc[:, :] = np.nan
    out = f.transform(test)
    # Eksik test değerleri EĞİTİM medyanıyla doldurulur
    assert np.allclose(out.iloc[0].to_numpy(), f.medians_.to_numpy())


def test_filter_cache_gives_identical_results():
    X, y = _toy(seed=3)
    _STATS_CACHE.clear()
    a = UnivariateFilter(k=4, direction="up").fit(X, y)
    b = UnivariateFilter(k=4, direction="up").fit(X, y)   # önbellekten
    assert a.selected_ == b.selected_
    assert a.auc_.equals(b.auc_)
    # Farklı alt küme önbellekten yanlışlıkla sonuç almamalı
    c = UnivariateFilter(k=4, direction="up").fit(X.iloc[:40], y[:40])
    assert not c.auc_.equals(a.auc_)


def test_auc_ci_contains_point_estimate():
    X, y = _toy()
    auc, lo, hi = auc_with_ci(y, X["up"].to_numpy(), n_boot=300)
    assert lo <= auc <= hi


def test_acinar_program_positive_correlation_only():
    rng = np.random.default_rng(0)
    frac = pd.Series(rng.uniform(0, 50, 40), index=[f"t{i}" for i in range(40)])
    X = pd.DataFrame({"acinar_gene": frac + rng.normal(size=40),
                      "tumor_gene": -frac + rng.normal(size=40),
                      "noise": rng.normal(size=40)}, index=frac.index)
    prog = acinar_program(X, frac)
    assert prog == frozenset({"acinar_gene"})           # negatif korelasyonlu tümör geni korunur
