"""Hasta düzeyinde, sızıntısız çapraz doğrulama altyapısı (WP3, K-003, K-004).

- `grouped_folds`: tekrarlı StratifiedGroupKFold; aynı hastanın örnekleri
  hiçbir zaman eğitim ve test katmanına birlikte düşmez.
- `UnivariateFilter`: yalnızca eğitim verisinde öğrenilen denetimli ön seçim
  (yön kısıtı + asiner program dışlaması, K-012).
- `auc_with_ci`: bootstrap güven aralıklı AUC.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold


def grouped_folds(y: np.ndarray, groups: np.ndarray, n_splits: int = 5, repeats: int = 1,
                  seed: int = 42) -> list[tuple[int, int, np.ndarray, np.ndarray]]:
    """(tekrar, katman, eğitim_idx, test_idx) listesi döndürür."""
    out = []
    for r in range(repeats):
        cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed + r)
        for k, (tr, te) in enumerate(cv.split(np.zeros(len(y)), y, groups)):
            assert not set(groups[tr]) & set(groups[te]), "hasta sızıntısı"
            out.append((r, k, tr, te))
    return out


def fast_auc(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Her sütun için AUC (Mann-Whitney U / n1*n0); NaN'lar sütun medyanıyla doldurulmuş olmalı."""
    pos = y == 1
    n1, n0 = pos.sum(), (~pos).sum()
    ranks = rankdata(x, axis=0)
    return (ranks[pos].sum(axis=0) - n1 * (n1 + 1) / 2) / (n1 * n0)


_STATS_CACHE: dict = {}


def _univariate_stats(X: pd.DataFrame, y: np.ndarray, max_missing: float,
                      exclude: frozenset) -> tuple[pd.Series, pd.Series]:
    """Özellik başına AUC ve medyan; aynı eğitim alt kümesi için önbelleklenir.

    Optuna denemeleri aynı iç eğitim katmanını tekrar tekrar kullandığından sonuç
    her seferinde aynıdır; önbellek yalnızca hesap süresini kısaltır.
    """
    key = (tuple(X.index), X.shape[1], y.tobytes(), max_missing, exclude)
    if key in _STATS_CACHE:
        return _STATS_CACHE[key]
    v = X.to_numpy(dtype=float)
    ok = (np.isnan(v).mean(axis=0) <= max_missing) & ~X.columns.isin(list(exclude))
    v = v[:, ok]
    cols = X.columns[ok]
    med = np.nanmedian(v, axis=0)
    v = np.where(np.isnan(v), med, v)
    res = (pd.Series(fast_auc(v, y), index=cols), pd.Series(med, index=cols))
    if len(_STATS_CACHE) > 256:
        _STATS_CACHE.clear()
    _STATS_CACHE[key] = res
    return res


class UnivariateFilter(BaseEstimator, TransformerMixin):
    """Eğitim katmanında tek değişkenli ön seçim.

    direction='up'   -> yalnızca tümörde artan özellikler (AUC > 0.5), kan testi yönelimi
    direction='both' -> |AUC - 0.5| büyüklüğüne göre
    exclude          -> dışlanacak özellik adları (ör. asiner program, K-012)
    k                -> tutulacak en fazla özellik sayısı
    Eksik değerler eğitim medyanıyla doldurulur (medyan yalnızca eğitimden öğrenilir).
    """

    def __init__(self, k: int = 200, direction: str = "up", exclude: frozenset = frozenset(),
                 max_missing: float = 0.3):
        self.k = k
        self.direction = direction
        self.exclude = exclude
        self.max_missing = max_missing

    def fit(self, X: pd.DataFrame, y):
        y = np.asarray(y)
        auc, medians = _univariate_stats(X, y, self.max_missing, self.exclude)
        self.medians_ = medians
        if self.direction == "up":
            score = auc[auc > 0.5] - 0.5
        else:
            score = (auc - 0.5).abs()
        self.auc_ = auc
        self.selected_ = list(score.sort_values(ascending=False).index[: self.k])
        self.medians_ = self.medians_[self.selected_]
        return self

    def transform(self, X: pd.DataFrame) -> pd.DataFrame:
        return X.reindex(columns=self.selected_).fillna(self.medians_)


def auc_with_ci(y: np.ndarray, score: np.ndarray, n_boot: int = 2000, seed: int = 0) -> tuple[float, float, float]:
    """AUC ve yüzdelik bootstrap %95 güven aralığı."""
    y, score = np.asarray(y), np.asarray(score)
    auc = roc_auc_score(y, score)
    rng = np.random.default_rng(seed)
    boots = []
    idx1, idx0 = np.where(y == 1)[0], np.where(y == 0)[0]
    for _ in range(n_boot):
        b = np.concatenate([rng.choice(idx1, len(idx1)), rng.choice(idx0, len(idx0))])  # tabakalı
        boots.append(roc_auc_score(y[b], score[b]))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return float(auc), float(lo), float(hi)


def acinar_program(X: pd.DataFrame, acinar_fraction: pd.Series, rho_min: float = 0.5) -> frozenset:
    """Tümör örneklerinde asiner fraksiyonla pozitif korelasyonlu (ρ > rho_min) özellikler.

    Sınıf etiketi değil patolojik doku bileşimi kullanıldığından tüm tümörlerde
    hesaplanması sızıntı yaratmaz (K-012).
    """
    ok = acinar_fraction.dropna().index.intersection(X.index)
    x = X.loc[ok]
    x = x.loc[:, x.notna().mean() >= 0.7]
    with np.errstate(divide="ignore", invalid="ignore"):  # sabit sütunlar -> NaN korelasyon
        r = x.rank().corrwith(acinar_fraction.loc[ok].rank())
    return frozenset(r.index[r > rho_min])
