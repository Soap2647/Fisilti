"""MOFA+ çoklu-omik faktör analizi (WP2).

Tasarım (K-011):
- Tamamen denetimsiz: özellik ön seçimi yalnızca varyansa göre yapılır, sınıf
  etiketi kullanılmaz. Faktörlerin klinik anlamı model kurulduktan *sonra*
  değerlendirilir.
- WP3'te faktörler sınıflandırıcıya girdi olarak kullanılırsa `fit_mofa`
  çapraz doğrulamanın her eğitim katmanında yeniden çağrılır (sızıntı önleme).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

# Görünüm başına varyansa göre tutulacak en fazla özellik sayısı
DEFAULT_TOP = {"mrna": 3000, "protein": 2000, "methylation": 3000, "mirna": 500}


def beta_to_m(beta: pd.DataFrame, eps: float = 1e-3) -> pd.DataFrame:
    """Metilasyon beta -> M değeri; Gauss olabilirliğine daha uygun ölçek."""
    b = beta.clip(eps, 1 - eps)
    return np.log2(b / (1 - b))


def select_top_variance(m: pd.DataFrame, n: int, max_missing: float = 0.5) -> pd.DataFrame:
    """Eksik oranı `max_missing` altındaki özelliklerden en değişken `n` tanesi."""
    m = m[m.isna().mean(axis=1) < max_missing]
    return m.loc[m.var(axis=1).sort_values(ascending=False).index[:n]]


def prepare_views(views: dict[str, pd.DataFrame], top: dict[str, int] | None = None,
                  samples: list[str] | None = None) -> dict[str, pd.DataFrame]:
    """Görünümleri örnek x özellik biçimine getirir, ortak örnek eksenine hizalar.

    Bir örnekte ölçülmemiş görünüm tamamen NaN satır olur (MOFA+ bunu destekler).
    Özellik adları görünüm önekiyle benzersizleştirilir (ör. 'mrna:KRAS').
    """
    top = top or DEFAULT_TOP
    if samples is None:
        samples = sorted(set().union(*[set(v.columns) for v in views.values()]))
    out = {}
    for name, m in views.items():
        m = m.loc[:, m.columns.intersection(samples)]
        if name == "methylation":
            m = beta_to_m(m)
        m = select_top_variance(m, top.get(name, 2000))
        x = m.T.reindex(samples)
        x.columns = [f"{name}:{f}" for f in x.columns]
        out[name] = x
    return out


@dataclass
class MofaResult:
    factors: pd.DataFrame                      # örnek x faktör
    weights: dict[str, pd.DataFrame]           # görünüm -> özellik x faktör
    r2: pd.DataFrame                           # faktör x görünüm (açıklanan varyans, %)
    views: dict[str, pd.DataFrame] = field(repr=False, default_factory=dict)


def fit_mofa(views: dict[str, pd.DataFrame], n_factors: int = 15, seed: int = 42,
             outfile: Path | None = None, verbose: bool = False,
             convergence_mode: str = "medium", max_iter: int = 1000) -> MofaResult:
    from mofapy2.run.entry_point import entry_point

    names = list(views)
    samples = list(views[names[0]].index)
    ent = entry_point()
    ent.set_data_options(scale_views=True, center_groups=True)
    ent.set_data_matrix(
        data=[[views[v].to_numpy(dtype=float)] for v in names],
        likelihoods=["gaussian"] * len(names),
        views_names=names,
        groups_names=["all"],
        samples_names=[samples],
        features_names=[list(views[v].columns) for v in names],
    )
    ent.set_model_options(factors=n_factors, spikeslab_weights=True, ard_weights=True)
    ent.set_train_options(iter=max_iter, convergence_mode=convergence_mode, dropR2=0.005,
                          seed=seed, verbose=verbose, quiet=False)
    ent.build()
    ent.run()
    if outfile is not None:
        outfile.parent.mkdir(parents=True, exist_ok=True)
        ent.save(str(outfile), save_data=False)

    z = ent.model.nodes["Z"].getExpectation()
    k = z.shape[1]
    fnames = [f"F{i + 1}" for i in range(k)]
    factors = pd.DataFrame(z, index=samples, columns=fnames)
    w_all = ent.model.nodes["W"].getExpectation()
    weights = {v: pd.DataFrame(w_all[i], index=views[v].columns, columns=fnames) for i, v in enumerate(names)}
    r2 = variance_explained(views, factors, weights)
    return MofaResult(factors, weights, r2, views)


def load_mofa(hdf5: Path, views: dict[str, pd.DataFrame]) -> tuple[MofaResult, pd.DataFrame]:
    """Kaydedilmiş MOFA+ modelini yükler (yeniden eğitim gerekmez).

    Dönüş: (sonuç, MOFA+'ın kendi hesapladığı faktör x görünüm R² tablosu, %).
    İkincisi, `variance_explained` hesabımızı bağımsız olarak doğrulamak içindir.
    """
    import h5py

    def s(a):
        return [x.decode() if isinstance(x, bytes) else str(x) for x in a]

    with h5py.File(hdf5, "r") as f:
        names = s(f["views/views"][()])
        samples = s(f["samples/all"][()])
        z = f["expectations/Z/all"][()].T
        fnames = [f"F{i + 1}" for i in range(z.shape[1])]
        factors = pd.DataFrame(z, index=samples, columns=fnames)
        weights = {v: pd.DataFrame(f[f"expectations/W/{v}"][()].T, index=s(f[f"features/{v}"][()]),
                                   columns=fnames) for v in names}
        r2_mofa = pd.DataFrame(f["variance_explained/r2_per_factor/all"][()].T, index=fnames, columns=names)
    views = {v: views[v].reindex(index=samples, columns=weights[v].index) for v in names}
    return MofaResult(factors, weights, variance_explained(views, factors, weights), views), r2_mofa


def variance_explained(views: dict[str, pd.DataFrame], factors: pd.DataFrame,
                       weights: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Her faktörün her görünümde açıkladığı varyans (%, MOFA+ tanımı: 1 - SS_res/SS_tot).

    Veri, modelin gördüğü biçimde (özellik başına merkezlenmiş, görünüm ölçekli) hesaplanır;
    yalnızca gözlenen (NaN olmayan) değerler kullanılır.
    """
    rows = {}
    for v, x in views.items():
        y = x.to_numpy(dtype=float)
        y = y - np.nanmean(y, axis=0)
        y = y / np.nanstd(y)
        mask = ~np.isnan(y)
        ss_tot = np.sum(y[mask] ** 2)
        z = factors.to_numpy()
        w = weights[v].to_numpy()
        rows[v] = {
            f: 100 * (1 - np.sum((y - np.outer(z[:, j], w[:, j]))[mask] ** 2) / ss_tot)
            for j, f in enumerate(factors.columns)
        }
        rows[v]["toplam"] = 100 * (1 - np.sum((y - z @ w.T)[mask] ** 2) / ss_tot)
    return pd.DataFrame(rows)
