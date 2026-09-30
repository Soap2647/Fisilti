"""Tablo verisi sınıflandırıcıları ve iç CV ile Optuna hiperparametre ayarı (WP3).

Her model bir sklearn Pipeline'dır: UnivariateFilter -> (ölçekleme) -> sınıflandırıcı.
Ön filtre de pipeline'ın parçası olduğundan iç CV'de her iç eğitim katmanında
yeniden öğrenilir (sızıntı yok).
"""

from __future__ import annotations

import warnings

import numpy as np
import optuna
import pandas as pd
from lightgbm import LGBMClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier

from fisilti.cv import UnivariateFilter, grouped_folds

optuna.logging.set_verbosity(optuna.logging.WARNING)
MODELS = ["logreg_l1", "xgboost", "lightgbm"]


def make_model(name: str, params: dict, filt: dict, pos_weight: float, seed: int = 42) -> Pipeline:
    f = UnivariateFilter(**filt)
    if name == "logreg_l1":
        clf = LogisticRegression(penalty="l1", solver="liblinear", C=params["C"],
                                 class_weight="balanced", max_iter=5000, random_state=seed)
        return Pipeline([("filter", f), ("scale", StandardScaler()), ("clf", clf)])
    if name == "xgboost":
        clf = XGBClassifier(
            n_estimators=params["n_estimators"], max_depth=params["max_depth"],
            learning_rate=params["learning_rate"], subsample=params["subsample"],
            colsample_bytree=params["colsample"], reg_lambda=params["reg_lambda"],
            min_child_weight=params["min_child_weight"], scale_pos_weight=pos_weight,
            tree_method="hist", n_jobs=1, random_state=seed, verbosity=0,
        )
        return Pipeline([("filter", f), ("clf", clf)])
    if name == "lightgbm":
        clf = LGBMClassifier(
            n_estimators=params["n_estimators"], num_leaves=params["num_leaves"],
            learning_rate=params["learning_rate"], min_child_samples=params["min_child_samples"],
            subsample=params["subsample"], subsample_freq=1, colsample_bytree=params["colsample"],
            reg_lambda=params["reg_lambda"], scale_pos_weight=pos_weight,
            n_jobs=1, random_state=seed, verbose=-1,
        )
        return Pipeline([("filter", f), ("clf", clf)])
    raise ValueError(name)


def suggest(trial: optuna.Trial, name: str) -> dict:
    if name == "logreg_l1":
        return {"C": trial.suggest_float("C", 1e-3, 10.0, log=True)}
    common = {
        "n_estimators": trial.suggest_int("n_estimators", 50, 400, step=50),
        "learning_rate": trial.suggest_float("learning_rate", 0.01, 0.2, log=True),
        "subsample": trial.suggest_float("subsample", 0.6, 1.0),
        "colsample": trial.suggest_float("colsample", 0.3, 1.0),
        "reg_lambda": trial.suggest_float("reg_lambda", 1e-3, 10.0, log=True),
    }
    if name == "xgboost":
        common["max_depth"] = trial.suggest_int("max_depth", 2, 5)
        common["min_child_weight"] = trial.suggest_float("min_child_weight", 0.5, 5.0)
    else:
        common["num_leaves"] = trial.suggest_int("num_leaves", 3, 15)
        common["min_child_samples"] = trial.suggest_int("min_child_samples", 3, 20)
    return common


def tune_and_fit(name: str, X: pd.DataFrame, y: np.ndarray, groups: np.ndarray, filt: dict,
                 n_trials: int = 15, inner_splits: int = 3, seed: int = 42) -> tuple[Pipeline, dict, float]:
    """İç grup-CV'de Optuna ile ayar; en iyi parametrelerle tüm eğitim verisinde son model."""
    pos_weight = float((y == 0).sum() / max((y == 1).sum(), 1))
    inner = grouped_folds(y, groups, n_splits=inner_splits, repeats=1, seed=seed)

    def objective(trial):
        params = suggest(trial, name)
        aucs = []
        for _, _, tr, te in inner:
            m = make_model(name, params, filt, pos_weight, seed)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                m.fit(X.iloc[tr], y[tr])
            aucs.append(roc_auc_score(y[te], m.predict_proba(X.iloc[te])[:, 1]))
        return float(np.mean(aucs))

    study = optuna.create_study(direction="maximize", sampler=optuna.samplers.TPESampler(seed=seed))
    study.optimize(objective, n_trials=n_trials)
    model = make_model(name, study.best_params, filt, pos_weight, seed)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model.fit(X, y)
    return model, study.best_params, study.best_value


def feature_importance(model: Pipeline) -> pd.Series:
    """Seçilen özelliklerin model içi önemi (lojistik: |katsayı|, ağaç: kazanç)."""
    feats = model.named_steps["filter"].selected_
    clf = model.named_steps["clf"]
    if isinstance(clf, LogisticRegression):
        imp = np.abs(clf.coef_[0])
    elif isinstance(clf, XGBClassifier):
        booster = clf.get_booster()
        gain = booster.get_score(importance_type="gain")
        # DataFrame ile eğitildiğinde anahtarlar özellik adıdır, aksi halde f0, f1, ...
        imp = np.array([gain.get(f, gain.get(f"f{i}", 0.0)) for i, f in enumerate(feats)])
    else:
        imp = clf.booster_.feature_importance(importance_type="gain").astype(float)
    s = pd.Series(imp, index=feats)
    return s / s.sum() if s.sum() > 0 else s
