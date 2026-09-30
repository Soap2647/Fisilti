"""Proje yapılandırmasını yükler ve proje köküne göre yolları çözer."""

from __future__ import annotations

import random
from functools import lru_cache
from pathlib import Path

import numpy as np
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"


@lru_cache(maxsize=1)
def load_config() -> dict:
    with open(CONFIG_PATH, encoding="utf-8") as f:
        return yaml.safe_load(f)


def path(key: str) -> Path:
    """config.yaml içindeki `paths.<key>` değerini mutlak yola çevirir ve klasörü oluşturur."""
    p = PROJECT_ROOT / load_config()["paths"][key]
    p.mkdir(parents=True, exist_ok=True)
    return p


def set_seed(seed: int | None = None) -> int:
    """Tüm rastgelelik kaynaklarını sabitler (tekrarlanabilirlik için)."""
    seed = load_config()["project"]["seed"] if seed is None else seed
    random.seed(seed)
    np.random.seed(seed)
    try:
        import torch

        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    except ImportError:
        pass
    return seed
