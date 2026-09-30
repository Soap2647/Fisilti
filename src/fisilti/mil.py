"""Dikkat mekanizmalı çoklu-örnek öğrenimi (Attention-MIL) ile çoklu-omik füzyon (WP3).

Ilse, Tomczak & Welling (ICML 2018) kapılı dikkat havuzlaması uyarlanmıştır:
  - Torba (bag)     = bir örnek (hasta dokusu)
  - Parça (instance)= o örnekte ölçülmüş her omik görünümü (mRNA, protein, miRNA, metilasyon)
Her görünüm kendi kodlayıcısıyla ortak d-boyutlu uzaya taşınır; dikkat ağırlıkları
örnek başına hangi omiğin karara katkı verdiğini gösterir. Ölçülmemiş görünümler
maskelenir — eksik veri doldurulmaz (MIL torbaları değişken boyutludur).
Eğitimde "görünüm düşürme" (modality dropout) ile tek omiğe bağımlılık azaltılır.
"""

from __future__ import annotations

import numpy as np
import torch
from torch import nn


class GatedAttentionMIL(nn.Module):
    def __init__(self, view_dims: dict[str, int], d: int = 32, attn_dim: int = 16, dropout: float = 0.3):
        super().__init__()
        self.views = list(view_dims)
        self.encoders = nn.ModuleDict({
            v: nn.Sequential(nn.Linear(n, d), nn.LayerNorm(d), nn.ReLU(), nn.Dropout(dropout))
            for v, n in view_dims.items()
        })
        self.att_v = nn.Linear(d, attn_dim)
        self.att_u = nn.Linear(d, attn_dim)
        self.att_w = nn.Linear(attn_dim, 1)
        self.head = nn.Sequential(nn.Dropout(dropout), nn.Linear(d, 1))

    def forward(self, xs: dict[str, torch.Tensor], mask: torch.Tensor):
        """xs[v]: (B, n_v) ; mask: (B, V) — 1 = görünüm ölçülmüş."""
        h = torch.stack([self.encoders[v](xs[v]) for v in self.views], dim=1)          # (B, V, d)
        a = self.att_w(torch.tanh(self.att_v(h)) * torch.sigmoid(self.att_u(h))).squeeze(-1)  # (B, V)
        a = a.masked_fill(mask == 0, float("-inf"))
        a = torch.softmax(a, dim=1)
        z = (a.unsqueeze(-1) * h).sum(dim=1)                                               # (B, d)
        return self.head(z).squeeze(-1), a


def _to_tensors(xs: dict[str, np.ndarray], mask: np.ndarray, idx: np.ndarray):
    return ({v: torch.as_tensor(x[idx], dtype=torch.float32) for v, x in xs.items()},
            torch.as_tensor(mask[idx], dtype=torch.float32))


def train_mil(xs: dict[str, np.ndarray], mask: np.ndarray, y: np.ndarray, train_idx: np.ndarray,
              val_idx: np.ndarray, seed: int = 42, epochs: int = 300, patience: int = 30,
              lr: float = 1e-3, weight_decay: float = 1e-2, view_dropout: float = 0.3) -> GatedAttentionMIL:
    """Erken durdurmalı eğitim; en iyi doğrulama kaybındaki ağırlıklar döndürülür."""
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    model = GatedAttentionMIL({v: x.shape[1] for v, x in xs.items()})
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    ytr = y[train_idx]
    pos_weight = torch.tensor((ytr == 0).sum() / max((ytr == 1).sum(), 1), dtype=torch.float32)
    loss_fn = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
    xv, mv = _to_tensors(xs, mask, val_idx)
    yv = torch.as_tensor(y[val_idx], dtype=torch.float32)
    best, best_state, wait = np.inf, None, 0
    for _ in range(epochs):
        model.train()
        order = rng.permutation(train_idx)
        for b in np.array_split(order, max(1, len(order) // 32)):
            xb, mb = _to_tensors(xs, mask, b)
            # Görünüm düşürme: her örnekte ölçülmüş görünümlerden en az biri korunur
            drop = torch.as_tensor(rng.random(mb.shape) < view_dropout, dtype=torch.float32)
            m2 = mb * (1 - drop)
            empty = m2.sum(1) == 0
            m2[empty] = mb[empty]
            logit, _ = model(xb, m2)
            loss = loss_fn(logit, torch.as_tensor(y[b], dtype=torch.float32))
            opt.zero_grad()
            loss.backward()
            opt.step()
        model.eval()
        with torch.no_grad():
            vl = loss_fn(model(xv, mv)[0], yv).item()
        if vl < best - 1e-4:
            best, wait = vl, 0
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
        else:
            wait += 1
            if wait >= patience:
                break
    model.load_state_dict(best_state)
    model.eval()
    return model


@torch.no_grad()
def predict_mil(model: GatedAttentionMIL, xs: dict[str, np.ndarray], mask: np.ndarray,
                idx: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x, m = _to_tensors(xs, mask, idx)
    logit, att = model(x, m)
    return torch.sigmoid(logit).numpy(), att.numpy()
