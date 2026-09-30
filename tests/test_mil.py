import numpy as np
import torch

from fisilti.mil import GatedAttentionMIL, predict_mil, train_mil


def test_masked_view_gets_zero_attention_and_no_influence():
    torch.manual_seed(0)
    m = GatedAttentionMIL({"a": 4, "b": 3}).eval()
    xa = torch.randn(2, 4)
    mask = torch.tensor([[1.0, 0.0], [1.0, 1.0]])
    logit1, att = m({"a": xa, "b": torch.randn(2, 3)}, mask)
    assert att[0, 1].item() == 0.0                       # ölçülmemiş görünüm -> dikkat 0
    assert torch.allclose(att.sum(1), torch.ones(2))
    # Maskeli görünümün değeri değişse de çıktı değişmemeli
    logit2, _ = m({"a": xa, "b": torch.randn(2, 3) * 100}, mask)
    assert torch.isclose(logit1[0], logit2[0])


def test_mil_learns_separable_signal():
    rng = np.random.default_rng(0)
    n = 120
    y = np.array([1] * 60 + [0] * 60)
    xs = {"a": (rng.normal(size=(n, 5)) + y[:, None] * 1.5).astype(np.float32),
          "b": rng.normal(size=(n, 3)).astype(np.float32)}      # b yalnızca gürültü
    mask = np.ones((n, 2), dtype=np.float32)
    mask[::3, 0] = 0                                             # bazı örneklerde a eksik
    mask[::3, 1] = 1
    idx = rng.permutation(n)
    tr, va, te = idx[:70], idx[70:90], idx[90:]
    model = train_mil(xs, mask, y, tr, va, seed=0, epochs=150, patience=30)
    s, att = predict_mil(model, xs, mask, te)
    has_a = mask[te, 0] == 1
    from sklearn.metrics import roc_auc_score
    assert roc_auc_score(y[te][has_a], s[has_a]) > 0.9
    assert att[has_a, 0].mean() > att[has_a, 1].mean()          # bilgili görünüme daha çok dikkat
