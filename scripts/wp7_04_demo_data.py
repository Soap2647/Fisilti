"""WP7-04: Canlı jüri demosu (sunum/demo/fisilti_demo.html) için veri dosyası üretir.

Kilitli 8 proteinlik panelin (results/panel/protein_panel_v2.json) üç bağımsız GEO doku kohortundaki
örnek bazlı değerlerini, panel skorunu ve rastgele 8 genlik panellerin AUC dağılımını dışa aktarır.
Skor, wp6_03 ile aynıdır: kohort içi gen bazlı z-skorların ortalaması. Karar eşiği etiket kullanılmadan
sabitlenir (skor > 0: kohort ortalamasının üstü). Hasta kimlikleri yayımlanmaz; örnekler sıra numarasıyla anılır.

Çıktı: sunum/demo/data.js  ->  window.FISILTI_DEMO = {...}
Çalıştırma:  .venv\\Scripts\\python.exe scripts/wp7_04_demo_data.py
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from fisilti.config import PROJECT_ROOT, path
from fisilti.cv import fast_auc

OUT = PROJECT_ROOT / "sunum" / "demo"
COHORTS = {"GSE15471": "Romanya · Badea 2008", "GSE62452": "ABD · Yang 2016", "GSE16515": "ABD (Mayo) · Pei 2009"}
N_RANDOM = 2000
rng = np.random.default_rng(42)


def main() -> None:
    lock = json.loads((PROJECT_ROOT / "results/panel/protein_panel_v2.json").read_text(encoding="utf-8"))
    genes = [m["gene"] for m in lock["members"]]
    out = {"genes": genes, "panel_locked": lock["locked_utc"], "cohorts": []}
    for gse, desc in COHORTS.items():
        e = pd.read_parquet(path("processed") / "geo" / f"{gse}_expr.parquet")
        l = pd.read_csv(path("processed") / "geo" / f"{gse}_labels.tsv", sep="\t", index_col=0).loc[e.columns]
        m = l.group.isin(["PDAC", "normal"]).to_numpy()
        e, l = e.loc[:, m], l[m]
        z = e.sub(e.mean(axis=1), axis=0).div(e.std(axis=1).replace(0, np.nan), axis=0)
        y = (l.group == "PDAC").to_numpy().astype(int)
        score = z.loc[genes].mean(axis=0).to_numpy()
        auc = float(fast_auc(score, y))
        # Rastgele 8 genlik paneller (aynı skor tanımı)
        zz = z.dropna().to_numpy()
        rand = np.array([fast_auc(zz[rng.choice(len(zz), len(genes), replace=False)].mean(axis=0), y) for _ in range(N_RANDOM)])
        pred = score > 0
        samples = [{"z": [round(float(v), 2) for v in z.loc[genes, s]], "s": round(float(sc), 3), "y": int(t),
                    "early": (None if pd.isna(l.at[s, "early"]) else int(l.at[s, "early"]))}
                   for s, sc, t in zip(l.index, score, y)]
        out["cohorts"].append({
            "id": gse, "desc": desc, "n_pdac": int(y.sum()), "n_normal": int((1 - y).sum()), "auc": round(auc, 3),
            "sens": round(float(pred[y == 1].mean()), 3), "spec": round(float((~pred)[y == 0].mean()), 3),
            "random_auc": [round(float(a), 3) for a in rand], "p_random": round(float((rand >= auc).mean()), 4),
            "samples": samples})
        print(f"[ok] {gse}: n={len(y)} AUC {auc:.3f}  duyarlılık {pred[y == 1].mean():.2f}  özgüllük {(~pred)[y == 0].mean():.2f}  "
              f"rastgele ≥ panel: {(rand >= auc).mean():.4f}")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "data.js").write_text("window.FISILTI_DEMO = " + json.dumps(out, ensure_ascii=False, separators=(",", ":")) + ";\n", encoding="utf-8")
    print(f"[ok] {OUT / 'data.js'}")


if __name__ == "__main__":
    main()
