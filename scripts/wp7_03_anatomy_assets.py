"""WP7-03: BodyParts3D anatomi verisinden kapanış filmi için 3B varlıklar üretir.

Kaynak: BodyParts3D 3.0 (99% indirgenmiş OBJ), © The Database Center for Life Science,
        CC Attribution-Share Alike 2.1 Japan.  (data/raw/bodyparts3d/)
Çıktı : sunum/assets/anatomy.js  ->  window.ANATOMY = {skin, bones, meshes}

Koordinatlar: BodyParts3D (mm; X = hastanın solu +, Y = arka +, Z = yukarı) ->
sahne (1 birim = 100 mm; x = X, y = Z, z = -Y; hasta kameraya bakar). Pankreasın
merkezi sahnede (0, 0.25, 0) noktasına taşınır (dalış sahnesiyle uyumlu).
Pankreas modeli 2 kez alt bölünür ve Taubin yöntemiyle yumuşatılır (yakın plan için).

Çalıştırma:  .venv\\Scripts\\python.exe scripts/wp7_03_anatomy_assets.py
"""

from __future__ import annotations

import base64
import json
import zipfile
from pathlib import Path

import numpy as np

from fisilti.config import PROJECT_ROOT, path

SRC = path("raw") / "bodyparts3d"
OUT = PROJECT_ROOT / "sunum" / "assets"
rng = np.random.default_rng(42)

# (FMA kimliği, ad, renk, saydamlık, tür)
MESHES = [
    ("FMA7198nsn", "pankreas", "#F08A5D", 1.0, "pancreas"),
    ("FMA10419", "pankreas kanalı", "#FFD08A", 0.55, "duct"),
    ("FMA7206", "onikiparmak bağırsağı", "#E8A27A", 0.42, "organ"),
    ("FMA7197", "karaciğer", "#C0604A", 0.20, "organ"),
    ("FMA7148", "mide", "#D98C98", 0.16, "organ"),
    ("FMA7196", "dalak", "#B0466A", 0.32, "organ"),
    ("FMA7202", "safra kesesi", "#7FB06A", 0.32, "organ"),
    ("FMA7204", "sağ böbrek", "#B86A50", 0.30, "organ"),
    ("FMA7205", "sol böbrek", "#B86A50", 0.30, "organ"),
    ("FMA3784", "inen aort", "#E04A40", 0.55, "vessel"),
    ("FMA50737", "çölyak arter", "#E04A40", 0.6, "vessel"),
    ("FMA14749", "üst mezenterik arter", "#E04A40", 0.6, "vessel"),
    ("FMA14773", "dalak arteri", "#E04A40", 0.6, "vessel"),
    ("FMA10951", "alt ana toplardamar", "#4A78D0", 0.5, "vessel"),
    ("FMA14332", "üst mezenterik ven", "#4A78D0", 0.6, "vessel"),
    ("FMA14331", "dalak veni", "#4A78D0", 0.6, "vessel"),
    ("FMA13072", "L1", "#E6DCC8", 0.30, "bone"),
    ("FMA13073", "L2", "#E6DCC8", 0.30, "bone"),
    ("FMA13074", "L3", "#E6DCC8", 0.30, "bone"),
]
RIBS = ["FMA7857", "FMA7882", "FMA7909", "FMA7957", "FMA7987", "FMA8012", "FMA8039", "FMA8066", "FMA8093", "FMA8148",
        "FMA8175", "FMA8202", "FMA8229", "FMA8256", "FMA8283", "FMA8310", "FMA8364", "FMA8391", "FMA8445", "FMA8472",
        "FMA8531", "FMA8532", "FMA8533", "FMA8534"]


def load_obj(fma: str) -> tuple[np.ndarray, np.ndarray]:
    p = SRC / "obj" / f"{fma}.obj"
    if not p.exists():
        p.parent.mkdir(exist_ok=True)
        with zipfile.ZipFile(SRC / "BodyParts3D_3.0_obj_99.zip") as z:
            p.write_bytes(z.read(f"BodyParts3D_3.0_obj_99/{fma}.obj"))
    v, f = [], []
    for line in p.read_text(errors="replace").splitlines():
        if line.startswith("v "):
            v.append([float(x) for x in line.split()[1:4]])
        elif line.startswith("f "):
            idx = [int(x.split("/")[0]) - 1 for x in line.split()[1:]]
            for k in range(1, len(idx) - 1):  # çokgenleri üçgenle
                f.append([idx[0], idx[k], idx[k + 1]])
    return np.asarray(v, dtype=np.float64), np.asarray(f, dtype=np.int64)


ORIGIN = None  # pankreas merkezi (mm), ilk yüklemede belirlenir


def to_scene(v: np.ndarray) -> np.ndarray:
    x, y, z = v[:, 0], v[:, 1], v[:, 2]
    o = ORIGIN
    return np.stack([(x - o[0]) / 100, (z - o[2]) / 100 + 0.25, -(y - o[1]) / 100], axis=1)


def subdivide(v: np.ndarray, f: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Orta nokta alt bölmesi: her üçgen 4 üçgene."""
    edges = {}
    verts = list(v)

    def mid(a, b):
        k = (a, b) if a < b else (b, a)
        if k not in edges:
            edges[k] = len(verts)
            verts.append((v[a] + v[b]) / 2)
        return edges[k]

    nf = []
    for a, b, c in f:
        ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
        nf += [[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]]
    return np.asarray(verts), np.asarray(nf)


def taubin(v: np.ndarray, f: np.ndarray, iters: int = 12, lam: float = 0.5, mu: float = -0.53) -> np.ndarray:
    """Hacmi koruyan Laplace yumuşatması (Taubin 1995)."""
    n = len(v)
    nb = [set() for _ in range(n)]
    for a, b, c in f:
        nb[a].update((b, c)); nb[b].update((a, c)); nb[c].update((a, b))
    nb = [np.fromiter(s, dtype=np.int64) for s in nb]
    v = v.copy()
    for _ in range(iters):
        for w in (lam, mu):
            avg = np.stack([v[s].mean(axis=0) if len(s) else v[i] for i, s in enumerate(nb)])
            v = v + w * (avg - v)
    return v


def sample_surface(v: np.ndarray, f: np.ndarray, n: int) -> tuple[np.ndarray, np.ndarray]:
    """Alanla orantılı yüzey örneklemesi; nokta + yüz normali."""
    a, b, c = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
    cr = np.cross(b - a, c - a)
    area = np.linalg.norm(cr, axis=1)
    keep = area > 0
    a, b, c, cr, area = a[keep], b[keep], c[keep], cr[keep], area[keep]
    idx = rng.choice(len(area), n, p=area / area.sum())
    r1, r2 = rng.random(n), rng.random(n)
    s = np.sqrt(r1)
    p = (1 - s)[:, None] * a[idx] + (s * (1 - r2))[:, None] * b[idx] + (s * r2)[:, None] * c[idx]
    nrm = cr[idx] / area[idx][:, None]
    return p, nrm


def b64(arr: np.ndarray, dtype) -> str:
    return base64.b64encode(np.ascontiguousarray(arr, dtype=dtype).tobytes()).decode()


def main() -> None:
    global ORIGIN
    pv, _ = load_obj("FMA7198nsn")
    ORIGIN = pv.mean(axis=0)
    out = {"source": "BodyParts3D 3.0 (99%), © DBCLS, CC BY-SA 2.1 JP", "meshes": []}

    # Deri: 60 bin yüzey noktası (normal: x-ışını kenar parlaklığı için)
    sv, sf = load_obj("FMA7163")
    sp, sn = sample_surface(sv, sf, 60000)
    sp = to_scene(sp)
    sn = np.stack([sn[:, 0], sn[:, 2], -sn[:, 1]], axis=1)
    out["skin"] = {"count": len(sp), "p": b64(sp, np.float32), "n": b64(sn, np.float32)}

    # Kaburgalar: nokta bulutu
    bp, bn = [], []
    for r in RIBS:
        v, f = load_obj(r)
        p, n = sample_surface(v, f, 450)
        bp.append(to_scene(p)); bn.append(np.stack([n[:, 0], n[:, 2], -n[:, 1]], axis=1))
    bp, bn = np.concatenate(bp), np.concatenate(bn)
    out["bones"] = {"count": len(bp), "p": b64(bp, np.float32), "n": b64(bn, np.float32)}

    for fma, name, color, opacity, kind in MESHES:
        v, f = load_obj(fma)
        if kind == "pancreas":
            for _ in range(2):
                v, f = subdivide(v, f)
            v = taubin(v, f)
        v = to_scene(v)
        out["meshes"].append({"id": fma, "name": name, "color": color, "opacity": opacity, "kind": kind,
                              "nv": len(v), "v": b64(v, np.float32), "i": b64(f, np.uint32)})
        print(f"[ok] {name:<22} köşe {len(v):>6}  üçgen {len(f):>6}  merkez {np.round(v.mean(axis=0), 2)}")

    OUT.mkdir(parents=True, exist_ok=True)
    js = "window.ANATOMY = " + json.dumps(out, separators=(",", ":")) + ";\n"
    (OUT / "anatomy.js").write_text(js, encoding="utf-8")
    print(f"[ok] deri {len(sp)} nokta, kaburga {len(bp)} nokta · anatomy.js {len(js) / 1e6:.1f} MB")
    print(f"     beden y aralığı {sp[:, 1].min():.2f} … {sp[:, 1].max():.2f}; ön karın yüzeyi z ≈ {np.percentile(sp[np.abs(sp[:, 1] - 0.25) < 0.3, 2], 98):.2f}")


if __name__ == "__main__":
    main()
