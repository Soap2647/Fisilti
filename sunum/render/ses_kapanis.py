"""Fısıltı kapanış filmi için prosedürel ses tasarımı (telifsiz, tamamen sentez).

Ses olayları görüntüyle aynı zaman çizelgesinden gelir: render_kapanis.mjs,
sayfadaki window.__film.cues listesini kapanis_cues.json olarak yazar.
Çıktı: kapanis_ses.wav (48 kHz, stereo, 24-bit)

Çalıştırma:  ..\\..\\.venv\\Scripts\\python.exe ses_kapanis.py
"""

from __future__ import annotations

import json
import wave
from pathlib import Path

import numpy as np
from scipy.signal import butter, fftconvolve, sosfilt

SR = 48000
HERE = Path(__file__).parent
rng = np.random.default_rng(7)


# --------------------------------------------------------------------------- yardımcılar
def t_axis(dur):
    return np.arange(int(dur * SR)) / SR


def noise(dur):
    return rng.standard_normal(int(dur * SR))


def brown(dur):
    x = np.cumsum(noise(dur))
    x -= np.linspace(x[0], x[-1], len(x))
    return x / (np.abs(x).max() + 1e-9)


def lp(x, fc, order=2):
    return sosfilt(butter(order, fc, "low", fs=SR, output="sos"), x)


def hp(x, fc, order=2):
    return sosfilt(butter(order, fc, "high", fs=SR, output="sos"), x)


def svf_sweep(x, f0, f1, q=0.7, curve="exp"):
    """Zamanla değişen bant geçiren filtre (Chamberlin SVF)."""
    n = len(x)
    u = np.linspace(0, 1, n)
    fc = f0 * (f1 / f0) ** u if curve == "exp" else f0 + (f1 - f0) * u
    f = 2 * np.sin(np.pi * np.minimum(fc, SR / 6) / SR)
    low = band = 0.0
    out = np.empty(n)
    damp = 1 / q
    for i in range(n):
        high = x[i] - low - damp * band
        band += f[i] * high
        low += f[i] * band
        out[i] = band
    return out


def env_exp(dur, attack=0.005, decay=0.3):
    t = t_axis(dur)
    a = np.minimum(1, t / max(attack, 1e-4))
    return a * np.exp(-np.maximum(0, t - attack) / decay)


def pan2(x, p):
    """Eşit güçlü pan: p = -1 (sol) ... +1 (sağ); p dizi olabilir."""
    ang = (np.asarray(p) + 1) * np.pi / 4
    return np.stack([x * np.cos(ang), x * np.sin(ang)])


def gliss(f0, f1, dur, curve=2.0):
    u = np.linspace(0, 1, int(dur * SR)) ** curve
    f = f0 + (f1 - f0) * u
    return np.sin(2 * np.pi * np.cumsum(f) / SR)


# --------------------------------------------------------------------------- ses nesneleri
def heart(c):
    g = c.get("gain", 1.0)
    out = []
    for off, f0, amp in [(0.0, 58, 1.0), (0.3, 70, 0.72)]:
        d = 0.45
        s = gliss(f0, f0 * 0.62, d, 0.6) * env_exp(d, 0.004, 0.11)
        s += lp(noise(d), 120) * env_exp(d, 0.002, 0.05) * 0.8
        out.append((off, s * amp * g))
    L = int(0.8 * SR)
    y = np.zeros(L)
    for off, s in out:
        i = int(off * SR)
        y[i:i + len(s)] += s[: L - i]
    return pan2(y * 1.1, 0.0), 0.12


def gather(c):
    d = c["dur"]
    x = svf_sweep(noise(d), 180, 4200, q=1.4) * (np.linspace(0, 1, int(d * SR)) ** 2.2)
    x[-int(0.04 * SR):] *= np.linspace(1, 0, int(0.04 * SR))
    return pan2(x * 0.5, np.linspace(-0.7, 0.5, len(x))), 0.35


def drone(c):
    d = c["dur"]
    t = t_axis(d)
    rise = np.clip(t / (d - 3.5), 0, 1) ** 1.6
    fade = np.clip((d - t) / 3.0, 0, 1)
    y = np.zeros_like(t)
    for f, a in [(55, 1.0), (82.41, 0.6), (110, 0.45), (164.81, 0.25), (220, 0.12)]:
        y += a * np.sin(2 * np.pi * f * t + 0.3 * np.sin(2 * np.pi * 0.07 * t))
    y += 0.35 * lp(brown(d), 300)
    y = lp(y, 900) * (0.25 + 0.75 * rise) * fade * np.minimum(1, t / 1.5)
    st = np.stack([y, np.roll(y, int(0.012 * SR))])
    return st * 0.28, 0.2


def glow(c):
    d = 4.0
    t = t_axis(d)
    f = c["freq"]
    y = sum(a * np.sin(2 * np.pi * f * m * t) for m, a in [(1, 1), (2, 0.5), (3, 0.25), (4.02, 0.12)])
    y *= np.minimum(1, t / 0.25) * np.exp(-t / 1.4)
    return pan2(y * 0.3, 0.1), 0.45


def riser(c):
    d = c["dur"]
    n = int(d * SR)
    u = np.linspace(0, 1, n)
    x = svf_sweep(noise(d), 200, 7000, q=2.0) * 0.6 + gliss(90, 900, d, 2.2) * 0.35
    x *= u ** 2.6
    return pan2(x * 0.55, 0.0), 0.3


def whoosh(c):
    d = c["dur"]
    n = int(d * SR)
    u = np.linspace(0, 1, n)
    x = svf_sweep(noise(d), 2600, 140, q=1.1)
    e = np.minimum(1, u / 0.18) * (1 - u) ** 1.3
    y = x * e
    if c.get("big"):
        y += lp(noise(d), 90) * e * 1.2 + gliss(70, 35, d, 1.0) * e * 0.5
    return pan2(y * (0.9 if c.get("big") else 0.6), np.linspace(-0.8, 0.8, n)), 0.35


def boom(c):
    g = c.get("gain", 1.0)
    d = 4.0
    t = t_axis(d)
    y = gliss(48, 26, d, 0.5) * np.exp(-t / 1.1) * np.minimum(1, t / 0.004)
    y += lp(noise(d), 380) * np.exp(-t / 0.18) * 0.9
    y += hp(noise(d), 2500) * np.exp(-t / 0.05) * 0.25
    return pan2(y * g * 1.2, 0.0), 0.55


def tissue(c):
    d = c["dur"]
    t = t_axis(d)
    base = lp(brown(d), 260) * (0.6 + 0.4 * np.sin(2 * np.pi * 0.18 * t))
    wet = svf_sweep(noise(d), 380, 520, q=6.0, curve="lin") * (0.5 + 0.5 * np.sin(2 * np.pi * 0.31 * t + 1)) * 0.25
    y = (base + wet) * np.minimum(1, t / 1.5) * np.minimum(1, (d - t) / 2.0)
    st = np.stack([y, lp(brown(d), 260) * 0.6 * np.minimum(1, t / 1.5) * np.minimum(1, (d - t) / 2.0) + wet])
    return st * 0.22, 0.25


def bubble(c):
    d = 0.14
    f = c["freq"]
    y = gliss(f, f * 1.7, d, 0.7) * env_exp(d, 0.003, 0.04)
    return pan2(y * 0.16, rng.uniform(-0.7, 0.7)), 0.35


def crackle(c):
    d = 0.12
    y = np.zeros(int(d * SR))
    for _ in range(rng.integers(3, 6)):
        i = int(rng.uniform(0, 0.05) * SR)
        b = hp(noise(0.004), 1800) * np.hanning(int(0.004 * SR))
        y[i:i + len(b)] += b * rng.uniform(0.4, 1)
    y += np.sin(2 * np.pi * 190 * t_axis(d)) * env_exp(d, 0.002, 0.03) * 0.5
    return pan2(y * 0.36, c.get("pan", 0)), 0.25


def dissonance(c):
    d = c["dur"]
    t = t_axis(d)
    sw = np.clip(t / (d * 0.8), 0, 1) ** 1.5
    trem = 1 - 0.35 * sw * (0.5 + 0.5 * np.sin(2 * np.pi * (2 + 4 * sw) * t))
    y = sum(a * np.sin(2 * np.pi * f * t) for f, a in [(110, 1), (155.56, 0.8), (116.5, 0.5), (233.1, 0.3)])
    y = lp(y, 700) * sw * trem * np.minimum(1, (d - t) / 1.5)
    return np.stack([y, np.roll(y, int(0.02 * SR))]) * 0.3, 0.3


def ping(c):
    d = 1.6
    t = t_axis(d)
    f = c["freq"]
    y = (np.sin(2 * np.pi * f * t) + 0.35 * np.sin(2 * np.pi * f * 2.76 * t) * np.exp(-t / 0.2)) * env_exp(d, 0.002, 0.45)
    return pan2(y * 0.2, c.get("pan", 0)), 0.45


def flow(c):
    d = c["dur"]
    t = t_axis(d)
    pulse = 0.55 + 0.45 * np.maximum(0, np.sin(2 * np.pi * 1.2 * t)) ** 3
    y = lp(brown(d), 420) * pulse + svf_sweep(noise(d), 250, 330, q=3, curve="lin") * 0.15 * pulse
    y *= np.minimum(1, t / 1.0) * np.minimum(1, (d - t) / 0.8)
    st = np.stack([y, lp(brown(d), 420) * pulse * np.minimum(1, t / 1.0) * np.minimum(1, (d - t) / 0.8)])
    beats = np.zeros_like(t)
    for k in np.arange(0.2, d, 1 / 1.2):
        s, _ = heart({"gain": 0.35})
        i = int(k * SR)
        m = min(s.shape[1], len(t) - i)
        beats[i:i + m] += s[0, :m]
    return st * 0.35 + np.stack([beats, beats]) * 0.6, 0.2


def reverse(c):
    d = c["dur"]
    u = np.linspace(0, 1, int(d * SR))
    y = svf_sweep(noise(d), 400, 5000, q=1.2) * u ** 3
    return pan2(y * 0.5, 0.0), 0.4


def chime(c):
    d = 3.2
    t = t_axis(d)
    f = c["freq"]
    parts = [(1, 1.0, 1.4), (2.0, 0.45, 0.8), (3.01, 0.25, 0.5), (4.2, 0.12, 0.3), (5.43, 0.06, 0.2)]
    y = sum(a * np.sin(2 * np.pi * f * m * t) * np.exp(-t / dec) for m, a, dec in parts)
    y *= np.minimum(1, t / 0.003)
    y += np.sin(2 * np.pi * f / 2 * t) * np.exp(-t / 0.9) * 0.2
    return pan2(y * 0.2, c.get("pan", 0)), 0.5


def pad(c):
    d = c["dur"]
    t = t_axis(d)
    chord = [110, 164.81, 220, 261.63, 329.63, 493.88]
    y = np.zeros_like(t)
    for f in chord:
        for h, a in [(1, 1), (2, 0.35), (3, 0.15)]:
            y += a * np.sin(2 * np.pi * f * h * t * (1 + 0.0015 * np.sin(2 * np.pi * 0.2 * t + f)))
    y = lp(y, 1400) * np.minimum(1, t / 1.6) * np.clip((d - t) / 4.5, 0, 1)
    return np.stack([y, np.roll(y, int(0.018 * SR))]) * 0.05, 0.6


SYNTHS = dict(heart=heart, gather=gather, drone=drone, glow=glow, riser=riser, whoosh=whoosh, boom=boom,
              tissue=tissue, bubble=bubble, crackle=crackle, dissonance=dissonance, ping=ping, flow=flow,
              reverse=reverse, chime=chime, pad=pad)


def reverb_ir(dur=3.2):
    t = t_axis(dur)
    irs = []
    for ch in range(2):
        n = rng.standard_normal(len(t)) * np.exp(-t / 0.85)
        n = lp(n, 5200) * 0.7 + lp(n, 1400) * 0.3
        n[: int(0.012 * SR)] *= np.linspace(0, 1, int(0.012 * SR))
        irs.append(n / np.sqrt(np.sum(n ** 2)))
    return irs


def main() -> None:
    spec = json.loads((HERE / "kapanis_cues.json").read_text(encoding="utf-8"))
    total, cues = spec["total"], spec["cues"]
    n = int(total * SR)
    dry = np.zeros((2, n))
    send = np.zeros((2, n))
    for c in sorted(cues, key=lambda c: c["t"]):
        sig, rv = SYNTHS[c["type"]](c)
        i = int(c["t"] * SR)
        m = min(sig.shape[1], n - i)
        if m <= 0:
            continue
        dry[:, i:i + m] += sig[:, :m]
        send[:, i:i + m] += sig[:, :m] * rv
    ir = reverb_ir()
    wet = np.stack([fftconvolve(send[k], ir[k])[:n] for k in range(2)]) * 0.9
    mix = dry + wet
    mix = np.tanh(mix * 1.6) / np.tanh(1.6)            # yumuşak sınırlayıcı
    mix /= np.abs(mix).max() / 0.89                    # tepe ≈ -1 dBFS
    fade = int(0.02 * SR)
    mix[:, :fade] *= np.linspace(0, 1, fade)
    mix[:, -int(0.8 * SR):] *= np.linspace(1, 0, int(0.8 * SR))
    pcm = (np.clip(mix.T, -1, 1) * (2 ** 23 - 1)).astype(np.int32)
    b = pcm.astype("<i4").tobytes()
    b24 = bytearray()
    raw = np.frombuffer(b, dtype=np.uint8).reshape(-1, 4)[:, :3]
    b24 = raw.tobytes()
    out = HERE / "kapanis_ses.wav"
    with wave.open(str(out), "wb") as w:
        w.setnchannels(2)
        w.setsampwidth(3)
        w.setframerate(SR)
        w.writeframes(b24)
    rms = np.sqrt(np.mean(mix ** 2))
    print(f"[ok] {out.name} · {total} sn · {len(cues)} ses olayı · RMS {20 * np.log10(rms):.1f} dBFS")


if __name__ == "__main__":
    main()
