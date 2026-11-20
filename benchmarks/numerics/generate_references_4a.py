"""Regenerate the literal reference tables from the repository root with
python benchmarks/numerics/generate_references_4a.py (mpmath 1.3, outside CI).

Numbers are the actual binary64 inputs lifted to 80 decimals. Bulk uses a
complex partial-fraction primitive in dispersed energy. For the hard-gap KK
and optical moments, integrate the real-axis kernel on width-resolved panels
and integrate the remaining infinite tail. Remove the PV pole algebraically.
The triangle uses its actual rounded knots. No xpsloss code is imported.
"""

from functools import lru_cache
import json
import math
from pathlib import Path

import mpmath as mp


mp.mp.dps = 80
ROOT = Path(__file__).resolve().parents[2]


@lru_cache(None)
def bulk(w, a, gamma, alpha, t=20):
    w, a, gamma, alpha = map(mp.mpf, (w, a, gamma, alpha))
    e = 1000
    plus = mp.sqrt(e) + mp.sqrt(e - t)
    lo, hi = (w + alpha * z**2 for z in (t / plus, plus))
    root = mp.sqrt(t**2 + 1j * gamma * t)
    poles = [w, root, -root]
    primitive = mp.fsum(
        (mp.log(hi - p) - mp.log(lo - p)) / mp.fprod(p - q for q in poles if q != p)
        for p in poles
    )
    return a * mp.im(primitive) / (2 * mp.pi * mp.mpf(0.529177210903) * e)


def response(width, gap, loss, kind):
    w, a, width, gap, t = map(mp.mpf, (15, 100, width, gap, loss))

    def denominator(u):
        return ((w - u) * (w + u))**2 + (width * u)**2

    def integrand(u):
        if kind == "first":
            return a * width * u**2 / denominator(u)
        if kind == "lhs" or t == 0:
            return 2 * a * width / (mp.pi * denominator(u))
        # Subtract t²/D(t) from u²/D(u) before dividing by u²-t².
        return (2 * a * width / mp.pi) * (w**4 - u**2 * t**2) / (denominator(u) * denominator(t))

    cutoff = 64
    nodes = {gap, cutoff}
    if gap < w:
        nodes.add(w)
    flank = width / 2
    while flank < 16 * cutoff:
        nodes.update(p for p in (w - flank, w + flank) if gap < p < cutoff)
        flank *= 16
    value = mp.quad(integrand, sorted(nodes)) + mp.quad(integrand, [cutoff, mp.inf])
    if kind == "kk":
        value -= 1
        if t:
            value += a * width * t / (mp.pi * denominator(t)) * mp.log(abs((gap + t) / (gap - t)))
    return float(value)


def triangle(width):
    lo, w, hi, height = map(mp.mpf, (15 - width, 15, 15 + width, 1 / width))

    def kernel(u):
        level = (u - lo) / (w - lo) if u <= w else (hi - u) / (hi - w)
        return 2 * height * level * u / (mp.pi * (u**2 - 400))

    return float(-1 + mp.quad(kernel, [lo, w, hi]))


if __name__ == "__main__":
    broad = bulk(15, 10, 20, 0.5)
    for width in (1, 0.1, 0.01, 0.001, 1e-4, 1e-5, 1e-6, 1e-7, 1e-9, 1e-12):
        k = bulk(15, 180, width, 0.5)
        print("width", width, float(k), float(k + broad))
    for strength in (0.001, 0.01, 0.1, 1, 10, 50, 100, 180, 1000):
        print("strength", strength, float(bulk(15, strength, 1e-5, 0.5) + broad))
    for alpha, w in ((0.5, 19.949493661166535), (0.001, 16.040101012677667)):
        for width in (1e-7, 1e-9):
            for n in (-5, -0.5, 0, 0.5, 5, 100):
                center = w + n * width
                print("edge", center, alpha, width, float(bulk(center, 180, width, alpha) + broad))
    for strength in (1e-8, 180):
        print("13a", strength, float(bulk(15, strength, 1e-24, 0.5) + broad))
    for path in sorted((ROOT / "benchmarks/inputs").glob("*.bulk.json")):
        rows = json.loads(path.read_text())["material"]["oscillators"]
        total = mp.fsum(bulk(o["energy_eV"], o["strength_eV2"], o["width_eV"], o["dispersion"])
                        for o in rows)
        print(path.name, float(total))
        if path.name.startswith("au_"):
            at_97 = mp.fsum(bulk(o["energy_eV"], o["strength_eV2"], o["width_eV"],
                                 o["dispersion"], 97) for o in rows)
            print("Au at 97 eV", float(at_97))
    for width in (1e-14, 1e-16):
        t = mp.mpf(math.nextafter(15.0, math.inf))
        print("neighbor", width, complex(-1 + 100 / (225 - t**2 - 1j * mp.mpf(width) * t)))
    for width in (1e-4, 1e-6, 1e-8, 1e-12):
        print("triangle", width, triangle(width))
    for width, gap, t in ((1e-5, 5, 20), (1e-5, 16, 20), (1e-5, 14.99995, 14.9999500001),
                          (1e-12, 15, 20), (1e-24, 15, 20)):
        print("kk", width, gap, t, response(width, gap, t, "kk"))
    for width, gap in ((1e-5, 5), (1e-10, 15.00000000005), (1e-24, 15)):
        print("optical", width, gap, *(response(width, gap, 0, k) for k in ("lhs", "first")))
