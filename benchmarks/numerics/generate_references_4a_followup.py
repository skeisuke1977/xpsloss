"""Print the follow-up literals with python benchmarks/numerics/generate_references_4a_followup.py.

Frozen with NumPy 2.2.6, SciPy 1.16.3, mpmath 1.3.0 and Python 3.12.10.
Other Gaussian-rule/reduction versions can differ by a few ulps. For the rounded Al samples
use direct 32-point Gaussian quadrature on every interpolation segment, with
the PV value subtracted and integrated separately. Integrate the optional
infinite tail in x=150/u. For gapped oscillators use an 80-decimal partial
fraction primitive of the retained real-axis kernel. No xpsloss imports.
"""
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre


def al_table(step, losses=None):
    root = Path(__file__).resolve().parents[2]
    source = json.loads((root / "benchmarks/inputs/al_elfdb_v1.0.json").read_text(encoding="utf-8"))
    osc = source["material"]["oscillators"][0]
    w, a, g = (osc[key] for key in ("energy_eV", "strength_eV2", "width_eV"))
    nodes = np.round(np.arange(0, 150 + step/2, step), 10)
    values = a * g * nodes / ((w*w - nodes**2)**2 + (g*nodes)**2)
    x, weight = roots_legendre(32)
    fraction = (x + 1) / 2
    widths = np.diff(nodes)[:, None]
    u = nodes[:-1, None] + widths * fraction
    level = values[:-1, None] + np.diff(values)[:, None] * fraction
    results = []
    if losses is None:
        losses = np.arange(0.5, 150, 0.5)
    for t in losses:
        at_t = values[np.searchsorted(nodes, t)]
        density = (level - at_t) / (u - t) + level / (u + t)
        pv = np.sum(density * (widths * (weight / 2))) + at_t * math.log((150 - t) / t)
        tail = 2 * values[-1] * quad(lambda y: y*y / (1 - (t*y/150)**2), 0, 1,
                                    epsabs=1e-13, epsrel=1e-13)[0]
        results.append((-1 + pv / math.pi, -1 + (pv + tail) / math.pi))
    return results


STATIC = ((7.878141405147874, 18.425072384979863, 0.11939570572947698, 0.05),
          (15.8296550770744, 243.3360003362837, 9.096397682802841, 0.0))
NEAR = ((5.884967692641011, 49.344009033807374, 0.3028168533899312, 0.09443850922792511),)
BROAD = ((2.302194719940513, 98.49158300581752, 38.787012565194765, 0.0),
         (4.657905346631697, 0.06897270027955346, 1.00925757845866, 0.0))
SCAN_28 = ((3.556680710239554, 71.07923360993777, 26.068995055720393, 0),
           (22.567507487432778, 58.6259905571535, 8.42521634641305, 0))
DISPERSED = ((4.017220767437674, 2.497700902738755, 9.185425063805724e-8, 0.03712174177269153),)


def kk(lines, gap, t, q):
    gap, t, q = map(mp.mpf, (gap, t, q))
    total = mp.mpf(-1)
    for w, a, g, alpha in lines:
        w, a, g, alpha = map(mp.mpf, (w, a, g, alpha))
        w += alpha * mp.mpf(3.80998212) * q*q
        if not gap:
            total += a / (w*w - t*t - 1j*g*t)
            continue
        h = mp.sqrt(w*w - (g/2)**2)
        poles = [s*h + 1j*r*g/2 for s in (-1, 1) for r in (-1, 1)]
        if t:
            poles.extend([t, -t])
        residues = [(p**2 if t else 1) / mp.fprod(p - k for k in poles if k != p) for p in poles]
        total -= 2*a*g/mp.pi * mp.re(mp.fsum(c * mp.log(gap-p) for c, p in zip(residues, poles)))
    return total


if __name__ == "__main__":
    for coarse, fine in zip(al_table(0.1), al_table(0.05)):
        print("    " + repr(tuple(float(x) for x in coarse + fine)) + ",")
    mp.mp.dps = 80
    for name, lines, gap, ts, q in (
            ("static", STATIC, 7.89193759455495, [0, 3], 0),
            ("near", NEAR, 5.846831168012581, [5.887660937002207], 0),
            ("overdamped", BROAD, 1.550890771276142, [0], 0.9009708343237725),
            ("scan_28", SCAN_28, 0.6671472392854397, [0.67], 0.5),
            ("dispersed", DISPERSED, 0, [4.02099068386561], 0.16326392005751011),
            ("dispersed_gap", DISPERSED, 1, [4.02099068386561], 0.16326392005751011)):
        for t in ts:
            print(name, t, mp.nstr(kk(lines, gap, t, q), 20))
