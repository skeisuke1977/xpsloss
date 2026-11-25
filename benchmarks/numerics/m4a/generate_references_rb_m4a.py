"""Independent references R-B and CE for M4a (docs/M4A_PREREG.md 6.2 and 6.3).

Written by the reference author R-B, separately from the implementation and from the trial
prototypes of work/consult/m4a/trial/.  The only package function that is called is the public
M3a kernel xpsloss.effective_cross_section_depth_averaged; the public input types
xpsloss.Material and xpsloss.Oscillator are used to build its argument.  No module of the M4a
implementation (xps_single_loss, _loss_panels, _fold) is imported.  The candidate roots, the
partitions, the Gauss-Legendre rules, the overlap of the primary with the output regions and all
sums are computed in this file.

Usage (the rows are computed into per-row cache files, then assembled):

    python benchmarks/numerics/m4a/generate_references_rb_m4a.py --cache DIR --rows RB1,RB1a
    python benchmarks/numerics/m4a/generate_references_rb_m4a.py --cache DIR --assemble
    python benchmarks/numerics/m4a/generate_references_rb_m4a.py --check-rows RB1   (recompute, compare)

The JSON contains no time stamps and no timings, so a new run on the same environment gives the
same bytes.  The timings are written to DIR/timing_*.json and recorded in the appendix of
docs/M4A_PREREG.md.

Quantities (docs/M4A_DESIGN.md 2.3-2.6, prereg 4.1 and 6.2).  K_j(T) is the M3a output at the
binary64 loss T (eV^-1 A^-1), the eleven columns in the order of COLUMNS; `total` is the M3a
column itself (it is integrated directly, never rebuilt from no_hole and hole).

    window_area.j          A_j = Int_a^b K_j dT                                    (A^-1)
    single_loss_weight     w = L A_total
    inverse interval       [1/(r+u_r), 1/(r-u_r)] of window_area.total (prereg 7.2 (i))
    cutoff_sensitivity.j   H_j = L a K_j(a)
    region D zero_loss     J0 = G_D(0), exact rational
    region D single_loss_c J1_c = L Int dE' F(E') Int_{W cap [E'-v, E'-u]} K_c(T) dT
    region D spectrum      J0 + J1_total
    region D H+-_c         L a K_c(a) G_D(a +- 0), exact rational G

The bins are [u, v) and the last bin [u, v]; D_- = (-inf, u_0), D_+ = (u_N, inf).  For a delta
line (E_i, s_i) the outer integral is the value at E_i.  For a binned primary the outer variable
E' runs over each primary bin [p, q] with density f/(q-p); the outer pieces are cut where an end
of the inner interval switches between moving and fixed (E' = a+u, a+v, b+u, b+v).

Partitions (own rule, more conservative than the design's):
  * candidates z: the roots enumerated below, and the extra points T = 0, T = E0 and, for a
    gapped material, T = g;
  * the inner window is first split uniformly so that no panel is wider than Gamma_min/2;
  * every panel P must satisfy dist(z, P) >= |P| for every candidate (bisection at the midpoint),
    which grades the panels geometrically around near candidates and towards 0 and E0;
  * the outer pieces obey the same cap and the same rule for the shifted candidates z+u (moving
    upper end E'-u) and z+v (moving lower end E'-v).
Gauss-Legendre nodes and weights are computed here at 80 digits; the inner nodes are rounded to
binary64 before the kernel call, the outer nodes are never rounded.  All value sums are formed in
mpmath at 50 digits.

Candidate roots (own code; design 3.2 (3) gives the idea).  With z = E0 t (2 - t) the endpoint
momentum is H q^2 = E0 t^2 (k_total) or E0 t^3 (4 - t)/4 (k_parallel) at both ends; Re t < 1 is the
lower end (t = 1 - s/sqrt(E0)) and Re t > 1 the upper end (t = 1 + s/sqrt(E0)), s = sqrt(E0 - z).
D_i = W_i^2 - z^2 - i G_i z, W_i = E_i + d_i H q^2, S = sum A_i/D_i.  The polynomials of the
families bulk (D_i), surface (2 prod D - sum A_i prod_{k!=i} D_k) and epsilon (prod D - sum ...)
are built with exact rational coefficients from the binary64 inputs; only exactly zero leading
coefficients are dropped.  If every active line is dispersionless the polynomials are in z (fixed
poles, end "both").  Gapped materials use the gapless coefficients.  The roots are found by
mpmath.polyroots at 60 digits (extra precision 600 bits) and confirmed on the original rational
expressions with the frozen criteria (prereg 3.1).

Uncertainty u_r of a key (prereg 6.2), every part stored separately:
  (1) inner_order   |V(inner Gauss-Legendre order N_IN_ALT) - V(base)|
  (2) outer_order   |V(outer order N_OUT_ALT) - V(base)| (binned primaries; 0 otherwise)
  (3) local_breaks  |V(shifted partition) - V(base)|: every interior inner breakpoint moved by
                    one quarter of its smaller neighbour panel and, for binned primaries, every
                    outer piece split at 2/5 of its length
  (4) upstream      sum over nodes of |weight| times the M3a per-column error (component_errors)
  (5) rounding      binary64 rounding of the inner nodes (|weight| |dK/dT| |fl(T)-T|, slope from
                    adjacent nodes, times 2), the 50-digit sums and the 40-digit decimal
The value r is the base computation.  f_r = u_r/|r| for r != 0.

The rule "no reference" (prereg 6.3) is applied here, before any implementation value exists.
"""
import os

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_name, "1")

import argparse
from bisect import bisect_right
from fractions import Fraction
import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy

import xpsloss
from xpsloss import Material, Oscillator, effective_cross_section_depth_averaged

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TABLE_FILE = HERE / "expected_states_m4a.json"
TABLE_SHA256 = "250b45f863ccbd1c4c94474f23769976380c5dad5c2671328712ef4acdc7f408"
OUTPUT = HERE / "references_rb_m4a.json"
SCHEMA = "xpsloss-m4a-references-rb-v1"

COLUMNS = ("term_13a", "term_13b", "term_14", "term_15_T1", "term_15_T2", "term_15_T3",
           "term_16_e", "term_16_h", "no_hole", "hole", "total")
COLIDX = {c: i for i, c in enumerate(COLUMNS)}
REGION_COLUMNS = ("total", "no_hole", "hole")
QUANTITY_COLUMN = {"single_loss_total": "total", "single_loss_no_hole": "no_hole",
                   "single_loss_hole": "hole"}

DPS = 50                     # working digits of all sums
GL_DPS = 80                  # digits of the Gauss-Legendre rules
ROOT_DPS = 60                # digits of the root finder
DIGITS = 40                  # significant digits of the stored decimals
N_IN, N_IN_ALT = 16, 24      # inner Gauss-Legendre orders (base, part (1))
N_OUT, N_OUT_ALT = 12, 20    # outer orders (base, part (2))
KAPPA = 1                    # dist(z, P) >= KAPPA |P|
CAP_DIVISOR = 2              # panel width <= Gamma_min / CAP_DIVISOR
SHIFT = Fraction(1, 4)       # part (3): breakpoint shift, fraction of the smaller neighbour
OUTER_SPLIT = Fraction(2, 5)  # part (3): split point of every outer piece
SLOPE_SAFETY = 2
CHUNK = 1500                 # losses per kernel call (plus the two window ends)
KERNEL_KW = ("kernel_epsabs", "kernel_epsrel", "kernel_limit")
# prereg 3.1 (frozen criteria for confirming and selecting candidates)
CONFIRM = mp.mpf("1e-8")
CANCEL = mp.mpf("1e-10")
DUPLICATE = mp.mpf("1e-9")
# prereg 6.3 and 7.2
FR_MAX = Fraction(125, 10 ** 7)      # f_r <= 1.25e-5, i.e. d_rel = 8 f_r <= 1e-4
CE_AREA_FLOOR = Fraction(1, 10 ** 18)
H_FLOOR = Fraction(1, 10 ** 12)
RB_ROWS = ("RB1", "RB2", "RB3", "RB4", "RB5", "RB6", "RB7", "RB8", "RB1a", "RB2a", "RB3a", "RB4a")
CE_ROWS = ("CE1", "CE2", "CE3", "CE4", "CE5", "CE6")
ALL_ROWS = RB_ROWS + CE_ROWS
LEDGER = {"R-B": 1048, "CE": 66}

mp.mp.dps = DPS


# ----------------------------------------------------------------------------- numbers

def mpf_of(x):
    if isinstance(x, Fraction):
        return mp.mpf(x.numerator) / x.denominator
    return mp.mpf(x)


def dec(x):
    """Decimal string with DIGITS significant digits (exact zero is "0")."""
    if isinstance(x, Fraction):
        with mp.workdps(DPS + 20):
            x = mpf_of(x)
    if x == 0:
        return "0"
    return mp.nstr(x, DIGITS, strip_zeros=False, min_fixed=1, max_fixed=0)


def dec_error(x, text):
    """|decimal(text) - x| (the rounding of the stored decimal)."""
    with mp.workdps(DPS + 20):
        exact = mpf_of(x) if isinstance(x, Fraction) else x
        return abs(mp.mpf(text) - exact) if text != "0" else abs(exact)


def up(x):
    """An mpf bound rounded up when written (used for u_r)."""
    return x * (1 + mp.mpf(10) ** (-(DIGITS - 5))) if x else x


def frac_str(x):
    return Fraction(x) if x != "0" else Fraction(0)


def hexval(item):
    v = float.fromhex(item["hex"])
    if v != item["value"]:
        raise ValueError(f"hex {item['hex']} and decimal {item['value']} differ")
    return v


# ----------------------------------------------------------------------------- Gauss-Legendre

_GL = {}


def gauss_legendre(n):
    """Nodes (ascending) and weights of the n-point rule on [-1, 1], own Newton iteration."""
    if n in _GL:
        return _GL[n]
    with mp.workdps(GL_DPS + 20):
        tol = mp.mpf(10) ** (-(GL_DPS + 10))
        pairs = []
        for i in range(1, n + 1):
            x = mp.cos(mp.pi * (4 * i - 1) / (4 * n + 2))
            for _ in range(100):
                p0, p1 = mp.mpf(1), x
                for k in range(2, n + 1):
                    p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
                dp = n * (x * p1 - p0) / (x * x - 1)
                dx = p1 / dp
                x -= dx
                if abs(dx) < tol:
                    break
            else:
                raise ArithmeticError("Gauss-Legendre Newton iteration did not converge")
            p0, p1 = mp.mpf(1), x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1)
            pairs.append((x, 2 / ((1 - x * x) * dp * dp)))
        pairs.sort()
        # self-check: sum of weights 2, exact for x^(2n-2) and x^(2n-1)
        s0 = mp.fsum(w for _, w in pairs)
        s2 = mp.fsum(w * x ** (2 * n - 2) for x, w in pairs)
        s1 = mp.fsum(w * x ** (2 * n - 1) for x, w in pairs)
        if abs(s0 - 2) > tol * 100 or abs(s2 - mp.mpf(2) / (2 * n - 1)) > tol * 100 or abs(s1) > tol * 100:
            raise ArithmeticError("Gauss-Legendre self-check failed")
    with mp.workdps(GL_DPS):
        rule = ([+x for x, _ in pairs], [+w for _, w in pairs])
    _GL[n] = rule
    return rule


def panel_rule(x0, x1, n):
    xs, ws = gauss_legendre(n)
    c, h = (x0 + x1) / 2, (x1 - x0) / 2
    return [c + h * x for x in xs], [h * w for w in ws]


# ----------------------------------------------------------------------------- inputs

def load_table():
    raw = TABLE_FILE.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != TABLE_SHA256:
        raise SystemExit(f"expected_states_m4a.json has SHA-256 {digest}, not the frozen one")
    return json.loads(raw.decode("ascii"))


def row_conditions(table, row):
    d = table["defaults"]
    inp = row["inputs"]

    def val(name):
        x = inp.get(name, d[name])
        return hexval(x) if isinstance(x, dict) else float(x)

    m = table["materials"][inp["material"]]
    osc = [tuple(float.fromhex(h) for h in o) for o in m["oscillators_hex"]]
    for o, ov in zip(osc, m["oscillators"]):
        if list(o) != [float(v) for v in ov]:
            raise ValueError(f"material {inp['material']}: hex and decimal oscillators differ")
    gap = float.fromhex(m["gap_hex"])
    if gap != m["gap_eV"]:
        raise ValueError("gap hex and decimal differ")
    settings = dict(d["settings"])
    settings.update(row["settings"])
    cond = {
        "row": row["id"], "family": row["family"], "entry": row["entry"],
        "material_name": inp["material"], "gap": gap, "oscillators": osc,
        "E0": val("electron_energy_eV"), "L": val("weight_length_A"), "Q": val("hole_charge"),
        "mode": inp.get("epsilon_momentum", d["epsilon_momentum"]), "theta": float(d["theta_deg"]),
        "window": tuple(hexval(x) for x in inp["loss_window_eV"]), "settings": settings,
        "primary": None, "edges": None,
    }
    if "primary" in inp:
        p = table["primaries"][inp["primary"]]
        kind = p["kind"]
        pos = [hexval({"hex": hx, "value": v}) for hx, v in
               zip(p["energies_eV" if kind == "delta" else "edges_eV"]["hex"],
                   p["energies_eV" if kind == "delta" else "edges_eV"]["values"])]
        inten = [hexval({"hex": hx, "value": v}) for hx, v in zip(p["intensities"]["hex"], p["intensities"]["values"])]
        cond["primary"] = {"name": inp["primary"], "kind": kind, "positions": pos, "intensities": inten}
        e = table["output_edges"][inp["output_edges_eV"]]
        cond["edges"] = [hexval({"hex": hx, "value": v}) for hx, v in zip(e["hex"], e["values"])]
        cond["edges_name"] = inp["output_edges_eV"]
    return cond


def make_material(cond):
    return Material(cond["material_name"], cond["gap"],
                    tuple(Oscillator(*o) for o in cond["oscillators"]))


# ----------------------------------------------------------------------------- candidate roots

Z0 = (Fraction(0), Fraction(0))


def padd(p, q):
    n = max(len(p), len(q))
    p = p + [Z0] * (n - len(p))
    q = q + [Z0] * (n - len(q))
    return [(a[0] + b[0], a[1] + b[1]) for a, b in zip(p, q)]


def pscale(p, c):
    cr, ci = c
    return [(a * cr - b * ci, a * ci + b * cr) for a, b in p]


def pmul(p, q):
    out = [Z0] * (len(p) + len(q) - 1)
    for i, (a, b) in enumerate(p):
        if a == 0 and b == 0:
            continue
        for j, (c, d) in enumerate(q):
            r = out[i + j]
            out[i + j] = (r[0] + a * c - b * d, r[1] + a * d + b * c)
    return out


def preal(coeffs):
    return [(Fraction(x), Fraction(0)) for x in coeffs]


def pstrip(p):
    while len(p) > 1 and p[-1] == Z0:
        p = p[:-1]
    return p


def to_mpc(c):
    return mp.mpc(mpf_of(c[0]), mpf_of(c[1]))


def candidate_roots(cond):
    """Enumerate the endpoint roots (or fixed poles) of the families bulk, surface, epsilon."""
    E0 = Fraction(cond["E0"])
    active = [o for o in cond["oscillators"] if o[1] > 0]
    if not active:
        return {"dispersive": False, "variable": None, "families": [], "roots": []}
    dispersive = any(o[3] != 0 for o in active)
    if dispersive:
        zpoly = preal([0, 2 * E0, -E0])
        hq2 = preal([0, 0, E0]) if cond["mode"] == "k_total" else preal([0, 0, 0, E0, -E0 / 4])
    else:
        zpoly = preal([0, 1])
        hq2 = preal([0])
    D = []
    for E, A, G, dd in active:
        W = padd(preal([Fraction(E)]), pscale(hq2, (Fraction(dd), Fraction(0))))
        Di = padd(pmul(W, W), pscale(pmul(zpoly, zpoly), (Fraction(-1), Fraction(0))))
        Di = padd(Di, pscale(zpoly, (Fraction(0), -Fraction(G))))
        D.append(Di)
    prod = [(Fraction(1), Fraction(0))]
    for Di in D:
        prod = pmul(prod, Di)
    ssum = [Z0]
    for i, (E, A, G, dd) in enumerate(active):
        term = [(Fraction(A), Fraction(0))]
        for k, Dk in enumerate(D):
            if k != i:
                term = pmul(term, Dk)
        ssum = padd(ssum, term)
    fams = [("bulk", i, pstrip(Di)) for i, Di in enumerate(D)]
    fams.append(("surface", None, pstrip(padd(pscale(prod, (Fraction(2), Fraction(0))),
                                               pscale(ssum, (Fraction(-1), Fraction(0)))))))
    fams.append(("epsilon", None, pstrip(padd(prod, pscale(ssum, (Fraction(-1), Fraction(0)))))))

    a, b = (mpf_of(Fraction(x)) for x in cond["window"])
    E0m = mpf_of(E0)
    out, fam_info = [], []
    with mp.workdps(ROOT_DPS):
        for family, index, p in fams:
            degree = len(p) - 1
            fam_info.append({"family": family, "line": index, "degree": degree})
            if degree < 1:
                continue
            coeffs = [to_mpc(c) for c in reversed(p)]
            roots, err = mp.polyroots(coeffs, maxsteps=4000, extraprec=600, error=True)
            if err > mp.mpf(10) ** (-(ROOT_DPS - 15)):
                raise ArithmeticError(f"polyroots error {err} for {family}")
            for t in roots:
                # Newton polish on the polynomial itself (well conditioned at 60 digits)
                for _ in range(5):
                    f, df = mp.polyval(coeffs, t, derivative=True)
                    if df == 0:
                        break
                    t -= f / df
                if dispersive:
                    z = E0m * t * (2 - t)
                    hq = E0m * t * t if cond["mode"] == "k_total" else E0m * t ** 3 * (4 - t) / 4
                    end = "lower" if mp.re(t) < 1 else "upper"
                else:
                    z, hq, end = t, mp.mpf(0), "both"
                Ws = [mp.mpf(E) + mp.mpf(dd) * hq for E, A, G, dd in active]
                Ds = [W * W - z * z - 1j * mp.mpf(G) * z for W, (E, A, G, dd) in zip(Ws, active)]
                sig = [abs(W) ** 2 + abs(z) ** 2 + mp.mpf(G) * abs(z) for W, (E, A, G, dd) in zip(Ws, active)]
                if family == "bulk":
                    residual = abs(Ds[index]) / sig[index]
                    cancelled = False
                else:
                    if any(Dk == 0 for Dk in Ds):
                        residual, cancelled = mp.inf, True
                    else:
                        S = mp.fsum(mp.mpf(A) / Dk for (E, A, G, dd), Dk in zip(active, Ds))
                        scale = 1 + mp.fsum(abs(mp.mpf(A) / Dk) for (E, A, G, dd), Dk in zip(active, Ds))
                        residual = abs((2 if family == "surface" else 1) - S) / scale
                        cancelled = any(abs(Dk) <= CANCEL * sk for Dk, sk in zip(Ds, sig))
                x = min(max(mp.re(z), a), b)
                out.append({"z": mp.mpc(z), "family": family, "line": index, "end": end,
                            "residual": residual, "confirmed": bool(residual <= CONFIRM),
                            "cancelled": bool(cancelled), "distance": abs(z - x),
                            "near_origin": bool(abs(z) < a / 2), "kind": "enumerated", "multiplicity": 1})
    # duplicates (same family and end)
    kept = []
    for r in out:
        for k in kept:
            if (k["family"], k["end"]) == (r["family"], r["end"]) and \
                    abs(k["z"] - r["z"]) <= DUPLICATE * (1 + abs(r["z"])):
                k["multiplicity"] += 1
                break
        else:
            kept.append(r)
    for r in kept:
        r["adopted"] = bool((not r["cancelled"]) and r["distance"] < b - a)
    kept.sort(key=lambda r: (float(r["distance"]), r["family"], float(mp.re(r["z"])), float(mp.im(r["z"]))))
    return {"dispersive": dispersive, "variable": "t" if dispersive else "z", "families": fam_info,
            "roots": kept}


def extra_points(cond):
    pts = [("origin", mp.mpc(0)), ("E0_endpoint", mp.mpc(mpf_of(Fraction(cond["E0"]))))]
    if cond["gap"] > 0:
        pts.append(("gap_edge", mp.mpc(mpf_of(Fraction(cond["gap"])))))
    return pts


def partition_candidates(cond, roots):
    a, b = (mpf_of(Fraction(x)) for x in cond["window"])
    zs = [r["z"] for r in roots["roots"] if not r["cancelled"]]
    zs += [z for _, z in extra_points(cond)]
    return zs


# ----------------------------------------------------------------------------- partitions

def dist(z, x0, x1):
    xr = mp.re(z)
    xc = x0 if xr < x0 else (x1 if xr > x1 else xr)
    return abs(z - xc)


def refine(edges, cands, cap, scale):
    """Split uniformly to width <= cap, then bisect until dist(z, P) >= KAPPA |P|."""
    out = []
    for x0, x1 in zip(edges[:-1], edges[1:]):
        n = int(mp.ceil((x1 - x0) / cap))
        pts = [x0 + (x1 - x0) * k / n for k in range(n)] + [x1]
        work = list(zip(pts[:-1], pts[1:]))[::-1]
        while work:
            p0, p1 = work.pop()
            w = p1 - p0
            near = [z for z in cands if abs(mp.re(z) - (p0 + p1) / 2) < 2 * w + abs(mp.im(z)) + w]
            if all(dist(z, p0, p1) >= KAPPA * w for z in near):
                out.append((p0, p1))
            else:
                if w < scale * mp.mpf(10) ** -12:
                    raise ArithmeticError("partition: a candidate lies on the integration path")
                m = (p0 + p1) / 2
                work.extend([(m, p1), (p0, m)])
    pts = [out[0][0]] + [p1 for _, p1 in out]
    return pts


def shifted(points):
    """Part (3): move every interior breakpoint by SHIFT of its smaller neighbour panel."""
    s = mpf_of(SHIFT)
    new = [points[0]]
    for k in range(1, len(points) - 1):
        left, right = points[k] - points[k - 1], points[k + 1] - points[k]
        new.append(points[k] + s * min(left, right))
    new.append(points[-1])
    return new


def gamma_cap(cond):
    widths = [o[2] for o in cond["oscillators"] if o[1] > 0]
    return mpf_of(Fraction(min(widths))) / CAP_DIVISOR


# ----------------------------------------------------------------------------- regions

def regions(cond):
    """[(name, u, v)] with u, v as Fractions or None for -inf / +inf."""
    e = [Fraction(x) for x in cond["edges"]]
    out = [("outside.low", None, e[0])]
    out += [(f"bin[{i}]", e[i], e[i + 1]) for i in range(len(e) - 1)]
    out.append(("outside.high", e[-1], None))
    return out


def zero_loss(cond, name, u, v):
    p = cond["primary"]
    nb = len(cond["edges"]) - 1
    last = name == f"bin[{nb - 1}]"
    total = Fraction(0)
    if p["kind"] == "delta":
        for E, s in zip(p["positions"], p["intensities"]):
            E, s = Fraction(E), Fraction(s)
            if u is None:
                inside = E < v
            elif v is None:
                inside = E > u
            else:
                inside = (u <= E < v) or (last and E == v)
            total += s if inside else 0
    else:
        for (lo, hi), f in zip(zip(p["positions"][:-1], p["positions"][1:]), p["intensities"]):
            lo, hi, f = Fraction(lo), Fraction(hi), Fraction(f)
            x0 = lo if u is None else max(lo, u)
            x1 = hi if v is None else min(hi, v)
            if x1 > x0:
                total += f * (x1 - x0) / (hi - lo)
    return total


def g_at_a(cond, name, u, v, side):
    """G_D(a+0) (side "plus") or G_D(a-0) (side "minus"), exact."""
    p = cond["primary"]
    a = Fraction(cond["window"][0])
    total = Fraction(0)
    if p["kind"] == "delta":
        for E, s in zip(p["positions"], p["intensities"]):
            x = Fraction(E) - a
            if side == "plus":       # detected energy x - 0
                inside = (x <= v) if u is None else ((x > u) if v is None else (u < x <= v))
            else:                    # detected energy x + 0
                inside = (x < v) if u is None else ((x >= u) if v is None else (u <= x < v))
            total += Fraction(s) if inside else 0
    else:                            # continuous in T: overlap of [p, q] with D + a
        for (lo, hi), f in zip(zip(p["positions"][:-1], p["positions"][1:]), p["intensities"]):
            lo, hi, f = Fraction(lo), Fraction(hi), Fraction(f)
            x0 = lo if u is None else max(lo, u + a)
            x1 = hi if v is None else min(hi, v + a)
            if x1 > x0:
                total += f * (x1 - x0) / (hi - lo)
    return total


# ----------------------------------------------------------------------------- integration plans

def inner_panels(lo, hi, points):
    """Panels of [lo, hi] cut at the interior points of the partition."""
    i0 = bisect_right(points, lo)
    i1 = bisect_right(points, hi)
    cuts = [lo] + [x for x in points[i0:i1] if lo < x < hi] + [hi]
    return list(zip(cuts[:-1], cuts[1:]))


def add_inner(acc, lo, hi, weight, points, n):
    for x0, x1 in inner_panels(lo, hi, points):
        T, W = panel_rule(x0, x1, n)
        acc.append((T, [weight * w for w in W]))


def outer_pieces(cond, u, v, p, q, cands, cap, split):
    """Outer pieces of [p, q] for the region (u, v): switch points, cap and tolerance rule."""
    a, b = (mpf_of(Fraction(x)) for x in cond["window"])
    sw = []
    for base in (u, v):
        if base is not None:
            sw += [mpf_of(Fraction(cond["window"][0]) + base), mpf_of(Fraction(cond["window"][1]) + base)]
    pm, qm = mpf_of(p), mpf_of(q)
    pts = sorted({pm, qm} | {s for s in sw if pm < s < qm})
    pieces = []
    for x0, x1 in zip(pts[:-1], pts[1:]):
        m = (x0 + x1) / 2
        lo = a if v is None else max(a, m - mpf_of(v))
        hi = b if u is None else min(b, m - mpf_of(u))
        if lo >= hi:
            continue
        shift = []
        if v is not None and m - mpf_of(v) > a:
            shift += [z + mpf_of(v) for z in cands]
        if u is not None and m - mpf_of(u) < b:
            shift += [z + mpf_of(u) for z in cands]
        sub = refine([x0, x1], shift, cap, x1 - x0) if shift else refine([x0, x1], [], cap, x1 - x0)
        for y0, y1 in zip(sub[:-1], sub[1:]):
            if split:
                c = y0 + (y1 - y0) * mpf_of(OUTER_SPLIT)
                pieces += [(y0, c), (c, y1)]
            else:
                pieces.append((y0, y1))
    return pieces


def plan(cond, points, n_in, n_out, split, cands, cap):
    """{integral name: [(T list, W list)]} for one variant."""
    a, b = (mpf_of(Fraction(x)) for x in cond["window"])
    L = mpf_of(Fraction(cond["L"]))
    out = {}
    acc = []
    add_inner(acc, a, b, mp.mpf(1), points, n_in)
    out["area"] = acc
    if cond["primary"] is None:
        return out
    prim = cond["primary"]
    for name, u, v in regions(cond):
        acc = []
        if prim["kind"] == "delta":
            for E, s in zip(prim["positions"], prim["intensities"]):
                if s == 0:
                    continue
                Em = mpf_of(Fraction(E))
                lo = a if v is None else max(a, Em - mpf_of(v))
                hi = b if u is None else min(b, Em - mpf_of(u))
                if lo < hi:
                    add_inner(acc, lo, hi, L * mpf_of(Fraction(s)), points, n_in)
        else:
            xs, ws = gauss_legendre(n_out)
            for (p, q), f in zip(zip(prim["positions"][:-1], prim["positions"][1:]), prim["intensities"]):
                if f == 0:
                    continue
                p, q, f = Fraction(p), Fraction(q), Fraction(f)
                rho = mpf_of(f / (q - p))
                for y0, y1 in outer_pieces(cond, u, v, p, q, cands, cap, split):
                    c, h = (y0 + y1) / 2, (y1 - y0) / 2
                    for x, w in zip(xs, ws):
                        Ep = c + h * x
                        lo = a if v is None else max(a, Ep - mpf_of(v))
                        hi = b if u is None else min(b, Ep - mpf_of(u))
                        if lo < hi:
                            add_inner(acc, lo, hi, L * rho * h * w, points, n_in)
        out[name] = acc
    return out


# ----------------------------------------------------------------------------- kernel

class Kernel:
    """Cached calls of the public M3a kernel; every batch contains the two window ends."""

    def __init__(self, cond):
        self.cond = cond
        self.material = make_material(cond)
        self.cache = {}
        self.calls = 0
        self.points = 0
        a, b = cond["window"]
        self.sentinels = (a, b)
        s = cond["settings"]
        self.kw = dict(weight_length_A=cond["L"], epsilon_momentum=cond["mode"], theta_deg=cond["theta"],
                       hole_charge=cond["Q"], epsabs=s["kernel_epsabs"], epsrel=s["kernel_epsrel"],
                       limit=s["kernel_limit"])

    def call(self, losses):
        r = effective_cross_section_depth_averaged(self.material, self.cond["E0"], np.asarray(losses, float), **self.kw)
        vals = np.stack([np.asarray(getattr(r, c), float) for c in COLUMNS], axis=1)
        errs = np.stack([np.asarray(r.component_errors[c], float) for c in COLUMNS], axis=1)
        self.calls += 1
        self.points += len(losses)
        return vals, errs

    def need(self, losses):
        todo = sorted({t for t in losses if t not in self.cache})
        for i in range(0, len(todo), CHUNK):
            batch = todo[i:i + CHUNK]
            full = list(self.sentinels) + [t for t in batch if t not in self.sentinels]
            vals, errs = self.call(full)
            for t, v, e in zip(full, vals, errs):
                if t in self.cache:
                    if not (np.array_equal(self.cache[t][0], v) and np.array_equal(self.cache[t][1], e)):
                        raise ArithmeticError(f"M3a value at {t!r} depends on the batch")
                else:
                    self.cache[t] = (v, e)

    def ensure_sentinels(self):
        if any(t not in self.cache for t in self.sentinels):
            vals, errs = self.call(list(self.sentinels))
            for t, v, e in zip(self.sentinels, vals, errs):
                self.cache.setdefault(t, (v, e))

    def batch_independence(self):
        """Self-check: a few cached losses again in a different batch, bitwise equal."""
        keys = sorted(self.cache)
        probe = keys[1:len(keys):max(1, len(keys) // 5)][:5]
        if not probe:
            return {"checked": 0, "bitwise_equal": True}
        vals, errs = self.call(probe[::-1] + list(self.sentinels))
        ok = all(np.array_equal(self.cache[t][0], v) and np.array_equal(self.cache[t][1], e)
                 for t, v, e in zip(probe[::-1], vals, errs))
        return {"checked": len(probe), "bitwise_equal": bool(ok)}


def integrate(panels, cols, kernel):
    """Values (mpf), upstream error and node-rounding/arithmetic rounding for the columns."""
    idx = [COLIDX[c] for c in cols]
    allW, allV = [], [[] for _ in cols]
    upstream = np.zeros(len(cols))
    node = np.zeros(len(cols))
    absvsum = np.zeros(len(cols))
    deferred = []
    gslope = np.zeros(len(cols))
    for T, W in panels:
        Tf = [float(t) for t in T]
        V = np.array([kernel.cache[t][0][idx] for t in Tf])
        E = np.array([kernel.cache[t][1][idx] for t in Tf])
        aw = np.array([float(abs(w)) for w in W])
        upstream += aw @ E
        absvsum += aw @ np.abs(V)
        dT = np.array([float(abs(mp.mpf(tf) - t)) for tf, t in zip(Tf, T)])
        x = np.array(Tf)
        dx = np.diff(x)
        ok = dx > 0
        if ok.any():
            q = np.zeros((len(x) - 1, len(cols)))
            q[ok] = np.abs(np.diff(V, axis=0)[ok]) / dx[ok, None]
            slope = np.zeros((len(x), len(cols)))
            slope[:-1] = q
            slope[1:] = np.maximum(slope[1:], q)
            slope *= SLOPE_SAFETY
            gslope = np.maximum(gslope, slope.max(axis=0))
            node += (aw * dT) @ slope
        else:
            deferred.append(float(np.sum(aw * dT)))
        allW.extend(W)
        for k in range(len(cols)):
            allV[k].extend(V[:, k].tolist())
    for wd in deferred:
        node += wd * gslope
    with mp.workdps(DPS):
        values = [mp.fdot(allW, allV[k]) if allW else mp.mpf(0) for k in range(len(cols))]
    arith = absvsum * 10.0 ** (-(DPS - 5))
    return values, upstream * (1 + 1e-12), node + arith


# ----------------------------------------------------------------------------- one row

def compute_row(table, row_id, kernels):
    row = next(r for r in table["rows"] if r["id"] == row_id)
    cond = row_conditions(table, row)
    t_wall, t_cpu = time.perf_counter(), time.process_time()
    ckey = (cond["material_name"], cond["E0"], cond["L"], cond["Q"], cond["mode"], cond["window"])
    kernel = kernels.setdefault(ckey, Kernel(cond))
    pts_before, calls_before = kernel.points, kernel.calls

    tr = time.process_time()
    roots = candidate_roots(cond)
    cands = partition_candidates(cond, roots)
    a, b = (mpf_of(Fraction(x)) for x in cond["window"])
    cap = gamma_cap(cond)
    base_points = refine([a, b], cands, cap, b - a)
    alt_points = shifted(base_points)
    t_roots = time.process_time() - tr

    binned = cond["primary"] is not None and cond["primary"]["kind"] == "binned"
    variants = {
        "base": (base_points, N_IN, N_OUT, False),
        "inner_order": (base_points, N_IN_ALT, N_OUT, False),
        "local_breaks": (alt_points, N_IN, N_OUT, True),
    }
    if binned:
        variants["outer_order"] = (base_points, N_IN, N_OUT_ALT, False)
    results, nodes = {}, {}
    failure = None
    try:
        for vname, (pts, nin, nout, split) in variants.items():
            pl = plan(cond, pts, nin, nout, split, cands, cap)
            losses = [float(t) for panels in pl.values() for T, _ in panels for t in T]
            nodes[vname] = len(losses)
            kernel.need(losses)
            kernel.ensure_sentinels()
            res = {}
            for name, panels in pl.items():
                cols = COLUMNS if name == "area" else REGION_COLUMNS
                res[name] = dict(zip(("values", "upstream", "rounding"), integrate(panels, cols, kernel)))
                res[name]["cols"] = cols
            results[vname] = res
        kernel.ensure_sentinels()
        indep = kernel.batch_independence()
    except Exception as exc:  # reference path failure (prereg 6.3 (5))
        failure = {"type": type(exc).__name__, "message": str(exc)}
        indep = None

    keys = build_keys(cond, row, results, kernel, failure)
    entry = {
        "row": row_id, "family": row["family"], "entry": row["entry"],
        "conditions": conditions_record(cond),
        "candidates": roots_record(roots, cond),
        "partition": {
            "cap_eV": dec(cap), "kappa": KAPPA, "inner_panels": len(base_points) - 1,
            "min_panel_width_eV": dec(min(y - x for x, y in zip(base_points[:-1], base_points[1:]))),
            "breakpoints_eV": [dec(x) for x in base_points],
        },
        "integration_nodes": nodes,
        "keys": keys,
        "reference_path_failure": failure,
        "self_checks": self_checks(cond, results, indep) if failure is None else None,
    }
    timing = {"row": row_id, "wall_s": time.perf_counter() - t_wall, "cpu_s": time.process_time() - t_cpu,
              "roots_and_partition_cpu_s": t_roots, "kernel_points_new": kernel.points - pts_before,
              "kernel_calls_new": kernel.calls - calls_before}
    return entry, timing


def conditions_record(cond):
    out = {
        "material": cond["material_name"], "gap_eV": cond["gap"].hex(),
        "oscillators_hex": [[float(x).hex() for x in o] for o in cond["oscillators"]],
        "electron_energy_eV": cond["E0"].hex(), "weight_length_A": cond["L"].hex(),
        "hole_charge": cond["Q"].hex(), "epsilon_momentum": cond["mode"], "theta_deg": cond["theta"].hex(),
        "loss_window_eV": [x.hex() for x in cond["window"]],
        "kernel_settings": {k: cond["settings"][k] for k in KERNEL_KW},
        "comparison_settings": {k: cond["settings"][k] for k in ("epsrel", "area_epsabs", "spectrum_epsabs_fraction")},
    }
    if cond["primary"] is not None:
        p = cond["primary"]
        out["primary"] = {"name": p["name"], "kind": p["kind"], "positions_eV": [x.hex() for x in p["positions"]],
                          "intensities": [x.hex() for x in p["intensities"]]}
        out["output_edges_eV"] = [x.hex() for x in cond["edges"]]
    return out


def roots_record(roots, cond):
    rec = []
    for r in roots["roots"]:
        if not (r["adopted"] or r["distance"] < 4 * (mpf_of(Fraction(cond["window"][1])) - mpf_of(Fraction(cond["window"][0])))):
            continue
        rec.append({"z_re": mp.nstr(mp.re(r["z"]), 20, min_fixed=1, max_fixed=0),
                    "z_im": mp.nstr(mp.im(r["z"]), 20, min_fixed=1, max_fixed=0),
                    "family": r["family"], "line": r["line"], "end": r["end"],
                    "residual": mp.nstr(r["residual"], 3), "confirmed": r["confirmed"], "cancelled": r["cancelled"],
                    "near_origin": r["near_origin"], "multiplicity": r["multiplicity"], "adopted": r["adopted"],
                    "distance_to_window_eV": mp.nstr(r["distance"], 15, min_fixed=1, max_fixed=0)})
    extras = [{"name": n, "z_re": mp.nstr(mp.re(z), 20), "z_im": "0"} for n, z in extra_points(cond)]
    return {"variable": roots["variable"], "dispersive": roots["dispersive"], "families": roots["families"],
            "method": "exact rational coefficients, mpmath.polyroots at 60 digits (extraprec 600), "
                      "Newton polish on the polynomial, residual on the original rational expressions",
            "near_roots": rec, "extra_partition_points": extras,
            "adopted_count": sum(1 for r in roots["roots"] if r["adopted"])}


# ----------------------------------------------------------------------------- keys

def comparison_floor(cond, key):
    s = cond["settings"]
    L, a = Fraction(cond["L"]), Fraction(cond["window"][0])
    if key.startswith("window_area."):
        return CE_AREA_FLOOR if cond["family"] == "CE" else Fraction(s["area_epsabs"])
    if key == "single_loss_weight":
        return L * Fraction(s["area_epsabs"])
    if key.startswith("cutoff_sensitivity."):
        return L * a * H_FLOOR
    sabs = sum(abs(Fraction(x)) for x in cond["primary"]["intensities"])
    if ".cutoff_sensitivity." in key:
        return L * a * H_FLOOR * sabs
    return Fraction(s["spectrum_epsabs_fraction"]) * sabs


def classify(cond, key, r, u):
    """prereg 6.3: status, condition, f_r, d_abs, d_rel, width (all exact)."""
    d_abs = comparison_floor(cond, key)
    epsrel = Fraction(cond["settings"]["epsrel"])
    if r == 0:
        if u <= d_abs / 4:
            return "value", "6.3(1): r = 0, u_r <= d_abs/4; width d_abs", None, d_abs, None, d_abs
        return "no_reference", "6.3(1): r = 0 and u_r > d_abs/4", None, d_abs, None, d_abs
    f = u / abs(r)
    d_rel = max(8 * f, epsrel)
    if f <= FR_MAX:
        return "value", "6.3(4): r != 0, f_r <= 1.25e-5; d_rel = max(8 f_r, epsrel)", f, d_abs, d_rel, d_abs + d_rel * abs(r)
    return "no_reference", "6.3(4): f_r > 1.25e-5 (d_rel would exceed 1e-4)", f, d_abs, d_rel, d_abs + d_rel * abs(r)


def entry_for(cond, key, value, parts, unit):
    """value: mpf or Fraction; parts: dict of mpf (or Fraction) -> stored entry."""
    text = dec(value)
    parts = dict(parts)
    parts["rounding"] = parts.get("rounding", mp.mpf(0)) + dec_error(value, text)
    u = mp.fsum(mpf_of(p) if isinstance(p, Fraction) else p for p in parts.values())
    parts_text = {k: dec(up(mpf_of(p) if isinstance(p, Fraction) else p)) for k, p in parts.items()}
    u_text = dec(up(u))
    r_f, u_f = frac_str(text), frac_str(u_text)
    status, condition, f, d_abs, d_rel, width = classify(cond, key, r_f, u_f)
    return {"kind": "direct_reference", "unit": unit, "r": text, "u_r": u_text,
            "u_r_parts": parts_text, "f_r": None if f is None else dec(f),
            "d_abs": dec(d_abs), "d_rel": None if d_rel is None else dec(d_rel), "width": dec(width),
            "status": status, "condition": condition}


ZERO_PARTS = ("inner_order", "outer_order", "local_breaks", "upstream", "rounding")


def integral_parts(results, name, col, scale=1):
    base = results["base"][name]
    k = base["cols"].index(col)
    v0 = base["values"][k]
    parts = {}
    for part in ("inner_order", "outer_order", "local_breaks"):
        if part in results:
            parts[part] = abs(results[part][name]["values"][k] - v0) * abs(scale)
        else:
            parts[part] = mp.mpf(0)
    parts["upstream"] = mp.mpf(float(base["upstream"][k])) * abs(scale)
    parts["rounding"] = mp.mpf(float(base["rounding"][k])) * abs(scale)
    return v0 * scale, parts


def build_keys(cond, row, results, kernel, failure):
    planned = [k.split("|") for k in row["check_plan"]]
    planned = [(k, kind) for k, kind, _ in planned if kind in ("direct_reference", "inverse_interval")]
    keys = {}
    if failure is not None:
        for k, kind in planned:
            keys[k] = {"kind": kind, "status": "reference_path_failure", "condition": "6.3(5)", "failure": failure}
        return keys
    L = mpf_of(Fraction(cond["L"]))
    a = cond["window"][0]
    am = mpf_of(Fraction(a))
    Ka, Ea = kernel.cache[a]
    reg = {n: (u, v) for n, u, v in regions(cond)} if cond["primary"] is not None else {}
    for k, kind in planned:
        if kind == "inverse_interval":
            continue
        if k.startswith("window_area."):
            col = k.split(".", 1)[1]
            v, parts = integral_parts(results, "area", col)
            keys[k] = entry_for(cond, k, v, parts, "A^-1")
        elif k == "single_loss_weight":
            v, parts = integral_parts(results, "area", "total", scale=L)
            keys[k] = entry_for(cond, k, v, parts, "1")
        elif k.startswith("cutoff_sensitivity."):
            col = k.split(".", 1)[1]
            j = COLIDX[col]
            v = L * am * mp.mpf(float(Ka[j]))
            parts = {p: mp.mpf(0) for p in ZERO_PARTS}
            parts["upstream"] = L * am * mp.mpf(float(Ea[j]))
            keys[k] = entry_for(cond, k, v, parts, "1")
        else:
            region, rest = split_region(k)
            u, v_ = reg[region]
            if rest == "zero_loss":
                g = zero_loss(cond, region, u, v_)
                keys[k] = entry_for(cond, k, g, {p: mp.mpf(0) for p in ZERO_PARTS}, "I")
            elif rest in QUANTITY_COLUMN:
                val, parts = integral_parts(results, region, QUANTITY_COLUMN[rest])
                keys[k] = entry_for(cond, k, val, parts, "I")
            elif rest == "spectrum":
                val, parts = integral_parts(results, region, "total")
                g = zero_loss(cond, region, u, v_)
                with mp.workdps(DPS + 20):
                    total = val + mpf_of(g)
                keys[k] = entry_for(cond, k, total, parts, "I")
            elif rest.startswith("cutoff_sensitivity."):
                _, side, col = rest.split(".")
                g = g_at_a(cond, region, u, v_, side)
                j = COLIDX[col]
                with mp.workdps(DPS + 20):
                    val = L * am * mp.mpf(float(Ka[j])) * mpf_of(g)
                parts = {p: mp.mpf(0) for p in ZERO_PARTS}
                parts["upstream"] = L * am * mp.mpf(float(Ea[j])) * abs(mpf_of(g))
                keys[k] = entry_for(cond, k, val, parts, "I")
            else:
                raise KeyError(k)
    tot = keys.get("window_area.total")
    for k, kind in planned:
        if kind != "inverse_interval":
            continue
        if tot is None:
            raise KeyError("inverse interval without window_area.total")
        r, u = frac_str(tot["r"]), frac_str(tot["u_r"])
        e = {"kind": "inverse_interval", "unit": "A", "derived_from": "window_area.total",
             "r_area": tot["r"], "u_area": tot["u_r"]}
        if tot["status"] != "value":
            e.update(status="no_reference", condition="window_area.total has no reference")
        elif r - u <= 0:
            e.update(status="structure", condition="6.3(2): r - u_r <= 0, compare value, interval and reason")
        else:
            with mp.workdps(DPS + 20):
                lo = 1 / mpf_of(r + u)
                hi = 1 / mpf_of(r - u)
            lo_t = mp.nstr(lo * (1 - mp.mpf(10) ** (-(DIGITS - 3))), DIGITS, strip_zeros=False, min_fixed=1, max_fixed=0)
            hi_t = mp.nstr(hi * (1 + mp.mpf(10) ** (-(DIGITS - 3))), DIGITS, strip_zeros=False, min_fixed=1, max_fixed=0)
            e.update(status="value", condition="prereg 7.2 (i): the reference interval [1/(r+u_r), 1/(r-u_r)] "
                     "(outward rounded) of window_area.total",
                     reference_interval_A=[lo_t, hi_t])
        keys[k] = e
    return keys


def split_region(key):
    if key.startswith("outside.low."):
        return "outside.low", key[len("outside.low."):]
    if key.startswith("outside.high."):
        return "outside.high", key[len("outside.high."):]
    i = key.index("].")
    return key[:i + 1], key[i + 2:]


def self_checks(cond, results, indep):
    base = results["base"]
    out = {"batch_independence": indep}
    area = base["area"]
    # component identity of the reference itself: total - (no_hole + Q hole), no_hole and hole sums
    v = dict(zip(area["cols"], area["values"]))
    Q = mpf_of(Fraction(cond["Q"]))
    ident = {
        "total_minus_no_hole_plus_Q_hole": dec(v["total"] - (v["no_hole"] + Q * v["hole"])),
        "no_hole_minus_terms": dec(v["no_hole"] - (v["term_13a"] + v["term_13b"] + v["term_15_T1"] + v["term_15_T3"] + v["term_16_e"])),
        "hole_minus_terms": dec(v["hole"] - (v["term_14"] + v["term_15_T2"] + v["term_16_h"])),
        "scale_sum_abs_upstream_total": dec(mp.mpf(float(area["upstream"][COLIDX["total"]]))),
    }
    out["area_component_identity"] = ident
    if cond["primary"] is not None:
        S = sum(Fraction(x) for x in cond["primary"]["intensities"])
        L = mpf_of(Fraction(cond["L"]))
        bal = {}
        for c in REGION_COLUMNS:
            tot = mp.fsum(base[n]["values"][REGION_COLUMNS.index(c)] for n, _, _ in regions(cond))
            ref = mpf_of(S) * L * v[c]
            bal[c] = {"sum_regions": dec(tot), "S_L_A": dec(ref), "difference": dec(tot - ref)}
        out["balance_not_used_for_passing"] = bal
    return out


# ----------------------------------------------------------------------------- assemble

def ledger_counts(table):
    per_row, fam = {}, {"R-B": 0, "CE": 0}
    for r in table["rows"]:
        if r["id"] not in ALL_ROWS:
            continue
        n = sum(1 for k in r["check_plan"] if k.split("|")[1] in ("direct_reference", "inverse_interval"))
        per_row[r["id"]] = n
        fam["CE" if r["family"] == "CE" else "R-B"] += n
    return per_row, fam


def width_of(entry):
    return frac_str(entry["width"])


def mutation_record(rows, row, key, y, uy):
    ref = rows[row]["keys"][key]
    r, u, w = frac_str(ref["r"]), frac_str(ref["u_r"]), width_of(ref)
    d = abs(y - r)
    return {"row": row, "key": key, "y_mutant": dec(y), "u_mutant": dec(uy), "r": ref["r"], "u_r": ref["u_r"],
            "width": ref["width"], "lhs_|y-r|+u_r": dec(d + u), "ratio_lhs_to_width": dec((d + u) / w),
            "conservative_ratio_(|y-r|-u_mutant)/width": dec((d - uy) / w),
            "exceeds_width": bool(d + u > w), "exceeds_width_conservatively": bool(d - uy > w)}


def pre_comparison_records(rows):
    rec = {}
    ident = []
    for pair, keys in ((("RB1", "RB2"), ("bin[1].single_loss_total", "bin[2].single_loss_total")),
                       (("RB3", "RB4"), ("bin[0].single_loss_total", "bin[2].single_loss_total"))):
        for k in keys:
            e1, e2 = rows[pair[0]]["keys"][k], rows[pair[1]]["keys"][k]
            r1, r2 = frac_str(e1["r"]), frac_str(e2["r"])
            u1, u2 = frac_str(e1["u_r"]), frac_str(e2["u_r"])
            w1, w2 = width_of(e1), width_of(e2)
            diff = abs(r1 - r2)
            ident.append({"rows": list(pair), "key": k, "r": [e1["r"], e2["r"]], "width": [e1["width"], e2["width"]],
                          "difference": dec(diff), "relative_difference": dec(diff / abs(r1)),
                          "ratio_difference_to_width_sum": dec(diff / (w1 + w2)),
                          "ratio_conservative_(diff-u1-u2)/(w1+w2)": dec((diff - u1 - u2) / (w1 + w2)),
                          "exceeds": bool(diff - u1 - u2 > w1 + w2)})
    rec["identification_bins"] = ident
    # M04: the M3a length stays 20 A, only the outer multiplier follows L (5 A in RB5)
    L5, L20 = Fraction(5), Fraction(20)
    m04 = []
    for k in ("bin[1].single_loss_total", "bin[2].single_loss_total"):
        e = rows["RB1"]["keys"][k]
        y = frac_str(e["r"]) * L5 / L20
        m04.append(mutation_record(rows, "RB5", k, y, frac_str(e["u_r"]) * L5 / L20))
    all_m04 = {"exceeding": 0, "compared": 0}
    for k, e in rows["RB5"]["keys"].items():
        if e["kind"] != "direct_reference" or e["status"] != "value":
            continue
        e1 = rows["RB1"]["keys"][k]
        r1 = frac_str(e1["r"])
        if k.startswith("window_area."):
            y = r1                                  # the area uses the M3a length 20
        elif k.endswith(".zero_loss"):
            y = r1
        elif k == "single_loss_weight":
            y = r1 * L5 / L20                       # 5 * A(20)
        elif k.endswith(".spectrum"):
            r0 = frac_str(rows["RB1"]["keys"][k.replace(".spectrum", ".zero_loss")]["r"])
            y = r0 + (r1 - r0) * L5 / L20
        else:
            y = r1 * L5 / L20                       # multiplier 5 on a 20 A kernel
        all_m04["compared"] += 1
        r, u, w = frac_str(e["r"]), frac_str(e["u_r"]), width_of(e)
        all_m04["exceeding"] += int(abs(y - r) + u > w)
    rec["M04"] = {"detected_by_row": "RB5", "keys": m04, "all_RB5_value_keys": all_m04,
                  "normal_control": {"row": "RB1", "note": "RB1 already uses 20 A, so under M04 the call and the "
                                     "values equal the normal version (no width excess is required)"}}
    m05 = []
    for row in ("RB6", "RB7"):
        k = "bin[1].single_loss_total"
        e = rows["RB1"]["keys"][k]               # M3a total with the default unit charge
        m05.append(mutation_record(rows, row, k, frac_str(e["r"]), frac_str(e["u_r"])))
    rec["M05a"] = {"keys": m05, "note": "without hole_charge upstream the M3a total is no_hole + 1*hole, i.e. "
                   "the RB1 value of the same key (same material, window, L, primary and edges)"}
    return rec


PREREG_CANDIDATES = {
    # prereg 5.2 (CE rows): family, end, distance to the window (eV) and the tolerance written there
    "CE1": [("surface", "lower", "0.0222", "0.05")],
    "CE2": [("bulk", "lower", "0.02536", "0.01")],
    "CE3": [("bulk", "lower", "0.00503", "0.01")],
}
PREREG_POSITIONS = {
    # prereg 5.2: surface candidates within 1e-3 eV (complex distance) of these positions
    "CE5": [("surface", "8.3811", "-0.5"), ("surface", "10.0130", "-0.5")],
    "CE6": [("surface", "8.3516", "-0.5"), ("surface", "9.9875", "-0.5")],
}


def candidate_checks(rows):
    """Own candidate roots against the values written in the pre-registration (record only)."""
    out = {}
    for rid in CE_ROWS:
        near = rows[rid]["candidates"]["near_roots"]
        rec = []
        for family, end, dist, tol in PREREG_CANDIDATES.get(rid, []):
            hits = [c for c in near if c["family"] == family and c["end"] == end]
            best = min(hits, key=lambda c: abs(Fraction(c["distance_to_window_eV"]) - Fraction(dist)))
            rel = abs(Fraction(best["distance_to_window_eV"]) - Fraction(dist)) / Fraction(dist)
            rec.append({"family": family, "end": end, "prereg_distance_eV": dist, "own_distance_eV":
                        best["distance_to_window_eV"], "relative_difference": dec(rel),
                        "within_prereg_tolerance": bool(rel <= Fraction(tol))})
        for family, re_, im_ in PREREG_POSITIONS.get(rid, []):
            target = mp.mpc(mp.mpf(re_), mp.mpf(im_))
            hits = [c for c in near if c["family"] == family]
            best = min(hits, key=lambda c: abs(mp.mpc(mp.mpf(c["z_re"]), mp.mpf(c["z_im"])) - target))
            d = abs(mp.mpc(mp.mpf(best["z_re"]), mp.mpf(best["z_im"])) - target)
            rec.append({"family": family, "prereg_position_eV": [re_, im_], "own_root_eV": [best["z_re"], best["z_im"]],
                        "complex_distance_eV": mp.nstr(d, 6), "within_1e-3": bool(d <= mp.mpf("1e-3"))})
        bulk_both = sum(1 for c in near if c["family"] == "bulk" and c["end"] == "both" and c["adopted"])
        rec.append({"adopted_bulk_poles_end_both": bulk_both})
        out[rid] = rec
    return out


def assemble(cache):
    table = load_table()
    rows = {}
    for rid in ALL_ROWS:
        f = cache / f"row_{rid}.json"
        if not f.exists():
            raise SystemExit(f"missing {f}")
        rows[rid] = json.loads(f.read_text(encoding="ascii"))
    per_row, fam = ledger_counts(table)
    if fam != LEDGER:
        raise SystemExit(f"ledger counts {fam} differ from {LEDGER}")
    counts = {"by_row": {}, "by_family": {"R-B": 0, "CE": 0}}
    no_ref, failures = [], []
    status_counts = {}
    for rid, entry in rows.items():
        n = len(entry["keys"])
        if n != per_row[rid]:
            raise SystemExit(f"{rid}: {n} keys, ledger {per_row[rid]}")
        counts["by_row"][rid] = n
        counts["by_family"]["CE" if entry["family"] == "CE" else "R-B"] += n
        for k, e in entry["keys"].items():
            status_counts[e["status"]] = status_counts.get(e["status"], 0) + 1
            if e["status"] == "no_reference":
                no_ref.append({"row": rid, "key": k, "condition": e["condition"], "f_r": e.get("f_r")})
            if e["status"] == "reference_path_failure":
                failures.append({"row": rid, "key": k, "failure": e.get("failure")})
    counts["by_status"] = dict(sorted(status_counts.items()))
    counts["ledger"] = LEDGER
    summary = maxima(rows)
    doc = {
        "schema": SCHEMA,
        "prereg": "docs/M4A_PREREG.md (6.2, 6.3, 7.2, appendix R-B and CE)",
        "table": {"file": "benchmarks/numerics/m4a/expected_states_m4a.json", "sha256": TABLE_SHA256},
        "generator": {"file": "benchmarks/numerics/m4a/generate_references_rb_m4a.py",
                      "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        "kernel": "xpsloss.effective_cross_section_depth_averaged (public M3a API), the only package function called",
        "environment": environment(),
        "method": {
            "inner_gauss_legendre": [N_IN, N_IN_ALT], "outer_gauss_legendre": [N_OUT, N_OUT_ALT],
            "kappa": KAPPA, "cap": "Gamma_min/2", "shift": str(SHIFT), "outer_split": str(OUTER_SPLIT),
            "digits": {"sums": DPS, "gauss_legendre": GL_DPS, "roots": ROOT_DPS, "stored": DIGITS},
            "u_r_parts": list(ZERO_PARTS),
            "value": "base variant",
            "partition": "uniform split to width <= Gamma_min/2, then bisection until dist(z, P) >= |P| for every "
                         "candidate z (enumerated roots and the points 0, E0 and the gap edge); outer pieces cut at "
                         "the switch points a+u, a+v, b+u, b+v with the same rule for z+u and z+v",
            "kernel_settings": "kernel_epsabs, kernel_epsrel, kernel_limit of the row (API defaults 1e-12, 1e-9, 200); "
                               "every call contains the two window ends",
        },
        "counts": counts,
        "maxima": summary,
        "no_reference": no_ref,
        "reference_path_failures": failures,
        "pre_comparison_records": pre_comparison_records(rows),
        "candidate_checks_against_prereg": candidate_checks(rows),
        "rows": rows,
    }
    text = json.dumps(doc, indent=1, sort_keys=True, ensure_ascii=True) + "\n"
    return text


def maxima(rows):
    best_f = (Fraction(-1), None, None)
    best = {p: (Fraction(-1), None, None) for p in ZERO_PARTS}
    for rid, entry in rows.items():
        for k, e in entry["keys"].items():
            if e["kind"] != "direct_reference" or e["status"] == "reference_path_failure":
                continue
            if e["f_r"] is not None:
                f = frac_str(e["f_r"])
                if f > best_f[0]:
                    best_f = (f, rid, k)
            r = frac_str(e["r"])
            for p in ZERO_PARTS:
                x = frac_str(e["u_r_parts"][p])
                rel = x / abs(r) if r else None
                if rel is not None and rel > best[p][0]:
                    best[p] = (rel, rid, k)
    out = {"f_r": {"value": dec(best_f[0]), "row": best_f[1], "key": best_f[2]}}
    out["u_r_part_relative"] = {p: {"value": dec(v), "row": rid, "key": k} for p, (v, rid, k) in best.items()}
    return out


def environment():
    return {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
            "mpmath": mp.__version__, "xpsloss": xpsloss.__version__,
            "platform": platform.platform(), "machine": platform.machine(),
            "base_commit": "38ab37a (main, M4a design decided and checks pre-registered)"}


# ----------------------------------------------------------------------------- main

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cache", type=Path)
    ap.add_argument("--rows")
    ap.add_argument("--assemble", action="store_true")
    ap.add_argument("--check-rows")
    args = ap.parse_args(argv)
    if args.rows:
        table = load_table()
        args.cache.mkdir(parents=True, exist_ok=True)
        kernels, timings = {}, []
        t0, c0 = time.perf_counter(), time.process_time()
        for rid in args.rows.split(","):
            entry, timing = compute_row(table, rid, kernels)
            (args.cache / f"row_{rid}.json").write_text(json.dumps(entry, indent=1, sort_keys=True, ensure_ascii=True) + "\n",
                                                        encoding="ascii", newline="\n")
            timings.append(timing)
            print(json.dumps(timing), flush=True)
        total = {"rows": args.rows, "wall_s": time.perf_counter() - t0, "cpu_s": time.process_time() - c0,
                 "kernel_points": sum(k.points for k in kernels.values()),
                 "kernel_calls": sum(k.calls for k in kernels.values()), "per_row": timings}
        (args.cache / f"timing_{args.rows.replace(',', '_')}.json").write_text(json.dumps(total, indent=1) + "\n",
                                                                               encoding="ascii", newline="\n")
        print(json.dumps({k: v for k, v in total.items() if k != "per_row"}), flush=True)
    if args.assemble:
        text = assemble(args.cache)
        OUTPUT.write_text(text, encoding="ascii", newline="\n")
        print("wrote", OUTPUT, hashlib.sha256(text.encode("ascii")).hexdigest())
    if args.check_rows:
        table = load_table()
        frozen = json.loads(OUTPUT.read_text(encoding="ascii"))
        kernels = {}
        for rid in args.check_rows.split(","):
            entry, _ = compute_row(table, rid, kernels)
            same = json.dumps(entry, sort_keys=True) == json.dumps(frozen["rows"][rid], sort_keys=True)
            print(rid, "identical" if same else "DIFFERENT")
            if not same:
                return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
