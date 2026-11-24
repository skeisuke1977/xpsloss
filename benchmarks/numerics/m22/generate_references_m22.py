"""Independent reference values for M2-2 (docs/M2_2_PREREG.md, section 4).

Run offline (about ten minutes on 12 processes):

    PYTHONPATH=src python benchmarks/numerics/m22/generate_references_m22.py

Written by the reference author, separately from the implementation.  Nothing from
xpsloss.dielectric, xpsloss.xps_fixed_depth or the internals of xpsloss.model is used.
The only package import happens inside the frozen tests/xps_route_b_reference.py
(two physical constants), which is called for the I3-type route.

Route "closed_form" (mpmath, all outputs, both momentum conventions)
  * Im F = L (E1, docs/EQUATIONS.md) with the hard gap, F = -1/eps.
  * Re F = -1 + (2/pi) PV Int_g^inf u L0(q,u)/(u^2-T^2) du (E4a).  The integral is done
    by partial fractions written here: u L0/(u^2-T^2) = Im[A u/((W^2-u^2-i G u)(u^2-T^2))]
    for real u, a rational function with simple poles u1, u2 (lower half plane) and +-T.
    A few points are compared with a direct principal-value quadrature (mp.quad).
  * Eqs. (13)-(16) of Simonsen, Yubero, Tougaard, PRB 56, 1612 (1997) (journal pp.
    1614-1615, PDF pp. 3-4), typed from the printed braces as transcribed in
    docs/M1_EQUATION_MAP.md section 4.1, at theta = 0, hbar w > 0.  The dielectric factors
    are built from F: 1/eps-1 = -S, (eps-1)/(eps+1) = R = S/(2-S),
    (eps-1)(eps+2)/(eps(eps+1)) = C = 2S-R, S = 1+F.  Eq. (17) gives the limits.
  * k_parallel integration: adaptive composite 24-point Gauss-Legendre in ln k.
    Run A: 50 digits, adaptive panels.  Run B: 80 digits, every panel of run A halved.
    Value = run B.  u_j = 2 |A - B| + (accepted adaptive residual of run A).

Route "i3" (binary64, k_parallel only)
  * The k_perp integrals of the frozen tests/xps_route_b_reference.py (KPerp.integrals)
    are combined with coefficients formed directly from h = -F (no reciprocal of 1/eps),
    F from the KK code of this file, and integrated over k_parallel by a separate
    composite Gauss-Legendre rule.  Refinements: k_parallel panels halved; k_perp panels
    halved with higher orders in the core and the tails.

V15: the uniform upper bound of the design (section 3.2) for Si, in exact rationals.
"""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from multiprocessing import Pool
import os
from pathlib import Path
import platform
import sys
import time

import mpmath as mp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]

H_DEC = "3.80998212"            # hbar^2/2m [eV A^2], docs/EQUATIONS.md E1
A0_DEC = "0.529177210903"       # Bohr radius [A], docs/EQUATIONS.md E2
E0_DEC = "1000"
DEPTHS_DEC = ("2", "10", "50")
SI_INPUT = "benchmarks/inputs/si_elfdb_v1.0.json"
SI_INPUT_SHA256 = "49020931fcb6d83f799b6a8b0cfd0a6eedb8c3a5dc446b88793bc9e4e5704079"

# (E_i, A_i, Gamma_i, alpha_i) as decimal strings: eV, eV^2, eV, dimensionless.
MATERIALS = {
    "Si": {"gap_eV": "1.12", "oscillators": [("10", "6.08458", "5", "0.5"),
                                               ("14", "30.4229", "5", "0.5"),
                                               ("16.8", "212.9603", "3.8", "0.5")]},
    "S0": {"gap_eV": "1.12", "oscillators": [("10", "6.08458", "5", "0"),
                                               ("14", "30.4229", "5", "0"),
                                               ("16.8", "212.9603", "3.8", "0")]},
    "N1": {"gap_eV": "2", "oscillators": [("3", "6", "0.02", "0.5")]},
}
# Losses of section 4 as decimals (g = 1.12): label, decimal string.
LOSSES = [("0.2", "0.2"), ("0.4", "0.4"), ("1.0", "1.0"),
          ("g-1e-3", "1.119"), ("g-1e-6", "1.119999"), ("g+1e-6", "1.120001"), ("g+1e-3", "1.121"),
          ("1.2", "1.2"), ("2", "2"), ("5", "5"), ("11", "11"), ("16.8", "16.8"), ("17", "17"),
          ("30", "30"), ("60", "60"), ("99.8", "99.8")]
N1_LOSS = ("4", "4")
N1_DEPTH = "10"
SPECS = ("k_parallel", "k_total")
TERMS = ("term_13a", "term_13b", "term_14", "term_15_T1", "term_15_T2",
         "term_15_T3", "term_16_e", "term_16_h")
NO_HOLE_TERMS = ("term_13a", "term_13b", "term_15_T1", "term_15_T3", "term_16_e")
ELECTRON_GROUP = ("term_13b", "term_15_T1", "term_15_T3", "term_16_e")
HOLE_TERMS = ("term_14", "term_15_T2", "term_16_h")
HOLE_CHARGES = (("total_q0", 0.0), ("total_q0.5", 0.5), ("total_q1", 1.0))
EPSABS, EPSREL, LIMIT = 1e-12, 1e-9, 200

DPS_A, DPS_B = 50, 80
GL_DEGREE = 4                    # mpmath Gauss-Legendre degree 4 = 24 nodes
PANEL_TOL = mp.mpf("1e-24")      # accepted |parent - children| per panel and output
MAX_LEVEL = 40


def b64(dec):
    return float(dec)


def number(dec):
    x = float(dec)
    return {"decimal": dec, "binary64_hex": x.hex(), "binary64_repr": repr(x)}


# ----------------------------------------------------------------------------- response

class Mat:
    """Material with the binary64 values of the decimal inputs, as exact mpf numbers."""

    def __init__(self, name):
        spec = MATERIALS[name]
        self.name = name
        self.g = mp.mpf(b64(spec["gap_eV"]))
        self.osc = [tuple(mp.mpf(b64(v)) for v in o) for o in spec["oscillators"]]
        self.dispersive = any(o[3] != 0 and o[1] != 0 for o in self.osc)
        self._cache = {}


def kk_line(E, A, G, alpha, g, q2, T, h):
    """Int_g^inf u L0_i(q,u)/(u^2-T^2) du (principal value) for one oscillator, T > 0, T != g.

    f(u) = A u/((W^2-u^2-iGu)(u^2-T^2)) = -A u/((u-u1)(u-u2)(u-T)(u+T)), Im f = u L0/(u^2-T^2)
    on the real axis.  With f = sum_p c_p/(u-p) and sum_p c_p = 0,
    PV Int_g^inf f du = -sum_p c_p log(g-p), principal log for u1, u2 (Im(g-p) > 0),
    log|g-p| for p = +-T.
    """
    W = E + alpha * h * q2
    root = mp.sqrt(4 * W * W - G * G)          # imaginary when overdamped
    j = mp.mpc(0, 1)
    u1, u2 = (-j * G + root) / 2, (-j * G - root) / 2
    poles = (u1, u2, T, -T)
    total = mp.mpc(0)
    for i, p in enumerate(poles):
        den = mp.mpc(1)
        for k, pp in enumerate(poles):
            if k != i:
                den *= p - pp
        c = -A * p / den
        total += c * (mp.log(g - p) if i < 2 else mp.log(abs(g - p)))
    return -mp.im(total)


def response(mat, q2, T, h):
    """(Re S, Im S) with S = 1 + F: Re S = (2/pi) sum_i PV Int (E4a), Im S = L (E1)."""
    key = None
    if not mat.dispersive:
        key = (mp.mp.dps, T)
        if key in mat._cache:
            return mat._cache[key]
    re = mp.mpf(0)
    im = mp.mpf(0)
    for E, A, G, alpha in mat.osc:
        if A == 0:
            continue
        re += kk_line(E, A, G, alpha, mat.g, q2, T, h)
        if T > mat.g:
            W = E + alpha * h * q2
            im += A * G * T / ((W * W - T * T) ** 2 + G * G * T * T)
    out = (2 * re / mp.pi, im)
    if key is not None:
        mat._cache[key] = out
    return out


def kk_direct(mat, q2, T, h):
    """Re S by direct quadrature of E4a (check of kk_line only)."""
    g = mat.g

    def B(u):   # u L0(q, u)
        s = mp.mpf(0)
        for E, A, G, alpha in mat.osc:
            W = E + alpha * h * q2
            s += A * G * u * u / ((W * W - u * u) ** 2 + G * G * u * u)
        return s

    scale = max(o[0] + o[3] * h * q2 for o in mat.osc)
    marks = sorted({g + (scale - g) * mp.mpf(f) for f in ("0.25", "0.5", "0.75", "1", "1.5", "2", "4", "16")})
    if T < g:
        pts = [g] + [m for m in marks if m > g] + [mp.inf]
        val = mp.quad(lambda u: B(u) / (u * u - T * T), pts)
    else:
        d = T - g
        hfun = lambda u: B(u) / (u + T)
        # PV over [g, 2T-g] by the antisymmetric subtraction, then the regular rest.
        core = mp.quad(lambda s: (hfun(T + s) - hfun(T - s)) / s, mp.linspace(0, d, 5))
        up = T + d
        pts = [up] + [m for m in marks if m > up] + [mp.inf]
        val = core + mp.quad(lambda u: B(u) / (u * u - T * T), pts)
    return 2 * val / mp.pi


# ----------------------------------------------------------------------------- Eqs. (13)-(17)

def kinematics(T, E0, h):
    """b = Omega0 = omega_perp = T/(hbar v), hbar v = 2 sqrt(H E0); Eq. (17) at theta = 0:
    k_par+- = sqrt(k+-^2 - b^2), k+- = (sqrt(E0) +- sqrt(E0 - T))/sqrt(H)."""
    b = T / (2 * mp.sqrt(h * E0))
    s, r = mp.sqrt(E0), mp.sqrt(E0 - T)
    lo = mp.sqrt(((s - r) ** 2) / h - b * b)
    hi = mp.sqrt(((s + r) ** 2) / h - b * b)
    return b, lo, hi


def term_densities(k, S, b, depths, E0, a0):
    """The eight terms per unit k_parallel (before the measure k dk), for each depth.

    Printed forms at theta = 0 (v_par = 0, Omega0 = omega_perp = b, E_perp = E0, x = a),
    hbar w > 0 so that pi delta(omega_perp) + i/omega_perp -> i/b.
    """
    j = mp.mpc(0, 1)
    R = S / (2 - S)
    C = 2 * S - R
    n = b * b + k * k
    p, m = mp.mpc(b, k), mp.mpc(b, -k)
    P1 = 1 / (mp.pi * a0 * E0)
    out = []
    for a in depths:
        P2 = 1 / (2 * mp.pi * E0 * a0 * a * b)
        E = mp.expj(b * a)               # exp(+i Omega0 a)
        Ei = 1 / E
        d = mp.exp(-k * a)
        # Eq. (13), first term: cos/(pi a0 E_perp) Re[i (1/eps-1)/(O^2+k^2)]
        t13a = P1 * mp.re(j * (-S) / n)
        # Eq. (13), second term (three pieces of the brace)
        brace = (j * p * p * (Ei * d - 1) + (-j) * m * m * (E * d - 1)
                 + j * (Ei - d) * (E - d) * n)
        t13b = -P2 * mp.re(j * (-S) / (n * n) * brace)
        # Eq. (14)
        brace = j * p * (1 - d * Ei) + (-j) * m * (d * d - d * Ei)
        t14 = -P2 * mp.re((j / b) * (-S) / n * brace)
        # Eq. (15): common factor (-i)(O - i k)(e^{-ka} - e^{-iOa})/(O^2+k^2)
        X = (-j) * m * (d - Ei) / n
        t15_1 = -P2 * mp.re(j * X * C * (E - d) * p / n)
        t15_2 = -P2 * mp.re(j * X * C * (-j) * (j / b) * d)
        t15_3 = -P2 * mp.re(j * X * R * m / n * E)
        # Eq. (16): (1-eps)/(1+eps) = -R
        Y = (-R) * p * Ei / n
        t16_e = P2 * mp.re(Y * (2 * j * k * E - p * d) / n)
        t16_h = P2 * mp.re(Y * (-j) * (j / b) * d)
        out.extend((t13a, t13b, t14, t15_1, t15_2, t15_3, t16_e, t16_h))
    return out


class Integrand:
    def __init__(self, material, spec, T_float, depths_float):
        self.mat = Mat(material)
        self.spec = spec
        self.T = mp.mpf(T_float)
        self.depths = [mp.mpf(a) for a in depths_float]
        self.h = mp.mpf(b64(H_DEC))
        self.a0 = mp.mpf(b64(A0_DEC))
        self.E0 = mp.mpf(b64(E0_DEC))
        self.b, self.lo, self.hi = kinematics(self.T, self.E0, self.h)
        self.evaluations = 0

    def q2(self, k):
        return k * k + (self.b * self.b if self.spec == "k_total" else 0)

    def S(self, k):
        re, im = response(self.mat, self.q2(k), self.T, self.h)
        return mp.mpc(re, im) if im != 0 else re

    def __call__(self, t):
        self.evaluations += 1
        k = mp.exp(t)
        dens = term_densities(k, self.S(k), self.b, self.depths, self.E0, self.a0)
        k2 = k * k                           # k dk = k^2 d(ln k)
        return [v * k2 for v in dens]

    def feature_points(self):
        """ln k of bulk lines W_i(q) = T and of the largest |1/(1-F)| on a scan."""
        t0, t1 = mp.log(self.lo), mp.log(self.hi)
        pts = []
        for E, A, G, alpha in self.mat.osc:
            if A != 0 and alpha != 0 and self.T > E:
                k2 = (self.T - E) / (alpha * self.h) - (self.b ** 2 if self.spec == "k_total" else 0)
                if k2 > 0:
                    pts.append(mp.log(k2) / 2)
        if self.mat.dispersive:
            n = 400
            grid = [t0 + (t1 - t0) * i / n for i in range(n + 1)]
            w = [abs(1 / (2 - self.S(mp.exp(t)))) for t in grid]
            for i in range(1, n):
                if w[i] > w[i - 1] and w[i] >= w[i + 1] and w[i] > 2:
                    pts.append(grid[i])
        return [t for t in pts if t0 < t < t1]


def gl_nodes(dps):
    return mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(GL_DEGREE, mp.mp.prec)


def panel(f, nodes, lo, hi):
    half, mid = (hi - lo) / 2, (hi + lo) / 2
    acc = None
    for x, w in nodes:
        v = f(mid + half * x)
        if acc is None:
            acc = [w * c for c in v]
        else:
            acc = [s + w * c for s, c in zip(acc, v)]
    return [half * s for s in acc]


def integrate_adaptive(f, edges):
    """Bisect until |I(parent) - I(left) - I(right)| <= PANEL_TOL for every output.
    Returns the sum over accepted children, the accepted children and the residual sum."""
    nodes = gl_nodes(mp.mp.dps)
    stack = [(lo, hi, panel(f, nodes, lo, hi), 0) for lo, hi in zip(edges[:-1], edges[1:])]
    total, leaves, resid = None, [], None
    while stack:
        lo, hi, whole, level = stack.pop()
        mid = (lo + hi) / 2
        left, right = panel(f, nodes, lo, mid), panel(f, nodes, mid, hi)
        diff = [abs(w - l - r) for w, l, r in zip(whole, left, right)]
        if max(diff) <= PANEL_TOL or level >= MAX_LEVEL:
            if level >= MAX_LEVEL and max(diff) > PANEL_TOL:
                raise RuntimeError("adaptive quadrature reached the level limit")
            part = [l + r for l, r in zip(left, right)]
            total = part if total is None else [s + c for s, c in zip(total, part)]
            resid = diff if resid is None else [s + c for s, c in zip(resid, diff)]
            leaves.extend(((lo, mid), (mid, hi)))
        else:
            stack.append((lo, mid, left, level + 1))
            stack.append((mid, hi, right, level + 1))
    return total, sorted(leaves), resid


def integrate_fixed(f, panels):
    nodes = gl_nodes(mp.mp.dps)
    total = None
    for lo, hi in panels:
        part = panel(f, nodes, lo, hi)
        total = part if total is None else [s + c for s, c in zip(total, part)]
    return total


def closed_form_task(args):
    """One (material, spec, loss): all depths.  Returns plain strings (80-digit values)."""
    material, spec, T_float, depths_float, brute = args
    start = time.perf_counter()
    mp.mp.dps = DPS_A
    f = Integrand(material, spec, T_float, depths_float)
    t0, t1 = mp.log(f.lo), mp.log(f.hi)
    n0 = int(mp.ceil((t1 - t0) / mp.mpf("0.5")))
    edges = sorted(set([t0 + (t1 - t0) * i / n0 for i in range(n0 + 1)] + f.feature_points()))
    run_a, leaves, resid = integrate_adaptive(f, edges)
    evals_a = f.evaluations

    mp.mp.dps = DPS_B
    fb = Integrand(material, spec, T_float, depths_float)
    tb0, tb1 = mp.log(fb.lo), mp.log(fb.hi)
    fine = []
    for lo, hi in leaves:
        lo, hi = mp.mpf(lo), mp.mpf(hi)
        fine.extend(((lo, (lo + hi) / 2), ((lo + hi) / 2, hi)))
    # The end points are recomputed at 80 digits (the 50-digit ones carry 1e-50 errors).
    fine[0] = (tb0, fine[0][1])
    fine[-1] = (fine[-1][0], tb1)
    run_b = integrate_fixed(fb, fine)

    brute_diff = None
    if brute:
        # Independent partition: uniform panels of width <= 0.002 in ln k, 50 digits.
        mp.mp.dps = DPS_A
        fc = Integrand(material, spec, T_float, depths_float)
        n = int(mp.ceil((t1 - t0) / mp.mpf("0.002")))
        uniform = [(t0 + (t1 - t0) * i / n, t0 + (t1 - t0) * (i + 1) / n) for i in range(n)]
        run_c = integrate_fixed(fc, uniform)
        mp.mp.dps = DPS_B
        brute_diff = [mp.nstr(abs(mp.mpf(c) - b), 5) for c, b in zip(run_c, run_b)]

    mp.mp.dps = DPS_B
    out = {"material": material, "spec": spec, "loss": T_float, "depths": depths_float,
           "values": [mp.nstr(v, 40) for v in run_b],
           "diff_ab": [mp.nstr(abs(mp.mpf(a) - b), 5) for a, b in zip(run_a, run_b)],
           "residual_a": [mp.nstr(r, 5) for r in resid],
           "panels_a": len(leaves), "evaluations_a": evals_a, "evaluations_b": fb.evaluations,
           "k_limits": [mp.nstr(fb.lo, 25), mp.nstr(fb.hi, 25)], "omega_perp": mp.nstr(fb.b, 25),
           "brute_diff": brute_diff, "seconds": time.perf_counter() - start}
    return out


# ----------------------------------------------------------------------------- I3-type route

def i3_task(args):
    """Binary64 I3-type values for Si and S0 at one (loss, depth, variant), k_parallel."""
    T_float, a_float, variant = args
    import numpy as np
    from numpy.polynomial.legendre import leggauss
    from numpy.polynomial.laguerre import laggauss
    sys.path.insert(0, str(ROOT / "tests"))
    import xps_route_b_reference as ref

    saved = (ref._XG, ref._WG, ref._XL, ref._WL, ref._XS, ref._WS)

    class FineKPerp(ref.KPerp):
        """Same integrands; core panels of half the width."""

        def _core(self, q, K, extra=()):
            hmax = 0.5 * np.pi / self.a
            e = [0.0]
            x = 0.0
            while x < K:
                x = min(x + min(0.125 * np.hypot(x, q), hmax), K)
                e.append(x)
            e = np.asarray(e)
            e = np.concatenate([-e[::-1], e[1:]])
            if extra:
                e = np.unique(np.concatenate([e, np.asarray(extra, dtype=float)]))
            lo, hi = e[:-1], e[1:]
            k = (0.5 * (hi - lo)[:, None] * ref._XG + 0.5 * (hi + lo)[:, None]).ravel()
            wts = (0.5 * (hi - lo)[:, None] * ref._WG).ravel()
            return k, wts

    try:
        if variant == "kperp_refined":
            # Higher orders in the frozen module's rules: core 28, Laguerre tails 120, 1/k tail 48.
            ref._XG, ref._WG = leggauss(28)
            ref._XL, ref._WL = laggauss(120)
            xs, ws = leggauss(48)
            ref._XS, ref._WS = 0.5 * (xs + 1.0), 0.5 * ws
            kp = FineKPerp(float(E0_DEC), T_float, a_float)
        else:
            kp = ref.KPerp(float(E0_DEC), T_float, a_float)
        width = 0.175 if variant == "kpar_refined" else 0.35
        mp.mp.dps = 40
        h = mp.mpf(ref.H)
        E0 = mp.mpf(float(E0_DEC))
        T = mp.mpf(T_float)
        # Eq. (17) at theta = 0 from the definitions (not the frozen module's helper).
        _, lo, hi = kinematics(T, E0, h)
        u0, u1 = float(mp.log(lo)), float(mp.log(hi))
        n = int(math.ceil((u1 - u0) / width))
        edges = list(np.linspace(u0, u1, n + 1))
        # Bulk lines W_i(q) = T of the dispersive material: poles of F at a distance
        # wt = G_i/(4 alpha_i H q_r^2) from the real ln k axis.  Geometric panels around them.
        for E, A, G, alpha in MATERIALS["Si"]["oscillators"]:
            E, G, alpha = float(E), float(G), float(alpha)
            if alpha > 0 and T_float > E:
                qr2 = (T_float - E) / (alpha * float(ref.H))
                ur, wt = 0.5 * math.log(qr2), G / (4 * alpha * float(ref.H) * qr2)
                step = 0.5 * wt * (0.5 if variant == "kpar_refined" else 1.0)
                offsets = [0.0]
                while offsets[-1] < width:
                    offsets.append(max(2 * offsets[-1], step))
                edges += [ur + sgn * o for o in offsets for sgn in (1, -1) if u0 < ur + sgn * o < u1]
        edges = np.unique(np.asarray(edges))
        xg, wg = leggauss(16)
        u = (0.5 * (edges[1:] - edges[:-1])[:, None] * xg + 0.5 * (edges[1:] + edges[:-1])[:, None]).ravel()
        wu = (0.5 * (edges[1:] - edges[:-1])[:, None] * wg).ravel()
        qs = np.exp(u)
        ws = wu * qs
        mats = {name: Mat(name) for name in ("Si", "S0")}
        acc = {name: {t: 0.0 + 0.0j for t in TERMS + ("M1",)} for name in mats}
        for q, wq in zip(qs, ws):
            I = kp.integrals(q)
            f = wq * 2 * np.pi * q
            for name, mat in mats.items():
                re, im = response(mat, mp.mpf(float(q)) ** 2, T, h)
                S = complex(re, im)
                inv1 = -S                       # 1/eps - 1 = h - 1 with h = -F
                R = S / (2.0 - S)
                C1 = 2.0 * S - R                # (eps-1)(eps+2)/(eps(eps+1))
                C2 = -R                         # (1-eps)/(1+eps)
                d = acc[name]
                d["M1"] += f * inv1 * I["I_ee"]
                d["term_14"] += f * inv1 * I["I_hh"]
                d["term_15_T1"] += f * I["I_S"] * C1 * I["F_ee"]
                d["term_15_T2"] += f * I["I_S"] * C1 * I["F_hh"]
                d["term_15_T3"] += f * I["I_S"] * C2 * I["F_V"]
                d["term_16_e"] += f * I["I_V"] * C2 * (I["F_ee"] + I["F_V"])
                d["term_16_h"] += f * I["I_V"] * C2 * I["F_hh"]
                d["term_13a"] += wq * q * inv1 / (kp.Q ** 2 + q ** 2)
        pref = 2.0 / ((2 * np.pi) ** 4 * a_float * T_float)
        pre13a = 2.0 * ref.E2 / (np.pi * kp.V ** 2)
        out = {}
        for name in mats:
            d = acc[name]
            r = {}
            m1 = pref * (1j * d["M1"]).real
            r["term_13a"] = pre13a * (1j * d["term_13a"]).real
            r["term_13b"] = m1 - r["term_13a"]
            for t in ("term_14", "term_15_T1", "term_15_T2", "term_15_T3", "term_16_e", "term_16_h"):
                r[t] = pref * (1j * d[t]).real
            r["no_hole"] = sum(r[t] for t in NO_HOLE_TERMS)
            r["hole"] = sum(r[t] for t in HOLE_TERMS)
            out[name] = {k: float(v) for k, v in r.items()}
        return {"loss": T_float, "depth": a_float, "variant": variant, "values": out,
                "k_parallel_nodes": int(qs.size)}
    finally:
        ref._XG, ref._WG, ref._XL, ref._WL, ref._XS, ref._WS = saved


# ----------------------------------------------------------------------------- V15

PI_LOWER = Fraction(314159265358979, 10 ** 14)      # pi = 3.14159265358979323... > this


def log_upper(x, terms=60):
    """Rational upper bound of log(x), x >= 1 a Fraction.
    x = 2^n r, 1 <= r < 2; log r = 2 atanh(y), y = (r-1)/(r+1) in [0, 1/3);
    atanh(y) <= sum_{k<N} y^(2k+1)/(2k+1) + y^(2N+1)/((2N+1)(1-y^2))."""
    def atanh_upper(y):
        s = sum(y ** (2 * k + 1) / (2 * k + 1) for k in range(terms))
        return s + y ** (2 * terms + 1) / ((2 * terms + 1) * (1 - y * y))
    n = 0
    r = Fraction(x)
    while r >= 2:
        r /= 2
        n += 1
    log2_upper = 2 * atanh_upper(Fraction(1, 3))
    return n * log2_upper + 2 * atanh_upper((r - 1) / (r + 1))


def delta_g_binary64(g):
    """D7 of docs/DIELECTRIC_DESIGN.md in binary64 arithmetic."""
    return max(math.ulp(0.0), min(g / 2, max(32 * math.ulp(g), 1e-12 * g)))


def si_uniform_bound(binary64):
    """Upper bound of S(q,T) for all q >= 0 and 0 < T <= g - delta_g with U = 3 eV:
    S <= P0 + c log(2g/(g-T)), P0 = B0/(1-(g/U)^2), c = M g/pi,
    B0 = sum A_i/E_i^2, M = sum A_i G_i/(E_i^2-U^2)^2 (design section 3.2).
    Exact rationals; pi from below and log from above."""
    conv = (lambda s: Fraction(float(s))) if binary64 else (lambda s: Fraction(s))
    spec = MATERIALS["Si"]
    g = conv(spec["gap_eV"])
    osc = [tuple(conv(v) for v in o) for o in spec["oscillators"]]
    U = Fraction(3)
    if binary64:
        delta = Fraction(delta_g_binary64(float(spec["gap_eV"])))
    else:
        delta = max(32 * Fraction(math.ulp(1.12)), Fraction(1, 10 ** 12) * g)
    assert g < U < min(o[0] for o in osc)
    B0 = sum(A / (E * E) for E, A, G, _ in osc)
    P0 = B0 / (1 - (g / U) ** 2)
    M = sum(A * G / (E * E - U * U) ** 2 for E, A, G, _ in osc)
    c_upper = M * g / PI_LOWER
    bound = P0 + c_upper * log_upper(2 * g / delta)
    return {"inputs": "binary64" if binary64 else "decimal", "U_eV": 3,
            "delta_g_eV": float(delta), "B0": float(B0), "B0_minus_1": float(B0 - 1),
            "P0": float(P0), "M_per_eV": float(M), "c_upper": float(c_upper),
            "S_upper_bound": float(bound), "F_upper_bound": float(bound - 1),
            "S_bound_below_1.31827": bound < Fraction("1.31827"),
            "static_F_bound_below_-0.0293991": B0 - 1 < Fraction("-0.0293991"),
            "exact": {"B0_minus_1": str(B0 - 1) if not binary64 else None}}


def v15_samples():
    """Sampled check of the bound against the multiprecision F (not part of the proof):
    max over q and T of S(q,T) - bound(T), and of F(q,0) - (B0-1) (must be <= 0)."""
    mp.mp.dps = 40
    mat = Mat("Si")
    h = mp.mpf(b64(H_DEC))
    g = mat.g
    B0 = sum(A / (E * E) for E, A, G, _ in mat.osc)
    P0 = B0 / (1 - (g / 3) ** 2)
    c = sum(A * G / (E * E - 9) ** 2 for E, A, G, _ in mat.osc) * g / mp.pi
    delta = mp.mpf(delta_g_binary64(float(g)))
    worst_s, worst_static, max_s = mp.mpf("-inf"), mp.mpf("-inf"), mp.mpf("-inf")
    for q in ("0", "0.001", "0.05", "0.3", "1", "3", "10", "32"):
        q2 = mp.mpf(q) ** 2
        for d in (delta, mp.mpf("1e-9"), mp.mpf("1e-6"), mp.mpf("1e-3"), mp.mpf("0.12"),
                  mp.mpf("0.72"), mp.mpf("0.92"), mp.mpf("1.1199")):
            T = g - d
            s = response(mat, q2, T, h)[0]
            max_s = max(max_s, s)
            worst_s = max(worst_s, s - (P0 + c * mp.log(2 * g / d)))
        # static value: Int_g^inf L0(u)/u du by quadrature
        val = mp.quad(lambda u: sum(A * G / ((((E + al * h * q2) ** 2 - u * u) ** 2) + G * G * u * u)
                                    for E, A, G, al in mat.osc), [g, 5, 10, 14, 17, 30, 100, mp.inf])
        worst_static = max(worst_static, 2 * val / mp.pi - 1 - (B0 - 1))
    return {"max_S_minus_bound": float(worst_s), "max_static_F_minus_B0m1": float(worst_static),
            "max_sampled_S": float(max_s),
            "note": "sampled with the multiprecision KK of this file; both maxima must be negative"}


# ----------------------------------------------------------------------------- assembly

def kk_cross_check():
    """kk_line against the direct principal-value quadrature of E4a at a few points."""
    mp.mp.dps = 40
    h = mp.mpf(b64(H_DEC))
    rows = []
    cases = [("Si", "0", "0.2"), ("Si", "0.3", "1.119999"), ("Si", "0.3", "1.120001"),
             ("Si", "1", "5"), ("Si", "2.5", "17"), ("Si", "30", "99.8"),
             ("S0", "0", "11"), ("N1", "0.7245", "4"), ("N1", "0.05", "4")]
    for name, q, T in cases:
        mat = Mat(name)
        q2 = mp.mpf(q) ** 2
        Tm = mp.mpf(b64(T))
        closed = response(mat, q2, Tm, h)[0]
        direct = kk_direct(mat, q2, Tm, h)
        rows.append({"material": name, "q_inv_A": q, "loss_eV": T, "re_S_partial_fractions": mp.nstr(closed, 25),
                     "abs_difference_direct_quadrature": float(abs(closed - direct))})
    return rows


def budget(value):
    return max(EPSABS, EPSREL * abs(value))


def build(processes):
    started = time.perf_counter()
    raw = (ROOT / SI_INPUT).read_bytes()
    if hashlib.sha256(raw).hexdigest() != SI_INPUT_SHA256:
        raise SystemExit("Si input record does not match the pre-registered SHA-256")
    record = json.loads(raw.decode("utf-8"))["material"]
    want = MATERIALS["Si"]
    got = [(o["energy_eV"], o["strength_eV2"], o["width_eV"], o["dispersion"]) for o in record["oscillators"]]
    if record["gap_eV"] != float(want["gap_eV"]) or got != [tuple(float(v) for v in o) for o in want["oscillators"]]:
        raise SystemExit("Si input record differs from the decimal inputs of this script")

    depths = [b64(d) for d in DEPTHS_DEC]
    tasks = [(m, s, b64(dec), depths, False) for m in ("Si", "S0") for s in SPECS for _, dec in LOSSES]
    tasks += [("N1", s, b64(N1_LOSS[1]), [b64(N1_DEPTH)], True) for s in SPECS]
    i3_tasks = [(b64(dec), a, v) for _, dec in LOSSES for a in depths
                for v in ("base", "kpar_refined", "kperp_refined")]
    with Pool(processes) as pool:
        closed_async = pool.map_async(closed_form_task, tasks, chunksize=1)
        i3_async = pool.map_async(i3_task, i3_tasks, chunksize=1)
        closed = closed_async.get()
        t_closed = time.perf_counter() - started
        i3 = i3_async.get()
    t_all = time.perf_counter() - started

    labels = {b64(dec): (label, dec) for label, dec in LOSSES + [N1_LOSS]}
    i3_index = {(r["loss"], r["depth"], r["variant"]): r for r in i3}
    points = []
    mp.mp.dps = DPS_B
    for res in closed:
        gap = b64(MATERIALS[res["material"]]["gap_eV"])
        below = res["loss"] < gap
        for i, a in enumerate(res["depths"]):
            sl = slice(8 * i, 8 * i + 8)
            val = dict(zip(TERMS, (mp.mpf(v) for v in res["values"][sl])))
            unc = dict(zip(TERMS, (2 * mp.mpf(d) + mp.mpf(r) for d, r in
                                   zip(res["diff_ab"][sl], res["residual_a"][sl]))))
            groups = {"no_hole": NO_HOLE_TERMS, "hole": HOLE_TERMS,
                      "electron_group": ELECTRON_GROUP, "hole_group": HOLE_TERMS}
            for name, members in groups.items():
                val[name] = mp.fsum(val[t] for t in members)
                unc[name] = mp.fsum(unc[t] for t in members)
            for name, charge in HOLE_CHARGES:
                val[name] = val["no_hole"] + mp.mpf(charge) * val["hole"]
                unc[name] = unc["no_hole"] + mp.mpf(charge) * unc["hole"]
            residual = {}
            zero = []
            if below:
                # Im S = 0: 13a vanishes; electron and hole groups cancel pointwise (design 3.3).
                for name in ("term_13a", "total_q1"):
                    residual[name] = float(abs(val[name]))
                    val[name], unc[name] = mp.mpf(0), mp.mpf(0)
                    zero.append(name)
                val["no_hole"] = val["electron_group"]
                unc["no_hole"] = unc["electron_group"]
            outputs = {}
            for name in val:
                v = float(val[name])
                u = float(unc[name])
                if u != 0:
                    u = math.nextafter(u, math.inf)
                # binary64 rounding of the stored value is part of the reference uncertainty
                u_total = u + (abs(float(val[name] - mp.mpf(v))) if name not in zero else 0.0)
                outputs[name] = {"value": v, "value_hex": v.hex(), "value_decimal": mp.nstr(val[name], 30),
                                 "u": u_total, "B": budget(v), "u_over_B": u_total / budget(v)}
            point = {"material": res["material"], "epsilon_momentum": res["spec"], "depth_A": a,
                     "loss_label": labels[res["loss"]][0], "loss_eV": number(labels[res["loss"]][1]),
                     "region": "below_gap" if below else "above_gap",
                     "identically_zero": zero, "numerical_residual_of_zero": residual,
                     "outputs": outputs,
                     "closed_form_run": {"panels_run_A": res["panels_a"], "nodes_run_A": res["evaluations_a"],
                                         "nodes_run_B": res["evaluations_b"],
                                         "k_parallel_limits_inv_A": res["k_limits"],
                                         "omega_perp_inv_A": res["omega_perp"]}}
            if res["brute_diff"] is not None:
                point["closed_form_run"]["uniform_partition_abs_difference"] = dict(
                    zip(TERMS, (float(v) for v in res["brute_diff"][sl])))
            if res["spec"] == "k_parallel" and res["material"] in ("Si", "S0"):
                base = i3_index[res["loss"], a, "base"]["values"][res["material"]]
                change = {}
                for name in base:
                    change[name] = max(abs(base[name] - i3_index[res["loss"], a, v]["values"][res["material"]][name])
                                       for v in ("kpar_refined", "kperp_refined"))
                point["i3"] = {"values": base, "refinement_change": change,
                               "k_parallel_nodes": i3_index[res["loss"], a, "base"]["k_parallel_nodes"]}
            points.append(point)

    # I3 scale (V3): max |total| at unit hole charge over all reference losses, per material and depth.
    scales = {}
    for p in points:
        if "i3" in p:
            key = f'{p["material"]}|{p["depth_A"]:g}'
            total = abs(p["i3"]["values"]["no_hole"] + p["i3"]["values"]["hole"])
            scales[key] = max(scales.get(key, 0.0), total)
    for p in points:
        if "i3" in p:
            scale = scales[f'{p["material"]}|{p["depth_A"]:g}']
            p["i3"]["v3_scale"] = scale
            worst = 0.0
            for name, v in p["i3"]["values"].items():
                width = 1e-6 * abs(v) + 1e-6 * scale
                worst = max(worst, p["i3"]["refinement_change"][name] / width)
            for name, charge in HOLE_CHARGES:
                v = p["i3"]["values"]["no_hole"] + charge * p["i3"]["values"]["hole"]
                dv = p["i3"]["refinement_change"]["no_hole"] + charge * p["i3"]["refinement_change"]["hole"]
                worst = max(worst, dv / (1e-6 * abs(v) + 1e-6 * scale))
            p["i3"]["refinement_change_over_v3_width_max"] = worst
            # Agreement of the two independent routes of this file (information only).
            agree = 0.0
            for name, v in p["i3"]["values"].items():
                agree = max(agree, abs(v - p["outputs"][name]["value"]) / (1e-6 * abs(v) + 1e-6 * scale))
            p["i3"]["closed_form_difference_over_v3_width_max"] = agree

    # V13: dispersive points where the two conventions differ by >= 1e3 B_j.
    index = {(p["material"], p["epsilon_momentum"], p["depth_A"], p["loss_label"]): p for p in points}
    v13 = []
    for p in points:
        if p["material"] != "Si" or p["epsilon_momentum"] != "k_parallel":
            continue
        other = index["Si", "k_total", p["depth_A"], p["loss_label"]]
        ratios = {}
        for name in p["outputs"]:
            a, b = p["outputs"][name], other["outputs"][name]
            ratios[name] = abs(a["value"] - b["value"]) / max(a["B"], b["B"])
        separated = sorted(n for n, r in ratios.items() if r >= 1e3)
        if all(n in separated for n in ("no_hole", "hole")) and (
                p["region"] == "below_gap" or "total_q1" in separated):
            v13.append({"depth_A": p["depth_A"], "loss_label": p["loss_label"],
                        "separated_outputs": separated,
                        "difference_over_B": {n: ratios[n] for n in separated}})

    max_ratio = max(o["u_over_B"] for p in points for o in p["outputs"].values())
    worst = max(((o["u_over_B"], p["material"], p["epsilon_momentum"], p["depth_A"], p["loss_label"], n)
                 for p in points for n, o in p["outputs"].items()))
    doc = {
        "schema": "xpsloss.m22.references.v1",
        "generated_by": "benchmarks/numerics/m22/generate_references_m22.py",
        "prereg": "docs/M2_2_PREREG.md section 4",
        "units": {"K": "eV^-1 A^-1 per solid path", "loss": "eV", "depth": "A", "momentum": "A^-1"},
        "constants": {"H_eV_A2": number(H_DEC), "a0_A": number(A0_DEC)},
        "conditions": {"electron_energy_eV": number(E0_DEC), "theta_deg": 0,
                       "epsabs": EPSABS, "epsrel": EPSREL, "limit": LIMIT,
                       "budget": "B_j = max(epsabs, epsrel*|K_j,ref|)",
                       "hole_charges": [c for _, c in HOLE_CHARGES]},
        "inputs": {"si_record": SI_INPUT, "si_record_sha256": SI_INPUT_SHA256},
        "materials": {name: {"gap_eV": number(spec["gap_eV"]),
                             "oscillators": [dict(zip(("energy_eV", "strength_eV2", "width_eV", "dispersion"),
                                                      (number(v) for v in o))) for o in spec["oscillators"]]}
                      for name, spec in MATERIALS.items()},
        "outputs": {"terms": list(TERMS), "no_hole": list(NO_HOLE_TERMS), "hole": list(HOLE_TERMS),
                    "electron_group": list(ELECTRON_GROUP), "hole_group": list(HOLE_TERMS),
                    "total_qX": "no_hole + X*hole (hole terms are for unit hole charge)"},
        "method": {
            "closed_form": {
                "library": "mpmath", "digits_run_A": DPS_A, "digits_run_B": DPS_B,
                "response": "Im F = L (E1) with hard gap; Re F by E4a, partial fractions of the retained band",
                "equations": "Simonsen et al. (1997) Eqs. (13)-(17), theta = 0, coefficients from F",
                "quadrature": "adaptive composite 24-point Gauss-Legendre in ln k_parallel; run B halves every panel of run A",
                "panel_tolerance": float(PANEL_TOL),
                "value": "run B", "u": "2*|run A - run B| + accepted residual of run A + binary64 rounding of the stored value",
                "zero": "below the gap term_13a and total_q1 are identically zero (value 0, u 0); numerical residuals are recorded"},
            "i3": {"library": "numpy (binary64), frozen tests/xps_route_b_reference.py KPerp.integrals",
                   "coefficients": "from h = -F of the KK code of this file (40 digits), no reciprocal of 1/eps",
                   "k_parallel": "composite 16-point Gauss-Legendre in ln k_parallel, panel width <= 0.35, geometric panels around the bulk lines W_i(q) = T",
                   "refinements": ["k_parallel panel width <= 0.175 and half the innermost panel at the bulk lines",
                                   "k_perp: core panels halved, orders core 20->28, Laguerre tails 80->120, 1/k tail 32->48"],
                   "scope": "k_parallel, Si and S0"}},
        "kk_cross_check": kk_cross_check(),
        "v13_points": v13,
        "v15": {"decimal": si_uniform_bound(False), "binary64": si_uniform_bound(True),
                "samples": v15_samples()},
        "summary": {"points": len(points), "max_u_over_B": max_ratio,
                    "max_u_over_B_at": list(worst[1:]),
                    "max_i3_refinement_change_over_v3_width": max(
                        p["i3"]["refinement_change_over_v3_width_max"] for p in points if "i3" in p),
                    "max_i3_closed_form_difference_over_v3_width": max(
                        p["i3"]["closed_form_difference_over_v3_width_max"] for p in points if "i3" in p),
                    "max_abs_numerical_residual_of_zero": max(
                        (v for p in points for v in p["numerical_residual_of_zero"].values()), default=0.0)},
        "environment": {"platform": platform.platform(), "python": sys.version.split()[0],
                        "mpmath": mp.__version__, "numpy": __import__("numpy").__version__,
                        "processes": processes},
        "timing_seconds": {"closed_form_wall": round(t_closed, 1), "all_wall": round(t_all, 1),
                           "closed_form_cpu_sum": round(sum(r["seconds"] for r in closed), 1)},
        "points": points,
    }
    return doc


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--processes", type=int, default=max(1, min(12, (os.cpu_count() or 2) - 2)))
    parser.add_argument("--output", default=str(HERE / "references_m22.json"))
    args = parser.parse_args()
    doc = build(args.processes)
    if doc["summary"]["max_u_over_B"] > 0.01:
        print("WARNING: u_j > 0.01 B_j at", doc["summary"]["max_u_over_B_at"], file=sys.stderr)
    text = json.dumps(doc, indent=1, ensure_ascii=True) + "\n"
    Path(args.output).write_bytes(text.encode("ascii"))
    print("wrote", args.output)
    print("sha256", hashlib.sha256(text.encode("ascii")).hexdigest())
    print("summary", json.dumps(doc["summary"]))
    print("timing", json.dumps(doc["timing_seconds"]))


if __name__ == "__main__":
    main()
