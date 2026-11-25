"""Independent multiprecision references for M3a (docs/M3A_PREREG.md, section 5, R1).

Written by the reference author R1, separately from the implementation.  Nothing from
the package xpsloss is imported.  Only mpmath and the standard library are used.

Usage (the results of the 164 groups (material, energy, convention, loss) are cached in
DIR, one file per group, so the run can be interrupted and resumed; several instances may
share DIR):

    python benchmarks/numerics/m3a/generate_references_m3a.py --cache DIR --stage compute
    python benchmarks/numerics/m3a/generate_references_m3a.py --cache DIR --stage assemble
    python benchmarks/numerics/m3a/generate_references_m3a.py --cache DIR --time-one 35

The frozen file was computed with 12 processes (two instances, 10 and 2) on a shared
12-core laptop in 59 minutes of wall time, 39364 s of CPU (timing_seconds in the JSON).  The JSON contains timings, so a
new run gives the same values but another SHA-256.

Physics (all typed here from docs/M1_EQUATION_MAP.md section 4.1, the transcription of
Simonsen, Yubero, Tougaard, PRB 56, 1612 (1997), journal pp. 1614-1615, PDF pp. 3-4):

  * Im F = L (E1, docs/EQUATIONS.md) with the hard gap, F = -1/eps, S = 1 + F.
    g = 0: S = sum_i A_i/(W_i^2 - T^2 - i G_i T) (E4b), W_i = E_i + alpha_i H q^2.
    g > 0: Re S = (2/pi) PV Int_g^inf u L0(q,u)/(u^2-T^2) du (E4a) by partial fractions of
    the retained band (kk_line); the two forms are compared at g = 0 (kk_cross_check).
  * Eqs. (13)-(16) at theta = 0 (v_par = 0, Omega0 = omega_perp = b, E_perp = E0, x = a),
    hbar w > 0 so that pi delta(omega_perp) + i/omega_perp -> i/b.  Eq. (17) gives the
    k_parallel limits.  R = S/(2-S), C = 2S - R.
  * Eq. (18) (journal p. 1615): K_sc = Int_0^inf (x/lambda^2) exp(-x/lambda) K(a=x) dx.

Every term except 13a is  (1/x) P2a k^2 Re[Z_j piece_j(k, a)]  per d(ln k), with
P2a = 1/(2 pi E0 a0 b), a dielectric coefficient Z_j(k) and a printed depth factor
piece_j.  Hence K_sc,j = (1/lambda^2) Int d(ln k) P2a k^2 Re[Z_j Int_0^inf e^{-pa} piece_j da],
p = 1/lambda.  The depth integral is done in three independent ways:

  (a)  closed form derived here by hand (compact_laplace), value of record;
  (b)  numerical depth quadrature of the printed depth factors, composite Gauss-Legendre
       on panels cut by the decay rates p+2k, p+k, p and by the phase b*a.  It does not use
       the weights of (a).  For the cost (about 1.5e8 depth nodes) the quadrature is done
       for the elementary depth factors of the printed braces, 1, e^{-iba-ka}, e^{-2ka} and
       (e^{-iba}-e^{-ka})(e^{iba}-e^{-ka}), evaluated with numbers at every depth node
       (fixed-point integers, 256 fractional bits), and the depth-independent factors of
       the braces are applied afterwards (pieces_from_factors).  This regrouping is checked
       pointwise against printed_pieces, the fixed-point arithmetic against mpf, and for
       the cases BLACKBOX the whole printed pieces are integrated node by node
       (formula_self_check, blackbox_task);
  (g)  the electron and hole groups are assembled from the full printed braces with a
       small symbolic "mode" algebra (products of e^{+-i b a} and e^{-k a}) and
       transformed mode by mode, 1/(p + n k - i m b).  This is neither (a) nor a sum of
       the eight terms.

k_parallel integration of (a): adaptive composite 24-point Gauss-Legendre in ln k.
Run A: 50 digits, adaptive.  Run B: 80 digits, every panel of run A halved.  Value = run B.
Route (b) uses the accepted parent panels of run A (24 points each) at 50 digits.
The dielectric response (KK) is evaluated once per k node and shared by all lengths; the
hole charge enters only in the test (linear combination).  Narrow lines and the surface
resonance: panel edges at W_i(q) = T and at the maxima of |1/(2-S)| on a 400-point scan,
then bisection.  Above the gap the principal value is taken in closed form (kk_line).

Uncertainty of an averaged term j (stored with its parts for every key and column):

    u_j = |a - b| + 2 |run A - run B| + residual(run A)
          + max(|order change|, |phase-split change|, |tail change|) of route (b)
          + rounding of the stored 40-digit decimal
    u(term_13a) = 2 |A - B| + residual + rounding               (no depth dependence)
    u(group)    = |group(g) - sum of its terms of route (b)| + 2 |A - B| + residual
                  + sum over members of the route (b) changes + rounding
    u(no_hole)  = u(term_13a) + u(electron_group),  u(hole) = u(hole_group)
    total(Q_h)  = K_e + Q_h K_h with u = u_e + |Q_h| u_h (formed by the test in rationals)

Below a hard gap (T < g) Im S = 0, so term_13a = 0 and the unit-charge total = 0
analytically (value "0", u "0"); the numerical residuals of the raw sums are stored.
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
from mpmath.libmp import to_fixed
from mpmath.libmp.libelefun import cos_sin_fixed, exp_fixed, ln2_fixed, pi_fixed

HERE = Path(__file__).resolve().parent
TABLE = HERE / "expected_states_m3a.json"
TABLE_SHA256 = "d87dad8f02a2fd211679229bd7c5205d908dbb2b3b241b30a27e94469f414ed5"

H_DEC = "3.80998212"            # hbar^2/2m [eV A^2], docs/EQUATIONS.md E1
A0_DEC = "0.529177210903"       # Bohr radius [A], docs/EQUATIONS.md E2

TERMS = ("term_13a", "term_13b", "term_14", "term_15_T1", "term_15_T2",
         "term_15_T3", "term_16_e", "term_16_h")
DEPTH_TERMS = TERMS[1:]
ELECTRON_IDX = (0, 2, 4, 5)       # 13b, 15-T1, 15-T3, 16-e in DEPTH_TERMS
HOLE_IDX = (1, 3, 6)              # 14, 15-T2, 16-h
OUTPUTS = TERMS + ("no_hole", "hole", "total")
LIMIT_CHARGES = (0.0, 0.5, 1.0)

DPS_A, DPS_B = 50, 80
GL_DEGREE = 4                    # mpmath degree 4 = 24 nodes, degree 5 = 48 nodes
REL_TOL = mp.mpf("1e-24")        # per panel and column, relative to the coarse integral
ENV_REL_TOL = mp.mpf("1e-13")    # envelopes (|.| has kinks); they are used as upper bounds
FLOOR_TOL = mp.mpf("1e-30")      # times the largest column
MAX_LEVEL = 60
PANEL_WIDTH = "0.5"

# Route (b): depth panels.  Decay rates of exp(-p a) times {1, u, |u|^2}: p, p + k, p + 2k.
TAIL = 52                        # modes are dropped where exp(-rate*a) < exp(-TAIL) = 2.6e-23
TAIL_EXT = 64
# (panel size c, GL degree).  n-point Gauss-Legendre on exp(z x), |z| * length/2 = c/2:
# error about 1e-27 for (20, 24 points), 1e-22 for (26, 24 points), 1e-54 for (40, 48 points).
B_VARIANTS = {"base": (20, 4), "order": (40, 5), "phase": (26, 4)}

FIXED_BITS = 256                 # fixed-point depth quadrature of route (b): absolute 2^-256
VALUE_DIGITS = 40
J = mp.mpc(0, 1)


def b64(dec):
    return float(dec)


def fhex(x):
    return float(x).hex()


# ----------------------------------------------------------------------------- response

class Mat:
    """Material with the binary64 inputs as exact mpf numbers."""

    def __init__(self, gap, oscillators):
        self.g = mp.mpf(gap)
        self.osc = [tuple(mp.mpf(v) for v in o) for o in oscillators]
        self.dispersive = any(o[3] != 0 and o[1] != 0 for o in self.osc)
        self._cache = {}


def kk_line(E, A, G, alpha, g, q2, T, h):
    """Int_g^inf u L0_i(q,u)/(u^2-T^2) du (principal value), T > 0, T != g.

    f(u) = A u/((W^2-u^2-iGu)(u^2-T^2)) = -A u/((u-u1)(u-u2)(u-T)(u+T)), Im f = u L0/(u^2-T^2)
    on the real axis.  With f = sum_p c_p/(u-p) and sum_p c_p = 0,
    PV Int_g^inf f du = -sum_p c_p log(g-p): principal log for u1, u2, log|g-p| for +-T.
    """
    W = E + alpha * h * q2
    root = mp.sqrt(4 * W * W - G * G)
    u1, u2 = (-J * G + root) / 2, (-J * G - root) / 2
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


def response(mat, q2, T, h, force_kk=False):
    """S = 1 + F (complex, or real below the gap)."""
    key = None
    if not mat.dispersive:
        key = (mp.mp.prec, T, force_kk)
        if key in mat._cache:
            return mat._cache[key]
    if mat.g == 0 and not force_kk:
        out = mp.mpc(0)
        for E, A, G, alpha in mat.osc:
            if A == 0:
                continue
            W = E + alpha * h * q2
            out += A / (W * W - T * T - J * G * T)
    else:
        re = mp.mpf(0)
        im = mp.mpf(0)
        for E, A, G, alpha in mat.osc:
            if A == 0:
                continue
            re += kk_line(E, A, G, alpha, mat.g, q2, T, h)
            if T > mat.g:
                W = E + alpha * h * q2
                im += A * G * T / ((W * W - T * T) ** 2 + G * G * T * T)
        out = mp.mpc(2 * re / mp.pi, im) if im != 0 else 2 * re / mp.pi
    if key is not None:
        mat._cache[key] = out
    return out


def loss_function(mat, q2, T, h):
    if T <= mat.g:
        return mp.mpf(0)
    s = mp.mpf(0)
    for E, A, G, alpha in mat.osc:
        W = E + alpha * h * q2
        s += A * G * T / ((W * W - T * T) ** 2 + G * G * T * T)
    return s


def kinematics(T, E0, h):
    """b = T/(hbar v), hbar v = 2 sqrt(H E0); Eq. (17) at theta = 0:
    k_par+- = sqrt(k+-^2 - b^2), k+- = (sqrt(E0) +- sqrt(E0 - T))/sqrt(H)."""
    b = T / (2 * mp.sqrt(h * E0))
    s, r = mp.sqrt(E0), mp.sqrt(E0 - T)
    qlo, qhi = (s - r) / mp.sqrt(h), (s + r) / mp.sqrt(h)
    return b, mp.sqrt(qlo * qlo - b * b), mp.sqrt(qhi * qhi - b * b), qlo, qhi


# ----------------------------------------------------------------------------- printed forms

class Modes:
    """Linear combinations of E^m d^n with E = exp(i b a), d = exp(-k a): {(m, n): coefficient}."""

    __slots__ = ("c",)

    def __init__(self, c):
        self.c = c

    @staticmethod
    def lift(x):
        return x if isinstance(x, Modes) else Modes({(0, 0): x})

    def __add__(self, other):
        other = Modes.lift(other)
        out = dict(self.c)
        for key, v in other.c.items():
            out[key] = out.get(key, 0) + v
        return Modes(out)

    __radd__ = __add__

    def __neg__(self):
        return Modes({key: -v for key, v in self.c.items()})

    def __sub__(self, other):
        return self + (-Modes.lift(other))

    def __rsub__(self, other):
        return Modes.lift(other) + (-self)

    def __mul__(self, other):
        if not isinstance(other, Modes):
            return Modes({key: v * other for key, v in self.c.items()})
        out = {}
        for (m1, n1), v1 in self.c.items():
            for (m2, n2), v2 in other.c.items():
                key = (m1 + m2, n1 + n2)
                out[key] = out.get(key, 0) + v1 * v2
        return Modes(out)

    __rmul__ = __mul__

    def __truediv__(self, scalar):
        return Modes({key: v / scalar for key, v in self.c.items()})

    def laplace(self, k, b, p):
        """Int_0^inf exp(-p a) (.) da, mode by mode."""
        return mp.fsum(v / (p + n * k - J * m * b) for (m, n), v in self.c.items())

    def coefficients(self):
        """(c_1, c_u, c_ubar, c_q) for the modes 1, u = E^-1 d, ubar = E d, q = d^2."""
        names = {(0, 0): 0, (-1, 1): 1, (1, 1): 2, (0, 2): 3}
        out = [mp.mpc(0)] * 4
        for key, v in self.c.items():
            if key not in names:
                raise RuntimeError(f"unexpected depth mode {key}")
            out[names[key]] = mp.mpc(v)
        return out


def mode_symbols():
    return Modes({(1, 0): mp.mpf(1)}), Modes({(-1, 0): mp.mpf(1)}), Modes({(0, 1): mp.mpf(1)})


def printed_pieces(k, b, E, Ei, d, one=1):
    """The depth factors of the seven depth-dependent terms as printed (theta = 0).

    E = exp(+i b a), Ei = exp(-i b a), d = exp(-k a); numbers or Modes.  Term j is
    (1/a) P2a Re[Z_j piece_j] per k dk with Z = coefficients(S, n)."""
    n = b * b + k * k
    pp, mm = mp.mpc(b, k), mp.mpc(b, -k)         # Omega0 + i k, Omega0 - i k
    hf = J / b                                   # pi delta(omega_perp) + i/omega_perp
    # Eq. (13), second term, the three pieces of the brace
    b13 = J * pp * pp * (Ei * d - one) + (-J) * mm * mm * (E * d - one) + J * ((Ei - d) * (E - d)) * n
    # Eq. (14)
    b14 = hf * (J * pp * (one - d * Ei) + (-J) * mm * (d * d - d * Ei))
    # Eq. (15): common factor (-i)(O - i k)(e^{-ka} - e^{-iOa})/(O^2+k^2)
    X = ((-J) * mm / n) * (d - Ei)
    b15_1 = J * (X * (E - d)) * (pp / n)
    b15_2 = J * (X * d) * ((-J) * hf)
    b15_3 = J * (X * E) * (mm / n)
    # Eq. (16): (O + i k) e^{-iOa}/(O^2+k^2) times the brace
    Y = (pp / n) * Ei
    b16_e = (Y * (2 * J * k * E - pp * d)) / n
    b16_h = (Y * d) * ((-J) * hf)
    return [b13, b14, b15_1, b15_2, b15_3, b16_e, b16_h]


class NumericPieces:
    """printed_pieces for numbers, with the depth-independent factors of each printed
    product computed once per k (same factors, same order of the depth factors)."""

    def __init__(self, k, b):
        n = b * b + k * k
        pp, mm = mp.mpc(b, k), mp.mpc(b, -k)
        hf = J / b
        self.a13, self.b13, self.c13 = J * pp * pp, (-J) * mm * mm, J * n
        self.a14, self.b14 = hf * J * pp, hf * (-J) * mm
        self.x = (-J) * mm / n
        self.c151, self.c152, self.c153 = J * pp / n, J * ((-J) * hf), J * mm / n
        self.y, self.two_jk, self.pp, self.inv_n = pp / n, 2 * J * k, pp, 1 / n
        self.c16h = (-J) * hf

    def at(self, E, Ei, d):
        u, ub = Ei * d, E * d
        e_d, ei_d = E - d, Ei - d
        b13 = self.a13 * (u - 1) + self.b13 * (ub - 1) + self.c13 * (ei_d * e_d)
        b14 = self.a14 * (1 - u) + self.b14 * (d * d - u)
        X = self.x * (d - Ei)
        Y = self.y * Ei
        return [b13, b14, (X * e_d) * self.c151, (X * d) * self.c152, (X * E) * self.c153,
                (Y * (self.two_jk * E - self.pp * d)) * self.inv_n, (Y * d) * self.c16h]


def coefficients(S, n):
    """Z_j with term_j = (1/a) P2a Re[Z_j piece_j]: signs and dielectric factors of Eqs. (13)-(16).
    1/eps - 1 = -S, (eps-1)(eps+2)/(eps(eps+1)) = C, (eps-1)/(eps+1) = R, (1-eps)/(1+eps) = -R."""
    R = S / (2 - S)
    C = 2 * S - R
    return [-(J * (-S)) / (n * n), -(-S) / n, -C, -C, -R, -R, -R]


def printed_groups(k, b, S, E, Ei, d, one=1):
    """a*K/P2a of the full Eqs. (13) second term, (14), (15), (16) with the hole factor hf
    kept symbolic: returns (electron part, hole part) before taking the real part."""
    n = b * b + k * k
    pp, mm = mp.mpc(b, k), mp.mpc(b, -k)
    R = S / (2 - S)
    C = 2 * S - R
    hf = J / b
    brace13 = (J * pp * pp * (Ei * d - one) + (-J) * mm * mm * (E * d - one)
               + J * ((Ei - d) * (E - d)) * n)
    k13 = -(J * (-S) / (n * n)) * brace13
    brace14 = J * pp * (one - d * Ei) + (-J) * mm * (d * d - d * Ei)
    k14 = -(hf * (-S) / n) * brace14
    X = ((-J) * mm / n) * (d - Ei)
    k15_e = -J * (X * (C * (pp / n) * (E - d) + R * (mm / n) * E))
    k15_h = -J * (X * (C * ((-J) * hf) * d))
    Y = ((-R) * pp / n) * Ei
    k16_e = Y * ((2 * J * k * E - pp * d) / n)
    k16_h = Y * (((-J) * hf) * d)
    return k13 + k15_e + k16_e, k14 + k15_h + k16_h


def compact_laplace(k, b, p):
    """Route (a): Int_0^inf exp(-p a) piece_j(a) da in closed form (derived by hand).

    s = k + i b, D = |s|^2, L = p + k, Jw = p + 2k, H = |p + s|^2.
    L[1-u] = s/(p(p+s)), L[u-|u|^2] = conj(s)/((p+s) Jw), L[|1-u|^2] = 2 L D/(p Jw H).
    With O+ik = i conj(s), O-ik = -i s:
      piece_13b -> D (2 b Jw + 2 i L D/Jw)/(p H)        piece_14   -> -2 k D hf/(p Jw (p+s))
      piece_15_1 -> -2 L/(p Jw H)                       piece_15_2 -> hf/(Jw (p+s))
      piece_15_3 -> 1/(conj(s) p (p+conj(s)))           piece_16_e -> -Jw/(D p (p+s))
      piece_16_h -> hf/(s (p+s))
    """
    s = mp.mpc(k, b)
    sc = mp.mpc(k, -b)
    D = k * k + b * b
    L = p + k
    Jw = p + 2 * k
    H = L * L + b * b
    hf = J / b
    ps = p + s
    return [D * (2 * b * Jw + 2 * J * L * D / Jw) / (p * H),
            -2 * k * D * hf / (p * Jw * ps),
            -2 * L / (p * Jw * H),
            hf / (Jw * ps),
            1 / (sc * p * (p + sc)),
            -Jw / (D * p * ps),
            hf / (s * ps)]


# ----------------------------------------------------------------------------- quadrature

_NODES = {}


def gl_nodes(degree):
    key = (degree, mp.mp.prec)
    if key not in _NODES:
        _NODES[key] = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(degree, mp.mp.prec)
    return _NODES[key]


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


def integrate_adaptive(f, edges, rel):
    """Bisect until |I(parent) - I(left) - I(right)| <= tol_j for every column.
    tol_j = rel_j |coarse I_j| + FLOOR_TOL max_j |coarse I_j|.
    Returns the sum over accepted children, the children, the accepted parents, the sum of
    the accepted differences and the tolerances."""
    nodes = gl_nodes(GL_DEGREE)
    stack = [(lo, hi, panel(f, nodes, lo, hi), 0) for lo, hi in zip(edges[:-1], edges[1:])]
    coarse = [mp.fsum(abs(s[2][j]) for s in stack) for j in range(len(stack[0][2]))]
    top = max(coarse)
    tol = [r * c + FLOOR_TOL * top for r, c in zip(rel, coarse)]
    total, leaves, parents, resid = None, [], [], None
    while stack:
        lo, hi, whole, level = stack.pop()
        mid = (lo + hi) / 2
        left, right = panel(f, nodes, lo, mid), panel(f, nodes, mid, hi)
        diff = [abs(w - l - r) for w, l, r in zip(whole, left, right)]
        if all(dv <= t for dv, t in zip(diff, tol)):
            part = [l + r for l, r in zip(left, right)]
            total = part if total is None else [s + c for s, c in zip(total, part)]
            resid = diff if resid is None else [s + c for s, c in zip(resid, diff)]
            leaves.extend(((lo, mid), (mid, hi)))
            parents.append((lo, hi))
        else:
            if level >= MAX_LEVEL:
                raise RuntimeError("adaptive quadrature reached the level limit")
            stack.append((lo, mid, left, level + 1))
            stack.append((mid, hi, right, level + 1))
    return total, sorted(leaves), sorted(parents), resid, tol


def integrate_fixed(f, panels):
    nodes = gl_nodes(GL_DEGREE)
    total = None
    for lo, hi in panels:
        part = panel(f, nodes, lo, hi)
        total = part if total is None else [s + c for s, c in zip(total, part)]
    return total


def halved(leaves, t0, t1):
    fine = []
    for lo, hi in leaves:
        lo, hi = mp.mpf(lo), mp.mpf(hi)
        fine.extend(((lo, (lo + hi) / 2), ((lo + hi) / 2, hi)))
    fine[0] = (t0, fine[0][1])           # end points recomputed at the working precision
    fine[-1] = (fine[-1][0], t1)
    return fine


# ----------------------------------------------------------------------------- group task

class Context:
    """One (material, E0, convention, loss) at the current working precision."""

    def __init__(self, spec):
        self.mat = Mat(spec["gap"], spec["oscillators"])
        self.convention = spec["convention"]
        self.T = mp.mpf(spec["loss"])
        self.E0 = mp.mpf(spec["energy"])
        self.h = mp.mpf(b64(H_DEC))
        self.a0 = mp.mpf(b64(A0_DEC))
        self.b, self.lo, self.hi, self.qlo, self.qhi = kinematics(self.T, self.E0, self.h)
        self.lams = [mp.mpf(v) for v in spec["lengths"]]
        self.ps = [1 / v for v in self.lams]
        self.P2a = 1 / (2 * mp.pi * self.E0 * self.a0 * self.b)
        self.P1 = 1 / (mp.pi * self.a0 * self.E0)
        self.t0, self.t1 = mp.log(self.lo), mp.log(self.hi)
        self.below = self.T < self.mat.g
        self.symbols = mode_symbols()

    def S(self, k):
        q2 = k * k + (self.b * self.b if self.convention == "k_total" else 0)
        return response(self.mat, q2, self.T, self.h)

    def edges(self):
        n0 = int(mp.ceil((self.t1 - self.t0) / mp.mpf(PANEL_WIDTH)))
        pts = [self.t0 + (self.t1 - self.t0) * i / n0 for i in range(n0 + 1)]
        extra = [mp.log(self.b)] + [mp.log(p) for p in self.ps]
        for E, A, G, alpha in self.mat.osc:
            if A != 0 and alpha != 0 and self.T > E:
                k2 = (self.T - E) / (alpha * self.h) - (self.b ** 2 if self.convention == "k_total" else 0)
                if k2 > 0:
                    extra.append(mp.log(k2) / 2)
        if self.mat.dispersive:
            n = 400
            grid = [self.t0 + (self.t1 - self.t0) * i / n for i in range(n + 1)]
            w = [abs(1 / (2 - self.S(mp.exp(t)))) for t in grid]
            for i in range(1, n):
                if w[i] > w[i - 1] and w[i] >= w[i + 1] and w[i] > 2:
                    extra.append(grid[i])
        pts += [t for t in extra if self.t0 < t < self.t1]
        return sorted(set(pts))


class MainIntegrand:
    """Columns: term_13a, then for each length the seven terms of route (a) and the two
    groups of route (g).  Per d(ln k)."""

    def __init__(self, ctx):
        self.c = ctx
        self.evaluations = 0

    def __call__(self, t):
        c = self.c
        self.evaluations += 1
        k = mp.exp(t)
        S = c.S(k)
        n = c.b * c.b + k * k
        k2 = k * k
        out = [c.P1 * k2 * mp.im(S) / n]            # Re[i (1/eps - 1)/n] = Im S / n
        Z = coefficients(S, n)
        E, Ei, d = c.symbols
        ge, gh = printed_groups(k, c.b, S, E, Ei, d, Modes.lift(mp.mpf(1)))
        w = c.P2a * k2
        for p in c.ps:
            wl = w * p * p
            Ls = compact_laplace(k, c.b, p)
            out.extend(wl * mp.re(z * l) for z, l in zip(Z, Ls))
            out.append(wl * mp.re(ge.laplace(k, c.b, p)))
            out.append(wl * mp.re(gh.laplace(k, c.b, p)))
        return out


def depth_regions(k, b, p, size, tail=TAIL):
    """Regions (start, stop, number of equal panels) of the depth integral for one k: the
    mode |u|^2 e^{-pa} lives to A1, u e^{-pa} (phase b a) to A2, e^{-pa} to A3.  Panel
    length <= size/(decay rate + b) of the fastest live mode, so that |rate|*length/2 <= size/2."""
    out = []
    a = mp.mpf(0)
    for stop, rate in ((tail / (p + 2 * k), p + 2 * k + b), (tail / (p + k), p + k + b), (tail / p, p)):
        out.append((a, stop, int(mp.ceil((stop - a) * rate / size))))
        a = stop
    return out


def depth_panels(k, b, p, size, tail=TAIL):
    return [(a + (stop - a) * i / n, a + (stop - a) * (i + 1) / n)
            for a, stop, n in depth_regions(k, b, p, size, tail) for i in range(n)]


_FIXED_NODES = {}


def fixed_nodes(degree):
    """Positive Gauss-Legendre nodes and their weights as fixed-point integers (the rule is
    symmetric and has an even number of points)."""
    key = (degree, mp.mp.prec)
    if key not in _FIXED_NODES:
        nodes = gl_nodes(degree)
        pos = [(to_fixed(x._mpf_, FIXED_BITS), to_fixed(w._mpf_, FIXED_BITS)) for x, w in nodes if x > 0]
        if 2 * len(pos) != len(nodes) or abs(mp.fsum(w for _, w in nodes) - 2) > mp.mpf(2) ** (10 - mp.mp.prec):
            raise RuntimeError("unexpected Gauss-Legendre rule")
        _FIXED_NODES[key] = pos
    return _FIXED_NODES[key]


def depth_factor_integrals(k, b, p, regions, degree):
    """Int exp(-p a) F(a) da over the regions for the elementary depth factors of the printed
    braces, evaluated with numbers at every depth node:
    1, u = e^{-i b a} e^{-k a}, e^{-2 k a} and (e^{-i b a} - e^{-k a})(e^{i b a} - e^{-k a}).

    Fixed-point integers with FIXED_BITS fractional bits.  In a region of equal panels the
    exponentials at the node a = mid +- da are products of the values at mid and at da
    (addition theorems; the values at -da are reciprocals and reflections), and the values at
    the panel centres follow from the first centre by the step h.  Where k*h is large,
    exp(-k a) is taken directly at every node."""
    P = FIXED_BITS
    pi2, ln2 = pi_fixed(P - 1), ln2_fixed(P)
    nodes = fixed_nodes(degree)
    kf, bf, pf = (to_fixed(v._mpf_, P) for v in (k, b, p))
    two = 1 << (2 * P)
    t1 = tur = tui = tq = t3 = 0
    for start, stop, n in regions:
        h = (stop - start) / n
        half_f = to_fixed((h / 2)._mpf_, P)
        h_f = to_fixed(h._mpf_, P)
        mid_f = to_fixed((start + h / 2)._mpf_, P)
        direct = k * h > 80
        offsets = []
        for xf, wf in nodes:
            da = (half_f * xf) >> P
            ep = exp_fixed((pf * da) >> P, P, ln2)              # exp(+p da) >= 1
            cn, sn = cos_sin_fixed((bf * da) >> P, P, pi2)
            if direct:
                d_plus, d_minus = da, -da
            else:
                d_minus = exp_fixed((kf * da) >> P, P, ln2)     # exp(-k (-da))
                d_plus = two // d_minus
            offsets.append(((wf << P) // ep, d_plus, cn, sn))
            offsets.append(((wf * ep) >> P, d_minus, cn, -sn))
        em = exp_fixed(-((pf * mid_f) >> P), P, ln2)
        e_step = exp_fixed(-((pf * h_f) >> P), P, ln2)
        cm, sm = cos_sin_fixed((bf * mid_f) >> P, P, pi2)
        c_step, s_step = cos_sin_fixed((bf * h_f) >> P, P, pi2)
        if direct:
            dm, d_step = mid_f, None
        else:
            dm = exp_fixed(-((kf * mid_f) >> P), P, ln2)
            d_step = exp_fixed(-((kf * h_f) >> P), P, ln2)
        for j in range(n):
            if j:
                em = (em * e_step) >> P
                cm, sm = (cm * c_step - sm * s_step) >> P, (sm * c_step + cm * s_step) >> P
                dm = dm + h_f if direct else (dm * d_step) >> P
            s1 = sur = sui = sq = s3 = 0
            for wen, dn, cn, sn in offsets:
                c = (cm * cn - sm * sn) >> P
                sv = (sm * cn + cm * sn) >> P
                if direct:
                    d = exp_fixed(-((kf * (dm + dn)) >> P), P, ln2)
                else:
                    d = (dm * dn) >> P
                s1 += wen
                sur += wen * ((d * c) >> P)
                sui -= wen * ((d * sv) >> P)
                sq += wen * ((d * d) >> P)
                cd = c - d
                s3 += wen * ((cd * cd + sv * sv) >> P)
            scale = (em * half_f) >> P
            t1 += (s1 * scale) >> P
            tur += ((sur >> P) * scale) >> P
            tui += ((sui >> P) * scale) >> P
            tq += ((sq >> P) * scale) >> P
            t3 += ((s3 >> P) * scale) >> P
    back = lambda v: mp.ldexp(mp.mpf(v), -P)
    return [back(t1), mp.mpc(back(tur), back(tui)), back(tq), back(t3)]


def depth_factor_integrals_mpf(k, b, p, panels, degree):
    """The same integrals as depth_factor_integrals with mpf numbers and a direct evaluation
    of every exponential (check of the fixed-point version)."""
    nodes = gl_nodes(degree)
    i1 = iur = iui = iq = i3 = mp.mpf(0)
    for lo, hi in panels:
        half, mid = (hi - lo) / 2, (hi + lo) / 2
        for x, w in nodes:
            a = mid + half * x
            c, sn = mp.cos_sin(b * a)
            d = mp.exp(-k * a)
            wa = w * half * mp.exp(-p * a)
            wd = wa * d
            i1 += wa
            iur += wd * c
            iui -= wd * sn
            iq += wd * d
            cd = c - d
            i3 += wa * (cd * cd + sn * sn)
    return [i1, mp.mpc(iur, iui), iq, i3]


def pieces_from_factors(k, b, i1, iu, iq, i3):
    """The printed pieces as linear combinations of the elementary depth factors.  With
    pointwise values (1, u, |u|^2, |e^{-iba} - e^{-ka}|^2) this is printed_pieces at one depth
    (formula_self_check); with their depth integrals it is the depth integral of the pieces,
    because the other factors do not depend on the depth."""
    n = b * b + k * k
    pp, mm = mp.mpc(b, k), mp.mpc(b, -k)
    hf = J / b
    f1 = iu - i1                    # e^{-i O a - k a} - 1
    f1c = mp.conj(f1)               # e^{+i O a - k a} - 1
    f5 = iq - iu                    # e^{-2 k a} - e^{-k a - i O a}
    b13 = J * pp * pp * f1 + (-J) * mm * mm * f1c + J * i3 * n
    b14 = hf * (J * pp * (-f1) + (-J) * mm * f5)
    xc = (-J) * mm / n
    b15_1 = J * (xc * (-i3)) * (pp / n)          # (d - Ei)(E - d) = -(Ei - d)(E - d)
    b15_2 = J * (xc * f5) * ((-J) * hf)          # (d - Ei) d = d d - d Ei
    b15_3 = J * (xc * f1c) * (mm / n)            # (d - Ei) E = d E - 1
    yc = pp / n
    b16_e = (yc * (2 * J * k * i1 - pp * iu)) / n    # Ei (2 i k E - (O + i k) d)
    b16_h = (yc * iu) * ((-J) * hf)
    return [b13, b14, b15_1, b15_2, b15_3, b16_e, b16_h]


def depth_quadrature(k, b, p, panels, degree):
    """Int exp(-p a) piece_j(a) da over the panels: the whole printed pieces evaluated with
    numbers at every depth node (used by the black-box check of route (b))."""
    nodes = gl_nodes(degree)
    acc = [mp.mpc(0)] * 7
    pieces_at = NumericPieces(k, b).at
    for lo, hi in panels:
        half, mid = (hi - lo) / 2, (hi + lo) / 2
        for x, w in nodes:
            a = mid + half * x
            E = mp.expj(b * a)
            pieces = pieces_at(E, mp.conj(E), mp.exp(-k * a))
            wa = w * half * mp.exp(-p * a)
            acc = [s + wa * v for s, v in zip(acc, pieces)]
    return acc


def route_b(ctx, parents, variants=("base", "order", "phase", "tail"), blackbox=False):
    """Route (b) on the parent panels of run A: {variant: [per length: seven terms]}.
    The depth quadrature is done for the elementary depth factors and the depth-independent
    factors of the printed braces are applied afterwards (pieces_from_factors).  With
    blackbox=True the whole printed pieces are evaluated at every depth node instead."""
    nodes = gl_nodes(GL_DEGREE)
    nl = len(ctx.ps)
    acc = {v: [[mp.mpf(0)] * 7 for _ in range(nl)] for v in variants}
    depth_nodes = 0
    for lo, hi in parents:
        half, mid = (hi - lo) / 2, (hi + lo) / 2
        for x, w in nodes:
            k = mp.exp(mid + half * x)
            Z = coefficients(ctx.S(k), ctx.b * ctx.b + k * k)
            wk = w * half * ctx.P2a * k * k
            for i, p in enumerate(ctx.ps):
                base = None
                for v in variants:
                    if v == "tail":
                        # base plus the piece [TAIL/p, TAIL_EXT/p]; the phase is kept where k < p
                        rate = p + (ctx.b if k < p else 0)
                        a0, a1 = TAIL / p, TAIL_EXT / p
                        n = int(mp.ceil((a1 - a0) * rate / B_VARIANTS["base"][0]))
                        ext = depth_factor_integrals(k, ctx.b, p, [(a0, a1, n)], B_VARIANTS["base"][1])
                        vals = pieces_from_factors(k, ctx.b, *[s + e for s, e in zip(base, ext)])
                    else:
                        size, degree = B_VARIANTS[v]
                        regs = depth_regions(k, ctx.b, p, size)
                        if blackbox:
                            vals = depth_quadrature(k, ctx.b, p, depth_panels(k, ctx.b, p, size), degree)
                        else:
                            facs = depth_factor_integrals(k, ctx.b, p, regs, degree)
                            vals = pieces_from_factors(k, ctx.b, *facs)
                        if v == "base":
                            base = None if blackbox else facs
                            depth_nodes += sum(r[2] for r in regs) * 24
                    wl = wk * p * p
                    acc[v][i] = [s + wl * mp.re(z * l) for s, z, l in zip(acc[v][i], Z, vals)]
    return acc, depth_nodes


class BulkIntegrand:
    """E2: L(q, T) per d(ln q) (times 1/(pi a0 E0) outside)."""

    def __init__(self, ctx):
        self.c = ctx

    def __call__(self, t):
        q = mp.exp(t)
        return [loss_function(self.c.mat, q * q, self.c.T, self.c.h)]


def bulk_13a(spec, k13a):
    """V3: (a) no dispersion: L(0,T) log(q+/q-)/(pi a0 E0); (b) dispersion, k_total: the
    bulk integral over the total momentum q (another variable and other panels than 13a)."""
    mp.mp.dps = DPS_B
    c = Context(spec)
    if c.below:
        return None
    if not c.mat.dispersive:
        v = c.P1 * loss_function(c.mat, mp.mpf(0), c.T, c.h) * mp.log(c.qhi / c.qlo)
        return {"kind": "analytic_no_dispersion", "value": mp.nstr(v, VALUE_DIGITS + 5), "quadrature_u": "0",
                "abs_difference_from_term_13a": mp.nstr(abs(v - mp.mpf(k13a)), 6)}
    if spec["convention"] != "k_total":
        return None
    mp.mp.dps = DPS_A
    ca = Context(spec)
    t0, t1 = mp.log(ca.qlo), mp.log(ca.qhi)
    n0 = int(mp.ceil((t1 - t0) / mp.mpf("0.37")))
    pts = [t0 + (t1 - t0) * i / n0 for i in range(n0 + 1)]
    for E, A, G, alpha in ca.mat.osc:
        if A != 0 and alpha != 0 and ca.T > E:
            tq = mp.log((ca.T - E) / (alpha * ca.h)) / 2
            if t0 < tq < t1:
                pts.append(tq)
    run_a, leaves, _, resid, _ = integrate_adaptive(BulkIntegrand(ca), sorted(set(pts)), [REL_TOL])
    mp.mp.dps = DPS_B
    run_b = integrate_fixed(BulkIntegrand(c), halved(leaves, mp.log(c.qlo), mp.log(c.qhi)))
    u = c.P1 * (2 * abs(mp.mpf(run_a[0]) - run_b[0]) + resid[0])
    return {"kind": "bulk_integral_total_momentum", "value": mp.nstr(c.P1 * run_b[0], VALUE_DIGITS + 5),
            "quadrature_u": mp.nstr(u, 6),
            "abs_difference_from_term_13a": mp.nstr(abs(c.P1 * run_b[0] - mp.mpf(k13a)), 6)}


class LimitIntegrand:
    """Mode coefficients of G_j(a) = Int [c0 + Re(cs e^{-s a}) + cd e^{-2ka}] d(ln k)
    (design 2.2), taken from the printed factors with the mode algebra.

    Smooth columns: G(0) (7), G_inf (7), A1 (7).  Envelope columns for the outputs
    13b..16h, no_hole, hole, total(0.5), total(1): env0 (11), then per long length E_tr (11)
    and the bound of the remainder after the first order (11)."""

    def __init__(self, ctx, long_lengths):
        self.c = ctx
        self.long = [mp.mpf(v) for v in long_lengths]

    def columns(self):
        return 21 + 11 * (1 + 2 * len(self.long))

    def rel(self):
        return [REL_TOL] * 21 + [ENV_REL_TOL] * (11 * (1 + 2 * len(self.long)))

    def modes(self, t):
        """k, s, c0 (7), and the pairs (cs, cd) of the 11 envelope outputs at ln k = t."""
        c = self.c
        k = mp.exp(t)
        s = mp.mpc(k, c.b)
        n = c.b * c.b + k * k
        Z = coefficients(c.S(k), n)
        E, Ei, d = c.symbols
        pieces = printed_pieces(k, c.b, E, Ei, d, Modes.lift(mp.mpf(1)))
        w = c.P2a * k * k
        c0, cs, cd = [], [], []
        for z, piece in zip(Z, pieces):
            c1, cu, cub, cq = piece.coefficients()
            c0.append(w * mp.re(z * c1))
            cs.append(w * (z * cu + mp.conj(z * cub)))
            cd.append(w * mp.re(z * cq))
        se, sh = mp.fsum(cs[i] for i in ELECTRON_IDX), mp.fsum(cs[i] for i in HOLE_IDX)
        de, dh = mp.fsum(cd[i] for i in ELECTRON_IDX), mp.fsum(cd[i] for i in HOLE_IDX)
        half = mp.mpf("0.5")
        pairs = list(zip(cs, cd)) + [(se, de), (sh, dh), (se + half * sh, de + half * dh), (se + sh, de + dh)]
        return k, s, c0, pairs

    def kinks(self, edges):
        """Sign changes of the real coefficients cd (kinks of |cd|), located by bisection.
        Changes below 1e-30 of the largest |cd| are rounding noise of an identical zero."""
        samples = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            samples += [lo + (hi - lo) * i / 12 for i in range(12)]
        samples.append(edges[-1])
        rows = [[cc for _, cc in self.modes(t)[3]] for t in samples]
        scale = max(abs(v) for row in rows for v in row)
        roots = []
        for j in range(11):
            for i in range(len(samples) - 1):
                a, bb = rows[i][j], rows[i + 1][j]
                if a * bb < 0 and max(abs(a), abs(bb)) > mp.mpf("1e-30") * scale:
                    lo, hi, flo = samples[i], samples[i + 1], a
                    for _ in range(mp.mp.prec + 10):
                        mid = (lo + hi) / 2
                        fm = self.modes(mid)[3][j][1]
                        if fm == 0:
                            lo = hi = mid
                            break
                        if (fm < 0) == (flo < 0):
                            lo = mid
                        else:
                            hi = mid
                    roots.append((lo + hi) / 2)
        return roots

    def __call__(self, t):
        k, s, c0, pairs = self.modes(t)
        cs, cd = [bb for bb, _ in pairs[:7]], [cc for _, cc in pairs[:7]]
        out = [a + mp.re(bb) + cc for a, bb, cc in zip(c0, cs, cd)]          # G(0)
        out += c0                                                           # G_inf
        out += [mp.re(bb / s) + cc / (2 * k) for bb, cc in zip(cs, cd)]     # A1
        out += [abs(s * bb) + 2 * k * abs(cc) for bb, cc in pairs]
        for lam in self.long:
            one_s, one_k = abs(1 + lam * s), 1 + 2 * lam * k
            out += [abs(bb) / one_s + abs(cc) / one_k for bb, cc in pairs]
            out += [abs(bb) / (lam * abs(s) * one_s) + abs(cc) / (2 * lam * k * one_k) for bb, cc in pairs]
        return out


def limits(spec):
    want = spec.get("limits")
    if not want:
        return None
    mp.mp.dps = DPS_A
    ca = Context(spec)
    fa = LimitIntegrand(ca, want["long_lengths"])
    edges = ca.edges()
    kinks = fa.kinks(edges)
    edges = sorted(set(edges) | {t for t in kinks if ca.t0 < t < ca.t1})
    run_a, leaves, _, resid, _ = integrate_adaptive(fa, edges, fa.rel())
    mp.mp.dps = DPS_B
    cb = Context(spec)
    run_b = integrate_fixed(LimitIntegrand(cb, want["long_lengths"]), halved(leaves, cb.t0, cb.t1))
    u = [2 * abs(mp.mpf(a) - b) + mp.mpf(r) for a, b, r in zip(run_a, run_b, resid)]
    return {"long_lengths": want["long_lengths"], "panels": len(leaves), "kinks_of_cd": len(kinks),
            "values": [mp.nstr(v, VALUE_DIGITS + 5) for v in run_b], "u": [mp.nstr(v, 6) for v in u]}


def group_task(spec):
    """All lengths of one (material, E0, convention, loss).  Returns strings."""
    cache = spec.get("cache")
    if cache and os.path.exists(cache):
        with open(cache, encoding="ascii") as handle:
            return json.load(handle)
    t_start = time.perf_counter()
    mp.mp.dps = DPS_A
    ca = Context(spec)
    fa = MainIntegrand(ca)
    ncol = 1 + 9 * len(ca.ps)
    run_a, leaves, parents, resid, tol = integrate_adaptive(fa, ca.edges(), [REL_TOL] * ncol)
    mp.mp.dps = DPS_B
    cb = Context(spec)
    fb = MainIntegrand(cb)
    run_b = integrate_fixed(fb, halved(leaves, cb.t0, cb.t1))
    t_a = time.perf_counter() - t_start

    mp.mp.dps = DPS_A
    t_b0 = time.perf_counter()
    rb, depth_nodes = route_b(ca, parents)
    t_b = time.perf_counter() - t_b0

    mp.mp.dps = DPS_B
    nstr = lambda v, digits=6: mp.nstr(v, digits)
    out = {"key": spec["key"], "lengths": spec["lengths"], "below_gap": bool(ca.below),
           "values": [nstr(v, VALUE_DIGITS + 5) for v in run_b],
           "diff_ab": [nstr(abs(mp.mpf(a) - b)) for a, b in zip(run_a, run_b)],
           "residual_a": [nstr(r) for r in resid],
           "route_b": {v: [[nstr(x, VALUE_DIGITS + 5) for x in row] for row in rows] for v, rows in rb.items()},
           "panels_a": len(leaves), "nodes_a": fa.evaluations, "nodes_b_run": fb.evaluations,
           "route_b_k_nodes": 24 * len(parents), "route_b_depth_nodes_base": depth_nodes,
           "k_limits": [nstr(cb.lo, 30), nstr(cb.hi, 30)], "omega_perp": nstr(cb.b, 30),
           "seconds_route_a": t_a, "seconds_route_b": t_b}
    t0 = time.perf_counter()
    out["bulk_13a"] = bulk_13a(spec, run_b[0])
    out["limits"] = limits(spec)
    out["seconds_extra"] = time.perf_counter() - t0
    out["seconds"] = time.perf_counter() - t_start
    if cache:
        tmp = cache + ".tmp"
        with open(tmp, "w", encoding="ascii", newline="\n") as handle:
            json.dump(out, handle)
        os.replace(tmp, cache)
    return out


# Black-box check of route (b): (material, convention, loss, length).  The whole printed
# pieces are evaluated at every depth node (base rule) and compared with route (b).
BLACKBOX = (("W4", "k_parallel", math.sqrt(135.0), 20.0), ("Al", "k_total", 15.6, 20.0),
            ("Si", "k_total", 0.2, 5.0), ("S0", "k_parallel", 5.0, 80.0))


def blackbox_task(spec):
    cache = spec.get("cache")
    if cache and os.path.exists(cache):
        with open(cache, encoding="ascii") as handle:
            return json.load(handle)
    start = time.perf_counter()
    mp.mp.dps = DPS_A
    ca = Context(spec)
    _, _, parents, _, _ = integrate_adaptive(MainIntegrand(ca), ca.edges(), [REL_TOL] * (1 + 9 * len(ca.ps)))
    t0 = time.perf_counter()
    fast, _ = route_b(ca, parents, variants=("base",))
    t_fast = time.perf_counter() - t0
    t0 = time.perf_counter()
    full, _ = route_b(ca, parents, variants=("base",), blackbox=True)
    t_full = time.perf_counter() - t0
    rel = [abs(x - y) / abs(y) for x, y in zip(fast["base"][0], full["base"][0])]
    out = {"key": spec["key"], "length_A_hex": fhex(spec["lengths"][0]),
           "max_relative_difference": mp.nstr(max(rel), 6),
           "relative_difference": dict(zip(DEPTH_TERMS, (mp.nstr(v, 6) for v in rel))),
           "seconds_factor_route": round(t_fast, 2), "seconds_whole_pieces": round(t_full, 2),
           "seconds": round(time.perf_counter() - start, 2)}
    if cache:
        with open(cache + ".tmp", "w", encoding="ascii") as handle:
            json.dump(out, handle)
        os.replace(cache + ".tmp", cache)
    return out


def any_task(spec):
    return blackbox_task(spec) if spec.get("blackbox") else group_task(spec)


# ----------------------------------------------------------------------------- table

def load_table():
    raw = TABLE.read_bytes()
    if hashlib.sha256(raw).hexdigest() != TABLE_SHA256:
        raise SystemExit("expected_states_m3a.json does not match the pre-registered SHA-256")
    return json.loads(raw.decode("utf-8"))


def group_key(material, energy, convention, loss):
    return f"{material}|{fhex(energy)}|{convention}|{fhex(loss)}"


def build_groups(table):
    groups = {}
    for r in table["reference_keys"]:
        gk = group_key(r["material"], r["electron_energy_eV"], r["epsilon_momentum"], r["loss_eV"])
        m = table["materials"][r["material"]]
        g = groups.setdefault(gk, {"key": gk, "material": r["material"], "gap": m["gap_eV"],
                                   "oscillators": m["oscillators"], "energy": r["electron_energy_eV"],
                                   "convention": r["epsilon_momentum"], "loss": r["loss_eV"], "lengths": []})
        g["lengths"].append(r["weight_length_A"])
    for r in table["limit_keys"]:
        gk = group_key(r["material"], r["electron_energy_eV"], r["epsilon_momentum"], r["loss_eV"])
        if gk not in groups:
            raise SystemExit(f"limit key without a reference group: {r['key']}")
        lim = groups[gk].setdefault("limits", {"long_lengths": []})
        if r["kind"] == "lambda_to_inf":
            lim["long_lengths"] = sorted(set(lim["long_lengths"]) | set(r["weight_lengths_A"]))
    out = []
    for g in groups.values():
        g["lengths"] = sorted(set(g["lengths"]))
        out.append(g)
    return out


# ----------------------------------------------------------------------------- V1 (exact rationals)

V1_COMPONENTS = ("W13", "N", "W14", "W15b", "W15c", "W16e", "W16h")


class CQ:
    """Complex rational."""

    __slots__ = ("re", "im")

    def __init__(self, re, im=0):
        self.re, self.im = Fraction(re), Fraction(im)

    def __mul__(self, o):
        o = o if isinstance(o, CQ) else CQ(o)
        return CQ(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)

    def __add__(self, o):
        o = o if isinstance(o, CQ) else CQ(o)
        return CQ(self.re + o.re, self.im + o.im)

    def __neg__(self):
        return CQ(-self.re, -self.im)

    def inv(self):
        m = self.re * self.re + self.im * self.im
        return CQ(self.re / m, -self.im / m)


def weights_complex(k, b, p, mutation=None):
    """Scaled weights p^2 W/(2b) from the complex forms of design 3.1 (exact rationals)."""
    s, sc = CQ(k, b), CQ(k, -b)
    if mutation == "K2":
        s, sc = sc, s
    Jw = p + 2 * k if mutation != "K1" else p + k
    D = k * k + b * b
    i = CQ(0, 1)
    ps, psc = s + p, sc + p
    f14 = 2 if mutation != "K3" else 1
    sign15b = 1 if mutation != "K4" else -1
    w = {"W13": -(s * ps * p).inv(),
         "N": CQ(2 * (p + k) / (p * Jw * (ps * psc).re)),
         "W14": i * (f14 * k) * (ps * (p * Jw)).inv(),
         "W15b": -(i * (ps * Jw).inv()) * sign15b,
         "W15c": -(sc * psc * p).inv(),
         "W16e": -(ps * (D * p)).inv() * Jw,
         "W16h": i * (s * ps).inv()}
    scale = p * p / (2 * b)
    return {name: CQ(v.re * scale, v.im * scale) for name, v in w.items()}


def weights_real(k, b, p, mutation=None):
    """Scaled weights from the real and imaginary parts of the table of design 3.1."""
    L = p + k
    Jw = L + k if mutation != "K1" else L
    D = k * k + b * b
    H = L * L + b * b
    n = k * L - b * b
    f14 = 2 if mutation != "K3" else 1
    sign15b = 1 if mutation != "K4" else -1
    sb = b if mutation != "K2" else -b              # s <-> conj(s) is b -> -b in the parts
    w = {"W13": (-n / (p * D * H), sb * Jw / (p * D * H)),
         "N": (2 * L / (p * Jw * H), Fraction(0)),
         "W14": (f14 * k * sb / (p * Jw * H), f14 * k * L / (p * Jw * H)),
         "W15b": (sign15b * -sb / (Jw * H), sign15b * -L / (Jw * H)),
         "W15c": (-n / (p * D * H), -sb * Jw / (p * D * H)),
         "W16e": (-Jw * L / (p * D * H), sb * Jw / (p * D * H)),
         "W16h": (sb * Jw / (D * H), n / (D * H))}
    scale = p * p / (2 * b)
    return {name: (v[0] * scale, v[1] * scale) for name, v in w.items()}, n


def v1_widths(k, b, p):
    """V1 value widths of the scaled weights: relative 1e-13, or the absolute widths of the
    three parts with a zero crossing."""
    L = p + k
    D = k * k + b * b
    H = L * L + b * b
    scale = p * p / (2 * b)
    tol = Fraction(1, 10 ** 13)
    return {("W13", 0): tol * (k * L + b * b) / (p * D * H) * scale,
            ("W15c", 0): tol * (k * L + b * b) / (p * D * H) * scale,
            ("W16h", 1): tol * (k * L + b * b) / (D * H) * scale}


def frac_json(x, digits=45):
    x = Fraction(x)
    mp.mp.dps = digits + 10
    return {"decimal": mp.nstr(mp.mpf(x.numerator) / x.denominator, digits),
            "numerator": str(x.numerator), "denominator": str(x.denominator)}


def v1_fixed_depth_factors(k, b, a):
    """The fixed-depth factors of the seven weights (design 3.1, column "fixed depth"):
    plus = i/s, minus = conj(plus), diff = 1-u, norm = |1-u|^2/D, mixed = u-|u|^2."""
    s = mp.mpc(k, b)
    u = mp.exp(-s * a)
    D = k * k + b * b
    plus, minus = J / s, -J / mp.conj(s)
    diff = 1 - u
    norm = (diff * mp.conj(diff)).real / D
    mixed = u - mp.exp(-2 * k * a)
    return [-diff / (s * s), mp.mpc(norm), plus * diff + minus * mixed, minus * mixed,
            minus * minus * mp.conj(diff), plus * (2 * J * k / D - plus * u), plus * u]


def v1_quadrature(kf, bf, pf):
    """Scaled weights by depth quadrature in t = p a (phase-split panels), with the changes
    from the order, the panel split and the tail."""
    extra = 0
    for v in (kf, bf):
        ratio = Fraction(v) / Fraction(pf)
        if ratio < 1:
            extra = max(extra, 2 * int(-math.log10(float(ratio)) + 1))
    mp.mp.dps = 60 + extra                       # (s a)^2 cancellations of |1-u|^2 for k, b << p
    k, b, p = mp.mpf(kf), mp.mpf(bf), mp.mpf(pf)
    scale = p * p / (2 * b)

    def run(size, degree, tail):
        nodes = gl_nodes(degree)
        acc = [mp.mpc(0)] * 7
        for lo, hi in depth_panels(k, b, p, size, tail):
            half, mid = (hi - lo) / 2, (hi + lo) / 2
            for x, w in nodes:
                a = mid + half * x
                wa = w * half * mp.exp(-p * a)
                acc = [s + wa * v for s, v in zip(acc, v1_fixed_depth_factors(k, b, a))]
        return [scale * v for v in acc]

    base = run(10, 4, 64)
    change = [mp.mpf(0)] * 7
    for size, degree, tail in ((20, 5, 64), (7, 4, 64), (10, 4, 80)):
        other = run(size, degree, tail)
        change = [max(c, abs(x - y)) for c, x, y in zip(change, base, other)]
    exact = weights_complex(Fraction(kf), Fraction(bf), Fraction(pf))
    rows = {}
    for name, val, ch in zip(V1_COMPONENTS, base, change):
        ex = mp.mpc(mp.mpf(exact[name].re.numerator) / exact[name].re.denominator,
                    mp.mpf(exact[name].im.numerator) / exact[name].im.denominator)
        rows[name] = {"re": mp.nstr(val.real, 40), "im": mp.nstr(val.imag, 40),
                      "variation": mp.nstr(ch, 6), "abs_difference_from_exact": mp.nstr(abs(val - ex), 6),
                      "relative_difference": mp.nstr(abs(val - ex) / abs(ex), 6)}
    return {"digits": mp.mp.dps, "components": rows}


def v1_points(table):
    mp.mp.dps = 60
    pts = [("k<<p", 1e-3, 0.5, 1.0), ("k~p", 0.05, 0.03, 0.05), ("k>>p", 5.0, 0.03, 0.05),
           ("b>>p", 0.01, 2.0, 0.05),
           ("zero_crossing_before", 0.5, 1.0, math.nextafter(1.5, 0.0)),
           ("zero_crossing", 0.5, 1.0, 1.5),
           ("zero_crossing_after", 0.5, 1.0, math.nextafter(1.5, math.inf)),
           ("underflow", 1.0, 1.0, 1e82)]
    h, E0 = mp.mpf(b64(H_DEC)), mp.mpf(1000.0)
    for T in (1.0, 15.6, 99.8):
        b, lo, hi, _, _ = kinematics(mp.mpf(T), E0, h)
        for lam in (5.0, 20.0, 80.0):
            for tag, kv in (("k_par_minus", lo), ("k_par_geometric_mean", mp.sqrt(lo * hi)), ("k_par_plus", hi)):
                pts.append((f"Al|T={T!r}|lambda={lam!r}|{tag}", float(kv), float(b), 1.0 / lam))
    return pts


def v1_block(table):
    pts = v1_points(table)
    points = []
    mutation_hits = {m: [] for m in ("K1", "K2", "K3", "K4")}
    for idx, (label, kf, bf, pf) in enumerate(pts):
        k, b, p = Fraction(kf), Fraction(bf), Fraction(pf)
        wc = weights_complex(k, b, p)
        wr, n = weights_real(k, b, p)
        for name in V1_COMPONENTS:
            if (wc[name].re, wc[name].im) != wr[name]:
                raise SystemExit(f"V1: complex and real forms differ at {label} {name}")
        widths = v1_widths(k, b, p)
        entry = {"label": label, "k_hex": kf.hex(), "b_hex": bf.hex(), "p_hex": pf.hex(),
                 "n": frac_json(n),
                 "weights": {name: {"re": frac_json(wr[name][0]), "im": frac_json(wr[name][1])}
                             for name in V1_COMPONENTS}}
        if idx < 8:
            entry["depth_quadrature"] = v1_quadrature(kf, bf, pf)
        for m in mutation_hits:
            wm, _ = weights_real(k, b, p, m)
            wmc = weights_complex(k, b, p, m)
            hit = []
            for name in V1_COMPONENTS:
                # K1 is defined on the parts of the table (J also appears there in W13, W15c, W16h)
                if m != "K1" and (wmc[name].re, wmc[name].im) != wm[name]:
                    raise SystemExit(f"V1: mutated forms differ at {label} {name} {m}")
                for part in (0, 1):
                    width = widths.get((name, part), Fraction(1, 10 ** 13) * abs(wr[name][part]))
                    if abs(wm[name][part] - wr[name][part]) > width:
                        hit.append(f"{name}.{'re' if part == 0 else 'im'}")
            mutation_hits[m].append(hit)
        if idx in (0, 1, 2, 3, 6):
            corners = []
            for sk in (-1, 1):
                for sb in (-1, 1):
                    for sp in (-1, 1):
                        kc = k + sk * 4 * Fraction(math.ulp(kf))
                        bc = b + sb * 4 * Fraction(math.ulp(bf))
                        pc = p + sp * 4 * Fraction(math.ulp(pf))
                        wcn, _ = weights_real(kc, bc, pc)
                        corners.append({"signs": [sk, sb, sp],
                                        "weights": {name: {"re": frac_json(wcn[name][0]), "im": frac_json(wcn[name][1])}
                                                    for name in V1_COMPONENTS}})
            entry["input_error_ulps"] = 4
            entry["box_corners"] = corners
        points.append(entry)

    # zero crossing just after: binary64 sequence of design 3.1 and the e_n floor
    kf, bf, pf = 0.5, 1.0, math.nextafter(1.5, math.inf)
    L = pf + kf
    n_hat = kf * L - bf * bf
    k, b, p = Fraction(kf), Fraction(bf), Fraction(pf)
    wr, n_exact = weights_real(k, b, p)
    ur = Fraction(1, 2 ** 53)
    e_L = ur * abs(Fraction(L))                                  # rounding of L = p + k
    e_n = abs(k) * e_L + ur * (abs(k * Fraction(L)) + b * b) + ur * abs(Fraction(n_hat))
    D = k * k + b * b
    Hx = (p + k) ** 2 + b * b
    scale = p * p / (2 * b)
    floor_re = e_n / (p * D * Hx) * scale
    floor_im = e_n / (D * Hx) * scale
    after = {"n_hat_binary64_hex": float(n_hat).hex(), "n_exact": frac_json(n_exact),
             "L_hat_hex": float(L).hex(),
             "e_n_with_zero_input_errors": frac_json(e_n),
             "e_n_formula": "|k| e_L + u_r (|k L| + b^2) + u_r |n_hat|, e_L = u_r |L_hat|, u_r = 2^-53 (design 3.1 with e_k = e_b = e_p = 0)",
             "e_n_covers_exact_n": e_n >= abs(n_exact),
             "terms_proportional_to_n_hat_vanish": n_hat == 0.0,
             "e_n_part_of_error": {"W13.re": frac_json(floor_re), "W15c.re": frac_json(floor_re),
                                   "W16h.im": frac_json(floor_im)},
             "exact_magnitude": {"W13.re": frac_json(abs(wr["W13"][0])), "W15c.re": frac_json(abs(wr["W15c"][0])),
                                 "W16h.im": frac_json(abs(wr["W16h"][1]))},
             "e_n_part_covers_exact_value": bool(floor_re >= abs(wr["W13"][0]) and floor_im >= abs(wr["W16h"][1]))}
    labels = [pt[0] for pt in pts]
    reg = {"K1": {"k~p": ["W13", "W14", "W15b", "W16e"], "k>>p": ["W13", "W14", "W15b", "W16e"]},
           "K2": {"k~p": ["W13.im", "W16e.im", "W16h.re"]},
           "K3": {lab: ["W14"] for lab in labels},
           "K4": {lab: ["W15b"] for lab in labels}}
    mutations = {}
    for m, hits in mutation_hits.items():
        registered = True
        for lab, comps in reg[m].items():
            got = hits[labels.index(lab)]
            for comp in comps:
                if not any(g == comp or g.startswith(comp + ".") for g in got):
                    registered = False
        mutations[m] = {"points_with_a_component_beyond_the_width": sum(1 for hset in hits if hset),
                        "components_beyond_the_width": {lab: hset for lab, hset in zip(labels, hits)},
                        "registered_detection_of_prereg_table_holds": registered}
    names = {"K1": "p+2k -> p+k", "K2": "s <-> conj(s)", "K3": "weight_14 factor 2 -> 1", "K4": "sign of W15b"}
    for m in mutations:
        mutations[m]["mutation"] = names[m]
    kf, bf, pf = 1.0, 1.0, 1e82
    w14 = weights_real(Fraction(kf), Fraction(bf), Fraction(pf))[0]["W14"][0]
    return {"components": list(V1_COMPONENTS),
            "definition": "scaled weight p^2 W/(2b); exact rationals from the binary64 inputs",
            "points": points, "zero_crossing_after": after, "mutations": mutations,
            "underflow_point_scaled_ReW14": frac_json(w14),
            "underflow_point_relative_difference_from_1.0000000000000002e-164":
                float(abs(w14 - Fraction("1.0000000000000002e-164")) / w14),
            "out_of_range_control": {"k_hex": (1.0).hex(), "b_hex": (1.0).hex(), "p_hex": (1e200).hex(),
                                     "expect": "FloatingPointError"}}


# ----------------------------------------------------------------------------- assembly

def kk_cross_check():
    """E4b against the partial fractions at g = 0, and the partial fractions against a direct
    principal-value quadrature of E4a for g > 0."""
    mp.mp.dps = 40
    h = mp.mpf(b64(H_DEC))
    rows = []
    al = Mat(0.0, [[15.6, 243.3598, 1.5, 1.0]])
    for q, T in (("0.3", 5.0), ("1.7", 15.6), ("0.02", 0.05), ("9", 99.8)):
        q2 = mp.mpf(q) ** 2
        a = response(al, q2, mp.mpf(T), h)
        bb = response(al, q2, mp.mpf(T), h, force_kk=True)
        rows.append({"material": "Al", "q_inv_A": q, "loss_eV": T, "check": "E4b vs partial fractions (g = 0)",
                     "abs_difference": float(abs(a - bb))})
    si = Mat(1.12, [[10.0, 6.08458, 5.0, 0.5], [14.0, 30.4229, 5.0, 0.5], [16.8, 212.9603, 3.8, 0.5]])
    n2 = Mat(2.0, [[3.0, 6.0, 0.2, 0.5]])
    for name, mat, q, T in (("Si", si, "0.3", 0.2), ("Si", si, "1", 5.0), ("Si", si, "0.3", 1.120001),
                            ("N2", n2, "0.7245", 4.0)):
        q2 = mp.mpf(q) ** 2
        Tm = mp.mpf(T)
        g = mat.g

        def B(u):
            return mp.fsum(A * G * u * u / (((E + al_ * h * q2) ** 2 - u * u) ** 2 + G * G * u * u)
                           for E, A, G, al_ in mat.osc)

        scale = max(o[0] + o[3] * h * q2 for o in mat.osc)
        marks = sorted({g + (scale - g) * mp.mpf(f) for f in ("0.25", "0.5", "0.75", "1", "1.5", "2", "4", "16")})
        if Tm < g:
            val = mp.quad(lambda u: B(u) / (u * u - Tm * Tm), [g] + [m for m in marks if m > g] + [mp.inf])
        else:
            dlt = Tm - g
            hfun = lambda u: B(u) / (u + Tm)
            core = mp.quad(lambda s: (hfun(Tm + s) - hfun(Tm - s)) / s, mp.linspace(0, dlt, 5))
            up = Tm + dlt
            val = core + mp.quad(lambda u: B(u) / (u * u - Tm * Tm), [up] + [m for m in marks if m > up] + [mp.inf])
        closed = mp.re(response(mat, q2, Tm, h))
        rows.append({"material": name, "q_inv_A": q, "loss_eV": T, "check": "partial fractions vs direct PV quadrature",
                     "abs_difference": float(abs(closed - 2 * val / mp.pi))})
    return rows


def formula_self_check():
    """compact_laplace (route a) against the mode-by-mode transform of the printed factors."""
    mp.mp.dps = 60
    worst = mp.mpf(0)
    E, Ei, d = mode_symbols()
    for k, b, p in (("0.3", "0.07", "0.05"), ("5", "0.03", "0.2"), ("0.002", "0.5", "0.0125"),
                    ("1.3", "0.81", "1e-6"), ("0.4", "0.1", "1e6"), ("31", "0.2", "1")):
        k, b, p = mp.mpf(k), mp.mpf(b), mp.mpf(p)
        pieces = printed_pieces(k, b, E, Ei, d, Modes.lift(mp.mpf(1)))
        for piece, compact in zip(pieces, compact_laplace(k, b, p)):
            ref = piece.laplace(k, b, p)
            worst = max(worst, abs(ref - compact) / abs(compact))
    # the fast numeric evaluator against the generic printed forms
    k, b, a = mp.mpf("0.37"), mp.mpf("0.11"), mp.mpf("2.3")
    Ev = mp.expj(b * a)
    dv = mp.exp(-k * a)
    fast = NumericPieces(k, b).at(Ev, mp.conj(Ev), dv)
    slow = printed_pieces(k, b, Ev, mp.conj(Ev), dv)
    worst_fast = max(abs(x - y) / abs(y) for x, y in zip(fast, slow))
    # pieces_from_factors with pointwise factors against the printed forms at single depths
    worst_fac = mp.mpf(0)
    for k, b, a in (("0.37", "0.11", "2.3"), ("4.1", "0.8", "0.6"), ("0.003", "0.09", "40"), ("12", "0.01", "1e-3")):
        k, b, a = mp.mpf(k), mp.mpf(b), mp.mpf(a)
        Ev, dv = mp.expj(b * a), mp.exp(-k * a)
        Eiv = mp.conj(Ev)
        point = pieces_from_factors(k, b, mp.mpf(1), Eiv * dv, dv * dv, ((Eiv - dv) * (Ev - dv)).real)
        slow = printed_pieces(k, b, Ev, Eiv, dv)
        worst_fac = max(worst_fac, max(abs(x - y) / abs(y) for x, y in zip(point, slow)))
    # fixed-point depth quadrature against the mpf version with direct exponentials
    mp.mp.dps = DPS_A
    worst_fx = mp.mpf(0)
    for k, b, p in (("0.3", "0.07", "0.05"), ("5", "0.03", "0.2"), ("0.002", "0.5", "0.0125"),
                    ("1.3", "0.81", "1e6"), ("0.4", "0.1", "1e-6"), ("31", "0.2", "0.0125")):
        k, b, p = mp.mpf(k), mp.mpf(b), mp.mpf(p)
        for size, degree in B_VARIANTS.values():
            fx = depth_factor_integrals(k, b, p, depth_regions(k, b, p, size), degree)
            fm = depth_factor_integrals_mpf(k, b, p, depth_panels(k, b, p, size), degree)
            worst_fx = max(worst_fx, max(abs(x - y) / abs(y) for x, y in zip(fx, fm)))
    return {"compact_laplace_vs_mode_transform_max_relative": float(worst),
            "numeric_pieces_vs_printed_pieces_max_relative": float(worst_fast),
            "pieces_from_factors_vs_printed_pieces_pointwise_max_relative": float(worst_fac),
            "fixed_point_vs_mpf_depth_quadrature_max_relative": float(worst_fx)}


def up(x):
    """Decimal string of an upper bound of the non-negative mpf x."""
    if x == 0:
        return "0"
    return mp.nstr(x * (1 + mp.mpf("1e-4")), 6)


def dec_round(v_str_long):
    """(stored 40-digit decimal, rounding width) of a longer decimal string."""
    v = mp.mpf(v_str_long)
    s = mp.nstr(v, VALUE_DIGITS)
    return s, abs(mp.mpf(s) - v)


def assemble(table, results, t_wall, processes, blackbox=()):
    mp.mp.dps = DPS_B
    by_key = {r["key"]: r for r in results}
    runs = {r["id"]: r for r in table["runs"]}
    users = {}
    for check in ("V7", "V11iii", "V11v", "V3d"):
        for run_id, _, key in table["check_plan"][check]:
            if key is None or key == "analytic_zero":
                continue
            run = runs[run_id]
            users.setdefault(key, set()).add((run["settings"], run["hole_charge"]["value"]))
    references = {}
    worst = (mp.mpf(0), None)
    undetermined = []
    max_resid_zero = mp.mpf(0)
    for r in table["reference_keys"]:
        gk = group_key(r["material"], r["electron_energy_eV"], r["epsilon_momentum"], r["loss_eV"])
        g = by_key[gk]
        i = g["lengths"].index(r["weight_length_A"])
        col = lambda arr, j: mp.mpf(arr[j])
        below = g["below_gap"]
        vals, us, parts = {}, {}, {}
        # term_13a
        raw13a = col(g["values"], 0)
        s13, dr = dec_round(g["values"][0])
        parts["term_13a"] = {"diff_50_80": up(col(g["diff_ab"], 0)), "residual": up(col(g["residual_a"], 0)),
                             "decimal": up(dr)}
        vals["term_13a"] = mp.mpf(s13)
        us["term_13a"] = 2 * col(g["diff_ab"], 0) + col(g["residual_a"], 0) + dr
        base = 1 + 9 * i
        rb = g["route_b"]
        var = []
        for j, name in enumerate(DEPTH_TERMS):
            a_val = col(g["values"], base + j)
            b_val = mp.mpf(rb["base"][i][j])
            s, dr = dec_round(g["values"][base + j])
            d_ab = abs(a_val - b_val)
            ch = {v: abs(mp.mpf(rb[v][i][j]) - b_val) for v in ("order", "phase", "tail")}
            var.append(max(ch.values()))
            vals[name] = mp.mpf(s)
            us[name] = d_ab + 2 * col(g["diff_ab"], base + j) + col(g["residual_a"], base + j) + var[-1] + dr
            parts[name] = {"a_minus_b": up(d_ab), "diff_50_80": up(col(g["diff_ab"], base + j)),
                           "residual": up(col(g["residual_a"], base + j)), "order": up(ch["order"]),
                           "phase": up(ch["phase"]), "tail": up(ch["tail"]), "decimal": up(dr)}
        for name, idx, off in (("electron_group", ELECTRON_IDX, 7), ("hole_group", HOLE_IDX, 8)):
            a_val = col(g["values"], base + off)
            s, dr = dec_round(g["values"][base + off])
            b_sum = mp.fsum(mp.mpf(rb["base"][i][j]) for j in idx)
            d_asm = abs(a_val - b_sum)
            members = mp.fsum(var[j] for j in idx)
            vals[name] = mp.mpf(s)
            us[name] = d_asm + 2 * col(g["diff_ab"], base + off) + col(g["residual_a"], base + off) + members + dr
            parts[name] = {"group_minus_sum_of_route_b_terms": up(d_asm), "diff_50_80": up(col(g["diff_ab"], base + off)),
                           "residual": up(col(g["residual_a"], base + off)), "route_b_changes_of_members": up(members),
                           "decimal": up(dr)}
        zero, residual = [], {}
        if below:
            residual["term_13a"] = up(abs(raw13a))
            raw_total = raw13a + mp.fsum(col(g["values"], base + j) for j in range(7))
            residual["total_unit_charge_from_raw_terms"] = up(abs(raw_total))
            residual["total_unit_charge_from_groups"] = up(abs(col(g["values"], base + 7) + col(g["values"], base + 8)))
            max_resid_zero = max(max_resid_zero, abs(raw_total), abs(raw13a))
            vals["term_13a"], us["term_13a"] = mp.mpf(0), mp.mpf(0)
            zero = ["term_13a", "total_unit_charge"]
        vals["no_hole"] = vals["term_13a"] + vals["electron_group"]
        us["no_hole"] = us["term_13a"] + us["electron_group"]
        vals["hole"], us["hole"] = vals["hole_group"], us["hole_group"]
        outputs = {}
        for name in TERMS + ("no_hole", "hole", "electron_group", "hole_group"):
            v = vals[name]
            if name == "no_hole":
                txt = mp.nstr(v, VALUE_DIGITS + 2)
                extra = abs(mp.mpf(txt) - v)
            else:
                txt, extra = (mp.nstr(v, VALUE_DIGITS) if v != 0 else "0"), mp.mpf(0)
            outputs[name] = {"value": txt, "u": up(us[name] + extra)}
        # tightest use
        ratio_max, at = mp.mpf(0), None
        for setting, charge in sorted(users.get(r["key"], ())):
            st = table["settings"][setting]
            ea, er = mp.mpf(st["epsabs"]), mp.mpf(st["epsrel"])
            q = mp.mpf(charge)
            cand = [(n, vals[n], us[n]) for n in TERMS + ("no_hole", "hole")]
            if not (below and charge == 1.0):
                cand.append(("total", vals["no_hole"] + q * vals["hole"], us["no_hole"] + abs(q) * us["hole"]))
            if below and charge != 1.0:
                cand.append(("total_as_(1-Q)K_e", (1 - q) * vals["electron_group"], abs(1 - q) * us["electron_group"]))
            for n, v, u in cand:
                ratio = u / max(ea, er * abs(v))
                if ratio > ratio_max:
                    ratio_max, at = ratio, [n, setting, charge]
        if r["key"] not in users:
            raise SystemExit(f"reference key not used by any check: {r['key']}")
        if ratio_max > worst[0]:
            worst = (ratio_max, [r["key"]] + at)
        if ratio_max > mp.mpf("0.01"):
            undetermined.append({"key": r["key"], "u_over_Bref": float(ratio_max), "at": at})
        entry = {"material": r["material"], "electron_energy_eV_hex": fhex(r["electron_energy_eV"]),
                 "epsilon_momentum": r["epsilon_momentum"], "weight_length_A_hex": fhex(r["weight_length_A"]),
                 "loss_eV_hex": fhex(r["loss_eV"]), "region": "below_gap" if below else "above_gap",
                 "group": gk, "analytic_zero": zero, "numerical_residual_of_zero": residual,
                 "outputs": outputs, "u_parts": parts,
                 "max_u_over_Bref": float(ratio_max), "max_u_over_Bref_at": at}
        references[r["key"]] = entry

    groups = {}
    for g in results:
        cols = g["values"]
        entry = {"lengths_A_hex": [fhex(v) for v in g["lengths"]], "below_gap": g["below_gap"],
                 "k_parallel_limits_inv_A": g["k_limits"], "omega_perp_inv_A": g["omega_perp"],
                 "panels_run_A": g["panels_a"], "nodes_run_A": g["nodes_a"], "nodes_run_B": g["nodes_b_run"],
                 "route_b_k_nodes": g["route_b_k_nodes"], "route_b_depth_nodes_base": g["route_b_depth_nodes_base"],
                 "seconds_route_a": round(g["seconds_route_a"], 2), "seconds_route_b": round(g["seconds_route_b"], 2),
                 "seconds": round(g["seconds"], 2)}
        if g["bulk_13a"]:
            bulk = g["bulk_13a"]
            s, dr = dec_round(bulk["value"])
            entry["term_13a_independent"] = {"kind": bulk["kind"], "value": s,
                                             "u": up(mp.mpf(bulk["quadrature_u"]) + dr),
                                             "abs_difference_from_term_13a": bulk["abs_difference_from_term_13a"]}
        groups[g["key"]] = entry

    limits_out = {}
    names7 = DEPTH_TERMS
    for r in table["limit_keys"]:
        gk = group_key(r["material"], r["electron_energy_eV"], r["epsilon_momentum"], r["loss_eV"])
        g = by_key[gk]
        lim = g["limits"]
        V = [mp.mpf(v) for v in lim["values"]]
        U = [mp.mpf(v) for v in lim["u"]]
        q = r["hole_charge"]
        qi = {0.0: None, 0.5: 9, 1.0: 10}[q]                # envelope column of total(Q); Q = 0 is no_hole
        k13 = mp.mpf(g["values"][0]) if not g["below_gap"] else mp.mpf(0)
        u13 = (2 * mp.mpf(g["diff_ab"][0]) + mp.mpf(g["residual_a"][0])) if not g["below_gap"] else mp.mpf(0)

        def linear(off):
            """values and u of the 11 outputs of a linear quantity stored as 7 columns at off."""
            v = {n: V[off + j] for j, n in enumerate(names7)}
            u = {n: U[off + j] for j, n in enumerate(names7)}
            v["term_13a"], u["term_13a"] = mp.mpf(0), mp.mpf(0)
            v["no_hole"] = mp.fsum(V[off + j] for j in ELECTRON_IDX)
            u["no_hole"] = mp.fsum(U[off + j] for j in ELECTRON_IDX)
            v["hole"] = mp.fsum(V[off + j] for j in HOLE_IDX)
            u["hole"] = mp.fsum(U[off + j] for j in HOLE_IDX)
            v["total"] = v["no_hole"] + mp.mpf(q) * v["hole"]
            u["total"] = u["no_hole"] + abs(mp.mpf(q)) * u["hole"]
            return v, u

        def envelope(off):
            v = {n: V[off + j] for j, n in enumerate(names7)}
            u = {n: U[off + j] for j, n in enumerate(names7)}
            v["term_13a"], u["term_13a"] = mp.mpf(0), mp.mpf(0)
            v["no_hole"], u["no_hole"] = V[off + 7], U[off + 7]
            v["hole"], u["hole"] = V[off + 8], U[off + 8]
            j = 7 if qi is None else qi
            v["total"], u["total"] = V[off + j], U[off + j]
            return v, u

        def pack(v, u, digits=VALUE_DIGITS):
            out = {}
            for n in OUTPUTS:
                txt = mp.nstr(v[n], digits) if v[n] != 0 else "0"
                out[n] = {"value": txt, "u": up(u[n] + abs(mp.mpf(txt) - v[n]))}
            return out

        entry = {"kind": r["kind"], "material": r["material"], "epsilon_momentum": r["epsilon_momentum"],
                 "loss_eV_hex": fhex(r["loss_eV"]), "hole_charge_hex": fhex(q),
                 "weight_lengths_A_hex": [fhex(v) for v in r["weight_lengths_A"]], "group": gk}
        if r["kind"] == "lambda_to_0":
            entry["G0"] = pack(*linear(0))
            entry["envelope"] = pack(*envelope(21), digits=20)
            kv = {n: (k13 if n in ("term_13a", "no_hole", "total") else mp.mpf(0)) for n in OUTPUTS}
            ku = {n: (u13 if n in ("term_13a", "no_hole", "total") else mp.mpf(0)) for n in OUTPUTS}
            entry["kappa"] = pack(kv, ku)
        else:
            entry["G_inf"] = pack(*linear(7))
            entry["A1"] = pack(*linear(14))
            entry["lengths"] = {}
            for lam in r["weight_lengths_A"]:
                m = lim["long_lengths"].index(lam)
                entry["lengths"][fhex(lam)] = {"E_tr": pack(*envelope(32 + 22 * m), digits=20),
                                               "remainder_bound_after_first_order": pack(*envelope(43 + 22 * m), digits=20)}
        limits_out[r["key"]] = entry

    # Internal consistency of the two code paths: the averaged references (compact forms) at
    # the extreme lengths against the limit quantities (mode coefficients).
    limit_checks = {"lambda_to_0": {"cases": 0, "max_deviation_over_bound": 0.0},
                    "lambda_to_inf_envelope": {"cases": 0, "max_deviation_over_bound": 0.0},
                    "lambda_to_inf_first_order": {"cases": 0, "max_deviation_over_bound": 0.0}}

    def note(name, lhs, bound):
        if bound == 0 and lhs == 0:
            return
        item = limit_checks[name]
        item["cases"] += 1
        item["max_deviation_over_bound"] = max(item["max_deviation_over_bound"],
                                               float(lhs / bound) if bound != 0 else math.inf)

    for r in table["limit_keys"]:
        entry = limits_out[r["key"]]
        q = mp.mpf(r["hole_charge"])
        val = lambda field, n, e=entry: (mp.mpf(e[field][n]["value"]), mp.mpf(e[field][n]["u"]))
        for lam in r["weight_lengths_A"]:
            rk = f'{r["material"]}|{r["electron_energy_eV"]!r}|{r["epsilon_momentum"]}|{fhex(lam)}|{fhex(r["loss_eV"])}'
            if rk not in references:
                continue
            out = references[rk]["outputs"]
            lm = mp.mpf(lam)
            K = {n: mp.mpf(out[n]["value"]) for n in TERMS + ("no_hole", "hole")}
            U = {n: mp.mpf(out[n]["u"]) for n in TERMS + ("no_hole", "hole")}
            K["total"], U["total"] = K["no_hole"] + q * K["hole"], U["no_hole"] + abs(q) * U["hole"]
            for n in OUTPUTS:
                if r["kind"] == "lambda_to_0":
                    (g0, ug), (env, ue), (kap, uk) = val("G0", n), val("envelope", n), val("kappa", n)
                    note("lambda_to_0", abs(lm * K[n] - g0), lm * (abs(kap) + env + ue + uk + U[n]) + ug)
                else:
                    per = entry["lengths"][fhex(lam)]
                    k13 = K["term_13a"] if n in ("term_13a", "no_hole", "total") else 0
                    htr = lm * (K[n] - k13)
                    (gi, ugi), (a1, ua1) = val("G_inf", n), val("A1", n)
                    etr, utr = val("E_tr", n, per)
                    rem, urem = val("remainder_bound_after_first_order", n, per)
                    slack = lm * (U[n] + U["term_13a"]) + ugi
                    note("lambda_to_inf_envelope", abs(htr - gi), etr + utr + slack)
                    note("lambda_to_inf_first_order", abs(htr - gi - a1 / lm), rem + urem + ua1 / lm + slack)

    rep = max(results, key=lambda g: g["seconds"])
    doc = {
        "schema": "xpsloss.m3a.references.v1",
        "generated_by": "benchmarks/numerics/m3a/generate_references_m3a.py",
        "prereg": "docs/M3A_PREREG.md section 5 (R1)",
        "table": {"file": "benchmarks/numerics/m3a/expected_states_m3a.json", "sha256": TABLE_SHA256},
        "units": {"K_sc": "A^-1 eV^-1 per solid path", "G, H": "eV^-1", "loss": "eV", "length": "A", "momentum": "A^-1"},
        "constants": {"H_eV_A2": {"decimal": H_DEC, "binary64_hex": b64(H_DEC).hex()},
                      "a0_A": {"decimal": A0_DEC, "binary64_hex": b64(A0_DEC).hex()}},
        "outputs": {"terms": list(TERMS), "electron_group": [DEPTH_TERMS[i] for i in ELECTRON_IDX],
                    "hole_group": [DEPTH_TERMS[i] for i in HOLE_IDX],
                    "no_hole": "term_13a + electron_group", "hole": "hole_group (unit hole charge)",
                    "total": "no_hole + Q_h*hole, u = u(no_hole) + |Q_h| u(hole); formed by the test in rationals",
                    "below_gap": "term_13a = 0 and total(Q_h = 1) = 0 analytically (u = 0); total(Q_h) = (1-Q_h)*electron_group"},
        "method": {
            "library": "mpmath", "digits_run_A": DPS_A, "digits_run_B": DPS_B,
            "response": "Im F = L (E1) with the hard gap; g = 0: E4b; g > 0: Re F by E4a, partial fractions of the retained band",
            "equations": "Simonsen et al. (1997) Eqs. (13)-(18), theta = 0, typed from docs/M1_EQUATION_MAP.md 4.1",
            "route_a": "closed-form depth average derived by hand (compact_laplace); adaptive composite 24-point Gauss-Legendre in ln k_parallel, run A 50 digits, run B 80 digits with every panel halved; value = run B",
            "route_b": "depth quadrature of the printed fixed-depth factors (the elementary depth factors 1, u, |u|^2 and (e^{-iba}-e^{-ka})(e^{iba}-e^{-ka}) are integrated numerically and the depth-independent factors of the printed braces are applied afterwards; the identity is checked pointwise in self_checks; the depth quadrature uses fixed-point integers with 256 fractional bits and the addition theorems of the exponentials inside a region of equal panels, checked against a direct mpf evaluation in self_checks; and the whole printed pieces are integrated for the cases of route_b_whole_pieces_check): composite Gauss-Legendre on panels cut by the decay rates p+2k, p+k, p and the phase b*a, up to exp(-52); on the accepted parent panels of run A, 50 digits",
            "route_b_variants": {"base": "panel size 20/(rate+b), 24 points, tail 52/p",
                                 "order": "panel size 40/(rate+b), 48 points",
                                 "phase": "panel size 26/(rate+b), 24 points (another node set)",
                                 "tail": "base plus [52/p, 64/p]"},
            "v1_depth_quadrature": "t = p*a panels of size 10/(rate+b), 24 points, tail 64/p; changes: 48 points with size 20, size 7, tail 80/p; digits raised by 2*log10(p/min(k,b)) for the cancellations",
            "groups": "full printed braces with the hole factor kept apart, transformed mode by mode (not a sum of the terms, not route a)",
            "u_term": "|a-b| + 2*|A-B| + residual(A) + max(order, phase, tail changes) + decimal rounding",
            "u_term_13a": "2*|A-B| + residual(A) + decimal rounding",
            "u_group": "|group - sum of route-b terms| + 2*|A-B| + residual(A) + sum of route-b changes of the members + decimal rounding",
            "u_no_hole": "u(term_13a) + u(electron_group)", "u_hole": "u(hole_group)",
            "panel_tolerance": "|parent - children| <= 1e-24*|coarse integral| + 1e-30*max column (envelopes: 1e-13)",
            "term_13a_independent": "V3: (a) L(0,T) log(q+/q-)/(pi a0 E0) without dispersion; (b) bulk integral over the total momentum for k_total with dispersion",
            "limits": "mode coefficients c0, cs, cd of G_j(a) from the printed factors (design 2.2); G0 = Int(c0+Re cs+cd), envelope = Int(|s cs|+2k|cd|), G_inf = Int c0, A1 = Int(Re(cs/s)+cd/(2k)), E_tr = Int(|cs|/|1+lam s|+|cd|/(1+2 lam k)), remainder bound = Int(|cs|/(lam|s||1+lam s|)+|cd|/(2 lam k(1+2 lam k))); u = 2*|A-B| + residual(A) + decimal rounding; aggregated outputs use the aggregated complex coefficients in the envelopes",
            "storage": "values: decimal strings of 40 significant digits (envelopes 20); u: decimal strings rounded up; the rounding of the stored decimal is inside u"},
        "self_checks": dict(formula_self_check(), kk_cross_check=kk_cross_check(),
                            route_b_whole_pieces_check=list(blackbox)),
        "v1": v1_block(table),
        "references": references,
        "groups": groups,
        "limits": limits_out,
        "summary": {"reference_keys": len(references), "limit_keys": len(limits_out), "groups": len(groups),
                    "max_u_over_Bref": float(worst[0]), "max_u_over_Bref_at": worst[1],
                    "undetermined": undetermined,
                    "internal_limit_checks": limit_checks,
                    "max_abs_numerical_residual_of_zero": float(max_resid_zero),
                    "max_abs_13a_difference_from_independent": max(
                        (float(e["term_13a_independent"]["abs_difference_from_term_13a"]) for e in groups.values()
                         if "term_13a_independent" in e), default=0.0)},
        "environment": {"platform": platform.platform(), "python": sys.version.split()[0],
                        "mpmath": mp.__version__, "mpmath_backend": mp.libmp.BACKEND, "processes": processes},
        "timing_seconds": {"wall_compute": round(t_wall, 1),
                           "cpu_sum": round(sum(g["seconds"] for g in results), 1),
                           "cpu_route_a": round(sum(g["seconds_route_a"] for g in results), 1),
                           "cpu_route_b": round(sum(g["seconds_route_b"] for g in results), 1),
                           "slowest_group": {"key": rep["key"], "seconds": round(rep["seconds"], 1)}},
    }
    return doc


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--processes", type=int, default=max(1, min(12, (os.cpu_count() or 2))))
    parser.add_argument("--cache", required=True, help="directory of the per-group results (resumable)")
    parser.add_argument("--stage", choices=("compute", "assemble", "all"), default="all")
    parser.add_argument("--output", default=str(HERE / "references_m3a.json"))
    parser.add_argument("--time-one", type=int, default=None, help="time one group (index) without the cache")
    parser.add_argument("--only", default=None, help="comma separated group indices (development)")
    parser.add_argument("--partial", action="store_true",
                        help="assemble only the cached groups (development; not a frozen reference)")
    args = parser.parse_args()
    table = load_table()
    groups = build_groups(table)
    os.makedirs(args.cache, exist_ok=True)
    if args.time_one is not None:
        spec = dict(groups[args.time_one])
        res = group_task(spec)
        print(json.dumps({k: res[k] for k in ("key", "lengths", "panels_a", "nodes_a", "route_b_k_nodes",
                                              "route_b_depth_nodes_base", "seconds_route_a", "seconds_route_b",
                                              "seconds_extra", "seconds")}, indent=1))
        return
    for i, g in enumerate(groups):
        g["cache"] = os.path.join(args.cache, f"group_{i:03d}.json")
    bb_specs = []
    for n, (material, convention, loss, length) in enumerate(BLACKBOX):
        g = next(g for g in groups if g["material"] == material and g["convention"] == convention and g["loss"] == loss)
        if length not in g["lengths"]:
            raise SystemExit("black-box case is not a reference key")
        bb_specs.append({**{k: v for k, v in g.items() if k != "limits"}, "lengths": [length], "blackbox": True,
                         "cache": os.path.join(args.cache, f"blackbox_{n}.json")})
    wall_file = os.path.join(args.cache, "wall_seconds.txt")
    if args.stage in ("compute", "all"):
        todo = [g for g in groups if not os.path.exists(g["cache"])]
        if args.only:
            wanted = {int(v) for v in args.only.split(",")}
            todo = [g for i, g in enumerate(groups) if i in wanted and not os.path.exists(g["cache"])]
        # long groups first (extreme lengths and narrow lines)
        todo.sort(key=lambda g: -len(g["lengths"]))
        started = time.perf_counter()
        if not args.only:
            todo += [g for g in bb_specs if not os.path.exists(g["cache"])]
        with Pool(args.processes) as pool:
            for n, res in enumerate(pool.imap_unordered(any_task, todo, chunksize=1)):
                print(f"{n + 1}/{len(todo)} {res['key']} {res['seconds']:.1f}s", flush=True)
        previous = float(open(wall_file).read()) if os.path.exists(wall_file) else 0.0
        with open(wall_file, "w") as handle:
            handle.write(repr(previous + time.perf_counter() - started))
    if args.stage in ("assemble", "all"):
        if args.partial:
            groups = [g for g in groups if os.path.exists(g["cache"])]
            have = {g["key"] for g in groups}
            gk = lambda r: group_key(r["material"], r["electron_energy_eV"], r["epsilon_momentum"], r["loss_eV"])
            table["reference_keys"] = [r for r in table["reference_keys"] if gk(r) in have]
            table["limit_keys"] = [r for r in table["limit_keys"] if gk(r) in have]
        results = [group_task(g) for g in groups]
        t_wall = float(open(wall_file).read()) if os.path.exists(wall_file) else 0.0
        bb = [blackbox_task(g) for g in bb_specs if os.path.exists(g["cache"]) or not args.partial]
        doc = assemble(table, results, t_wall, args.processes, bb)
        if args.partial:
            doc["partial"] = True
        text = json.dumps(doc, indent=1, ensure_ascii=True) + "\n"
        Path(args.output).write_bytes(text.encode("ascii"))
        print("wrote", args.output)
        print("sha256", hashlib.sha256(text.encode("ascii")).hexdigest())
        print("summary", json.dumps({k: v for k, v in doc["summary"].items()}))
        print("timing", json.dumps(doc["timing_seconds"]))


if __name__ == "__main__":
    main()
