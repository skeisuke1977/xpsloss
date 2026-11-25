"""Independent depth paths for M3a (docs/M3A_PREREG.md section 5, R2; docs/M3A_DESIGN.md section 3.2).

Written by R2 without the depth-averaged implementation.  Two paths average the
fixed-depth cross section over the depth weight of Eq. (18) of Simonsen, Yubero and
Tougaard, Phys. Rev. B 56, 1612 (1997):

    K_sc,j(T; lambda) = Int_0^inf w(x) K_j(T; a=x) dx,   w(x) = x exp(-x/lambda) / lambda^2,

at normal emission (x = a).  With t = a/lambda, w(x) dx = t exp(-t) dt.

Path 1 (V8): the values come only from the public fixed-depth function
xpsloss.effective_cross_section_fixed_depth.  One call per depth node evaluates every
loss of the registered set P (material, E0, momentum, lambda).  The public dielectric
function (xpsloss.negative_inverse_dielectric_function) is used only for the tail
envelope C_j, which enters multiplied by exp(-r).

Path 3 (V9): the frozen tests/xps_route_b_reference.py (k_perp kernel of Eq. (2), the
I3 route of M2) at the same depth nodes as the base rule of path 1, k_parallel only.

Depth rule (prereg section 5, written here as it is):
  * t in [0, r], r = 40.  Edges: the registered 0, 0.01, 0.04, 0.12, 0.3, 0.7, 1.4, 2.5,
    4.5, 8, 14, 25, 40, plus (i) geometric edges r 2^-m (m = 1, 2, ...) while
    r 2^-m >= 1/(8 k_par+,max lambda), k_par+,max the largest Eq. (17) upper limit of P,
    and (ii) an equal split of each interval until the phase increment
    b_eff(t) lambda dt <= pi/4.  b_eff(t) is the largest b(T) = T/(hbar v) over the losses
    of P with t < t_c(T) = 40/(k_par-(T) lambda) (no split if there is none).  b_eff does
    not increase with t, so an interval [t0, t1] is cut into the smallest number n of equal
    pieces with b_eff(t0) lambda (t1 - t0)/n <= pi/4.  Gauss-Legendre with 10 nodes on each
    piece.  Node weight alpha_i = (GL weight) * t_i exp(-t_i).
  * Variations, one at a time: (A) 16 nodes per piece; (B) phase increment <= pi/8 and
    geometric ratio sqrt(2) (edges r 2^(-m/2)); (C) r = 50, read as the base pieces on
    [0, 40] plus the interval [40, 50] cut by the same phase rule (geometric edges stay
    those of r = 40).
  * Tail t > r: |kappa_j| (1 + r) e^-r + C_j e^-r / lambda, C_j = Int (|c0j| + |csj| + |cdj|)
    d ln k (upper side of its quadrature), from the mode coefficients of design section 2.2.
  * u_j = max(|Delta_A|, |Delta_B|, |Delta_C|) + sum_i alpha_i e_ij + rounding + tail.
  * Allocation: alpha_i e_ij <= f_i eta D_j with eta = 1/9 and f_i = 1/N (N base nodes).
    Gapless nodes publish one common quadrature_error, given to every column.  Each hole
    charge of P is its own fixed-depth call (hole_charge = Q_h), so the total of every
    charge carries the estimate the public function published for it; no total is rebuilt
    from unit-charge components and the small direct-sum error of a sub-gap unit-charge
    total is never reused for another charge.
"""
from dataclasses import dataclass, field
import math
from pathlib import Path
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss

from xpsloss import (
    FixedDepthConvergenceError, Material, NumericalConvergenceError, Oscillator, effective_cross_section_fixed_depth,
    negative_inverse_dielectric_function,
)
from xpsloss.model import BOHR_RADIUS_A, HBAR2_OVER_2ME_EV_A2

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
H = HBAR2_OVER_2ME_EV_A2
A0 = BOHR_RADIUS_A
EPS = np.finfo(float).eps

TERMS = ("term_13a", "term_13b", "term_14", "term_15_T1", "term_15_T2", "term_15_T3", "term_16_e", "term_16_h")
NO_HOLE_TERMS = ("term_13a", "term_13b", "term_15_T1", "term_15_T3", "term_16_e")
HOLE_TERMS = ("term_14", "term_15_T2", "term_16_h")
OUTPUTS = TERMS + ("no_hole", "hole", "total")
UNIT = OUTPUTS[:10]  # charge-free columns; total is formed per hole charge

REGISTERED_EDGES = (0.0, 0.01, 0.04, 0.12, 0.3, 0.7, 1.4, 2.5, 4.5, 8.0, 14.0, 25.0, 40.0)
BASE_R = 40.0
DECAY = 40.0          # t_c(T) = DECAY / (k_par-(T) lambda): exp(-k a) < e^-40 beyond
ETA = 1.0 / 9.0
V8_REL, V8_ABS = 1e-7, 1e-8
V9_REL, V9_ABS = 1e-6, 1e-6
ALLOCATION_SAFETY = 0.5   # requested node accuracy is this fraction of the allotted share
PRELIM_SETTINGS = (1e-9, 1e-7, 200)
NODE_EPSREL = 1e-15       # effectively absolute budgets: max(epsabs_i, 1e-15 |K|)
NODE_LIMIT = 200

VARIANTS_PATH1 = {
    "base": dict(order=10, phase=math.pi / 4, ratio=2.0, upper=None),
    "A": dict(order=16, phase=math.pi / 4, ratio=2.0, upper=None),
    "B": dict(order=10, phase=math.pi / 8, ratio=math.sqrt(2.0), upper=None),
    "C": dict(order=10, phase=math.pi / 4, ratio=2.0, upper=50.0),
}
VARIANTS_PATH3 = ("base", "A", "B", "C", "D", "E")


# ----------------------------------------------------------------------------- kinematics (layout and envelope only)

def hbar_v(E0):
    """hbar v in eV A for the non-relativistic E0 = m v^2 / 2."""
    return 2.0 * math.sqrt(H * E0)


def kpar_limits(E0, T):
    """Eq. (17) at theta = 0 without cancellation:
    k-^2 = (s-r)^3 (3s+r) / (4 s^2 H), k+^2 = (s+r)^3 (3s-r) / (4 s^2 H), s = sqrt(E0), r = sqrt(E0-T)."""
    T = np.asarray(T, dtype=float)
    s, r = math.sqrt(E0), np.sqrt(E0 - T)
    d = T / (s + r)
    lo = np.sqrt(d ** 3 * (3 * s + r) / (4 * s * s * H))
    hi = np.sqrt((s + r) ** 3 * (3 * s - r) / (4 * s * s * H))
    return lo, hi


# ----------------------------------------------------------------------------- materials of the table

def material_from_table(name, table):
    entry = table["materials"][name]
    return Material(name, float(entry["gap_eV"]), tuple(Oscillator(*map(float, line)) for line in entry["oscillators"]))


# ----------------------------------------------------------------------------- depth rule

@dataclass
class DepthRule:
    lam: float
    order: int
    phase: float
    ratio: float
    r: float
    upper: float
    edges: list
    pieces: list
    t: np.ndarray = field(repr=False)
    gl_weight: np.ndarray = field(repr=False)
    alpha: np.ndarray = field(repr=False)

    @property
    def n(self):
        return int(self.t.size)


def _geometric_edges(r, ratio, floor):
    """r ratio^-m for m = 1, 2, ... while >= floor (ratio 2: exact halving; sqrt 2: r 2^(-m/2))."""
    out = []
    m = 1
    while True:
        if ratio == 2.0:
            e = math.ldexp(r, -m)
        elif ratio == math.sqrt(2.0):
            e = r * 2.0 ** (-m / 2.0)
        else:
            e = r * ratio ** (-m)
        if e < floor:
            return out
        if 0.0 < e < r:
            out.append(e)
        m += 1


def depth_rule(losses_eV, E0, lam, *, order=10, phase=math.pi / 4, ratio=2.0, r=BASE_R, upper=None,
               extra_edges=()):
    """Pieces and nodes of the path-1 depth quadrature (see the module docstring).

    upper: (C) adds the interval [r, upper], split by the phase rule only.
    extra_edges: test-only additional edges (V2 (ii)).
    """
    losses = np.atleast_1d(np.asarray(losses_eV, dtype=float))
    k_lo, k_hi = kpar_limits(E0, losses)
    b = losses / hbar_v(E0)
    t_c = DECAY / (k_lo * lam)
    edges = {e for e in REGISTERED_EDGES if e <= r}
    edges.add(0.0)
    edges.add(float(r))
    edges.update(_geometric_edges(float(r), ratio, 1.0 / (8.0 * float(k_hi.max()) * lam)))
    edges.update(float(e) for e in extra_edges if 0.0 < e < r)
    if upper is not None:
        if not upper > r:
            raise ValueError("upper must exceed r")
        edges.add(float(upper))
    edges = sorted(edges)
    pieces = []
    for lo, hi in zip(edges, edges[1:]):
        active = b[lo < t_c]
        beff = float(active.max()) if active.size else 0.0
        n = 1
        if beff > 0.0:
            n = max(1, math.ceil(beff * lam * (hi - lo) / phase))
            while n > 1 and beff * lam * ((hi - lo) / (n - 1)) <= phase:
                n -= 1
            while beff * lam * ((hi - lo) / n) > phase:
                n += 1
        for i in range(n):
            a = lo + (hi - lo) * i / n
            c = hi if i == n - 1 else lo + (hi - lo) * (i + 1) / n
            pieces.append((a, c))
    x, w = leggauss(order)
    lo = np.array([p[0] for p in pieces])
    hi = np.array([p[1] for p in pieces])
    t = (0.5 * (hi - lo)[:, None] * x + 0.5 * (hi + lo)[:, None]).ravel()
    gw = (0.5 * (hi - lo)[:, None] * w).ravel()
    alpha = gw * t * np.exp(-t)
    return DepthRule(float(lam), order, phase, ratio, float(r), upper, edges, pieces, t, gw, alpha)


def variant_rule(losses_eV, E0, lam, variant):
    """Node set of a path-1 variant.  For C only the nodes on [40, 50] are new."""
    spec = VARIANTS_PATH1["base" if variant in ("D", "E") else variant]
    return depth_rule(losses_eV, E0, lam, order=spec["order"], phase=spec["phase"], ratio=spec["ratio"],
                      upper=spec["upper"])


def extension_mask(rule):
    """Nodes of rule C beyond the base upper limit r = 40."""
    return rule.t > BASE_R


# ----------------------------------------------------------------------------- path 1 nodes

def _pack(result, gapped):
    values = np.array([np.asarray(getattr(result, name), dtype=float) for name in OUTPUTS])
    if gapped:
        errors = np.array([np.asarray(result.component_errors[name], dtype=float) for name in OUTPUTS])
    else:
        # The gapless contract publishes one common estimate that covers all eleven outputs.
        errors = np.broadcast_to(np.asarray(result.quadrature_error, dtype=float), values.shape).copy()
    return values, errors


def fixed_depth_node(material, E0, losses, depth, spec, epsabs, epsrel=NODE_EPSREL, limit=NODE_LIMIT, *,
                     hole_charge=1.0):
    """One depth node of path 1: values and errors (11, n) for hole_charge, and a status per loss.

    The whole loss array goes to one call.  If that call fails, each loss is replayed alone so
    that the failure is attributed to its loss; a failing loss is a failure of the path there.
    """
    losses = np.atleast_1d(np.asarray(losses, dtype=float))
    n = losses.size
    gapped = material.gap_eV > 0
    values = np.full((11, n), np.nan)
    errors = np.full((11, n), np.nan)
    status = [""] * n

    def call(grid):
        return effective_cross_section_fixed_depth(material, E0, grid, depth, theta_deg=0.0, hole_charge=hole_charge,
                                                   epsilon_momentum=spec, epsabs=epsabs, epsrel=epsrel, limit=limit)
    try:
        values[:], errors[:] = _pack(call(losses), gapped)
        return values, errors, status
    except (NumericalConvergenceError, ValueError, ArithmeticError) as exc:
        if n == 1:
            status[0] = f"{type(exc).__name__}: {exc}"
            return values, errors, status
    for i in range(n):
        try:
            v, e = _pack(call(losses[i:i + 1]), gapped)
            values[:, i], errors[:, i] = v[:, 0], e[:, 0]
        except (NumericalConvergenceError, ValueError, ArithmeticError) as exc:
            status[i] = f"{type(exc).__name__}: {exc}"
    return values, errors, status


def node_epsabs(rule, d_floor, charge_factor, n_alloc=None, safety=ALLOCATION_SAFETY):
    """Requested absolute accuracy per node: safety * f_i eta D_floor / (alpha_i * charge_factor)."""
    n_alloc = rule.n if n_alloc is None else n_alloc
    with np.errstate(divide="ignore"):
        eps = safety * ETA * (1.0 / n_alloc) * d_floor / (rule.alpha * charge_factor)
    return np.minimum(eps, 1.0)


def charge_factor(charges):
    return max(1.0 if q == 1.0 else 1.0 + abs(q) for q in charges)


# ----------------------------------------------------------------------------- assembly of path values

def weighted_sum(alpha, values):
    """sum_i alpha_i values[i, ...] with math.fsum, and its rounding bound."""
    flat = values.reshape(values.shape[0], -1)
    prod = alpha[:, None] * flat
    out = np.array([math.fsum(col) for col in prod.T]).reshape(values.shape[1:])
    mag = np.array([math.fsum(col) for col in np.abs(prod).T]).reshape(values.shape[1:])
    # alpha_i (GL weight, t, exp) ~ 4 roundings, product 1, fsum exact then 1.
    return out, 8 * EPS * mag + EPS * np.abs(out)


# ----------------------------------------------------------------------------- tail envelope (design section 2.2)

def _resonance_edges(material, T, spec, E0, lo, hi):
    """ln k edges near the bulk and single-oscillator surface lines of dispersive inputs (layout only)."""
    b = T / hbar_v(E0)
    out = []
    for o in material.oscillators:
        if o.dispersion <= 0 or o.strength_eV2 == 0:
            continue
        for target in (T, math.sqrt(T * T + o.strength_eV2 / 2.0)):
            if target <= o.energy_eV:
                continue
            q2 = (target - o.energy_eV) / (o.dispersion * H)
            k2 = q2 - b * b if spec == "k_total" else q2
            if k2 <= 0:
                continue
            kr = math.sqrt(k2)
            hw = o.width_eV / (4.0 * o.dispersion * H * kr)
            for f in (0.0, 0.25, 0.5, 1, 2, 4, 8, 16, 32):
                for sgn in (1, -1):
                    kk = kr + sgn * f * hw
                    if lo < kk < hi:
                        out.append(math.log(kk))
    return out


def mode_coefficients(material, E0, T, k, spec):
    """kappa density and the mode coefficients c0, cs (complex), cd of G_j for the eleven outputs
    (unit hole charge; total of charge Q is no_hole + Q hole) in d ln k, from the fixed-depth
    integrand with u = exp(-s a), 1 - u, |1-u|^2 = 1 - 2 Re u + e^-2ka and u - e^-2ka:

      13b:  (M/2b)   [c0 = 2 ImS Im(1/s^2) - ReS/D,  cs = 2i ImS/s^2 + 2 ReS/D,  cd = -ReS/D]
      14:   (M/2b^2) [c0 = Re(-iS/s),  cs = iS(1/s + 1/conj s),  cd = Re(-iS/conj s)]
      15T1: (M/2b)   [ReC/D, -2 ReC/D, ReC/D]
      15T2: (M/2b^2) [0, -iC/conj s, Re(iC/conj s)]
      15T3: (M/2b)   [Re(-R/conj s^2), conj(R)/s^2, 0]
      16e:  (M/2b)   [Re(2kR/(sD)), -R/s^2, 0]
      16h:  (M/2b^2) [0, -iR/s, 0]
      13a:  kappa density (M/D) ImS (G = 0)
    with s = k + ib, D = |s|^2, M = k^2/(pi a0 E0), S = 1 + F, R = S/(2-S), C = 2S - R."""
    k = np.asarray(k, dtype=float)
    b = T / hbar_v(E0)
    q = np.hypot(k, b) if spec == "k_total" else k
    F = np.asarray(negative_inverse_dielectric_function(material, np.full(k.shape, T), q), dtype=complex)
    S = 1.0 + F
    R = S / (2.0 - S)
    C = 2.0 * S - R
    s = k + 1j * b
    sc = np.conj(s)
    D = k * k + b * b
    M = k * k / (math.pi * A0 * E0)
    p1, p2 = M / (2 * b), M / (2 * b * b)
    zero = np.zeros_like(k)
    c0 = {
        "term_13b": p1 * (2 * S.imag * (1 / s ** 2).imag - S.real / D),
        "term_14": p2 * (-1j * S / s).real,
        "term_15_T1": p1 * C.real / D,
        "term_15_T2": zero,
        "term_15_T3": p1 * (-R / sc ** 2).real,
        "term_16_e": p1 * (2 * k * R / (s * D)).real,
        "term_16_h": zero,
    }
    cs = {
        "term_13b": p1 * (2j * S.imag / s ** 2 + 2 * S.real / D),
        "term_14": p2 * 1j * S * (1 / s + 1 / sc),
        "term_15_T1": p1 * (-2 * C.real / D) + 0j,
        "term_15_T2": p2 * (-1j * C / sc),
        "term_15_T3": p1 * np.conj(R) / s ** 2,
        "term_16_e": -p1 * R / s ** 2,
        "term_16_h": -p2 * 1j * R / s,
    }
    cd = {
        "term_13b": -p1 * S.real / D,
        "term_14": p2 * (-1j * S / sc).real,
        "term_15_T1": p1 * C.real / D,
        "term_15_T2": p2 * (1j * C / sc).real,
        "term_15_T3": zero,
        "term_16_e": zero,
        "term_16_h": zero,
    }
    kappa = M / D * S.imag
    out = {}
    for name in TERMS[1:]:
        out[name] = (c0[name], cs[name], cd[name])
    out["term_13a"] = (zero, zero + 0j, zero)
    for group, members in (("no_hole", NO_HOLE_TERMS[1:]), ("hole", HOLE_TERMS)):
        out[group] = tuple(sum(out[m][i] for m in members) for i in range(3))
    return kappa, out


def _ln_k_nodes(material, E0, T, spec, order, span=0.1):
    lo, hi = kpar_limits(E0, T)
    lo, hi = float(lo), float(hi)
    u0, u1 = math.log(lo), math.log(hi)
    edges = list(np.linspace(u0, u1, max(2, math.ceil((u1 - u0) / span) + 1)))
    edges += _resonance_edges(material, T, spec, E0, lo, hi)
    edges = np.unique(np.asarray(edges))
    x, w = leggauss(order)
    a, c = edges[:-1], edges[1:]
    u = (0.5 * (c - a)[:, None] * x + 0.5 * (c + a)[:, None]).ravel()
    wu = (0.5 * (c - a)[:, None] * w).ravel()
    return np.exp(u), wu


def envelope(material, E0, T, spec, charges=(0.0, 0.5, 1.0)):
    """Upper bounds |kappa_j| and C_j = Int (|c0| + |cs| + |cd|) d ln k for each output and charge.

    Two Gauss-Legendre orders (16 and 8) on the same ln k panels; the bound is the 16-node value
    plus the difference of the two (the upper side of the quadrature).  Returns None if a bound
    is not finite (the path comparison is then held)."""
    res = {}
    for order in (16, 8):
        k, w = _ln_k_nodes(material, E0, T, spec, order)
        kappa, coeff = mode_coefficients(material, E0, T, k, spec)
        row = {"kappa": float(np.sum(w * kappa))}
        for name in UNIT:
            c0, cs, cd = coeff[name]
            row["C:" + name] = float(np.sum(w * (np.abs(c0) + np.abs(cs) + np.abs(cd))))
        for qh in charges:
            c0 = coeff["no_hole"][0] + qh * coeff["hole"][0]
            cs = coeff["no_hole"][1] + qh * coeff["hole"][1]
            cd = coeff["no_hole"][2] + qh * coeff["hole"][2]
            row[f"C:total:{qh!r}"] = float(np.sum(w * (np.abs(c0) + np.abs(cs) + np.abs(cd))))
        res[order] = row
    out = {}
    for key, v16 in res[16].items():
        v8 = res[8][key]
        bound = abs(v16) + abs(v16 - v8)
        if not math.isfinite(bound):
            return None
        out[key] = bound
    return out


def signed_G(material, E0, T, spec, a, order=16):
    """G_j(a) = Int [c0 + Re(cs e^{-sa}) + cd e^{-2ka}] d ln k (self-check of the mode coefficients)."""
    k, w = _ln_k_nodes(material, E0, T, spec, order)
    kappa, coeff = mode_coefficients(material, E0, T, k, spec)
    b = T / hbar_v(E0)
    u = np.exp(-(k + 1j * b) * a)
    out = {"kappa": float(np.sum(w * kappa))}
    for name in UNIT:
        c0, cs, cd = coeff[name]
        out[name] = float(np.sum(w * (c0 + (cs * u).real + cd * np.exp(-2 * k * a))))
    return out


def tail_bound(env, column, lam, r=BASE_R, charge=None):
    """|kappa_j|(1+r)e^-r + C_j e^-r / lambda for output column (with hole charge for total)."""
    has_kappa = column in ("term_13a", "no_hole", "total")
    key = f"C:total:{charge!r}" if column == "total" else "C:" + column
    kap = env["kappa"] if has_kappa else 0.0
    return kap * (1 + r) * math.exp(-r) + env[key] * math.exp(-r) / lam


# ----------------------------------------------------------------------------- path 3 (I3 type)

def import_route_b():
    if str(HERE) not in sys.path:
        sys.path.insert(0, str(HERE))
    import xps_route_b_reference as ref  # noqa: E402  (frozen; never edited)
    return ref


def configure_route_b(ref, variant):
    """(D) kpar_nodes with n_base 24; (E) k_perp core and tail node arrays twice as fine.
    Call once per process: the frozen module is changed in memory only."""
    if variant == "D":
        import functools
        original = ref.kpar_nodes
        ref.kpar_nodes = functools.partial(original, n_base=24)
    elif variant == "E":
        from numpy.polynomial.laguerre import laggauss
        ref._XG, ref._WG = leggauss(40)
        ref._XL, ref._WL = laggauss(160)
        xs, ws = leggauss(64)
        ref._XS, ref._WS = 0.5 * (xs + 1.0), 0.5 * ws


def i3_node(ref, materials, E0, losses, depth):
    """Ten unit-charge columns (8 terms, no_hole, hole) per material at one depth: shape (n_mat, 10, n)."""
    res = ref.reference_terms(materials, E0, losses, depth)
    return np.array([[np.asarray(r[name], dtype=float) for name in UNIT] for r in res])


# ----------------------------------------------------------------------------- V2: the quadrature on test functions

def v2_functions(lam, kappa=2.0, intercept=3.0, s=(0.1 + 0.05j, 1 + 0.5j)):
    """Test functions K(x) and their exact averages over w(x) on [0, inf)."""
    fs = {
        "one": (lambda x: np.ones_like(x), 1.0, 1.0),
        "x": (lambda x: x, 2 * lam, 1.0 / lam),
        "inv_x": (lambda x: 1.0 / x, 1.0 / lam, lam),
        "kappa_plus_I_over_x": (lambda x: kappa + intercept / x, kappa + intercept / lam, 1.0),
    }
    for z in s:
        fs[f"exp_sx_over_x({z!r})"] = (lambda x, z=z: np.exp(-z * x) / x, 1.0 / (lam * (1 + lam * z)), lam)
    return fs


def quadrature_average(rule, func):
    """Sum alpha_i K(lambda t_i) for a test function (complex allowed)."""
    x = rule.lam * rule.t
    vals = func(x)
    return complex(math.fsum((rule.alpha * np.real(vals)).tolist()), math.fsum((rule.alpha * np.imag(vals)).tolist()))
