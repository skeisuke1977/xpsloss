"""Independent references R-B for M4b-1 (docs/M4B1_PREREG.md 6.2, 6.3, 7.2; design 4.1, 10, 14 T11).

Written by the reference author R-B, separately from the implementation.  The only package function
that is called is the public M3a kernel xpsloss.effective_cross_section_depth_averaged (API default
accuracy 1e-12, 1e-9, 200); the public input types xpsloss.Material and xpsloss.Oscillator build its
argument.  No module of the M4a or M4b implementation (xps_single_loss, _loss_panels, _fold,
xps_multiple_loss, _loss_orders) is imported: no panels, no panel convolution, no cells.  What is
shared with the implementation is the upstream point value K_j(T) itself; an error of that point
value is outside the scope of this comparison.

Usage (point values are computed into per-configuration cache files, then assembled):

    python benchmarks/numerics/m4b1/generate_references_rb_m4b1.py --cache DIR --points Al_k_total_Q1
    python benchmarks/numerics/m4b1/generate_references_rb_m4b1.py --cache DIR --assemble
    python benchmarks/numerics/m4b1/generate_references_rb_m4b1.py --cache DIR --assemble --check

The JSON contains no time stamps and no timings; a new run on the same environment gives the same
bytes (the cache holds the binary64 point values and errors, bit for bit).  Timings are written to
DIR/timing_*.json and recorded in the appendix of docs/M4B1_PREREG.md.

Model B1 first_loss_13a_chain_v0 (design 4.1, prereg 4).  K_j(T) is the M3a output at the binary64
loss T (eV^-1 A^-1; eleven columns; the hole columns are for a unit hole; `total` is the M3a column
itself, never rebuilt).  With the window W = [a, b], the multiplier m (A) and R = K_term_13a,

    C_{2,j}(s) = m^2 Int K_j(t) R(s-t) dt,   C_{3,j} = m^3 K_j * (R * R),   supp C_n = [n a, n b]
    order[n].<region>.<c>      J = Int dE' F(E') Int_{E'-E in region} C_{n,c}(E'-E) dE
    order[n].area.<column>     w_{n,j} = m^n A_j A_13a^(n-1),  A_j = Int_a^b K_j dT
    order[n].cutoff_sensitivity.<column>
                               h_1 = m a K_j(a),  h_n = h_{1,j} r^(n-1) + (n-1) w_{1,j} h_r r^(n-2)
    order_area_ratio           r = m A_13a

Numerical route (prereg 6.2).  Uniform grids h = 2^-k eV, k = 4..8, T_i = a + i h (exact binary64;
a, b, the bin edges and the primary edges are grid points for k >= 1).  At every level the
convolution is the trapezoid rule on the exact overlap interval [max(a, s-b'), min(b, s-a')]:
full discrete convolution minus half of the two end terms.  This treats the end sums separately:
at s = 2a and 2b the overlap is one point and the value is exactly 0, at s = a+b the two end terms
are K(a)R(b) and K(b)R(a); for order 3 the same holds at 3a, 2a+b, a+2b, 3b, and the kink of R*R at
a+b falls on a grid node of the inner variable (full weight).  The region integral is the trapezoid
rule with a positive weight vector (difference of the cumulative distribution written as one local
sum, so that no large cumulative sums are subtracted); for a binned primary the outer integral over
E' is the trapezoid rule on the same grid (integral of the cumulative distribution).  All these
discrete values have an expansion in h^2, h^4, ... with break points only at grid nodes, and are
Richardson-extrapolated three times from the levels k = 5, 6, 7, 8 (the value r).  The window
areas are extrapolated separately over the whole window and are not taken from the convolutions.

Uncertainty u_r (every part stored per key, then added; prereg 6.2):
  (1) finest_two     |T3(levels 5..8) - T3(levels 4..7)|: the extrapolations ending at the two finest h
  (2) extrapolations |T3(levels 5..8) - T2(levels 6..8)|: one extrapolation less
  (3) upstream       the M3a per-column errors (component_errors) carried through the absolute values
                     of the convolution and region weights, on the level-6 grid with the errors and
                     |K| replaced by their maximum over the neighbouring finest nodes, times the sum
                     of the absolute Richardson coefficients
  (4) primary_bins   binned primary only: |main - alternative|; the alternative integrates over E'
                     by a 16-interval Romberg rule (h/32 = 1/32 eV nodes) of the level values
  (5) rounding       (8 M + 64) 2^-53 times the absolute sums (same carrier as (3)), times the sum of
                     the absolute Richardson coefficients, plus the 40-digit decimal
The rule "no reference" (prereg 6.3: f_r = u_r/|r| > 1.25e-5 for r != 0) is applied here, before
any implementation value exists.  A region that no order-n loss can reach is an exact zero
(value 0, u_r 0); this is decided with exact rationals.
"""
import os

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_name, "1")

import argparse
from fractions import Fraction
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy

import xpsloss
from xpsloss import Material, Oscillator, effective_cross_section_depth_averaged

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TABLE_FILE = HERE / "expected_states_m4b1.json"
TABLE_SHA256 = "d8d0f892a49c3ed9933b08d6a7c92e5b3aaf05ee5d6c7575050fc164991035e2"
OUTPUT = HERE / "references_rb_m4b1.json"
M4A_REFERENCES = HERE.parent / "m4a" / "references_rb_m4a.json"
SCHEMA = "xpsloss-m4b1-references-rb-v1"

COLUMNS = ("term_13a", "term_13b", "term_14", "term_15_T1", "term_15_T2", "term_15_T3",
           "term_16_e", "term_16_h", "no_hole", "hole", "total")
COLIDX = {c: i for i, c in enumerate(COLUMNS)}
COMPONENTS = ("total", "no_hole", "hole")
RB_ROWS = ("RB1", "RB2", "RB3", "RB4", "RB5", "RB6", "RB7", "RB8", "RB1a")
LEDGER = 6096

KMAX = 8                     # finest level, h = 2^-8 eV
LEVELS = (4, 5, 6, 7, 8)     # main: 5..8, part (1): 4..7
KE = 6                       # level of the carriers of parts (3) and (5)
POOL = 2 ** (KMAX - KE)      # half width (finest nodes) of the maximum taken for the carriers
OUTER_INTERVALS = 16         # part (4): Romberg intervals per primary bin
CHUNK = 1024                 # losses per kernel call (plus the two window ends)
DIGITS = 40                  # significant digits of the stored decimals
EPS53 = Fraction(1, 2 ** 53)
FR_MAX = Fraction(125, 10 ** 7)      # prereg 6.3: f_r <= 1.25e-5
DREL_CAP = Fraction(1, 10 ** 4)      # prereg 7.2: d_rel = max(8 f_r, n eps + (n-1) eps_o), at most 1e-4
H_FLOOR = Fraction(1, 10 ** 12)      # prereg 7.2: sensitivity floor n m a 1e-12 max(1, r_abs)^(n-1)
PARTS = ("finest_two", "extrapolations", "upstream", "primary_bins", "rounding")


# ----------------------------------------------------------------------------- small helpers

def fr(x):
    return Fraction(float(x))


def from_hex(text):
    return float.fromhex(text)


def dec(x, mode="nearest", digits=DIGITS):
    """Decimal string of a Fraction with `digits` significant digits (nearest, up or down in magnitude)."""
    x = Fraction(x)
    if x == 0:
        return "0"
    sign = "-" if x < 0 else ""
    x = abs(x)
    e = len(str(x.numerator)) - len(str(x.denominator))
    while Fraction(10) ** e > x:
        e -= 1
    while Fraction(10) ** (e + 1) <= x:
        e += 1
    scaled = x / Fraction(10) ** (e - digits + 1)
    q, rem = divmod(scaled.numerator, scaled.denominator)
    if mode == "up" and rem:
        q += 1
    elif mode == "nearest" and 2 * rem >= scaled.denominator:
        q += 1
    if q == 10 ** digits:
        q //= 10
        e += 1
    s = str(q)
    return f"{sign}{s[0]}.{s[1:]}e{e:+d}"


def romberg_coefficients(levels):
    """Exact coefficients of the fully extrapolated trapezoid value from `levels` values (coarse to fine)."""
    table = [[Fraction(int(i == j)) for j in range(levels)] for i in range(levels)]
    for m in range(1, levels):
        f = Fraction(4) ** m
        table = [[(f * table[i + 1][j] - table[i][j]) / (f - 1) for j in range(levels)]
                 for i in range(len(table) - 1)]
    return table[0]


C4 = romberg_coefficients(4)
C3 = romberg_coefficients(3)
AMP = sum(abs(c) for c in C4)
OUTER_W = romberg_coefficients(5)       # levels 1, 2, 4, 8, 16 intervals


def outer_node_weights():
    """Weights on the 17 nodes of [0, 1] of the 16-interval Romberg rule (sum 1)."""
    w = [Fraction(0)] * (OUTER_INTERVALS + 1)
    for level, c in enumerate(OUTER_W):
        n = 2 ** level
        step = OUTER_INTERVALS // n
        for i in range(n + 1):
            w[i * step] += c * Fraction(1, n) * (Fraction(1, 2) if i in (0, n) else 1)
    assert sum(w) == 1
    return w


OUTER_NODES = outer_node_weights()


def combine(coeffs, values):
    return sum(c * fr(v) for c, v in zip(coeffs, values))


def tw(n):
    """Trapezoid weights of n intervals in units of the step."""
    w = np.ones(n + 1)
    w[0] = w[-1] = 0.5
    return w


def conv_trap(a, b, h):
    """Trapezoid convolution of two node arrays (steps h) on the exact overlap interval of every node."""
    ma, mb = len(a) - 1, len(b) - 1
    full = np.convolve(a, b)
    p = np.arange(ma + mb + 1)
    i0 = np.maximum(0, p - mb)
    i1 = np.minimum(ma, p)
    out = h * (full - 0.5 * (a[i0] * b[p - i0] + a[i1] * b[p - i1]))
    out[0] = out[-1] = 0.0          # one-point overlap: exactly 0
    return out


def pooled_max(x, half, step):
    """Maximum of |x| over the finest nodes within `half` of every `step`-th node."""
    ax = np.abs(x)
    n = len(ax)
    idx = np.arange(0, n, step)
    out = np.zeros(len(idx))
    for d in range(-half, half + 1):
        j = np.clip(idx + d, 0, n - 1)
        out = np.maximum(out, ax[j])
    return out


def mono_delta(xs, ds):
    """prod(|x| + d) - prod(|x|) for exact rationals."""
    full = Fraction(1)
    base = Fraction(1)
    for x, d in zip(xs, ds):
        full *= abs(x) + d
        base *= abs(x)
    return full - base


# ----------------------------------------------------------------------------- table and conditions

def load_table():
    raw = TABLE_FILE.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != TABLE_SHA256:
        raise SystemExit(f"expected_states_m4b1.json has SHA-256 {digest}, not the frozen one")
    doc = json.loads(raw.decode("utf-8"))
    doc["_rows"] = {r["id"]: r for r in doc["rows"]}
    return doc


def config_id(material, momentum, charge):
    return f"{material}_{momentum}_Q{charge:g}"


def row_conditions(table, rid):
    row = table["_rows"][rid]
    inp = row["inputs"]
    pub = table["defaults"]["public"]
    first = dict(table["defaults"]["first_order_settings"])
    order = dict(table["defaults"]["order_settings"])
    for k, v in row.get("settings", {}).items():
        if k in first:
            first[k] = v
        elif k in order:
            order[k] = v
        else:
            raise SystemExit(f"{rid}: unknown setting {k}")
    mat = table["materials"][inp["material"]]
    window = [from_hex(x["hex"]) for x in inp["loss_window_eV"]]
    cond = {
        "row": rid, "entry": row["entry"], "material": inp["material"],
        "gap": from_hex(mat["gap_hex"]),
        "oscillators": [[from_hex(x) for x in o] for o in mat["oscillators_hex"]],
        "E0": float(pub["electron_energy_eV"]),
        "momentum": inp.get("epsilon_momentum", pub["epsilon_momentum"]),
        "charge": float(inp.get("hole_charge", pub["hole_charge"])),
        "L": float(inp.get("weight_length_A", pub["weight_length_A"])),
        "theta": float(pub["theta_deg"]),
        "window": window, "max_order": inp["max_order"],
        "multiplier": inp.get("multiplier", pub["multiplier"]),
        "epsrel": float(first["epsrel"]), "area_epsabs": float(first["area_epsabs"]),
        "spectrum_epsabs_fraction": float(first["spectrum_epsabs_fraction"]),
        "kernel": {"epsabs": float(first["kernel_epsabs"]), "epsrel": float(first["kernel_epsrel"]),
                   "limit": int(first["kernel_limit"])},
        "order_epsrel": float(order["order_epsrel"]),
        "order_epsabs_fraction": float(order["order_epsabs_fraction"]),
    }
    if cond["multiplier"] == "weight_length":
        cond["m"] = cond["L"]
    elif cond["multiplier"] == "external":
        cond["m"] = float(inp["multiplier_A"])
    else:
        raise SystemExit(f"{rid}: multiplier {cond['multiplier']} is not handled by R-B")
    if cond["max_order"] != 3 or inp["model"] != table["model_id"]:
        raise SystemExit(f"{rid}: unexpected model or order")
    if cond["kernel"] != {"epsabs": 1e-12, "epsrel": 1e-9, "limit": 200}:
        raise SystemExit(f"{rid}: kernel settings are not the API defaults")
    cond["config"] = config_id(cond["material"], cond["momentum"], cond["charge"])
    if row["entry"] == "api_spectrum":
        prim = table["primaries"][inp["primary"]]
        cond["primary_name"] = inp["primary"]
        cond["primary_kind"] = prim["kind"]
        cond["intensities"] = [from_hex(x) for x in prim["intensities"]["hex"]]
        key = "energies_eV" if prim["kind"] == "delta" else "edges_eV"
        cond["primary_points"] = [from_hex(x) for x in prim[key]["hex"]]
        cond["edges_name"] = inp["output_edges_eV"]
        cond["edges"] = [from_hex(x) for x in table["output_edges"][inp["output_edges_eV"]]["hex"]]
    elif row["entry"] != "api_areas":
        raise SystemExit(f"{rid}: unexpected entry {row['entry']}")
    return cond


def all_conditions(table):
    return {rid: row_conditions(table, rid) for rid in RB_ROWS}


def make_material(cond):
    return Material(cond["material"], cond["gap"], tuple(Oscillator(*o) for o in cond["oscillators"]))


# ----------------------------------------------------------------------------- point values (cache)

def grid(cond, k=KMAX):
    a, b = cond["window"]
    n = (Fraction(b) - Fraction(a)) * 2 ** k
    assert n.denominator == 1
    t = a + np.arange(int(n) + 1) / 2.0 ** k
    assert t[-1] == b and all(Fraction(float(x)) == Fraction(a) + Fraction(i, 2 ** k)
                              for i, x in ((0, t[0]), (1, t[1]), (int(n), t[-1])))
    return t


def compute_points(cond, cache):
    material = make_material(cond)
    kw = dict(weight_length_A=cond["L"], epsilon_momentum=cond["momentum"], theta_deg=cond["theta"],
              hole_charge=cond["charge"], epsabs=cond["kernel"]["epsabs"], epsrel=cond["kernel"]["epsrel"],
              limit=cond["kernel"]["limit"])
    t = grid(cond)
    a, b = cond["window"]
    vals = np.empty((len(COLUMNS), len(t)))
    errs = np.empty_like(vals)
    sentinel = None
    calls = 0
    start = time.perf_counter()
    cpu = time.process_time()
    interior = t[1:-1]
    for i in range(0, len(interior), CHUNK):
        batch = interior[i:i + CHUNK]
        losses = np.concatenate(([a, b], batch))
        r = effective_cross_section_depth_averaged(material, cond["E0"], losses, **kw)
        v = np.stack([np.asarray(getattr(r, c), float) for c in COLUMNS])
        e = np.stack([np.asarray(r.component_errors[c], float) for c in COLUMNS])
        calls += 1
        if sentinel is None:
            sentinel = (v[:, :2].copy(), e[:, :2].copy())
        elif not (np.array_equal(sentinel[0], v[:, :2]) and np.array_equal(sentinel[1], e[:, :2])):
            raise ArithmeticError("the M3a value at a window end depends on the batch")
        vals[:, 1 + i:1 + i + len(batch)] = v[:, 2:]
        errs[:, 1 + i:1 + i + len(batch)] = e[:, 2:]
    vals[:, 0], vals[:, -1] = sentinel[0][:, 0], sentinel[0][:, 1]
    errs[:, 0], errs[:, -1] = sentinel[1][:, 0], sentinel[1][:, 1]
    if not (np.isfinite(vals).all() and np.isfinite(errs).all() and (errs >= 0).all()):
        raise ArithmeticError("non-finite point value or error")
    np.savez(cache / f"points_{cond['config']}.npz", T=t, V=vals, E=errs)
    timing = {"config": cond["config"], "points": int(len(t)), "calls": calls,
              "elapsed_s": time.perf_counter() - start, "cpu_s": time.process_time() - cpu}
    (cache / f"timing_points_{cond['config']}.json").write_text(json.dumps(timing, indent=1) + "\n", encoding="ascii")
    print(json.dumps(timing))


def load_points(cond, cache):
    path = cache / f"points_{cond['config']}.npz"
    if not path.exists():
        raise SystemExit(f"missing cache file {path.name}: run --points {cond['config']}")
    z = np.load(path)
    if not np.array_equal(z["T"], grid(cond)):
        raise SystemExit(f"{path.name}: the grid differs")
    return z["V"], z["E"]


# ----------------------------------------------------------------------------- per-configuration arrays

class Config:
    """Level arrays of one (material, momentum, hole charge): K, C2, C3 without the multiplier."""

    def __init__(self, cond, cache, with_x01=False):
        self.a, self.b = Fraction(cond["window"][0]), Fraction(cond["window"][1])
        self.V, self.E = load_points(cond, cache)
        self.m_fine = self.V.shape[1] - 1
        self.gamma = (8 * self.m_fine + 64) * EPS53
        comp = [COLIDX[c] for c in COMPONENTS]
        self.K, self.C2, self.C3, self.X2, self.X3, self.area = {}, {}, {}, {}, {}, {}
        for k in LEVELS:
            step = 2 ** (KMAX - k)
            h = 2.0 ** -k
            kk = self.V[:, ::step]
            self.K[k] = kk
            r = kk[0]
            d2 = conv_trap(r, r, h)
            self.C2[k] = np.stack([conv_trap(kk[j], r, h) for j in comp])
            self.C3[k] = np.stack([conv_trap(kk[j], d2, h) for j in comp])
            self.area[k] = h * (kk.sum(axis=1) - 0.5 * (kk[:, 0] + kk[:, -1]))
            if with_x01:        # known wrong answer X01: the repeated kernel is `total`
                tot = kk[COLIDX["total"]]
                x2 = conv_trap(tot, tot, h)
                self.X2[k] = x2[None, :]
                self.X3[k] = conv_trap(tot, x2, h)[None, :]
        # carriers of parts (3) and (5) on the level-KE grid (upper envelopes, full end weights)
        step = 2 ** (KMAX - KE)
        h = 2.0 ** -KE
        kabs = np.stack([pooled_max(self.V[j], POOL, step) for j in range(len(COLUMNS))])
        kerr = np.stack([pooled_max(self.E[j], POOL, step) for j in range(len(COLUMNS))])
        rabs, rerr = kabs[0], kerr[0]
        rp = rabs + rerr
        p2 = h * np.convolve(rp, rp)
        q2 = h * np.convolve(2 * rabs + rerr, rerr)
        r2 = h * np.convolve(rabs, rabs)
        self.U2 = np.stack([h * (np.convolve(kerr[j], rp) + np.convolve(kabs[j], rerr)) for j in comp])
        self.U3 = np.stack([h * (np.convolve(kerr[j], p2) + np.convolve(kabs[j], q2)) for j in comp])
        self.A2 = np.stack([h * np.convolve(kabs[j], rabs) for j in comp])
        self.A3 = np.stack([h * np.convolve(kabs[j], r2) for j in comp])
        # window areas A_j (A^-1), exact rationals of the extrapolation, and their uncertainty parts
        hf = Fraction(1, 2 ** KMAX)
        self.A, self.dA = [], []
        for j in range(len(COLUMNS)):
            lv = [self.area[k][j] for k in LEVELS]
            main = combine(C4, lv[1:])
            e = self.E[j]
            a_abs = np.abs(self.V[j])
            up = hf * (fr(e.sum()) - (fr(e[0]) + fr(e[-1])) / 2)
            l1 = hf * (fr(a_abs.sum()) - (fr(a_abs[0]) + fr(a_abs[-1])) / 2)
            self.A.append(main)
            self.dA.append({"finest_two": abs(main - combine(C4, lv[:4])),
                            "extrapolations": abs(main - combine(C3, lv[2:])),
                            "upstream": AMP * up * (1 + self.gamma), "primary_bins": Fraction(0),
                            "rounding": AMP * self.gamma * l1})
        # norms without the multiplier (finest trapezoid of |K|; used only in the absolute floors)
        self.l1 = [hf * (fr(np.abs(self.V[j]).sum()) - (fr(abs(self.V[j][0])) + fr(abs(self.V[j][-1]))) / 2)
                   for j in range(len(COLUMNS))]
        self.Ka = [fr(self.V[j][0]) for j in range(len(COLUMNS))]
        self.dKa = [fr(self.E[j][0]) for j in range(len(COLUMNS))]


# ----------------------------------------------------------------------------- regions and weights

def regions_of(edges):
    e = [Fraction(x) for x in edges]
    out = [("outside.low", None, e[0])]
    out += [(f"bin[{i}]", e[i], e[i + 1]) for i in range(len(e) - 1)]
    out.append(("outside.high", e[-1], None))
    return out


def primary_support(cond):
    pts = [Fraction(x) for x in cond["primary_points"]]
    f = cond["intensities"]
    if cond["primary_kind"] == "delta":
        used = [p for p, s in zip(pts, f) if s > 0]
    else:
        used = [p for i, s in enumerate(f) if s > 0 for p in (pts[i], pts[i + 1])]
    if any(s < 0 for s in f) or not used:
        raise SystemExit("R-B handles primaries with non-negative intensities and a non-empty support")
    return min(used), max(used)


def reached(region, n, a, b, pmin, pmax):
    """True if losses in (n a, n b) of positive measure reach the region (exact rationals)."""
    _, u, v = region
    lo = None if v is None else pmin - v
    hi = None if u is None else pmax - u
    return (hi is None or hi > n * a) and (lo is None or lo < n * b)


def primary_sequences(cond, k, alternative=False):
    """List of (first node E' as Fraction, weights on the level-k nodes) of the primary."""
    pts = [Fraction(x) for x in cond["primary_points"]]
    f = cond["intensities"]
    if cond["primary_kind"] == "delta":
        return [(p, np.array([s])) for p, s in zip(pts, f)]
    out = []
    for i, s in enumerate(f):
        nseg = (pts[i + 1] - pts[i]) * 2 ** k
        assert nseg.denominator == 1
        nseg = int(nseg)
        if alternative:
            stride, rest = divmod(nseg, OUTER_INTERVALS)
            assert rest == 0 and stride >= 1
            seq = np.zeros(nseg + 1)
            seq[::stride] = [s * float(w) for w in OUTER_NODES]
        else:
            seq = s / nseg * tw(nseg)
        out.append((pts[i], seq))
    return out


def weight_matrix(cond, k, n, alternative=False):
    """Rows: regions; columns: nodes s = n a + p 2^-k.  J = h * (W @ C)."""
    a, b = Fraction(cond["window"][0]), Fraction(cond["window"][1])
    nm = n * (b - a) * 2 ** k
    assert nm.denominator == 1
    nm = int(nm)

    def idx(s):
        v = (s - n * a) * 2 ** k
        assert v.denominator == 1
        return int(v)

    seqs = primary_sequences(cond, k, alternative)
    regions = regions_of(cond["edges"])
    w = np.zeros((len(regions), nm + 1))
    for row, (_, u, v) in enumerate(regions):
        for x0, seq in seqs:
            if u is not None and v is not None:
                length = (v - u) * 2 ** k
                assert length.denominator == 1
                start = idx(x0 - v)
                s = np.convolve(seq, tw(int(length)))
            elif u is None:                      # (-inf, v): losses >= E' - v
                start = idx(x0 - v)
                length = nm - start + 1
                if length <= 0:
                    continue
                theta = np.ones(length)
                theta[0] = 0.5
                s = np.convolve(seq, theta)
            else:                                # (u, inf): losses <= E' - u
                pad = len(seq)
                length = idx(x0 - u) + 1 + pad
                if length <= 0:
                    continue
                theta = np.ones(length)
                theta[-1] = 0.5
                start = -pad
                s = np.convolve(seq, theta)
            lo, hi = max(start, 0), min(start + len(s), nm + 1)
            if hi > lo:
                w[row, lo:hi] += s[lo - start:hi - start]
    return w


# ----------------------------------------------------------------------------- records

def widths_bins(cond, cfg, n, column, s_abs):
    m = Fraction(cond["m"])
    eps, eps_o = Fraction(cond["epsrel"]), Fraction(cond["order_epsrel"])
    eta, eta_o = Fraction(cond["spectrum_epsabs_fraction"]), Fraction(cond["order_epsabs_fraction"])
    r_abs = m * cfg.l1[0]
    w_abs = m * cfg.l1[COLIDX[column]]
    phi = r_abs ** (n - 1) + (n - 1) * w_abs * r_abs ** (n - 2)
    psi = sum(r_abs ** (n - k) for k in range(2, n + 1))
    return (phi * eta + psi * eta_o) * s_abs, n * eps + (n - 1) * eps_o


def widths_area(cond, cfg, n, column):
    m = Fraction(cond["m"])
    eps, eps_o = Fraction(cond["epsrel"]), Fraction(cond["order_epsrel"])
    eta_a, eta_o = Fraction(cond["area_epsabs"]), Fraction(cond["order_epsabs_fraction"])
    r_abs = m * cfg.l1[0]
    w_abs = m * cfg.l1[COLIDX[column]]
    phi = r_abs ** (n - 1) + (n - 1) * w_abs * r_abs ** (n - 2)
    psi = sum(r_abs ** (n - k) for k in range(2, n + 1))
    return phi * m * eta_a + psi * eta_o, n * eps + (n - 1) * eps_o


def widths_sensitivity(cond, cfg, n):
    m, a = Fraction(cond["m"]), Fraction(cond["window"][0])
    eps, eps_o = Fraction(cond["epsrel"]), Fraction(cond["order_epsrel"])
    r_abs = m * cfg.l1[0]
    return n * m * a * H_FLOOR * max(Fraction(1), r_abs) ** (n - 1), n * eps + (n - 1) * eps_o


def record(key, kind, value, parts, d_abs, d_rel_base, exact_zero=False):
    """One reference record; applies prereg 6.3 (no reference) and 7.2 (d_rel of R-B)."""
    rec = {"key": key, "kind": kind, "oracle": "R-B", "interval": None, "no_value_reason": None}
    if exact_zero:
        rec.update(value="0", u_r="0", f_r=None, u_r_parts={p: "0" for p in PARTS}, zero_kind="exact_support",
                   width=None if d_abs is None else {"d_abs": dec(d_abs, "down"), "d_rel": dec(d_rel_base, "down")})
        return rec, Fraction(0), Fraction(0)
    parts = dict(parts)
    parts["rounding"] += abs(value) / Fraction(10) ** (DIGITS - 1)
    u = sum(parts[p] for p in PARTS)
    shown_u = Fraction(dec(u, "up"))
    rec["u_r_parts"] = {p: dec(parts[p], "up") for p in PARTS}
    rec["zero_kind"] = None
    if value == 0:
        ok = d_abs is not None and d_abs > 0 and shown_u <= d_abs / 4
        f = None
        d_rel = d_rel_base
        reason = None if ok else "zero_value_u_r_exceeds_quarter_floor"
    else:
        f = shown_u / abs(value)
        ok = f <= FR_MAX
        d_rel = None if d_rel_base is None else min(max(8 * f, d_rel_base), DREL_CAP)
        reason = None if ok else "f_r_exceeds_1.25e-5"
    rec["f_r"] = None if f is None else dec(f, "up")
    width = None if d_abs is None else {"d_abs": dec(d_abs, "down"), "d_rel": dec(d_rel, "down")}
    if ok:
        rec.update(value=dec(value), u_r=dec(u, "up"), width=width)
        if kind == "interval_contains":
            rec["interval"] = [dec(value - shown_u, "down" if value - shown_u > 0 else "up"),
                               dec(value + shown_u, "up" if value + shown_u > 0 else "down")]
    else:
        rec.update(value=None, u_r=None, width=None, no_value_reason=reason,
                   rejected={"value": dec(value), "u_r": dec(u, "up")})
    return rec, value, shown_u


def extrapolate_regions(cfg, cond, arrays, n, alternative=False, levels=LEVELS):
    """Per level J = h W C for the three components; returns {level: array (3, regions)}."""
    out = {}
    for k in levels:
        w = weight_matrix(cond, k, n, alternative)
        out[k] = 2.0 ** -k * (arrays[k] @ w.T)
    return out


def spectrum_records(cond, cfg):
    a, b = cfg.a, cfg.b
    m = Fraction(cond["m"])
    regions = regions_of(cond["edges"])
    pmin, pmax = primary_support(cond)
    s_total = sum(Fraction(x) for x in cond["intensities"])
    s_abs = sum(abs(Fraction(x)) for x in cond["intensities"])
    binned = cond["primary_kind"] == "binned"
    records, values = [], {}
    checks = {}
    for n, arrays, ucar, acar in ((2, cfg.C2, cfg.U2, cfg.A2), (3, cfg.C3, cfg.U3, cfg.A3)):
        mn = m ** n
        lev = extrapolate_regions(cfg, cond, arrays, n)
        alt = extrapolate_regions(cfg, cond, arrays, n, True, LEVELS[1:]) if binned else None
        we = np.abs(weight_matrix(cond, KE, n))
        up = 2.0 ** -KE * (ucar @ we.T)
        ab = 2.0 ** -KE * (acar @ we.T)
        sums = {c: Fraction(0) for c in COMPONENTS}
        for ci, c in enumerate(COMPONENTS):
            for ri, region in enumerate(regions):
                key = f"order[{n}].{region[0]}.{c}"
                d_abs, d_rel = widths_bins(cond, cfg, n, c, s_abs)
                if not reached(region, n, a, b, pmin, pmax):
                    if any(lev[k][ci, ri] != 0.0 for k in LEVELS):
                        raise ArithmeticError(f"{cond['row']} {key}: non-zero weight outside the support")
                    rec, val, u = record(key, "direct_reference", Fraction(0), None, d_abs, d_rel, exact_zero=True)
                else:
                    ys = [lev[k][ci, ri] for k in LEVELS]
                    main = mn * combine(C4, ys[1:])
                    parts = {"finest_two": abs(main - mn * combine(C4, ys[:4])),
                             "extrapolations": abs(main - mn * combine(C3, ys[2:])),
                             "upstream": mn * AMP * fr(up[ci, ri]) * (1 + cfg.gamma),
                             "primary_bins": (abs(main - mn * combine(C4, [alt[k][ci, ri] for k in LEVELS[1:]]))
                                              if binned else Fraction(0)),
                             "rounding": mn * AMP * cfg.gamma * fr(ab[ci, ri])}
                    rec, val, u = record(key, "direct_reference", main, parts, d_abs, d_rel)
                records.append(rec)
                values[key] = (val, u)
                sums[c] += val
        # self check (not used for any value): sum over all regions against S w_n from the window areas
        for c in COMPONENTS:
            w_n = mn * cfg.A[COLIDX[c]] * cfg.A[0] ** (n - 1)
            target = s_total * w_n
            checks[f"order[{n}].{c}"] = dec(abs(sums[c] - target) / abs(target), "up", 6) if target else None
    return records, values, checks


def area_parts(cfg, m, n, j):
    xs = [cfg.A[j]] + [cfg.A[0]] * (n - 1)
    parts = {p: m ** n * mono_delta(xs, [cfg.dA[j][p]] + [cfg.dA[0][p]] * (n - 1)) for p in PARTS}
    total = m ** n * mono_delta(xs, [sum(cfg.dA[j].values())] + [sum(cfg.dA[0].values())] * (n - 1))
    parts["rounding"] += total - sum(parts.values())      # cross terms between the parts
    return m ** n * cfg.A[j] * cfg.A[0] ** (n - 1), parts


def area_records(cond, cfg):
    m = Fraction(cond["m"])
    records, values = [], {}
    for n in (2, 3):
        for j, column in enumerate(COLUMNS):
            key = f"order[{n}].area.{column}"
            value, parts = area_parts(cfg, m, n, j)
            d_abs, d_rel = widths_area(cond, cfg, n, column)
            rec, val, u = record(key, "direct_reference", value, parts, d_abs, d_rel)
            records.append(rec)
            values[key] = (val, u)
    return records, values


def ratio_record(cond, cfg):
    m = Fraction(cond["m"])
    parts = {p: m * cfg.dA[0][p] for p in PARTS}
    rec, val, u = record("order_area_ratio", "interval_contains", m * cfg.A[0], parts, None, None)
    return rec, {"order_area_ratio": (val, u)}


def sensitivity_records(cond, cfg):
    m, a = Fraction(cond["m"]), Fraction(cond["window"][0])
    records, values = [], {}
    zero = {p: Fraction(0) for p in PARTS}

    def split(xs_list, coeffs, deltas_of):
        """Parts of sum_i coeff_i prod(xs_i); deltas_of(p) gives the uncertainties of every factor."""
        parts = {p: sum(c * mono_delta(xs, ds) for c, xs, ds in zip(coeffs, xs_list, deltas_of(p))) for p in PARTS}
        total = sum(c * mono_delta(xs, ds) for c, xs, ds in zip(coeffs, xs_list, deltas_of(None)))
        parts["rounding"] += total - sum(parts.values())
        return parts

    r, hr = m * cfg.A[0], m * a * cfg.Ka[0]

    def d_r(p):
        return m * (cfg.dA[0][p] if p else sum(cfg.dA[0].values()))

    def d_h(j, p):
        return m * a * cfg.dKa[j] if p in ("upstream", None) else Fraction(0)

    def d_w(j, p):
        return m * (cfg.dA[j][p] if p else sum(cfg.dA[j].values()))

    for n in (1, 2, 3):
        for j, column in enumerate(COLUMNS):
            key = f"order[{n}].cutoff_sensitivity.{column}"
            h1, w1 = m * a * cfg.Ka[j], m * cfg.A[j]
            if n == 1:
                value = h1
                parts = dict(zero)
                parts["upstream"] = d_h(j, "upstream")
            else:
                value = h1 * r ** (n - 1) + (n - 1) * w1 * hr * r ** (n - 2)
                xs_list = [[h1] + [r] * (n - 1), [w1, hr] + [r] * (n - 2)]
                parts = split(xs_list, [1, n - 1],
                              lambda p: [[d_h(j, p)] + [d_r(p)] * (n - 1), [d_w(j, p), d_h(0, p)] + [d_r(p)] * (n - 2)])
            d_abs, d_rel = widths_sensitivity(cond, cfg, n)
            rec, val, u = record(key, "direct_reference", value, parts, d_abs, d_rel)
            records.append(rec)
            values[key] = (val, u)
    return records, values


def expand_ledger(row):
    """Counts of the R-B key groups of the frozen table for this row."""
    out = {}
    for item in row["check_plan"]:
        if item["oracle"] == "R-B":
            out[(item["keys"], item["kind"])] = item["count"]
    return out


def identification(values, nbins):
    """Largest order-2 total bin of the reference and its two neighbours (prereg 6.3 item 5)."""
    vals = [values[f"order[2].bin[{i}].total"][0] for i in range(nbins)]
    best = max(range(nbins), key=lambda i: (vals[i], -i))
    return best, [i for i in (best - 1, best, best + 1) if 0 <= i < nbins]


def width_of(rec):
    """Comparison width d_abs + d_rel |r| of a record with a value (exact rationals of the stored decimals)."""
    r = Fraction(rec["value"])
    return Fraction(rec["width"]["d_abs"]) + Fraction(rec["width"]["d_rel"]) * abs(r)


def x01_record(cond, cfg, by_key, ident):
    """Known wrong answer X01 (the repeated kernel is `total`) against the RB1 widths."""
    m = Fraction(cond["m"])
    out = {"mutation": "X01", "row": cond["row"],
           "wrong_model": "C_n = K_total * K_total^(*(n-1)) (the repeated kernel is total instead of term_13a)",
           "orders": {}}
    regions = regions_of(cond["edges"])
    for n, arrays in ((2, cfg.X2), (3, cfg.X3)):
        lev = extrapolate_regions(cfg, cond, arrays, n)
        compared, exceeding, detail = 0, [], []
        for ri, region in enumerate(regions):
            key = f"order[{n}].{region[0]}.total"
            rec = by_key[key]
            if rec["value"] is None or rec.get("zero_kind") == "exact_support":
                continue
            ys = [lev[k][0, ri] for k in LEVELS]
            wrong = m ** n * combine(C4, ys[1:])
            u_wrong = abs(wrong - m ** n * combine(C4, ys[:4])) + abs(wrong - m ** n * combine(C3, ys[2:]))
            r, u_r = Fraction(rec["value"]), Fraction(rec["u_r"])
            width = width_of(rec)
            compared += 1
            beyond = abs(wrong - r) - u_wrong - u_r > width
            if beyond and region[0].startswith("bin["):
                exceeding.append(int(region[0][4:-1]))
            if n == 2 and region[0] in [f"bin[{i}]" for i in ident]:
                detail.append({"key": key, "r": dec(r), "y_wrong": dec(wrong), "u_wrong": dec(u_wrong, "up"),
                               "u_r": dec(u_r, "up"), "width": dec(width, "down"),
                               "ratio_(|y-r|-u_wrong-u_r)/width": dec((abs(wrong - r) - u_wrong - u_r) / width, "down", 12),
                               "exceeds_width": beyond})
        out["orders"][str(n)] = {"regions_compared": compared, "bins_exceeding_count": len(exceeding),
                                 "bins_exceeding": exceeding, "identification_bins": detail}
    return out


def momentum_record(name, ra, rb, rows):
    """Difference of the references of the two momentum prescriptions in the main order-2 bins."""
    a, b = rows[ra], rows[rb]
    ka = {r["key"]: r for r in a["records"]}
    kb = {r["key"]: r for r in b["records"]}
    out = {"rows": [ra, rb], "material": name, "bins": [], "all_order2_total_bins": None}
    compared = exceeding = 0
    nb = len(a["conditions"]["output_edges_hex"]) - 1
    for i in range(nb):
        key = f"order[2].bin[{i}].total"
        x, y = ka[key], kb[key]
        if x["value"] is None or y["value"] is None or x.get("zero_kind") or y.get("zero_kind"):
            continue
        diff = abs(Fraction(x["value"]) - Fraction(y["value"]))
        total = width_of(x) + width_of(y) + Fraction(x["u_r"]) + Fraction(y["u_r"])
        compared += 1
        exceeding += diff > total
        if i in a["identification"]["bins"]:
            out["bins"].append({"key": key, f"r_{ra}": x["value"], f"r_{rb}": y["value"], "difference": dec(diff),
                                "sum_of_widths_and_u_r": dec(total, "up"),
                                "ratio": dec(diff / total, "down", 12), "exceeds": diff > total})
    out["all_order2_total_bins"] = {"compared": compared, "exceeding": exceeding}
    return out


def m4a_area_check(table, conds, cfgs):
    """Self check: the window areas against the frozen M4a R-B references (other route: Gauss-Legendre)."""
    if not M4A_REFERENCES.exists():
        return {"available": False}
    doc = json.loads(M4A_REFERENCES.read_text(encoding="utf-8"))
    out = {"available": True, "file": "benchmarks/numerics/m4a/references_rb_m4a.json",
           "sha256": hashlib.sha256(M4A_REFERENCES.read_bytes()).hexdigest(), "matches": []}
    for rid, cond in conds.items():
        cfg = cfgs[cond["config"]]
        for mid, mrow in sorted(doc["rows"].items()):
            c = mrow.get("conditions", {})
            try:
                same = (from_hex(c["electron_energy_eV"]) == cond["E0"] and c["epsilon_momentum"] == cond["momentum"]
                        and from_hex(c["hole_charge"]) == cond["charge"] and from_hex(c["weight_length_A"]) == cond["L"]
                        and [from_hex(x) for x in c["loss_window_eV"]] == cond["window"]
                        and from_hex(c["gap_eV"]) == cond["gap"]
                        and [[from_hex(x) for x in o] for o in c["oscillators_hex"]] == cond["oscillators"])
            except (KeyError, TypeError, ValueError):
                continue
            if not same:
                continue
            worst = None
            for key, rec in sorted(mrow.get("keys", {}).items()):
                if key.startswith("window_area.") and rec.get("r") not in (None, "0"):
                    col = key.split(".", 1)[1]
                    ref = Fraction(rec["r"])
                    rel = abs(cfg.A[COLIDX[col]] - ref) / abs(ref)
                    if worst is None or rel > worst[0]:
                        worst = (rel, col)
            if worst:
                out["matches"].append({"row": rid, "m4a_row": mid, "max_relative_difference_of_window_areas":
                                       dec(worst[0], "up", 6), "column": worst[1]})
            break
    return out


def environment():
    return {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
            "xpsloss": xpsloss.__version__, "platform": platform.platform(), "machine": platform.machine(),
            "base_commit": "db67d12 (m4b1-prereg, checks pre-registered and frozen)"}


def assemble(cache):
    table = load_table()
    conds = all_conditions(table)
    started = time.perf_counter()
    cpu = time.process_time()
    cfgs = {}
    for rid, cond in conds.items():
        if cond["config"] not in cfgs:
            cfgs[cond["config"]] = Config(cond, cache, with_x01=(cond["config"] == conds["RB1"]["config"]))
    rows, counts_row, no_reference, failures, insufficient = {}, {}, [], [], []
    maxima_f = (Fraction(0), None, None)
    maxima_part = {p: (Fraction(0), None, None) for p in PARTS}
    pre = {}
    kinds = {}
    for rid in RB_ROWS:
        cond, cfg = conds[rid], cfgs[conds[rid]["config"]]
        ledger = expand_ledger(table["_rows"][rid])
        records, values, checks, ident = [], {}, None, None
        if cond["entry"] == "api_spectrum":
            recs, vals, checks = spectrum_records(cond, cfg)
            assert len(recs) == ledger[("order[2..3].<outside.low|bin[i]|outside.high>.<total|no_hole|hole>", "direct_reference")]
            records += recs
            values.update(vals)
            recs, vals = area_records(cond, cfg)
            assert len(recs) == ledger[("order[2..3].area.<11 columns>", "direct_reference")]
            records += recs
            values.update(vals)
            best, near = identification(values, len(cond["edges"]) - 1)
            ident = {"largest_order2_total_bin": best, "bins": near,
                     "keys": [f"order[2].bin[{i}].total" for i in near]}
        else:
            recs, vals = area_records(cond, cfg)
            assert len(recs) == ledger[("order[2..3].area.<11 columns>", "direct_reference")]
            records += recs
            values.update(vals)
            recs, vals = sensitivity_records(cond, cfg)
            assert len(recs) == ledger[("order[1..3].cutoff_sensitivity.<11 columns>", "direct_reference")]
            records += recs
            values.update(vals)
        rec, vals = ratio_record(cond, cfg)
        assert ledger[("order_area_ratio", "interval_contains")] == 1
        records.append(rec)
        values.update(vals)
        assert len(records) == sum(ledger.values()) and len({(r["key"], r["kind"]) for r in records}) == len(records)
        by_key = {r["key"]: r for r in records}
        missing = [r["key"] for r in records if r["value"] is None]
        for r in records:
            kinds[r["kind"]] = kinds.get(r["kind"], 0) + 1
            if r["value"] is None:
                no_reference.append({"row": rid, "key": r["key"], "reason": r["no_value_reason"],
                                     "rejected": r["rejected"]})
            elif r["f_r"] is not None:
                f = Fraction(r["f_r"])
                if f > maxima_f[0]:
                    maxima_f = (f, rid, r["key"])
                for p in PARTS:
                    rel = Fraction(r["u_r_parts"][p]) / abs(Fraction(r["value"]))
                    if rel > maxima_part[p][0]:
                        maxima_part[p] = (rel, rid, r["key"])
        reasons = []
        if 10 * len(missing) > len(records):
            reasons.append("more than one tenth of the external reference keys have no value")
        if ident and any(by_key[k]["value"] is None for k in ident["keys"]):
            reasons.append("an identification key has no value")
        if cond["entry"] == "api_areas":
            must = [f"order[{n}].area.total" for n in (2, 3)] + [f"order[{n}].cutoff_sensitivity.total" for n in (1, 2, 3)]
            if any(by_key[k]["value"] is None for k in must):
                reasons.append("a required area or sensitivity key has no value")
        if reasons:
            insufficient.append({"row": rid, "reasons": reasons})
        if ident:
            ident["all_have_values"] = all(by_key[k]["value"] is not None for k in ident["keys"])
        m = Fraction(cond["m"])
        rows[rid] = {
            "entry": cond["entry"],
            "conditions": {
                "material": cond["material"], "gap_eV": cond["gap"].hex(),
                "oscillators_hex": [[x.hex() for x in o] for o in cond["oscillators"]],
                "electron_energy_eV": cond["E0"].hex(), "epsilon_momentum": cond["momentum"],
                "hole_charge": cond["charge"].hex(), "weight_length_A": cond["L"].hex(),
                "multiplier": cond["multiplier"], "multiplier_A": cond["m"].hex(), "theta_deg": cond["theta"].hex(),
                "loss_window_eV": [x.hex() for x in cond["window"]], "max_order": cond["max_order"],
                "primary": cond.get("primary_name"), "output_edges": cond.get("edges_name"),
                "output_edges_hex": [x.hex() for x in cond.get("edges", [])],
                "kernel_settings": cond["kernel"],
                "comparison_settings": {k: cond[k] for k in ("epsrel", "area_epsabs", "spectrum_epsabs_fraction",
                                                             "order_epsrel", "order_epsabs_fraction")},
                "point_value_config": cond["config"],
            },
            "norms": {"note": "finest-grid trapezoid of |K|, multiplier included (dimensionless); used in the floors d_abs only",
                      "r_abs": dec(m * cfg.l1[0]), "w_abs": {c: dec(m * cfg.l1[j]) for j, c in enumerate(COLUMNS)}},
            "primary_abs_integral": dec(sum(abs(Fraction(x)) for x in cond["intensities"])) if ident else None,
            "identification": ident,
            "self_check_sum_over_regions_vs_S_w_n_relative": checks,
            "counts": {"keys": len(records), "with_value": len(records) - len(missing), "no_reference": len(missing),
                       "exact_zero": sum(r.get("zero_kind") == "exact_support" for r in records)},
            "records": records,
        }
        counts_row[rid] = len(records)
        if rid == "RB1":
            pre["X01_wrong_answer_vs_RB1_width"] = x01_record(cond, cfg, by_key, ident["bins"])
        print(f"{rid}: {len(records)} keys, {len(missing)} without value", file=sys.stderr)
    pre["momentum_prescriptions"] = [momentum_record("Al", "RB1", "RB2", rows), momentum_record("Si", "RB4", "RB5", rows)]
    total = sum(counts_row.values())
    if total != LEDGER or total != table["counts"]["planned_numeric_keys_by_group"]["RM-B"]:
        raise SystemExit(f"{total} keys, the ledger has {LEDGER}")
    doc = {
        "schema": SCHEMA,
        "prereg": "docs/M4B1_PREREG.md (1, 4, 5, 6.2, 6.3, 7, appendix R-B)",
        "table": {"file": "benchmarks/numerics/m4b1/expected_states_m4b1.json", "sha256": TABLE_SHA256},
        "generator": {"file": "benchmarks/numerics/m4b1/generate_references_rb_m4b1.py",
                      "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        "kernel": "xpsloss.effective_cross_section_depth_averaged (public M3a API), the only package function called; "
                  "the shared point value itself is outside the scope of the comparison",
        "environment": environment(),
        "method": {
            "grid": "T = a + i 2^-k eV, window ends, bin edges and primary edges on the grid",
            "levels": list(LEVELS), "main_levels": list(LEVELS[1:]), "extrapolations": 3,
            "richardson_coefficients": [str(c) for c in C4], "amplification": str(AMP),
            "convolution": "trapezoid rule on the exact overlap interval (full discrete convolution minus half of the "
                           "two end terms; the one-point overlaps at n a and n b are exactly 0)",
            "regions": "trapezoid rule with a positive weight vector per region (difference of the cumulative "
                       "distribution as one local sum); binned primary: trapezoid rule over E' on the same grid",
            "areas": "window areas extrapolated separately over the whole window; w_n = m^n A_j A_13a^(n-1)",
            "u_r_parts": list(PARTS),
            "carrier_level": KE, "carrier_pool_half_width_finest_nodes": POOL,
            "rounding_factor": "(8 M + 64) 2^-53, M = finest intervals of the window",
            "outer_alternative": f"{OUTER_INTERVALS}-interval Romberg rule per primary bin",
            "points_per_call": CHUNK, "digits_stored": DIGITS,
            "no_reference_rule": "r != 0: f_r > 1.25e-5; unreachable regions are exact zeros with u_r = 0",
            "d_rel": "min(max(8 f_r, n epsrel + (n-1) order_epsrel), 1e-4)",
        },
        "counts": {"by_row": counts_row, "by_kind": kinds, "total": total, "ledger": LEDGER,
                   "with_value": total - len(no_reference), "no_reference": len(no_reference),
                   "reference_path_failures": len(failures),
                   "exact_zero": sum(r["counts"]["exact_zero"] for r in rows.values())},
        "maxima": {"f_r": {"value": dec(maxima_f[0], "up"), "row": maxima_f[1], "key": maxima_f[2]},
                   "u_r_part_relative": {p: {"value": dec(v[0], "up"), "row": v[1], "key": v[2]}
                                         for p, v in maxima_part.items()}},
        "no_reference": no_reference,
        "reference_path_failures": failures,
        "insufficient_rows": insufficient,
        "pre_comparison_records": pre,
        "self_checks": {"window_areas_vs_m4a_references": m4a_area_check(table, conds, cfgs),
                        "batch_independence": "every kernel call contains the two window ends; their values and "
                                              "errors are bit-identical in all calls (checked when the cache is written)"},
        "rows": rows,
    }
    text = json.dumps(doc, indent=1, sort_keys=True, ensure_ascii=True) + "\n"
    timing = {"assemble_elapsed_s": time.perf_counter() - started, "assemble_cpu_s": time.process_time() - cpu}
    (cache / "timing_assemble.json").write_text(json.dumps(timing, indent=1) + "\n", encoding="ascii")
    print(json.dumps(timing), file=sys.stderr)
    return text


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cache", required=True, type=Path)
    parser.add_argument("--points", help="configuration id, comma separated, or 'all'")
    parser.add_argument("--assemble", action="store_true")
    parser.add_argument("--check", action="store_true", help="with --assemble: compare with the stored file")
    parser.add_argument("--list", action="store_true")
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    table = load_table()
    conds = all_conditions(table)
    configs = {}
    for cond in conds.values():
        configs.setdefault(cond["config"], cond)
    if args.list:
        print("\n".join(configs))
    if args.points:
        names = list(configs) if args.points == "all" else args.points.split(",")
        for name in names:
            compute_points(configs[name], args.cache)
    if args.assemble:
        text = assemble(args.cache)
        if args.check:
            same = OUTPUT.read_bytes() == text.encode("ascii")
            print("identical" if same else "DIFFERENT")
            raise SystemExit(0 if same else 1)
        with open(OUTPUT, "w", encoding="ascii", newline="\n") as f:
            f.write(text)
        print(f"wrote {OUTPUT.name} sha256 {hashlib.sha256(text.encode('ascii')).hexdigest()}")


if __name__ == "__main__":
    main()
