"""M4b-1, R-A: analytic references of the multiple-loss rows (prereg 6.1, 6.3, 7.2, 9).

    python benchmarks/numerics/m4b1/generate_references_ra_m4b1.py            # writes references_ra_m4b1.json
    python benchmarks/numerics/m4b1/generate_references_ra_m4b1.py --check    # regenerates and compares the bytes

The generator is deterministic and does not import the package under test.  The
binary64 inputs of the frozen table are read as exact rationals.  Box kernels
are evaluated with exact rational truncated powers (the cumulative distribution
of a sum of uniform variables of different widths: one term per subset of the
shifted ends, sign by the parity of the subset, denominator n!); polynomial
kernels with exact rational piecewise polynomials; the tilted box, the
Lorentzian and the truncation bounds with mpmath at two working precisions
(40 digits or more), the difference entering u_r.

Every record is addressed by (row, key, kind, oracle) of the frozen ledger and
holds the value (40 significant digits), u_r (which contains the decimal
rounding), the comparison width of prereg 7.2 and the status of prereg 6.3,
applied here before any implementation value exists.
"""
from __future__ import annotations

import contextlib
import hashlib
import io
import json
import math
import pathlib
import platform
import re
import sys
import time
from fractions import Fraction as Fr

import mpmath
from mpmath import mp
from mpmath.libmp import to_rational

HERE = pathlib.Path(__file__).resolve().parent
TABLE_JSON = HERE / "expected_states_m4b1.json"
TABLE_PY = HERE / "expected_states_m4b1.py"
OUT = HERE / "references_ra_m4b1.json"
TABLE_SHA256 = "d8d0f892a49c3ed9933b08d6a7c92e5b3aaf05ee5d6c7575050fc164991035e2"
TABLE_PY_SHA256 = "7c06f0bab21b5f195696fcfe6fbd4bdb0232a7bcf7d9b5f994bd3e2bfbbc248b"
SCHEMA = "xpsloss-m4b1-references-ra-v1"
DIGITS = 40
DPS_LOW, DPS_HIGH = 45, 65
EPS = Fr(1, 2**52)
COLUMNS = ("term_13a", "term_13b", "term_14", "term_15_T1", "term_15_T2", "term_15_T3", "term_16_e", "term_16_h",
           "no_hole", "hole", "total")
NO_HOLE = ("term_13a", "term_13b", "term_15_T1", "term_15_T3", "term_16_e")
HOLE = ("term_14", "term_15_T2", "term_16_h")
COMPONENTS = ("total", "no_hole", "hole")
EXTERNAL = ("direct_reference", "bound_reference", "error_part", "interval_contains")


class ReferenceError_(RuntimeError):
    pass


def need(cond, message):
    if not cond:
        raise ReferenceError_(message)


# ---------------------------------------------------------------- the frozen table and its ledger
def load_table():
    raw = TABLE_JSON.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == TABLE_SHA256, "the expected-state table is not the frozen one")
    return json.loads(raw.decode("utf-8"))


def load_ledger():
    """Concrete (key, kind, oracle) records per row, expanded by the frozen table script itself."""
    src = TABLE_PY.read_bytes()
    need(hashlib.sha256(src).hexdigest() == TABLE_PY_SHA256, "the table script is not the frozen one")
    scope = {"__file__": str(TABLE_PY), "__name__": "m4b1_expected_states_table"}
    saved = sys.argv
    sys.argv = [str(TABLE_PY), "--check"]
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            exec(compile(src, str(TABLE_PY), "exec"), scope)
    except SystemExit as stop:
        need(stop.code == 0, "the table script does not reproduce the frozen table")
    finally:
        sys.argv = saved
    return {rid: [tuple(rec) for rec in records] for rid, records in scope["LEDGER"].items()}


def hexfr(obj):
    """Exact rational of a table number {"hex", "value"} (or of a plain number)."""
    if isinstance(obj, dict):
        x = float.fromhex(obj["hex"])
        need(x == float(obj["value"]), f"hex and decimal differ: {obj}")
        return Fr(x)
    return Fr(float(obj))


# ---------------------------------------------------------------- decimal output
def dec(x):
    """(string with DIGITS significant digits, bound of the rounding error) of an exact rational."""
    x = Fr(x)
    if x == 0:
        return "0", Fr(0)
    ax = abs(x)
    e = len(str(ax.numerator)) - len(str(ax.denominator))
    while Fr(10)**e > ax:
        e -= 1
    while Fr(10)**(e + 1) <= ax:
        e += 1
    scaled = ax / Fr(10)**(e - DIGITS + 1)
    q = (2*scaled.numerator + scaled.denominator) // (2*scaled.denominator)
    exact = Fr(q) == scaled
    if q == 10**DIGITS:
        q //= 10
        e += 1
    digits = str(q)
    text = f"{'-' if x < 0 else ''}{digits[0]}.{digits[1:]}e{e:+d}"
    return text, (Fr(0) if exact else Fr(10)**(e - DIGITS + 1)/2)


def dec_up(x, digits=3):
    """Decimal string of a non-negative rational rounded up to a few digits."""
    x = Fr(x)
    need(x >= 0, "negative uncertainty")
    if x == 0:
        return "0"
    e = len(str(x.numerator)) - len(str(x.denominator))
    while Fr(10)**e > x:
        e -= 1
    while Fr(10)**(e + 1) <= x:
        e += 1
    scaled = x / Fr(10)**(e - digits + 1)
    q = -((-scaled.numerator) // scaled.denominator)
    return f"{q}e{e - digits + 1:+d}"


def dec_down(x, digits=DIGITS):
    """Decimal string of a non-negative rational rounded down (used for widths: never wider than registered)."""
    x = Fr(x)
    need(x >= 0, "negative width")
    if x == 0:
        return "0"
    e = len(str(x.numerator)) - len(str(x.denominator))
    while Fr(10)**e > x:
        e -= 1
    while Fr(10)**(e + 1) <= x:
        e += 1
    scaled = x / Fr(10)**(e - digits + 1)
    q = scaled.numerator // scaled.denominator
    return f"{q}e{e - digits + 1:+d}"


def mpfr(x):
    """Exact rational of an mpf."""
    p, q = to_rational(mp.mpf(x)._mpf_)
    return Fr(p, q)


def tomp(x):
    x = Fr(x)
    return mp.mpf(x.numerator) / mp.mpf(x.denominator)


def two_precisions(function):
    """Evaluate at two working precisions; (value at the higher one, uncertainty)."""
    out = []
    for dps in (DPS_LOW, DPS_HIGH):
        with mp.workdps(dps):
            out.append(mpfr(function()))
    low, high = out
    return high, 4*abs(high - low) + abs(high)*Fr(1, 10**(DPS_LOW - 2))


# ---------------------------------------------------------------- truncated powers (box kernels)
_FACT = [1]


def fact(n):
    while len(_FACT) <= n:
        _FACT.append(_FACT[-1]*len(_FACT))
    return _FACT[n]


class TP:
    """sum of coef*(s - c)_+^p/p!  (exact).  A box c*1[a,b] is c*(s-a)_+^0 - c*(s-b)_+^0."""
    __slots__ = ("t",)

    def __init__(self, terms=None):
        self.t = {k: v for k, v in (terms or {}).items() if v != 0}

    @staticmethod
    def box(value, lo, hi):
        out = TP()
        out._add((lo, 0), value)
        out._add((hi, 0), -value)
        return out

    def _add(self, key, value):
        v = self.t.get(key, 0) + value
        if v == 0:
            self.t.pop(key, None)
        else:
            self.t[key] = v

    def __add__(self, other):
        out = TP(self.t)
        for k, v in other.t.items():
            out._add(k, v)
        return out

    def scale(self, factor):
        return TP({k: v*factor for k, v in self.t.items()})

    def __sub__(self, other):
        return self + other.scale(-1)

    def conv(self, other):
        out = TP()
        for (c1, p1), v1 in self.t.items():
            for (c2, p2), v2 in other.t.items():
                out._add((c1 + c2, p1 + p2 + 1), v1*v2)
        return out

    def antider(self):
        return TP({(c, p + 1): v for (c, p), v in self.t.items()})

    def __call__(self, s):
        """Right-continuous value at s."""
        total = Fr(0)
        for (c, p), v in self.t.items():
            if s >= c:
                total += v*(s - c)**p/fact(p) if p else v
        return total

    def breaks(self):
        return sorted({c for c, _ in self.t})

    def total(self):
        """Value to the right of the last break (the total of a cumulative of a compact density)."""
        b = self.breaks()
        return self(b[-1] + 1) if b else Fr(0)

    def l1_constant(self):
        need(all(p == 0 for _, p in self.t), "l1_constant needs a piecewise constant function")
        b = self.breaks()
        return sum((abs(self(lo))*(hi - lo) for lo, hi in zip(b, b[1:])), Fr(0))

    def first_moment_constant(self):
        b = self.breaks()
        return sum((self(lo)*(hi*hi - lo*lo)/2 for lo, hi in zip(b, b[1:])), Fr(0))


def tp_power(f, n):
    out = f
    for _ in range(n - 1):
        out = out.conv(f)
    return out


# ---------------------------------------------------------------- exact piecewise polynomials
def ptrim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p


def padd(p, q):
    n = max(len(p), len(q))
    return ptrim([(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)])


def pscale(p, k):
    return ptrim([c*k for c in p])


def pmul(p, q):
    if not p or not q:
        return []
    out = [Fr(0)]*(len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a*b
    return ptrim(out)


def pshift(p, c):
    """Coefficients of p(s - c)."""
    out = [Fr(0)]*len(p)
    for k, a in enumerate(p):
        for j in range(k + 1):
            out[j] += a*math.comb(k, j)*(-c)**(k - j)
    return ptrim(out)


def pint(p):
    return ptrim([Fr(0)] + [c/(k + 1) for k, c in enumerate(p)])


def pval(p, x):
    total = Fr(0)
    for c in reversed(p):
        total = total*x + c
    return total


class PP:
    """Piecewise polynomial: 0 left of the first break, ``tail`` right of the last one."""

    def __init__(self, pieces, tail=()):
        self.p = sorted(((lo, hi, ptrim(c)) for lo, hi, c in pieces if lo < hi), key=lambda x: x[0])
        for (_, h0, _), (l1, _, _) in zip(self.p, self.p[1:]):
            need(h0 <= l1, "overlapping pieces")
        self.tail = ptrim(tail)

    @staticmethod
    def from_tp(tp, grid=()):
        b = sorted(set(tp.breaks()) | set(grid))
        terms = sorted(tp.t.items())

        def poly_right_of(x):
            out = []
            for (c, p), v in terms:
                if c <= x:
                    out = padd(out, pshift([Fr(0)]*p + [Fr(v, fact(p))], c))
            return out
        return PP([(lo, hi, poly_right_of(lo)) for lo, hi in zip(b, b[1:])], poly_right_of(b[-1]) if b else ())

    def end(self):
        return self.p[-1][1] if self.p else Fr(0)

    def __call__(self, s):
        """Right-continuous value."""
        for lo, hi, c in self.p:
            if lo <= s < hi:
                return pval(c, s)
        if self.p and s >= self.end():
            return pval(self.tail, s)
        return Fr(0)

    def antider(self):
        """Cumulative from minus infinity (continuous; the tail continues the last value)."""
        out, value, last = [], Fr(0), None
        for lo, hi, c in self.p:
            if last is not None and lo > last:
                out.append((last, lo, [value]))
            q = pint(c)
            q = padd(q, [value - pval(q, lo)])
            out.append((lo, hi, q))
            value, last = pval(q, hi), hi
        t = pint(self.tail)
        return PP(out, padd(t, [value - pval(t, last)]) if last is not None else ())

    def scale(self, k):
        return PP([(lo, hi, pscale(c, k)) for lo, hi, c in self.p], pscale(self.tail, k))

    def breaks(self):
        return sorted({x for lo, hi, _ in self.p for x in (lo, hi)})

    def l1_linear(self):
        """Exact integral of |f| for pieces of degree <= 1 (split at the zero)."""
        need(not self.tail, "l1 of a non-compact function")
        total = Fr(0)
        for lo, hi, c in self.p:
            need(len(c) <= 2, "l1_linear needs degree <= 1")
            cuts = [lo, hi]
            if len(c) == 2:
                z = -c[0]/c[1]
                if lo < z < hi:
                    cuts = [lo, z, hi]
            q = pint(c)
            total += sum((abs(pval(q, y) - pval(q, x)) for x, y in zip(cuts, cuts[1:])), Fr(0))
        return total


def merge_pieces(contributions):
    """Sum overlapping (lo, hi, poly) contributions on the union of their breaks."""
    b = sorted({x for lo, hi, _ in contributions for x in (lo, hi)})
    out = []
    for lo, hi in zip(b, b[1:]):
        total = []
        for l0, h0, c in contributions:
            if l0 <= lo and hi <= h0:
                total = padd(total, c)
        out.append((lo, hi, total))
    return out


def _conv_piece(p, a1, b1, q, a2, b2):
    """int p(t) q(s - t) dt over t in [max(a1, s - b2), min(b1, s - a2)] as pieces in s."""
    H = {}
    for i, pi in enumerate(p):
        if pi == 0:
            continue
        for k, qk in enumerate(q):
            if qk == 0:
                continue
            for j in range(k + 1):
                key = (k - j, i + j + 1)
                H[key] = H.get(key, 0) + pi*qk*math.comb(k, j)*(-1)**j/Fr(i + j + 1)

    def at_constant(c):
        out = []
        for (es, et), v in H.items():
            out = padd(out, [Fr(0)]*es + [v*c**et])
        return out

    def at_shift(c):
        out = []
        for (es, et), v in H.items():
            out = padd(out, pmul([Fr(0)]*es + [v], pshift([Fr(0)]*et + [Fr(1)], c)))
        return out
    cuts = sorted({a1 + a2, a1 + b2, b1 + a2, b1 + b2})
    out = []
    for lo, hi in zip(cuts, cuts[1:]):
        mid = (lo + hi)/2
        upper = at_constant(b1) if b1 <= mid - a2 else at_shift(a2)
        lower = at_constant(a1) if a1 >= mid - b2 else at_shift(b2)
        out.append((lo, hi, padd(upper, pscale(lower, -1))))
    return out


def pp_conv(f, g):
    need(not f.tail and not g.tail, "convolution of non-compact functions")
    contributions = []
    for a1, b1, p in f.p:
        for a2, b2, q in g.p:
            if p and q:
                contributions += _conv_piece(p, a1, b1, q, a2, b2)
    return PP(merge_pieces(contributions))


# ---------------------------------------------------------------- inputs of a row
class Config:
    """Exact inputs of a row: kernel columns, window, multiplier, primary, edges and the settings of the widths."""

    def __init__(self, table, row, **override):
        self.table, self.row = table, row
        inp = dict(row["inputs"])
        inp.update(override)
        self.inp = inp
        d = table["defaults"]
        settings = {**d["first_order_settings"], **d["order_settings"], **row["settings"]}
        self.eps = Fr(float(settings["epsrel"]))
        self.eta = Fr(float(settings["spectrum_epsabs_fraction"]))
        self.eta_a = Fr(float(settings["area_epsabs"]))
        self.eps_o = Fr(float(settings["order_epsrel"]))
        self.eta_o = Fr(float(settings["order_epsabs_fraction"]))
        self.charge = Fr(float(inp.get("hole_charge", d["analytic"]["hole_charge"])))
        self.length = Fr(float(inp.get("weight_length_A", d["analytic"]["weight_length_A"])))
        self.kname = inp["kernel"]
        spec = table["kernels"][self.kname]
        self.spec = spec
        window = inp.get("window_override_eV", spec["window_eV"])
        self.a, self.b = hexfr(window[0]), hexfr(window[1])
        kinds = {p["type"] for pieces in spec["columns"].values() for p in pieces}
        need(len(kinds) == 1, "mixed kernel")
        self.kind = kinds.pop()
        self.N = inp["max_order"] if isinstance(inp["max_order"], int) else row["expect"]["computed"]
        if "primary" in inp:
            p = table["primaries"][inp["primary"]]
            self.primary_kind = p["kind"]
            key = "energies_eV" if p["kind"] == "delta" else "edges_eV"
            self.px = [Fr(float.fromhex(h)) for h in p[key]["hex"]]
            self.ps = [Fr(float.fromhex(h)) for h in p["intensities"]["hex"]]
            self.s_abs = sum(abs(s) for s in self.ps)
            self.edges = [Fr(float.fromhex(h)) for h in table["output_edges"][inp["output_edges_eV"]]["hex"]]
            self.regions = ([("outside.low", None, self.edges[0])]
                            + [(f"bin[{i}]", lo, hi) for i, (lo, hi) in enumerate(zip(self.edges, self.edges[1:]))]
                            + [("outside.high", self.edges[-1], None)])
        if self.kind == "box":
            self.K = self.box_columns()
            self.areas = {c: self.K[c].antider().total() for c in COLUMNS}
            self.norms = {c: self.K[c].l1_constant() for c in COLUMNS}
            self.m = self.multiplier(self.areas["total"])

    def box_columns(self):
        K = {c: TP() for c in COLUMNS}
        for name, pieces in self.spec["columns"].items():
            for p in pieces:
                lo, hi = (hexfr(x) for x in p["interval_eV"])
                need(self.a <= lo < hi <= self.b, "box outside the window")
                K[name] = K[name] + TP.box(hexfr(p["value"]), lo, hi)
        self.derive(K, lambda x, y: x + y, lambda x, k: x.scale(k), TP())
        return K

    def derive(self, K, add, scale, zero):
        nh, h = zero, zero
        for c in NO_HOLE:
            nh = add(nh, K[c])
        for c in HOLE:
            h = add(h, K[c])
        K["no_hole"], K["hole"], K["total"] = nh, h, add(nh, scale(h, self.charge))

    def multiplier(self, area_total):
        kind = self.inp.get("multiplier", "weight_length")
        if kind == "weight_length":
            return self.length
        if kind == "external":
            return Fr(float(self.inp["multiplier_A"]))
        need(kind == "inverse_window_area", kind)
        return 1/area_total

    # ---- widths of prereg 7.2
    def phi(self, n, wabs, rabs):
        if n == 1:
            return Fr(1)
        if n == 2:
            return rabs + wabs
        return rabs**(n - 1) + (n - 1)*wabs*rabs**(n - 2)

    def psi(self, n, rabs):
        return sum((rabs**(n - k) for k in range(2, n + 1)), Fr(0))

    def width(self, quantity, n, column, exact_polynomial=False):
        """(d_abs, d_rel) of prereg 7.2; the norms carry the multiplier."""
        wabs, rabs = self.m*self.norms[column], self.m*self.norms["term_13a"]
        d_rel = n*self.eps + (n - 1)*self.eps_o
        if quantity == "bin":
            if exact_polynomial:
                return Fr(1, 10**13)*self.s_abs, Fr(1, 10**12)
            return (self.phi(n, wabs, rabs)*self.eta + self.psi(n, rabs)*self.eta_o)*self.s_abs, d_rel
        if quantity == "area":
            if exact_polynomial:
                return Fr(1, 10**15), Fr(1, 10**12)
            return self.phi(n, wabs, rabs)*self.m*self.eta_a + self.psi(n, rabs)*self.eta_o, d_rel
        if quantity == "sensitivity":
            return n*self.m*self.a*Fr(1, 10**12)*max(Fr(1), rabs)**(n - 1), d_rel
        raise ReferenceError_(quantity)

    def spectrum_width(self):
        N = self.N
        d_abs = self.eta*self.s_abs + sum((self.width("bin", n, "total")[0] for n in range(1, N + 1)), Fr(0))
        return d_abs, N*self.eps + (N - 1)*self.eps_o

    # ---- masses of a density given by its cumulative
    def loss_interval(self, p, lo, hi):
        return (None if hi is None else p - hi), (None if lo is None else p - lo)

    def masses(self, cumulative, total, second=None):
        """Region masses (exact) of a loss density with cumulative ``cumulative`` (and its antiderivative)."""
        out = {}
        for name, lo, hi in self.regions:
            value = Fr(0)
            if self.primary_kind == "delta":
                for p, s in zip(self.px, self.ps):
                    if s == 0:
                        continue
                    l1, l2 = self.loss_interval(p, lo, hi)
                    upper = total if l2 is None else cumulative(l2)
                    lower = Fr(0) if l1 is None else cumulative(l1)
                    value += s*(upper - lower)
            else:
                for (e0, e1), s in zip(zip(self.px, self.px[1:]), self.ps):
                    if s == 0:
                        continue
                    rho = s/(e1 - e0)
                    upper = total*(e1 - e0) if lo is None else second(e1 - lo) - second(e0 - lo)
                    lower = Fr(0) if hi is None else second(e1 - hi) - second(e0 - hi)
                    value += rho*(upper - lower)
            out[name] = value
        return out

    def tp_masses(self, density):
        g = density.antider()
        return self.masses(g, g.total(), g.antider() if self.primary_kind != "delta" else None)

    def pp_masses(self, density):
        g = density.antider()
        need(len(g.tail) <= 1, "cumulative tail")
        need(self.primary_kind == "delta", "binned primaries use the truncated powers")
        return self.masses(g, pval(g.tail, 0))


def b1_orders(cfg, N=None, K=None, R=None, m=None):
    """C_{n,j} = m^n K_j * R^{*(n-1)} for the 11 columns (box kernels, truncated powers)."""
    K = cfg.K if K is None else K
    R = K["term_13a"] if R is None else R
    m = cfg.m if m is None else m
    orders, rep = {}, None
    for n in range(1, (cfg.N if N is None else N) + 1):
        orders[n] = {c: (K[c] if rep is None else K[c].conv(rep)).scale(m**n) for c in COLUMNS}
        rep = R if rep is None else rep.conv(R)
    return orders


# ---------------------------------------------------------------- records
class Book:
    def __init__(self, ledger):
        self.ledger = ledger
        self.rows, self.no_reference, self.wrong, self.self_checks = {}, [], [], {}
        self.max_ratio = {}

    def planned(self, rid):
        return {(k, kind): oracle for k, kind, oracle in self.ledger[rid]
                if kind in EXTERNAL and (oracle.startswith("R-A") or oracle.startswith("R-INV"))}

    def put(self, rid, key, kind, value, u=Fr(0), *, d_abs=None, d_rel=None, klass, interval=None, **extra):
        planned = self.planned(rid)
        need((key, kind) in planned, f"{rid}: {key} [{kind}] is not in the ledger")
        row = self.rows.setdefault(rid, {"keys": {}})
        need(key not in row["keys"], f"{rid}: duplicate reference key {key}")
        value = Fr(value)
        text, rounding = dec(value)
        u_total = Fr(u) + rounding
        rec = {"kind": kind, "oracle": planned[(key, kind)], "class": klass, "value": text,
               "u_r": dec_up(u_total), "zero_reference": value == 0, "interval": interval,
               "no_value_reason": None, "status": "reference"}
        if value == 0:
            need(u_total == 0, f"{rid} {key}: a zero reference must be exact")
        if d_abs is not None:
            rec["width"] = {"d_abs": dec_down(d_abs), "d_rel": dec_down(d_rel)}
            width = Fr(rec["width"]["d_abs"]) + (Fr(rec["width"]["d_rel"])*abs(value) if value else 0)
            u_stored = Fr(rec["u_r"])
            if width > 0:
                ratio = u_stored/width
                if ratio > self.max_ratio.get(rid, (Fr(-1), None))[0]:
                    self.max_ratio[rid] = (ratio, key)
            if 4*u_stored > width:
                rec["status"], rec["no_value_reason"] = "no_reference", "u_r > width/4 (prereg 6.3)"
                self.no_reference.append({"row": rid, "key": key, "kind": kind,
                                          "condition": "prereg 6.3 item 1" if value == 0 else "prereg 6.3 item 2"})
        else:
            rec["width"] = None
        rec.update(extra)
        row["keys"][key] = rec
        return rec

    def note(self, rid, name, value):
        self.rows.setdefault(rid, {"keys": {}}).setdefault("notes", {})[name] = value

    def wrong_answer(self, mutation, rid, key, wrong, reference, d_abs, d_rel, text, u=Fr(0), rule="outside the width"):
        wrong, reference = Fr(wrong), Fr(reference)
        width = Fr(d_abs) + Fr(d_rel)*abs(reference)
        outside = abs(wrong - reference) - Fr(u) > width
        self.wrong.append({"mutation": mutation, "row": rid, "key": key, "wrong_value": dec(wrong)[0],
                           "reference": dec(reference)[0], "difference": dec(wrong - reference)[0],
                           "width": dec_down(width), "outside_width": bool(outside), "rule": rule, "note": text})
        return outside


REGION_KEY = re.compile(r"^order\[(\d+)\]\.(outside\.low|bin\[\d+\]|outside\.high)\.(total|no_hole|hole)$")
AREA_KEY = re.compile(r"^order\[(\d+)\]\.area\.(\w+)$")
SENS_KEY = re.compile(r"^order\[(\d+)\]\.cutoff_sensitivity\.(\w+)$")
SPECTRUM_KEY = re.compile(r"^spectrum\.(outside\.low|bin\[\d+\]|outside\.high)$")


def put_standard(book, rid, cfg, masses, areas, sens=None, exact_polynomial=False, u=None):
    """Fill the region, area and sensitivity keys of a row from exact (or mp) values.

    masses[n][component][region], areas[n][column], sens[n][column]; u has the same shape (optional)."""
    def unc(kind, n, name, region=None):
        if u is None:
            return Fr(0)
        box = u[kind][n][name]
        return box[region] if region is not None else box
    for (key, kind), _ in sorted(book.planned(rid).items()):
        if kind != "direct_reference":
            continue
        m = REGION_KEY.match(key)
        if m and masses is not None:
            n, region, comp = int(m.group(1)), m.group(2), m.group(3)
            d_abs, d_rel = cfg.width("bin", n, comp, exact_polynomial)
            book.put(rid, key, kind, masses[n][comp][region], unc("masses", n, comp, region), d_abs=d_abs,
                     d_rel=d_rel, klass="region", order=n, region=region, component=comp)
            continue
        m = AREA_KEY.match(key)
        if m and areas is not None:
            n, col = int(m.group(1)), m.group(2)
            d_abs, d_rel = cfg.width("area", n, col, exact_polynomial)
            book.put(rid, key, kind, areas[n][col], unc("areas", n, col), d_abs=d_abs, d_rel=d_rel, klass="area",
                     order=n, column=col)
            continue
        m = SENS_KEY.match(key)
        if m and sens is not None:
            n, col = int(m.group(1)), m.group(2)
            d_abs, d_rel = cfg.width("sensitivity", n, col)
            book.put(rid, key, kind, sens[n][col], unc("sens", n, col), d_abs=d_abs, d_rel=d_rel,
                     klass="sensitivity", order=n, column=col)


def box_values(cfg, orders=None):
    orders = b1_orders(cfg) if orders is None else orders
    masses = None
    if hasattr(cfg, "regions"):
        masses = {n: {c: cfg.tp_masses(cols[c]) for c in COMPONENTS} for n, cols in orders.items()}
    areas = {n: {c: cols[c].antider().total() for c in COLUMNS} for n, cols in orders.items()}
    return orders, masses, areas


def box_sensitivity(cfg, N):
    """h_{n,j} = h_{1,j} r^(n-1) + (n-1) w_{1,j} h_r r^(n-2), h_{1,j} = m a K_j(a+0)."""
    m, a = cfg.m, cfg.a
    h1 = {c: m*a*cfg.K[c](a) for c in COLUMNS}
    w1 = {c: m*cfg.areas[c] for c in COLUMNS}
    r, hr = w1["term_13a"], h1["term_13a"]
    out = {}
    for n in range(1, N + 1):
        out[n] = {c: h1[c] if n == 1 else h1[c]*r**(n - 1) + (n - 1)*w1[c]*hr*r**(n - 2) for c in COLUMNS}
    return out


# ---------------------------------------------------------------- box rows
def check_box_against_pp(book, cfg, orders, rid):
    """Self-check: the truncated-power masses equal the piecewise-polynomial convolution (orders <= 3)."""
    K = {c: PP.from_tp(cfg.K[c]) for c in ("total", "term_13a")}
    c2 = pp_conv(K["total"], K["term_13a"])
    c3a = pp_conv(c2, K["term_13a"])
    c3b = pp_conv(K["total"], pp_conv(K["term_13a"], K["term_13a"]))
    same = c3a.p == c3b.p
    ok = same
    for n, density in ((2, c2), (3, c3a)):
        if n > cfg.N or cfg.primary_kind != "delta":
            continue
        reference = cfg.tp_masses(orders[n]["total"])
        other = cfg.pp_masses(density.scale(cfg.m**n))
        ok = ok and reference == other
    book.self_checks[f"{rid}: associativity (K1*R)*R == K1*(R*R) in exact rationals"] = bool(same)
    book.self_checks[f"{rid}: truncated powers == piecewise-polynomial convolution (orders 2, 3)"] = bool(ok)
    need(ok, f"{rid}: self-check of the box convolution failed")


def row_box(book, table, rows, rid, exact_polynomial=False, sens=False, check=False):
    cfg = Config(table, rows[rid])
    orders, masses, areas = box_values(cfg)
    if check:
        check_box_against_pp(book, cfg, orders, rid)
    put_standard(book, rid, cfg, masses, areas, box_sensitivity(cfg, cfg.N) if sens else None, exact_polynomial)
    return cfg, orders, masses, areas


def put_ratio(book, rid, cfg):
    r = cfg.m*cfg.areas["term_13a"]
    book.put(rid, "order_area_ratio", "interval_contains", r, klass="order_area_ratio",
             rule="the returned interval contains value (within u_r); disjoint is a failure; otherwise unresolved")


# ---------------------------------------------------------------- mp helpers for the tilted box (K_A2)
class Tilted:
    """13a = kappa*exp(-mu*T) on [a, b]; C_n(s) = (m*kappa)^n exp(-mu*s) IH_n(s)."""

    def __init__(self, cfg):
        piece, = cfg.spec["columns"]["term_13a"]
        self.cfg = cfg
        self.kappa, self.mu = hexfr(piece["amplitude"]), hexfr(piece["decay_per_eV"])
        self.a, self.b, self.m = cfg.a, cfg.b, cfg.length
        need(cfg.inp.get("multiplier", "weight_length") == "weight_length", "tilted box: multiplier")
        unit = TP.box(Fr(1), self.a, self.b)
        self.ih = {n: PP.from_tp(tp_power(unit, n)) for n in range(1, cfg.N + 1)}

    def area(self):
        mu = tomp(self.mu)
        return tomp(self.kappa)*(mp.exp(-mu*tomp(self.a)) - mp.exp(-mu*tomp(self.b)))/mu

    def first_moment(self):
        mu, a, b = tomp(self.mu), tomp(self.a), tomp(self.b)
        prim = lambda t: -mp.exp(-mu*t)*(t/mu + 1/mu**2)
        return tomp(self.kappa)*(prim(b) - prim(a))

    def laplace(self, tau):
        q = tomp(self.mu) + tomp(tau)
        return tomp(self.kappa)*(mp.exp(-q*tomp(self.a)) - mp.exp(-q*tomp(self.b)))/q

    @staticmethod
    def exp_poly(coeffs, mu, x0, x1):
        """int_{x0}^{x1} exp(-mu s) sum_k c_k s^k ds."""
        def prim(x):
            total = mp.mpf(0)
            for k, c in enumerate(coeffs):
                if c == 0:
                    continue
                inner = sum(mp.mpf(fact(k))/fact(k - j)*x**(k - j)/mu**(j + 1) for j in range(k + 1))
                total += tomp(c)*inner
            return -mp.exp(-mu*x)*total
        return prim(tomp(x1)) - prim(tomp(x0))

    def mass(self, n, l1, l2, pieces=None, factor=None):
        """Mass of C_n on the loss interval [l1, l2] (None is infinite)."""
        mu = tomp(self.mu)
        total = mp.mpf(0)
        for lo, hi, c in (self.ih[n].p if pieces is None else pieces):
            x0 = lo if l1 is None else max(lo, l1)
            x1 = hi if l2 is None else min(hi, l2)
            if x0 < x1:
                total += self.exp_poly(c, mu, x0, x1)
        return total*(tomp(self.m*self.kappa)**n if factor is None else factor)


def tilted_values(cfg, t):
    """masses, areas, sens (values and uncertainties) of the 13a chain of the tilted box."""
    val = {"masses": {}, "areas": {}, "sens": {}}
    unc = {"masses": {}, "areas": {}, "sens": {}}
    a1, ua1 = two_precisions(t.area)
    h1, uh1 = two_precisions(lambda: tomp(t.m*t.a*t.kappa)*mp.exp(-tomp(t.mu)*tomp(t.a)))
    for n in range(1, cfg.N + 1):
        w, uw = two_precisions(lambda: (tomp(t.m)*t.area())**n)
        h, uh = two_precisions(lambda: n*tomp(t.m*t.a*t.kappa)*mp.exp(-tomp(t.mu)*tomp(t.a))
                               * (tomp(t.m)*t.area())**(n - 1))
        val["areas"][n] = {c: (w if c in ("term_13a", "no_hole", "total") else Fr(0)) for c in COLUMNS}
        unc["areas"][n] = {c: (uw if c in ("term_13a", "no_hole", "total") else Fr(0)) for c in COLUMNS}
        val["sens"][n] = {c: (h if c in ("term_13a", "no_hole", "total") else Fr(0)) for c in COLUMNS}
        unc["sens"][n] = {c: (uh if c in ("term_13a", "no_hole", "total") else Fr(0)) for c in COLUMNS}
        if hasattr(cfg, "regions"):
            mv, mu_ = {}, {}
            for name, lo, hi in cfg.regions:
                value, u = Fr(0), Fr(0)
                for p, s in zip(cfg.px, cfg.ps):
                    l1, l2 = cfg.loss_interval(p, lo, hi)
                    if l2 is not None and l2 <= n*t.a or l1 is not None and l1 >= n*t.b:
                        continue
                    v, du = two_precisions(lambda: t.mass(n, l1, l2))
                    value, u = value + s*v, u + abs(s)*du
                mv[name], mu_[name] = value, u
            val["masses"][n] = {"total": mv, "no_hole": mv, "hole": {k: Fr(0) for k in mv}}
            unc["masses"][n] = {"total": mu_, "no_hole": mu_, "hole": {k: Fr(0) for k in mu_}}
    return val, unc, (a1, ua1, h1, uh1)


def tilted_norms(cfg, t):
    """Norms for the widths: the kernel is positive, so the L1 norms are the areas (rounded down)."""
    a1, _ = two_precisions(t.area)
    a1 = Fr(dec_down(a1))
    cfg.norms = {c: (a1 if c in ("term_13a", "no_hole", "total") else Fr(0)) for c in COLUMNS}
    cfg.m = t.m


def row_a2(book, table, rows):
    rid = "A2"
    cfg = Config(table, rows[rid])
    t = Tilted(cfg)
    tilted_norms(cfg, t)
    val, unc, _ = tilted_values(cfg, t)
    put_standard(book, rid, cfg, val["masses"], val["areas"], val["sens"], u=unc)
    # sum rule of the reference itself: all regions add up to S*w_n
    for n in range(1, cfg.N + 1):
        total = sum(val["masses"][n]["total"].values())
        need(abs(total - cfg.s_abs*val["areas"][n]["total"]) < Fr(1, 10**35), "A2: region sum")
    book.self_checks["A2: sum over the regions equals S*w_n within 1e-35 (orders 1..4)"] = True
    # Laplace transform and moments of the unmultiplied panels (column term_13a), prereg 4.6 item 9 and 7.2
    for n in (2, 3, 4):
        d_rel = n*cfg.eps + (n - 1)*cfg.eps_o
        floor = cfg.width("area", n, "term_13a")[0]/cfg.m**n
        for tau in ("0.05", "0.3"):
            v, u = two_precisions(lambda: t.laplace(Fr(float(tau)))**n)
            book.put(rid, f"panels.laplace[{n}][{tau}]", "direct_reference", v, u, d_abs=floor, d_rel=d_rel,
                     klass="panel_laplace", order=n, column="term_13a", tau_per_eV=float(tau),
                     unit="unmultiplied panel units, angstrom^-n")
        v, u = two_precisions(lambda: t.area()**n)
        book.put(rid, f"panels.moment[{n}][0]", "direct_reference", v, u, d_abs=floor, d_rel=d_rel,
                 klass="panel_moment", order=n, column="term_13a", moment=0)
        v, u = two_precisions(lambda: n*t.area()**(n - 1)*t.first_moment())
        book.put(rid, f"panels.moment[{n}][1]", "direct_reference", v, u, d_abs=floor*n*cfg.b, d_rel=d_rel,
                 klass="panel_moment", order=n, column="term_13a", moment=1)
    return cfg, t, val


def row_ct7(book, table, rows):
    rid = "CT7"
    cfg = Config(table, rows[rid])
    t = Tilted(cfg)
    tilted_norms(cfg, t)
    val, unc, _ = tilted_values(cfg, t)
    put_standard(book, rid, cfg, None, val["areas"], None, u=unc)


# ---------------------------------------------------------------- polynomial kernels (K_A4)
def legendre(d):
    p0, p1 = [Fr(1)], [Fr(0), Fr(1)]
    if d == 0:
        return p0
    for n in range(1, d):
        p0, p1 = p1, padd(pscale(pmul([Fr(0), Fr(1)], p1), Fr(2*n + 1, n + 1)), pscale(p0, Fr(-n, n + 1)))
    return p1


def a4_kernel(table, name):
    piece, = table["kernels"][name]["columns"]["term_13a"]
    lo, hi = (hexfr(x) for x in piece["interval_eV"])
    scale, d = hexfr(piece["scale"]), piece["degree"]
    need((lo, hi) == (2, 5) and piece["form"] == "scale*(1 + P_degree(u)), u = (2*T - 7)/3", "K_A4 form")
    # compose P_d(u) with u = (2 T - 7)/3
    u = [Fr(-7, 3), Fr(2, 3)]
    poly, power = [], [Fr(1)]
    for c in legendre(d):
        poly = padd(poly, pscale(power, c))
        power = pmul(power, u)
    return PP([(lo, hi, pscale(padd([Fr(1)], poly), scale))]), scale, d


def row_a4_point(book, table, rows, rid, name):
    k, scale, d = a4_kernel(table, name)
    c2 = pp_conv(k, k)
    value = c2(Fr(7))
    left = pval(c2.p[0][2], Fr(7))
    need(c2.breaks() == [4, 7, 10] and left == value, f"{rid}: order-2 pieces")
    expected = scale**2*3*(Fr(18, 19) if d == 9 else Fr(34, 33))
    need(value == expected, f"{rid}: the value at 7 eV is not the registered closed form")
    book.self_checks[f"{rid}: C_2(7 eV) == scale^2*3*{'18/19' if d == 9 else '34/33'} exactly"] = True
    for col in ("term_13a", "total"):
        book.put(rid, f"order_panels[2].value_at_7eV.{col}", "direct_reference", value, d_abs=Fr(1, 10**17),
                 d_rel=Fr(1, 10**13), klass="panel_point", order=2, column=col, at_eV=7.0,
                 unit="unmultiplied panel units, angstrom^-2 eV^-1; read at the end nodes of both panels at 7 eV")
    return value, scale, d


def row_a4b(book, table, rows):
    rid = "A4b"
    cfg = Config(table, rows[rid])
    k, scale, d = a4_kernel(table, cfg.kname)
    cfg.m = cfg.length
    area = pval(k.antider().tail, 0)
    # 1 + P_d >= 0: the L1 norm is the area
    cfg.norms = {c: (area if c in ("term_13a", "no_hole", "total") else Fr(0)) for c in COLUMNS}
    dens = {1: k}
    for n in range(2, cfg.N + 1):
        dens[n] = pp_conv(dens[n - 1], k)
    c3_other = pp_conv(k, pp_conv(k, k))
    need(dens[3].p == c3_other.p, "A4b: associativity")
    book.self_checks["A4b: (K*K)*K == K*(K*K) in exact rationals (degree 29)"] = True
    masses, areas = {}, {}
    for n, f in dens.items():
        mv = cfg.pp_masses(f.scale(cfg.m**n))
        zero = {name: Fr(0) for name in mv}
        masses[n] = {"total": mv, "no_hole": mv, "hole": zero}
        w = (cfg.m*area)**n
        need(sum(mv.values()) == cfg.s_abs*w, "A4b: region sum")
        areas[n] = {c: (w if c in ("term_13a", "no_hole", "total") else Fr(0)) for c in COLUMNS}
    book.self_checks["A4b: sum over the regions == S*w_n exactly (orders 1..3)"] = True
    put_standard(book, rid, cfg, masses, areas)


# ---------------------------------------------------------------- Lorentzian (K_A5)
class Lorentz:
    def __init__(self, cfg):
        piece, = cfg.spec["columns"]["term_13a"]
        self.kappa, self.c, self.g = (hexfr(piece[k]) for k in ("amplitude", "center_eV", "half_width_eV"))
        self.a, self.b, self.m = cfg.a, cfg.b, cfg.length

    def setup(self, level, step):
        """Nodes of the composite Gauss-Legendre rule on [2a, 2b] with panels of width ``step``."""
        from mpmath.calculus.quadrature import GaussLegendre
        rule = GaussLegendre(mp).calc_nodes(level, mp.prec)
        self.kap, self.cc, self.gg = tomp(self.kappa), tomp(self.c), tomp(self.g)
        self.aa, self.bb = tomp(self.a), tomp(self.b)
        self.atan_a = mp.atan((self.aa - self.cc)/self.gg)
        self.area = self.cum(self.bb)
        lo, hi = 2*self.a, 2*self.b
        count = int((hi - lo)/step)
        need(lo + count*step == hi, "panel grid")
        self.panels = []
        for i in range(count):
            x0, x1 = lo + i*step, lo + (i + 1)*step
            mid, half = (tomp(x0) + tomp(x1))/2, (tomp(x1) - tomp(x0))/2
            nodes = [(mid + half*x, half*w) for x, w in rule]
            self.panels.append((x0, x1, [(x, w, self.c2(x)) for x, w in nodes]))

    def cum(self, y):
        """Cumulative of the unmultiplied kernel (mp)."""
        y = min(max(y, self.aa), self.bb)
        return self.kap*self.gg*(mp.atan((y - self.cc)/self.gg) - self.atan_a)

    def kernel(self, t):
        return self.kap/(1 + ((t - self.cc)/self.gg)**2) if self.aa <= t <= self.bb else mp.mpf(0)

    def c2(self, x):
        """(K*K)(x) in closed form (partial fractions with complex logarithms)."""
        lo, hi = max(self.aa, x - self.bb), min(self.bb, x - self.aa)
        if hi <= lo:
            return mp.mpf(0)
        z0 = mp.mpc(self.cc, self.gg)
        w = mp.mpc(x - self.cc, -self.gg)
        wbar = mp.conj(w)
        d1 = 2*self.cc - x
        if d1 == 0:
            i1 = 1/(hi - z0) - 1/(lo - z0)
        else:
            i1 = -((mp.log(hi - z0) - mp.log(hi - wbar)) - (mp.log(lo - z0) - mp.log(lo - wbar)))/d1
        i2 = -((mp.log(hi - z0) - mp.log(hi - w)) - (mp.log(lo - z0) - mp.log(lo - w)))/(z0 - w)
        return (self.kap*self.gg)**2*mp.re(i1 - i2)/2

    def mass2(self, l1, l2):
        total = mp.mpf(0)
        for x0, x1, nodes in self.panels:
            if (l1 is None or x0 >= l1) and (l2 is None or x1 <= l2):
                total += sum(w*v for _, w, v in nodes)
            else:
                need((l2 is not None and x0 >= l2) or (l1 is not None and x1 <= l1), "bin edge inside a panel")
        return total

    def mass3(self, l1, l2):
        """int C2(x) [A(l2 - x) - A(l1 - x)] dx; the kinks x = l - a, l - b are panel ends."""
        total = mp.mpf(0)
        for x0, x1, nodes in self.panels:
            for l in (l1, l2):
                if l is not None:
                    for kink in (l - self.a, l - self.b):
                        need(not x0 < kink < x1, "kink inside a panel")
            for x, w, v in nodes:
                upper = self.area if l2 is None else self.cum(tomp(l2) - x)
                lower = mp.mpf(0) if l1 is None else self.cum(tomp(l1) - x)
                total += w*v*(upper - lower)
        return total


def row_a5(book, table, rows):
    rid = "A5"
    cfg = Config(table, rows[rid])
    lz = Lorentz(cfg)
    cfg.m = lz.m
    settings = ((DPS_LOW, 4, Fr(1, 4)), (DPS_HIGH, 5, Fr(1, 4)), (DPS_LOW, 4, Fr(1, 8)))
    runs = []
    for dps, level, step in settings:
        with mp.workdps(dps):
            lz.setup(level, step)
            run = {"area": mpfr(lz.area), "masses": {}}
            for n in (1, 2, 3):
                factor = tomp(lz.m)**n
                mv = {}
                for name, lo, hi in cfg.regions:
                    value = Fr(0)
                    for p, s in zip(cfg.px, cfg.ps):
                        l1, l2 = cfg.loss_interval(p, lo, hi)
                        if l2 is not None and l2 <= n*lz.a or l1 is not None and l1 >= n*lz.b:
                            continue
                        if n == 1:
                            x = (lz.area if l2 is None else lz.cum(tomp(l2))) - (0 if l1 is None else lz.cum(tomp(l1)))
                        else:
                            x = (lz.mass2 if n == 2 else lz.mass3)(l1, l2)
                        value += s*mpfr(factor*x)
                    mv[name] = value
                run["masses"][n] = mv
            if dps == DPS_HIGH:
                # the closed form of K*K against direct quadrature at a few points (independent check)
                worst = mp.mpf(0)
                for x in ("5.1", "7.3", "10.6", "13.9", "17.3", "22.25"):
                    x = mp.mpf(x)
                    lo, hi = max(lz.aa, x - lz.bb), min(lz.bb, x - lz.aa)
                    cuts = sorted({lo, hi} | {p for p in (lz.cc, x - lz.cc) if lo < p < hi})
                    direct = mp.quad(lambda t: lz.kernel(t)*lz.kernel(x - t), cuts)
                    worst = max(worst, abs(direct - lz.c2(x))/direct)
                need(worst < mp.mpf(10)**-25, f"A5: closed form of K*K ({worst})")
                book.self_checks["A5: closed form of K*K equals direct quadrature at 6 points (relative 1e-25)"] = True
            runs.append(run)
    base, high, fine = runs
    area = high["area"]
    cfg.norms = {c: (Fr(dec_down(area)) if c in ("term_13a", "no_hole", "total") else Fr(0)) for c in COLUMNS}
    masses, um, areas, ua = {}, {}, {}, {}
    worst_sum = Fr(0)
    for n in (1, 2, 3):
        mv = high["masses"][n]
        mu_ = {k: 4*max(abs(mv[k] - base["masses"][n][k]), abs(mv[k] - fine["masses"][n][k]))
               + abs(mv[k])*Fr(1, 10**(DPS_LOW - 2)) if mv[k] != 0 else Fr(0) for k in mv}
        zero = {k: Fr(0) for k in mv}
        masses[n] = {"total": mv, "no_hole": mv, "hole": zero}
        um[n] = {"total": mu_, "no_hole": mu_, "hole": zero}
        w = (lz.m*area)**n
        uw = abs(w)*Fr(1, 10**(DPS_LOW - 2)) + 4*n*abs(w)*abs(area - base["area"])/area
        areas[n] = {c: (w if c in ("term_13a", "no_hole", "total") else Fr(0)) for c in COLUMNS}
        ua[n] = {c: (uw if c in ("term_13a", "no_hole", "total") else Fr(0)) for c in COLUMNS}
        worst_sum = max(worst_sum, abs(sum(mv.values()) - cfg.s_abs*w)/w)
    need(worst_sum < Fr(1, 10**30), f"A5: region sum {float(worst_sum)}")
    book.self_checks["A5: sum over the regions equals S*w_n (relative 1e-30, orders 1..3)"] = True
    book.note(rid, "method", "order 1: arctangent; order 2: closed form of K*K (complex logarithms) integrated by a "
              "composite Gauss-Legendre rule; order 3: int C2(x)[A(l2-x) - A(l1-x)] dx with the kinks on panel ends; "
              f"settings (digits, nodes per panel, panel width eV): ({DPS_LOW}, 24, 0.25), ({DPS_HIGH}, 48, 0.25), "
              f"({DPS_LOW}, 24, 0.125); u_r = 4*max difference + precision floor")
    put_standard(book, rid, cfg, masses, areas, u={"masses": um, "areas": ua})


# ---------------------------------------------------------------- R-INV (renewal equation for the cumulative)
def inv_solution(cfg, N, with_tail=True, s_max=40):
    """H(s) = int_{-inf}^{s} J^(N) in the loss coordinate for a delta primary of strength S at loss 0.

    (I - A_R) J = (I + A_gbar) F - A_{K1*R^N} F with R = c_R*1[alpha, beta]:
    H(s) = S*(1 + Gbar(s) - T(s)) + c_R*(H2(s - alpha) - H2(s - beta)), H2' = H, solved from the known side."""
    m = cfg.m
    R = cfg.K["term_13a"]
    need(sorted(R.t) == sorted(TP.box(R(R.breaks()[0]), *R.breaks()).t), "R-INV needs a single box R")
    alpha, beta = R.breaks()
    c_r = m*R(alpha)
    S = cfg.s_abs
    gbar = (cfg.K["total"] - R).scale(m)
    rhs_tp = gbar.antider()
    if with_tail:
        tail = cfg.K["total"].scale(m).conv(tp_power(R.scale(m), N))
        rhs_tp = rhs_tp - tail.antider()
    step = Fr(1, 2)
    grid = [step*i for i in range(int(s_max/step) + 1)]
    need(all(b in grid for b in rhs_tp.breaks()), "R-INV: breaks off the grid")
    need(alpha/step == int(alpha/step) and beta/step == int(beta/step), "R-INV: delay off the grid")
    rhs = PP.from_tp(rhs_tp, grid)
    rhs_poly = {lo: c for lo, hi, c in rhs.p}
    H, H2 = {}, {}        # polynomial on [x, x + step)
    h2_value = Fr(0)      # H2 at the left end of the current interval
    for x in grid[:-1]:
        poly = padd([S], pscale(rhs_poly.get(x, rhs.tail if x >= rhs.end() else []), S))
        for delay, sign in ((alpha, 1), (beta, -1)):
            src = x - delay
            if src >= 0:
                poly = padd(poly, pscale(pshift(H2[src], delay), sign*c_r))
        H[x] = poly
        q = pint(poly)
        H2[x] = padd(q, [h2_value - pval(q, x)])
        h2_value = pval(H2[x], x + step)

    def cumulative(s):
        if s < 0:
            return Fr(0)
        need(s < s_max, "R-INV: beyond the solved range")
        x = step*int(s/step)
        return pval(H[x], s)
    return cumulative, H


def row_inv(book, table, rows, rid):
    cfg = Config(table, rows[rid])
    N = cfg.N
    need(cfg.primary_kind == "delta" and cfg.px == [1000] and cfg.ps == [1], "INV primary")
    cumulative, H = inv_solution(cfg, N)
    last = sorted(H)[-1]
    need(len(H[last]) <= 1 and H[last] == H[sorted(H)[-8]], "R-INV: the cumulative is not constant at the end")
    total = pval(H[last], 0)
    solution = cfg.masses(cumulative, total)
    # self-check with a different function: the order sum by truncated powers
    orders = b1_orders(cfg)
    by_orders = {name: Fr(0) for name, _, _ in cfg.regions}
    for n in orders:
        for name, v in cfg.tp_masses(orders[n]["total"]).items():
            by_orders[name] += v
    zero_loss = cfg.masses(lambda s: cfg.s_abs if s >= 0 else Fr(0), cfg.s_abs)
    for name in by_orders:
        by_orders[name] += zero_loss[name]
    need(solution == by_orders, f"{rid}: renewal solution differs from the order sum")
    book.self_checks[f"{rid}: renewal solution == J0 + sum of the orders (separate function), exactly"] = True
    d_abs, d_rel = cfg.spectrum_width()
    for name, _, _ in cfg.regions:
        book.put(rid, f"spectrum.{name}", "direct_reference", solution[name], d_abs=d_abs, d_rel=d_rel,
                 klass="spectrum", region=name)
    # without the tail term: the full series; equal below the threshold, outside the width on the reaching side
    wrong_cumulative, Hw = inv_solution(cfg, N, with_tail=False)
    threshold = cfg.a + N*(cfg.K["term_13a"].breaks()[0])
    reaching, controls, outside = [], [], []
    for name, lo, hi in cfg.regions:
        l1, l2 = cfg.loss_interval(cfg.px[0], lo, hi)
        s1 = Fr(-1) if l1 is None else l1
        if l2 is None:
            wrong = None            # the full series has no finite total on the solved range
        else:
            wrong = wrong_cumulative(l2) - (wrong_cumulative(s1) if s1 >= 0 else Fr(0))
        reaches = l2 is None or l2 > threshold
        if not reaches:
            need(wrong == solution[name], f"{rid}: the tail term changes a bin below the threshold")
            controls.append(name)
            continue
        reaching.append(name)
        if wrong is not None and abs(wrong - solution[name]) > d_abs + d_rel*abs(solution[name]):
            outside.append({"region": name, "without_tail": dec(wrong)[0], "reference": dec(solution[name])[0],
                            "difference": dec(wrong - solution[name])[0],
                            "width": dec_down(d_abs + d_rel*abs(solution[name]))})
    need(outside, f"{rid}: no reaching-side bin is outside the width without the tail term")
    book.note(rid, "reference: reaching-side bin outside the width without the tail term",
              {"kind": "diagnostic", "oracle": "R-INV", "threshold_loss_eV": float(threshold),
               "normal_controls_exactly_equal": controls, "reaching_side": reaching,
               "outside_width_without_tail": outside, "passes": True})
    return cfg, solution, orders


# ---------------------------------------------------------------- reference-only models
def ref_models(cfg, model, N):
    """Densities of the three reference-only models for total and no_hole (hole = total - no_hole)."""
    m = cfg.m
    R = cfg.K["term_13a"].scale(m)
    out = {n: {} for n in range(1, N + 1)}
    for comp in ("total", "no_hole"):
        K = cfg.K[comp].scale(m)
        if model in ("full_kernel_geometric_ref", "full_kernel_factorial_ref"):
            f = K
            for n in range(1, N + 1):
                out[n][comp] = f.scale(Fr(1, fact(n))) if model.endswith("factorial_ref") else f
                f = f.conv(K)
        else:
            need(model == "excess_poisson_ref" and N <= 3, model)
            g = K - R
            g2, r2 = g.conv(g), R.conv(R)
            out[1][comp] = g + R
            out[2][comp] = g2.scale(Fr(1, 2)) + g.conv(R) + r2
            if N >= 3:
                out[3][comp] = (g2.conv(g).scale(Fr(1, 6)) + g2.conv(R).scale(Fr(1, 2)) + g.conv(r2) + r2.conv(R))
    for n in out:
        out[n]["hole"] = out[n]["total"] - out[n]["no_hole"]
    return out


def row_ref(book, table, rows, rid):
    cfg = Config(table, rows[rid])
    need(cfg.charge == 1, "the reference entry takes hole charge 1 only")
    dens = ref_models(cfg, cfg.inp["model"], cfg.N)
    masses = {n: {c: cfg.tp_masses(f[c]) for c in COMPONENTS} for n, f in dens.items()}
    put_standard(book, rid, cfg, masses, None)
    return cfg, dens, masses


# ---------------------------------------------------------------- truncation bounds
def tau_grid(a):
    """tau_i = 10^(-3 + i/20) per eV while tau_i <= 40/a (prereg 3.3)."""
    out, i = [], 0
    limit = mp.mpf(40)/tomp(a)
    while True:
        tau = mp.power(10, mp.mpf(-3) + mp.mpf(i)/20)
        if tau > limit*(1 + mp.mpf(10)**-30):
            return out
        out.append((i, tau))
        i += 1


def box_lambda(c, a, b, tau):
    return c*(mp.exp(-tau*a) - mp.exp(-tau*b))/tau


def bound_grid(s_abs, U, N, k1, rr, a, panels=None, variant=None):
    """Grid minimum of S e^(tau U) Lambda_K1 Lambda_R^N/(1 - Lambda_R) for boxes on [a, b].

    k1, rr: (height with the multiplier, b, error per eV with the multiplier); ``panels``: first-order panel
    edges for the error terms exp(-tau t_lo) e_panel.  Returns (value, tau or 0 for the limit, grid index)."""
    S, Um = tomp(s_abs), tomp(U)
    best = None

    def lam(spec, tau):
        c, b, e = spec
        value = box_lambda(tomp(c), tomp(a), tomp(b), tau) if tau else tomp(c)*(tomp(b) - tomp(a))
        if panels and e:
            for lo, hi in zip(panels, panels[1:]):
                value += (mp.exp(-tau*tomp(lo)) if tau else 1)*tomp(e)*(tomp(hi) - tomp(lo))
        return value
    candidates = list(tau_grid(a))
    if lam(rr, 0) < 1:
        candidates.append((None, mp.mpf(0)))
    power = N + (1 if variant == "N+1" else -1 if variant == "N-1" else 0)
    for index, tau in candidates:
        lr = lam(rr, tau)
        if lr >= 1:
            continue
        growth = 1 if variant == "no_exp" else mp.exp(tau*Um)
        value = S*growth*lam(k1, tau)*lr**power/(1 - lr)
        if best is None or value < best[0]:
            best = (value, tau, index)
    return best


def ih_cdf(n, x):
    """P(sum of n uniform(0,1) <= x)."""
    if x <= 0:
        return Fr(0)
    if x >= n:
        return Fr(1)
    return sum(((-1)**k*math.comb(n, k)*(x - k)**n for k in range(int(x) + 1)), Fr(0))/fact(n)


def t2_bounds(book, table, rows, rid, ledger_done):
    cfg = Config(table, rows[rid])
    N, m, a, b = cfg.N, cfg.m, cfg.a, cfg.b
    ctot, c13 = cfg.K["total"](a), cfg.K["term_13a"](a)
    need(cfg.K["total"].breaks() == [a, b] and cfg.K["term_13a"].breaks() == [a, b], "T2 kernel")
    p_max = max(p for p, s in zip(cfg.px, cfg.ps) if s > 0)
    out = {}
    for i, (lo, hi) in enumerate(zip(cfg.edges, cfg.edges[1:])):
        U = p_max - lo
        if U < (N + 1)*a:
            out[i] = None
            continue
        (v, u) = two_precisions(lambda: bound_grid(cfg.s_abs, U, N, (m*ctot, b, 0), (m*c13, b, 0), a)[0])
        with mp.workdps(DPS_HIGH):
            _, tau, index = bound_grid(cfg.s_abs, U, N, (m*ctot, b, 0), (m*c13, b, 0), a)
        need(index is not None, f"{rid}: the minimum is the tau -> 0 limit")
        # exact cumulative neglected quantity: orders n > N with loss <= U (identical boxes: Irwin-Hall)
        h = b - a
        exact = sum((cfg.s_abs*m**n*ctot*c13**(n - 1)*h**n*ih_cdf(n, (U - n*a)/h)
                     for n in range(N + 1, int(U/a) + 1)), Fr(0))
        out[i] = dict(U=U, grid=v, u=u, tau=tau, index=index, exact=exact)
        book.put(rid, f"truncation.bound.bin[{i}].grid_minimum", "bound_reference", v, u, klass="bound_grid_minimum",
                 bin=i, U_D_eV=float(U), tau_index=index, tau_per_eV=mp.nstr(tau, 20), lower_factor="1 - 1e-9",
                 upper_factor="1 + 1e-6",
                 rule="two verdicts: bound >= value*(1 - 1e-9) - u_r and bound <= value*(1 + 1e-6) + u_r; tau > 0")
        book.put(rid, f"truncation.bound.bin[{i}].covers_exact", "bound_reference", exact,
                 klass="bound_covers_exact", bin=i, U_D_eV=float(U),
                 rule="bound >= value (exact cumulative neglected quantity at or above the region)")
        if i in (0, 1):
            book.put(rid, f"truncation.bound.bin[{i}].ratio_le_10", "bound_reference", exact,
                     klass="bound_ratio", bin=i, factor=10, rule="bound <= 10*value")
            need(mpfr(mp.mpf(1)) and v <= 10*exact and v >= exact, f"{rid}: the reference bound itself breaks the rule")
    return cfg, out


# ---------------------------------------------------------------- main
def build():
    table = load_table()
    ledger = load_ledger()
    rows = {r["id"]: r for r in table["rows"]}
    book = Book(ledger)
    keep = {}

    # ---- RM-A1 and the other exact box rows
    cfg, orders, masses, areas = row_box(book, table, rows, "A1", exact_polynomial=True)
    keep["A1"] = (cfg, orders, masses, areas)
    need(any(masses[n]["total"]["outside.high"] > 0 for n in masses)
         and any(masses[n]["total"]["outside.low"] > 0 for n in masses), "A1: outside regions")
    book.note("A1", "outside_regions", {"order[1].outside.high.total > 0": masses[1]["total"]["outside.high"] > 0,
                                        "order[6].outside.low.total > 0": masses[6]["total"]["outside.low"] > 0})
    keep["A1s"] = row_box(book, table, rows, "A1s", exact_polynomial=True, check=True)
    row_box(book, table, rows, "A1p", exact_polynomial=True)
    keep["T2c"] = row_box(book, table, rows, "T2c", exact_polynomial=True)
    keep["DS1"] = row_box(book, table, rows, "DS1", exact_polynomial=True)
    for rid in ("A3a", "A3b", "A3c", "A3d", "A3e"):
        keep[rid] = row_box(book, table, rows, rid, check=rid == "A3a")
        put_ratio(book, rid, keep[rid][0])
    keep["A3z"] = row_box(book, table, rows, "A3z")
    row_box(book, table, rows, "OI1")
    row_box(book, table, rows, "OI2")
    row_box(book, table, rows, "REF4a")
    row_box(book, table, rows, "CT1")
    row_box(book, table, rows, "I-areas", sens=True)

    # A1s: the support statement of the row (7, 6, 6 regions without overlap, three components)
    cfg, orders, masses, areas = keep["A1s"]
    zero_regions = {n: [name for name in masses[n]["total"] if masses[n]["total"][name] == 0] for n in masses}
    need([len(zero_regions[n]) for n in (1, 2, 3)] == [7, 6, 6], "A1s: regions outside the support")
    need(masses[2]["total"]["bin[2]"] > 0 and masses[1]["total"]["bin[2]"] == 0 == masses[3]["total"]["bin[2]"],
         "A1s: bin[2]")
    book.note("A1s", "regions_without_overlap", zero_regions)

    # A3z: moments of the hole column of the unmultiplied panels
    cfg = keep["A3z"][0]
    a_h, a_r = cfg.areas["hole"], cfg.areas["term_13a"]
    m1_h, m1_r = cfg.K["hole"].first_moment_constant(), cfg.K["term_13a"].first_moment_constant()
    for n in (2, 3):
        d_rel = n*cfg.eps + (n - 1)*cfg.eps_o
        floor = cfg.width("area", n, "hole")[0]/cfg.m**n
        book.put("A3z", f"panels.moment[{n}][0].hole", "direct_reference", a_h*a_r**(n - 1), d_abs=floor, d_rel=d_rel,
                 klass="panel_moment", order=n, column="hole", moment=0)
        book.put("A3z", f"panels.moment[{n}][1].hole", "direct_reference",
                 m1_h*a_r**(n - 1) + (n - 1)*a_h*m1_r*a_r**(n - 2), d_abs=floor*n*cfg.b, d_rel=d_rel,
                 klass="panel_moment", order=n, column="hole", moment=1)
    book.note("A3z", "hole_area_of_the_binary_inputs", dec(a_h)[0])
    need(any(v != 0 for v in keep["A3z"][2][2]["hole"].values()), "A3z: hole bins")

    # A3a: negative hole bins exist
    masses = keep["A3a"][2]
    negative = sorted(f"order[{n}].{name}.hole" for n in masses for name, v in masses[n]["hole"].items() if v < 0)
    need(negative, "A3a: no negative hole bin")
    book.note("A3a", "negative_hole_keys", negative)
    cfg = keep["A3a"][0]
    book.note("A3a", "areas_per_angstrom", {"A_13a": dec(cfg.areas["term_13a"])[0], "A_total": dec(cfg.areas["total"])[0],
                                            "l1_total": dec(cfg.norms["total"])[0]})

    # ---- T2 bounds, neglected totals
    for rid in ("T2a", "T2c"):
        cfg, bounds = t2_bounds(book, table, rows, rid, None)
        keep[rid + ".bounds"] = (cfg, bounds)
    cfg = keep["T2c"][0]
    w, r = cfg.m*cfg.areas["total"], cfg.m*cfg.areas["term_13a"]
    book.put("T2c", "neglected_total_area value S*w*r^4/(1-r)", "interval_contains", cfg.s_abs*w*r**4/(1 - r),
             klass="neglected_total_area", rule="the returned interval contains value")
    cfg4 = Config(table, rows["T4v"])
    w, r = cfg4.m*cfg4.areas["total"], cfg4.m*cfg4.areas["term_13a"]
    value = cfg4.s_abs*w*r**3/(1 - r)
    need(value == 1024*(1 - Fr(1, 1024))**4, "T4v closed form")
    book.put("T4v", "neglected_total_area interval contains S*w*r^3/(1-r) = 1024*(1-2^-10)^4", "interval_contains",
             value, klass="neglected_total_area", rule="the returned interval contains value")

    # ---- T2e: bounds with the error terms of the first-order panels
    row = rows["T2e"]
    cfg = Config(table, row)
    err = {k: hexfr(v) for k, v in row["inputs"]["inject"]["errors"].items()}
    panels = [Fr(x) for x in (1, 2, 4, 8, 16, 32, 50)]
    a, b, m, N = cfg.a, cfg.b, cfg.m, cfg.N
    ctot, c13 = cfg.K["total"](a), cfg.K["term_13a"](a)
    differ = []
    for i, (lo, hi) in enumerate(zip(cfg.edges, cfg.edges[1:])):
        U = max(cfg.px) - lo
        if U < (N + 1)*a:
            continue
        args = (cfg.s_abs, U, N, (m*ctot, b, m*err["total"]), (m*c13, b, m*err["term_13a"]), a)
        v, u = two_precisions(lambda: bound_grid(*args, panels=panels)[0])
        plain, _ = two_precisions(lambda: bound_grid(*args)[0])
        with mp.workdps(DPS_HIGH):
            _, tau, index = bound_grid(*args, panels=panels)
        book.put("T2e", f"truncation.bound.bin[{i}].grid_minimum_with_error_terms", "bound_reference", v, u,
                 klass="bound_grid_minimum", bin=i, U_D_eV=float(U), tau_index=index,
                 tau_per_eV=mp.nstr(tau, 20), lower_factor="1 - 1e-9", upper_factor="1 + 1e-6",
                 rule="two verdicts: bound >= value*(1 - 1e-9) - u_r and bound <= value*(1 + 1e-6) + u_r",
                 panel_edges_eV=[float(x) for x in panels],
                 error_per_eV={"total": float(err["total"]), "term_13a": float(err["term_13a"])})
        outside = not (v*(1 - Fr(1, 10**9)) <= plain <= v*(1 + Fr(1, 10**6)))
        differ.append({"bin": i, "with_error_terms": dec(v)[0], "without_error_terms": dec(plain)[0],
                       "ratio_without_over_with": dec(plain/v)[0], "outside_width": outside})
        book.wrong_answer("X29", "T2e", f"truncation.bound.bin[{i}].grid_minimum_with_error_terms", plain, v,
                          0, Fr(1, 10**6), "Lambda from point values only (no error terms)",
                          rule="below value*(1 - 1e-9)")
    need(all(d["outside_width"] for d in differ), "T2e: the error terms do not move the bound")
    book.note("T2e", "the same bins without the error terms differ by more than the width (recorded by the reference)",
              {"kind": "diagnostic", "oracle": "R-A", "bins": differ, "passes": True})

    # ---- X18: wrong bounds of T2c
    cfg, bounds = keep["T2c.bounds"]
    a, b, m, N = cfg.a, cfg.b, cfg.m, cfg.N
    ctot, c13 = cfg.K["total"](a), cfg.K["term_13a"](a)
    for i in (0, 1):
        for variant, text in (("no_exp", "exp(tau*U) dropped"), ("N+1", "power N -> N+1"), ("N-1", "power N -> N-1")):
            with mp.workdps(DPS_HIGH):
                wrong = mpfr(bound_grid(cfg.s_abs, bounds[i]["U"], N, (m*ctot, b, 0), (m*c13, b, 0), a,
                                        variant=variant)[0])
            ref = bounds[i]["grid"]
            outside = not (ref*(1 - Fr(1, 10**9)) <= wrong <= ref*(1 + Fr(1, 10**6)))
            book.wrong.append({"mutation": "X18", "row": "T2c", "key": f"truncation.bound.bin[{i}].grid_minimum",
                               "wrong_value": dec(wrong)[0], "reference": dec(ref)[0],
                               "difference": dec(wrong - ref)[0], "width": "factors 1 - 1e-9 and 1 + 1e-6",
                               "outside_width": bool(outside), "rule": "outside [value*(1 - 1e-9), value*(1 + 1e-6)]",
                               "note": text})

    # ---- REF4b, REF4c: exact neglected quantity of the model
    for rid in ("REF4b", "REF4c"):
        cfg = Config(table, rows[rid])
        N, a = cfg.N, cfg.a
        K = cfg.K["total"].scale(cfg.m)
        factorial = cfg.inp["model"].endswith("factorial_ref")
        p_max = max(cfg.px)
        for i, (lo, hi) in enumerate(zip(cfg.edges, cfg.edges[1:])):
            U = p_max - lo
            reached = U >= (N + 1)*a
            exact, n, f = Fr(0), N + 1, tp_power(K, N + 1)
            while n*a <= U:
                exact += cfg.s_abs*f.antider()(U)/(fact(n) if factorial else 1)
                n, f = n + 1, f.conv(K)
            book.put(rid, f"truncation.bound.bin[{i}].positive_where_reached_and_covers_exact", "bound_reference",
                     exact, klass="bound_model", bin=i, U_D_eV=float(U), reached=reached,
                     rule="reached (U_D >= (N+1)a, equality on the reached side): bound > 0 and bound >= value; "
                          "not reached: (0.0, 'complete', None)")

    # ---- A2, CT7 (tilted box), A4 (polynomials), A5 (Lorentzian)
    keep["A2"] = row_a2(book, table, rows)
    row_ct7(book, table, rows)
    keep["A4a"] = row_a4_point(book, table, rows, "A4a", "K_A4_9")
    keep["A4c"] = row_a4_point(book, table, rows, "A4c", "K_A4_16")
    row_a4b(book, table, rows)
    row_a5(book, table, rows)

    # ---- INV
    keep["INV1"] = row_inv(book, table, rows, "INV1")
    keep["INV2"] = row_inv(book, table, rows, "INV2")

    # ---- reference-only models
    for rid in ("REF1", "REF2", "REF3"):
        keep[rid] = row_ref(book, table, rows, rid)
    cfg, dens, masses = keep["REF3"]
    R = cfg.K["term_13a"].scale(cfg.m)
    b1 = keep["A3a"][2]
    need(Config(table, rows["A3a"]).m == cfg.m, "REF3 and A3a use the same multiplier")
    for comp in ("total", "no_hole"):
        g = cfg.K[comp].scale(cfg.m) - R
        half = cfg.tp_masses(g.conv(g).scale(Fr(1, 2)))
        d_abs, d_rel = cfg.width("bin", 2, comp)
        for name, _, _ in cfg.regions:
            need(half[name] == masses[2][comp][name] - b1[2][comp][name], "REF3 - A3a")
            book.put("REF3", f"order[2].{name}.{comp} minus A3a.result", "direct_reference", half[name], d_abs=d_abs,
                     d_rel=d_rel, klass="region_difference", order=2, region=name, component=comp,
                     minus="A3a.result", rule="width of prereg 7.2 for order 2 applied to the difference")
    book.self_checks["REF3: order 2 minus B1 order 2 == gbar*gbar/2 exactly (total, no_hole)"] = True
    for n in (1, 2, 3):
        for comp in COMPONENTS:
            for name in keep["REF1"][2][n][comp]:
                need(keep["REF2"][2][n][comp][name]*fact(n) == keep["REF1"][2][n][comp][name], "REF2 = REF1/n!")

    # ---- error parts (A3f, AE1)
    error_parts(book, table, rows)

    # ---- known wrong answers
    known_wrong(book, table, rows, keep)
    return table, ledger, book


def error_parts(book, table, rows):
    # A3f: Young right-hand side of order 2 with exact L1 norms (S2, K_A3, error 1e-3 on all 11 columns)
    row = rows["A3f"]
    cfg = Config(table, row)
    e_point = hexfr(row["inputs"]["inject"]["errors"]["all 11 columns"])
    width = cfg.b - cfg.a
    e_area = e_point*width                      # per angstrom, interpolation-only injection: (b - a)*e
    m = cfg.m
    n_r = cfg.norms["term_13a"]
    c2 = pp_conv(PP.from_tp(cfg.K["total"]), PP.from_tp(cfg.K["term_13a"]))
    book.put("A3f", "single_loss.window_area.total", "direct_reference", cfg.areas["total"], d_abs=0,
             d_rel=Fr(1, 10**12), klass="window_area", unit="per angstrom")
    book.put("A3f", "single_loss.window_area_error.total", "direct_reference", e_area, klass="window_area_error",
             interval=[dec(e_area)[0], dec(e_area*(1 + Fr(1, 10**9)))[0]], unit="per angstrom",
             rule="value <= reported error <= value*(1 + 1e-9)")
    for col in ("total", "no_hole", "hole", "term_13a"):
        nominal = m*m*(n_r*e_area + cfg.norms[col]*e_area + e_area*e_area)
        book.put("A3f", f"order[2].area_error_parts[inherited][{col}]", "error_part", nominal, klass="young_area",
                 column=col, norm_R_per_A=dec(n_r)[0], norm_column_per_A=dec(cfg.norms[col])[0],
                 nominal_area_error_per_A=dec(e_area)[0], nominal_multiplier_A=dec(m)[0],
                 rule="rebuilt = m*nR*e_j + m*nj*e_R + e_j*e_R with the reported binary64 multiplier m and the "
                      "reported seeds e_j, e_R = first_order.fixed_multiplier_area_error (dimensionless), exact "
                      "rationals; the product value is not below rebuilt and not above rebuilt*(1 + 1e-9) + "
                      "8*eps*rebuilt; value is the nominal one for m = 1/A_total and e = (b - a)*1e-3 exactly")
    book.note("A3f", "norms_for_order_3", {
        "norm_R_per_A": dec(n_r)[0], "norm_C2_total_per_A2": dec(c2.l1_linear())[0],
        "norm_K_total_per_A": dec(cfg.norms["total"])[0],
        "text": "unmultiplied exact L1 norms; ||C_2||_1 is the integral of |K_total * K_13a| (piecewise linear, "
                "split at its zeros)"})
    need(c2.l1_linear() >= abs(pval(c2.antider().tail, 0)), "A3f: L1 norm below the algebraic area")
    # control of the point path: K_A1s, e = 1e-3 on 13a, no_hole and total, s = 7 eV, unmultiplied
    control = row["inputs"]["controls"]["control_point_path"]
    e = hexfr(control["inject"]["errors"]["total"])
    kc = Config(table, row, kernel=control["kernel"], multiplier="weight_length")
    value_k = kc.K["total"](kc.a)
    width_k = kc.b - kc.a
    need((kc.a, kc.b) == (2, 5), "control kernel")
    point = width_k*(value_k*e + e*value_k + e*e)
    need(abs(point - Fr(603, 10**6)) < Fr(1, 10**15), "A3f control value")
    book.put("A3f", "control_point_path.order[2].node_parts[inherited][total]@7eV", "error_part", point,
             d_abs=Fr(1, 10**12), d_rel=0, klass="node_part", at_eV=7.0,
             terms={"|C_1|*eps_R": dec(width_k*value_k*e)[0], "eps_1*|R|": dec(width_k*e*value_k)[0],
                    "mixed": dec(width_k*e*e)[0]},
             rule="|product - value| <= 1e-12 (absolute), unmultiplied panel units")
    for name, missing in (("point path: drop |C_1|*eps_R", width_k*value_k*e),
                          ("point path: drop eps_1*|R|", width_k*e*value_k), ("point path: drop the mixed term", width_k*e*e)):
        book.wrong_answer("X13", "A3f", "control_point_path.order[2].node_parts[inherited][total]@7eV",
                          point - missing, point, Fr(1, 10**12), 0, name)
    nominal = m*m*(n_r*e_area + cfg.norms["total"]*e_area + e_area*e_area)
    for name, missing in (("area path: drop ||R||*e_1", m*m*n_r*e_area),
                          ("area path: drop ||C_1||*e_R", m*m*cfg.norms["total"]*e_area),
                          ("area path: drop the mixed term", m*m*e_area*e_area)):
        book.wrong_answer("X13", "A3f", "order[2].area_error_parts[inherited][total]", nominal - missing, nominal,
                          0, 0, name, rule="below the rebuilt value")

    # AE1: K_LIN, error 1e-3 on 13a, no_hole, total (interpolation only); hole = 0.02*(T - 3)
    row = rows["AE1"]
    spec = table["kernels"]["K_LIN"]
    a, b = hexfr(spec["window_eV"][0]), hexfr(spec["window_eV"][1])
    lin, = spec["columns"]["term_14"]
    box, = spec["columns"]["term_13a"]
    slope, zero = hexfr(lin["slope_per_eV"]), hexfr(lin["zero_eV"])
    hole = PP([(a, b, [-slope*zero, slope])])
    n_hole = hole.l1_linear()
    n_r = hexfr(box["value"])*(b - a)
    e_r = hexfr(row["inputs"]["inject"]["errors"]["term_13a"])*(b - a)
    m = Fr(float(table["defaults"]["analytic"]["weight_length_A"]))
    E = m*m*n_hole*e_r
    need(abs(E - Fr(6, 10**4)) < Fr(1, 10**18), "AE1: E")
    book.put("AE1", "order[2].area_error_parts[inherited][hole]", "error_part", E, klass="young_area", column="hole",
             norm_R_per_A=dec(n_r)[0], norm_column_per_A=dec(n_hole)[0], nominal_area_error_R_per_A=dec(e_r)[0],
             nominal_multiplier_A=dec(m)[0],
             rule="rebuilt = m*nR*e_hole + m*n_hole*e_R + e_hole*e_R with the reported seeds (dimensionless); the "
                  "product value is not below max(value, rebuilt) and not above rebuilt*(1 + 1e-9) + 8*eps*rebuilt; "
                  "value = m^2*||hole||_1*(b - a)*1e-3 is the exact minimum (|hole| split at T = 3)")
    # X33: 17-point Gauss rule on |hole| over the unsplit panel [2, 4]
    with mp.workdps(DPS_HIGH):
        nodes = _gauss_nodes(17)
        unsplit = sum(w*abs(x) for x, w in nodes)          # int_{-1}^{1} |u| du by the rule
        wrong_norm = tomp(slope)*unsplit + tomp(slope)*(tomp(b) - 4)*(tomp(b) + 4 - 2*tomp(zero))/2
        wrong = mpfr(tomp(m*m*e_r)*wrong_norm)
        book.note("AE1", "gauss17_of_abs_u_on_[-1,1]", mp.nstr(unsplit, 20))
    book.wrong.append({"mutation": "X33", "row": "AE1", "key": "order[2].area_error_parts[inherited][hole]",
                       "wrong_value": dec(wrong)[0], "reference": dec(E)[0], "difference": dec(wrong - E)[0],
                       "ratio": dec(wrong/E)[0], "width": "not below the reference", "outside_width": bool(wrong < E),
                       "rule": "below the reference", "note": "17-point Gauss rule on |hole| over the unsplit panel [2,4]"})


def _gauss_nodes(n):
    """Gauss-Legendre nodes and weights on [-1, 1] at the current precision (Newton on P_n)."""
    out = []
    for k in range(1, n + 1):
        x = mp.cos(mp.pi*(k - mp.mpf(1)/4)/(n + mp.mpf(1)/2))
        for _ in range(100):
            p0, p1 = mp.mpf(1), x
            for j in range(1, n):
                p0, p1 = p1, ((2*j + 1)*x*p1 - j*p0)/(j + 1)
            dp = n*(x*p1 - p0)/(x*x - 1)
            step = p1/dp
            x -= step
            if abs(step) < mp.mpf(10)**(-mp.dps + 5):
                break
        p0, p1 = mp.mpf(1), x
        for j in range(1, n):
            p0, p1 = p1, ((2*j + 1)*x*p1 - j*p0)/(j + 1)
        dp = n*(x*p1 - p0)/(x*x - 1)
        out.append((x, 2/((1 - x*x)*dp*dp)))
    return out


def known_wrong(book, table, rows, keep):
    """Wrong values of the mutations of prereg 9 that are computable by the reference; all outside the width."""
    def ref(rid, key):
        rec = book.rows[rid]["keys"][key]
        return Fr(rec["value"]), Fr(rec["width"]["d_abs"]), Fr(rec["width"]["d_rel"])

    def add(mutation, rid, key, wrong, text, u=Fr(0)):
        r, d_abs, d_rel = ref(rid, key)
        return book.wrong_answer(mutation, rid, key, wrong, r, d_abs, d_rel, text, u=u + Fr(1, 10**30)*abs(r))

    # ---- A3a family
    cfg = keep["A3a"][0]
    m, K = cfg.m, cfg.K
    R = K["term_13a"]
    add("X01", "A3a", "order[2].bin[9].total", cfg.tp_masses(K["total"].conv(K["total"]).scale(m*m))["bin[9]"],
        "later losses use the total column: C_2 = K_total * K_total")
    cb = keep["A3b"][0]
    add("X02", "A3b", "order[2].bin[9].total",
        cb.tp_masses(cb.K["total"].conv(cb.K["no_hole"]).scale(cb.m**2))["bin[9]"],
        "later losses use the no_hole column (hole charge 0): C_2 = K_total * K_no_hole")
    g = K["total"] - R
    x04 = (R + g.scale(2)).conv(R).scale(m*m)
    add("X04", "A3a", "order[2].bin[9].total", cfg.tp_masses(x04)["bin[9]"], "the excess is multiplied by n = 2")
    add("X04", "A3a", "order[2].area.total", x04.antider().total(), "the excess is multiplied by n = 2")
    masses = keep["A3a"][2]
    add("X20", "A3a", "order[1].bin[12].hole", max(Fr(0), masses[1]["hole"]["bin[12]"]), "negative values clipped to 0")
    add("X20", "A3a", "order[2].bin[10].hole", max(Fr(0), masses[2]["hole"]["bin[10]"]), "negative values clipped to 0")
    need(masses[1]["hole"]["bin[12]"] < 0 and masses[2]["hole"]["bin[10]"] < 0, "X20: the reference is not negative")
    mc = keep["A3c"][2]
    need(mc[1]["total"]["bin[12]"] < 0, "X20: A3c reference is not negative")
    add("X20", "A3c", "order[1].bin[12].total", Fr(0), "negative values clipped to 0")
    cd = keep["A3d"][0]
    c2 = cd.K["total"].conv(cd.K["term_13a"])
    add("X06", "A3d", "order[2].bin[9].total", cd.tp_masses(c2.scale(cd.m*cd.length))["bin[9]"],
        "the repeated factor uses L_w = 2 instead of m = 3")
    add("X07", "A3d", "order[2].area.total", c2.scale(cd.m).antider().total(), "each order carries m only once")
    add("X07", "A3d", "order[2].bin[10].total", cd.tp_masses(c2.scale(cd.m))["bin[10]"], "each order carries m only once")
    ce = keep["A3e"][0]
    add("X08", "A3e", "order[1].area.no_hole", Fr(1), "S2 per component: m_no_hole = 1/A_no_hole")
    add("X08", "A3e", "order[1].area.hole", Fr(1), "S2 per component: m_hole = 1/A_hole")
    r_true, r_wrong = ce.m*ce.areas["term_13a"], ce.areas["term_13a"]/ce.areas["no_hole"]
    book.wrong.append({"mutation": "X08", "row": "A3e", "key": "order_area_ratio", "wrong_value": dec(r_wrong)[0],
                       "reference": dec(r_true)[0], "difference": dec(r_wrong - r_true)[0],
                       "width": "interval_contains", "outside_width": abs(r_wrong - r_true) > Fr(1, 10**3),
                       "rule": "an interval of the registered accuracy around the wrong value cannot contain the "
                               "reference (difference > 1e-3)", "note": "S2 per component: r = A_13a/A_no_hole"})

    # ---- A1 family
    cfg, orders, masses, areas = keep["A1"]
    add("X05", "A1", "order[2].bin[5].total", masses[2]["total"]["bin[5]"]/2, "each order divided by n!")
    add("X05", "A1", "order[2].area.total", areas[2]["total"]/2, "each order divided by n!")
    add("X15", "A1", "order[1].outside.high.total", Fr(0), "outside.high dropped")
    add("X15", "A1", "order[6].outside.low.total", masses[6]["total"]["outside.low"] + masses[6]["total"]["outside.high"],
        "outside.low taken as the balance residual S*w_6 - sum of the bins")
    w, r = cfg.m*cfg.areas["total"], cfg.m*cfg.areas["term_13a"]
    add("X16", "A1", "order[6].outside.low.total",
        masses[6]["total"]["outside.low"] + cfg.s_abs*w*r**6/(1 - r), "the neglected total area added to outside.low")
    cfg, orders, masses, areas = keep["T2c"]
    w, r = cfg.m*cfg.areas["total"], cfg.m*cfg.areas["term_13a"]
    add("X16", "T2c", "order[4].outside.low.total",
        masses[4]["total"]["outside.low"] + cfg.s_abs*w*r**4/(1 - r), "the neglected total area added to outside.low")
    cfg, orders, masses, areas = keep["A1s"]
    need(masses[2]["total"]["bin[2]"] > 0, "X09 reference")
    add("X09", "A1s", "order[2].bin[2].total", Fr(0), "second order set to zero for s > b (bin[2] is loss 5..6 eV)")

    # ---- X19: the sum omits the last order (INV1 bin[9])
    cfg, solution, orders = keep["INV1"]
    last = cfg.tp_masses(orders[cfg.N]["total"])
    rec = book.rows["INV1"]["keys"]["spectrum.bin[9]"]
    book.wrong_answer("X19", "INV1", "spectrum.bin[9]", solution["bin[9]"] - last["bin[9]"], solution["bin[9]"],
                      Fr(rec["width"]["d_abs"]), Fr(rec["width"]["d_rel"]), "the sum omits order N = 2",
                      u=Fr(1, 10**30))

    # ---- A2: X10, X14, X23
    cfg, t, val = keep["A2"]
    kap, mu_, a, b, m = t.kappa, t.mu, t.a, t.b, t.m
    lo, hi = None, None
    for name, e0, e1 in cfg.regions:
        if name == "bin[5]":
            lo, hi = cfg.px[0] - e1, cfg.px[0] - e0

    def correlation_mass():
        # C_wrong(s) = int K(s + t) K(t) dt = kappa^2 e^(-mu s)(e^(-2 mu a) - e^(-2 mu (b - s)))/(2 mu), 0 <= s <= b - a
        mu, aa, bb, k2 = tomp(mu_), tomp(a), tomp(b), tomp(kap)**2
        x0, x1 = max(tomp(lo), mp.mpf(0)), min(tomp(hi), bb - aa)
        prim = lambda s: k2/(2*mu)*(-mp.exp(-mu*s - 2*mu*aa)/mu - mp.exp(-2*mu*bb + mu*s)/mu)
        return tomp(m)**2*(prim(x1) - prim(x0))
    wrong, u = two_precisions(correlation_mass)
    add("X10", "A2", "order[2].bin[5].total", wrong, "C_1(s - t) replaced by C_1(s + t) (correlation)", u=u)
    wrong, u = two_precisions(lambda: t.laplace(Fr(float("0.3")))*t.laplace(-Fr(float("0.3"))))
    add("X10", "A2", "panels.laplace[2][0.3]", wrong, "correlation: L(tau)*L(-tau)", u=u)
    add("X10", "A2", "panels.moment[2][1]", Fr(0), "correlation: the first moment vanishes")
    h1, _ = two_precisions(lambda: tomp(m*a*kap)*mp.exp(-tomp(mu_)*tomp(a)))
    w1, _ = two_precisions(lambda: tomp(m)*t.area())
    add("X23", "A2", "order[3].cutoff_sensitivity.total", h1*w1**2 + w1*h1*w1, "(n - 1) dropped from h_n",
        u=Fr(1, 10**30))
    # X14: first order averaged over the display cells on the loss axis, then convolved
    cells = sorted((cfg.px[0] - e1, cfg.px[0] - e0) for _, e0, e1 in cfg.regions if e0 is not None and e1 is not None)
    coarse = TP()
    for c0, c1 in cells:
        x0, x1 = max(c0, a), min(c1, b)
        if x0 < x1:
            mean, _ = two_precisions(lambda: tomp(kap)*(mp.exp(-tomp(mu_)*tomp(x0)) - mp.exp(-tomp(mu_)*tomp(x1)))
                                     / tomp(mu_)/tomp(c1 - c0))
            coarse = coarse + TP.box(mean, c0, c1)
    both = cfg.tp_masses(coarse.conv(coarse).scale(m*m))
    outside_bins = []
    for name in both:
        key = f"order[2].{name}.total"
        r, d_abs, d_rel = ref("A2", key)
        if abs(both[name] - r) > d_abs + d_rel*abs(r) + Fr(1, 10**25):
            outside_bins.append(name)
    need("bin[5]" in outside_bins, "X14: bin[5] is not outside the width")
    add("X14", "A2", "order[2].bin[5].total", both["bin[5]"],
        "first order replaced by its cell means (both factors) before the second convolution; bins outside the "
        f"width: {outside_bins}", u=Fr(1, 10**25))

    def mixed_mass():
        # coarse first order convolved with the exact repeated kernel kappa*exp(-mu*t) on [a, b]
        mu, k = tomp(mu_), tomp(kap)

        def g2(y):          # antiderivative of the cumulative of the exact kernel
            y = tomp(y)
            aa, bb = tomp(a), tomp(b)
            if y <= aa:
                return mp.mpf(0)
            z = min(y, bb)
            inner = k/mu*(mp.exp(-mu*aa)*(z - aa) + (mp.exp(-mu*z) - mp.exp(-mu*aa))/mu)
            return inner + (k/mu*(mp.exp(-mu*aa) - mp.exp(-mu*bb))*(y - bb) if y > bb else 0)
        total = mp.mpf(0)
        for c0, c1 in zip(coarse.breaks(), coarse.breaks()[1:]):
            height = tomp(coarse(c0))
            total += height*(g2(hi - c0) - g2(hi - c1) - g2(lo - c0) + g2(lo - c1))
        return tomp(m)**2*total
    wrong, u = two_precisions(mixed_mass)
    add("X14", "A2", "order[2].bin[5].total", wrong,
        "first order replaced by its cell means, repeated kernel exact (second reading of the mutation)", u=u)

    # ---- X12: 9 Gauss points in the inner rule (A4a, A4c)
    for rid, name in (("A4a", "K_A4_9"), ("A4c", "K_A4_16")):
        value, scale, d = keep[rid]
        leg = legendre(d)
        with mp.workdps(DPS_HIGH):
            total = mp.mpf(0)
            for x, w in _gauss_nodes(9):
                pd = sum(tomp(c)*x**k for k, c in enumerate(leg))
                pm = sum(tomp(c)*(-x)**k for k, c in enumerate(leg))
                total += w*(1 + pd)*(1 + pm)
            wrong = mpfr(tomp(scale)**2*mp.mpf(3)/2*total)
        key = f"order_panels[2].value_at_7eV.term_13a"
        add("X12", rid, key, wrong, "9-point inner Gauss rule" + (" (all nodes are zeros of P9: scale^2*3)" if d == 9 else ""),
            u=Fr(1, 10**40))


# ---------------------------------------------------------------- output
def document():
    table, ledger, book = build()
    rows = {r["id"]: r for r in table["rows"]}
    counts, total, missing = {}, 0, []
    for rid in table["execution_order"]:
        planned = book.planned(rid)
        have = {(k, v["kind"]) for k, v in book.rows.get(rid, {"keys": {}})["keys"].items()}
        if set(planned) != have:
            missing.append((rid, sorted(set(planned) - have)[:5], sorted(have - set(planned))[:5]))
        if planned:
            by_kind = {}
            for _, kind in planned:
                by_kind[kind] = by_kind.get(kind, 0) + 1
            counts[rid] = dict(sorted(by_kind.items()), total=len(planned))
            total += len(planned)
        # the expansion agrees with the counts of the frozen table
        for p in rows[rid]["check_plan"]:
            if p["kind"] in EXTERNAL and (p["oracle"].startswith("R-A") or p["oracle"].startswith("R-INV")):
                n = sum(1 for k, kind, oracle in ledger[rid] if kind == p["kind"] and oracle == p["oracle"])
                same = sum(q["count"] for q in rows[rid]["check_plan"] if q["kind"] == p["kind"] and q["oracle"] == p["oracle"])
                need(n == same, f"{rid}: the expansion of {p['keys']} differs from the count of the table")
    need(not missing, f"reference keys differ from the ledger: {missing[:3]}")
    planned_rb = table["counts"]["planned_numeric_keys_by_group"]["RM-B"]
    rb_inside = sum(1 for rid in ledger for k, kind, oracle in ledger[rid] if kind in EXTERNAL and oracle.startswith("R-B"))
    need(total + rb_inside == table["counts"]["planned_external_reference_keys"], "external key count")
    for item in book.wrong:
        need(item["outside_width"], f"a known wrong answer is inside the width: {item}")
    need(all(book.self_checks.values()), "self-check")
    by_mutation = {}
    for item in book.wrong:
        by_mutation[item["mutation"]] = by_mutation.get(item["mutation"], 0) + 1
    shortage = []
    for rid, row in book.rows.items():
        keys = row["keys"]
        bad = sum(v["status"] != "reference" for v in keys.values())
        if keys and (10*bad > len(keys) or (bad and not any(v["kind"] == "direct_reference" for v in keys.values()))):
            shortage.append(rid)
    doc = {
        "schema": SCHEMA,
        "registration_table_sha256": TABLE_SHA256,
        "registration_table_script_sha256": TABLE_PY_SHA256,
        "generator": "benchmarks/numerics/m4b1/generate_references_ra_m4b1.py",
        "imports_package_under_test": False,
        "arithmetic": {"exact": "fractions.Fraction on the binary64 inputs read as rationals",
                       "multiprecision": f"mpmath {mpmath.__version__}, {DPS_LOW} and {DPS_HIGH} digits",
                       "decimal_digits": DIGITS,
                       "u_r": "decimal rounding (half a unit of the last digit) plus, for multiprecision values, "
                              "4 times the difference between the two working precisions (and quadrature settings) "
                              "and a floor of 1e-43 relative"},
        "record_format": {"value": "decimal string", "u_r": "decimal string, rounded up",
                          "width": "d_abs and d_rel of prereg 7.2 (rounded down); null where the rule is not "
                                   "|y - r| + u_r <= d_abs + d_rel*|r| (see the field rule)",
                          "status": "reference or no_reference (prereg 6.3, applied before any implementation value)"},
        "no_reference_rule": {"zero_reference": "u_r = 0 exactly (prereg 6.3 item 1)",
                              "nonzero_reference": "u_r <= width/4, otherwise no_reference (prereg 6.3 item 2)"},
        "counts": {"keys": total, "by_row": counts, "rows": len(counts),
                   "by_kind": {k: sum(c.get(k, 0) for c in counts.values()) for k in EXTERNAL},
                   "table_external_reference_keys": table["counts"]["planned_external_reference_keys"],
                   "table_keys_of_R-B": rb_inside},
        "max_u_r_over_width": {rid: {"ratio": dec_up(ratio), "key": key}
                               for rid, (ratio, key) in sorted(book.max_ratio.items())},
        "no_reference": book.no_reference,
        "rows_short_of_references": shortage,
        "self_checks": dict(sorted(book.self_checks.items())),
        "known_wrong_answers": book.wrong,
        "known_wrong_answers_by_mutation": dict(sorted(by_mutation.items())),
        "rows": {rid: book.rows[rid] for rid in table["execution_order"] if rid in book.rows},
    }
    text = json.dumps(doc, ensure_ascii=True, indent=1, sort_keys=True) + "\n"
    return text.encode("ascii"), doc


def main(argv):
    start = time.perf_counter()
    data, doc = document()
    digest = hashlib.sha256(data).hexdigest()
    elapsed = time.perf_counter() - start
    summary = {"sha256": digest, "keys": doc["counts"]["keys"], "rows": doc["counts"]["rows"],
               "no_reference": len(doc["no_reference"]), "known_wrong_answers": len(doc["known_wrong_answers"]),
               "seconds": round(elapsed, 2), "python": platform.python_version(), "mpmath": mpmath.__version__}
    if argv == ["--check"]:
        same = OUT.is_file() and OUT.read_bytes() == data
        print(json.dumps({"identical": same, **summary}))
        return 0 if same else 1
    need(not argv, "usage: generate_references_ra_m4b1.py [--check]")
    OUT.write_bytes(data)
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
