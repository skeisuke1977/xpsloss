"""R-A: analytic references of M4a (prereg 6.1). Does not import xpsloss.

The binary64 inputs of the registration table are read as exact rationals.  The
closed forms are evaluated with rationals or with mpmath at two working
precisions; u_r is the difference of the two plus the decimal rounding.

    python benchmarks/numerics/m4a/generate_references_ra_m4a.py          # write
    python benchmarks/numerics/m4a/generate_references_ra_m4a.py --check  # compare bytes
"""
from decimal import Decimal, getcontext
from fractions import Fraction as Fr
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import mpmath as mp

HERE = Path(__file__).resolve().parent
TABLE = HERE/"expected_states_m4a.json"
OUT = HERE/"references_ra_m4a.json"
TABLE_SHA256 = "250b45f863ccbd1c4c94474f23769976380c5dad5c2671328712ef4acdc7f408"
DPS = (50, 90)
DIGITS = 40
COLUMNS = ("term_13a", "term_13b", "term_14", "term_15_T1", "term_15_T2", "term_15_T3",
           "term_16_e", "term_16_h", "no_hole", "hole", "total")
NO_HOLE = ("term_13a", "term_13b", "term_15_T1", "term_15_T3", "term_16_e")
HOLE = ("term_14", "term_15_T2", "term_16_h")
QUANTITIES = ("zero_loss", "single_loss_total", "single_loss_no_hole", "single_loss_hole", "spectrum")
EPS = Fr(1, 2**52)
TINY = Fr(1, 2**400)  # smaller than any gap between two different break points of binary64 inputs


def fd(x):
    """A decimal literal of the table, as the binary64 that is passed to the code."""
    return Fr(*float(x).as_integer_ratio())


def fx(obj):
    if isinstance(obj, dict):
        value = float.fromhex(obj["hex"])
        if value != float(obj["value"]):
            raise AssertionError(obj)
        return Fr(*value.as_integer_ratio())
    return fd(obj)


def M(x):
    return mp.mpf(x.numerator)/mp.mpf(x.denominator) if isinstance(x, Fr) else x


def add(x, y):
    return x + y if isinstance(x, Fr) and isinstance(y, Fr) else M(x) + M(y)


def mul(x, y):
    return x*y if isinstance(x, Fr) and isinstance(y, Fr) else M(x)*M(y)


# ---------------------------------------------------------------- kernel forms

class Poly:
    """sum_k c_k T^k with rational coefficients (exact)."""
    def __init__(self, coefficients):
        self.c = [Fr(c) for c in coefficients]

    def val(self, t):
        return sum((c*t**k for k, c in enumerate(self.c)), Fr(0))

    def integ_lin(self, t1, t2, al, be):
        out = Fr(0)
        for k, c in enumerate(self.c):
            out += c*(al*(t2**(k+1) - t1**(k+1))/(k+1) + be*(t2**(k+2) - t1**(k+2))/(k+2))
        return out


class InvT:
    """c/T."""
    def __init__(self, c):
        self.c = c

    def val(self, t):
        return self.c/t

    def integ_lin(self, t1, t2, al, be):
        return add(mul(self.c*al, mp.log(M(t2)/M(t1))), self.c*be*(t2 - t1))


class Lorentz:
    """c/(1 + ((T - t0)/g)^2)."""
    def __init__(self, c, t0, g):
        self.c, self.t0, self.g = c, t0, g

    def val(self, t):
        return self.c/(1 + ((t - self.t0)/self.g)**2)

    def integ_lin(self, t1, t2, al, be):
        x1, x2 = M((t1 - self.t0)/self.g), M((t2 - self.t0)/self.g)
        out = mul(self.c*self.g*(al + be*self.t0), mp.atan(x2) - mp.atan(x1))
        if be:
            out = add(out, mul(self.c*self.g*self.g*be/2, mp.log1p(x2*x2) - mp.log1p(x1*x1)))
        return out


def poly_mul(p, q):
    out = [Fr(0)]*(len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i+j] += x*y
    return out


def poly_add(p, q, scale=1):
    out = [Fr(0)]*max(len(p), len(q))
    for i, x in enumerate(p):
        out[i] += x
    for i, y in enumerate(q):
        out[i] += scale*y
    return out


def chebyshev_in_T(n, x_of_t):
    """T_n(x(T)) as a polynomial in T; x_of_t are the coefficients of x(T)."""
    previous, current = [Fr(1)], list(x_of_t)
    for _ in range(n - 1):
        previous, current = current, poly_add(poly_mul([2*c for c in x_of_t], current), previous, -1)
    return current


def kernel_pieces(row):
    """[(low, high, {term: form})]; bounds None are infinite; zero outside the pieces."""
    k = row["inputs"]["kernel"]
    form = k["form"]
    if form == "rectangle":
        lo, hi = map(fd, k["support_eV"])
        terms = k["terms"] if "terms" in k else {"term_13a": k["term_13a"]}
        return [(lo, hi, {name: Poly([fd(v)]) for name, v in terms.items()})]
    if form == "piecewise_linear":
        if (k["term_13a"], k["term_14"], k["support_eV"]) != ("0.1*max(0, min(T-2, (5-T)/2))", "0.02*(T-3)", [2, 5]):
            raise AssertionError(k)
        c1, c2 = fd(0.1), fd(0.02)
        hole = Poly([-3*c2, c2])
        return [(Fr(2), Fr(3), {"term_13a": Poly([-2*c1, c1]), "term_14": hole}),
                (Fr(3), Fr(5), {"term_13a": Poly([c1*Fr(5, 2), -c1/2]), "term_14": hole})]
    if form == "c_over_T":
        if k["term_13a"] != "0.03/T":
            raise AssertionError(k)
        return [(None, None, {"term_13a": InvT(fd(k["c"]))})]
    if form == "lorentzian":
        if "c" in k:
            if k["expression"] != "c/(1+((T-t0)/gamma)^2)":
                raise AssertionError(k)
            form = Lorentz(fd(k["c"]), fx(k["t0_eV"]), fx(k["gamma_eV"]))
        elif k["expression"] == "1/(1+((T-2.25)/0.3)^2)":
            form = Lorentz(Fr(1), fd(2.25), fd(0.3))
        else:
            raise AssertionError(k)
        return [(None, None, {k["column"]: form})]
    if form == "chebyshev":
        if k["expression"] != "kappa*(1 - T_32(2*T-5))":
            raise AssertionError(k)
        kappa = fx(k["kappa"])
        lo, hi = map(fd, k["support_eV"])
        t32 = chebyshev_in_T(32, [Fr(-5), Fr(2)])
        return [(lo, hi, {k["column"]: Poly([kappa*c for c in poly_add([Fr(1)], t32, -1)])})]
    raise AssertionError(form)


def column_terms(column, charge):
    """The column as a list of (coefficient, term)."""
    if column in NO_HOLE or column in HOLE:
        return [(Fr(1), column)]
    if column == "no_hole":
        return [(Fr(1), t) for t in NO_HOLE]
    if column == "hole":
        return [(Fr(1), t) for t in HOLE]
    return [(Fr(1), t) for t in NO_HOLE] + [(charge, t) for t in HOLE]


def inside(piece, t):
    return (piece[0] is None or piece[0] < t) and (piece[1] is None or t < piece[1])


def kernel_value(pieces, column, charge, t, side):
    """K_column(t + side*0), side = +1 or -1 (exact rational)."""
    probe = t + side*TINY
    out = Fr(0)
    for piece in pieces:
        if inside(piece, probe):
            for coefficient, term in column_terms(column, charge):
                if term in piece[2]:
                    out += coefficient*piece[2][term].val(t)
    return out


def integrate(pieces, column, charge, a, b, weight_cuts=(), weight=None):
    """int_a^b K_column(T) g(T) dT, g piecewise linear with breaks in weight_cuts (g = 1 if weight is None)."""
    cuts = {a, b}
    for lo, hi, _ in pieces:
        cuts.update(x for x in (lo, hi) if x is not None and a < x < b)
    cuts.update(x for x in weight_cuts if a < x < b)
    cuts = sorted(cuts)
    total = Fr(0)
    for t1, t2 in zip(cuts, cuts[1:]):
        if weight is None:
            al, be = Fr(1), Fr(0)
        else:
            m1, m2 = t1 + (t2 - t1)/3, t1 + 2*(t2 - t1)/3
            g1, g2 = weight(m1), weight(m2)
            be = (g2 - g1)/(m2 - m1)
            al = g1 - be*m1
            if weight((t1 + t2)/2) != al + be*(t1 + t2)/2:
                raise AssertionError("the weight is not linear on a cell")
            if al == 0 and be == 0:
                continue
        middle = (t1 + t2)/2
        for piece in pieces:
            if inside(piece, middle):
                for coefficient, term in column_terms(column, charge):
                    if term in piece[2]:
                        total = add(total, mul(coefficient, piece[2][term].integ_lin(t1, t2, al, be)))
    return total


# ---------------------------------------------------------------- primary and output regions

def region_list(edges):
    out = [("outside.low", None, edges[0], False)]
    out += [(f"bin[{i}]", u, v, i == len(edges) - 2) for i, (u, v) in enumerate(zip(edges, edges[1:]))]
    return out + [("outside.high", edges[-1], None, False)]


def member(energy, region):
    _, u, v, last = region
    if u is None:
        return energy < v
    if v is None:
        return energy > u
    return u <= energy < v or (last and energy == v)


def g_value(primary, region, t):
    """G_D(T) = int 1_D(E' - T) F(dE') (exact)."""
    _, u, v, _ = region
    if primary["kind"] == "delta":
        return sum((s for e, s in zip(primary["x"], primary["s"]) if member(e - t, region)), Fr(0))
    out = Fr(0)
    for p, q, f in zip(primary["x"], primary["x"][1:], primary["s"]):
        lo = p if u is None else max(p, u + t)
        hi = q if v is None else min(q, v + t)
        if hi > lo:
            out += f*(hi - lo)/(q - p)
    return out


def g_limit(primary, region, t, side):
    return g_value(primary, region, t + side*TINY) if primary["kind"] == "delta" else g_value(primary, region, t)


def read_primary(table, name):
    p = table["primaries"][name]
    key = "energies_eV" if p["kind"] == "delta" else "edges_eV"
    return {"kind": p["kind"], "x": [Fr(*float.fromhex(h).as_integer_ratio()) for h in p[key]["hex"]],
            "s": [Fr(*float.fromhex(h).as_integer_ratio()) for h in p["intensities"]["hex"]]}


def read_edges(table, name):
    return [Fr(*float.fromhex(h).as_integer_ratio()) for h in table["output_edges"][name]["hex"]]


# ---------------------------------------------------------------- one row

def settings_of(table, row):
    return {**table["defaults"]["settings"], **row["settings"]}


def row_conditions(table, row):
    i = row["inputs"]
    a, b = (fx(x) for x in i["loss_window_eV"])
    return dict(a=a, b=b, L=fx(i["weight_length_A"]), charge=fx(i.get("hole_charge", 1.0)),
                settings=settings_of(table, row))


def compute_layer(table, row):
    """{key: exact rational or mpf} for a row of the analytic layers A (and B)."""
    c = row_conditions(table, row)
    a, b, L, charge = c["a"], c["b"], c["L"], c["charge"]
    pieces = kernel_pieces(row)
    out = {}
    for column in COLUMNS:
        out[f"window_area.{column}"] = integrate(pieces, column, charge, a, b)
        out[f"cutoff_sensitivity.{column}"] = L*a*kernel_value(pieces, column, charge, a, +1)
    out["single_loss_weight"] = mul(L, out["window_area.total"])
    if "primary" not in row["inputs"]:
        return out
    primary = read_primary(table, row["inputs"]["primary"])
    edges = read_edges(table, row["inputs"]["output_edges_eV"])
    breaks = sorted({e - d for e in primary["x"] for d in edges})
    for region in region_list(edges):
        name = region[0]
        weight = lambda t, region=region: g_value(primary, region, t)
        out[f"{name}.zero_loss"] = g_value(primary, region, Fr(0))
        for column in ("total", "no_hole", "hole"):
            out[f"{name}.single_loss_{column}"] = mul(L, integrate(pieces, column, charge, a, b, breaks, weight))
            for label, side in (("plus", +1), ("minus", -1)):
                out[f"{name}.cutoff_sensitivity.{label}.{column}"] = (
                    L*a*kernel_value(pieces, column, charge, a, side)*g_limit(primary, region, a, side))
        out[f"{name}.spectrum"] = add(out[f"{name}.zero_loss"], out[f"{name}.single_loss_total"])
    return out


def balance(table, row, values):
    """Sums over the regions against S and S*L*A (a self check of the reference; not a criterion)."""
    if "primary" not in row["inputs"]:
        return None
    c = row_conditions(table, row)
    primary = read_primary(table, row["inputs"]["primary"])
    names = [r[0] for r in region_list(read_edges(table, row["inputs"]["output_edges_eV"]))]
    s = sum(primary["s"], Fr(0))
    out = {"zero_loss": M(sum((values[f"{n}.zero_loss"] for n in names), Fr(0))) - M(s)}
    for column in ("total", "no_hole", "hole"):
        total = mp.mpf(0)
        for n in names:
            total += M(values[f"{n}.single_loss_{column}"])
        out[f"single_loss_{column}"] = total - M(s)*M(c["L"])*M(values[f"window_area.{column}"])
    return {k: mp.nstr(abs(v), 3) for k, v in out.items()}


def lobatto_nodes(left, right):
    """17 Chebyshev-Lobatto nodes, index 0 at the left end (mpf)."""
    mid, half = (M(left) + M(right))/2, (M(right) - M(left))/2
    return [mid - half*mp.cos(mp.pi*k/16) for k in range(17)]


def basis_integral(nodes, index, lo, hi):
    """int_lo^hi l_index(T) dT for the Lagrange basis on the nodes (expanded around lo)."""
    coefficients = [mp.mpf(1)]
    for k, node in enumerate(nodes):
        if k != index:
            scale = nodes[index] - node
            shift = M(lo) - node
            new = [mp.mpf(0)]*(len(coefficients) + 1)
            for i, c in enumerate(coefficients):
                new[i] += c*shift/scale
                new[i+1] += c/scale
            coefficients = new
    h = M(hi) - M(lo)
    return mp.fsum(c*h**(i+1)/(i+1) for i, c in enumerate(coefficients))


def compute_ra9(table, row):
    c = row_conditions(table, row)
    primary = read_primary(table, row["inputs"]["primary"])
    edges = read_edges(table, row["inputs"]["output_edges_eV"])
    left, right = map(fd, row["inputs"]["declared_features"]["panels_eV"])
    e = fx(row["inputs"]["injection"]["point_error"])
    lo, hi = max(c["a"], primary["x"][0] - edges[1]), min(c["b"], primary["x"][0] - edges[0])
    if (lo, left, right) != (Fr(2), Fr(2), Fr(5)) or primary["s"] != [Fr(1)]:
        raise AssertionError("RA9 cell")
    omega = basis_integral(lobatto_nodes(left, right), 2, lo, hi)
    # sensitivity to the binary64 rounding of the nodes (the implementation rounds its nodes)
    rounded = [mp.mpf(3.5 - 1.5*math.cos(math.pi*k/16)) for k in range(17)]
    rounded[0], rounded[16] = mp.mpf(2), mp.mpf(5)
    node_rounding = abs(basis_integral(rounded, 2, lo, hi) - omega)
    value = abs(omega)*M(c["L"])*M(e)
    keys = {f"bin[0].single_loss_{name}.kernel_point_error": value for name in ("total", "no_hole", "hole")}
    return keys, dict(omega_2=omega, cell=(lo, hi), node_rounding=node_rounding*M(c["L"])*M(e))


def compute_ra10(table, row):
    c = row_conditions(table, row)
    k = row["inputs"]["kernel"]
    primary = read_primary(table, row["inputs"]["primary"])
    edges = read_edges(table, row["inputs"]["output_edges_eV"])
    lo, hi = map(fd, k["support_eV"])
    region = region_list(edges)[1]
    breaks = sorted({e - d for e in primary["x"] for d in edges})
    out = {}
    for name in ("total", "no_hole", "hole"):
        pieces = [(lo, hi, {"term_13a": Poly([fd(k[name])])})]
        out[f"bin[0].single_loss_{name}"] = c["L"]*integrate(
            pieces, "term_13a", Fr(1), c["a"], c["b"], breaks, lambda t: g_value(primary, region, t))
    return out


def compute_ra14(table, row):
    """Area estimate of the total column: verification part + point part (interpolation, edge are 0)."""
    c = row_conditions(table, row)
    k = row["inputs"]["kernel"]
    s = row["inputs"]["declared_features"]["singularities"][0]
    form = Lorentz(Fr(*float.fromhex(k["k_hex"]).as_integer_ratio()), fx(s["center_eV"]), fx(s["half_width_eV"]))
    if k["expression"] != "k*1/(1+((T-2.4)/0.5)^2)" or float.fromhex(k["k_hex"]) != 1e-10/0.7753994963:
        raise AssertionError(k)
    a, b = c["a"], c["b"]
    mid, half = (M(a) + M(b))/2, (M(b) - M(a))/2

    def value(t):
        return M(form.c)/(1 + ((t - M(form.t0))/M(form.g))**2)
    nodes = lobatto_nodes(a, b)
    point = mp.fsum(basis_integral(nodes, q, a, b)*value(nodes[q]) for q in range(17))
    checks = [mid - half*mp.cos(mp.pi*(2*m + 1)/32) for m in range(16)]
    verification = (M(b) - M(a))*max(2*value(t) for t in checks)
    return ({"exception.diagnostics.error_budget.outputs[window_area.total].estimate": verification + point},
            dict(verification=verification, kernel_point_error=point,
                 true_area=form.integ_lin(a, b, Fr(1), Fr(0))))


# ---------------------------------------------------------------- widths (prereg 7.2)

def replacement_strength(row_id, region, total):
    if row_id not in ("RA6", "RA6b"):
        return total
    table = {"outside.low": Fr(1), "bin[0]": Fr(2)**60, "outside.high": Fr(2)}
    if row_id == "RA6":
        table["bin[1]"] = Fr(2)
    return table[region]


def width_of(table, row, key):
    """(d_abs, d_rel, unit) of a direct_reference key."""
    c = row_conditions(table, row)
    s = c["settings"]
    rel = fd(s["epsrel"])
    floor = fd(s["area_epsabs"])
    kernel_floor = Fr(1, 10**12)
    if key.startswith("window_area."):
        return floor, rel, "angstrom^-1"
    if key == "single_loss_weight":
        return c["L"]*floor, rel, "1"
    if key.startswith("cutoff_sensitivity."):
        return c["L"]*c["a"]*kernel_floor, rel, "1"
    region, quantity = key.split(".", 1) if key.startswith("bin[") else (key[:key.index(".", 8)], key[key.index(".", 8)+1:])
    primary = read_primary(table, row["inputs"]["primary"])
    strength = replacement_strength(row["id"], region, sum((abs(x) for x in primary["s"]), Fr(0)))
    if quantity in QUANTITIES:
        return fd(s["spectrum_epsabs_fraction"])*strength, rel, "I"
    if quantity.startswith("cutoff_sensitivity."):
        return c["L"]*c["a"]*kernel_floor*strength, rel, "I"
    raise AssertionError(key)


# ---------------------------------------------------------------- decimal output

def exact(x):
    if isinstance(x, Fr):
        return x
    sign, man, exp, _ = x._mpf_
    return Fr(-int(man) if sign else int(man))*Fr(2)**int(exp)


def decimal_string(x):
    getcontext().prec = DIGITS
    if x == 0:
        return "0"
    return format(Decimal(x.numerator)/Decimal(x.denominator), f".{DIGITS-1}e")


def bound_string(x):
    """x >= 0 rounded up to 4 significant digits."""
    if x == 0:
        return "0"
    exponent = math.floor(mp.log10(M(x)))
    mantissa = int(mp.ceil(M(x)/mp.mpf(10)**(exponent - 3))) + 1
    return f"{mantissa}e{exponent-3}"


def at_precisions(function):
    results = []
    for dps in DPS:
        mp.mp.dps = dps
        results.append(function())
    mp.mp.dps = DPS[-1]
    return results


def record(low, high, extra=Fr(0)):
    """Value string, u_r string and the parsed rationals."""
    low, high = exact(low), exact(high)
    text = decimal_string(high)
    stored = Fr(text)
    u = abs(low - high) + abs(stored - high) + exact(extra)
    u_text = bound_string(u)
    if Fr(u_text) < u:
        raise AssertionError("u_r not rounded up")
    return text, u_text, stored, Fr(u_text)


# ---------------------------------------------------------------- RA13 prototype (prereg 6.1, M14a)

def ra13_prototype(table, row, reference, width):
    """An own binary64 prototype: 17-point Clenshaw-Curtis panels on the registered initial mesh.

    'truncated' stops after the six initial panels; 'normal' bisects until the sum of the
    nested estimates |CC17 - CC9| is within epsrel*|A| (the success rule of prereg 2).
    """
    c = row_conditions(table, row)
    a, b = float(c["a"]), float(c["b"])
    half_width = float(fx(row["inputs"]["declared_features"]["min_half_width_eV"]))
    mp.mp.dps = 40
    unit = [mp.cos(mp.pi*(16 - k)/16) for k in range(17)]
    w17 = [float(basis_integral(unit, q, -1, 1)) for q in range(17)]
    w9 = [float(basis_integral(unit[::2], q, -1, 1)) for q in range(9)]
    x17 = [float(u) for u in unit]
    mp.mp.dps = DPS[-1]

    def kernel(t):
        return 1.0/(1.0 + ((t - 2.25)/0.3)**2)

    def panel(lo, hi):
        mid, half = 0.5*(lo + hi), 0.5*(hi - lo)
        values = [kernel(lo)] + [kernel(mid + half*x) for x in x17[1:-1]] + [kernel(hi)]
        i17 = half*math.fsum(w*v for w, v in zip(w17, values))
        i9 = half*math.fsum(w*v for w, v in zip(w9, values[::2]))
        return i17, abs(i17 - i9)

    # initial mesh: grading a*2^k, then equal parts no wider than twice the minimal half width
    grading = [a]
    while grading[-1]*2 < b:
        grading.append(grading[-1]*2)
    grading.append(b)
    edges = [a]
    for lo, hi in zip(grading, grading[1:]):
        n = math.ceil(Fr(hi - lo)/Fr(2*half_width))
        edges += [lo + (hi - lo)*i/n for i in range(1, n + 1)]
    panels = [(lo, hi, *panel(lo, hi)) for lo, hi in zip(edges, edges[1:])]
    initial = len(panels)
    truncated = math.fsum(p[2] for p in panels)
    stages = 0
    rel = c["settings"]["epsrel"]
    while True:
        area = math.fsum(p[2] for p in panels)
        estimate = math.fsum(p[3] for p in panels)
        budget = max(c["settings"]["area_epsabs"], rel*abs(area))
        if estimate <= budget:
            break
        share = budget/(b - a)
        new = []
        for lo, hi, value, est in panels:
            if est > 0.5*share*(hi - lo):
                mid = 0.5*(lo + hi)
                new += [(lo, mid, *panel(lo, mid)), (mid, hi, *panel(mid, hi))]
            else:
                new.append((lo, hi, value, est))
        if len(new) == len(panels) or len(new) > 400:
            raise AssertionError("the RA13 prototype did not converge")
        panels, stages = new, stages + 1
    r, u = reference

    def judge(y):
        d = abs(Fr(y) - r)
        return dict(value=float(y).hex(), abs_difference=bound_string(d), relative_difference=bound_string(d/abs(r)),
                    passes=bool(d + u <= width))
    return dict(key="window_area.term_13a", initial_panels=initial, final_panels=len(panels), stages=stages,
                final_estimate=float(estimate), normal=judge(area), truncated_after_initial_panels=judge(truncated),
                width=bound_string(width))


# ---------------------------------------------------------------- known wrong answers (prereg 6.1, 9)

def known_wrong(table, rows, values):
    """Values that the registered mutations would give; all must lie outside the width."""
    def ref(row_id, key):
        return exact(values[row_id][key][1])
    out = []

    def case(mutation, row_id, key, wrong, note):
        r = ref(row_id, key)
        d_abs, d_rel, _ = width_of(table, rows[row_id], key) if not key.endswith("kernel_point_error") else (Fr(0), Fr(1, 10**9), "I")
        width = d_abs + d_rel*abs(r) if r else d_abs
        wrong = exact(wrong)
        out.append(dict(mutation=mutation, row=row_id, key=key, wrong_value=decimal_string(wrong), note=note,
                        reference=decimal_string(r), width=bound_string(width),
                        outside_width=bool(abs(wrong - r) > width)))
    c1 = row_conditions(table, rows["RA1"])
    area1 = ref("RA1", "window_area.total")
    total1 = sum((ref("RA1", f"{n}.single_loss_total") for n in ("outside.low", "outside.high", *(f"bin[{i}]" for i in range(7)))), Fr(0))
    for b in (2, 3):
        key = f"bin[{b}].single_loss_total"
        case("M01", "RA1", key, Fr(0), "axis reversed: the displayed satellite is 0")
        case("M02a", "RA1", key, ref("RA1", key)/c1["L"], "L dropped (sum 0.30)")
        case("M03", "RA1", key, ref("RA1", key)/c1["L"]/area1, "L replaced by 1/A (sum 1.0)")
    case("M01", "RA1", "outside.high.single_loss_total", total1, "axis reversed: 0.60 on the high side")
    case("M12", "RA1", "bin[5].zero_loss", Fr(1), "delta at the inner edge counted on both sides")
    case("M10", "RA2", "outside.low.single_loss_total", ref("RA2", "outside.high.single_loss_total"), "outsides exchanged")
    case("M10", "RA2", "outside.high.single_loss_total", ref("RA2", "outside.low.single_loss_total"), "outsides exchanged")
    case("M05b", "RA3", "window_area.total", ref("RA3", "window_area.no_hole") + ref("RA3", "window_area.hole"),
         "total rebuilt with the coefficient 1 (0.18)")
    case("M07", "RA3", "bin[2].single_loss_hole", 3*ref("RA3", "bin[2].single_loss_hole"),
         "primary bin width 3 multiplied once more")
    case("M13", "RA3", "bin[2].single_loss_hole", Fr(0), "negative hole clipped")
    # M09: G fixed at the midpoint of each cell; cells [2,3], [3,4] (the coarsest possible) and 400 equal cells
    c3 = row_conditions(table, rows["RA3"])
    primary3 = read_primary(table, "B_RA3")
    region3 = region_list(read_edges(table, "E_RA3"))[3]
    pieces3 = kernel_pieces(rows["RA3"])
    for cells, note in ((2, "G at the cell midpoint, cells [2,3] and [3,4]"),
                        (400, "G at the cell midpoint, 400 equal cells on [2,4] (at the panel limit)")):
        wrong = Fr(0)
        for i in range(cells):
            lo, hi = 2 + Fr(2*i, cells), 2 + Fr(2*i + 2, cells)
            wrong += g_value(primary3, region3, (lo + hi)/2)*integrate(pieces3, "hole", c3["charge"], lo, hi)
        case("M09", "RA3", "bin[2].single_loss_hole", c3["L"]*wrong, note)
    # M08a: half-width factor dropped on the initial graded mesh a*2^k of RA4a (finer meshes are further away)
    c4 = row_conditions(table, rows["RA4a"])
    pieces4 = kernel_pieces(rows["RA4a"])
    edges, wrong = [c4["a"]], mp.mpf(0)
    while edges[-1]*2 < c4["b"]:
        edges.append(edges[-1]*2)
    edges.append(c4["b"])
    for lo, hi in zip(edges, edges[1:]):
        wrong += M(integrate(pieces4, "total", Fr(1), lo, hi))*2/M(hi - lo)
    for key in ("window_area.term_13a", "window_area.total"):
        case("M08a", "RA4a", key, wrong, "half-width factor dropped on the initial graded mesh (7 panels)")
    # M11a, M11b: residuals in binary64
    low6 = float(ref("RA6", "outside.low.single_loss_total"))
    b0, b1 = float(ref("RA6", "bin[0].single_loss_total")), float(ref("RA6", "bin[1].single_loss_total"))
    whole = (2.0**60 + 3.0)*2.0*float(ref("RA6", "window_area.total"))
    case("M11a", "RA6", "outside.low.single_loss_total", Fr(whole - b0 - b1 - 0.0),
         "outside as the binary64 residual S*L*A - bins - other outside")
    case("M11b", "RA6", "bin[1].single_loss_total", Fr((low6 + b0 + b1) - (low6 + b0)),
         "inner bin as a binary64 difference of cumulative sums")
    high6b = float(ref("RA6b", "outside.high.single_loss_total"))
    b0b = float(ref("RA6b", "bin[0].single_loss_total"))
    low6b = float(ref("RA6b", "outside.low.single_loss_total"))
    case("M11a", "RA6b", "outside.high.single_loss_total", Fr(whole - b0b - low6b),
         "outside as the binary64 residual S*L*A - bins - other outside")
    del high6b
    case("M15", "RA7", "bin[1].cutoff_sensitivity.plus.total", ref("RA7", "bin[1].cutoff_sensitivity.minus.total"),
         "one-sided sensitivities exchanged")
    case("M15", "RA7", "bin[2].cutoff_sensitivity.minus.total", ref("RA7", "bin[2].cutoff_sensitivity.plus.total"),
         "one-sided sensitivities exchanged")
    case("M05c", "RA10", "bin[0].single_loss_total",
         ref("RA10", "bin[0].single_loss_no_hole") + ref("RA10", "bin[0].single_loss_hole"),
         "total rebuilt as no_hole + Q_h*hole (0.72)")
    case("M17", "RA12", "bin[0].single_loss_total", Fr(0), "verification residual not propagated: satellite 0")
    case("M18", "RA9", "bin[0].single_loss_total.kernel_point_error", -ref("RA9", "bin[0].single_loss_total.kernel_point_error"),
         "signed weight: omega_2*L*e is negative (relative width 1e-9)")
    case("M18", "RA9", "bin[0].single_loss_total.kernel_point_error", Fr(0),
         "signed weights clipped at 0")
    return out


# ---------------------------------------------------------------- main

def build():
    raw = TABLE.read_bytes()
    if hashlib.sha256(raw).hexdigest() != TABLE_SHA256:
        raise AssertionError("the registration table changed")
    table = json.loads(raw)
    rows = {r["id"]: r for r in table["rows"]}
    out_rows, values, counts, no_reference = {}, {}, {}, []
    max_ratio = dict(value=Fr(0), row=None, key=None)
    notes = {}
    for row in table["rows"]:
        plan = [p.split("|") for p in row["check_plan"]]
        wanted = [(k, kind, oracle) for k, kind, oracle in plan
                  if kind in ("direct_reference", "inverse_interval", "error_part") and oracle.startswith("R-A")]
        if not wanted:
            continue
        rid = row["id"]
        extra = {}
        if row["entry"] in ("layer_ab", "layer_a") and rid != "RA14":
            low, high = at_precisions(lambda: compute_layer(table, row))
            notes[rid] = dict(balance=balance(table, row, high))
        elif rid == "RA9":
            (low, info_low), (high, info) = at_precisions(lambda: compute_ra9(table, row))
            extra = {k: info["node_rounding"] for k in high}
            notes[rid] = dict(omega_2_eV=decimal_string(exact(info["omega_2"])),
                              omega_2_difference_of_precisions=bound_string(abs(exact(info["omega_2"]) - exact(info_low["omega_2"]))),
                              cell_eV=[decimal_string(x) for x in info["cell"]],
                              cell_upper_hex=float(info["cell"][1]).hex(),
                              node_rounding_part_of_u_r=bound_string(exact(info["node_rounding"])))
        elif rid == "RA10":
            low, high = at_precisions(lambda: compute_ra10(table, row))
        elif rid == "RA14":
            (low, _), (high, info) = at_precisions(lambda: compute_ra14(table, row))
            notes[rid] = {k: decimal_string(exact(v)) for k, v in info.items()}
        else:
            raise AssertionError(rid)
        values[rid] = {k: (low[k], high[k]) for k in high}
        keys = {}
        area = None
        for key, kind, oracle in wanted:
            if kind == "inverse_interval":
                continue
            text, u_text, r, u = record(low[key], high[key], extra.get(key, Fr(0)))
            if kind == "error_part":
                d_abs, d_rel, unit = Fr(0), (Fr(1, 100) if rid == "RA14" else Fr(1, 10**9)), ("angstrom^-1" if rid == "RA14" else "I")
            else:
                d_abs, d_rel, unit = width_of(table, row, key)
            width = d_abs + d_rel*abs(r) if r else d_abs
            # prereg 6.3: rule 1 for r = 0, rule 3 for r != 0
            rule = 1 if r == 0 else 3
            valid = 4*u <= width
            if not valid:
                no_reference.append(dict(row=rid, key=key, rule=rule))
            if width and u/width > max_ratio["value"]:
                max_ratio = dict(value=u/width, row=rid, key=key)
            keys[key] = dict(kind=kind, status="reference" if valid else "no_reference", value=text, u_r=u_text,
                             d_abs=decimal_string(d_abs), d_rel=decimal_string(d_rel), width=bound_string(width),
                             unit=unit, zero_reference=r == 0, no_reference_rule_applied=rule)
            if key == "window_area.total":
                area = (r, u)
        for key, kind, oracle in wanted:
            if kind != "inverse_interval":
                continue
            r, u = area
            if r - u <= 0:
                raise AssertionError("inverse not defined")
            end = key.rsplit(".", 1)[1]
            bound = 1/(r + u) if end == "lower" else 1/(r - u)
            getcontext().prec = DIGITS
            keys[key] = dict(kind=kind, status="reference", reference_key="window_area.total", end=end,
                             value=decimal_string(bound), unit="angstrom", no_reference_rule_applied=None,
                             definition="1/(r + u_r) for lower, 1/(r - u_r) for upper, r and u_r of window_area.total "
                                        "of the same row; the test rebuilds the exact rational (prereg 7.2)")
        if list(keys) != [k for k, _, _ in wanted] and set(keys) != {k for k, _, _ in wanted}:
            raise AssertionError(rid)
        c = row_conditions(table, row)
        out_rows[rid] = dict(entry=row["entry"], keys={k: keys[k] for k, _, _ in wanted},
                             conditions=dict(loss_window_eV=[float(c["a"]).hex(), float(c["b"]).hex()],
                                             weight_length_A=float(c["L"]).hex(), hole_charge=float(c["charge"]).hex(),
                                             settings=c["settings"]),
                             notes=notes.get(rid))
        counts[rid] = {kind: sum(1 for _, k2, _ in wanted if k2 == kind)
                       for kind in ("direct_reference", "inverse_interval", "error_part")}
    numeric = sum(c["direct_reference"] + c["inverse_interval"] for c in counts.values())
    by_row = {rid: c["direct_reference"] + c["inverse_interval"] for rid, c in counts.items()
              if c["direct_reference"] + c["inverse_interval"]}
    ledger = {k: v for k, v in table["counts"]["planned_numeric_keys_by_row"].items() if k.startswith("RA")}
    if numeric != table["counts"]["planned_numeric_keys_by_family"]["RA"] or by_row != ledger:
        raise AssertionError((numeric, by_row, ledger))
    # RA13 check before the freeze
    k13 = out_rows["RA13"]["keys"]["window_area.term_13a"]
    r13, u13 = Fr(k13["value"]), Fr(k13["u_r"])
    width13 = Fr(k13["d_abs"]) + Fr(k13["d_rel"])*abs(r13)
    ra13 = ra13_prototype(table, rows["RA13"], (r13, u13), width13)
    if not ra13["normal"]["passes"] or ra13["truncated_after_initial_panels"]["passes"] or ra13["initial_panels"] != 6:
        raise AssertionError(ra13)
    wrong = known_wrong(table, rows, values)
    if not all(w["outside_width"] for w in wrong):
        raise AssertionError([w for w in wrong if not w["outside_width"]])
    # auxiliary analytic values used by the state test (not R-A ledger keys)
    mp.mp.dps = DPS[-1]
    wa = exact(values["RA4a"]["single_loss_weight"][1])
    wb = exact(values["RA4b"]["single_loss_weight"][1])
    n27 = rows["N27"]
    auxiliary = {
        "RA4b.single_loss_weight@difference": dict(
            value=decimal_string(wa - wb), definition="w(RA4a) - w(RA4b) with the binary64 inputs (0.06*ln(10) nominal)",
            width=bound_string(Fr(out_rows["RA4a"]["keys"]["single_loss_weight"]["width"])
                               + Fr(out_rows["RA4b"]["keys"]["single_loss_weight"]["width"]))),
        "N27.window_area.total.kernel_point_error": dict(
            value=decimal_string(sum((fd(y) - fd(x) for x, y in zip(n27["inputs"]["initial_edges_eV"],
                                                                       n27["inputs"]["initial_edges_eV"][1:])), Fr(0))*fd(1e-4)),
            definition="sum over the panels [2,4], [4,5] of |P| times the injected point error fl(1e-4)", d_rel="1e-12"),
        "RA1.satellite_sum": dict(
            value=decimal_string(sum((exact(values["RA1"][f"{n}.single_loss_total"][1])
                                      for n in ("outside.low", "outside.high", *(f"bin[{i}]" for i in range(7)))), Fr(0))),
            S_L_A_total=decimal_string(exact(values["RA1"]["single_loss_weight"][1])),
            definition="sum of single_loss_total over all regions, and S*L*A_total computed separately (0.60 nominal)"),
    }
    return {
        "schema": "xpsloss-m4a-references-ra-v1",
        "registration_table_sha256": TABLE_SHA256,
        "method": "exact rationals of the binary64 inputs; closed forms with mpmath at %d and %d digits; "
                  "u_r = difference of the two precisions + decimal rounding (+ node rounding for RA9)" % DPS,
        "decimal_digits": DIGITS,
        "conventions": {
            "cutoff_sensitivity": "H_j = L*a*K_j(a+0), the limit from inside the window",
            "one_sided": "plus: T -> a+0, minus: T -> a-0; H = L*a*K_j(a+-0)*G_D(a+-0)",
            "bins": "left closed, right open; the last bin contains its right end",
            "inverse_interval": "no reference of its own; taken from window_area.total of the row",
            "width": "d_abs + d_rel*|r| for r != 0, d_abs for r = 0 (prereg 7.2); RA6, RA6b use the replacement strengths",
        },
        "counts": {"by_row": counts, "numeric_keys": numeric, "numeric_keys_by_row": by_row,
                   "error_part_keys": sum(c["error_part"] for c in counts.values()),
                   "ledger_numeric_keys_RA": table["counts"]["planned_numeric_keys_by_family"]["RA"],
                   "zero_references": sum(1 for r in out_rows.values() for k in r["keys"].values() if k.get("zero_reference"))},
        "max_u_r_over_width": dict(value=bound_string(max_ratio["value"]), row=max_ratio["row"], key=max_ratio["key"]),
        "no_reference": no_reference,
        "reference_path_failures": [],
        "ra13_prefreeze_check": ra13,
        "known_wrong_answers": wrong,
        "auxiliary": auxiliary,
        "rows": out_rows,
    }


def main(argv):
    start = time.perf_counter()
    text = json.dumps(build(), indent=1, sort_keys=True, ensure_ascii=True) + "\n"
    data = text.encode("ascii")
    if "--check" in argv:
        same = OUT.read_bytes() == data
        print("check:", "identical" if same else "DIFFERENT")
        return 0 if same else 1
    OUT.write_bytes(data)
    print("wrote", OUT.name, "sha256", hashlib.sha256(data).hexdigest(), "seconds %.2f" % (time.perf_counter() - start))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
