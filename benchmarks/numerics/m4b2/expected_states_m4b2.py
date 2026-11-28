"""M4b-2 pre-registration: the expected-state table (docs/M4B2_PREREG.md section 5).  Fourth draft.

Writes expected_states_m4b2.json next to this file (ASCII, LF) and prints its SHA-256 and the counts.
One row is one new entry point call (entry, inputs, settings) with its own expected terminal state.
"in_row_checks" are control runs and comparisons made inside a row; "reuse" names artifacts of other
rows ("ROW.result").  Neither adds a row.

Standard library only (not xpsloss, not NumPy).  Kernels, the KE primary, the KE output edges and the
materials are copied from the frozen M4b-1 table after its SHA-256 is checked.

What this script proves with exact rationals (fractions) for every binding-energy row, and writes under
"geometry": the conversion C - B with one rounding, the union injectivity test with the first colliding
pair and its roles, the support orders S_B, S_K and S = max, the per-bin loss widths of both sides and
their maximum rounded up, the bins whose bound is "complete" (both sides below (N+1)a), the counts of
disagreeing comparisons, and the zero loss per region of the rounded-KE problem.  What it cannot prove
(a finite bound against a tolerance, convergence of a public run) is written under "expect" with a
"basis" that names where the expectation comes from.

The ledger is expanded here into concrete (row, key, kind, oracle) records.  The JSON keeps the groups
with their counts and three example keys; the expansion is checked for duplicates, every mutation
target must be a ledger record of its row, and every CLI argv is parsed against the option arities.

Usage: python benchmarks/numerics/m4b2/expected_states_m4b2.py          (write the JSON)
       python benchmarks/numerics/m4b2/expected_states_m4b2.py --check  (compare byte for byte)
"""
import copy
import hashlib
import json
import math
import pathlib
import sys
from fractions import Fraction as Q

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "expected_states_m4b2.json"
M4B1 = HERE.parent / "m4b1" / "expected_states_m4b1.json"
M4B1_SHA256 = "d8d0f892a49c3ed9933b08d6a7c92e5b3aaf05ee5d6c7575050fc164991035e2"
FILES = {"benchmarks/inputs/al_elfdb_v1.0.json": "17ded0465a5e42124eebf44d9573369a5cbea87df59b5df11216ea2f7432d1c4",
         "benchmarks/inputs/al_elfdb_v1.0.bulk.json":
             "978931c1651881e6c7dfa739f4d626e4496225bf24bceb815466ff651934eae2"}
BASELINE = "8c33d2939dab866decc37fdf9663f5ac2e29b78d"
MODEL = "first_loss_13a_chain_v0"
SCHEMA_KE = "xpsloss-single-loss-spectrum-v0"
SCHEMA_BE = "xpsloss-measured-spectrum-v0"
N_ROWS, N_CI, N_FULL = 82, 70, 12
STATES = {"success": 40, "complete": 8, "error": 32, "not_converged": 2}
NR = "not_recorded"
COLUMNS = ["term_13a", "term_13b", "term_14", "term_15_T1", "term_15_T2", "term_15_T3", "term_16_e", "term_16_h",
           "no_hole", "hole", "total"]
COMPONENTS = ["total", "no_hole", "hole"]
PARTS = ["inherited", "propagated", "local", "multiplier_interval"]


class TableError(RuntimeError):
    pass


def need(cond, message):
    if not cond:
        raise TableError(message)


if sys.argv[1:] not in ([], ["--check"]):
    raise SystemExit("usage: expected_states_m4b2.py [--check]")

ROOT = next((p for p in HERE.parents if (p / "benchmarks/inputs/al_elfdb_v1.0.json").is_file()), None)
need(ROOT is not None, "repository root with benchmarks/inputs not found")
for name, digest in FILES.items():
    need(hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest, f"{name} changed")
need(hashlib.sha256(M4B1.read_bytes()).hexdigest() == M4B1_SHA256, "the frozen M4b-1 table changed")
OLD = json.loads(M4B1.read_text(encoding="utf-8"))


def num(x):
    x = float(x)
    return {"value": x, "hex": x.hex()}


def nums(xs):
    xs = [float(x) for x in xs]
    return {"values": xs, "hex": [x.hex() for x in xs]}


def cmp(x, y):
    return (x > y) - (x < y)


REL = {-1: "<", 0: "=", 1: ">"}


def up(x):
    """Smallest binary64 that is >= the rational x."""
    f = float(x)
    return math.nextafter(f, math.inf) if Q(f) < x else f


def frac(x):
    return f"{x.numerator}/{x.denominator}"


U40 = 2.0**-40
A_EPS = 1.0 + 2.0**-45
need(float("1.0000000000000284") == A_EPS, "decimal form of 1 + 2^-45")

# ---------------------------------------------------------------- inputs
KERNELS = {k: OLD["kernels"][k] for k in ("K_A1s", "K_A3", "K_T1", "K_T2")}
KERNELS["K_T1e"] = {
    "definition": "K_T1 with the lower window edge 1 + 2^-45 eV (design section 5 counterexample)",
    "columns": {"term_13a": [{"type": "box", "value": num(0.2), "interval_eV": [num(A_EPS), num(2.5)]}]},
    "features": {"jumps_eV": [A_EPS, 2.5], "kinks_eV": [], "singularities": []},
    "window_eV": [num(A_EPS), num(2.5)],
    "derived_columns": KERNELS["K_T1"]["derived_columns"],
}
WINDOW = {k: tuple(w["value"] for w in v["window_eV"]) for k, v in KERNELS.items()}
MATERIALS = {k: OLD["materials"][k] for k in ("Al", "synthetic_wide")}
KE_PRIMARIES = {"D1000": OLD["primaries"]["D1000"]}
KE_EDGES = {k: OLD["output_edges"][k] for k in ("E_M", "E_A3")}
KE_EDGES["E_Ms"] = {"values": [995.0, 999.0, 1000.0, 1001.0],
                    "hex": [x.hex() for x in (995.0, 999.0, 1000.0, 1001.0)]}

REC_FULL = {
    "source": {"name": "synthetic monochromatic source", "monochromatic": True,
               "satellite_treatment": "none (synthetic)"},
    "calibration": {"line": "synthetic line", "reference_binding_eV": 84.0, "measured_binding_eV": 84.0,
                    "applied_shift_eV": 0.0},
    "charge_correction": {"applied": False, "method": None, "shift_eV": None},
    "core_level": "synthetic level", "emission_angle_deg": 0.0, "acceptance_half_angle_deg": 5.0}
RECORD_FIELDS = list(REC_FULL)


def record(**changes):
    out = copy.deepcopy(REC_FULL)
    for key, value in changes.items():
        parent, _, child = key.partition("__")
        if child:
            out[parent][child] = value
        else:
            out[parent] = value
    return out


RECORDS = {  # all values are synthetic; these dictionaries are exactly the six record fields
    "REC_FULL": REC_FULL,
    "REC_NR": {k: NR for k in RECORD_FIELDS},
    "REC_TILT": record(emission_angle_deg=30.0),
    "REC_NAMED": record(source__name="Al Kalpha monochromatic"),
    "REC_CORR": record(calibration={"line": "synthetic line", "reference_binding_eV": 84.0,
                                    "measured_binding_eV": 84.25, "applied_shift_eV": -0.25},
                       charge_correction={"applied": True, "method": "synthetic flood-gun reference",
                                          "shift_eV": 1.5}),
    "REC_NEST": record(source__monochromatic=NR),
}
NOT_RECORDED_ORDER = ["reference.line", "reference.citation", "primary_reference_choice",
                      "source", "source.name", "source.monochromatic", "source.satellite_treatment",
                      "calibration", "calibration.line", "calibration.reference_binding_eV",
                      "calibration.measured_binding_eV", "calibration.applied_shift_eV",
                      "charge_correction", "charge_correction.applied", "charge_correction.method",
                      "charge_correction.shift_eV", "core_level", "emission_angle_deg",
                      "acceptance_half_angle_deg"]


def po(hv, phi):
    return {"photon_energy_eV": num(hv),
            "reference": {"kind": "photon_and_offset", "kinetic_energy_offset_eV": num(phi)}}


def pa(hv, b, e):
    return {"photon_energy_eV": num(hv),
            "reference": {"kind": "paired_anchor", "binding_anchor_eV": num(b), "kinetic_anchor_eV": num(e),
                          "line": "synthetic anchor line", "citation": "synthetic"}}


AXES = {  # all synthetic
    "AX_PO": po(1486.6, 4.5),            # C = Q(1486.6) - 9/2, a binary64
    "AX_PA": pa(1486.6, 84.0, 1398.25),  # C = 1482.25, implied offset Q(1486.6) - 1482.25
    "AX_R": po(1486.6, 4.6),             # C is not a binary64 (double rounding differs)
    "AX_D": po(1253.5, 4.25),            # C = 1249.25 (KE images equal M4b-1 inputs)
    "AX_DA": pa(1253.5, 249.25, 1000.0),  # C = 1249.25, implied offset 4.25
    "AX_S": po(130.5, 2.5),              # C = 128 (one-ulp edges stay exact)
    "AX_T2": po(1069.0, 4.5),            # C = 1064.5 (KE image equals the M4b-1 row T2c)
    "AX_H": po(6000.0, 0.0),             # synthetic HAXPES
    "AX_8192": po(8192.0, 0.0),          # design counterexamples
    "AX_16": po(16.0, 0.0),              # exact control of the counterexamples
}


def delta(positions, intensities):
    return {"kind": "delta", "binding_energies_eV": nums(positions), "intensities": nums(intensities)}


def binned(edges, intensities):
    return {"kind": "binned", "binding_edges_eV": nums(edges), "intensities": nums(intensities)}


BE_PRIMARIES = {
    "BD_2": delta([72.5, 60.5], [1.0, 0.5]),
    "BD_2r": delta([60.5, 72.5], [0.5, 1.0]),
    "BB_neg": binned([1.5, 0.5, 0.0, -0.5, -1.5], [0.1, 0.4, 0.3, 0.2]),
    "BB_neg_inc": binned([-1.5, -0.5, 0.0, 0.5, 1.5], [0.2, 0.3, 0.4, 0.1]),
    "BD_ends": delta([70.0, 74.0, 76.0], [0.25, 0.5, 0.25]),
    "BD_dup": delta([74.0, 74.0, 70.0], [0.25, 0.25, 0.5]),
    "BD_731": delta([73.1], [1.0]),
    "BD_723": delta([72.3], [1.0]),
    "BD_249": delta([249.25], [1.0]),
    "BD_249zero": delta([249.25], [0.0]),
    "BD_645": delta([64.5], [1.0]),
    "BD_0": delta([0.0], [1.0]),
    "BD_eps": delta([1.0 + 2.0**-45], [1.0]),
    "BD_eps0": delta([1.0 + 2.0**-45], [0.0]),
    "BB_u": binned([0.0, 2 * U40], [1.0]),
    "BD_m1": delta([-1.0], [1.0]),
}
E_B120 = [880.5 + i for i in range(121)]
UP1 = math.nextafter(67.5, math.inf)
DN1 = math.nextafter(67.5, -math.inf)
BE_EDGES = {
    "BE_E1": nums([70, 72.5, 73, 74.5, 75.5, 77.5, 80, 81]),
    "BE_E1r": nums([81, 80, 77.5, 75.5, 74.5, 73, 72.5, 70]),
    "BE_Eneg": nums([14, 10, 7, 4, 2.5, 1, 0.25, -0.75]),
    "BE_E4": nums([70, 72, 74, 76]),
    "BE_E5a": nums([70.3, 73.1, 75.2, 77.7, 80.9, 85.4]),
    "BE_E5": nums([113.7, 93.4, 89.9, 85.2, 78.3, 74.6, 73.1, 72.4]),
    "BE_E8": nums([70.1, 72.3, 74.6, 77.8, 80.3, 85.6]),
    "BE_EM": nums([248.25, 249.25, 250.25, 254.25, 261.25, 265.25, 269.25, 289.25]),
    "BE_EMr": nums([289.25, 269.25, 265.25, 261.25, 254.25, 250.25, 249.25, 248.25]),
    "BE_EMs": nums([248.25, 249.25, 250.25, 254.25]),
    "BE_B120": {**nums([1249.25 - x for x in E_B120]),
                "note": "1249.25 - (880.5 + i), i = 0..120: decreasing BE, the KE image is the M4b-1 E_B120"},
    "BE_T1_eq": nums([63.5, 64.5, 65.5, 66.5, 67.5]),
    "BE_T1_dn": nums([63.5, 64.5, 65.5, 66.5, DN1]),
    "BE_T1_up": nums([63.5, 64.5, 65.5, 66.5, UP1]),
    "BE_T2": nums([63.5, 64.5, 65.5, 74.5, 84.5, 94.5, 114.5]),
    "BE_c7": nums([0, 1, 2]),
    "BE_c5": nums([0, 2.0 + 2.0**-43]),
    "BE_c5b": nums([0, 2.0 + 2.0**-43, 9]),
    "BE_c5rev": nums([0, 2.0 - 2.0**-43]),
    "BE_c5q": nums([0, 2.25 + 2.0**-43]),
    "BE_u": nums([0, 0.75 * U40, 2 * U40]),
    "BE_grid": {**nums([-2 + 0.5 * i for i in range(9)]),
                "compact": {"min": -2, "max": 2, "step": 0.5},
                "note": "the CLI row gives the compact form; expansion as the existing _edges (exact decimal "
                        "interval count, then numpy.linspace)"},
    "BE_neg_ke": nums([126, 130]),
    "BE_bad": nums([70, 72, 71]),
}


# ---------------------------------------------------------------- exact reading of a BE input
def constant(axis):
    ax = AXES[axis]
    ref = ax["reference"]
    if ref["kind"] == "photon_and_offset":
        return Q(ax["photon_energy_eV"]["value"]) - Q(ref["kinetic_energy_offset_eV"]["value"])
    return Q(ref["binding_anchor_eV"]["value"]) + Q(ref["kinetic_anchor_eV"]["value"])


def analyse(axis, primary, edges, window, max_order, ref=None, e0=None):
    """Exact-rational reading of one BE input (no product code).  Rules: prereg sections 2.1-2.5."""
    C = constant(axis)
    p, e = BE_PRIMARIES[primary], BE_EDGES[edges]["values"]
    field = "binding_energies_eV" if p["kind"] == "delta" else "binding_edges_eV"
    pos, inten = p[field]["values"], p["intensities"]["values"]
    role = {}
    for k, x in enumerate(pos):
        role.setdefault(x, f"primary.{field}[{k}]")
    for k, x in enumerate(e):
        role.setdefault(x, f"output_edges_eV[{k}]")
    coords = sorted(role)
    ke = {x: float(C - Q(x)) for x in coords}
    out = {"C": C, "ke": ke, "collision": None,
           "per_array_injective": all(len({ke[x] for x in arr}) == len(set(arr)) for arr in (pos, e)),
           "negative": min(ke.values()) < 0}
    for x, y in zip(coords, coords[1:]):
        if ke[x] <= ke[y]:
            out["collision"] = {"binding_eV": nums([x, y]), "kinetic_eV": nums([ke[x], ke[y]]),
                                "roles": [role[x], role[y]]}
            return out
    rounding = max(abs(Q(ke[x]) - (C - Q(x))) for x in coords)
    out.update(rounding=rounding, exact=rounding == 0,
               double_rounding_differs=sum(float(Q(float(C)) - Q(x)) != ke[x] for x in coords))
    if p["kind"] == "delta":
        support = sorted({x for x, f in zip(pos, inten) if f})
    else:
        support = sorted({x for a, b, f in zip(pos, pos[1:], inten) if f for x in (a, b)})
    a, b = Q(window[0]), Q(window[1])
    nb = len(e) - 1
    increasing = e[0] < e[-1]
    high = [max(e[i], e[i + 1]) for i in range(nb)]   # high-BE edge of input bin i
    if support:
        ub = [Q(d) - Q(min(support)) for d in high]
        uk = [Q(ke[min(support)]) - Q(ke[d]) for d in high]
        need(ke[min(support)] == max(ke[x] for x in support), "rounding is monotone")
        sb, sk = (max(0, math.floor(max(u) / a)) for u in (ub, uk))
    else:
        ub = uk = None
        sb = sk = 0
    s = max(sb, sk)
    n = s if max_order == "display_support" else max_order
    geometry = {"support_order_binding": sb, "support_order_kinetic": sk, "support_order": s, "computed": n,
                "conversion_exact": rounding == 0, "conversion_constant_exact": frac(C)}
    if support:
        nxt = (n + 1) * a
        geometry["region_loss_upper_eV"] = [num(up(max(x, y))) for x, y in zip(ub, uk)]
        geometry["widths_equal_on_both_sides"] = ub == uk
        geometry["bound_class"] = ["complete" if x < nxt and y < nxt else "reached" for x, y in zip(ub, uk)]
        geometry["kinetic_side_reaches_next_order"] = [y >= nxt for y in uk]
        geometry["width_excess_over_kinetic_eV"] = num(up(max(max(x, y) - y for x, y in zip(ub, uk))))
        geometry["width_equals_kinetic"] = [up(max(x, y)) == up(y) for x, y in zip(ub, uk)]

        def example(kind, order, name, x, d, rb, rk):
            return {"comparison": kind, "order": order, "bound": name,
                    "primary_binding_eV": num(x), "primary_kinetic_eV": num(ke[x]), "primary_index": pos.index(x),
                    "edge_binding_eV": num(d), "edge_kinetic_eV": num(ke[d]), "edge_index": e.index(d),
                    "binding_relation": REL[rb], "kinetic_relation": REL[rk]}

        # exhaustive search (the product counts by bisection); order: order, lower/upper, primary BE, edge BE
        examples = []
        for order in range(1, n + 1):
            for name, c in (("lower", order * a), ("upper", order * b)):
                for x in support:
                    for d in sorted(e):
                        rb, rk = cmp(Q(d) - Q(x), c), cmp(Q(ke[x]) - Q(ke[d]), c)
                        if rb != rk:
                            examples.append(example("computed_order", order, name, x, d, rb, rk))
        n_computed = len(examples)
        next_rel, next_reach = [], []
        for d, i in sorted((high[i], i) for i in range(nb)):
            rb, rk = cmp(ub[i], nxt), cmp(uk[i], nxt)
            if rb != rk:
                next_rel.append(i)
                examples.append({**example("next_order", n + 1, "lower", min(support), d, rb, rk), "bin_index": i})
            if (ub[i] >= nxt) != (uk[i] >= nxt):
                next_reach.append(i)
        geometry["support_comparison"] = {
            "computed_relation_mismatch_count": n_computed, "next_relation_mismatch_count": len(next_rel),
            "next_reachability_mismatch_bins": sorted(next_reach), "first_mismatches": examples[:8]}
        disagree = sb != sk or bool(n_computed) or bool(next_rel)
    else:
        geometry["region_loss_upper_eV"] = [None] * nb
        geometry["widths_equal_on_both_sides"] = True
        geometry["bound_class"] = ["complete"] * nb
        geometry["kinetic_side_reaches_next_order"] = [False] * nb
        geometry["width_equals_kinetic"] = [True] * nb
        geometry["width_excess_over_kinetic_eV"] = num(0.0)
        geometry["support_comparison"] = {"computed_relation_mismatch_count": 0, "next_relation_mismatch_count": 0,
                                          "next_reachability_mismatch_bins": [], "first_mismatches": []}
        disagree = False
    geometry["support_disagreement_warning"] = disagree
    geometry["truncation_equals_ke_control"] = not disagree and geometry["widths_equal_on_both_sides"]
    geometry["state_without_tolerance"] = "complete_on_display" if n >= s else "truncated"
    # zero loss of the rounded-KE problem (KE bins [e_i, e_i+1), the last one closed), keyed in input order
    asc = sorted(Q(ke[d]) for d in e)
    bounds = [None, *asc, None]
    names = ["outside.high_binding",
             *[f"bin[{nb - 1 - k if increasing else k}]" for k in range(nb)], "outside.low_binding"]
    mass = dict.fromkeys(names, Q(0))
    overlap = dict.fromkeys(names, Q(0))
    if p["kind"] == "delta":
        for x, f in zip(pos, inten):
            t = Q(ke[x])
            k = 0 if t < asc[0] else nb + 1 if t > asc[-1] else nb if t == asc[-1] else \
                next(j for j in range(nb) if t < asc[j + 1]) + 1
            mass[names[k]] += Q(f)
            overlap[names[k]] += Q(f)
    else:
        basc = sorted(map(Q, e))
        bbounds = [None, *basc, None]
        bnames = ["outside.low_binding",
                  *[f"bin[{j if increasing else nb - 1 - j}]" for j in range(nb)], "outside.high_binding"]
        for s0, s1, f in zip(pos, pos[1:], inten):
            for table, grid, labels, lo, hi in (
                    (mass, bounds, names, *sorted((Q(ke[s0]), Q(ke[s1])))),
                    (overlap, bbounds, bnames, *sorted((Q(s0), Q(s1))))):
                for name, left, right in zip(labels, grid, grid[1:]):
                    start = lo if left is None else max(lo, left)
                    stop = hi if right is None else min(hi, right)
                    if start < stop:
                        table[name] += Q(f) * (stop - start) / (hi - lo)
    geometry["zero_loss_kinetic_widths"] = {k: frac(v) for k, v in sorted(mass.items())}
    if overlap != mass:
        geometry["zero_loss_exact_binding_overlap_diagnostic_only"] = {k: frac(v) for k, v in sorted(overlap.items())}
    if ref is not None:
        geometry["implied_offset_is_none"] = AXES[axis]["reference"]["kind"] == "photon_and_offset"
        if e0 is not None:
            geometry["electron_energy_minus_reference_exact"] = frac(Q(e0) - (C - Q(ref)))
    out["geometry"] = geometry
    out["n_coordinates"] = {"edges": len(e), "primary": len(pos)}
    return out


# ---------------------------------------------------------------- ledger
ROWS = []
LEDGER = {}


def add(rid, key, kind, oracle):
    record_ = (key, kind, oracle)
    need(record_ not in LEDGER.setdefault(rid, []), f"{rid}: duplicate ledger record {record_}")
    LEDGER[rid].append(record_)


def regions(nb):
    return ["outside.low_binding", *[f"bin[{i}]" for i in range(nb)], "outside.high_binding"]


def physical_ledger(rid, nb, orders, oracle):
    """Keys compared bit for bit with the physical control (prereg 7.1)."""
    for n in range(1, orders + 1):
        for region in regions(nb):
            for c in COMPONENTS:
                add(rid, f"order[{n}].{region}.{c}", "bit_identity", oracle)
                add(rid, f"order[{n}].{region}.{c}_error", "bit_identity", oracle)
                for part in PARTS:
                    add(rid, f"order[{n}].{region}.{c}_error_parts[{part}]", "bit_identity", oracle)
        for col in COLUMNS:
            for name in ("area", "area_error", "cutoff_sensitivity", "cutoff_sensitivity_error", "area_zero_kind"):
                add(rid, f"order[{n}].{name}.{col}", "bit_identity", oracle)
            for part in (*PARTS, "convention"):
                add(rid, f"order[{n}].area_error_parts[{part}][{col}]", "bit_identity", oracle)
        add(rid, f"order[{n}].mesh", "record_identity", oracle)
    for region in regions(nb):
        for name in ("zero_loss", "zero_loss_error", "spectrum", "spectrum_error"):
            add(rid, f"{name}.{region}", "bit_identity", oracle)
    for name in ("order_area_ratio.value", "order_area_ratio.error", "order_area_ratio.interval[0]",
                 "order_area_ratio.interval[1]", "convention_ratio"):
        add(rid, name, "bit_identity", oracle)
    for name in ("single_loss (every field; the KE result inside)", "truncation.neglected_total_area",
                 "diagnostics (parts that do not depend on the axis; the tolerance entry of "
                 "effective_settings is checked against the input, not against the control)",
                 "warnings (identifiers and order without the BE identifiers and without 'truncated')"):
        add(rid, name, "record_identity", oracle)


def truncation_ledger(rid, nb, geometry, tolerance):
    rc = "R-C (both-sided rule)"
    # with a tolerance the state depends on the numerical bound, which R-C does not produce: oracle "expect"
    add(rid, "truncation.state", "truncation", rc if tolerance is None else "expect")
    for name in ("truncation.support_order", "truncation.computed", "truncation.requested",
                 "energy_axis.support_order_binding", "energy_axis.support_order_kinetic",
                 "energy_axis.support_comparison.computed_relation_mismatch_count",
                 "energy_axis.support_comparison.next_relation_mismatch_count",
                 "energy_axis.support_comparison.next_reachability_mismatch_bins",
                 "warning axis_conversion_support_disagreement present or absent"):
        add(rid, name, "truncation", rc)
    add(rid, "energy_axis.support_comparison.first_mismatches", "record_identity", rc)
    add(rid, "warning 'truncated' present if and only if truncation.state == 'truncated'", "truncation",
        "same result")
    if tolerance is not None:
        add(rid, "truncation.neglect_tolerance_fraction and diagnostics.effective_settings entry == input",
            "truncation", "inputs")
    for i in range(nb):
        add(rid, f"energy_axis.region_loss_upper_eV[bin[{i}]]", "truncation", rc)
        # judged from the three returned fields: complete <=> exactly (0.0, 'complete', None)
        add(rid, f"truncation.bound_class[bin[{i}]]", "truncation", rc)
        if geometry["bound_class"][i] != "reached":
            continue
        if not geometry["kinetic_side_reaches_next_order"][i]:
            add(rid, f"truncation three fields of bin[{i}]: finite and > 0, reason None, finite tau",
                "truncation", "expect")
        elif geometry["width_equals_kinetic"][i]:
            for name in ("neglected_at_or_above_region_bound", "bound_reason", "tau_per_eV"):
                add(rid, f"truncation.{name}[bin[{i}]]", "bit_identity", "KE control")
        else:
            add(rid, f"truncation.neglected_at_or_above_region_bound[bin[{i}]]", "bound_bracket",
                "KE control: b_K <= b <= the candidate at tau_K with the R-C width (prereg 7.2)")


def conversion_ledger(rid, facts, public):
    rc = "R-C"
    for k in range(facts["n_coordinates"]["edges"]):
        add(rid, f"energy_axis.kinetic_output_edges_eV[{k}]", "exact_conversion", rc)
    for k in range(facts["n_coordinates"]["primary"]):
        add(rid, f"energy_axis.kinetic_primary.coordinates[{k}]", "exact_conversion", rc)
    for name in ("conversion_constant_eV", "conversion_constant_exact", "implied_offset_eV", "rounding_max_eV",
                 "primary_reference_kinetic_eV"):
        add(rid, f"energy_axis.{name}", "exact_conversion", rc)
    add(rid, "energy_axis.electron_energy_minus_reference_eV", "exact_conversion", rc if public else "None (analytic)")
    for name in ("input_binding_edges_eV (input bits and order)", "input_primary (input bits and order)",
                 "direction", "record (deep copy of the six fields)", "not_recorded (names and order)",
                 "axis, binding_reference, reference_kind and the copied reference fields"):
        add(rid, f"energy_axis.{name}", "record_identity", "inputs")
    add(rid, "warnings (BE identifiers, their order, source and message pairing)", "diagnostic", "prereg 4.4")


def membership_ledger(rid, nb):
    for region in regions(nb):
        add(rid, f"zero_loss.{region}", "bin_membership", "R-C (KE widths)")


# ---------------------------------------------------------------- rows
PUBLIC = {"electron_energy_eV": 1000.0, "weight_length_A": 20.0, "epsilon_momentum": "k_total",
          "loss_window_eV": [num(1.0), num(50.0)]}
ANALYTIC = {"weight_length_A": 2.0, "hole_charge": 1.0}
DEFAULT_CHOICE = "delta position (synthetic)"


def row(rid, group, entry, inputs, expect, *, ci=True, settings=None, checks=(), reuse=(), provides=(),
        denominator="state", geometry=None):
    add(rid, "terminal_state", "state", "expect")
    item = {"id": rid, "group": group, "entry": entry, "inputs": inputs, "settings": settings or {},
            "expect": expect, "ci": ci, "reuse": list(reuse), "provides": list(provides),
            "in_row_checks": list(checks), "denominator": denominator, "mutations": []}
    if geometry is not None:
        item["geometry"] = geometry
    ROWS.append(item)


def axis_in(axis, rec, ref, choice=DEFAULT_CHOICE, defects=()):
    out = {"axis": axis, "record": rec, "primary_reference_binding_eV": num(ref),
           "primary_reference_choice": choice}
    if defects:
        out["defects"] = list(defects)
    return out


def n_bins(edges):
    return len(BE_EDGES[edges]["values"]) - 1


def an(kernel, **extra):
    return {"kernel": kernel, **ANALYTIC, **extra}


def be_ok(rid, group, axis, rec, primary, edges, max_order, *, ref, kernel=None, public=None, state,
          support, settings=None, checks=(), ci=True, reuse=(), warnings=(), basis="geometry", expect_extra=None,
          control_order=None, choice=DEFAULT_CHOICE):
    """A successful BE row.  `support` is the registered (S_B, S_K); `state` the registered state."""
    entry = "api_spectrum" if public else "analytic_spectrum"
    window = WINDOW[kernel] if kernel else tuple(w["value"] for w in public["loss_window_eV"])
    e0 = public["electron_energy_eV"] if public else None
    facts = analyse(axis, primary, edges, window, max_order, ref=ref, e0=e0)
    need(facts["collision"] is None and not facts["negative"], f"{rid}: not a valid BE input")
    g = facts["geometry"]
    need((g["support_order_binding"], g["support_order_kinetic"]) == tuple(support),
         f"{rid}: support orders {g['support_order_binding']}, {g['support_order_kinetic']}")
    tolerance = (settings or {}).get("neglect_tolerance_fraction")
    if tolerance is None:
        need(state == g["state_without_tolerance"], f"{rid}: state {g['state_without_tolerance']}")
    else:
        need(basis != "geometry" and state in ("bounded_on_display",), f"{rid}: a tolerance row needs a basis")
    nb = n_bins(edges)
    computed = g["computed"]
    inputs = {"energy_axis": axis_in(axis, rec, ref, choice), "primary": primary, "output_edges_eV": edges,
              "model": MODEL, "max_order": max_order}
    inputs.update(public or an(kernel))
    codes = ["binding_axis_assumed_calibrated", *warnings]
    if g["support_disagreement_warning"]:
        codes.append("axis_conversion_support_disagreement")
    expect = {"status": "success", "truncation_state": state, "support_order": g["support_order"],
              "computed": computed, "basis": basis, "warnings_include": codes,
              "warnings_exclude": [] if g["support_disagreement_warning"] else ["axis_conversion_support_disagreement"]}
    expect.update(expect_extra or {})
    if control_order is None:
        control_order = "display_support" if (max_order == "display_support" and g["support_order"] == 0) \
            else computed
    control = (f"physical control: the same entry on the KE coordinates of the reference with max_order = "
               f"{control_order!r} and no tolerance; the keys of prereg 7.1 bit for bit")
    conversion_ledger(rid, facts, bool(public))
    membership_ledger(rid, nb)
    physical_ledger(rid, nb, computed, "KE control")
    truncation_ledger(rid, nb, g, tolerance)
    row(rid, group, entry, inputs, expect, ci=ci, settings=settings, checks=[control, *checks], reuse=reuse,
        provides=["result"], denominator="bit_identity_and_exact", geometry=g)
    return facts


# --- BE: conversion, orientation, membership
be_ok("BE1", "BE", "AX_PO", "REC_FULL", "BD_2", "BE_E1", 3, ref=72.5, kernel="K_A1s", support=(10, 10),
      state="truncated",
      checks=["both outside regions are nonzero: outside.low_binding holds the zero loss 1/2 of the 60.5 eV line "
              "and its first order; outside.high_binding holds higher orders of the 72.5 eV line",
              "the delta at 72.5 eV lies on an interior edge: its zero loss is in bin[0] (lower BE side)"])
be_ok("BE1r", "BE", "AX_PO", "REC_FULL", "BD_2r", "BE_E1r", 3, ref=72.5, kernel="K_A1s", support=(10, 10),
      state="truncated", reuse=["BE1.result"],
      checks=["decreasing input: only the BE display-bin index is reversed against BE1.result; kinetic_output_edges_eV "
              "and kinetic_primary are in the same ascending KE order; the outside names do not move"])
add("BE1r", "physical keys of prereg 7.1 of BE1.result with BE display bin[i] -> bin[6 - i]; KE coordinates "
    "in the same ascending order", "record_identity", "BE1.result")
be_ok("BE2", "BE", "AX_PA", "REC_FULL", "BD_2", "BE_E1", 3, ref=72.5, kernel="K_A1s", support=(10, 10),
      state="truncated", reuse=["BE1.result"],
      checks=["paired anchor: C = 1482.25 exactly; implied offset = RN(Q(1486.6) - 1482.25)",
              "the 90 registered values (3 orders x 9 regions x 3 components, and spectrum in 9 regions) are bit "
              "for bit those of BE1.result; this is registered for this pair only, not for common shifts in general"])
for region_ in regions(7):
    add("BE2", f"spectrum.{region_}", "bit_identity", "BE1.result")
    for n_ in (1, 2, 3):
        for c_ in COMPONENTS:
            add("BE2", f"order[{n_}].{region_}.{c_}", "bit_identity", "BE1.result")
be_ok("BE3", "BE", "AX_PA", "REC_FULL", "BB_neg", "BE_Eneg", 2, ref=0.0, kernel="K_A3", support=(15, 15),
      state="truncated", choice="centre of the synthetic band",
      checks=["binned primary given in decreasing BE with negative BE; unequal bins; the part of the primary "
              "below -0.75 eV is zero loss of outside.low_binding; negative hole bins exist and are not clipped"])
be_ok("BE3r", "BE", "AX_PA", "REC_FULL", "BB_neg_inc", "BE_Eneg", 2, ref=0.0, kernel="K_A3", support=(15, 15),
      state="truncated", reuse=["BE3.result"], choice="centre of the synthetic band",
      checks=["the same primary given in increasing BE with the intensities reversed together: every numeric key "
              "equals BE3.result; energy_axis.input_primary keeps the input order"])
add("BE3r", "physical keys of prereg 7.1, truncation lists and KE coordinates of BE3.result (the input BE edge "
    "array of the primary differs on purpose)", "record_identity", "BE3.result")
be_ok("BE4", "BE", "AX_PO", "REC_FULL", "BD_ends", "BE_E4", 2, ref=74.0, kernel="K_A1s", support=(3, 3),
      state="truncated",
      checks=["deltas on the lowest edge, an interior edge and the highest edge: zero loss 1/4, 1/2, 1/4 in "
              "bin[0], bin[1], bin[2]; both outside regions hold no zero loss"])
f5a = be_ok("BE5a", "BE", "AX_H", "REC_FULL", "BD_731", "BE_E5a", 2, ref=73.1, kernel="K_A1s", support=(6, 6),
            state="truncated",
            checks=["decimal BE at C = 6000: the conversion rounds; rounding_max_eV > 0 equals the reference",
                    "the delta and the edge 73.1 are the same binary64 (coordinates equal from the start are merged)"])
need(not f5a["exact"], "BE5a must round")
f5 = be_ok("BE5", "BE", "AX_H", "REC_FULL", "BD_731", "BE_E5", 2, ref=73.1, ci=False,
           public={"material": "synthetic_wide", **PUBLIC, "electron_energy_eV": 5926.9}, support=(40, 40),
           state="truncated", warnings=["E0_above_5keV"],
           checks=["synthetic HAXPES (design section 8, U9): comparison stays not_matched; decreasing input"])
need(not f5["exact"], "BE5 must round")
be_ok("BE6", "BE", "AX_D", "REC_FULL", "BD_249", "BE_B120", 3, ref=249.25, ci=False,
      public={"material": "Al", **PUBLIC}, settings={"epsrel": 1e-7}, support=(119, 119), state="truncated",
      checks=["real material; the KE image is the input of the M4b-1 row RB1 (D1000, E_B120)"])
be_ok("BE7", "BE", "AX_D", "REC_FULL", "BD_249", "BE_EM", 3, ref=250.25, choice="synthetic centroid, not the line",
      public={"material": "synthetic_wide", **PUBLIC}, support=(40, 40), state="truncated",
      warnings=["cumulative_loss_over_E0_gt_0p1"],
      checks=["the KE image is the M4b-1 input (D1000, E_M); electron_energy_eV = 1000 is given explicitly while "
              "C - B_ref = 999: electron_energy_minus_reference_eV = 1 exactly and the kernel is called with "
              "E0 = 1000 (the test wraps _call_kernel and reads the energy argument)"])
add("BE7", "kernel energy argument of every _call_kernel call", "diagnostic", "inputs")
need(ROWS[-1]["geometry"]["electron_energy_minus_reference_exact"] == "1/1", "BE7 difference")
f8 = be_ok("BE8", "BE", "AX_R", "REC_FULL", "BD_723", "BE_E8", 2, ref=72.3, kernel="K_A1s", support=(6, 6),
           state="truncated",
           checks=["C = Q(1486.6) - Q(4.6) is not a binary64 (RN(C) = 1482): rounding C first changes every one "
                   "of the six KE edges"])
need(f8["double_rounding_differs"] == 6 and Q(float(f8["C"])) != f8["C"], "BE8 must separate double rounding")
be_ok("BE9", "BE", "AX_PO", "REC_FULL", "BD_dup", "BE_E4", 2, ref=74.0, kernel="K_A1s", support=(3, 3),
      state="truncated",
      checks=["two deltas at the same BE and coordinates shared with edges are kept (merged only for the "
              "injectivity test); zero loss 1/2, 1/2, 0; kinetic_primary keeps both lines"])
fr = be_ok("BR1", "BE", "AX_8192", "REC_FULL", "BB_u", "BE_u", 1, ref=0.0, kernel="K_A1s", support=(0, 0),
           state="complete_on_display", choice="lower edge of the synthetic bin",
           checks=["a binned primary whose middle display edge rounds (u = 2^-40 eV): the returned zero loss is "
                   "that of the rounded-KE problem, 1/2 and 1/2; the exact BE overlap 3/8 and 5/8 is recorded "
                   "by the reference as a diagnostic and is not the expected value"])
need(fr["geometry"]["zero_loss_kinetic_widths"]["bin[0]"] == "1/2"
     and fr["geometry"]["zero_loss_exact_binding_overlap_diagnostic_only"]["bin[0]"] == "3/8", "BR1 masses")
be_ok("BZ0", "BE", "AX_D", "REC_FULL", "BD_249zero", "BE_EM", "display_support", ref=249.25,
      public={"material": "synthetic_wide", **PUBLIC}, support=(0, 0), state="complete_on_display",
      basis="geometry; the KE image succeeded on the current code (state only)",
      expect_extra={"orders": "empty tuple"},
      checks=["all intensities zero: the support is empty, both support orders are 0, display_support computes "
              "0 orders; region_loss_upper_eV is None for every bin; the KE control is called with display_support"])

# --- BT: truncation on the BE axis (both sides agree)
for rid_, edges_, sup_, state_ in (("BT1a", "BE_T1_eq", 3, "truncated"), ("BT1b", "BE_T1_dn", 2, "complete_on_display"),
                                   ("BT1c", "BE_T1_up", 3, "truncated")):
    f_ = be_ok(rid_, "BT", "AX_S", "REC_FULL", "BD_645", edges_, 2, ref=64.5, kernel="K_T1", support=(sup_, sup_),
               state=state_,
               checks=["high-BE end minus p = 3a (equality reaches) / one ulp below / one ulp above; C = 128 keeps "
                       "the one-ulp edge exact"])
    need(f_["exact"] and f_["geometry"]["truncation_equals_ke_control"], f"{rid_} must convert exactly")
be_ok("BDS", "BT", "AX_S", "REC_FULL", "BD_645", "BE_T1_eq", "display_support", ref=64.5, kernel="K_T1",
      support=(3, 3), state="complete_on_display",
      checks=["display_support on the BE axis computes the order with equality (3)"])
be_ok("BT2c", "BT", "AX_T2", "REC_FULL", "BD_645", "BE_T2", 4, ref=64.5, kernel="K_T2", support=(50, 50),
      state="truncated",
      expect_extra={"bound_detail": {"basis": "image of the M4b-1 row T2c",
                                     "finite_with_tau_gt_0": ["bin[2]", "bin[3]", "bin[4]", "bin[5]"],
                                     "largest": "bin[5]"}},
      checks=["the KE image is the M4b-1 row T2c; the bound of a bin uses its high-BE edge: bin[5] (94.5-114.5 eV) "
              "has the largest bound; bin[0], bin[1] are (0.0, 'complete', None)"])
need(ROWS[-1]["geometry"]["bound_class"] == ["complete", "complete", "reached", "reached", "reached", "reached"],
     "BT2c bound classes")
be_ok("BT2a", "BT", "AX_T2", "REC_FULL", "BD_645", "BE_T2", 4, ref=64.5, kernel="K_T2", support=(50, 50),
      state="bounded_on_display", settings={"neglect_tolerance_fraction": 0.05},
      basis="image of the M4b-1 row T2a (the bound against the tolerance is not proven here)",
      checks=["in-row: the three per-bin lists equal those of BT2c.result bit for bit; required_order_estimate "
              "equals that of the KE control run with the same tolerance"],
      reuse=["BT2c.result"])
for i_ in range(6):
    for name_ in ("neglected_at_or_above_region_bound", "bound_reason", "tau_per_eV"):
        add("BT2a", f"truncation.{name_}[bin[{i_}]] == BT2c.result", "bit_identity", "BT2c.result")
add("BT2a", "truncation.required_order_estimate", "truncation", "KE control with the same tolerance")
BT2B = {"energy_axis": axis_in("AX_T2", "REC_FULL", 64.5), "primary": "BD_645", "output_edges_eV": "BE_T2",
        "model": MODEL, "max_order": 4, **an("K_T2")}
row("BT2b", "BT", "analytic_spectrum", BT2B,
    {"status": "not_converged", "exception": "MultipleLossConvergenceError",
     "reason": "truncation_bound_exceeds_tolerance", "stage": "truncation", "required_order_estimate_at_least": 5,
     "basis": "image of the M4b-1 row T2b"},
    settings={"neglect_tolerance_fraction": 2.5e-3}, reuse=["BT2c.result"],
    checks=["tolerance between the bounds of bin[5] and bin[4]; the three per-bin lists of the exception are in "
            "the input BE order; R-C gives S_B = S_K = 50, so the KE control with the same tolerance must raise "
            "the same reason with the same estimate"])
for key_ in ("exception_type", "reason", "stage", "required_order_estimate >= 5"):
    add("BT2b", key_, "state", "expect")
add("BT2b", "required_order_estimate", "truncation", "KE control with the same tolerance")
add("BT2b", "the three per-bin lists of the exception (input BE order) == BT2c.result", "record_identity",
    "BT2c.result")
add("BT2b", "R-C support orders (50, 50) are the premise of the KE comparison", "truncation",
    "R-C (both-sided rule)")

# --- BS: the two sides disagree (design section 5, 6th edition)
BS_CHECK = "the KE control alone cannot give the truncation fields of this row (prereg 7.1)"
be_ok("BS1", "BS", "AX_8192", "REC_FULL", "BD_0", "BE_c5", 1, ref=0.0, kernel="K_T1e", support=(2, 1),
      state="truncated", reuse=[],
      expect_extra={"bound_detail": {"basis": "hand estimate about 1e-4 (0.16/tau^2 at the top of the tau grid)",
                                     "bin[0]": "finite, > 0, reason None, finite tau"}},
      checks=["design section 5: (S_B, S_K, S, N) = (2, 1, 2, 1); the rounded-KE run alone is complete with a "
              "zero bound", BS_CHECK])
be_ok("BS4", "BS", "AX_16", "REC_FULL", "BD_0", "BE_c5", 1, ref=0.0, kernel="K_T1e", support=(2, 2),
      state="truncated", reuse=["BS1.result"],
      checks=["the same kernel, primary and edges at C = 16: every coordinate converts exactly, so the KE control "
              "gives the truncation fields too", "in-row: the bound, reason and tau of bin[0] equal BS1.result bit "
              "for bit (the same rounded-up width 2 + 2^-43)"])
need(ROWS[-1]["geometry"]["truncation_equals_ke_control"]
     and ROWS[-1]["geometry"]["region_loss_upper_eV"] == ROWS[-2]["geometry"]["region_loss_upper_eV"], "BS1/BS4")
for name_ in ("neglected_at_or_above_region_bound", "bound_reason", "tau_per_eV"):
    add("BS4", f"truncation.{name_}[bin[0]] == BS1.result", "bit_identity", "BS1.result")
be_ok("BS1t", "BS", "AX_8192", "REC_FULL", "BD_0", "BE_c5", 1, ref=0.0, kernel="K_T1e", support=(2, 1),
      state="bounded_on_display", settings={"neglect_tolerance_fraction": 1e-2}, reuse=["BS4.result"],
      basis="hand estimate of the bound (about 1e-4; Codex computed 1.0095e-4 from the current KE code); the "
            "tolerance is about 99 times larger",
      expect_extra={"required_order_estimate": 1})
for name_ in ("neglected_at_or_above_region_bound", "bound_reason", "tau_per_eV"):
    add("BS1t", f"truncation.{name_}[bin[0]] == BS4.result", "bit_identity", "BS4.result")
add("BS1t", "truncation.required_order_estimate == 1", "truncation", "expect")
BS1_IN = {"energy_axis": axis_in("AX_8192", "REC_FULL", 0.0), "primary": "BD_0", "output_edges_eV": "BE_c5",
          "model": MODEL, **an("K_T1e")}
row("BS1x", "BS", "analytic_spectrum", {**BS1_IN, "max_order": 1},
    {"status": "not_converged", "exception": "MultipleLossConvergenceError",
     "reason": "truncation_bound_exceeds_tolerance", "stage": "truncation", "required_order_estimate": 2,
     "basis": "hand estimate of the bound (about 1e-4); the tolerance is 100 times smaller; at order 2 = S the "
              "display is complete"},
    settings={"neglect_tolerance_fraction": 1e-6}, reuse=["BS4.result"],
    checks=["the rounded-KE run alone would succeed (complete, bound 0)"])
for key_ in ("exception_type", "reason", "stage", "required_order_estimate == 2"):
    add("BS1x", key_, "state", "expect")
add("BS1x", "required_order_estimate == R-C support order S = 2", "truncation", "R-C (both-sided rule)")
add("BS1x", "the three per-bin lists of the exception == BS4.result", "record_identity", "BS4.result")
be_ok("BS1d", "BS", "AX_8192", "REC_FULL", "BD_0", "BE_c5", "display_support", ref=0.0, kernel="K_T1e",
      support=(2, 1), state="complete_on_display", control_order=2,
      checks=["display_support takes N = S = 2; the KE control is called with the integer 2", BS_CHECK])
row("BS1dl", "BS", "analytic_spectrum", {**BS1_IN, "max_order": "display_support"},
    {"status": "error", "exception": "MultipleLossInputError", "reason": "order_limit_exceeded",
     "requested_order": 2, "order_limit": 1, "kernel_evaluations": 0},
    settings={"order_limit": 1}, checks=["the order needed by the weaker side (2) exceeds order_limit = 1"])
for key_ in ("exception_type", "reason", "requested_order", "order_limit", "kernel_evaluations"):
    add("BS1dl", key_, "state", "expect")
add("BS1dl", "requested_order == R-C max(S_B, S_K) (the two orders are not read from the exception)",
    "truncation", "R-C (both-sided rule)")
be_ok("BS2", "BS", "AX_8192", "REC_FULL", "BD_0", "BE_c5rev", 1, ref=0.0, kernel="K_T1", support=(1, 2),
      state="truncated",
      expect_extra={"bound_detail": {"basis": "hand estimate", "bin[0]": "finite, > 0"}},
      checks=["the reverse counterexample (a = 1, top edge 2 - 2^-43): S_B = 1, S_K = 2"])
be_ok("BS3", "BS", "AX_8192", "REC_FULL", "BD_0", "BE_c5b", 1, ref=0.0, kernel="K_T1e", support=(8, 8),
      state="truncated", reuse=["BS4.result"],
      expect_extra={"bound_detail": {"basis": "hand estimate", "bin[0]": "finite, > 0, reason None, finite tau"}},
      checks=["both support orders are 8 and the state is truncated on either side; only bin[0] differs: the "
              "edge 2 + 2^-43 is above 2a in BE and below it in KE, so bin[0] is reached",
              "in-row: the bound of bin[0] equals that of BS4.result bit for bit (same width, kernel and order)"])
g3 = ROWS[-1]["geometry"]
need(g3["bound_class"] == ["reached", "reached"] and g3["support_comparison"]["next_reachability_mismatch_bins"] == [0]
     and g3["support_comparison"]["computed_relation_mismatch_count"] == 0, "BS3 geometry")
for name_ in ("neglected_at_or_above_region_bound", "bound_reason", "tau_per_eV"):
    add("BS3", f"truncation.{name_}[bin[0]] == BS4.result", "bit_identity", "BS4.result")
be_ok("BS3n2", "BS", "AX_8192", "REC_FULL", "BD_0", "BE_c5b", 2, ref=0.0, kernel="K_T1e", support=(8, 8),
      state="truncated",
      checks=["with N = 2 the same difference is a disagreement of a computed order (count 1, order 2, lower); "
              "bin[0] is complete on both sides"])
g3n = ROWS[-1]["geometry"]
need(g3n["bound_class"] == ["complete", "reached"]
     and g3n["support_comparison"]["computed_relation_mismatch_count"] == 1
     and g3n["support_comparison"]["first_mismatches"][0]["order"] == 2
     and g3n["support_comparison"]["first_mismatches"][0]["bound"] == "lower", "BS3n2 geometry")
be_ok("BS5", "BS", "AX_8192", "REC_FULL", "BD_0", "BE_c5q", 1, ref=0.0, kernel="K_T1e", support=(2, 2),
      state="truncated",
      checks=["top edge 2.25 + 2^-43: orders, state and every three-valued comparison agree, yet the width "
              "differs; the larger one (BE) is used even without a disagreement"])
g5 = ROWS[-1]["geometry"]
need(not g5["support_disagreement_warning"] and not g5["widths_equal_on_both_sides"], "BS5 geometry")
be_ok("BS6", "BS", "AX_16", "REC_FULL", "BD_0", "BE_c5q", 1, ref=0.0, kernel="K_T1e", support=(2, 2),
      state="truncated", reuse=["BS5.result"],
      checks=["exact control of BS5 at C = 16", "in-row: the bound, reason and tau of bin[0] equal BS5.result bit "
              "for bit; the KE control of BS5 (width 2.25) has a smaller bound"])
need(ROWS[-1]["geometry"]["truncation_equals_ke_control"]
     and ROWS[-1]["geometry"]["region_loss_upper_eV"] == g5["region_loss_upper_eV"], "BS5/BS6")
for name_ in ("neglected_at_or_above_region_bound", "bound_reason", "tau_per_eV"):
    add("BS6", f"truncation.{name_}[bin[0]] == BS5.result", "bit_identity", "BS5.result")
add("BS6", "bound of bin[0] > bound of the KE control of BS5", "inequality", "BS5 KE control")

# --- BW: the record
BW = dict(primary="BD_ends", edges="BE_E4", max_order=2, ref=74.0, kernel="K_A1s", support=(3, 3),
          state="truncated", reuse=["BE4.result"])
be_ok("BW1", "BW", "AX_PO", "REC_NR", **BW, warnings=["energy_axis_field_not_recorded"],
      expect_extra={"not_recorded": RECORD_FIELDS},
      checks=["all six record fields are 'not_recorded': success with one warning; a parent that is not "
              "recorded is listed once; numbers equal BE4.result"])
be_ok("BW2", "BW", "AX_PO", "REC_TILT", **BW, warnings=["recorded_emission_angle_differs_from_theta"],
      checks=["recorded emission angle 30 deg with theta_deg = 0: a warning, no correction"])
be_ok("BW3", "BW", "AX_PO", "REC_CORR", **BW,
      checks=["a nonzero calibration shift (-0.25 eV) and charge correction (1.5 eV) are recorded: the package "
              "does not shift the axis; every number and KE coordinate equals BE4.result"])
be_ok("BW4", "BW", "AX_PO", "REC_NEST", **BW, warnings=["energy_axis_field_not_recorded"],
      expect_extra={"not_recorded": ["source.monochromatic"]},
      checks=["one nested value is 'not_recorded': the list holds source.monochromatic only"])
for rid_ in ("BW1", "BW2", "BW3", "BW4"):
    add(rid_, "physical keys of prereg 7.1, truncation lists and KE coordinates of BE4.result (the record "
        "differs on purpose)", "record_identity", "BE4.result")
for rid_ in ("BW1", "BW4"):
    add(rid_, "energy_axis.not_recorded", "diagnostic", "expect")


# --- BN: rejections (kernel evaluations 0 in every row)
def reject(rid, entry, inputs, reason, *, exception="MultipleLossInputError", settings=None, checks=(),
           group="BN", **attrs):
    row(rid, group, entry, inputs,
        {"status": "error", "exception": exception, "reason": reason, "kernel_evaluations": 0, **attrs},
        settings=settings, checks=checks)
    for key in ("exception_type", "reason", "kernel_evaluations", *attrs):
        add(rid, key, "state", "expect")


def be_in(axis, rec, primary, edges, max_order, ref, model=MODEL, defects=(), **extra):
    return {"energy_axis": axis_in(axis, rec, ref, defects=defects), "primary": primary,
            "output_edges_eV": edges, "model": model, "max_order": max_order, **extra}


def op(path, value="<unset>", how="set"):
    out = {"op": how, "path": path}
    if how == "set":
        out["value"] = value
    return out


c7 = analyse("AX_8192", "BD_eps", "BE_c7", WINDOW["K_A1s"], 1)
need(c7["per_array_injective"] and c7["collision"]["binding_eV"]["values"] == [1.0, 1.0 + 2.0**-45]
     and c7["collision"]["kinetic_eV"]["values"] == [8191.0, 8191.0]
     and c7["collision"]["roles"] == ["output_edges_eV[1]", "primary.binding_energies_eV[0]"], "design section 7")
reject("BN01", "analytic_spectrum", be_in("AX_8192", "REC_FULL", "BD_eps", "BE_c7", 1, 1.0, **an("K_A1s")),
       "axis_conversion_not_injective", **c7["collision"],
       checks=["design section 7: each array alone is injective; the union is not"])
c70 = analyse("AX_8192", "BD_eps0", "BE_c7", WINDOW["K_A1s"], 1)
need(c70["collision"] == c7["collision"], "BN23")
reject("BN23", "analytic_spectrum", be_in("AX_8192", "REC_FULL", "BD_eps0", "BE_c7", 1, 1.0, **an("K_A1s")),
       "axis_conversion_not_injective", **c70["collision"],
       checks=["the colliding line has intensity 0: lines of zero intensity belong to the coordinate set"])
BASE_BN = dict(primary="BD_ends", edges="BE_E4", max_order=2, ref=74.0)
NAN = {"float": "nan"}
for rid_, axis_, rec_, defects_, field_, detail_, note_ in (
        ("BN04", "AX_PO", "REC_FULL", [op("reference.binding_anchor_eV", 84.0), op("reference.kinetic_anchor_eV", 1398.25)],
         "reference", "exclusive", "fields of both reference kinds"),
        ("BN05", "AX_PO", "REC_FULL", [op("reference.kinetic_energy_offset_eV", None)],
         "reference.kinetic_energy_offset_eV", "missing", "no value for the chosen kind"),
        ("BN06", "AX_PO", "REC_FULL", [op("reference.kind", "vacuum_level")], "reference.kind", "unknown_kind",
         "a vacuum-level reference is not supported"),
        ("BN07", "AX_PO", "REC_FULL", [op("photon_energy_eV", NR)], "photon_energy_eV", "not_recorded_not_allowed",
         "photon energy 'not_recorded'"),
        ("BN08", "AX_PO", "REC_NAMED", [op("photon_energy_eV", None)], "photon_energy_eV", "missing",
         "the source name says Al Kalpha; the number is never filled in from the name"),
        ("BN09", "AX_PA", "REC_FULL", [op("reference.kinetic_anchor_eV", NR)], "reference.kinetic_anchor_eV",
         "not_recorded_not_allowed", "anchor value 'not_recorded'"),
        ("BN12", "AX_PO", "REC_FULL", [{"op": "pass_plain_dict", "path": ""}], "energy_axis", "type",
         "a plain dict with the same fields, not the dataclass"),
        ("BN13", "AX_PO", "REC_FULL", [op("photon_energy_eV", NAN)], "photon_energy_eV", "not_finite", "nan"),
        ("BN14", "AX_PO", "REC_FULL", [op("source", 5)], "source", "type", "the record field is an integer"),
        ("BN24", "AX_PO", "REC_FULL", [op("photon_energy_eV", -1.0)], "photon_energy_eV", "out_of_range",
         "photon energy <= 0"),
        ("BN25", "AX_PO", "REC_FULL", [op("core_level", "  ")], "core_level", "empty",
         "a string of blanks")):
    reject(rid_, "analytic_spectrum", be_in(axis_, rec_, defects=defects_, **BASE_BN, **an("K_A1s")),
           "invalid_energy_axis", field=field_, detail=detail_, checks=[note_])
reject("BN10", "analytic_spectrum",
       {"energy_axis": axis_in("AX_PO", "REC_FULL", 74.0), "primary": {"kinetic_primary": "D1000"},
        "output_edges_eV": "BE_E4", "model": MODEL, "max_order": 2, **an("K_A1s")},
       "invalid_energy_axis", field="primary", detail="requires_binding_primary")
need(analyse("AX_S", "BD_645", "BE_neg_ke", WINDOW["K_A1s"], 2)["negative"], "BN11 needs a negative KE")
reject("BN11", "analytic_spectrum", be_in("AX_S", "REC_FULL", "BD_645", "BE_neg_ke", 2, 64.5, **an("K_A1s")),
       "invalid_energy_axis", field="output_edges_eV", detail="negative_kinetic_energy",
       checks=["C = 128 and the edge 130 eV: KE -2 eV"])
reject("BN15", "analytic_spectrum", be_in("AX_PO", "REC_FULL", "BD_ends", "BE_bad", 2, 74.0, **an("K_A1s")),
       "invalid_output_edges", exception="SingleLossInputError",
       checks=["BE edges that are neither increasing nor decreasing keep the M4a reason"])
reject("BN16", "api_spectrum",
       {"energy_axis": None, "primary": "BD_ends", "output_edges_eV": "BE_E4", "model": MODEL, "max_order": 2,
        "material": "synthetic_wide", **PUBLIC},
       "invalid_primary", exception="SingleLossInputError",
       checks=["a binding primary without energy_axis: the KE path is untouched and rejects it"])
PUB8192 = {"material": "synthetic_wide", **PUBLIC, "electron_energy_eV": 8192.0}
reject("BN17", "api_spectrum",
       {**be_in("AX_PO", "REC_FULL", defects=[op("reference.kind", "vacuum_level")], **BASE_BN),
        "material": "synthetic_wide", **PUBLIC, "loss_window_eV": [num(0.0), num(50.0)]},
       "invalid_energy_axis", field="reference.kind", detail="unknown_kind",
       checks=["order: the axis is checked before the M4a checks (the window [0, 50] is invalid too)"])
reject("BN18", "api_spectrum", {**be_in("AX_8192", "REC_FULL", "BD_eps", "BE_c7", 1, 1.0, model="b1"), **PUB8192},
       "axis_conversion_not_injective", **c7["collision"], checks=["order: injectivity before the model check"])
WIN_EPS = {"loss_window_eV": [num(A_EPS), num(50.0)]}
need(analyse("AX_8192", "BD_0", "BE_c5", (A_EPS, 50.0), 1)["geometry"]["support_order"] == 2, "BN19-BN21 window")
reject("BN19", "api_spectrum",
       {**be_in("AX_8192", "REC_FULL", "BD_0", "BE_c5", 1, 0.0, model="b1"), **PUB8192, **WIN_EPS},
       "invalid_model", checks=["the M4b-1 input checks are unchanged by a support disagreement"])
reject("BN20", "api_spectrum", {**be_in("AX_8192", "REC_FULL", "BD_0", "BE_c5", 4, 0.0), **PUB8192, **WIN_EPS},
       "work_lower_bound_exceeded", settings={"max_convolution_points": 4500},
       resource="max_convolution_points", lower_bound=4743, limit=4500,
       checks=["a support disagreement is not a rejection: the resource check of stage (i) answers"])
reject("BN21", "api_spectrum", {**be_in("AX_8192", "REC_FULL", "BD_0", "BE_c5", 65, 0.0), **PUB8192, **WIN_EPS},
       "order_limit_exceeded", requested_order=65, order_limit=64)
reject("BN22", "api_spectrum", {**be_in("AX_8192", "REC_FULL", "BD_0", "BE_c5b", 65, 0.0), **PUB8192, **WIN_EPS},
       "order_limit_exceeded", requested_order=65, order_limit=64,
       checks=["no loop over orders runs before the limit check"])

# --- BC: CLI
COMMON = ["--loss-window-eV", "1", "50", "--weight-length-A", "20", "--epsilon-momentum", "k_total"]
ARGV = [*COMMON, "--model", MODEL, "--max-order", "3"]
ARITY = {"--loss-window-eV": 2, "--weight-length-A": 1, "--epsilon-momentum": 1, "--model": 1, "--max-order": 1,
         "--neglect-tolerance-fraction": 1, "--order-limit": 1, "--max-convolution-points": 1}
REQUIRED = {"xps-multiple-loss": ["--loss-window-eV", "--weight-length-A", "--epsilon-momentum", "--model",
                                  "--max-order"],
            "xps-single-loss": ["--loss-window-eV", "--weight-length-A", "--epsilon-momentum"],
            "capabilities": []}
FIXED_VALUES = {"--loss-window-eV": {("1", "50"), ("1", "8"), ("1.0000000000000284", "50")},
                "--weight-length-A": {("20",)}, "--epsilon-momentum": {("k_total",)}, "--model": {(MODEL,)},
                "--max-order": {("1",), ("2",), ("3",)}}
CONFIG = {"config_material": "synthetic_wide", "config_electron_energy_eV": 1000.0,
          "config_integration": {"epsabs": 1e-12, "epsrel": 1e-9, "limit": 200}}


def check_argv(rid, command, argv):
    seen, i = [], 0
    while i < len(argv):
        need(argv[i] in ARITY, f"{rid}: unknown option {argv[i]}")
        values = argv[i + 1:i + 1 + ARITY[argv[i]]]
        need(len(values) == ARITY[argv[i]] and not any(v.startswith("--") for v in values),
             f"{rid}: option {argv[i]} lacks its value")
        seen.append(argv[i])
        i += 1 + ARITY[argv[i]]
    need(sorted(seen) == sorted(set(seen)) and set(REQUIRED[command]) <= set(seen),
         f"{rid}: required options of {command}")
    for option, allowed in FIXED_VALUES.items():
        if option in argv:
            k = argv.index(option)
            need(tuple(argv[k + 1:k + 1 + ARITY[option]]) in allowed, f"{rid}: value of {option}")
    if command == "xps-single-loss":
        need(not {"--model", "--max-order"} & set(seen), f"{rid}: options of the other command")


def measured(name, rec, primary, edges, ref, defects=(), choice=DEFAULT_CHOICE):
    spec = {"schema": SCHEMA_BE, "axis": "binding_energy_eV",
            "measurement": axis_in(name, rec, ref, choice), "primary": primary, "output_binding_edges_eV": edges}
    if defects:
        spec["defects"] = list(defects)
    return spec


def cli(rid, expect, *, argv=ARGV, spectrum, entry="cli", ci=True, checks=(), command="xps-multiple-loss",
        group="BC", geometry=None):
    check_argv(rid, command, list(argv))
    inputs = {"command": command, "argv": list(argv), **CONFIG, "spectrum": spectrum}
    row(rid, group, entry, inputs, expect, ci=ci, checks=checks, denominator="replay_and_state", geometry=geometry)
    add(rid, "exit_code", "state", "expect")
    if expect["status"] == "error":
        for key in ("error_type", "reason", "output_dir_after == absent", "stderr is one JSON object",
                    *[k for k in expect if k in ("field", "detail", "binding_eV", "kinetic_eV", "roles")]):
            add(rid, key, "state", "expect")
        if expect["reason"] == "invalid_schema" and command == "xps-multiple-loss":
            add(rid, "detail attribute is a non-empty ASCII string", "state", "expect")


def cli_ok(rid, axis, rec, primary, edges, order, ref, *, support, state, argv, checks=(), warnings=(),
           window=(1.0, 50.0), ci=True):
    facts = analyse(axis, primary, edges, window, order, ref=ref, e0=1000.0)
    g = facts["geometry"]
    need(facts["collision"] is None and (g["support_order_binding"], g["support_order_kinetic"]) == tuple(support)
         and g["state_without_tolerance"] == state, f"{rid}: geometry")
    codes = ["binding_axis_assumed_calibrated", *warnings]
    if g["support_disagreement_warning"]:
        codes.append("axis_conversion_support_disagreement")
    cli(rid, {"status": "complete", "exit_code": 0, "truncation_state": state, "support_order": g["support_order"],
              "basis": "geometry; the KE image succeeded on the current code (state only)",
              "warnings_include": codes},
        argv=argv, spectrum=measured(axis, rec, primary, edges, ref), checks=checks, geometry=g, ci=ci)
    nb = n_bins(edges)
    add(rid, "replay: API result rebuilt from run.json alone equals the CLI results", "record_identity", "replay")
    for name in ("multiple_loss_spectrum.csv", "order_areas.csv"):
        add(rid, f"replay: every cell of {name} through an independent column map", "record_identity", "replay")
    for name in ("header of multiple_loss_spectrum.csv", "row order of multiple_loss_spectrum.csv",
                 "the four edge columns (blank cells on the infinite sides)", "files == the three names"):
        add(rid, name, "diagnostic", "prereg 4.6")
    for name in ("energy_axis == results.energy_axis", "measurement == the file", "spectrum_input.sha256",
                 "spectrum_input.schema", "spectrum_input.binding_output_edges_sha256_float64_le",
                 "spectrum_input.kinetic_output_edges_sha256_float64_le", "spectrum_definition.axis",
                 "units.binding_energy", "model_id", "primary and output_edges_eV are the BE input record",
                 "results.single_loss.primary and edges are KE", "warnings", "truncation",
                 "effective_settings has the 13 names"):
        add(rid, f"run.json: {name}", "diagnostic", "prereg 4.6")
    for k in range(nb + 1):
        add(rid, f"run.json: energy_axis.kinetic_output_edges_eV[{k}]", "exact_conversion", "R-C")
    for name in ("conversion_constant_eV", "conversion_constant_exact", "implied_offset_eV", "rounding_max_eV",
                 "electron_energy_minus_reference_eV", "support_order_binding", "support_order_kinetic"):
        add(rid, f"run.json: energy_axis.{name}", "exact_conversion", "R-C")
    for region in regions(nb):
        add(rid, f"CSV zero_loss.{region}", "bin_membership", "R-C (KE widths)")
    return g


ROW_ORDER = "rows of the CSV: low_binding, 0.., high_binding for increasing input; reversed for decreasing input"
cli_ok("BC01s", "AX_D", "REC_FULL", "BD_249", "BE_EMs", 2, 249.25, support=(5, 5), state="truncated",
       argv=["--loss-window-eV", "1", "8", *COMMON[3:], "--model", MODEL, "--max-order", "2"], window=(1.0, 8.0),
       checks=["the small CI form of BC01 (window [1, 8] eV, N = 2, three bins): the same replay and CSV checks; "
               "no cumulative-loss warning (16/1000)", ROW_ORDER])
ROWS[-1]["expect"]["warnings_exclude"] = ["cumulative_loss_over_E0_gt_0p1"]
cli_ok("BC01", "AX_D", "REC_FULL", "BD_249", "BE_EM", 3, 249.25, support=(40, 40), state="truncated", argv=ARGV,
       warnings=["cumulative_loss_over_E0_gt_0p1"], ci=False,
       checks=["replay: rebuild the axis and the binding primary from run.json alone (prereg 4.6), call the API, "
               "compare values, errors, state, warnings and every numeric CSV cell bit for bit", ROW_ORDER])
cli_ok("BC02", "AX_DA", "REC_NR", "BD_249", "BE_EMr", 1, 249.25, support=(40, 40), state="truncated",
       argv=[*COMMON, "--model", MODEL, "--max-order", "1"], warnings=["energy_axis_field_not_recorded"],
       ci=False,
       checks=["paired anchor and 'not_recorded' records through the CLI; decreasing input: rows high_binding, "
               "0..6, low_binding; implied offset 4.25 exactly", ROW_ORDER])
add("BC02", "run.json: energy_axis.not_recorded", "diagnostic", "expect")
gc7 = cli_ok("BC07", "AX_8192", "REC_FULL", "BD_0", "BE_c5", 1, 0.0, support=(2, 1), state="truncated",
             argv=["--loss-window-eV", "1.0000000000000284", "50", *COMMON[3:], "--model", MODEL, "--max-order", "1"],
             window=(A_EPS, 50.0),
             checks=["design section 5 through the CLI: complete with exit 0, truncated, a finite positive bound of "
                     "bin[0] and the disagreement warning in the API replay and in run.json; E0 = 1000 stays in CONFIG",
                     "1.0000000000000284 parses to 1 + 2^-45"])
need(gc7["bound_class"] == ["reached"] and gc7["support_disagreement_warning"], "BC07 geometry")
for where_ in ("run.json", "API replay"):
    for name_ in ("truncation three fields of bin[0]: finite and > 0, reason None, finite tau",
                  "warning axis_conversion_support_disagreement"):
        add("BC07", f"{where_}: {name_}", "truncation", "expect")
cli_ok("BC12", "AX_D", "REC_FULL", "BD_m1", "BE_grid", 1, -1.0, support=(3, 3), state="truncated",
       argv=[*COMMON, "--model", MODEL, "--max-order", "1"], warnings=["E0_outside_primary_support"], ci=False,
       checks=["output_binding_edges_eV is the compact form {min: -2, max: 2, step: 0.5} with a negative min; the "
               "expanded edges equal the registered binary64 values; a delta at BE -1 eV"])
ROWS[-1]["inputs"]["spectrum"]["output_binding_edges_eV"] = {"compact_form_of": "BE_grid"}
CLI_ERR = {"status": "error", "exit_code": 2, "output_dir_after": "absent", "stderr": "one JSON object"}
GOOD = dict(name="AX_D", rec="REC_NAMED", primary="BD_249", edges="BE_EM", ref=249.25)
cli("BC03", {**CLI_ERR, "error_type": "MultipleLossInputError", "reason": "invalid_schema"},
    spectrum=measured(**GOOD, defects=[op("measurement.photon_energy_eV", how="delete_key")]),
    checks=["a missing key is a structure error; the source name says Al Kalpha and nothing is filled in"])
cli("BC04", {**CLI_ERR, "error_type": "MultipleLossInputError", "reason": "invalid_energy_axis",
             "field": "photon_energy_eV", "detail": "not_recorded_not_allowed"},
    spectrum=measured(**GOOD, defects=[op("measurement.photon_energy_eV", NR)]))
cli("BC05", {**CLI_ERR, "error_type": "MultipleLossInputError", "reason": "invalid_schema"},
    spectrum=measured(**GOOD, defects=[op("observed_counts", [1, 2, 3])]),
    checks=["an unknown top-level key: observed data are not part of this schema"])
cli("BC06", {**CLI_ERR, "error_type": "MultipleLossInputError", "reason": "axis_conversion_not_injective",
             **c7["collision"]},
    argv=[*COMMON, "--model", MODEL, "--max-order", "1"],
    spectrum=measured("AX_8192", "REC_FULL", "BD_eps", "BE_c7", 1.0),
    checks=["design section 7 through the CLI; the JSON numbers read back to the same binary64"])
add("BC06", "kernel_evaluations == 0", "state", "expect")
cli("BC08", {**CLI_ERR, "error_type": "SingleLossInputError", "reason": "invalid_schema"},
    argv=COMMON, command="xps-single-loss", spectrum=measured(**GOOD),
    checks=["the measured-spectrum schema is accepted by xps-multiple-loss only; the arguments are complete, so "
            "the rejection comes from reading the JSON"])
cli("BC09", {**CLI_ERR, "error_type": "MultipleLossInputError", "reason": "invalid_schema"},
    spectrum=measured(**GOOD, defects=[op("axis", "kinetic_energy_eV")]),
    checks=["the measured-spectrum schema has the axis binding_energy_eV only"])
cli("BC10", {**CLI_ERR, "error_type": "MultipleLossInputError", "reason": "invalid_energy_axis",
             "field": "reference", "detail": "exclusive"},
    spectrum=measured(**GOOD, defects=[op("measurement.reference.binding_anchor_eV", None)]),
    checks=["a key of the other reference kind is present (its value is null): exclusive, not a structure error"])
cli("BC11", {"status": "complete", "exit_code": 0}, argv=[], spectrum=None, entry="cli_capabilities",
    command="capabilities")
for key_ in ("multiple_scattering.unsupported_physics == the nine names in order (binding_energy_axis removed)",
             "multiple_scattering.energy_axes", "multiple_scattering.measured_spectrum_schema",
             "baseline capabilities with one removal and two additions == current (not_run without git history)",
             "supported == baseline", "each side reports its own version"):
    add("BC11", key_, "diagnostic", "prereg 4.6")


# --- BASE: energy_axis=None and the existing entries against the baseline code
def base(rid, entry, inputs, *, ci=True, checks=(), settings=None, command=None, expect=None):
    if entry == "cli":
        check_argv(rid, command, inputs["argv"])
        inputs = {"command": command, **inputs}
    expect = expect or ({"status": "complete", "exit_code": 0} if entry == "cli" else {"status": "success"})
    row(rid, "BASE", entry, inputs, {**expect, "without_git_history": "not_run/sdist_without_git_history"},
        ci=ci, settings=settings, denominator="bit_identity",
        checks=["baseline child process (git archive of the full baseline SHA): names, types, dictionary keys, "
                "warning order, array dtype, shape and bytes, scalar bits", *checks])
    if entry == "cli":
        add(rid, "exit_code == 0 on both sides", "state", "expect")
        for name in ("multiple_loss_spectrum.csv or single_loss_spectrum.csv bytes", "second CSV bytes",
                     "run.json without time, elapsed time, paths and version (no energy_axis key anywhere)"):
            add(rid, name, "record_identity", "baseline 8c33d29")
        for name in ("created_utc is valid", "elapsed_seconds is valid", "paths are valid",
                     "each side reports its own version"):
            add(rid, name, "diagnostic", "baseline 8c33d29")
    else:
        add(rid, "every field of the result or exception", "record_identity", "baseline 8c33d29")


KE_PUB = {"material": "synthetic_wide", **PUBLIC, "primary": "D1000", "output_edges_eV": "E_M", "model": MODEL,
          "max_order": 3}
KE_PUB_S = {**KE_PUB, "loss_window_eV": [num(1.0), num(8.0)], "output_edges_eV": "E_Ms", "max_order": 2}
for rid_, inputs_, ci_ in (("BASE-ms-s", KE_PUB_S, True), ("BASE-ms", KE_PUB, False)):
    base(rid_, "api_spectrum", inputs_, ci=ci_,
         checks=["in-row: the call with an explicit energy_axis=None returns the same bytes",
                 "the only field the baseline lacks is energy_axis directly under the spectrum result, and it is "
                 "None", "BASE-ms-s is the small CI form (window [1, 8] eV, N = 2, three bins)"])
    add(rid_, "explicit energy_axis=None equals the omitted call", "record_identity", "same row")
    add(rid_, "result.energy_axis is None", "diagnostic", "expect")
base("BASE-ms-s2", "api_spectrum", {**KE_PUB, "multiplier": "inverse_window_area"}, ci=False)
base("BASE-areas", "api_areas", {k: v for k, v in KE_PUB.items() if k not in ("primary", "output_edges_eV")},
     ci=False, checks=["in-row: multiple_loss_order_areas(..., energy_axis=...) raises TypeError (no such parameter)"])
add("BASE-areas", "energy_axis keyword raises TypeError", "diagnostic", "expect")
base("BASE-analytic", "analytic_spectrum",
     {"primary": "D1000", "output_edges_eV": "E_A3", "model": MODEL, "max_order": 3, **an("K_A3")},
     checks=["the second element (panels_by_order) is compared too"])
base("BASE-single", "api_single_spectrum", {k: v for k, v in KE_PUB.items() if k not in ("model", "max_order")},
     ci=False)
base("BASE-area1", "api_single_area",
     {k: v for k, v in KE_PUB.items() if k not in ("model", "max_order", "primary", "output_edges_eV")}, ci=False)
base("BASE-n13", "api_spectrum", {**KE_PUB, "max_order": 4}, settings={"max_convolution_points": 4500},
     expect={"status": "error", "exception": "MultipleLossInputError", "reason": "work_lower_bound_exceeded",
             "resource": "max_convolution_points", "lower_bound": 4743, "limit": 4500, "kernel_evaluations": 0},
     checks=["the M4b-1 row N13 against the baseline: type, reason, vars(exc) with every old attribute, message; "
             "omitted and explicit energy_axis=None; the new attributes exist only for the BE reasons"])
for key_ in ("exception_type", "reason", "resource", "lower_bound", "limit", "kernel_evaluations"):
    add("BASE-n13", key_, "state", "expect")
add("BASE-n13", "explicit energy_axis=None raises the same exception (vars and message)", "record_identity",
    "same row")
KE_SPEC = {"schema": SCHEMA_KE, "axis": "kinetic_energy_eV", "primary": "D1000", "output_edges_eV": "E_M"}
base("BASE-cli", "cli", {"argv": ARGV, **CONFIG, "spectrum": KE_SPEC}, command="xps-multiple-loss", ci=False)
base("BASE-cli-single", "cli", {"argv": COMMON, **CONFIG, "spectrum": KE_SPEC}, command="xps-single-loss", ci=False)

# ---------------------------------------------------------------- mutations (14 required)
MUTATIONS = []


def mutation(mid, origin, text, branches):
    """branches: list of (name, [(row, key, kind, change)], predicted) ; predicted names a terminal change."""
    out = {"id": mid, "class": "required", "origin": origin, "mutation": text, "branches": []}
    for name, records, predicted in branches:
        for r, key, kind, _ in records:
            need(any(k == key and kd == kind for k, kd, _o in LEDGER.get(r, [])),
                 f"{mid}: ({r}, {key}, {kind}) is not in the ledger")
        out["branches"].append({"branch": name, "predicted_terminal_change": predicted,
                                "detected_by": [{"row": r, "key": k, "kind": kd, "change": ch}
                                                for r, k, kd, ch in records]})
    MUTATIONS.append(out)


EC, BM, BI, TR, ST, DG = "exact_conversion", "bin_membership", "bit_identity", "truncation", "state", "diagnostic"
mutation("X25", "design", "drop the sign of BE or drop phi", [
    ("E = B - C", [("BE1", "terminal_state", ST, "success becomes invalid_energy_axis")],
     "MultipleLossInputError invalid_energy_axis, detail negative_kinetic_energy"),
    ("E = C + B", [("BE1", "energy_axis.kinetic_output_edges_eV[0]", EC, "KE edge differs"),
                   ("BE1", "zero_loss.outside.low_binding", BM, "1/2 becomes 0"),
                   ("BE1", "order[1].bin[4].total", BI, "the loss moves to lower BE")],
     "or axis_conversion_not_injective if the strict decrease of section 2.2 is tested"),
    ("C = h nu (phi dropped)", [("BE7", "energy_axis.conversion_constant_eV", EC, "1253.5 instead of 1249.25"),
                                ("BE7", "energy_axis.kinetic_output_edges_eV[0]", EC, "KE edge differs"),
                                ("BE7", "energy_axis.electron_energy_minus_reference_eV", EC, "1 becomes -3.25")],
     None)])
mutation("X26", "design", "do not restore the input order / do not swap the outside names", [
    ("arrays not reordered", [("BE1", "zero_loss.bin[0]", BM, "1 becomes 0"),
                              ("BE1", "order[1].bin[4].total", BI, "bin 4 holds the value of bin 2 (no loss there)"),
                              ("BE1r", "physical keys of prereg 7.1 of BE1.result with BE display bin[i] -> bin[6 - i]; KE coordinates in "
                                       "the same ascending order",
                               "record_identity", "a decreasing input must not be reversed twice")], None),
    ("outside.low_binding := internal outside.low",
     [("BE1", "zero_loss.outside.low_binding", BM, "1/2 becomes 0"),
      ("BE1", "order[3].outside.high_binding.total", BI, "the third order moves to the other side")], None)])
mutation("X27", "design", "a delta on an interior BE edge goes to the higher-BE bin", [
    ("", [("BE4", "zero_loss.bin[1]", BM, "1/2 becomes 0"), ("BE4", "zero_loss.bin[2]", BM, "1/4 becomes 3/4"),
          ("BE1", "zero_loss.bin[0]", BM, "1 becomes 0")], None)])
mutation("X32", "design", "judge the support from the rounded KE differences only", [
    ("", [("BS1", "truncation.support_order", TR, "2 becomes 1"),
          ("BS1", "truncation.state", TR, "truncated becomes complete_on_display"),
          ("BS1", "truncation.bound_class[bin[0]]", TR, "reached becomes complete"),
          ("BS1", "truncation three fields of bin[0]: finite and > 0, reason None, finite tau", TR,
           "(0.0, 'complete', None)"),
          ("BS1", "energy_axis.support_order_binding", TR, "2 is missing or 1")],
     "both versions succeed; the terminal state alone does not detect it")])
mutation("X34", "added", "round twice: RN(RN(C) - B)", [
    ("", [("BE8", f"energy_axis.kinetic_output_edges_eV[{k}]", EC, "differs") for k in range(6)]
     + [("BE8", "energy_axis.rounding_max_eV", EC, "differs")], None)])
mutation("X35", "added", "test injectivity per array only", [
    ("only the union test is weakened",
     [("BN01", "terminal_state", ST, "the registered rejection disappears: success"),
      ("BN01", "roles", ST, "no colliding pair"),
      ("BC06", "terminal_state", ST, "complete, exit 0")],
     "BN01 succeeds and BC06 completes (the rounded-KE images succeed); another exception is not a detection")])
mutation("X36", "added", "paired anchor with the wrong sign (C = E* - B*)", [
    ("", [("BE2", "energy_axis.conversion_constant_eV", EC, "1314.25 instead of 1482.25"),
          ("BE2", "energy_axis.kinetic_output_edges_eV[0]", EC, "KE edge differs")], None)])
mutation("X37", "added", "fill a missing photon energy with 1486.6 eV from the source name (applied before the "
         "schema check)", [
    ("API: explicit None", [("BN08", "terminal_state", ST, "the rejection disappears"),
                            ("BN08", "detail", ST, "no 'missing'")], "success"),
    ("CLI: key absent", [("BC03", "terminal_state", ST, "the rejection disappears"),
                         ("BC03", "reason", ST, "no invalid_schema")], "complete")])
mutation("X38", "added", "accept 'not_recorded' silently", [
    ("no warning (the list of names may stay)",
     [("BW1", "warnings (BE identifiers, their order, source and message pairing)", DG,
       "energy_axis_field_not_recorded is missing")], None),
    ("photon energy 'not_recorded' replaced by the finite value 1486.6",
     [("BN07", "terminal_state", ST, "success"), ("BN07", "detail", ST, "no not_recorded_not_allowed")],
     "success; an accidental TypeError is not a detection")])
mutation("X39", "added", "reorder the values but not the per-bin truncation lists", [
    ("bounds", [("BT2c", "truncation.neglected_at_or_above_region_bound[bin[5]]", BI, "largest bound moves")], None),
    ("reasons", [("BT2c", "truncation.bound_reason[bin[5]]", BI, "None becomes 'complete'")], None),
    ("tau", [("BT2c", "truncation.tau_per_eV[bin[5]]", BI, "tau moves")], None)])
mutation("X40", "added", "take the support order from both sides but the per-bin widths from KE only", [
    ("per-bin bound", [("BS3", "truncation.bound_class[bin[0]]", TR, "reached becomes complete"),
                       ("BS3", "truncation three fields of bin[0]: finite and > 0, reason None, finite tau", TR,
                        "(0.0, 'complete', None)"),
                       ("BS3", "energy_axis.region_loss_upper_eV[bin[0]]", TR, "2 + 2^-43 becomes 2")],
     "both versions return truncated"),
    ("diagnosis of computed orders omitted",
     [("BS3n2", "energy_axis.support_comparison.computed_relation_mismatch_count", TR, "1 becomes 0"),
      ("BS3n2", "energy_axis.support_comparison.first_mismatches", "record_identity", "empty")], None),
    ("larger width only when a comparison disagrees",
     [("BS6", "truncation.neglected_at_or_above_region_bound[bin[0]] == BS5.result", BI,
       "BS5 falls back to the KE width 2.25")], None)])
mutation("X41", "added", "apply the implied offset of a paired anchor (C = h nu)", [
    ("", [("BE2", "energy_axis.conversion_constant_eV", EC, "Q(1486.6)"),
          ("BE2", "energy_axis.kinetic_output_edges_eV[0]", EC, "KE edge differs"),
          ("BC02", "run.json: energy_axis.implied_offset_eV", EC, "4.25 with wrong KE edges, or 0"),
          ("BC02", "run.json: energy_axis.kinetic_output_edges_eV[0]", EC, "KE edge differs")], None)])
mutation("X42", "added", "take E0 from C - B_ref instead of electron_energy_eV", [
    ("", [("BE7", "order[1].bin[5].total", BI, "the kernel is evaluated at 999 eV"),
          ("BE7", "kernel energy argument of every _call_kernel call", DG, "999 instead of 1000")], None)])
mutation("X43", "added", "put a delta on a global end of the output into the outside region", [
    ("low end", [("BE4", "zero_loss.bin[0]", BM, "1/4 becomes 0"),
                 ("BE4", "zero_loss.outside.low_binding", BM, "0 becomes 1/4")], None),
    ("high end", [("BE4", "zero_loss.bin[2]", BM, "1/4 becomes 0"),
                  ("BE4", "zero_loss.outside.high_binding", BM, "0 becomes 1/4")], None)])

# ---------------------------------------------------------------- checks of the table
ENTRIES = {"analytic_spectrum", "api_spectrum", "api_areas", "api_single_spectrum", "api_single_area", "cli",
           "cli_capabilities"}
GROUPS = ["BE", "BT", "BS", "BW", "BN", "BC", "BASE"]
KINDS = {"state", "exact_conversion", "bin_membership", "bit_identity", "record_identity", "truncation",
         "bound_bracket", "inequality", "diagnostic"}
RECORD_KINDS = {"record_identity"}
NEW_REASONS = {"invalid_energy_axis", "axis_conversion_not_injective", "invalid_schema"}
ids = [r["id"] for r in ROWS]
need(len(ids) == len(set(ids)) and set(ids) == set(LEDGER), "row ids and ledger")
by_id = {r["id"]: r for r in ROWS}


def alias(x):
    if isinstance(x, dict):
        return {k: alias(v) for k, v in x.items()}
    if isinstance(x, list):
        return [alias(v) for v in x]
    if isinstance(x, str):
        for table in (BE_PRIMARIES, BE_EDGES, KE_PRIMARIES, KE_EDGES):
            if x in table:
                return {k: v for k, v in table[x].items() if k not in ("note", "compact")}
    return x


def build_axis(spec):
    """The final axis object: axis numbers, the six record fields and the two reference-position fields."""
    out = copy.deepcopy(AXES[spec["axis"]])
    out.update(copy.deepcopy(RECORDS[spec["record"]]))
    out["primary_reference_binding_eV"] = spec["primary_reference_binding_eV"]
    out["primary_reference_choice"] = spec["primary_reference_choice"]
    return out


def apply(obj, ops):
    for o in ops:
        if o["op"] == "pass_plain_dict":
            obj = {"plain_dict": obj}
            continue
        *parents, last = o["path"].split(".")
        node = obj
        for name in parents:
            node = node[name]
        if o["op"] == "delete_key":
            del node[last]
        else:
            node[last] = o["value"]     # an explicit null is kept as a key and differs from a deleted key
    return obj


def resolved(r):
    """Signature built from the final arguments of the entry: aliases, defects, kernel and defaults expanded."""
    final = copy.deepcopy(r["inputs"])
    axis = final.get("energy_axis")
    if isinstance(axis, dict):
        final["energy_axis"] = apply(build_axis(axis), axis.get("defects", ()))
    if "kernel" in final:
        final["kernel"] = {"name": final["kernel"], "window": WINDOW[final["kernel"]]}
    spectrum = final.get("spectrum")
    if isinstance(spectrum, dict) and "measurement" in spectrum:
        body = {k: v for k, v in spectrum.items() if k != "defects"}
        body["measurement"] = build_axis(spectrum["measurement"])
        final["spectrum"] = apply(body, spectrum.get("defects", ()))
    defaults = dict(OLD["defaults"]["first_order_settings"], **OLD["defaults"]["order_settings"])
    return json.dumps(binary64([r["entry"], alias(final), {**defaults, **r["settings"]}]), sort_keys=True)


def binary64(x):
    """Every number as the final binary64: a {value, hex} pair, a {values, hex} list and a bare number agree."""
    if isinstance(x, dict):
        if set(x) == {"value", "hex"}:
            return {"binary64": float(x["value"]).hex()}
        if {"values", "hex"} <= set(x):
            rest = {k: binary64(v) for k, v in x.items() if k not in ("values", "hex")}
            return {**rest, "binary64": [float(v).hex() for v in x["values"]]}
        return {k: binary64(v) for k, v in x.items()}
    if isinstance(x, list):
        return [binary64(v) for v in x]
    if isinstance(x, (int, float)) and not isinstance(x, bool):
        return {"binary64": float(x).hex()}
    return x


def _probe(defects):
    spec = {"axis": "AX_PO", "record": "REC_FULL", "primary_reference_binding_eV": 1,
            "primary_reference_choice": "c"}
    if defects:
        spec["defects"] = defects
    return resolved({"entry": "x", "settings": {}, "inputs": {"energy_axis": spec}})


need(_probe([{"op": "set", "path": "core_level", "value": "synthetic level"}]) == _probe([]),
     "a set of the same value must not change the signature")
need(_probe([{"op": "set", "path": "photon_energy_eV", "value": 1486.6}]) == _probe([])
     and _probe([{"op": "set", "path": "photon_energy_eV", "value": 1486.7}]) != _probe([]),
     "a numeric set of the same binary64 must not change the signature, another number must")
need(_probe([{"op": "set", "path": "core_level", "value": None}]) != _probe([{"op": "delete_key", "path": "core_level"}]),
     "an explicit null and a deleted key must differ")
signatures = [resolved(r) for r in ROWS]
need(len(signatures) == len(set(signatures)), "two rows share entry, inputs and settings after expansion")
for r in ROWS:
    need(r["entry"] in ENTRIES and r["group"] in GROUPS, f"{r['id']}: vocabulary")
    need(all(kind in KINDS for _k, kind, _o in LEDGER[r["id"]]), f"{r['id']}: check kind")
    for target in r["reuse"]:
        src, _, artifact = target.partition(".")
        need(src in by_id and artifact in by_id[src]["provides"], f"{r['id']}: reuse {target}")
        need(ids.index(src) < ids.index(r["id"]), f"{r['id']}: reuse must point to an earlier row")
        need(by_id[src]["ci"] or not r["ci"], f"{r['id']}: a CI row reuses a full-only row")
    for _key, _kind, oracle in LEDGER[r["id"]]:
        src = oracle.partition(".")[0]
        if oracle.endswith(".result"):
            need(oracle in r["reuse"], f"{r['id']}: oracle {oracle} is not in reuse")
    groups = {}
    for key, kind, oracle in LEDGER[r["id"]]:
        groups.setdefault((kind, oracle), []).append(key)
    r["check_plan"] = [{"kind": kind, "oracle": oracle, "count": len(keys), "examples": keys[:3]}
                       for (kind, oracle), keys in groups.items()]
need(NEW_REASONS <= {r["expect"].get("reason") for r in ROWS}, "a new reason has no row")
need("axis_conversion_changes_support" not in json.dumps(ROWS), "the removed reason is still used")
for m in MUTATIONS:
    for branch in m["branches"]:
        for d in branch["detected_by"]:
            if m["id"] not in by_id[d["row"]]["mutations"]:
                by_id[d["row"]]["mutations"].append(m["id"])
need([m["id"] for m in MUTATIONS] == ["X25", "X26", "X27", "X32", *[f"X{k}" for k in range(34, 44)]], "mutations")
CI_ORDER = ["BZ0", "BE7", "BC01s"]
ci_counts, moved = [], set()
for step in [None, *CI_ORDER]:
    if step:
        need(by_id[step]["ci"], f"{step} is not a CI row")
        moved.add(step)
    ci_now = {r["id"] for r in ROWS if r["ci"]} - moved
    for r in ROWS:
        need(r["id"] not in ci_now or all(t.partition(".")[0] in ci_now for t in r["reuse"]),
             f"after moving {sorted(moved)}: {r['id']} reuses a full-only row")
    ci_counts.append({"ci": len(ci_now), "full_only": len(ROWS) - len(ci_now)})
states, kinds, by_group = {}, {}, {}
for r in ROWS:
    states[r["expect"]["status"]] = states.get(r["expect"]["status"], 0) + 1
    by_group[r["group"]] = by_group.get(r["group"], 0) + 1
    for _key, kind, _oracle in LEDGER[r["id"]]:
        kinds[kind] = kinds.get(kind, 0) + 1
need((len(ROWS), sum(r["ci"] for r in ROWS), sum(not r["ci"] for r in ROWS)) == (N_ROWS, N_CI, N_FULL),
     f"row counts {len(ROWS)}, {sum(r['ci'] for r in ROWS)}")
need(states == STATES, f"terminal states {states}")
scalar_keys = sum(v for k, v in kinds.items() if k not in RECORD_KINDS)

TABLE = {
    "schema": "xpsloss-m4b2-expected-states-v1-draft4",
    "status": "fourth draft before the freeze (after the third Codex review)",
    "prereg": "docs/M4B2_PREREG.md",
    "design": "docs/M4B_DESIGN.md (6th edition), sections 5, 7, 9, 10; scope M4b-2 (U11 (a))",
    "model_id": MODEL,
    "baseline": {"commit": BASELINE, "package_version": "0.4.0.dev0",
                 "note": "main at the completion of M4b-1; src/ equals e7fa6a3; no tag exists for it",
                 "without_git_history": "not_run/sdist_without_git_history (not counted as passed)"},
    "copied_from_m4b1": {"file": "benchmarks/numerics/m4b1/expected_states_m4b1.json", "sha256": M4B1_SHA256,
                         "items": ["kernels K_A1s, K_A3, K_T1, K_T2", "primaries D1000",
                                   "output_edges E_M, E_A3", "materials Al, synthetic_wide",
                                   "defaults of the first-order and order settings"]},
    "input_files": FILES,
    "schemas": {"kinetic": SCHEMA_KE, "measured": SCHEMA_BE},
    "defaults": {"analytic": ANALYTIC, "public": PUBLIC, "theta_deg": 0.0,
                 "primary_reference_choice": DEFAULT_CHOICE,
                 "note": "every argument not listed in a row keeps the M4b-1 default (prereg M4b-1 section 2)"},
    "materials": MATERIALS, "kernels": KERNELS, "axes": AXES, "records": RECORDS,
    "inputs_are_synthetic": "every axis, record, binding primary and binding edge list of this table",
    "record_fields": RECORD_FIELDS, "not_recorded_order": NOT_RECORDED_ORDER,
    "defect_operations": {"set": "set the value at the dotted path (API: the dataclass field; CLI: the JSON key, "
                                 "created if absent); {'float': 'nan'} stands for the float nan",
                          "delete_key": "CLI only: remove the JSON key",
                          "pass_plain_dict": "API only: pass a plain dict with the same fields"},
    "binding_primaries": BE_PRIMARIES, "binding_output_edges": BE_EDGES,
    "kinetic_primaries": KE_PRIMARIES, "kinetic_output_edges": KE_EDGES,
    "cli": {"common_argv": COMMON, "option_arity": ARITY, "required_options": REQUIRED,
            "registered_option_values": {k: sorted(map(list, v)) for k, v in FIXED_VALUES.items()}},
    "ci_time_basis": {"note": "Codex, this PC, current KE code as a stand-in: the public CI rows of draft 2 took "
                              "128.64 s; with this CI set 11 public-equivalent runs took 38.99 s including the "
                              "baseline child, and the analytic rows about 1.83 s; not a measurement of M4b-2",
                      "seconds_per_row_draft2_measurement": {"BE7": 13.16, "BZ0": 10.83, "BC07": 7.00,
                                                             "BE5": 9.49, "BC01": 14.26, "BC02": 8.40,
                                                             "BC12": 9.57, "BASE-ms": 20.20, "BASE-areas": 10.95,
                                                             "BASE-cli": 16.65}},
    "key_names": {"regions": ["outside.low_binding", "bin[i] (input order)", "outside.high_binding"],
                  "components": COMPONENTS, "columns": COLUMNS, "error_parts": PARTS,
                  "read_from": {"zero_loss.bin[i]": "result.zero_loss[i]",
                                "zero_loss.outside.<side>": "result.outside[side]['zero_loss']",
                                "spectrum.outside.<side>": "result.outside[side]['spectrum']",
                                "order[n].bin[i].<c>": "getattr(result.orders[n-1], c)[i]",
                                "order[n].outside.<side>.<c>": "result.orders[n-1].outside[side][c]",
                                "KE control, outside zero loss": "control.single_loss.outside[side]['zero_loss']"},
                  "map_to_kinetic": "bin[i] (BE, input order) = KE bin[nb-1-i] for increasing BE input and KE "
                                    "bin[i] for decreasing BE input; outside.high_binding = KE outside.low; "
                                    "outside.low_binding = KE outside.high"},
    "vocabulary": {"entries": sorted(ENTRIES), "groups": GROUPS, "check_kinds": sorted(KINDS),
                   "new_reasons": sorted(NEW_REASONS), "removed_reason": "axis_conversion_changes_support",
                   "invalid_energy_axis_details": ["empty", "exclusive", "missing", "negative_kinetic_energy",
                                                   "not_finite", "not_recorded_not_allowed", "out_of_range",
                                                   "requires_binding_primary", "type", "unknown_kind"],
                   "warning_codes_in_order": ["binding_axis_assumed_calibrated", "energy_axis_field_not_recorded",
                                              "recorded_emission_angle_differs_from_theta",
                                              "axis_conversion_support_disagreement"]},
    "counts": {"rows": len(ROWS), "ci": N_CI, "full_only": N_FULL, "rows_by_group": by_group,
               "terminal_states": states, "ledger_records_by_kind": dict(sorted(kinds.items())),
               "ledger_fixed_schema_values": scalar_keys,
               "ledger_recursive_records": sum(kinds.get(k, 0) for k in RECORD_KINDS),
               "mutations_required": len(MUTATIONS), "mutations_from_design": 4, "mutations_added": 10},
    "limits_seconds": {"ci_increment": 60, "full": 1800, "row": 120, "reference": 300, "mutations": 3600},
    "ci_reduction_order": CI_ORDER, "ci_reduction_counts": ci_counts,
    "execution_order": ids,
    "rows": ROWS, "mutations": MUTATIONS,
}
TEXT = json.dumps(TABLE, indent=1, sort_keys=True, ensure_ascii=True, allow_nan=False) + "\n"
DATA = TEXT.encode("ascii")
if sys.argv[1:] == ["--check"]:
    need(OUT.is_file() and OUT.read_bytes() == DATA, "expected_states_m4b2.json differs from the generator")
    print("check passed")
else:
    OUT.write_bytes(DATA)
print("sha256", hashlib.sha256(DATA).hexdigest())
print("rows", len(ROWS), "ci", N_CI, "full-only", N_FULL, "groups", by_group, "states", states)
print("ledger", dict(sorted(kinds.items())), "values", scalar_keys)
print("mutations required", len(MUTATIONS))
