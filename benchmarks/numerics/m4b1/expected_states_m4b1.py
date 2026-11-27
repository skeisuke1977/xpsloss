"""M4b-1 pre-registration: the expected-state table (docs/M4B1_PREREG.md section 5).

Writes expected_states_m4b1.json next to this file (ASCII, LF) and prints its SHA-256 and the counts.
One row is one new entry point call (entry, inputs, settings) with its own expected terminal state.
"in_row_checks" are control runs and comparisons made inside a row; "reuse" names the artifacts of
other rows that are really read ("ROW.artifact").  Neither adds a row.  "check_plan" is the ledger
of planned keys, written before any reference value exists, as compact groups
{"keys", "kind", "oracle", "count"}; the counts are the denominators of the pre-registration.

The script checks, with explicit exceptions, that ids and normalized entries are unique, that every
reason, entry, group, state and check kind is in the registered vocabulary, that every public reason
of M4b-1 has at least one row, that reuse targets exist, form no cycle and never point from a CI row
to a full-only row, and that every mutation names existing rows.  It imports only the standard
library (not xpsloss, not NumPy).  The total number of rows and the CI/full split are asserted against
registered constants.  Third draft (after the second Codex review): key groups are expanded to concrete keys
at generation time (never written to the JSON) and checked for empty and duplicate expansions; entries
are compared by a typed signature after expanding the defaults; reused artifacts must be provided by
the source row; every mutation is bound to concrete (row, key, kind, oracle) records or to the terminal
state of a row.  Fourth draft (after the third Codex review): the signature expands every default of
prereg 4.1 (an explicit multiplier="weight_length" equals the omitted one); every ledger key and every
mutation target is checked against the orders, bins and regions allowed by the row (or by the in-row
control it names); three counterexamples are run as self-tests that must raise TableError.

Usage: python benchmarks/numerics/m4b1/expected_states_m4b1.py          (write the JSON)
       python benchmarks/numerics/m4b1/expected_states_m4b1.py --check  (compare byte for byte)
"""
import hashlib
import json
import math
import pathlib
import re
import sys

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "expected_states_m4b1.json"
INPUTS = {"Al": ("benchmarks/inputs/al_elfdb_v1.0.json", "benchmarks/inputs/al_elfdb_v1.0.bulk.json"),
          "Si": ("benchmarks/inputs/si_elfdb_v1.0.json", "benchmarks/inputs/si_elfdb_v1.0.bulk.json")}


class TableError(RuntimeError):
    pass


def need(cond, message):
    if not cond:
        raise TableError(message)


if sys.argv[1:] not in ([], ["--check"]):
    raise SystemExit("usage: expected_states_m4b1.py [--check]")

ROOT = next((p for p in HERE.parents if (p / INPUTS["Al"][0]).is_file()), None)
need(ROOT is not None, "repository root with benchmarks/inputs not found")


def hx(x):
    return float(x).hex()


def num(x):
    x = float(x)
    if not math.isfinite(x):
        return {"invalid": repr(x)}
    return {"value": x, "hex": x.hex()}


def nums(seq):
    return {"values": [float(x) for x in seq], "hex": [hx(x) for x in seq]}


# ---------------------------------------------------------------- materials
def _sha(rel):
    return hashlib.sha256((ROOT / rel).read_bytes()).hexdigest()


def _material(name):
    record, bulk = INPUTS[name]
    cfg = json.loads((ROOT / record).read_bytes().decode("utf-8"))["material"]
    run = json.loads((ROOT / bulk).read_bytes().decode("utf-8"))["material"]
    lines = [[float(o[k]) for k in ("energy_eV", "strength_eV2", "width_eV", "dispersion")]
             for o in cfg["oscillators"]]
    run_lines = [[float(o[k]) for k in ("energy_eV", "strength_eV2", "width_eV", "dispersion")]
                 for o in run["oscillators"]]
    need(lines == run_lines and float(cfg["gap_eV"]) == float(run["gap_eV"]),
         f"{name}: record and bulk run config disagree")
    return lines, float(cfg["gap_eV"]), record, bulk


def mat(gap, lines, origin, synthetic=True, **more):
    lines = [[float(v) for v in line] for line in lines]
    return {"gap_eV": float(gap), "gap_hex": hx(gap), "oscillators": lines,
            "oscillators_hex": [[hx(v) for v in line] for line in lines],
            "columns": ["energy_eV", "strength_eV2", "width_eV", "dispersion"],
            "synthetic": synthetic, "from": origin, **more}


AL, G_AL, AL_REC, AL_BULK = _material("Al")
SI, G_SI, SI_REC, SI_BULK = _material("Si")
need(G_AL == 0 and G_SI == 1.12 and AL == [[15.6, 243.3598, 1.5, 1.0]] and len(SI) == 3,
     "unexpected Al or Si input")
MATERIALS = {
    "Al": mat(G_AL, AL, "ELF database v1.0 (CC BY 4.0)", synthetic=False, file=AL_REC, sha256=_sha(AL_REC),
              run_config=AL_BULK, run_config_sha256=_sha(AL_BULK)),
    "Si": mat(G_SI, SI, "ELF database v1.0 (CC BY 4.0)", synthetic=False, file=SI_REC, sha256=_sha(SI_REC),
              run_config=SI_BULK, run_config_sha256=_sha(SI_BULK)),
    "AlZ": mat(G_AL, [[e, 0.0, w, d] for e, _, w, d in AL], "Al, strengths 0 (structural zero)"),
    "synthetic_wide": mat(0.0, [[20, 200, 8, 0]], "M4a pre-registration, cheap CI material, no dispersion"),
}

# ---------------------------------------------------------------- vocabulary and defaults
COLUMNS = ["term_13a", "term_13b", "term_14", "term_15_T1", "term_15_T2", "term_15_T3", "term_16_e",
           "term_16_h", "no_hole", "hole", "total"]
COMPONENTS = ["total", "no_hole", "hole"]
MODEL_ID = "first_loss_13a_chain_v0"
REFERENCE_MODELS = ["full_kernel_geometric_ref", "full_kernel_factorial_ref", "excess_poisson_ref"]
INPUT_REASONS = ["invalid_model", "invalid_max_order", "invalid_neglect_tolerance", "order_limit_exceeded",
                 "work_lower_bound_exceeded", "invalid_schema", "invalid_settings"]
DEFERRED_REASONS = ["invalid_energy_axis", "axis_conversion_not_injective", "axis_conversion_changes_support"]
CONVERGENCE_REASONS = ["order_budget_not_met", "inherited_error_exceeds_order_budget", "order_mesh_unresolved",
                       "convolution_work_limit_exceeded", "convolution_work_limit_reached",
                       "truncation_bound_exceeds_tolerance", "truncation_bound_unavailable", "nonfinite_output"]
INHERITED_REASONS = ["invalid_window", "multiplier_undefined", "weight_length_undefined", "loss_budget_not_met"]
TRUNCATION_STATES = ["complete_on_display", "bounded_on_display", "truncated"]
NEGLECTED_REASONS = ["diverges", "undetermined", "repeated_kernel_not_nonnegative"]
BOUND_REASONS = ["complete", "no_usable_tau", "not_representable"]
CONVENTION_RATIOS = ["exactly_one", "below_one", "not_determined_by_convention"]
ENTRIES = {
    "analytic_spectrum": "xpsloss._loss_orders.integrate_orders(kernel, window, primary=..., output_edges_eV=...)",
    "analytic_areas": "xpsloss._loss_orders.integrate_orders(kernel, window) without primary and edges",
    "layer_c": "xpsloss._loss_orders.convolve_panels(first_panels, window, max_order=..., features=...)",
    "reference_model": "xpsloss._loss_reference_models.integrate_orders_reference(kernel, window, model=...)",
    "api_spectrum": "xpsloss.multiple_loss_spectrum",
    "api_areas": "xpsloss.multiple_loss_order_areas",
    "cli": "xpsloss xps-multiple-loss (called as xpsloss.cli.main(argv) in the same process)",
    "cli_capabilities": "xpsloss capabilities",
}
GROUPS = ["RM-A1", "RM-A2", "RM-A3", "RM-A4", "RM-A5", "RM-B", "I", "INV", "REF", "M4A", "T", "N", "C", "CT"]
DENOMINATORS = ["numeric", "state", "bitwise", "replay_and_state", "identity"]
CHECK_KINDS = ["state", "direct_reference", "bound_reference", "error_part", "identity", "exact_zero_or_null",
               "bitwise", "replay", "diagnostic", "interval_contains", "recomputed_from_report", "inequality"]
EXTERNAL_REFERENCE_KINDS = ("direct_reference", "bound_reference", "error_part", "interval_contains")
STATUSES = ["success", "error", "not_converged", "complete"]
DEFAULTS = {
    "public": {"electron_energy_eV": 1000.0, "weight_length_A": 20.0, "hole_charge": 1.0,
               "epsilon_momentum": "k_total", "theta_deg": 0.0, "weight_length_kind": "external_constant",
               "extra_inverse_length_per_A": 0.0, "weight_length_seed_A": None, "multiplier": "weight_length",
               "multiplier_A": None},
    "analytic": {"weight_length_A": 2.0, "hole_charge": 1.0, "multiplier": "weight_length", "multiplier_A": None},
    "first_order_settings": {"epsrel": 1e-6, "area_epsabs": 1e-10, "spectrum_epsabs_fraction": 1e-9,
                             "panel_limit": 400, "max_cells": 2000000, "kernel_epsabs": 1e-12,
                             "kernel_epsrel": 1e-9, "kernel_limit": 200},
    "order_settings": {"order_epsrel": 1e-6, "order_epsabs_fraction": 1e-9, "order_limit": 64,
                       "max_convolution_points": 200000000, "neglect_tolerance_fraction": None},
    "test_only_settings": {"order_panel_limit": "same as panel_limit", "order_min_panel_width_eV": None},
    "model": MODEL_ID,
}
RESOURCE_FORMULAS = {
    "tau_grid_per_eV": "10**(-3 + i/20) for i = 0, 1, ... while the value is <= 40/a (a = lower window edge in eV); "
                       "93 points for a = 1 eV; tau -> 0 is a further candidate only when the upper estimate of "
                       "||R||_1 is below 1",
    "work_lower_bound_points": "527 * sum(n for n in 2..N) (31 interior points per panel, n panels, one interval)",
    "work_lower_bound_cells": "n * (output bins + 2) for each folded order n in 2..N, compared with max_cells per "
                              "order; not applied to the area entry",
    "work_estimate_after_first_order": "280.5 * sum(n * (P1 + n - 1) for n in 2..N), P1 = final first-order panels",
    "max_convolution_points_unit": "17 * (number of positive-width intersection intervals really integrated), "
                                   "cumulative over orders 2..N",
    "note": "the design table (section 14) proposed 561*sum(n^2); neither it nor 561*sum(n) is a lower bound of "
            "the points really evaluated (Codex review, finding 5); equality with the limit passes",
}
BASELINE = {"tag": "v0.3.0", "commit": "c40a38940d52f5e637ee987b9c6aad911adca4dc",
            "equivalent_src": "5a9f9b3742a2049a43da887d80161c707314ff74 (src differs only in __version__)",
            "version_fields": "each side reports its own __version__; never compared by value"}

# ---------------------------------------------------------------- analytic kernels (A^-1 eV^-1)
P2M10 = 2.0 ** -10


def box(value, lo, hi):
    return {"type": "box", "value": num(value), "interval_eV": [num(lo), num(hi)]}


def kernel(window, columns, features, text):
    feats = {"jumps_eV": [], "kinks_eV": [], "singularities": []}
    feats.update(features)
    return {"window_eV": [num(window[0]), num(window[1])], "columns": columns, "features": feats,
            "definition": text,
            "derived_columns": "no_hole = 13a+13b+15_T1+15_T3+16_e, hole = 14+15_T2+16_h, total = no_hole + Q_h*hole; "
                               "columns not listed are 0"}


A3_COLUMNS = {"term_13a": [box(0.04, 2, 5)], "term_15_T1": [box(0.03, 1, 3)], "term_14": [box(0.02, 3, 6)],
              "term_15_T2": [box(-0.05, 1.5, 2.5)]}
A3_FEATURES = {"jumps_eV": [1.5, 2, 2.5, 3, 5, 6]}
KERNELS = {
    "K_A1": kernel([1.3, 7.7], {"term_13a": [box(0.05, 1.3, 7.7)]}, {"jumps_eV": [1.3, 7.7]},
                   "box, non-aligned edges; K1 = R"),
    "K_A1s": kernel([2, 5], {"term_13a": [box(0.1, 2, 5)]}, {"jumps_eV": [2, 5]}, "box with dyadic edges; K1 = R"),
    "K_A2": kernel([1.3, 7.7], {"term_13a": [{"type": "exponential", "amplitude": num(0.08),
                                               "decay_per_eV": num(0.25), "interval_eV": [num(1.3), num(7.7)]}]},
                   {"jumps_eV": [1.3, 7.7], "min_variation_scale_eV": 4.0}, "0.08*exp(-0.25*T) on the window"),
    "K_A3": kernel([1, 8], A3_COLUMNS, A3_FEATURES, "boxes of different widths; negative hole box"),
    "K_A3z": kernel([1, 8], {**A3_COLUMNS, "term_15_T2": [box(-0.06, 1.5, 2.5)]}, A3_FEATURES,
                    "as K_A3 with a hole column of zero area (0.02*3 - 0.06*1)"),
    "K_A4_9": kernel([2, 5], {"term_13a": [{"type": "legendre", "scale": num(0.05), "degree": 9,
                                             "form": "scale*(1 + P_degree(u)), u = (2*T - 7)/3",
                                             "interval_eV": [num(2), num(5)]}]},
                     {"jumps_eV": [2, 5]}, "0.05*(1 + P9(u)) on [2,5]"),
    "K_A4_16": kernel([2, 5], {"term_13a": [{"type": "legendre", "scale": num(0.05), "degree": 16,
                                              "form": "scale*(1 + P_degree(u)), u = (2*T - 7)/3",
                                              "interval_eV": [num(2), num(5)]}]},
                      {"jumps_eV": [2, 5]}, "0.05*(1 + P16(u)) on [2,5]"),
    "K_A5": kernel([2, 12], {"term_13a": [{"type": "lorentzian", "amplitude": num(0.05), "center_eV": num(5.3),
                                            "half_width_eV": num(0.4), "interval_eV": [num(2), num(12)]}]},
                   {"singularities": [[5.3, 0.4]]}, "0.05/(1 + ((T-5.3)/0.4)^2) cut by the window"),
    "K_T1": kernel([1, 2.5], {"term_13a": [box(0.2, 1, 2.5)]}, {"jumps_eV": [1, 2.5]}, "box; a = 1 eV exactly"),
    "K_T2": kernel([1, 50], {"term_13a": [box(0.0072, 1, 50)], "term_15_T1": [box(0.005, 1, 50)]},
                   {"jumps_eV": [1, 50]}, "boxes; with m = 2: r = 0.7056, w = 1.1956"),
    "K_T4_one": kernel([1, 5], {"term_13a": [box(0.125, 1, 5)]}, {"jumps_eV": [1, 5]}, "with m = 2: r = 1 exactly"),
    "K_T4_above": kernel([1, 5], {"term_13a": [box(0.125 * (1 + P2M10), 1, 5)]}, {"jumps_eV": [1, 5]},
                         "with m = 2: r = 1 + 2^-10"),
    "K_T4_below": kernel([1, 5], {"term_13a": [box(0.125 * (1 - P2M10), 1, 5)]}, {"jumps_eV": [1, 5]},
                         "with m = 2: r = 1 - 2^-10"),
    "K_T5": kernel([1, 50], {"term_13a": [box(0.01, 1, 50)]}, {"jumps_eV": [1, 50]}, "box for huge multipliers"),
    "K_SGN": kernel([2, 5], {"term_13a": [box(0.10, 2, 5)],
                             "term_14": [box(0.04, 2, 3.5), box(-0.04, 3.5, 5)]}, {"jumps_eV": [2, 3.5, 5]},
                    "positive 13a box and a hole box that changes sign at 3.5 eV (zero area)"),
    "K_LIN": kernel([2, 5], {"term_13a": [box(0.10, 2, 5)],
                             "term_14": [{"type": "linear", "slope_per_eV": num(0.02), "zero_eV": num(3),
                                          "interval_eV": [num(2), num(5)]}]}, {"jumps_eV": [2, 5]},
                    "hole = 0.02*(T-3): changes sign at the centre of the first-order panel [2,4]"),
    "K_STEP": kernel([2, 5], {"term_13a": [box(k / 16, 2 + 0.5 * (k - 1), 2 + 0.5 * k) for k in range(1, 7)]},
                     {"jumps_eV": [2, 2.5, 3, 3.5, 4, 4.5, 5]},
                     "six steps k/16 of width 0.5 eV: 6 first-order panels and 12 mandatory order-2 panels"),
    "K_R0": kernel([2, 5], {"term_15_T1": [box(0.05, 2, 5)]}, {"jumps_eV": [2, 5]}, "R = 0 and K1 > 0"),
}

# ---------------------------------------------------------------- primaries and output edges
PRIMARIES = {
    "D1000": {"kind": "delta", "energies_eV": nums([1000]), "intensities": nums([1])},
    "D1000zero": {"kind": "delta", "energies_eV": nums([1000]), "intensities": nums([0])},
    "D1008": {"kind": "delta", "energies_eV": nums([1008]), "intensities": nums([1])},
    "D998": {"kind": "delta", "energies_eV": nums([998]), "intensities": nums([1])},
    "P_A1": {"kind": "delta", "energies_eV": nums([1000, 1012]), "intensities": nums([1, 0.5])},
    "B_A1p": {"kind": "binned", "edges_eV": nums([999.5, 1000, 1000 + P2M10]), "intensities": nums([0.6, 0.4])},
    "B5": {"kind": "binned", "edges_eV": nums([998, 998.5, 999, 999.5, 1000, 1000.5]),
           "intensities": nums([0.1, 0.2, 0.4, 0.2, 0.1])},
}
E_M = [960, 980, 984, 988, 995, 999, 1000, 1001]
E_A1S = [988, 990, 994, 995, 996, 998, 1000, 1001]
E_T1 = [997, 998, 999, 1000, 1001]
EDGES = {
    "E_A1": nums([960.3, 970.1, 980.7, 985.2, 990.9, 993.4, 996, 998.5, 1000, 1001]),
    "E_A1s": nums(E_A1S),
    "E_A1s_plus8": nums([e + 8 for e in E_A1S]),
    "E_A3": nums([968, 972, 976, 980, 984, 988, 990, 992, 994, 995, 996, 997, 998, 999, 1000, 1001]),
    "E_A5": nums([962, 970, 976, 980, 983, 985, 987, 989, 991, 993, 995, 997, 999, 1000, 1001]),
    "E_B120": {"grid": {"min": 880.5, "max": 1000.5, "step": 1}, **nums([880.5 + i for i in range(121)])},
    "E_M": nums(E_M),
    "E_T1_eq": nums(E_T1),
    "E_T1_up": nums([math.nextafter(997.0, math.inf)] + E_T1[1:]),
    "E_T1_down": nums([math.nextafter(997.0, -math.inf)] + E_T1[1:]),
    "E_T2": nums([950, 970, 980, 990, 999, 1000, 1001]),
    "E_SGN_rel": nums([989, 990, 992, 993, 993.125, 994, 996, 1000, 1001]),
    "E_SGN_abs": nums([989, 990, 992.5, 993.5, 996, 1000, 1001]),
    "E_D100": nums([900, 1000, 1001]),
    "E_DS2": nums([997.5, 998.5, 999.5, 1000.5]),
    "E_DS0": nums([1000, 1001, 1002]),
}


def nbins(name):
    return len(EDGES[name]["values"]) - 1


# ---------------------------------------------------------------- check plans
def plan(keys, kind, oracle, count=None, summary=None):
    """A key group.  A list is a list of concrete keys (stored as the summary when one is given); a string with
    a count is a free-text group whose members are enumerated by index."""
    need(kind in CHECK_KINDS, f"unknown check kind {kind}")
    if isinstance(keys, str):
        need(isinstance(count, int) and count > 0, f"empty expansion: {keys}")
        concrete = [keys] if count == 1 else [f"{keys} [{i + 1}/{count}]" for i in range(count)]
        stored = keys
    else:
        concrete = list(keys)
        need(concrete, f"empty expansion: {summary}")
        need(count is None or count == len(concrete), f"count differs from the expansion: {summary}")
        stored = summary if summary is not None else concrete
    need(len(set(concrete)) == len(concrete), f"duplicate keys in one group: {stored}")
    return {"keys": stored, "kind": kind, "oracle": oracle, "count": len(concrete), "_concrete": concrete}


def region_names(edges):
    return ["outside.low"] + [f"bin[{i}]" for i in range(nbins(edges))] + ["outside.high"]


STATE = plan("terminal_state", "state", "expect", 1)


def bins_plan(edges, orders, oracle, kind="direct_reference", comps=3):
    lo, hi = orders
    return plan([f"order[{n}].{reg}.{c}" for n in range(lo, hi + 1) for reg in region_names(edges)
                 for c in COMPONENTS[:comps]], kind, oracle,
                summary=f"order[{lo}..{hi}].<outside.low|bin[i]|outside.high>.<total|no_hole|hole>")


def areas_plan(orders, oracle, kind="direct_reference"):
    lo, hi = orders
    return plan([f"order[{n}].area.{c}" for n in range(lo, hi + 1) for c in COLUMNS], kind, oracle,
                summary=f"order[{lo}..{hi}].area.<11 columns>")


def sens_plan(orders, oracle):
    lo, hi = orders
    return plan([f"order[{n}].cutoff_sensitivity.{c}" for n in range(lo, hi + 1) for c in COLUMNS],
                "direct_reference", oracle, summary=f"order[{lo}..{hi}].cutoff_sensitivity.<11 columns>")


def sum_plan(edges, oracle, kind="identity"):
    return plan([f"spectrum.{reg}" for reg in region_names(edges)], kind, oracle,
                summary="spectrum.<outside.low|bin[i]|outside.high>")


def identity_plans(edges, orders):
    lo, hi = orders
    n = hi - lo + 1
    return [plan(f"order[{lo}..{hi}]: total = no_hole + Q_h*hole in every region and area", "identity",
                 "prereg 7.4", n * (nbins(edges) + 3)),
            plan(f"order[{lo}..{hi}]: sum over all regions = S*w_n (3 components)", "identity", "prereg 7.4", n * 3),
            sum_plan(edges, "prereg 7.4: J0 + sum of orders, returned errors added")]


rows = []


def ok(state=None, **attributes):
    e = {"status": "success", **attributes}
    if state is not None:
        need(state in TRUNCATION_STATES, state)
        e["truncation_state"] = state
    return e


def rejected(reason, exception="MultipleLossInputError", **attributes):
    return {"status": "error", "exception": exception, "reason": reason, "kernel_evaluations": 0, **attributes}


def unconverged(reason, exception="MultipleLossConvergenceError", **attributes):
    return {"status": "not_converged", "exception": exception, "reason": reason, **attributes}


def row(rid, group, entry, inputs, expect, *, ci, denominator, settings=None, reuse=(), checks=(), plans=()):
    need(group in GROUPS and entry in ENTRIES and denominator in DENOMINATORS, f"{rid}: vocabulary")
    need(expect["status"] in STATUSES, f"{rid}: status")
    provides = (["result"] if expect["status"] == "success" else []) + (["baseline"] if group == "M4A" else [])
    rows.append({"id": rid, "group": group, "entry": entry, "inputs": inputs, "settings": dict(settings or {}),
                 "expect": expect, "ci": bool(ci), "reuse": list(reuse), "provides": provides,
                 "in_row_checks": list(checks),
                 "denominator": denominator, "check_plan": [STATE] + list(plans)})


def an(kernel_name, max_order, primary="D1000", edges=None, **extra):
    d = {"kernel": kernel_name, "model": MODEL_ID, "max_order": max_order}
    if edges is not None:
        d.update(primary=primary, output_edges_eV=edges)
    d.update(extra)
    return d


LOOSE = {"epsrel": 0.5, "area_epsabs": 1.0, "spectrum_epsabs_fraction": 1.0}

# ---------------------------------------------------------------- RM-A1 box kernels
row("A1", "RM-A1", "analytic_spectrum", an("K_A1", 6, "P_A1", "E_A1"), ok("truncated", support_order=39),
    ci=True, denominator="numeric",
    checks=["rounding-level comparison (prereg 7.2, class exact_polynomial)",
            "order[n].mesh.knots_eV equals the mandatory knots k*a+(n-k)*b rounded once (diagnostic)",
            "order[3]: RN(2a+b) = RN(10.3) is an edge of the initial mesh and of the final mesh "
            "(mesh.initial_panel_edges_eV, mesh.panel_edges_eV and the returned panels), bit for bit",
            "both outside regions are non-zero in the reference for some order"],
    plans=[bins_plan("E_A1", (1, 6), "R-A exact_polynomial"), areas_plan((1, 6), "R-A exact_polynomial"),
           plan("order[2..6].mesh.knots_eV", "diagnostic", "rational knots", 5),
           plan(["order[3].mesh.initial_panel_edges_eV[RN(2a+b)]", "order[3].mesh.panel_edges_eV[RN(2a+b)]"],
                "diagnostic", "rational knots")]
    + identity_plans("E_A1", (1, 6)))
row("A1s", "RM-A1", "analytic_spectrum", an("K_A1s", 3, "D1000", "E_A1s"), ok("truncated", support_order=6),
    ci=True, denominator="numeric",
    checks=["bin[2] (loss 5..6 eV) holds second order only: orders 1 and 3 are exactly 0 there",
            "regions outside [n*a, n*b] are exactly 0 for each order (values and errors)"],
    plans=[bins_plan("E_A1s", (1, 3), "R-A exact_polynomial"), areas_plan((1, 3), "R-A exact_polynomial"),
           plan("order[n] regions with no positive-width overlap with [n*a, n*b]: 7, 6 and 6 regions for n = 1, 2, 3, "
                "three components each; includes order[2] in bin[6] (zero-loss bin)", "exact_zero_or_null",
                "support", (7 + 6 + 6) * 3)]
    + identity_plans("E_A1s", (1, 3)))
row("A1p", "RM-A1", "analytic_spectrum", an("K_A1s", 3, "B_A1p", "E_A1s"), ok("truncated"),
    ci=True, denominator="numeric",
    checks=["binned primary with a wide bin (0.5 eV) and a narrow bin (2^-10 eV)"],
    plans=[bins_plan("E_A1s", (1, 3), "R-A exact_polynomial"), areas_plan((1, 3), "R-A exact_polynomial")]
    + identity_plans("E_A1s", (1, 3)))

# ---------------------------------------------------------------- RM-A2 tilted box
row("A2", "RM-A2", "analytic_spectrum", an("K_A2", 4, "D1000", "E_A1"), ok("truncated", support_order=30),
    ci=True, denominator="numeric",
    checks=["Laplace transform of the returned unmultiplied panels (column term_13a) at tau = 0.05 and 0.3 per eV "
            "against the unmultiplied closed form; accuracy only, no coverage check",
            "moments M0 and M1 of the returned unmultiplied panels (column term_13a) against the closed form; "
            "accuracy only; the test quadrature error (24-point Gauss per panel) is added to u_r"],
    plans=[bins_plan("E_A1", (1, 4), "R-A"), areas_plan((1, 4), "R-A"), sens_plan((1, 4), "R-A"),
           plan([f"panels.laplace[{n}][{tau}]" for n in (2, 3, 4) for tau in ("0.05", "0.3")], "direct_reference",
                "R-A"),
           plan([f"panels.moment[{n}][{k}]" for n in (2, 3, 4) for k in (0, 1)], "direct_reference", "R-A")]
    + identity_plans("E_A1", (1, 4)))

# ---------------------------------------------------------------- RM-A3 two kinds of kernels
A3_NOTES = {
    "A3a": (1.0, {}, "base case; negative hole bins exist in the reference"),
    "A3b": (0.0, {}, "hole charge 0: the surface term 15_T1 must not repeat"),
    "A3c": (2.0, {}, "hole charge 2: the total kernel is negative on [1.5, 2); warning negative_total_kernel"),
    "A3d": (1.0, {"multiplier": "external", "multiplier_A": 3.0}, "m = 3 differs from L_w = 2: order n scales as m^n"),
    "A3e": (1.0, {"multiplier": "inverse_window_area"}, "S2: w = 1 and r = A_13a/A_total = 12/19"),
}
for rid, (charge, extra, note) in A3_NOTES.items():
    row(rid, "RM-A3", "analytic_spectrum", an("K_A3", 4, "D1000", "E_A3", hole_charge=charge, **extra),
        ok("truncated"), ci=True, denominator="numeric", checks=[note],
        plans=[bins_plan("E_A3", (1, 4), "R-A"), areas_plan((1, 4), "R-A"),
               plan(["order_area_ratio"], "interval_contains", "R-A")]
        + identity_plans("E_A3", (1, 4)))
row("A3f", "RM-A3", "analytic_areas",
    an("K_A3", 3, multiplier="inverse_window_area",
       controls={"control_external": {"same_as_row_except": {"multiplier": "external",
                                                               "multiplier_A": "the reported m of the row"},
                                      "max_order": 3},
                 "control_point_path": {"kernel": "K_A1s", "max_order": 2, "primary": "D1000",
                                        "output_edges_eV": "E_A1s", "initial_edges": [2, 5],
                                        "inject": {"entry": "inject", "roles": ["interpolation"],
                                                   "errors": {"term_13a": num(1e-3), "no_hole": num(1e-3),
                                                              "total": num(1e-3)}, "values": "unchanged"},
                                        "settings": "the loose settings of the row"}},
       inject={"entry": "inject", "roles": ["interpolation"],
               "errors": {"all 11 columns": num(1e-3)}, "values": "unchanged",
               "verification": "returned unchanged (values and errors)"}),
    ok(), ci=True, denominator="numeric", settings=LOOSE,
    checks=["first-order area and error: A_total = 0.19 (relative 1e-12) and e in [0.007, 0.007*(1 + 1e-9)] per angstrom",
            "the multiplier interval ends are rebuilt outside the product from the reported A and e: 1/(A+e) rounded "
            "down and 1/(A-e) rounded up, bit for bit (not only containment)",
            "order 1: area.total and its error are single_loss_weight and single_loss_weight_error bit for bit (S2: "
            "the error is |w_1 - 1| rounded up); the other ten columns follow the mapping of prereg 4.2",
            "orders 2 and 3: returned area error = e_f + rho_n*(|w| + e_f) with rho_n rebuilt from the reported "
            "interval and m by exact rationals, rounded up",
            "fixed-multiplier inherited part of order 2 (total) = m^2*(0.12*0.007 + 0.21*0.007 + 0.007^2) with the "
            "reported binary64 m and first-order errors: not below, and above by at most relative 1e-9",
            "order 3: the propagated part is rebuilt from the reported order-2 local part (converted to the same "
            "unit) and ||R||_1; the mixed term of that local part with e_R is added once to the inherited part; the "
            "local part is rebuilt from the per-panel local terms (mesh.local_area_terms); no independent "
            "reference exists for the local residuals",
            "in-row control (external): the same run with multiplier='external' and multiplier_A = the reported m "
            "gives the same fixed-multiplier parts for orders 2 and 3 (relative 4 eps) and zero interval parts",
            "in-row control (point path, X13): K_A1s, initial_edges=[2,5], interpolation-only error 1e-3 on term_13a, "
            "no_hole and total, the same loose budgets, D1000 and E_A1s, N = 2; the node parts handed to _fold_order "
            "are captured: at the mandatory knot s = 7 eV the unmultiplied inherited part of total is "
            "3*(0.1e + 0.1e + e^2) = 0.000603 within the absolute width 1e-12"],
    plans=[plan([f"order[{n}].area_error.{c}" for n in (2, 3) for c in COLUMNS], "recomputed_from_report",
                "in-row exact rational (rho_n rule)", summary="order[2..3].area_error.<11 columns>"),
           plan("order[1] area and error mapping (11 columns; total: convention part) and the interval ends",
                "recomputed_from_report", "in-row exact rational", 13),
           plan([f"order[2].area_error_parts[inherited][{c}]" for c in ("total", "no_hole", "hole", "term_13a")],
                "error_part", "R-A"),
           plan(["order[3].area_error_parts[inherited][total]", "order[3].area_error_parts[propagated][total]",
                 "order[3].area_error_parts[local][total]"], "recomputed_from_report",
                "reported order-2 parts and per-panel local terms with the R-A norms ||R||_1 and ||C_2||_1"),
           plan([f"control_external.order[{n}].area_error_parts[{p}][total]" for n in (2, 3)
                 for p in ("inherited", "propagated", "local")], "identity", "in-row control"),
           plan(["control_point_path.order[2].node_parts[inherited][total]@7eV"], "error_part",
                "R-A: 3*(0.1e + 0.1e + e^2) = 0.000603, absolute width 1e-12"),
           plan(["single_loss.window_area.total", "single_loss.window_area_error.total"], "direct_reference", "R-A")])
row("A3z", "RM-A3", "analytic_spectrum", an("K_A3z", 3, "D1000", "E_A3"), ok("truncated"),
    ci=True, denominator="numeric",
    checks=["hole has zero algebraic area and non-zero bins; first moment of orders 2 and 3 against the closed form"],
    plans=[bins_plan("E_A3", (1, 3), "R-A"), areas_plan((1, 3), "R-A"),
           plan([f"panels.moment[{n}][{k}].hole" for n in (2, 3) for k in (0, 1)], "direct_reference", "R-A")]
    + identity_plans("E_A3", (1, 3)))
row("AE1", "RM-A3", "analytic_areas",
    an("K_LIN", 2, initial_edges=[2, 4, 5],
       inject={"entry": "inject", "roles": ["interpolation"],
               "errors": {"term_13a": num(1e-3), "no_hole": num(1e-3), "total": num(1e-3)},
               "values": "unchanged", "verification": "returned unchanged (values and errors)"}),
    ok(), ci=True, denominator="numeric", settings=LOOSE,
    checks=["order[2].area_error_parts[inherited][hole] is not below E and above by at most relative 1e-9 plus the "
            "stated rounding, E = m^2*0.05*0.003 = 6e-4 (|hole| on [2,4] is split at its zero T = 3)",
            "a 17-point Gauss rule applied to |hole| on the unsplit panel [2,4] gives about 0.9978*E (below)",
            "the local part of the same key is rebuilt separately and is not counted twice",
            "precondition: the first-order panels are [2,4] and [4,5] (initial_edges; otherwise the row is not_run)"],
    plans=[plan(["order[2].area_error_parts[inherited][hole]"], "error_part", "R-A"),
           plan(["order[2].area_error_parts[local][hole]"], "recomputed_from_report",
                "sum of the per-panel local terms (mesh.local_area_terms), rounded outward"),
           plan(["order[2].area_error_parts[local][hole] <= local budget"], "inequality", "prereg 5"),
           plan("single_loss panels", "diagnostic", "precondition", 1)])
row("OI1", "RM-A3", "analytic_spectrum",
    an("K_A1s", 2, "D1000", "E_A1s",
       order_inject={"entry": "order_inject", "order": 2, "roles": ["interpolation"],
                     "columns": ["total", "no_hole", "term_13a"],
                     "local_errors": "1e-3 on panels of width >= 3 eV (the two initial panels [4,7] and [7,10]); "
                                     "1e-12 on panels of width <= 1.5 eV; unmultiplied panel units",
                     "values": "unchanged"}),
    ok("truncated"), ci=True, denominator="numeric",
    checks=["the initial order-2 mesh [4,7],[7,10] exceeds the local budget (iii) (recorded in mesh.history)",
            "order[2].mesh.stages >= 1 and every final panel has width <= 1.5 eV",
            "the final formal judgement passes; in bin[1] (KE 990..994, loss 6..10 eV, contains s = 7 eV) the local "
            "part is positive, the propagated part is exactly 0 and the inherited part is below 1e-10*|y| + 1e-13",
            "expected observations under X30: either inherited_error_exceeds_order_budget at order 2 (key in bin[1], "
            "bin[2] or bin[3], column total) or success with stages = 0; both break this row"],
    plans=[bins_plan("E_A1s", (1, 2), "R-A"),
           plan(["order[2].mesh.stages", "order[2].mesh.history", "order[2].mesh.final_panel_widths"], "diagnostic",
                "prereg 5"),
           plan([f"order[2].bin[1].total_error_parts[{p}]" for p in ("local", "propagated", "inherited")],
                "inequality", "prereg 5")])
row("OI2", "RM-A3", "analytic_spectrum",
    an("K_A1s", 3, "D1000", "E_A1s",
       order_inject={"entry": "order_inject", "order": 2, "roles": ["interpolation"],
                     "columns": ["total", "no_hole", "term_13a"],
                     "local_errors": "1e-12 on every order-2 panel; unmultiplied panel units", "values": "unchanged"}),
    ok("truncated", support_order=6), ci=True, denominator="numeric", reuse=["A1s.result"],
    checks=["the four inequalities read the total component of bin[1] (reached by orders 2 and 3; order 3 also "
            "reaches bin[0] and outside.low); nothing is required of hole, which is not injected",
            "order 3, bin[1].total: the propagated part (ii) is positive; the inherited part (i) agrees with that of "
            "A1s.result within relative 1e-6 plus 1e-18 (it carries the tiny mixed term of the injected local error)",
            "order 2, bin[1].total: the local part is positive and the propagated part is exactly 0",
            "order 2: mesh.panel_edges_eV equals that of A1s.result and mesh.builds = 1 (order 2 is entered once; no "
            "return to an earlier mesh)"],
    plans=[bins_plan("E_A1s", (1, 3), "R-A"),
           plan(["order[3].bin[1].total_error_parts[propagated] > 0",
                 "order[3].bin[1].total_error_parts[inherited] ~ A1s.result",
                 "order[2].bin[1].total_error_parts[local] > 0",
                 "order[2].bin[1].total_error_parts[propagated] = 0"], "inequality", "prereg 5"),
           plan(["order[2].mesh.panel_edges_eV"], "bitwise", "A1s.result"),
           plan(["order[2].mesh.builds"], "diagnostic", "prereg 4.2")])

# ---------------------------------------------------------------- RM-A4 polynomial kernels
for rid, kname, degree, value in (("A4a", "K_A4_9", 9, "0.05^2 * 3 * 18/19"), ("A4c", "K_A4_16", 16, "0.05^2 * 3 * 34/33")):
    row(rid, "RM-A4", "layer_c",
        {"panels": {"constructor": "Panel.from_lobatto(2, 5, values, errors)", "kernel": kname,
                    "values": "kernel at the 17 Lobatto nodes in term_13a, no_hole and total; errors 0"},
         "window": [num(2), num(5)], "max_order": 2, "model": MODEL_ID},
        ok(), ci=True, denominator="numeric",
        checks=[f"order-2 panels have a boundary at s = 7 eV and the end node value there (term_13a) is {value}; "
                "relative 1e-13 (the inner Gauss rule is exact for degree 32)",
                "a 9-point inner rule gives 0.05^2*3 for degree 9 (all nodes are zeros of P9)"],
        plans=[plan(["order_panels[2].value_at_7eV.term_13a", "order_panels[2].value_at_7eV.total"],
                    "direct_reference", "R-A exact rational")])
row("A4b", "RM-A4", "analytic_spectrum", an("K_A4_9", 3, "D1000", "E_A1s"), ok("truncated"),
    ci=True, denominator="numeric",
    checks=["order 2 is a piecewise polynomial of degree 19: the outer re-interpolation is not exact",
            "order[2].mesh.stages >= 1"],
    plans=[bins_plan("E_A1s", (1, 3), "R-A"), areas_plan((1, 3), "R-A"),
           plan("order[2].mesh.stages", "diagnostic", "prereg 5", 1)] + identity_plans("E_A1s", (1, 3)))

# ---------------------------------------------------------------- RM-A5 sharp peak
row("A5", "RM-A5", "analytic_spectrum", an("K_A5", 3, "D1000", "E_A5"), ok("truncated"),
    ci=False, denominator="numeric",
    checks=["order[2].mesh.stages >= 1 and order[3].mesh.stages >= 1"],
    plans=[bins_plan("E_A5", (1, 3), "R-A nested quadrature"), areas_plan((1, 3), "R-A"),
           plan("order[2..3].mesh.stages", "diagnostic", "prereg 5", 2)] + identity_plans("E_A5", (1, 3)))

# ---------------------------------------------------------------- RM-B real materials
RB_SET = {"epsrel": 1e-7}


def pub(material, max_order, primary="D1000", edges="E_B120", window=(1, 50), **extra):
    d = {"material": material, "loss_window_eV": [num(window[0]), num(window[1])], "model": MODEL_ID,
         "max_order": max_order}
    if edges is not None:
        d.update(primary=primary, output_edges_eV=edges)
    d.update(extra)
    return d


RB = (("RB1", "Al", (1, 50), "D1000", {}, True, []),
      ("RB2", "Al", (1, 50), "D1000", {"epsilon_momentum": "k_parallel"}, False, ["dispersive_k_parallel_13a"]),
      ("RB3", "Al", (1, 50), "B5", {}, False, []),
      ("RB4", "Si", (5, 50), "D1000", {}, False, []),
      ("RB5", "Si", (5, 50), "D1000", {"epsilon_momentum": "k_parallel"}, False, ["dispersive_k_parallel_13a"]),
      ("RB6", "Si", (5, 50), "B5", {}, False, []),
      ("RB7", "Al", (1, 50), "D1000", {"multiplier": "external", "multiplier_A": 40.0}, False,
       ["multiplier_differs_from_weight_length", "order_area_ratio_not_below_one"]),
      ("RB8", "Al", (1, 50), "D1000", {"hole_charge": 0.0}, False, []))
for rid, material, window, primary, extra, ci, warn in RB:
    expect = ok("truncated", warnings_include=["truncated"] + warn
                + (["resolution_basis_computed_candidates_only"] if material == "Si" else []))
    if rid == "RB7":
        expect["neglected_total_area"] = {"value": None, "reason": "diverges"}
    row(rid, "RM-B", "api_spectrum", pub(material, 3, primary, "E_B120", window, **extra), expect, ci=ci,
        denominator="numeric", settings=RB_SET,
        checks=["order 1 and the single_loss field are bitwise equal to single_loss_spectrum called in the same "
                "process with the same first-order arguments",
                "the arguments seen by _call_kernel (weight_length_A, hole_charge, epsilon_momentum)"],
        plans=[bins_plan("E_B120", (2, 3), "R-B"), areas_plan((2, 3), "R-B"),
               plan("order[1] and single_loss", "bitwise", "single_loss_spectrum, same process", 1),
               plan("order_area_ratio", "interval_contains", "R-B", 1)] + identity_plans("E_B120", (1, 3)))
row("RB1a", "RM-B", "api_areas", pub("Al", 3, edges=None), ok(), ci=False, denominator="numeric",
    settings=RB_SET, reuse=["RB1.result"],
    checks=["w_n, h_n and r agree with RB1.result within the sum of the two returned errors (not bitwise: the "
            "order meshes may differ)",
            "the single_loss field and order 1 are bitwise equal to loss_window_area called in the same process",
            "mesh.cells and mesh.formal_folds are 0 for every order (area entry)"],
    plans=[areas_plan((2, 3), "R-B"), sens_plan((1, 3), "R-B"),
           plan("order[1..3].area, cutoff_sensitivity and ratio against RB1.result", "identity", "prereg 7.4", 67),
           plan("order[1] and single_loss", "bitwise", "loss_window_area, same process", 1),
           plan("order_area_ratio", "interval_contains", "R-B", 1)])

# ---------------------------------------------------------------- identities
row("I-shift", "I", "analytic_spectrum", an("K_A1s", 3, "D1008", "E_A1s_plus8"), ok("truncated", support_order=6),
    ci=True, denominator="bitwise", reuse=["A1s.result"],
    checks=["all order values and errors are bitwise equal to A1s.result (exact rational differences)"],
    plans=[plan("orders, spectrum and errors against A1s.result", "bitwise", "A1s.result", 1)])
row("I-areas", "I", "analytic_areas", an("K_A3", 4), ok(), ci=True, denominator="identity", reuse=["A3a.result"],
    checks=["w_n agrees with A3a.result within the two returned errors and with the exact values"],
    plans=[areas_plan((1, 4), "R-A"), sens_plan((1, 4), "R-A"),
           plan("order[1..4].area against A3a.result", "identity", "prereg 7.4", 44)])

# ---------------------------------------------------------------- inverse relation (B1)
for rid, order in (("INV1", 2), ("INV2", 5)):
    row(rid, "INV", "analytic_spectrum", an("K_A3", order, "D1000", "E_A3"), ok("truncated"),
        ci=True, denominator="numeric",
        checks=["the product spectrum J^(N) is compared with the renewal-equation solution of the reference "
                "generator: (I - A_R) J = (I + A_gbar) F - A_{K1*R^N} F; the tail term reaches losses >= 1 + 2N eV",
                f"threshold {1 + 2 * order} eV: bins entirely below it are normal controls (the tail term is exactly "
                "0 there; only agreement is required); bins on the reaching side include the straddling bin, and the "
                "reference records at least one of them where the value without the tail term is outside the width",
                "the comparison with the product sum covers every bin and both outside regions"],
        plans=[sum_plan("E_A3", "R-INV", kind="direct_reference"),
               plan(["reference: reaching-side bin outside the width without the tail term"], "diagnostic", "R-INV")])

# ---------------------------------------------------------------- reference-only models
def ref(model, kname, order, **extra):
    return {"kernel": kname, "model": model, "max_order": order, "primary": "D1000", "output_edges_eV": "E_A3", **extra}


row("REF1", "REF", "reference_model", ref("full_kernel_geometric_ref", "K_A3", 3), ok(reference_only=True),
    ci=True, denominator="numeric", plans=[bins_plan("E_A3", (1, 3), "R-A (K1^{*n})")])
row("REF2", "REF", "reference_model", ref("full_kernel_factorial_ref", "K_A3", 3), ok(reference_only=True),
    ci=True, denominator="numeric", reuse=["REF1.result"],
    checks=["order n equals REF1 order n divided by n! within the returned errors"],
    plans=[bins_plan("E_A3", (1, 3), "R-A (K1^{*n}/n!)"),
           plan("order[2..3] against REF1.result/n!", "identity", "prereg 7.4", 2 * 17 * 3)])
row("REF3", "REF", "reference_model", ref("excess_poisson_ref", "K_A3", 3), ok(reference_only=True),
    ci=True, denominator="numeric", reuse=["A3a.result"],
    checks=["order 2 minus A3a order 2 equals gbar*gbar/2 of the reference (total and no_hole)",
            "hole is defined as total - no_hole"],
    plans=[bins_plan("E_A3", (1, 3), "R-A (model C)"),
           plan([f"order[2].{reg}.{c} minus A3a.result" for reg in region_names("E_A3")
                 for c in ("total", "no_hole")], "direct_reference", "R-A (gbar^{*2}/2)",
                summary="order[2].<outside.low|bin[i]|outside.high>.<total|no_hole> minus A3a.result")])
row("REF4a", "REF", "analytic_spectrum", an("K_R0", 1, "D1000", "E_A1s"),
    ok("truncated", support_order=6, neglected_total_area={"value": 0.0, "reason": None},
       neglected_bound="exactly 0 in every display bin"),
    ci=True, denominator="numeric",
    plans=[bins_plan("E_A1s", (1, 1), "R-A"),
           plan("truncation.neglected_at_or_above_region_bound[all bins] and neglected_total_area", "exact_zero_or_null",
                "R = 0", 8)])
for rid, model in (("REF4b", "full_kernel_geometric_ref"), ("REF4c", "full_kernel_factorial_ref")):
    row(rid, "REF", "reference_model", {**ref(model, "K_R0", 1), "output_edges_eV": "E_A1s"},
        ok(reference_only=True), ci=True, denominator="numeric",
        checks=["the model's own bound is positive in the bins reached by K1*K1 (loss 4..10 eV) and is not below "
                "the exact neglected quantity of that model"],
        plans=[plan([f"truncation.bound.bin[{i}].positive_where_reached_and_covers_exact" for i in range(7)],
                    "bound_reference", "R-A")])

# ---------------------------------------------------------------- M4A bitwise identity with v0.3.0
W13 = "inverse_13a_window_area"
CONVENTIONS = (
    ("c1", {}, True),
    ("c2", {"multiplier": "external", "multiplier_A": 20.0}, False),
    ("c3", {"multiplier": "external", "multiplier_A": 30.0}, False),
    ("c4", {"multiplier": "inverse_window_area"}, True),
    ("c5", {"weight_length_kind": W13, "weight_length_A": None, "extra_inverse_length_per_A": 0.01}, True),
    ("c6", {"weight_length_kind": W13, "weight_length_A": None, "extra_inverse_length_per_A": 0.01,
            "multiplier": "external", "multiplier_A": 30.0}, False),
    ("c7", {"weight_length_kind": W13, "weight_length_A": None, "extra_inverse_length_per_A": 0.01,
            "multiplier": "inverse_window_area"}, False),
)
M4A_CHECK = ["the single_loss field and order 1 are bitwise equal to single_loss_spectrum of the baseline child "
             "process (tag v0.3.0) with the same first-order arguments: all legacy fields, diagnostics and warnings",
             "the calls of _call_kernel (count, arrays, arguments) equal those of the baseline (the old entry of the "
             "baseline child has no delivery hook)",
             "on_final_panels is called exactly once in the current code (observed by wrapping _run) and "
             "diagnostics.first_order.beta_deliveries == 1"]
M4A_PLAN = [plan(["single_loss.all_legacy_fields", "order[1].values_and_errors"], "bitwise",
                 "baseline child process v0.3.0"),
            plan(["_call_kernel.call_sequence"], "diagnostic", "baseline child process v0.3.0"),
            plan(["on_final_panels.deliveries"], "diagnostic", "wrapped _run of the current code")]
for cid, extra, ci in CONVENTIONS:
    row(f"M4A-{cid}", "M4A", "api_spectrum", pub("synthetic_wide", 1, "D1000", "E_M", **extra), ok("truncated"),
        ci=ci, denominator="bitwise", checks=M4A_CHECK + [f"M4a-2 design 2.4 row {cid[1]}"], plans=M4A_PLAN)
for rid, cid, extra in (("M4A-o3", "c1", {}), ("M4A-o3s", "c4", {"multiplier": "inverse_window_area"})):
    row(rid, "M4A", "api_spectrum", pub("synthetic_wide", 3, "D1000", "E_M", **extra), ok("truncated"),
        ci=True, denominator="bitwise", reuse=[f"M4A-{cid}.result", f"M4A-{cid}.baseline"],
        checks=["the single_loss field and order 1 are bitwise equal to those of the max_order = 1 row and of the "
                "baseline; the _call_kernel calls are the same; at least one formal fold per higher order "
                "(diagnostic)",
                "the panel objects received by on_final_panels are captured directly and are the first-order input of "
                "order 2 (same objects or bitwise equal arrays as panels_by_order[1])"],
        plans=M4A_PLAN + [plan(["order[2].mesh.formal_folds >= 1", "order[3].mesh.formal_folds >= 1"], "diagnostic",
                               "prereg 4.6"),
                          plan(["on_final_panels.received_panels_are_order2_input"], "diagnostic",
                               "captured objects")]
        + identity_plans("E_M", (1, 3)))
row("M4A-Al", "M4A", "api_spectrum", pub("Al", 1, "D1000", "E_M"), ok("truncated"), ci=False,
    denominator="bitwise", checks=M4A_CHECK + ["conditions of the M4a row RB1"], plans=M4A_PLAN)
row("M4A-Z", "M4A", "api_spectrum", pub("AlZ", 2, "D1000", "E_M"),
    ok("truncated", kernel_evaluations=0, neglected_total_area={"value": 0.0, "reason": None}),
    ci=True, denominator="bitwise",
    checks=M4A_CHECK[:1] + ["orders 1 and 2: values, errors and areas are exactly 0; spectrum equals zero_loss"],
    plans=M4A_PLAN[:1] + [plan("order[1..2] values, errors, areas", "exact_zero_or_null", "structural zero",
                               2 * (9 * 3 + 11))])
for rid, extra, ci in (("M4A-area-c1", {}, True), ("M4A-area-c4", {"multiplier": "inverse_window_area"}, False)):
    row(rid, "M4A", "api_areas", pub("synthetic_wide", 1, edges=None, **extra), ok(), ci=ci, denominator="bitwise",
        checks=["the legacy area fields and order[1].area are bitwise equal to loss_window_area of the baseline; "
                "the _call_kernel calls are the same"],
        plans=[plan(["single_loss.all_legacy_fields", "order[1].area"], "bitwise", "baseline child process v0.3.0"),
               plan(["_call_kernel.call_sequence"], "diagnostic", "baseline child process v0.3.0")])

# ---------------------------------------------------------------- truncation
T1 = (("T1a", "E_T1_eq", "truncated", 3, "p_max - d_min = 3a exactly: the order with equality is kept"),
      ("T1b", "E_T1_up", "complete_on_display", 2, "d_min one ulp above 997"),
      ("T1c", "E_T1_down", "truncated", 3, "d_min one ulp below 997"))
for rid, edges, state, support, note in T1:
    e = ok(state, support_order=support)
    if state == "complete_on_display":
        e["neglected_bound"] = "exactly 0 in every display bin"
    row(rid, "T", "analytic_spectrum", an("K_T1", 2, "D1000", edges), e, ci=True, denominator="state",
        checks=[note], plans=[plan(["truncation.support_order", "truncation.state"], "diagnostic", "exact rational")])
T2_ZERO = plan("truncation bound of bin[4], bin[5] (U_D < (N+1)a): exactly 0, bound_reason complete, tau None",
               "exact_zero_or_null", "support", 6)
T2_BOUND = [plan([f"truncation.bound.bin[{i}].grid_minimum" for i in range(4)], "bound_reference",
                 "R-A: grid minimum of the design formula; lower side 1 - 1e-9, upper side 1 + 1e-6 (two verdicts)"),
            plan([f"truncation.bound.bin[{i}].covers_exact" for i in range(4)]
                 + ["truncation.bound.bin[0].ratio_le_10", "truncation.bound.bin[1].ratio_le_10"],
                 "bound_reference", "R-A Irwin-Hall sums"), T2_ZERO]
row("T2a", "T", "analytic_spectrum", an("K_T2", 4, "D1000", "E_T2", neglect_tolerance_fraction=0.05),
    ok("bounded_on_display", support_order=50), ci=True, denominator="state",
    checks=["the bound of bin[0] (U = 50 eV) is about 1.02e-2 (own exact estimate with exact Lambda; not a reference)"],
    plans=T2_BOUND)
row("T2b", "T", "analytic_spectrum", an("K_T2", 4, "D1000", "E_T2", neglect_tolerance_fraction=2.5e-3),
    unconverged("truncation_bound_exceeds_tolerance", required_order_estimate=">= 5", kernel_evaluations=">= 1"),
    ci=True, denominator="state",
    checks=["the tolerance lies between the bound of bin[1] (U = 30 eV, about 5.2e-4) and of bin[0] (about 1.02e-2)"])
row("T2c", "T", "analytic_spectrum", an("K_T2", 4, "D1000", "E_T2"), ok("truncated", support_order=50),
    ci=True, denominator="numeric",
    checks=["outside.low holds only computed orders (reference exact)",
            "bin[0..3]: finite bound with tau > 0 (||R||_1 < 1 allows tau -> 0 but it is not the minimum); "
            "bin[4..5] (U_D < 5a): (0.0, 'complete', None)"],
    plans=T2_BOUND + [bins_plan("E_T2", (1, 4), "R-A exact_polynomial"), areas_plan((1, 4), "R-A exact_polynomial"),
                      plan("neglected_total_area value S*w*r^4/(1-r)", "interval_contains", "R-A", 1)])
row("T2e", "T", "analytic_spectrum",
    an("K_T2", 4, "D1000", "E_T2",
       inject={"entry": "inject", "roles": ["interpolation"],
               "errors": {"term_13a": num(1.44e-4), "term_15_T1": num(1e-4), "no_hole": num(2.44e-4),
                          "total": num(2.44e-4)}, "values": "unchanged",
               "verification": "returned unchanged (values and errors)"}),
    ok("truncated"), ci=True, denominator="numeric", settings=LOOSE,
    checks=["Lambda uses the upper estimate: value part plus exp(-tau*t_lo)*e_panel on the first-order panels "
            "[1,2,4,8,16,32,50] (precondition: these panels; otherwise not_run)"],
    plans=[plan([f"truncation.bound.bin[{i}].grid_minimum_with_error_terms" for i in range(4)], "bound_reference",
                "R-A with the error terms (optimal tau chosen again), lower 1 - 1e-9, upper 1 + 1e-6"),
           plan("the same bins without the error terms differ by more than the width (recorded by the reference)",
                "diagnostic", "R-A", 1), T2_ZERO,
           plan("single_loss panels", "diagnostic", "precondition", 1)])
row("T4u", "T", "analytic_spectrum",
    an("K_T4_one", 3, "D1000", "E_A1s",
       inject={"entry": "inject", "roles": ["interpolation", "verification"],
               "errors": {"term_13a": num(1e-9), "no_hole": num(1e-9), "total": num(1e-9)}, "values": "unchanged"}),
    ok("truncated", neglected_total_area={"value": None, "reason": "undetermined"},
       warnings_include=["order_area_ratio_not_below_one"]),
    ci=True, denominator="state",
    checks=["the interval of r straddles 1; tau -> 0 is not a candidate",
            "U_D = 12, 10, 6, 5, 4, 2, 0 eV and (N+1)a = 4 eV: bin[0..4] have a finite bound with tau > 0 (the "
            "equality bin[4] is on the reached side); bin[5..6] are (0.0, 'complete', None)"],
    plans=[plan(["order_area_ratio.interval straddles 1"], "diagnostic", "prereg 5"),
           plan([f"truncation.bound.bin[{i}] finite with tau > 0" for i in range(5)], "diagnostic", "prereg 3.3"),
           plan([f"truncation.bound.bin[{i}] complete" for i in (5, 6)], "exact_zero_or_null", "support")])
row("T4d", "T", "analytic_spectrum", an("K_T4_above", 3, "D1000", "E_A1s"),
    ok("truncated", neglected_total_area={"value": None, "reason": "diverges"},
       warnings_include=["order_area_ratio_not_below_one"]),
    ci=True, denominator="state", plans=[plan("order_area_ratio.interval lower end >= 1", "diagnostic", "prereg 5", 1)])
row("T4v", "T", "analytic_spectrum", an("K_T4_below", 3, "D1000", "E_A1s"),
    ok("truncated", neglected_total_area={"value": "finite", "reason": None}),
    ci=True, denominator="numeric",
    plans=[plan("neglected_total_area interval contains S*w*r^3/(1-r) = 1024*(1-2^-10)^4", "interval_contains", "R-A", 1)])
row("T5", "T", "analytic_spectrum",
    an("K_T5", 2, "D1000", "E_T2", multiplier="external", multiplier_A=1e25, neglect_tolerance_fraction=1e-3),
    unconverged("truncation_bound_unavailable", kernel_evaluations=">= 1"), ci=True, denominator="state",
    checks=["Lambda_R >= 1 at every grid tau (at tau = 39.8 per eV it is about 1.3e4)"])
row("T5n", "T", "analytic_spectrum", an("K_T5", 2, "D1000", "E_T2", multiplier="external", multiplier_A=1e25),
    ok("truncated", neglected_bound="None with bound_reason no_usable_tau and tau None in bin[0..3]; exactly 0 with "
                                    "bound_reason complete in bin[4], bin[5] (U_D < 3a)",
       neglected_total_area={"value": None, "reason": "diverges"}),
    ci=True, denominator="state")
row("DS0", "T", "api_spectrum",
    {"material": "synthetic_wide", "loss_window_eV": [num(1), num(50)], "model": MODEL_ID,
     "max_order": "display_support", "primary": "D998", "output_edges_eV": "E_DS0"},
    ok("complete_on_display", support_order=0, computed=0, orders_empty=True, beta_deliveries=1,
       kernel_evaluations=">= 1", neglected_total_area={"value": "finite", "reason": None}),
    ci=True, denominator="state",
    checks=["the positive primary lies below the display: N_D = 0; M4a runs once (beta_deliveries = 1); orders = ()",
            "spectrum equals zero_loss in every region (outside.low holds J0 = 1); every display bound is exactly 0 "
            "with bound_reason complete",
            "neglected_total_area follows the N = 0 rule of the design (S*w/(1-r), r+ < 1); its interval is rebuilt "
            "from the w interval [m_-(A_total - e), m_+(A_total + e)] and the r interval by exact rational corners "
            "with outward rounding last",
            "order_area_ratio and single_loss are returned"],
    plans=[plan([f"spectrum.{reg} = zero_loss" for reg in ("outside.low", "bin[0]", "bin[1]", "outside.high")]
                + ["truncation.bound.bin[0] complete", "truncation.bound.bin[1] complete"], "exact_zero_or_null",
                "N_D = 0"),
           plan("neglected_total_area interval against the reported w and r intervals", "recomputed_from_report",
                "in-row exact rational", 1)])
row("DS1", "T", "analytic_spectrum", an("K_T1", "display_support", "D1000", "E_T1_eq"),
    ok("complete_on_display", support_order=3, computed=3), ci=True, denominator="numeric",
    checks=["display_support computes N = N_D = 3; the order with equality is kept and computed"],
    plans=[bins_plan("E_T1_eq", (1, 3), "R-A exact_polynomial"), areas_plan((1, 3), "R-A exact_polynomial"),
           plan("truncation bounds exactly 0 with bound_reason complete", "exact_zero_or_null", "support", 4)])

# ---------------------------------------------------------------- rejections and non-convergence
SW = dict(material="synthetic_wide", edges="E_M")


def n_api(rid, expect, ci=True, entry="api_spectrum", settings=None, checks=(), **kw):
    base = dict(SW)
    base.update(kw)
    material = base.pop("material")
    order = base.pop("max_order", 3)
    row(rid, "N", entry, pub(material, order, **base), expect, ci=ci, denominator="state", settings=settings,
        checks=checks)


n_api("N01", rejected("invalid_model"), model="full_kernel_geometric_ref",
      checks=["a reference-only model id is not accepted by the public entry"])
n_api("N02", rejected("invalid_model"), model="b1")
n_api("N03", rejected("invalid_max_order"), max_order=0)
n_api("N04", rejected("invalid_max_order"), max_order=2.5)
n_api("N05", rejected("invalid_max_order"), max_order=True)
n_api("N06", rejected("invalid_max_order"), max_order="auto")
n_api("N07", rejected("invalid_max_order"), entry="api_areas", edges=None, max_order="display_support",
      checks=["the area entry accepts integers only"])
n_api("N08", rejected("invalid_neglect_tolerance"), neglect_tolerance_fraction=-1.0)
n_api("N09", rejected("invalid_neglect_tolerance"), neglect_tolerance_fraction={"invalid": "nan"})
n_api("N10", rejected("invalid_neglect_tolerance"), neglect_tolerance_fraction=0.0)
n_api("N11", rejected("order_limit_exceeded", requested_order=65, order_limit=64), max_order=65)
n_api("N12", rejected("order_limit_exceeded", requested_order=100, order_limit=64), max_order="display_support",
      edges="E_D100", checks=["N_D = floor((1000 - 900)/1) = 100"])
n_api("N13", rejected("work_lower_bound_exceeded", resource="max_convolution_points", lower_bound=4743, limit=4500),
      max_order=4, settings={"max_convolution_points": 4500}, checks=["527*(2+3+4) = 4743 > 4500"])
n_api("N14", rejected("work_lower_bound_exceeded", resource="max_cells", order=34, lower_bound=306),
      max_order=40, settings={"max_cells": 300}, checks=["34*(7+2) = 306 > 300; orders up to 33 need at most 297"])
n_api("N15", unconverged("convolution_work_limit_exceeded", stage="after_first_order", kernel_evaluations=">= 1"),
      max_order=6, settings={"max_convolution_points": 50000},
      checks=["the input bound 527*20 = 10540 passes; the estimate 280.5*(20*P1 + 70) is at least 70125 for P1 >= 9",
              "the reported estimate equals the formula with the reported first-order panel count"])
row("N16", "N", "analytic_spectrum", an("K_STEP", 2, "D1000", "E_A1s"),
    unconverged("convolution_work_limit_reached", stage="order_refinement", order=2, limit=4500,
                kernel_evaluations=">= 1"),
    ci=True, denominator="state", settings={"max_convolution_points": 4500},
    checks=["the input bound 527*2 = 1054 passes; with P1 = 6 the estimate 280.5*2*7 = 3927 passes",
            "the 12 mandatory order-2 panels need at least 12*31*17 = 6324 points: the limit is reached during order 2",
            "precondition: 6 first-order panels (otherwise not_run); the measured count at the stop is reported"])
row("N18", "N", "analytic_spectrum", an("K_A5", 2, "D1000", "E_A5"),
    unconverged("order_budget_not_met", order=2), ci=True, denominator="state", settings={"order_panel_limit": 4},
    checks=["the order-2 mesh needs many more than 4 panels"])
SGN_INJECT = {"entry": "inject", "roles": ["interpolation", "verification"], "values": "unchanged"}
row("N19", "N", "analytic_spectrum",
    an("K_SGN", 2, "D1000", "E_SGN_rel",
       inject={**SGN_INJECT, "errors": {"term_13a": num(1.5e-8), "no_hole": num(1.5e-8), "total": num(1.5e-8)}}),
    unconverged("inherited_error_exceeds_order_budget", order=2,
                exceeded_includes={"key": "bin[3].hole", "side": "relative", "suggestion": "epsrel"}),
    ci=True, denominator="state", settings={"spectrum_epsabs_fraction": 1e-300, "area_epsabs": 1e-6},
    checks=["area_epsabs = 1e-6 keeps the order-2 area budget (0.84*2*1e-6) above the Young term (about 4.3e-8)",
            "the exception lists every exceeding key with inherited, both budgets and the side; no area key is listed",
            "bin[3] (loss 6.875..7 eV): |y| = m^2*0.004*0.0078125, inherited envelope integral = m^2*0.04*1.5e-8*0.3672, "
            "ratio 7.05e-6 against the budget 2e-6; the first order passes (relative point error 0.15*epsrel)"])
row("N20", "N", "analytic_spectrum",
    an("K_SGN", 2, "D1000", "E_SGN_abs",
       inject={**SGN_INJECT, "errors": {"term_13a": num(1e-8), "no_hole": num(1e-8), "total": num(1e-8)}}),
    unconverged("inherited_error_exceeds_order_budget", order=2,
                exceeded_includes={"key": "bin[2].hole", "side": "absolute",
                                   "suggestion": "spectrum_epsabs_fraction"}),
    ci=True, denominator="state", settings={"area_epsabs": 1e-6},
    checks=["area_epsabs = 1e-6 keeps the order-2 area budget (0.84*2*1e-6) above the Young term (about 2.9e-8)",
            "the exception lists every exceeding key with inherited, both budgets and the side; no area key is listed",
            "bin[2] (loss 6.5..7.5 eV) has y = 0 by antisymmetry; inherited envelope integral = m^2*0.04*1e-8*2.75 = "
            "4.4e-9 against the absolute budget phi_2*1e-9 = 8.4e-10"])
row("N21", "N", "analytic_spectrum", an("K_A5", 2, "D1000", "E_A5"),
    unconverged("order_mesh_unresolved", order=2), ci=True, denominator="state",
    settings={"order_min_panel_width_eV": 4.0},
    checks=["a panel that still needs a bisection cannot be split below the test-only minimum width"])
row("N26", "N", "analytic_spectrum", an("K_A1s", 2, "D1000", "E_A1s", multiplier="external", multiplier_A=1e200),
    unconverged("nonfinite_output", order=2), ci=True, denominator="state",
    checks=["m^2 overflows binary64 at order 2; order 1 is finite"])
n_api("N27", rejected("invalid_window", exception="SingleLossInputError"), window=(0, 50),
      checks=["the M4a exception passes through with its type, reason and attributes",
              "the old window check comes before the division by a in N_D"])
n_api("N28", rejected("multiplier_undefined", exception="SingleLossInputError"), material="AlZ",
      multiplier="inverse_window_area")
n_api("N29", rejected("weight_length_undefined", exception="SingleLossInputError"), material="AlZ",
      weight_length_kind=W13, weight_length_A=None, multiplier="inverse_window_area",
      checks=["weight_length_undefined comes before multiplier_undefined"])
n_api("N30", unconverged("loss_budget_not_met", exception="SingleLossConvergenceError", dominant_part="interpolation"),
      settings={"panel_limit": 9}, checks=["the first order fails as in the M4a row N25; no higher order starts"])
n_api("N31", rejected("invalid_settings", argument="order_epsrel"), settings={"order_epsrel": {"invalid": "nan"}})
n_api("N32", rejected("invalid_settings", argument="max_convolution_points"),
      settings={"max_convolution_points": True}, checks=["a bool is not accepted as a resource limit"])
n_api("N33", rejected("invalid_settings", argument="order_limit"), settings={"order_limit": 0})

# ---------------------------------------------------------------- CLI
CLI_OK = {"status": "complete", "exit_code": 0, "output_dir_after": "created"}
CLI_ERR = {"exit_code": 2, "output_dir_after": "absent"}
CLI_BASE = {"config_material": "synthetic_wide", "config_electron_energy_eV": 1000.0,
            "config_integration": {"epsabs": 1e-12, "epsrel": 1e-9, "limit": 200},
            "argv": ["--loss-window-eV", "1", "50", "--weight-length-A", "20", "--epsilon-momentum", "k_total",
                     "--model", MODEL_ID, "--max-order", "3"],
            "spectrum": {"schema": "xpsloss-single-loss-spectrum-v0", "axis": "kinetic_energy_eV",
                         "primary": "D1000", "output_edges_eV": "E_M"}}


def cli(rid, expect, checks=(), plans=(), ci=True, **changes):
    inputs = json.loads(json.dumps(CLI_BASE))
    for key, value in changes.items():
        if key == "argv_replace":
            for flag, new in value.items():
                i = inputs["argv"].index(flag)
                inputs["argv"][i + 1:i + 2] = [] if new is None else [new]
                if new is None:
                    inputs["argv"].remove(flag)
        elif key == "argv_add":
            inputs["argv"] += value
        elif key == "spectrum":
            inputs["spectrum"].update(value)
        else:
            inputs[key] = value
    row(rid, "C", "cli", inputs, expect, ci=ci, denominator="replay_and_state", checks=checks, plans=plans)


cli("C01", {**CLI_OK, "truncation_state": "truncated"},
    checks=["API replay from run.json effective_settings: all order values, errors and both CSV files are bitwise equal",
            "run.json: model_id, depends_on_model_id, repeated_kernel term_13a, hole_treatment first_event_only, "
            "kernel_energy_policy fixed_E0, truncation, order_area_ratio, comparison not_matched, "
            "physical_validity not_established, runtime, warnings, input SHA-256",
            "replaying run.json reproduces the input hashes, effective settings, state and warnings"],
    plans=[plan("CSV columns and run.json results against the API replay", "replay", "API", 3),
           plan("run.json record fields", "diagnostic", "prereg 4.3", 12)])
cli("C02", {**CLI_OK, "truncation_state": "complete_on_display"},
    checks=["N_D = floor((1000 - 997.5)/1) = 2; run.json truncation.computed = 2 and support_order = 2"],
    argv_replace={"--max-order": "display_support"}, spectrum={"output_edges_eV": "E_DS2"},
    plans=[plan("run.json truncation", "diagnostic", "prereg 4.3", 3)])
cli("C03", {**CLI_ERR, "status": "error", "error_type": "ValueError", "reason": None},
    argv_replace={"--model": None}, checks=["a missing required option is a usage error as in M4a C10c"])
cli("C04", {**CLI_ERR, "status": "error", "error_type": "ValueError", "reason": None},
    argv_replace={"--max-order": None})
cli("C05", {**CLI_ERR, "status": "not_converged", "error_type": "MultipleLossConvergenceError",
            "reason": "truncation_bound_exceeds_tolerance"},
    argv_replace={"--max-order": "1"}, argv_add=["--neglect-tolerance-fraction", "1e-30"],
    checks=["stderr JSON carries reason, order, stage, inner and required_order_estimate"])
cli("C06", {**CLI_ERR, "status": "not_converged", "error_type": "MultipleLossConvergenceError",
            "reason": "truncation_bound_unavailable"},
    argv_replace={"--max-order": "1"},
    argv_add=["--multiplier-A", "1e27", "--neglect-tolerance-fraction", "1e-3"],
    checks=["the reason differs from C05 and is compared separately"])
cli("C07", {**CLI_ERR, "status": "error", "error_type": "FileExistsError", "output_dir_after": "unchanged"},
    output_dir="existing",
    checks=["an existing output directory is refused; the bytes of the existing files are unchanged"])
cli("C08", {**CLI_ERR, "status": "error", "error_type": "MultipleLossInputError", "reason": "order_limit_exceeded"},
    argv_replace={"--max-order": "65"})
cli("C09", {**CLI_ERR, "status": "error", "error_type": "MultipleLossInputError", "reason": "invalid_schema"},
    spectrum={"axis": "binding_energy_eV"}, checks=["binding-energy input belongs to M4b-2"])
cli("C10", {**CLI_ERR, "status": "not_converged", "error_type": "MultipleLossConvergenceError",
            "reason": "convolution_work_limit_exceeded"},
    argv_replace={"--max-order": "6"}, argv_add=["--max-convolution-points", "50000"],
    checks=["stderr JSON: stage after_first_order, order null, inner null, estimate and limit; no file is created"])
row("C11", "C", "cli_capabilities", {"argv": ["capabilities"]}, {"status": "complete", "exit_code": 0},
    ci=True, denominator="replay_and_state",
    checks=["supported keeps the baseline elements in order and ends with xps_multiple_loss; the blocks of the "
            "baseline (including the per-command unsupported lists) are equal to the baseline; version follows the "
            "own-__version__ rule on both sides; the new block multiple_scattering names the model id, normal "
            "emission only, kernel fixed at E0, approximate, and the fixed unsupported list of prereg 4.3; the "
            "reference model ids do not appear anywhere"],
    plans=[plan("capabilities blocks", "diagnostic", "baseline child process v0.3.0", 4)])

# ---------------------------------------------------------------- controls
row("CT1", "CT", "analytic_spectrum", an("K_A3", 3, "D1000zero", "E_A3"),
    ok("complete_on_display", support_order=0), ci=True, denominator="numeric",
    checks=["all bins and both outside regions are exactly 0 with zero errors; the areas w_n are computed"],
    plans=[plan("order[1..3] regions and spectrum", "exact_zero_or_null", "F = 0", 3 * 17 * 3 + 17),
           areas_plan((1, 3), "R-A")])
row("CT3", "CT", "api_spectrum", pub("synthetic_wide", 3, "D1000", "E_M", epsilon_momentum="k_parallel"),
    ok("truncated", warnings_exclude=["dispersive_k_parallel_13a"]), ci=False, denominator="identity",
    reuse=["M4A-o3.result"],
    checks=["no dispersion: every order value agrees with M4A-o3.result (k_total) within the two returned errors"],
    plans=[plan("order[1..3] regions and areas against M4A-o3.result", "identity", "prereg 7.4", 3 * (9 * 3 + 11))])
CT_C = (("CT4", 0.0, "exactly_one", {"value": None, "reason": "undetermined"}),
        ("CT5", 1e-2, "below_one", {"value": "finite", "reason": None}),
        ("CT6", 1e-10, "below_one", {"value": None, "reason": "undetermined"}))
for rid, c, convention, neglected in CT_C:
    row(rid, "CT", "api_spectrum",
        pub("Al", 1, "D1000", "E_M", weight_length_kind=W13, weight_length_A=None, extra_inverse_length_per_A=c),
        ok("truncated", convention_ratio=convention, neglected_total_area=neglected), ci=False, denominator="state",
        checks=["the convention field and the computed interval are reported separately"])
row("CT7", "CT", "analytic_areas",
    {"kernel": "K_A2", "model": MODEL_ID, "max_order": 3,
     "window_override_eV": [num(1.3 * (1 + 2.0 ** -8)), num(7.7)],
     "note": "the kernel formula is evaluated on the shifted window"},
    ok(), ci=True, denominator="numeric", reuse=["A2.result"],
    checks=["h_n of A2.result agrees with (w_n(a) - w_n(a'))/ln(a'/a) within 2 percent for n = 2, 3 (11 columns "
            "with non-zero value)"],
    plans=[areas_plan((1, 3), "R-A"), plan(["fd.h[2]", "fd.h[3]"], "identity", "prereg 7.4")])

# ---------------------------------------------------------------- mutations
def mut(mid, text, detect, note=None, branches=None):
    """detect: (row, [concrete keys or 'terminal_state']).  Kinds and oracles are resolved from the ledger below."""
    m = {"id": mid, "class": "required", "mutation": text,
         "detected_by": [{"row": r, "keys": list(k)} for r, k in detect], "note": note,
         "state_change": "detected when any record of detected_by is broken; an unrelated exception, a crash of "
                         "the test harness or a failure of the normal control is not a detection"}
    if branches:
        m["branches"] = branches
        m["branch_rule"] = "every listed branch is run and reported; the id is detected when all required branches fail"
    return m


TS = "terminal_state"
MUTATIONS = [
    mut("X01", "later losses use the total column", [("A3a", ["order[2].bin[9].total"]), ("RB1", ["order[2].bin[89].total"])],
        "A3a bin[9] is loss 5..6 eV; RB1 bin[89] is loss 30..31 eV (auxiliary); the reference confirms the width"),
    mut("X02", "later losses use the no_hole column", [("A3b", ["order[2].bin[9].total"])],
        "with hole charge 0 the term 15_T1 would repeat"),
    mut("X03", "the first loss uses the 13a column",
        [("M4A-c1", ["single_loss.all_legacy_fields"]), ("A3a", ["order[1].bin[12].total", "order[2].bin[9].total"])]),
    mut("X04", "the excess is multiplied by n", [("A3a", ["order[2].bin[9].total", "order[2].area.total"])]),
    mut("X05", "each order is divided by n!", [("A1", ["order[2].bin[5].total", "order[2].area.total"])],
        "the positive second order halves"),
    mut("X06", "the repeated factor uses L_w instead of m",
        [("A3d", ["order[2].bin[9].total"]), ("RB7", ["order[2].bin[89].total"])],
        "A3d: the repeated factor changes from 3 to 2"),
    mut("X07", "each order carries m only once", [("A3d", ["order[2].area.total", "order[2].bin[10].total"])],
        "A3d: m^2 becomes m; registered on keys different from X06"),
    mut("X08", "S2 is applied per component",
        [("A3e", ["order[1].area.no_hole", "order[1].area.hole", "order_area_ratio"])],
        "the unmultiplied areas are 0.18 and 0.01 per angstrom; an exception raised by a per-component reciprocal is "
        "not a detection unless it breaks these keys"),
    mut("X09", "second order is set to zero for s > b", [("A1s", ["order[2].bin[2].total"])]),
    mut("X10", "C_{n-1}(s-t) becomes C_{n-1}(s+t)",
        [("A2", ["order[2].bin[5].total", "panels.laplace[2][0.3]", "panels.moment[2][1]"])]),
    mut("X11", "a mandatory interior knot is dropped: RN(2a+b) of order 3 in A1",
        [("A1", ["order[3].mesh.initial_panel_edges_eV[RN(2a+b)]", "order[3].mesh.panel_edges_eV[RN(2a+b)]"])],
        "the real panel edges are compared, not only the self-reported knots_eV or the numerical accuracy"),
    mut("X12", "17 Gauss points become 9",
        [("A4a", ["order_panels[2].value_at_7eV.term_13a"]), ("A4c", ["order_panels[2].value_at_7eV.term_13a"])],
        "A4a is decisive for the 9-point rule"),
    mut("X13", "one propagation term of the error is dropped",
        [("A3f", ["order[2].area_error_parts[inherited][total]",
                  "control_point_path.order[2].node_parts[inherited][total]@7eV"])],
        "two paths, one id: the area path (Young right-hand side) and the point path (node parts handed to "
        "_fold_order); each branch is run on both paths and recorded",
        branches=["area path: drop ||R||*e_1 (0.00084*m^2)", "area path: drop ||C_1||*e_R (0.00147*m^2)",
                  "area path: drop the mixed term (0.000049*m^2)", "point path: drop |C_1|*eps_R (0.0003)",
                  "point path: drop eps_1*|R| (0.0003)", "point path: drop the mixed term (0.000003)"]),
    mut("X14", "the first order of A2 is averaged over the cells [1000 - d_(i+1), 1000 - d_i] of E_A1 on the loss axis "
               "(mean over the full cell width, zero outside the window) before the second convolution",
        [("A2", ["order[2].bin[5].total"])],
        "the first-order satellite mass of A2 in both outside regions is 0; the reference records at least one "
        "second-order bin outside the width for this coarse-graining (bin[5] is the registered candidate)"),
    mut("X15", "outside.low is the balance residual and outside.high is dropped",
        [("A1", ["order[1].outside.high.total", "order[6].outside.low.total"])],
        "order 1 reaches outside.high (from the 1012 eV line) and order 6 reaches outside.low; both are compared "
        "directly, not through the balance"),
    mut("X16", "the neglected total area is added to the computed outside.low",
        [("A1", ["order[6].outside.low.total", "spectrum.outside.low"]), ("T2c", ["order[4].outside.low.total"])],
        "adding the neglected area only to the final outside.low sum leaves the per-order keys unchanged; it breaks "
        "the identity spectrum.outside.low = J0 + sum of the order totals of A1"),
    mut("X17", "the support boundary uses a strict comparison",
        [("T1a", [TS, "truncation.support_order", "truncation.state"])],
        "T1b (N_D = 2, complete_on_display) is the normal control of the pair"),
    mut("X18", "the bound loses exp(tau*U) or shifts the power of Lambda_R",
        [("T2c", ["truncation.bound.bin[0].grid_minimum", "truncation.bound.bin[1].grid_minimum"])],
        branches=["exp(tau*U) dropped", "power N -> N+1", "power N -> N-1"]),
    mut("X19", "the sum omits order N", [("A1s", ["spectrum.bin[0]", "spectrum.bin[1]"]), ("INV1", ["spectrum.bin[9]"])],
        "bins reached by the last order; INV1 bin[9] (loss 5..6 eV) is on the reaching side"),
    mut("X20", "negative values are clipped to 0",
        [("A3a", ["order[1].bin[12].hole", "order[2].bin[10].hole"]), ("A3c", ["order[1].bin[12].total"])],
        "bin[12] is loss 1..2 eV and bin[10] is loss 3..4 eV, where the exact values are negative"),
    mut("X21", "the first order is recomputed on a new path or computed twice",
        [("M4A-c1", ["single_loss.all_legacy_fields", "_call_kernel.call_sequence", "on_final_panels.deliveries"]),
         ("M4A-o3", ["on_final_panels.received_panels_are_order2_input"])],
        "M4A-c1 (max_order = 1) checks the old values, the call sequence and the delivery count; the hand-over to "
        "order 2 is observed in M4A-o3 by capturing the objects",
        branches=["new arithmetic path that copies the values", "the old path run twice"]),
    mut("X22", "rho_n uses the exponent 1 for every order",
        [("A3f", ["order[2].area_error.total", "order[3].area_error.total"])]),
    mut("X23", "(n-1) is dropped from h_n", [("A2", ["order[3].cutoff_sensitivity.total"]), ("CT7", ["fd.h[3]"])]),
    mut("X24", "the interval-corrected error is fed back into the recursion",
        [("A3f", ["order[2].area_error_parts[inherited][total]", "order[3].area_error_parts[inherited][total]",
                  "control_external.order[2].area_error_parts[inherited][total]",
                  "control_external.order[3].area_error_parts[inherited][total]"])],
        branches=["feedback 1 -> 2", "feedback 2 -> 3"]),
    mut("X28", "J0 is added to the formal higher-order total (reading the spectrum key of the fold, or an equivalent)",
        [("A1s", ["order[2].bin[6].total", "spectrum.bin[6]"])],
        "bin[6] is the zero-loss bin: order 2 is exactly 0 and the sum is 1; an access to a missing key that raises "
        "is not a detection"),
    mut("X29", "Lambda is built from point values only",
        [("T2e", ["truncation.bound.bin[0].grid_minimum_with_error_terms"])],
        "a violated panel precondition makes the test not_run"),
    mut("X30", "a local error is returned as inherited and no bisection is made",
        [("OI1", [TS, "order[2].mesh.stages", "order[2].bin[1].total_error_parts[local]",
                  "order[2].bin[1].total_error_parts[inherited]"])],
        "observed either as inherited_error_exceeds_order_budget at order 2 or as success with stages = 0"),
    mut("X31", "bounded_on_display includes the outside regions, or drops the lowest-KE bin",
        [("T2a", [TS]), ("T2b", [TS])],
        branches=["outside regions included (T2a fails: even the tau -> 0 limit 1.1956*0.7056^4/(1-0.7056) exceeds 0.05)",
                  "bin[0] dropped (T2b succeeds)"]),
    mut("X33", "the absolute-value envelope is integrated by the Gauss rule directly on the unsplit panel",
        [("AE1", ["order[2].area_error_parts[inherited][hole]"])],
        "the first-order panels are fixed by initial_edges = [2, 4, 5]"),
]
DEFERRED_MUTATIONS = {"X25": "M4b-2 (binding-energy axis)", "X26": "M4b-2", "X27": "M4b-2", "X32": "M4b-2"}

# ---------------------------------------------------------------- self checks
ids = [r["id"] for r in rows]
need(len(ids) == len(set(ids)), "duplicate row id")
by_id = {r["id"]: r for r in rows}


DESCRIPTIVE_FIELDS = ("note",)


def typed(v):
    if isinstance(v, bool):
        return ["bool", v]
    if isinstance(v, int):
        return ["int", v]
    if isinstance(v, float):
        return ["float", v.hex()]
    if v is None:
        return ["none"]
    if isinstance(v, str):
        return ["str", v]
    if isinstance(v, dict):
        return ["dict", {k: typed(x) for k, x in v.items() if k not in DESCRIPTIVE_FIELDS}]
    return ["list", [typed(x) for x in v]]


def signature(entry, inputs, settings):
    """Typed signature of a main entry after expanding every default of prereg 4.1 (an omitted argument equals
    its explicit default, also for multiplier and multiplier_A); bool and int stay distinct; descriptive fields
    are dropped; CLI rows keep their raw argv."""
    inputs, settings = dict(inputs), dict(settings)
    if entry in ("api_spectrum", "api_areas", "analytic_spectrum", "analytic_areas", "reference_model"):
        base = DEFAULTS["public" if entry.startswith("api") else "analytic"]
        for k, v in base.items():
            inputs.setdefault(k, v)
        for k, v in {**DEFAULTS["first_order_settings"], **DEFAULTS["order_settings"]}.items():
            settings.setdefault(k, v)
    return json.dumps([entry, typed(inputs), typed(settings)], sort_keys=True)


def find_duplicate(table):
    seen_sig = {}
    for r in table:
        k = signature(r["entry"], r["inputs"], r["settings"])
        if k in seen_sig:
            return seen_sig[k], r["id"]
        seen_sig[k] = r["id"]
    return None


def check_duplicates(table):
    need(find_duplicate(table) is None, f"rows with the same effective entry signature: {find_duplicate(table)}")


def must_raise(action, label):
    try:
        action()
    except TableError:
        return
    raise TableError(f"self-test did not raise: {label}")


check_duplicates(rows)
# self tests of the duplicate check: variants that must, and must not, collide
_o3 = by_id["M4A-o3"]
need(find_duplicate([_o3, {**_o3, "id": "x", "settings": {"epsrel": 1e-6}}]) is not None,
     "an explicit default setting must collide with the omitted one")
need(find_duplicate([_o3, {**_o3, "id": "x", "inputs": {**_o3["inputs"], "hole_charge": 1.0, "note": "text"}}])
     is not None, "an explicit default input and a descriptive field must collide")
need(find_duplicate([by_id["M4A-c1"], by_id["M4A-c2"]]) is None, "c1 and the real c2 (external, 20 A) differ")
_c2_default = {**by_id["M4A-c2"], "id": "x",
               "inputs": {**by_id["M4A-c2"]["inputs"], "multiplier": "weight_length", "multiplier_A": None}}
must_raise(lambda: check_duplicates([by_id["M4A-c1"], _c2_default]),
           "c2 with the explicit default multiplier must collide with c1")
need(find_duplicate([_o3, {**_o3, "id": "x", "settings": {"max_convolution_points": True}}]) is None,
     "bool and int settings differ")

# expansion of the ledger: no duplicate (key, kind, oracle) inside a row, and every key inside the range of
# orders and bins given by the row inputs (or by the in-row control that the key names)
ORDER_PATTERNS = (re.compile(r"order\[(\d+)(?:\.\.(\d+))?\]"), re.compile(r"order_panels\[(\d+)\]()"),
                  re.compile(r"panels\.(?:laplace|moment)\[(\d+)\]()"))
BIN_PATTERN = re.compile(r"bin\[(\d+)(?:\.\.(\d+))?\]")


def key_scope(r, key):
    """(max order, number of bins or None) that the key may address."""
    inputs = r["inputs"]
    for name, control in inputs.get("controls", {}).items():
        if key.startswith(name + "."):
            edges = control.get("output_edges_eV", inputs.get("output_edges_eV"))
            return control["max_order"], (nbins(edges) if edges else None)
    order = inputs.get("max_order")
    if order == "display_support":
        order = r["expect"].get("computed")
    if isinstance(order, bool) or not isinstance(order, int):
        order = None
    edges = inputs.get("output_edges_eV")
    return order, (nbins(edges) if isinstance(edges, str) and edges in EDGES else None)


def check_key_range(r, key):
    max_order, bins = key_scope(r, key)
    for pattern in ORDER_PATTERNS:
        for m in pattern.finditer(key):
            hi = int(m.group(2) or m.group(1))
            need(max_order is not None, f"{r['id']}: {key!r} names an order but the row has no integer max_order")
            need(1 <= int(m.group(1)) <= hi <= max_order,
                 f"{r['id']}: {key!r} is outside the orders 1..{max_order} of the row")
    for m in BIN_PATTERN.finditer(key):
        hi = int(m.group(2) or m.group(1))
        need(bins is not None, f"{r['id']}: {key!r} names a bin but the row has no output edges")
        need(int(m.group(1)) <= hi < bins, f"{r['id']}: {key!r} is outside the bins 0..{bins - 1} of the row")


def ledger_of(r):
    records = {}
    for p in r["check_plan"]:
        for key in p["_concrete"]:
            rec = (key, p["kind"], p["oracle"])
            need(rec not in records, f"{r['id']}: duplicate ledger record {rec}")
            check_key_range(r, key)
            records[rec] = True
    return list(records)


LEDGER = {r["id"]: ledger_of(r) for r in rows}


def bind(mutation, ledger, table):
    """Resolve the keys of a mutation to ledger records; the range of every target is checked again."""
    bound = []
    for d in mutation["detected_by"]:
        need(d["row"] in table, f"{mutation['id']}: unknown row {d['row']}")
        need(d["keys"], f"{mutation['id']}: no key for row {d['row']}")
        need(len(set(d["keys"])) == len(d["keys"]), f"{mutation['id']}: duplicate key for row {d['row']}")
        resolved = []
        for key in d["keys"]:
            check_key_range(table[d["row"]], key)
            hits = [{"key": k, "kind": kind, "oracle": oracle} for k, kind, oracle in ledger[d["row"]] if k == key]
            need(hits, f"{mutation['id']}: {d['row']} has no ledger record for {key!r} (wrong order, region or name)")
            resolved += hits
        bound.append({"row": d["row"], "records": resolved})
    need(len({(d["row"], rec["key"], rec["kind"]) for d in bound for rec in d["records"]})
         == sum(len(d["records"]) for d in bound), f"{mutation['id']}: duplicate record")
    return bound


# self-tests of the range checks (counterexamples of the third Codex review); each must raise TableError
_a1s = by_id["A1s"]
must_raise(lambda: ledger_of({**_a1s, "inputs": {**_a1s["inputs"], "max_order": 2}}),
           "A1s with N = 2 that keeps its order-3 keys")
_a1s_99 = {**_a1s, "check_plan": _a1s["check_plan"] + [plan(["order[99].bin[2].total"], "direct_reference", "R-A")]}
must_raise(lambda: ledger_of(_a1s_99), "an order-99 key in the ledger of A1s")


def _x09_with_order_99():
    ledger = {"A1s": [(key, p["kind"], p["oracle"]) for p in _a1s_99["check_plan"] for key in p["_concrete"]]}
    bind({"id": "X09", "detected_by": [{"row": "A1s", "keys": ["order[99].bin[2].total"]}]}, ledger, {"A1s": _a1s_99})


must_raise(_x09_with_order_99, "an order-99 key present in both the ledger and X09")
must_raise(lambda: ledger_of({**_a1s, "check_plan": [plan(["order[2].bin[7].total"], "direct_reference", "R-A")]}),
           "a bin index beyond the output edges of A1s")
PUBLIC_REASONS = INPUT_REASONS + CONVERGENCE_REASONS
seen = {r["expect"].get("reason") for r in rows}
for reason in PUBLIC_REASONS + INHERITED_REASONS:
    need(reason in seen, f"no row for reason {reason}")
for r in rows:
    reason = r["expect"].get("reason")
    need(reason is None or reason in INPUT_REASONS + CONVERGENCE_REASONS + INHERITED_REASONS,
         f"{r['id']}: unknown reason {reason}")
    need(reason not in DEFERRED_REASONS, f"{r['id']}: deferred reason")
    for target in r["reuse"]:
        other = by_id.get(target.split(".")[0])
        need(other is not None, f"{r['id']}: reuse target {target} missing")
        need(target.split(".", 1)[1] in other["provides"], f"{r['id']}: {target} is not provided by its row")
        need(other["ci"] or not r["ci"], f"{r['id']}: CI row reuses the full-only row {target}")
        need(ids.index(other["id"]) < ids.index(r["id"]), f"{r['id']}: reuse must point to an earlier row")
    inp = r["inputs"]
    for field, table in (("kernel", KERNELS), ("primary", PRIMARIES), ("output_edges_eV", EDGES),
                         ("material", MATERIALS)):
        if isinstance(inp.get(field), str):
            need(inp[field] in table, f"{r['id']}: unknown {field} {inp[field]}")
EXPECTED_MUTATIONS = ([f"X{i:02d}" for i in range(1, 25)] + ["X28", "X29", "X30", "X31", "X33"])
need([m["id"] for m in MUTATIONS] == EXPECTED_MUTATIONS, "mutation ids")
for m in MUTATIONS:
    m["detected_by"] = bind(m, LEDGER, by_id)
REGISTERED_COUNTS = {"rows": 109, "ci_rows": 90, "full_only_rows": 19}
REDUCTION_STEPS = [(90, 19), (89, 20), (88, 21), (87, 22)]
need(len(rows) == REGISTERED_COUNTS["rows"], f"row count {len(rows)} is not the registered one")
need(sum(r["ci"] for r in rows) == REGISTERED_COUNTS["ci_rows"], "CI row count is not the registered one")
need(sum(not r["ci"] for r in rows) == REGISTERED_COUNTS["full_only_rows"], "full-only row count")
_moved = set()
for step, name in enumerate(("M4A-o3s", "C02", "RB1"), start=1):
    need(by_id[name]["ci"], f"{name} of ci_reduction_order is not a CI row")
    _moved.add(name)
    ci_now = sum(r["ci"] and r["id"] not in _moved for r in rows)
    need((ci_now, len(rows) - ci_now) == REDUCTION_STEPS[step], f"CI/full counts after moving {name}")
    for r in rows:
        if r["ci"] and r["id"] not in _moved:
            need(all(t.split(".")[0] not in _moved and by_id[t.split(".")[0]]["ci"] for t in r["reuse"]),
                 f"{r['id']} would reuse a full-only row after moving {name}")
for r in rows:
    if r["expect"]["status"] == "success" and r["group"] in ("RM-A3", "T") and not any(
            p["kind"] in EXTERNAL_REFERENCE_KINDS + ("recomputed_from_report", "inequality", "exact_zero_or_null",
                                                     "diagnostic") for p in r["check_plan"]):
        need(r["denominator"] == "state", f"{r['id']}: a successful row without any planned key")

group_counts = {g: sum(r["group"] == g for r in rows) for g in GROUPS}
ci_counts = {g: sum(r["group"] == g and r["ci"] for r in rows) for g in GROUPS}
state_counts = {}
for r in rows:
    state_counts[r["expect"]["status"]] = state_counts.get(r["expect"]["status"], 0) + 1
kind_counts = {k: 0 for k in CHECK_KINDS}
numeric_by_group = {g: 0 for g in GROUPS}
for r in rows:
    for p in r["check_plan"]:
        kind_counts[p["kind"]] += p["count"]
        if p["kind"] in EXTERNAL_REFERENCE_KINDS:
            numeric_by_group[r["group"]] += p["count"]

for r in rows:
    r["check_plan"] = [{k: v for k, v in p.items() if k != "_concrete"} for p in r["check_plan"]]

doc = {
    "schema": "xpsloss-m4b1-expected-states-v1-draft4",
    "prereg": "docs/M4B1_PREREG.md",
    "design": "docs/M4B_DESIGN.md (5th edition), scope M4b-1 (U11 (a))",
    "status": "fourth draft before the freeze (after the third Codex review)",
    "reference_record_example": {
        "row": "A3a", "key": "order[2].bin[7].total", "kind": "direct_reference", "oracle": "R-A",
        "value": "<decimal string, 30 digits or more>", "u_r": "<decimal string>",
        "interval": None, "width": {"d_abs": "<decimal string>", "d_rel": "<decimal string>"},
        "no_value_reason": None,
        "note": "format sample without values; every (row, key, kind, oracle) is enumerable; accuracy and error "
                "coverage read the same record once and are stored as two verdicts"},
    "baseline": BASELINE,
    "model_id": MODEL_ID,
    "reference_only_models": REFERENCE_MODELS,
    "columns": COLUMNS,
    "components": COMPONENTS,
    "defaults": DEFAULTS,
    "resource_formulas": RESOURCE_FORMULAS,
    "vocabulary": {"input_reasons": INPUT_REASONS, "convergence_reasons": CONVERGENCE_REASONS,
                   "inherited_reasons_with_rows": INHERITED_REASONS, "deferred_to_m4b2": DEFERRED_REASONS,
                   "truncation_states": TRUNCATION_STATES, "neglected_total_area_reasons": NEGLECTED_REASONS,
                   "bound_reasons": BOUND_REASONS, "external_reference_kinds": list(EXTERNAL_REFERENCE_KINDS),
                   "convention_ratios": CONVENTION_RATIOS, "entries": ENTRIES, "groups": GROUPS,
                   "denominators": DENOMINATORS, "check_kinds": CHECK_KINDS, "statuses": STATUSES},
    "key_names": {"regions": "outside.low, bin[0..n-1], outside.high",
                  "order": "order[n] is the element of result.orders whose .order equals n (n >= 1)",
                  "components": "total, no_hole, hole", "areas": "order[n].area.<column> = w_{n,j}",
                  "sensitivity": "order[n].cutoff_sensitivity.<column> = h_{n,j}",
                  "error_parts": "order[n].<component>_error_parts[part][i], order[n].outside[side]"
                                 "['<component>_error_parts'][part], order[n].area_error_parts[part][column]; "
                                 "parts: inherited, propagated, local (fixed multiplier), multiplier_interval and, "
                                 "for areas only, convention (non-zero only for order 1, total, S2)",
                  "expansion": "key groups are expanded to concrete keys at generation time and checked; the file "
                               "stores a summary and the count; free-text groups are enumerated by index"},
    "kernel_evaluation_counters": {
        "public_entries": "calls of xpsloss.xps_single_loss._call_kernel",
        "analytic_entries": "calls of the kernel(T, side) passed to the entry",
        "layer_c": "no kernel is called"},
    "materials": MATERIALS,
    "kernels": KERNELS,
    "primaries": PRIMARIES,
    "output_edges": EDGES,
    "limits_seconds": {"ci_increment": 60, "full": 1800, "row": 120},
    "ci_reduction_order": ["M4A-o3s", "C02", "RB1"],
    "ci_reduction_counts": [list(x) for x in REDUCTION_STEPS],
    "ci_reduction_note": "A5 and CT3 are full-only from the start; a moved row stays in the full run and in the "
                         "mutation runs; DS1 keeps a light display_support success in CI if C02 is moved",
    "execution_order": ids,
    "counts": {"rows": len(rows), "ci_rows": sum(r["ci"] for r in rows),
               "full_only_rows": sum(not r["ci"] for r in rows), "by_group": group_counts,
               "ci_by_group": ci_counts, "by_expected_status": state_counts,
               "planned_keys_by_kind": kind_counts, "planned_numeric_keys_by_group": numeric_by_group,
               "planned_external_reference_keys": sum(numeric_by_group.values()),
               "note": "planned_numeric_keys_by_group counts the keys that read an external reference (kinds "
                       "direct_reference, bound_reference, error_part, interval_contains); different keys may read "
                       "the same reference quantity, so this is not a number of independent conditions",
               "registered": REGISTERED_COUNTS,
               "mutations_required": len(MUTATIONS)},
    "rows": rows,
    "mutations": MUTATIONS,
    "deferred_mutations": DEFERRED_MUTATIONS,
}
text = json.dumps(doc, ensure_ascii=True, indent=1, sort_keys=True) + "\n"
data = text.encode("ascii")
if sys.argv[1:] == ["--check"]:
    same = OUT.is_file() and OUT.read_bytes() == data
    print(json.dumps({"identical": same, "sha256": hashlib.sha256(data).hexdigest()}))
    sys.exit(0 if same else 1)
OUT.write_bytes(data)
print(json.dumps({"sha256": hashlib.sha256(data).hexdigest(), **doc["counts"]}, indent=1))
