"""M4a pre-registration: the 100-row expected-state table (docs/M4A_PREREG.md section 5).

Writes expected_states_m4a.json next to this file (ASCII, LF) and prints its SHA-256 and the counts.
One row is one new entry point call (entry, inputs, settings) with its own expected terminal state.
"in_row_checks" are control runs and comparisons made inside a row; "reuse" names the artifacts of
other rows that are really read ("ROW.artifact").  Neither adds a row.  Every row carries a
"check_plan": the ledger of planned keys as strings "key|check_kind|oracle", expanded now, before
any reference value exists.  The script checks, with explicit exceptions (not assert), that ids and
normalized entries (entry, inputs, settings) are unique, that there are exactly 100 rows with the
registered split by family and by CI/full, that every reason is in the registered vocabulary, that
reuse targets exist and form no cycle, that a CI row never reuses a full-only row, that the planned
key counts are the registered ones (CE 66, RB 1048), and that the mutation ids and their detecting
rows and keys exist.  It imports only the standard library (not xpsloss, not NumPy).

Usage: python benchmarks/numerics/m4a/expected_states_m4a.py          (write the JSON)
       python benchmarks/numerics/m4a/expected_states_m4a.py --check  (compare byte for byte)
"""
import hashlib
import json
import math
import pathlib
import struct
import sys

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "expected_states_m4a.json"
INPUT_AL = "benchmarks/inputs/al_elfdb_v1.0.bulk.json"
INPUT_SI = "benchmarks/inputs/si_elfdb_v1.0.bulk.json"


class TableError(RuntimeError):
    pass


def need(cond, message):
    if not cond:
        raise TableError(message)


if sys.argv[1:] not in ([], ["--check"]):
    raise SystemExit("usage: expected_states_m4a.py [--check]")

# The script lives in work/consult/m4a/prereg/ before the freeze and in benchmarks/numerics/m4a/ after it.
ROOT = next((p for p in HERE.parents if (p / INPUT_AL).is_file()), None)
need(ROOT is not None, "repository root with benchmarks/inputs not found")


def hx(x):
    return float(x).hex()


def num(x):
    """A binary64 input as decimal and hex; non-finite values are strings."""
    x = float(x)
    if not math.isfinite(x):
        return {"invalid": repr(x)}
    return {"value": x, "hex": x.hex()}


def nums(seq):
    return {"values": [float(x) for x in seq], "hex": [hx(x) for x in seq]}


# ---------------------------------------------------------------- materials
def _read(rel):
    data = (ROOT / rel).read_bytes()
    return json.loads(data.decode("utf-8")), hashlib.sha256(data).hexdigest()


AL_CFG, AL_SHA = _read(INPUT_AL)
SI_CFG, SI_SHA = _read(INPUT_SI)


def _lines(cfg):
    return [[o["energy_eV"], o["strength_eV2"], o["width_eV"], o["dispersion"]]
            for o in cfg["material"]["oscillators"]]


AL = _lines(AL_CFG)
SI = _lines(SI_CFG)
G_SI = float(SI_CFG["material"]["gap_eV"])
need(AL_CFG["material"]["gap_eV"] == 0 and G_SI == 1.12 and AL == [[15.6, 243.3598, 1.5, 1.0]],
     "unexpected Al or Si input")


def mat(gap, lines, origin, synthetic=True):
    lines = [[float(v) for v in line] for line in lines]
    return {"gap_eV": float(gap), "gap_hex": hx(gap), "oscillators": lines,
            "oscillators_hex": [[hx(v) for v in line] for line in lines],
            "columns": ["energy_eV", "strength_eV2", "width_eV", "dispersion"],
            "synthetic": synthetic, "from": origin}


MATERIALS = {
    "Al": {**mat(0.0, AL, INPUT_AL, synthetic=False), "file": INPUT_AL, "sha256": AL_SHA},
    "Si": {**mat(G_SI, SI, INPUT_SI, synthetic=False), "file": INPUT_SI, "sha256": SI_SHA},
    "Al0": mat(0.0, [[e, a, w, 0.0] for e, a, w, _ in AL], "Al, dispersion 0 (CE4)"),
    "SiZ": mat(G_SI, [[e, 0.0, w, d] for e, _, w, d in SI], "Si, strengths 0 (row Z1)"),
    "synthetic_wide": mat(0.0, [[20, 200, 8, 0]], "design 6.1, cheap CI material"),
    "N1": mat(2.0, [[3, 6, 0.02, 0.5]], "M2-2 and M3a pre-registrations"),
    "N2": mat(2.0, [[3, 6, 0.2, 0.5]], "N1 width x10; listed for the panel count only, no row uses it"),
    "CE1": mat(0.0, [[3, 3450, 0.1, 96]], "design 6.3 (review counterexample)"),
    "CE2": mat(0.0, [[1, 0.5, 0.1, 3037.62221442588]], "design 6.3 (review counterexample)"),
    "CE3": mat(0.0, [[5, 0.5, 0.1, 1.21604938271605]], "design 6.3 (review counterexample)"),
    "CE5": mat(0.0, [[10, 30, 1, 0], [10.05, 30, 1, 0]], "design 6.3 (two dispersionless lines)"),
    "CE6": mat(0.0, [[10, 30, 1, 0], [10.0000001, 30, 1, 0]], "design 6.3 (nearly double line)"),
    "P1": mat(2.0, [[3, 8.1, 0.5, 0.5]], "M2-2 pre-registration (V07 failure control)"),
}

# ---------------------------------------------------------------- defaults and vocabularies
AREA_SETTINGS = ["epsrel", "area_epsabs", "panel_limit", "kernel_epsabs", "kernel_epsrel", "kernel_limit"]
SPECTRUM_SETTINGS = ["epsrel", "area_epsabs", "spectrum_epsabs_fraction", "panel_limit", "max_cells",
                     "kernel_epsabs", "kernel_epsrel", "kernel_limit"]
DEFAULTS = {
    "electron_energy_eV": 1000.0, "weight_length_A": 20.0, "hole_charge": 1.0, "epsilon_momentum": "k_total",
    "theta_deg": 0.0,
    "area_settings": AREA_SETTINGS, "spectrum_settings": SPECTRUM_SETTINGS,
    "settings": {"epsrel": 1e-6, "area_epsabs": 1e-10, "spectrum_epsabs_fraction": 1e-9, "panel_limit": 400,
                 "max_cells": 2000000, "kernel_epsabs": 1e-12, "kernel_epsrel": 1e-9, "kernel_limit": 200},
}
COLUMNS = ["term_13a", "term_13b", "term_14", "term_15_T1", "term_15_T2", "term_15_T3", "term_16_e",
           "term_16_h", "no_hole", "hole", "total"]
INPUT_REASONS = ["invalid_window", "window_not_above_gap", "invalid_primary", "invalid_output_edges",
                 "invalid_settings", "invalid_schema", "cell_limit_exceeded", "resolution_exceeds_panel_limit",
                 "panel_width_unresolved"]
CONVERGENCE_REASONS = ["kernel_not_converged", "kernel_error_exceeds_budget", "loss_budget_not_met",
                       "cell_limit_reached", "panel_width_unresolved", "nonfinite_output"]
DOMINANT_PARTS = ["interpolation", "verification", "kernel_point_error", "edge", "rounding"]
BASES = ["fixed_poles_enumerated", "endpoint_roots_enumerated", "computed_candidates_only"]
ENTRIES = ["layer_ab", "layer_a", "layer_b", "api_area", "api_spectrum", "api_spectrum_pointwise",
           "cli_area", "cli_spectrum", "existing_api", "existing_cli"]
DENOMINATORS = ["numeric", "state", "numeric_relative_and_diagnostics", "replay_and_state", "invariance",
                "none"]

CHECK_KINDS = ["state", "direct_reference", "inverse_interval", "identity", "exact_zero_or_null", "replay",
               "error_part", "diagnostic", "record"]
AREA_KEYS = ([f"window_area.{c}" for c in COLUMNS] + ["single_loss_weight"]
             + ["inverse_window_area_interval_A.lower", "inverse_window_area_interval_A.upper"]
             + [f"cutoff_sensitivity.{c}" for c in COLUMNS])
QUANTITIES = ["zero_loss", "single_loss_total", "single_loss_no_hole", "single_loss_hole", "spectrum"]
SATELLITES = QUANTITIES[1:4]
SENSITIVITIES = [f"cutoff_sensitivity.{side}.{c}" for side in ("plus", "minus")
                 for c in ("total", "no_hole", "hole")]
need(len(AREA_KEYS) == 25, "area key count")


def regions(nbins):
    return ["outside.low"] + [f"bin[{i}]" for i in range(nbins)] + ["outside.high"]


def spectrum_keys(nbins, fields=None):
    return [f"{r}.{q}" for r in regions(nbins) for q in (fields or QUANTITIES + SENSITIVITIES)]


def plan(keys, kind, oracle):
    need(kind in CHECK_KINDS and "|" not in oracle, f"bad plan entry {kind} {oracle}")
    return [f"{k}|{kind}|{oracle}" for k in keys]


def direct(keys, oracle):
    return [f"{k}|{'inverse_interval' if k.startswith('inverse_') else 'direct_reference'}|{oracle}"
            for k in keys]


def _sha_f8(values):
    return hashlib.sha256(struct.pack(f"<{len(values)}d", *values)).hexdigest()


def _linspace(start, stop, num):
    step = (stop - start) / (num - 1)
    out = [i * step + start for i in range(num)]
    out[-1] = stop
    return out


GA_FULL = [0.2 * k for k in range(1, 500)]
GC_FULL = _linspace(0.2, 99.8, 499)
need(_sha_f8(GA_FULL) == "6214a4e35ef20fa9f1c9f5f19a376d6888847f081e5942672eba0ae9b9abf33e"
     and _sha_f8(GC_FULL) == "e69297cdd297d6432c4976f071938e792fb93f29f8681cba47da75677a842524",
     "GA or GC differs from the arrays of the M3a pre-registration")
GA_CI = GA_FULL[::4]
GC_CI = _linspace(0.2, 99.4, 125)
V_PROFILES = {
    "ci": {"GA": {"expression": "(0.2*np.arange(1, 500))[::4]", "length": len(GA_CI),
                  "sha256_float64_le": _sha_f8(GA_CI)},
           "GC": {"config_loss_grid_eV": {"min": 0.2, "max": 99.4, "step": 0.8}, "length": len(GC_CI),
                  "sha256_float64_le": _sha_f8(GC_CI)}},
    "full": {"GA": {"expression": "0.2*np.arange(1, 500)", "length": len(GA_FULL),
                    "sha256_float64_le": _sha_f8(GA_FULL)},
             "GC": {"config_loss_grid_eV": {"min": 0.2, "max": 99.8, "step": 0.2}, "length": len(GC_FULL),
                    "sha256_float64_le": _sha_f8(GC_FULL)}},
}
TAG = {"name": "v0.2.0", "commit": "023820cfbc95547a1cc90c11d1f2bcdf19400806"}
CONFIG_INTEGRATION = {"epsabs": 1e-12, "epsrel": 1e-9, "limit": 200}

H2 = math.ulp(2.0)
A_N08 = float.fromhex("0x1.1eb851eb85ab8p+0")
A_N09 = float.fromhex("0x1.1eb851eb8751bp+0")
need(A_N08 == G_SI + 5e-13 and A_N09 == G_SI + 2e-12 and H2 == 2.0 ** -51, "unexpected binary64 inputs")

PRIMARIES = {
    "D1000": {"kind": "delta", "energies_eV": nums([1000]), "intensities": nums([1])},
    "D1000x2": {"kind": "delta", "energies_eV": nums([1000]), "intensities": nums([2])},
    "D1000zero": {"kind": "delta", "energies_eV": nums([1000]), "intensities": nums([0])},
    "D998": {"kind": "delta", "energies_eV": nums([998]), "intensities": nums([1])},
    "D998_1000": {"kind": "delta", "energies_eV": nums([998, 1000]), "intensities": nums([1, 2])},
    "D1008": {"kind": "delta", "energies_eV": nums([1008]), "intensities": nums([1])},
    "D1200": {"kind": "delta", "energies_eV": nums([1200]), "intensities": nums([1])},
    "D3big": {"kind": "delta", "energies_eV": nums([980, 1000, 1030]), "intensities": nums([1, 2.0 ** 60, 2])},
    "B_RB3": {"kind": "binned", "edges_eV": nums([998, 999, 1001]), "intensities": nums([0.3, 0.7])},
    "B_RA3": {"kind": "binned", "edges_eV": nums([998, 999, 1002]), "intensities": nums([2, 1])},
}
E_RB1 = [960, 980, 984, 988, 995, 999, 1000, 1001]
EDGES = {
    "E_RB1": nums(E_RB1),
    "E_RB1_plus8": nums([e + 8 for e in E_RB1]),
    "E_RB1_plus200": nums([e + 200 for e in E_RB1]),
    "E_RB3": nums([960, 980, 984, 990, 998, 1000]),
    "E_RA1": nums([990, 994, 995, 996.5, 998, 999, 1000, 1001]),
    "E_RA2": nums([996, 997, 1000]),
    "E_RA3": nums([995, 996, 998, 1000]),
    "E_RA4": nums([997, 998.5, 999.5, 1000, 1001]),
    "E_RA5": nums([995, 996.85, 996.87, 999]),
    "E_RA6": nums([990, 1002, 1029]),
    "E_RA6b": nums([990, 1002]),
    "E_RA7": nums([994, 995, 998, 1001]),
    "E_RA7b": nums([0, 0.1, 0.2]),
    "E_RA9": nums([997.98, 998]),
    "E_RA10": nums([990, 1001]),
    "E_RA12": nums([997, 998]),
    "E_Z1": nums([990, 999, 1000, 1001]),
    "E_L": nums([999, 999.8, 999.95, 999.99, 1000, 1001]),
    "E_C04": {"grid": {"min": 960, "max": 1000, "step": 8}, **nums([960, 968, 976, 984, 992, 1000])},
}

rows = []


def window(a, b):
    return {"loss_window_eV": [num(a), num(b)]}


def ok(**attributes):
    return {"status": "success", **attributes}


def rejected(reason, exception="SingleLossInputError", **attributes):
    return {"status": "error", "exception": exception, "reason": reason, **attributes}


def unconverged(reason, **attributes):
    return {"status": "not_converged", "exception": "SingleLossConvergenceError", "reason": reason, **attributes}


def nb(edges):
    return len(EDGES[edges]["values"]) - 1


def row(rid, family, entry, inputs, expect, *, ci, denominator, settings=None, reuse=(), checks=(),
        reference=None, note=None, keys=(), auxiliary=None):
    rows.append({"id": rid, "family": family, "entry": entry, "inputs": inputs, "settings": dict(settings or {}),
                 "expect": expect, "ci": bool(ci), "reuse": list(reuse), "in_row_checks": list(checks),
                 "denominator": denominator, "reference": reference, "note": note,
                 "check_plan": ["terminal_state|state|expect"] + list(keys),
                 "auxiliary_controls": auxiliary})


# ---------------------------------------------------------------- R-A: analytic kernels, package not imported by the reference
def analytic(kernel, a, b, *, length=2.0, features, primary=None, edges=None, **extra):
    d = {"kernel": kernel, **window(a, b), "weight_length_A": num(length), "declared_features": features}
    if primary is not None:
        d["primary"] = primary
        d["output_edges_eV"] = edges
    d.update(extra)
    return d


RECT25 = {"form": "rectangle", "support_eV": [2, 5], "term_13a": 0.10, "other_terms": 0.0}
JUMPS25 = {"jumps_eV": [2, 5]}
row("RA1", "RA", "layer_ab",
    analytic({"form": "rectangle", "support_eV": [2, 5],
              "terms": {"term_13a": 0.04, "term_13b": 0.03, "term_15_T1": 0.02, "term_15_T3": -0.01,
                        "term_16_e": 0.005, "term_14": 0.02, "term_15_T2": -0.007, "term_16_h": 0.002}},
             1, 8, features=JUMPS25, primary="D1000", edges="E_RA1"),
    ok(), ci=True, denominator="numeric", reference="R-A",
    keys=direct(AREA_KEYS + spectrum_keys(nb("E_RA1")), "R-A")
    + plan([f"{r}.{q}@multiplier3" for r in regions(nb("E_RA1")) for q in SATELLITES], "identity",
           "abs(J3 - 1.5*J2) <= 4*eps*abs(1.5*J2), fold on the same panels"),
    checks=["layer B re-evaluated on the same panels with the multiplier 3 instead of 2 (areas unchanged, "
            "satellite sum 0.90)",
            "satellite sum 0.60 and S*L*A_total = 0.60 are compared separately",
            "zero loss: 0 in [999,1000), 1 in [1000,1001]"])
row("RA2", "RA", "layer_ab", analytic(RECT25, 2, 5, features=JUMPS25, primary="D1000", edges="E_RA2"),
    ok(), ci=True, denominator="numeric", reference="R-A",
    keys=direct(AREA_KEYS + spectrum_keys(nb("E_RA2")), "R-A"))
row("RA3", "RA", "layer_ab",
    analytic({"form": "piecewise_linear", "support_eV": [2, 5],
              "term_13a": "0.1*max(0, min(T-2, (5-T)/2))", "term_14": "0.02*(T-3)", "other_terms": 0.0},
             1, 6, features={"jumps_eV": [2, 5], "kinks_eV": [3]}, primary="B_RA3", edges="E_RA3",
             hole_charge=num(0.5)),
    ok(), ci=True, denominator="numeric", reference="R-A",
    keys=direct(AREA_KEYS + spectrum_keys(nb("E_RA3")), "R-A"),
    note="areas 0.15 (no_hole), 0.03 (hole), 0.165 (total); hole in [998,1000] is -2/225")
RA4K = {"form": "c_over_T", "c": 0.03, "term_13a": "0.03/T", "other_terms": 0.0}
row("RA4a", "RA", "layer_ab",
    analytic(RA4K, 0.02, 2, features={"declared": "smooth on the window; low-loss grading a*2^k applies"},
             primary="D1000", edges="E_RA4"),
    ok(), ci=True, denominator="numeric", reference="R-A",
    keys=direct(AREA_KEYS + spectrum_keys(nb("E_RA4")), "R-A"))
row("RA4b", "RA", "layer_ab",
    analytic(RA4K, 0.2, 2, features={"declared": "smooth on the window; low-loss grading a*2^k applies"},
             primary="B_RA3", edges="E_RA4"),
    ok(), ci=True, denominator="numeric", reference="R-A", reuse=["RA4a.result"],
    keys=direct(AREA_KEYS + spectrum_keys(nb("E_RA4")), "R-A")
    + plan(["single_loss_weight@difference"], "identity", "w(RA4a) - w(RA4b) = 0.06*ln(10)"),
    checks=["w(RA4a) - w(RA4b) = 0.06*ln(10) (area difference 0.03*ln(10) A^-1), within the sum of the two widths"])
LOR5 = {"form": "lorentzian", "expression": "c/(1+((T-t0)/gamma)^2)", "c": 0.03, "t0_eV": num(3.137),
        "gamma_eV": num(1e-5), "column": "term_13a", "other_terms": 0.0}
row("RA5", "RA", "layer_ab",
    analytic(LOR5, 2, 5, features={"singularities": [{"center_eV": num(3.137), "half_width_eV": num(1e-5)}]},
             primary="D1000", edges="E_RA5"),
    ok(), ci=True, denominator="numeric", reference="R-A",
    keys=direct(AREA_KEYS + spectrum_keys(nb("E_RA5")), "R-A"))
row("RA6", "RA", "layer_ab", analytic(RECT25, 2, 5, features=JUMPS25, primary="D3big", edges="E_RA6"),
    ok(), ci=True, denominator="numeric", reference="R-A",
    keys=direct(AREA_KEYS + spectrum_keys(nb("E_RA6")), "R-A"),
    note="replacement widths for the small tails (prereg 7.2)")
row("RA6b", "RA", "layer_ab", analytic(RECT25, 2, 5, features=JUMPS25, primary="D3big", edges="E_RA6b"),
    ok(), ci=True, denominator="numeric", reference="R-A",
    keys=direct(AREA_KEYS + spectrum_keys(nb("E_RA6b")), "R-A"),
    note="replacement widths for the small tails (prereg 7.2)")
row("RA7", "RA", "layer_ab",
    analytic({"form": "rectangle", "support_eV": [1, 5], "term_13a": 0.10, "other_terms": 0.0}, 2, 5,
             features={"jumps_eV": [1, 5]}, primary="D1000", edges="E_RA7"),
    ok(), ci=True, denominator="numeric", reference="R-A",
    keys=direct(AREA_KEYS + spectrum_keys(nb("E_RA7")), "R-A"),
    note="H+ 0.40 in [995,998), H- 0.40 in [998,1001]")
row("RA7b", "RA", "layer_ab",
    analytic({"form": "rectangle", "support_eV": [999, 1001], "term_13a": 0.10, "other_terms": 0.0},
             999.9, 1000.5, features={"jumps_eV": [999, 1001]}, primary="D1000", edges="E_RA7b"),
    ok(), ci=True, denominator="numeric", reference="R-A",
    keys=direct(AREA_KEYS + spectrum_keys(nb("E_RA7b")), "R-A"),
    note="both one-sided sensitivities fall in [0.1,0.2): 1000 - fl(999.9) > fl(0.1)")
LOR8 = {**LOR5, "gamma_eV": num(1e-6)}
row("RA8a", "RA", "layer_a",
    analytic(LOR8, 2, 5, features={"singularities": [{"center_eV": num(3.137), "half_width_eV": num(1e-6)}]}),
    ok(attributes={"required_panels": 41}), ci=True, denominator="numeric", reference="R-A",
    keys=direct(AREA_KEYS, "R-A") + plan(["diagnostics.required_panels"], "diagnostic", "41"))
row("RA8b", "RA", "layer_a",
    analytic(LOR8, 2, 5, features={"min_half_width_eV": num(1e-6)}),
    rejected("resolution_exceeds_panel_limit", attributes={"required_panels": 1500002}, kernel_evaluations=0),
    ci=True, denominator="state",
    keys=plan(["exception.required_panels"], "diagnostic", "1500002")
    + plan(["kernel_evaluations"], "diagnostic", "0 calls of the kernel passed to integrate_window"))
row("RA9", "RA", "layer_b",
    analytic(RECT25, 2, 5, features={"panels_eV": [2, 5]}, primary="D1000", edges="E_RA9",
             panel="Panel.from_lobatto(2, 5, values, errors)",
             read_from="folded.outputs[key] with value, the five error parts and estimate",
             injection={"point_error": num(1e-6), "at": "Lobatto node with 0-based index 2 (from the left) of "
                                                         "the panel [2,5] only", "columns": "all eleven"}),
    {"status": "success", "success_definition": "layer_b_returned_on_fixed_panels",
     "driver_success_criterion_applied": False, "component_identity_applied": True},
    ci=True, denominator="numeric", reference="R-A",
    keys=plan([f"bin[0].{q}.kernel_point_error" for q in SATELLITES], "error_part",
              "R-A: abs(omega_2)*L*e, omega_2 computed independently (negative, about -0.00160240 eV)"),
    note="layer B only: no budget judgement; the point-error part of the bin [997.98, 998] is compared")
row("RA10", "RA", "layer_b",
    analytic({"form": "constant_columns", "support_eV": [2, 5], "term_13a": 0.10, "term_14": 0.02,
              "no_hole": 0.10, "hole": 0.02, "total": 0.50, "other_terms": 0.0},
             2, 5, features={"panels_eV": [2, 5]}, primary="D1000", edges="E_RA10",
             panel="Panel.from_lobatto(2, 5, values, errors) with zero errors",
             read_from="folded.outputs[key]['value']"),
    {"status": "success", "success_definition": "layer_b_returned_on_fixed_panels",
     "driver_success_criterion_applied": False, "component_identity_applied": False},
    ci=True, denominator="numeric", reference="R-A",
    keys=direct([f"bin[0].{q}" for q in SATELLITES], "R-A: bin integral of the injected column (3.0, 0.6, 0.12)"),
    note="single_loss_total comes from the injected total column (3.0), not from no_hole + Q_h*hole (0.72); "
         "the component identity total = no_hole + Q_h*hole is not applied to this row")
row("RA12", "RA", "layer_ab",
    analytic({"form": "chebyshev", "expression": "kappa*(1 - T_32(2*T-5))", "kappa": num(4e-11),
              "support_eV": [2, 3], "column": "term_13a", "other_terms": 0.0},
             2, 3, length=100.0, features={"min_variation_scale_eV": num(0.5)}, primary="D1000", edges="E_RA12"),
    ok(attributes={"refinement_stages_min": 1}), ci=True, denominator="numeric", reference="R-A",
    keys=direct(AREA_KEYS + spectrum_keys(nb("E_RA12")), "R-A")
    + plan(["diagnostics.stages"], "diagnostic", ">= 1"),
    note="fault injection that breaks the declaration contract; the area alone would pass on the first panel")
row("RA13", "RA", "layer_a",
    analytic({"form": "lorentzian", "expression": "1/(1+((T-2.25)/0.3)^2)", "column": "term_13a",
              "other_terms": 0.0}, 2, 5, features={"min_half_width_eV": num(0.3)}),
    ok(attributes={"required_panels": 6}), ci=True, denominator="numeric", reference="R-A",
    keys=direct(AREA_KEYS, "R-A") + plan(["diagnostics.required_panels"], "diagnostic", "6"),
    settings={"epsrel": 1e-13, "area_epsabs": 1e-30},
    note="compared with d_rel 1e-13 including u_r; detects a stop after the initial panels (M14a)")
row("RA14", "RA", "layer_a",
    analytic({"form": "lorentzian", "expression": "k*1/(1+((T-2.4)/0.5)^2)", "k": "1e-10/0.7753994963",
              "k_hex": hx(1e-10 / 0.7753994963), "column": "term_13a", "other_terms": 0.0},
             2, 3, features={"singularities": [{"center_eV": num(2.4), "half_width_eV": num(0.5)}]},
             injection={"interpolation_nodes": "(value 0, error K)", "verification_nodes": "(value K, error K)"}),
    unconverged("loss_budget_not_met", attributes={"dominant_part": "verification"}),
    ci=True, denominator="state", settings={"panel_limit": 1},
    keys=plan(["exception.diagnostics.error_budget.outputs[window_area.total].estimate"], "error_part",
              "R-A: about 3.5585e-10 A^-1, within 1 percent")
    + plan(["exception.dominant_part"], "diagnostic", "verification"),
    checks=["the area estimate of the total column, read from the exception snapshot, is about 3.5585e-10 A^-1 "
            "(> budget 1e-10) within 1 percent"])

# ---------------------------------------------------------------- CE: review counterexamples, public area function
CE = (("CE1", "CE1", 49.5, 50.5, "k_total", 14, "endpoint_roots_enumerated", []),
      ("CE2", "CE2", 49.5, 50.5, "k_parallel", 12, "endpoint_roots_enumerated", []),
      ("CE3", "CE3", 989.5, 990.5, "k_total", 18, None, ["window_over_E0_gt_0p1"]),
      ("CE4", "Al0", 900, 999.9, "k_total", 74, "fixed_poles_enumerated", ["window_over_E0_gt_0p1"]),
      ("CE5", "CE5", 7, 11, "k_total", 5, "fixed_poles_enumerated", []),
      ("CE6", "CE6", 7, 11, "k_total", 4, "computed_candidates_only",
       ["resolution_basis_computed_candidates_only"]))
CE_CHECKS = {
    "CE1": ["a candidate with family surface and end lower at distance 0.0222 eV (5 percent) from the window"],
    "CE2": ["a candidate with family bulk and end lower at distance 0.02536 eV (1 percent) from the window"],
    "CE3": ["a candidate with family bulk and end lower at distance 0.00503 eV (1 percent)",
            "resolution_basis is not frozen; if it is computed_candidates_only the warning must be present"],
    "CE4": ["required_panels 74 detects a missing upper grading E0-(E0-b)*2^k (67 without it)"],
    "CE5": ["two different surface candidates within 1e-3 eV (complex distance) of 8.3811-0.5i and 10.0130-0.5i"],
    "CE6": ["two different surface candidates within 1e-3 eV of 8.3516-0.5i and 9.9875-0.5i, the second with "
            "confirmed false; the bulk poles with end both; the basis is demoted and the warning is present"],
}
CE_DIAG = {"CE1": ["diagnostics.candidates[surface,lower]"], "CE2": ["diagnostics.candidates[bulk,lower]"],
           "CE3": ["diagnostics.candidates[bulk,lower]", "warnings[conditional on the basis]"],
           "CE4": [], "CE5": ["diagnostics.candidates[surface,root1]", "diagnostics.candidates[surface,root2]"],
           "CE6": ["diagnostics.candidates[surface,root1]", "diagnostics.candidates[surface,root2,unconfirmed]",
                   "diagnostics.candidates[bulk,both]"]}
for rid, m, a, b, mom, panels, basis, warns in CE:
    attrs = {"required_panels": panels}
    if basis is not None:
        attrs["resolution_basis"] = basis
    row(rid, "CE", "api_area", {"material": m, **window(a, b), "epsilon_momentum": mom},
        ok(attributes=attrs, warnings_include=warns), ci=False, denominator="numeric_relative_and_diagnostics",
        settings={"area_epsabs": 1e-20}, reference="CE (R-B type)",
        keys=direct(AREA_KEYS[:11], "CE")
        + plan(["diagnostics.required_panels", "warnings"]
               + (["diagnostics.resolution_basis"] if basis else []) + CE_DIAG[rid], "diagnostic", "prereg 5.2"),
        checks=["comparison of the eleven areas with the R-B type reference in relative width, comparison floor "
                "1e-18 A^-1 (the run uses area_epsabs 1e-20)"] + CE_CHECKS[rid])

# ---------------------------------------------------------------- R-B: real materials
SI_WARN = ["resolution_basis_computed_candidates_only"]
RB = (("RB1", "Al", 1, 50, {}, "D1000", "E_RB1", []),
      ("RB2", "Al", 1, 50, {"epsilon_momentum": "k_parallel"}, "D1000", "E_RB1", []),
      ("RB3", "Si", 5, 50, {}, "B_RB3", "E_RB3", SI_WARN),
      ("RB4", "Si", 5, 50, {"epsilon_momentum": "k_parallel"}, "B_RB3", "E_RB3", SI_WARN),
      ("RB5", "Al", 1, 50, {"weight_length_A": num(5.0)}, "D1000", "E_RB1", []),
      ("RB6", "Al", 1, 50, {"hole_charge": num(0.0)}, "D1000", "E_RB1", ["hole_charge_not_unity"]),
      ("RB7", "Al", 1, 50, {"hole_charge": num(0.5)}, "D1000", "E_RB1", ["hole_charge_not_unity"]),
      ("RB8", "Al", 1, 50, {"electron_energy_eV": num(1200.0)}, "D1200", "E_RB1_plus200", []))
for rid, m, a, b, extra, prim, edges, warns in RB:
    checks = []
    if rid == "RB1":
        checks = ["full run only: the matched display of RB1 against the completed comparison table"]
    row(rid, "RB", "api_spectrum", {"material": m, **window(a, b), **extra, "primary": prim,
                                    "output_edges_eV": edges},
        ok(warnings_include=warns), ci=(rid == "RB1"), denominator="numeric", reference="R-B", checks=checks,
        keys=direct(AREA_KEYS + spectrum_keys(nb(edges)), "R-B")
        + plan(["_call_kernel.weight_length_A", "_call_kernel.hole_charge", "_call_kernel.epsilon_momentum"],
               "diagnostic", "actual arguments passed upstream"))
for rid, m, a, b, extra, prim, edges, warns in RB[:4]:
    row(rid + "a", "RB", "api_area", {"material": m, **window(a, b), **extra},
        ok(warnings_include=warns), ci=False, denominator="numeric", reference="R-B",
        keys=direct(AREA_KEYS, "R-B"))

# ---------------------------------------------------------------- identities, zero, invariance
AL150 = {"material": "Al", **window(1, 50)}
row("I1a", "I", "api_spectrum", {**AL150, "primary": "D998", "output_edges_eV": "E_RB1"},
    ok(warnings_include=["E0_outside_primary_support"]), ci=False, denominator="numeric",
    keys=plan(["warnings"], "diagnostic", "E0_outside_primary_support"),
    note="no independent reference: its keys are not counted as verified; I1b reads its result")
row("I1b", "I", "api_spectrum", {**AL150, "primary": "D998_1000", "output_edges_eV": "E_RB1"},
    ok(), ci=False, denominator="numeric", reuse=["RB1.result", "I1a.result"],
    keys=plan(spectrum_keys(nb("E_RB1"), QUANTITIES), "identity", "y(I1a) + 2*y(RB1), prereg 7.4"),
    checks=["linearity J(F1+2*F2) = J(F1) + 2*J(F2) for every bin and both outsides (prereg 7.4)"])
row("I3", "I", "api_spectrum", {**AL150, "primary": "D1008", "output_edges_eV": "E_RB1_plus8"},
    ok(warnings_include=["E0_outside_primary_support"]), ci=False, denominator="numeric", reuse=["RB1.result"],
    keys=plan(AREA_KEYS + spectrum_keys(nb("E_RB1")), "replay", "RB1.result (values and errors, bitwise)")
    + plan(["warnings"], "diagnostic", "E0_outside_primary_support")
    + plan(["diagnostics.primary_support"], "diagnostic", "min_minus_E0_eV and max_minus_E0_eV are both +8"),
    checks=["numeric arrays and errors bitwise equal to RB1; shifted fields, difference fields and warnings "
            "change as registered"])
row("I4a", "I", "api_spectrum", {"material": "Al", **window(1, 17), "primary": "D1000",
                                 "output_edges_eV": "E_RB1"},
    ok(), ci=False, denominator="numeric",
    note="no independent reference: its keys are not counted as verified; I4b reads its result")
row("I4b", "I", "api_spectrum", {"material": "Al", **window(17, 50), "primary": "D1000",
                                 "output_edges_eV": "E_RB1"},
    ok(), ci=False, denominator="numeric", reuse=["RB1.result", "I4a.result"],
    keys=plan(AREA_KEYS[:11] + spectrum_keys(nb("E_RB1"), SATELLITES), "identity",
              "y(I4a) + y(I4b) = y(RB1), prereg 7.4")
    + plan(spectrum_keys(nb("E_RB1"), ["zero_loss"]), "replay", "RB1.result and I4a.result (bitwise)"),
    checks=["additivity of [1,17] and [17,50] against RB1 for areas, bins and outsides (prereg 7.4)"])
SW150 = {"material": "synthetic_wide", **window(1, 50)}
row("I8", "I", "api_spectrum_pointwise", {**SW150, "primary": "D1000", "output_edges_eV": "E_RB1",
                                          "internal_setting": "_single_loss_spectrum(..., kernel_batch_size=1)",
                                          "effective_settings": "those of C02 (kernel precision 1e-12, 1e-9, 200)"},
    ok(), ci=False, denominator="numeric", reuse=["C02.api_replay"],
    keys=plan(AREA_KEYS + spectrum_keys(nb("E_RB1")), "replay", "C02.api_replay (values and errors, bitwise)")
    + plan(["diagnostics.kernel_calls"], "diagnostic", "one new interior point per call plus the sentinels"),
    checks=["values and errors bitwise equal to the API replay of C02; call counts and batch indices change "
            "as registered"])
row("Z1", "Z", "api_spectrum", {"material": "SiZ", **window(0.5, 2), "primary": "D1000",
                                "output_edges_eV": "E_Z1"},
    ok(attributes={"kernel_derived_exact_zero": True, "inverse_window_area_A": None,
                   "inverse_reason": "nonpositive_or_unresolved_area", "spectrum_equals_zero_loss": True}),
    ci=True, denominator="numeric",
    keys=plan(AREA_KEYS, "exact_zero_or_null", "exact 0; inverse: value, interval None and the reason")
    + plan(spectrum_keys(nb("E_Z1"), SATELLITES + SENSITIVITIES), "exact_zero_or_null", "exact 0")
    + plan(spectrum_keys(nb("E_Z1"), ["zero_loss"]), "identity", "G_D(0) of the primary (exact)")
    + plan(spectrum_keys(nb("E_Z1"), ["spectrum"]), "identity", "equal to zero_loss (bitwise)")
    + plan(["diagnostics.resolution_basis", "diagnostics.candidates", "diagnostics.kernel_calls"], "diagnostic",
           "null with the reason zero_response; empty list; 0 calls"),
    note="zero-strength material with a gap of 1.12 eV; the window straddles the gap and is accepted")
row("Z2", "Z", "api_spectrum", {**AL150, "primary": "D1000zero", "output_edges_eV": "E_RB1"},
    ok(attributes={"all_spectrum_outputs_exact_zero": True, "primary_support_fields": None,
                   "budget_ratio_of_zero_keys": None, "budget_flag": "zero_budget_zero_error",
                   "area_keys_judged_as_usual": True},
       warnings_exclude=["E0_outside_primary_support", "primary_support_below_window_upper"]),
    ci=False, denominator="numeric", reuse=["RB1a.result"],
    keys=plan(spectrum_keys(nb("E_RB1")), "exact_zero_or_null", "exact 0")
    + plan(AREA_KEYS[:11], "identity", "RB1a.result within the sum of the two reported errors")
    + plan(["diagnostics.error_budget.outputs[*].ratio", "diagnostics.primary_support"], "diagnostic",
           "null with the flag zero_budget_zero_error only for keys with e = B = 0; all four values null"),
    checks=["areas within the sum of the reported errors of the areas of RB1a"])

# ---------------------------------------------------------------- N: rejections and non-convergence
AREA_AL = {"material": "Al"}
N_WINDOWS = (("N01", num(0.0), num(50.0), "a <= 0"), ("N02", num(50.0), num(50.0), "a >= b"),
             ("N03", num(1.0), num(1000.0), "b >= E0"), ("N04", num(float("nan")), num(50.0), "non-finite"),
             ("N05", num(2.0), num(2.0 + 100 * H2), "100 ulp wide: node separation fails on one panel"))
for rid, a, b, note in N_WINDOWS:
    row(rid, "N", "api_area", {**AREA_AL, "loss_window_eV": [a, b]},
        rejected("invalid_window", kernel_evaluations=0), ci=True, denominator="state", note=note)
row("N06", "N", "api_area", {"material": "Si", **window(0.5, 1.0)},
    rejected("window_not_above_gap", kernel_evaluations=0), ci=True, denominator="state")
row("N07", "N", "api_area", {"material": "Si", **window(1.0, 5)},
    rejected("window_not_above_gap", kernel_evaluations=0), ci=True, denominator="state")
row("N08", "N", "api_area", {"material": "Si", **window(A_N08, 5)},
    {"status": "error", "exception": "DielectricGapEdgeError", "reason": None,
     "raised_in": "first M3a call", "public_attributes_preserved": True},
    ci=True, denominator="state",
    keys=plan(["_call_kernel.calls"], "diagnostic", "1 (the first batch raises before any panel)"), note="a = 1.12 + 5e-13, inside the upstream exclusion width 1.12e-12 eV")
row("N09", "N", "api_area", {"material": "Si", **window(A_N09, 5)},
    ok(warnings_include=SI_WARN), ci=False, denominator="state",
    note="a = 1.12 + 2e-12, just outside the exclusion width; full run only (CI time)")
SPEC_AL = {**AL150, "output_edges_eV": "E_RB1"}
N_PRIMARY = (
    ("N10", {"kind": "delta", "energies_eV": [], "intensities": []}, "empty"),
    ("N11", {"kind": "both", "delta": "D1000", "binned": "B_RB3"}, "delta and bins together"),
    ("N12", {"kind": "delta", "energies_eV": [999.0, 1000.0], "intensities": [1.0]}, "length mismatch"),
    ("N13", {"kind": "delta", "energies_eV": [1000.0], "intensities": [-1.0]}, "negative intensity"),
    ("N14", {"kind": "delta", "energies_eV": [1000.0], "intensities": ["inf"]}, "non-finite intensity"),
    ("N15", {"kind": "binned", "edges_eV": [998.0, 1001.0, 999.0], "intensities": [0.3, 0.7]}, "non-monotonic edges"),
    ("N16", {"kind": "delta", "energies_eV": [1000.0], "intensities": [hx(2.0 ** 901)]},
     "intensity 2^901 above the range [2^-900, 2^900] (given as binary64 hex)"),
)
for rid, prim, note in N_PRIMARY:
    row(rid, "N", "api_spectrum", {**SPEC_AL, "primary": prim},
        rejected("invalid_primary", kernel_evaluations=0), ci=True, denominator="state", note=note)
for rid, edges, note in (("N17", [1000.0], "one point"), ("N18", [960.0, 980.0, 980.0], "not strictly increasing"),
                         ("N19", [-1.0, 0.0, 1.0], "negative edge")):
    row(rid, "N", "api_spectrum", {**AL150, "primary": "D1000", "output_edges_eV": {"values": edges}},
        rejected("invalid_output_edges", kernel_evaluations=0), ci=True, denominator="state", note=note)
row("N23", "N", "api_area", dict(SW150),
    {"status": "error", "exception": "ValueError", "reason": None,
     "message_contains": ["integration limit too small"], "raised_in": "first M3a call"},
    ci=True, denominator="state", settings={"kernel_limit": 1})
row("N24", "N", "api_area", {"material": "N1", **window(2.1, 10)},
    rejected("resolution_exceeds_panel_limit", attributes={"required_panels": 401}, kernel_evaluations=0),
    ci=True, denominator="state",
    keys=plan(["exception.required_panels"], "diagnostic", "401") + plan(["kernel_evaluations"], "diagnostic", "0"),
    note="k_parallel gives the same count; only k_total is a row")
row("N25", "N", "api_area", dict(SW150),
    unconverged("loss_budget_not_met", attributes={"dominant_part": "interpolation"}),
    ci=True, denominator="state", settings={"panel_limit": 9},
    keys=plan(["exception.dominant_part", "exception.diagnostics.required_panels"], "diagnostic",
              "interpolation; 9"),
    note="the 9 initial panels are evaluated; the check before the next bisection stops the run")
INJ26 = {"interpolation_nodes": "error 1e-4 in every column", "verification_nodes": "error 0"}
for rid, epsabs, expect, checks in (
        ("N26", 2e-4, unconverged("kernel_error_exceeds_budget"), []),
        ("N27", 3.5e-4, ok(),
         ["the diagnostic sum of the point errors of the area is 3e-4 A^-1 within rounding (4e-4 if the "
          "half-width factor is dropped)"])):
    row(rid, "N", "layer_a",
        analytic(RECT25, 2, 5, features={"jumps_eV": [2, 5]}, initial_edges_eV=[2, 4, 5], injection=INJ26),
        expect, ci=True, denominator="state", settings={"panel_limit": 2, "area_epsabs": epsabs}, checks=checks,
        keys=(plan(["diagnostics.error_budget.outputs[window_area.total].kernel_point_error"], "error_part",
                   "analytic: (2+1)*1e-4 A^-1, relative 1e-12") if rid == "N27" else []))
row("N28", "N", "api_area", {"material": "synthetic_wide", **window(5, 6)},
    unconverged("kernel_not_converged",
                attributes={"node_kind": "sentinel", "loss_eV": 5.0, "index": 0, "stage": 0,
                            "diagnostics.error_budget.outputs": {}, "diagnostics.error_budget.max_ratio": None,
                            "inner_exception": "DepthAverageConvergenceError",
                            "inner_reason": "panel_floor_exceeds_budget"}),
    ci=True, denominator="state",
    keys=plan(["exception.node_kind", "exception.loss_eV", "exception.index", "exception.inner_exception",
               "exception.inner_reason"], "diagnostic", "prereg 5.2"), settings={"kernel_epsabs": 1e-30, "kernel_epsrel": 1e-15, "kernel_limit": 200})
row("N29", "N", "cli_spectrum", {"config_material": "synthetic_wide",
                                 "config_loss_grid_eV": {"min": 1, "max": 100, "step": 1},
                                 "config_electron_energy_eV": 1000, "config_integration": CONFIG_INTEGRATION,
                                 "argv": {"--loss-window-eV": [1, 50], "--weight-length-A": 20,
                                          "--epsilon-momentum": "k_total"},
                                 "spectrum": {"primary": "D1000", "output_edges_eV": "E_RB1",
                                              "edges_form": "array"},
                                 "injection": "the wrapped _call_kernel takes the first batch, picks the first "
                                              "index with a < T < b and raises a real "
                                              "DepthAverageConvergenceError(reason='budget_not_met') carrying that "
                                              "loss, index, a column, an estimate and a budget"},
    {"status": "not_converged", "exception": "SingleLossConvergenceError", "reason": "kernel_not_converged",
     "attributes": {"node_kind": "interior", "inner_exception": "DepthAverageConvergenceError",
                    "inner_reason": "budget_not_met", "stage": 0,
                    "diagnostics.error_budget.outputs": {}, "diagnostics.error_budget.max_ratio": None,
                    "index_loss_column_and_inner_attributes": "those of the injected exception"},
     "exit_code": 2, "output_dir_after": "absent"},
    ci=True, denominator="state",
    keys=plan(["stderr.reason", "stderr.node_kind", "stderr.index", "stderr.loss_eV", "stderr.column",
               "stderr.inner_reason", "stderr.inner_attributes", "exit_code", "output_dir_after",
               "injection_fired"], "diagnostic",
              "prereg 5.2; a normal C02 passes through the same _call_kernel at least once"))
row("N30", "N", "api_spectrum", {**SW150, "primary": "D1000", "output_edges_eV": "E_RB1"},
    rejected("cell_limit_exceeded", kernel_evaluations=0), ci=True, denominator="state",
    settings={"max_cells": 88},
    keys=plan(["kernel_evaluations"], "diagnostic", "0"),
    note="the initial cell count is 89: 10, 11, 10, 10, 11, 10, 9, 9, 9 from the low side (prereg 3.3)")
row("N31", "N", "api_area", dict(AL150), rejected("invalid_settings", kernel_evaluations=0),
    ci=True, denominator="state", settings={"epsrel": 1e-17})
row("N32a", "N", "layer_a",
    analytic({"form": "rectangle", "support_eV": [2, 3], "term_13a": 0.10, "other_terms": 0.0}, 2, 3,
             features={"kinks_eV": [num(2.5), num(math.nextafter(2.5, 3.0))]}),
    rejected("panel_width_unresolved", kernel_evaluations=0), ci=True, denominator="state",
    note="two declared boundaries one ulp apart")
C32 = 2.0 + 64 * H2
row("N32b", "N", "layer_a",
    analytic({"form": "bump", "expression": "4*K0*u*(1-u), u=(T-c)/h on [c,c+h], 0 elsewhere", "K0": 1e6,
              "c_eV": num(C32), "h_eV": num(H2), "column": "term_13a", "other_terms": 0.0},
             2.0, 2.0 + 256 * H2,
             features={"kinks_eV": [num(C32), num(C32 + H2)], "min_variation_scale_eV": num(H2 / 2)}),
    rejected("panel_width_unresolved", attributes={"panel_width_ulp": 1, "panel_eV": "a panel of the window"},
             kernel_evaluations=0),
    ci=True, denominator="state",
    keys=plan(["exception.panel_width_ulp", "exception.panel_eV", "kernel_evaluations"], "diagnostic",
              "1; a 1-ulp panel inside the window; 0 calls of the kernel passed to integrate_window"),
    note="256 initial panels of 1 ulp; boundaries are distinct floats, inner nodes coincide")
row("N33", "N", "layer_a",
    analytic({"form": "constant", "value": "inf", "support_eV": [2, 3], "column": "term_13a"}, 2, 3,
             features={"min_variation_scale_eV": num(0.5)}),
    unconverged("nonfinite_output"), ci=True, denominator="state")
row("N34", "N", "layer_a",
    analytic({"form": "constant", "value": 1.0, "support_eV": [2, 3], "column": "term_13a",
              "other_terms": 0.0}, 2, 3, features={"min_variation_scale_eV": num(0.5)},
             injection={"all_nodes": "error exp(-((T-2.5)/0.001)^2) in every column"}),
    ok(attributes={"refinement_stages_min": 1}), ci=True, denominator="state", settings={"area_epsabs": 0.01},
    keys=plan(["diagnostics.stages"], "diagnostic", ">= 1"),
    note="must not stop on the parent panel; succeeds after the refinement the rule requires")
for rid, limit, expect in (
        ("N35a", 23, rejected("resolution_exceeds_panel_limit", attributes={"required_panels": 24},
                              kernel_evaluations=0)),
        ("N35b", 24, ok(attributes={"required_panels": 24}))):
    row(rid, "N", "layer_a", analytic(RECT25, 2, 5, features={"jumps_eV": [2, 5],
                                                             "min_half_width_eV": num(0.0625)}),
        expect, ci=True, denominator="state", settings={"panel_limit": limit},
        keys=plan(["required_panels"], "diagnostic", "24"))
row("N36", "N", "api_spectrum", {**SW150, "primary": "D1000", "output_edges_eV": "E_RB1"},
    unconverged("cell_limit_reached"), ci=True, denominator="state", settings={"max_cells": 89},
    keys=plan(["exception.diagnostics.initial_cells"], "diagnostic", "89"),
    note="equal to the initial cell count; stops at the check before the first refinement")

# ---------------------------------------------------------------- C: CLI on synthetic_wide (same process, main())
CLI_OK = {"status": "complete", "exit_code": 0, "output_dir_after": "created"}
CLI_ERR = {"status": "error", "exit_code": 2, "output_dir_after": "absent"}
CLI_BASE = {"config_material": "synthetic_wide", "config_loss_grid_eV": {"min": 1, "max": 100, "step": 1},
            "config_electron_energy_eV": 1000, "config_integration": CONFIG_INTEGRATION,
            "argv": {"--loss-window-eV": [1, 50], "--weight-length-A": 20, "--epsilon-momentum": "k_total"}}


def cli(extra_argv=None, **more):
    d = json.loads(json.dumps(CLI_BASE))
    d["argv"].update(extra_argv or {})
    d.update(more)
    return d


SPEC_C02 = {"primary": "D1000", "output_edges_eV": "E_RB1", "edges_form": "array"}
row("C01", "C", "cli_area", cli({"--weight-length-source": "design example"}),
    {**CLI_OK, "files": ["window_area.csv", "run.json"], "run_json": {"primary": None,
                                                                      "spectrum_definition.config_loss_grid_used": False}},
    ci=True, denominator="replay_and_state",
    keys=plan(AREA_KEYS, "replay", "API replay with diagnostics.effective_settings (bitwise)")
    + plan(["run.json paths (prereg 4.4)", "run.json results", "window_area.csv header and order",
            "window_area.csv numeric columns"], "diagnostic",
           "prereg 4.3, 4.4; results and every numeric CSV column (values, errors, budgets) equal the API replay"),
    checks=["API replay with the effective arguments of run.json: all columns bitwise equal"])
row("C02", "C", "cli_spectrum", cli(spectrum=SPEC_C02, comparison_table="test_table_C02"),
    {**CLI_OK, "files": ["window_area.csv", "single_loss_spectrum.csv", "run.json"],
     "comparison": {"areas": "matched", "spectrum": "matched"}},
    ci=True, denominator="replay_and_state",
    keys=plan(AREA_KEYS + spectrum_keys(nb("E_RB1")), "replay",
              "API replay with diagnostics.effective_settings (bitwise)")
    + plan(["comparison.areas", "comparison.spectrum", "run.json paths (prereg 4.4)",
            "run.json results", "window_area.csv numeric columns", "single_loss_spectrum.csv header and order",
            "single_loss_spectrum.csv numeric columns", "_call_kernel.calls >= 1"], "diagnostic",
           "matched; matched; prereg 4.3, 4.4; results and every numeric CSV column equal the API replay"),
    checks=["API replay with the effective arguments of run.json: columns and outsides bitwise equal",
            "matched against the test version of the completed comparison table"])
C_DISPLAY = (
    ("C04", cli(spectrum={"primary": "D1000", "output_edges_eV": "E_C04", "edges_form": "min_max_step"}),
     {"areas": "kernel_conditions_checked", "spectrum": "kernel_conditions_checked"},
     "same kernel conditions as C02, other output edges (given as min/max/step)"),
    ("C05", cli({"--hole-charge": 0.5}, spectrum=SPEC_C02),
     {"areas": "not_matched", "spectrum": "not_matched"}, "only Q_h differs"),
    ("C06", cli({"--weight-length-source": "another source text"},
                spectrum={**SPEC_C02, "json_whitespace_and_key_order": "changed",
                          "description": "another description", "included_broadening": "another text"}),
     {"areas": "matched", "spectrum": "matched"}, "only whitespace and provenance strings differ"),
    ("C07", cli(spectrum={**SPEC_C02, "primary": "D1000x2"}),
     {"areas": "kernel_conditions_checked", "spectrum": "kernel_conditions_checked"},
     "only the intensity differs"),
    ("C08", cli(spectrum=SPEC_C02),
     {"areas": "matched", "spectrum": "not_matched", "spectrum_missing_keys_judgement": "no_reference"},
     "registration lacking one mandatory spectrum key"),
)
for rid, inputs, comparison, note in C_DISPLAY:
    inputs["comparison_table"] = "test_table_C08_missing_key" if rid == "C08" else "test_table_C02"
    inputs["input_template"] = "C02 (template only; each row is its own CLI run and reads no C02 artifact)"
    row(rid, "C", "cli_spectrum", inputs, {**CLI_OK, "comparison": comparison}, ci=True, denominator="state",
        keys=plan(["comparison.areas", "comparison.spectrum", "warnings[outside_matched_scope]"]
                  + (["comparison.key_status[spectrum][outside.low.single_loss_total]"] if rid == "C08" else []),
                  "diagnostic", "prereg 4.6"),
        checks=["comparison display in run.json and the matching warning code"]
        + (["the table test_table_C08_no_reference (same key present with status no_reference) gives the same "
            "display through the same judgement after loading"] if rid == "C08" else []), note=note)
row("C09", "C", "cli_area", cli({"--max-cells": 10}),
    {**CLI_ERR, "stderr_json": {"status": "error", "error_type": "ValueError", "reason": "absent"}},
    ci=True, denominator="state",
    keys=plan(["exit_code", "stderr.error_type", "output_dir_after"], "diagnostic", "2; ValueError; absent"),
    checks=["exit code, stderr JSON, no output directory"], note="spectrum-only argument given to the area command")
row("C10a", "C", "cli_spectrum", cli(spectrum={**SPEC_C02, "defect": "axis is 'binding_energy_eV'"}),
    {**CLI_ERR, "stderr_json": {"status": "error", "error_type": "SingleLossInputError",
                                "reason": "invalid_schema"}},
    ci=True, denominator="state", checks=["exit code, stderr JSON, no output directory"])
row("C10b", "C", "cli_spectrum", cli({"--loss-window-eV": [1, "nan"]}, spectrum=SPEC_C02),
    {**CLI_ERR, "stderr_json": {"status": "error", "error_type": "SingleLossInputError",
                                "reason": "invalid_window", "loss_window_eV": [1.0, "nan"]}},
    ci=True, denominator="state", checks=["exit code, stderr JSON, no output directory"])
c10c = cli(spectrum=SPEC_C02)
del c10c["argv"]["--epsilon-momentum"]
row("C10c", "C", "cli_spectrum", c10c,
    {**CLI_ERR, "stderr_json": {"status": "error", "error_type": "ValueError", "reason": "absent"}},
    ci=True, denominator="state", checks=["exit code, stderr JSON, no output directory"],
    note="required argument --epsilon-momentum missing")

# ---------------------------------------------------------------- V: existing API and CLI unchanged from tag v0.2.0
V_CHECK = ["bitwise comparison with the run of tag v0.2.0 (commit fixed in the table) extracted by git archive "
           "and run in one isolated child process for all V rows; the loaded module path is checked; time "
           "stamps, durations and absolute paths are excluded, and the version is checked against each "
           "__version__", "CI uses the 125-point profile, the full run the 499-point profile (v_profiles)"]
V_COMMON = {"epsilon_momentum": "k_total", "hole_charge": 1.0,
            "precision": {"api": [1e-12, 1e-9, 200], "config_integration": CONFIG_INTEGRATION}}
V = (("V01", "existing_api", {"function": "effective_cross_section_fixed_depth", "material": "Al",
                              "depth_A": 10.0, "loss_eV": "GA"}, ok(), None),
     ("V02", "existing_api", {"function": "effective_cross_section_fixed_depth", "material": "Si",
                              "depth_A": 10.0, "loss_eV": "GA"}, ok(), None),
     ("V03", "existing_api", {"function": "effective_cross_section_depth_averaged", "material": "Al",
                              "weight_length_A": 20.0, "loss_eV": "GA"}, ok(), None),
     ("V04", "existing_api", {"function": "effective_cross_section_depth_averaged", "material": "Si",
                              "weight_length_A": 20.0, "loss_eV": "GA"}, ok(), None),
     ("V05", "existing_cli", {"command": "xps-fixed-depth", "material": "Al", "depth_A": 10.0,
                              "config_loss_grid_eV": "GC"}, CLI_OK, None),
     ("V06", "existing_cli", {"command": "xps-depth-averaged", "material": "Al", "weight_length_A": 20.0,
                              "config_loss_grid_eV": "GC"}, CLI_OK, None),
     ("V07", "existing_api", {"function": "effective_cross_section_depth_averaged", "material": "P1",
                              "weight_length_A": 20.0, "loss_eV": [1.999]},
      {"status": "error", "exception": "FixedDepthGapDomainError", "reason": "surface_pole_not_excluded",
       "same_type_message_and_public_attributes_as_tag": True},
      [{"function": "effective_cross_section_fixed_depth", "material": "P1", "depth_A": 10.0, "loss_eV": [1.999],
        "exception": "FixedDepthGapDomainError", "reason": "surface_pole_not_excluded",
        "existing_regression": "tests/test_m3a_migration.py (M3a V14(b) failure control fail_p1)"},
       {"function": "effective_cross_section_depth_averaged", "material": "Si", "weight_length_A": 20.0,
        "loss_eV": "np.r_[GA_full[:10], 1.12]", "exception": "DielectricGapEdgeError", "offending_indices": [10],
        "existing_regression": "tests/test_m3a_expected_states.py (M3a RJ1)"}]),
     ("V08", "existing_cli", {"command": "capabilities"},
      {"status": "complete", "exit_code": 0, "output_dir_after": None, "stdout_json_status": "experimental",
       "supported_prefix_unchanged": True, "existing_blocks_equal": True,
       "appended_supported": ["xps_loss_window_area", "xps_single_loss"]}, None))
V_NOTES = {"V07": "the two auxiliary controls are reported separately and are not rows",
           "V08": "status complete is the harness reading of exit code 0; the JSON itself says experimental"}
_OUT = ("tag v0.2.0, same environment; stdout output_dir is checked against the real output directory of each "
        "version and then excluded")
V_ORACLE = {"V05": _OUT, "V06": _OUT,
            "V08": "tag v0.2.0, same environment; version is checked against each __version__; the old supported "
                   "prefix, the existing blocks and the two appended names are checked separately"}
for rid, entry, inputs, expect, aux in V:
    row(rid, "V", entry, {**V_COMMON, **inputs}, expect, ci=True, denominator="invariance", checks=V_CHECK,
        keys=plan(["all_fields"], "replay", V_ORACLE.get(rid, "tag v0.2.0, same environment")),
        auxiliary=aux, note=V_NOTES.get(rid))

# ---------------------------------------------------------------- L: low-loss record (no pass/fail)
for rid, a in (("L01", 0.2), ("L02", 0.05), ("L03", 0.01)):
    row(rid, "L", "api_spectrum", {"material": "Al", **window(a, 50), "primary": "D1000",
                                   "output_edges_eV": "E_L"},
        {"status": "record_only"}, ci=False, denominator="none",
        keys=plan(["terminal state as observed", "cutoff_sensitivity (11 columns)", "error_budget parts"],
                  "record", "no pass or fail"),
        note="recorded whatever the terminal state; a = 1 eV is RB1")

# ---------------------------------------------------------------- checks
ids = [r["id"] for r in rows]
need(len(ids) == len(set(ids)), "duplicate row id")
need(len(rows) == 100, f"expected 100 rows, got {len(rows)}")
FAMILY_COUNTS = {"RA": 17, "CE": 6, "RB": 12, "I": 6, "Z": 2, "N": 35, "C": 11, "V": 8, "L": 3}
family_counts = {f: sum(1 for r in rows if r["family"] == f) for f in FAMILY_COUNTS}
need(family_counts == FAMILY_COUNTS, f"family counts differ: {family_counts}")
CI_COUNTS = {"RA": 17, "CE": 0, "RB": 1, "I": 0, "Z": 1, "N": 34, "C": 11, "V": 8, "L": 0}
ci_counts = {f: sum(1 for r in rows if r["family"] == f and r["ci"]) for f in FAMILY_COUNTS}
need(ci_counts == CI_COUNTS, f"CI counts differ: {ci_counts}")
EXPECTED_N = ([f"N{i:02d}" for i in range(1, 20)] + [f"N{i}" for i in range(23, 32)]
              + ["N32a", "N32b", "N33", "N34", "N35a", "N35b", "N36"])
need(sorted(i for i in ids if i.startswith("N")) == sorted(EXPECTED_N), "N ids differ from the registered list")
by_id = {r["id"]: r for r in rows}
for r in rows:
    need(r["entry"] in ENTRIES, f"{r['id']}: unknown entry")
    need(r["denominator"] in DENOMINATORS, f"{r['id']}: unknown denominator")
    e = r["expect"]
    need(e["status"] in ("success", "complete", "error", "not_converged", "record_only"), f"{r['id']}: status")
    if e.get("exception") == "SingleLossInputError":
        need(e["reason"] in INPUT_REASONS, f"{r['id']}: unknown input reason")
    if e.get("exception") == "SingleLossConvergenceError":
        need(e["reason"] in CONVERGENCE_REASONS, f"{r['id']}: unknown convergence reason")
    part = (e.get("attributes") or {}).get("dominant_part")
    need(part is None or part in DOMINANT_PARTS, f"{r['id']}: unknown dominant part")
    basis = (e.get("attributes") or {}).get("resolution_basis")
    need(basis is None or basis in BASES, f"{r['id']}: unknown resolution basis")
    for artifact in r["reuse"]:
        target, _, name = artifact.partition(".")
        need(target in by_id and name in ("result", "api_replay"), f"{r['id']}: bad reuse {artifact}")
        need(not r["ci"] or by_id[target]["ci"], f"{r['id']}: CI row reuses the full-only row {target}")
    for entry_ in r["check_plan"]:
        parts = entry_.split("|")
        need(len(parts) == 3 and parts[1] in CHECK_KINDS, f"{r['id']}: bad plan entry {entry_}")
    need(len(set(r["check_plan"])) == len(r["check_plan"]), f"{r['id']}: duplicate plan entry")
    inp = r["inputs"]
    m = inp.get("material", inp.get("config_material"))
    need(m is None or m in MATERIALS, f"{r['id']}: unknown material {m}")
    for case in r["auxiliary_controls"] or []:
        need(case["material"] in MATERIALS, f"{r['id']}: unknown material in an auxiliary control")
    spec = inp.get("spectrum", inp)
    p, ed = spec.get("primary"), spec.get("output_edges_eV")
    need(not isinstance(p, str) or p in PRIMARIES, f"{r['id']}: unknown primary {p}")
    need(not isinstance(ed, str) or ed in EDGES, f"{r['id']}: unknown edges {ed}")
    unknown = set(r["settings"]) - set(DEFAULTS["settings"])
    need(not unknown, f"{r['id']}: unknown settings {unknown}")
    need(r["denominator"] != "none" or e["status"] == "record_only", f"{r['id']}: only records are unjudged")
need(by_id["N09"]["ci"] is False and by_id["N28"]["inputs"]["material"] == "synthetic_wide",
     "CI rules for N09 and N28 are not met")
DESCRIPTIVE_INPUT_FIELDS = ("input_template", "read_from")
PHYSICAL_DEFAULTS = ("electron_energy_eV", "weight_length_A", "hole_charge", "epsilon_momentum", "theta_deg")
SETTINGS_OF_ENTRY = {"api_area": AREA_SETTINGS, "api_spectrum": SPECTRUM_SETTINGS,
                     "api_spectrum_pointwise": SPECTRUM_SETTINGS}


def _norm(v):
    """Binary64 identity of a number given as a float, an int or a {"value", "hex"} record."""
    if isinstance(v, dict) and "hex" in v and "value" in v:
        return v["hex"]
    if isinstance(v, (int, float)) and not isinstance(v, bool):
        return float(v).hex()
    return v


def entry_key(r):
    """Identity of the run: the call, the effective settings and the replaced table; no descriptive field.
    API entries are compared after expanding the default physical arguments and numeric settings.
    CLI entries keep the raw inputs (replacement table, raw differences of C06)."""
    call = {k: v for k, v in r["inputs"].items() if k not in DESCRIPTIVE_INPUT_FIELDS}
    settings = dict(r["settings"])
    if r["entry"] in SETTINGS_OF_ENTRY:
        for name in PHYSICAL_DEFAULTS:
            call[name] = _norm(call.get(name, DEFAULTS[name]))
        settings = {k: _norm(settings.get(k, DEFAULTS["settings"][k])) for k in SETTINGS_OF_ENTRY[r["entry"]]}
        extra = set(r["settings"]) - set(SETTINGS_OF_ENTRY[r["entry"]])
        need(not extra, f"{r['id']}: settings {extra} are not arguments of {r['entry']}")
    return json.dumps([r["entry"], call, settings], sort_keys=True)


def find_duplicate(table):
    seen = {}
    for r in table:
        key = entry_key(r)
        if key in seen:
            return seen[key], r["id"]
        seen[key] = r["id"]
    return None


def variant(rid, **changes):
    """A copy of the table with one row changed (self-test of the duplicate check)."""
    table = json.loads(json.dumps(rows))
    target = next(r for r in table if r["id"] == rid)
    for field, value in changes.items():
        target[field].update(value)
    return table


need(find_duplicate(rows) is None, f"two rows are the same entry: {find_duplicate(rows)}")
# Counterexamples that must be rejected as duplicates (Codex review 3):
need(find_duplicate(variant("RB5", inputs={"weight_length_A": num(20.0)})) == ("RB1", "RB5"),
     "RB5 with an explicit length of 20 A was not rejected as a duplicate of RB1")
need(find_duplicate(variant("N36", settings={"max_cells": 88, "epsrel": 1e-6})) == ("N30", "N36"),
     "N36 with max_cells 88 and an explicit default epsrel was not rejected as a duplicate of N30")
need(find_duplicate(variant("C08", inputs={"comparison_table": "test_table_C02"})) == ("C02", "C08"),
     "C08 with the table of C02 was not rejected as a duplicate of C02")
need(find_duplicate(variant("C06", inputs={"comparison_table": "test_table_C02"})) is None,
     "the raw input differences of C06 must stay in the key")
# execution order: a row runs after the rows whose artifacts it reads (no cycle)
order, placed = [], set()
while len(order) < len(rows):
    ready = [r["id"] for r in rows if r["id"] not in placed
             and all(a.partition(".")[0] in placed for a in r["reuse"])]
    need(ready, "reuse cycle")
    order.extend(ready)
    placed.update(ready)


def kind_count(r, kinds):
    return sum(1 for e in r["check_plan"] if e.split("|")[1] in kinds)


NUMERIC_KINDS = ("direct_reference", "inverse_interval")
key_counts = {k: sum(kind_count(r, (k,)) for r in rows) for k in CHECK_KINDS}
numeric_by_family = {f: sum(kind_count(r, NUMERIC_KINDS) for r in rows if r["family"] == f) for f in FAMILY_COUNTS}
numeric_by_row = {r["id"]: kind_count(r, NUMERIC_KINDS) for r in rows if kind_count(r, NUMERIC_KINDS)}
need(numeric_by_family["CE"] == 66 and numeric_by_family["RB"] == 6 * 124 + 2 * 102 + 4 * 25 == 1048,
     f"planned numeric keys differ: {numeric_by_family}")
need(numeric_by_family["RA"] == 922 and all(numeric_by_family[f] == 0 for f in ("I", "Z", "N", "C", "V", "L")),
     f"planned numeric keys differ: {numeric_by_family}")
need(key_counts["state"] == 100, "every row has one terminal-state key")
all_keys = {r["id"]: {e.split("|")[0] for e in r["check_plan"]} for r in rows}
state_counts = {}
for r in rows:
    state_counts[r["expect"]["status"]] = state_counts.get(r["expect"]["status"], 0) + 1
den_counts = {d: sum(1 for r in rows if r["denominator"] == d) for d in DENOMINATORS}

STATE_CHANGE_RULE = ("detected when any key of detected_by is broken; an unrelated exception, a crash of the test "
                     "harness or a failure of the normal control is not a detection")


def mut(mid, cls, text, detect, note=None, controls=()):
    return {"id": mid, "class": cls, "mutation": text, "normal_control": "the same rows pass without the mutation",
            "detected_by": [{"row": r_, "keys": k_} for r_, k_ in detect], "state_change": STATE_CHANGE_RULE,
            "note": note,
            "unaffected_controls": [{"row": r_, "keys": k_, "condition": "equal to the normal version after the "
                                     "mutation"} for r_, k_ in controls]}


J1T = "single_loss_total"
MUTATIONS = [
    mut("M01", "required", "axis reversed",
        [("RA1", [f"bin[2].{J1T}", f"bin[3].{J1T}", f"outside.high.{J1T}"])]),
    mut("M02a", "required", "L dropped from the satellite (sum 0.30 instead of 0.60; w unchanged)",
        [("RA1", [f"bin[2].{J1T}", f"bin[3].{J1T}"])]),
    mut("M03", "required", "L replaced by 1/A in the satellite only (sum 1.0; S*L*A stays 0.60)",
        [("RA1", [f"bin[2].{J1T}", f"bin[3].{J1T}"])]),
    mut("M04", "required", "only the outer multiplier follows L; the length passed to M3a stays 20",
        [("RB5", ["_call_kernel.weight_length_A", f"bin[1].{J1T}", f"bin[2].{J1T}"])],
        controls=[("RB1", ["_call_kernel.weight_length_A"])]),
    mut("M05a", "required", "hole_charge not passed upstream (default 1 used)",
        [("RB6", ["_call_kernel.hole_charge", f"bin[1].{J1T}"]),
         ("RB7", ["_call_kernel.hole_charge", f"bin[1].{J1T}"])]),
    mut("M06", "required", "momentum convention swapped",
        [("RB1", [f"bin[1].{J1T}", f"bin[2].{J1T}"]), ("RB2", [f"bin[1].{J1T}", f"bin[2].{J1T}"]),
         ("RB3", [f"bin[0].{J1T}", f"bin[2].{J1T}"]), ("RB4", [f"bin[0].{J1T}", f"bin[2].{J1T}"])]),
    mut("M07", "required", "primary bin width multiplied once more in the density (hole in [998,1000] x3)",
        [("RA3", ["bin[2].single_loss_hole"])]),
    mut("M08a", "required", "panel half-width factor dropped (value path)",
        [("RA4a", ["window_area.term_13a", "window_area.total"])]),
    mut("M10", "required", "outside added to the opposite end",
        [("RA2", [f"outside.low.{J1T}", f"outside.high.{J1T}"])]),
    mut("M11a", "required", "outside obtained as a residual of the total",
        [("RA6", [f"outside.low.{J1T}"]), ("RA6b", [f"outside.high.{J1T}"])]),
    mut("M12", "required", "delta on an inner edge counted on both sides",
        [("RA1", ["bin[5].zero_loss", "bin[6].zero_loss"])]),
    mut("M13", "required", "negative hole clipped after integration (-2/225 becomes 0)",
        [("RA3", ["bin[2].single_loss_hole"])]),
    mut("M14a", "required", "stop after the 6 initial panels",
        [("RA13", ["window_area.term_13a", "window_area.total"])]),
    mut("M14b", "required", "failure of an interior node replaced by 0",
        [("N29", ["terminal_state", "stderr.reason", "stderr.node_kind", "stderr.index", "stderr.loss_eV",
                  "stderr.inner_reason", "exit_code", "output_dir_after"])],
        None),
    mut("M16", "required", "pre-splitting by the resolution scale, the panel-limit check and the node-separation "
        "check bypassed (the review states which of the three)",
        [("RA8b", ["terminal_state", "exception.required_panels", "kernel_evaluations"]),
         ("N32b", ["terminal_state", "exception.panel_width_ulp", "exception.panel_eV", "kernel_evaluations"])],
        "not only when success is returned; a failed allocation of a huge array is not a detection"),
    mut("M17", "required", "verification residual not propagated to outputs (initial satellite estimate 0)",
        [("RA12", [f"bin[0].{J1T}", "diagnostics.stages"])]),
    mut("M05b", "additional", "total rebuilt with coefficient 1 (0.18 instead of 0.165)",
        [("RA3", ["window_area.total"])]),
    mut("M05c", "additional", "total rebuilt with the correct Q_h (0.72 instead of 3.0)",
        [("RA10", [f"bin[0].{J1T}"])]),
    mut("M08b", "additional", "panel half-width factor dropped (error path): 4e-4 instead of 3e-4",
        [("N27", ["terminal_state",
                  "diagnostics.error_budget.outputs[window_area.total].kernel_point_error"])],
        "announced state change: N27 goes from success to kernel_error_exceeds_budget"),
    mut("M09", "additional", "G fixed at its midpoint value inside each cell",
        [("RA3", ["bin[2].single_loss_hole"])]),
    mut("M11b", "additional", "inner bin built from a cumulative difference",
        [("RA6", [f"bin[1].{J1T}"])]),
    mut("M15", "additional", "one-sided sensitivities swapped",
        [("RA7", ["bin[1].cutoff_sensitivity.plus.total", "bin[2].cutoff_sensitivity.minus.total"])]),
    mut("M18", "additional", "absolute value of the weight missing in the point error",
        [("RA9", [f"bin[0].{q}.kernel_point_error" for q in SATELLITES])]),
]
need(sum(1 for m in MUTATIONS if m["class"] == "required") == 16 and len(MUTATIONS) == 23, "mutation counts differ")
need(sorted(m["id"] for m in MUTATIONS) == sorted(
    ["M01", "M02a", "M03", "M04", "M05a", "M06", "M07", "M08a", "M10", "M11a", "M12", "M13", "M14a", "M14b",
     "M16", "M17", "M05b", "M05c", "M08b", "M09", "M11b", "M15", "M18"]), "mutation id set differs")
for m in MUTATIONS:
    for d in m["detected_by"] + m["unaffected_controls"]:
        need(d["row"] in by_id, f"{m['id']}: detection row missing")
        missing = [k for k in d["keys"] if k not in all_keys[d["row"]]]
        need(not missing, f"{m['id']}: {d['row']} has no planned key {missing}")

CONDITIONS_C02 = {
    "material": {"gap_eV": 0.0, "oscillators": [[20.0, 200.0, 8.0, 0.0]]},
    "electron_energy_eV": 1000.0, "weight_length_A": 20.0, "theta_deg": 0.0, "epsilon_momentum": "k_total",
    "hole_charge": 1.0, "loss_window_eV": [1.0, 50.0],
    "effective_settings": {k: DEFAULTS["settings"][k] for k in SPECTRUM_SETTINGS},
    "primary": {"kind": "delta", "energies_eV": [1000.0], "intensities": [1.0]},
    "output_edges_eV": [float(e) for e in E_RB1],
}
doc = {
    "schema": "xpsloss-m4a-expected-states-v4",
    "prereg": "docs/M4A_PREREG.md",
    "design": "work/consult/m4a/merged_design.md (6th edition)",
    "defaults": DEFAULTS,
    "columns": COLUMNS,
    "vocabulary": {"input_reasons": INPUT_REASONS, "convergence_reasons": CONVERGENCE_REASONS,
                   "dominant_parts": DOMINANT_PARTS, "resolution_bases": BASES, "entries": ENTRIES,
                   "denominators": DENOMINATORS, "check_kinds": CHECK_KINDS},
    "check_plan_format": "key|check_kind|oracle; numeric mandatory keys are direct_reference and "
                         "inverse_interval; a planned key stays in the ledger whatever its result",
    "key_groups": {"areas": AREA_KEYS, "spectrum_per_region": QUANTITIES + SENSITIVITIES,
                   "regions": "outside.low, bin[0..n-1], outside.high"},
    "v_profiles": V_PROFILES,
    "tag": TAG,
    "kernel_evaluation_counters": {
        "public_entries": "calls of xpsloss.xps_single_loss._call_kernel",
        "analytic_layer": "calls of the kernel(T, side) passed to integrate_window",
        "layer_b": "fold never calls a kernel"},
    "comparison_table_shape": {
        "version": "str", "registration_sha256": "str", "reference_sha256": "dict or null",
        "records": [{"entry": "loss_window_area or single_loss_spectrum",
                     "conditions": "normalized numeric conditions; see conditions_example_C02 (prereg 4.6)",
                     "key_status": {"areas": {"<key>": {"status": "pass|fail|no_reference|unresolved|not_run",
                                                        "reason": "str or null"}},
                                    "spectrum": "same shape; absent for the area entry"}}]},
    "conditions_example_C02": CONDITIONS_C02,
    "conditions_rules": "keys, nesting and array order as in the example; numbers are compared as binary64 values; "
                        "oscillators are [energy_eV, strength_eV2, width_eV, dispersion] in the order of the "
                        "config; the area entry has no primary and no output_edges_eV and its effective_settings "
                        "has the six area settings; binned primaries use edges_eV instead of energies_eV",
    "replacement_intensity_RA6": {"driver_S_abs": "2^60 + 3", "outside.low": 1, "bin[0]": "2^60",
                                  "bin[1] (RA6)": 2, "outside.high": 2},
    "ci_reduction_order": ["C05", "C07"],
    "execution_order": order,
    "materials": MATERIALS,
    "primaries": PRIMARIES,
    "output_edges": EDGES,
    "counts": {"rows": len(rows), "ci_rows": sum(r["ci"] for r in rows),
               "full_only_rows": sum(not r["ci"] for r in rows), "by_family": family_counts,
               "ci_by_family": ci_counts, "by_expected_status": state_counts, "by_denominator": den_counts,
               "planned_keys_by_kind": key_counts,
               "planned_numeric_keys": sum(numeric_by_family.values()),
               "planned_numeric_keys_by_family": numeric_by_family,
               "planned_numeric_keys_by_row": numeric_by_row},
    "rows": rows,
    "mutations": MUTATIONS,
}
text = json.dumps(doc, ensure_ascii=True, indent=1, sort_keys=True) + "\n"
data = text.encode("ascii")
if sys.argv[1:] == ["--check"]:
    same = OUT.is_file() and OUT.read_bytes() == data
    print(json.dumps({"identical": same, "sha256": hashlib.sha256(data).hexdigest()}))
    sys.exit(0 if same else 1)
OUT.write_bytes(data)
print(json.dumps({"sha256": hashlib.sha256(data).hexdigest(), **doc["counts"]}, indent=1))
