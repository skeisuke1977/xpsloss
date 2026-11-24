"""M3a pre-registration: expand the expected-state table (docs/M3A_PREREG.md section 3).

Writes expected_states_m3a.json next to this file and prints its SHA-256 and the counts.
Each run is one call of the API or the CLI. Each point key is
(material, E0, momentum, hole charge, weight length, loss as binary64 hex, precision set, theta,
entry, CLI source string, output-directory prior state and name). A run key adds the ordered
loss list. The script checks (with explicit exceptions, not assert) that
- run keys are unique, and no point key gets two categories;
- a non-rejection run never mixes categories (checked against a table built from all runs);
- rejection indices exist, are in range and are not repeated;
- every "either" (EI) run has one loss and a non-empty reason set from the registered vocabulary;
- every reference key and every 13a cross-length key that a run points to is defined.
It imports only the standard library and NumPy (not xpsloss).

Usage: python benchmarks/numerics/m3a/expected_states_m3a.py          (write the JSON)
       python benchmarks/numerics/m3a/expected_states_m3a.py --check  (compare byte for byte)
"""
import hashlib
import json
import math
import pathlib
import sys
from itertools import product

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / "expected_states_m3a.json"


class TableError(RuntimeError):
    pass


def need(cond, message):
    if not cond:
        raise TableError(message)


if sys.argv[1:] not in ([], ["--check"]):
    raise SystemExit("usage: expected_states_m3a.py [--check]")


# ---------------------------------------------------------------- materials
def _read(rel):
    path = ROOT / rel
    data = path.read_bytes()
    return json.loads(data.decode("utf-8")), hashlib.sha256(data).hexdigest()


AL_CFG, AL_SHA = _read("benchmarks/inputs/al_elfdb_v1.0.bulk.json")
SI_CFG, SI_SHA = _read("benchmarks/inputs/si_elfdb_v1.0.bulk.json")


def _lines(cfg):
    return [[o["energy_eV"], o["strength_eV2"], o["width_eV"], o["dispersion"]]
            for o in cfg["material"]["oscillators"]]


AL = _lines(AL_CFG)
SI = _lines(SI_CFG)
G = float(SI_CFG["material"]["gap_eV"])
need(AL_CFG["material"]["gap_eV"] == 0 and G == 1.12, "unexpected Al or Si input")

MATERIALS = {
    "Al": {"gap_eV": 0.0, "oscillators": AL, "from": f"benchmarks/inputs/al_elfdb_v1.0.bulk.json sha256 {AL_SHA}"},
    "Al0": {"gap_eV": 0.0, "oscillators": [[e, a, w, 0.0] for e, a, w, _ in AL], "from": "Al, dispersion 0"},
    "Y0": {"gap_eV": 0.0, "oscillators": [[18.0, 200.0, 3.0, 0.0]], "from": "synthetic"},
    "Y5": {"gap_eV": 0.0, "oscillators": [[18.0, 200.0, 3.0, 0.5]], "from": "synthetic"},
    "Z": {"gap_eV": 0.0, "oscillators": [[e, 0.0, w, d] for e, _, w, d in AL], "from": "Al, strengths 0"},
    "W4": {"gap_eV": 0.0, "oscillators": [[15.0, 180.0, 3e-4, 0.0]], "from": "synthetic (4f line)"},
    "W6": {"gap_eV": 0.0, "oscillators": [[15.0, 180.0, 3e-6, 0.0]], "from": "synthetic (4f line)"},
    "W7": {"gap_eV": 0.0, "oscillators": [[15.0, 180.0, 3e-7, 0.0]], "from": "synthetic (4f line)"},
    "W8": {"gap_eV": 0.0, "oscillators": [[15.0, 180.0, 3e-8, 0.0]], "from": "synthetic (4f line)"},
    "W7a": {"gap_eV": 0.0, "oscillators": [[15.0, 180.0, 3e-7, 2e-8]], "from": "synthetic (4f line)"},
    "W2": {"gap_eV": 0.0, "oscillators": [[15.0, 81.0, 4e-7, 0.0], [18.0, 180.0, 2e-7, 0.0]],
           "from": "synthetic (tests/test_numerics_4f.py)"},
    "Si": {"gap_eV": G, "oscillators": SI, "from": f"benchmarks/inputs/si_elfdb_v1.0.bulk.json sha256 {SI_SHA}"},
    "S0": {"gap_eV": G, "oscillators": [[e, a, w, 0.0] for e, a, w, _ in SI], "from": "Si, dispersion 0"},
    "Z0": {"gap_eV": G, "oscillators": [[e, 0.0, w, d] for e, _, w, d in SI], "from": "Si, strengths 0"},
    "Z1": {"gap_eV": G, "oscillators": SI + [[0.5, 0.0, 0.3, 0.5]], "from": "Si plus a zero-strength line"},
    "N1": {"gap_eV": 2.0, "oscillators": [[3.0, 6.0, 0.02, 0.5]], "from": "synthetic (M2-2)"},
    "N2": {"gap_eV": 2.0, "oscillators": [[3.0, 6.0, 0.2, 0.5]], "from": "synthetic (N1 width x10)"},
    "P1": {"gap_eV": 2.0, "oscillators": [[3.0, 8.1, 0.5, 0.5]], "from": "synthetic (M2-2)"},
    "P0": {"gap_eV": 2.0, "oscillators": [[3.0, 8.1, 0.5, 0.0]], "from": "synthetic (M2-2)"},
    "M1": {"gap_eV": 2.0, "oscillators": [[3.0, 8.1, 1.0, 0.0]], "from": "synthetic (M2-2)"},
    "UE": {"gap_eV": 2.0, "oscillators": [[1.9, 1.0, 0.05, 20.0]], "from": "synthetic (M2-2)"},
    "SP": {"gap_eV": 1.0, "oscillators": [[5.0, 40.0, 1.0, 0.0]], "from": "synthetic (M2-2)"},
}

# ---------------------------------------------------------------- arrays and settings
GA = 0.2 * np.arange(1, 500)
GC = np.linspace(0.2, 99.8, 499)
F = 0.2 / 4 * np.arange(4, 1997)
ARRAYS = {
    "GA": ("0.2*np.arange(1, 500)", GA),
    "GC": ("np.linspace(0.2, 99.8, 499)", GC),
    "F": ("0.2/4*np.arange(4, 1997)", F),
    "F2": ("F[::2]", F[::2]),
    "FRD": ("np.r_[F[::-4], F[:40:4]]", np.r_[F[::-4], F[:40:4]]),
    "GA10G": ("np.r_[GA[:10], 1.12]", np.r_[GA[:10], G]),
}
need(bool(np.all(F[::4] == GA)) and int(np.sum(GA != GC)) == 340, "unexpected grids")
# APIDEFAULT: the API is called without epsabs/epsrel/limit (its defaults are the DEF values).
# CLIDEFAULT: the CLI config has no "integration" block (the config defaults apply).
# CLIPREC: the API with the CLI default values passed explicitly (the API control of CLIDEFAULT).
SETTINGS = {"DEF": (1e-12, 1e-9, 200), "TIGHT": (1e-14, 1e-11, 200), "IMP": (1e-30, 1e-15, 200),
            "LOOSE": (1e-6, 1e-9, 200), "LIM1": (1e-12, 1e-9, 1),
            "APIDEFAULT": (1e-12, 1e-9, 200), "CLIDEFAULT": (1e-9, 1e-7, 200),
            "CLIPREC": (1e-9, 1e-7, 200)}
BOTH = ("k_parallel", "k_total")
SQ = math.sqrt(135.0)
Q3 = (0.0, 0.5, 1.0)
Q_NEXT = math.nextafter(1.0, 2.0)
GAPLESS = ("Al", "Al0", "Y0", "Y5")
V7_GAPLESS_T = (0.05, 0.2, 1.0, 5.0, 11.0, 15.6, 30.0, 60.0, 99.8)
V7_GAPPED_T = (*GA[:5], G - 1e-3, G - 1e-6, G + 1e-6, G + 1e-3, 1.2, 2.0, 5.0, 11.0, 16.8, 17.0, 30.0, 60.0, 99.8)
BAD_LENGTHS = ("0.0", "-1.0", "nan", "inf", "1+0j", "True")
JA_SOURCE = "文献値（例）: 合成の長さ"
JA_DIR = "深さ平均の出力"
GRID_GC = {"min": 0.2, "max": 99.8, "step": 0.2}

# Reason codes of DepthAverageConvergenceError (docs/M3A_PREREG.md section 4).
REASONS = ("panel_floor_exceeds_budget", "direct_sum_floor_exceeds_budget", "budget_not_met",
           "surface_denominator_unresolved", "oscillator_denominator_unresolved", "surface_resonance_unresolved",
           "d3_ill_conditioned", "coefficient_unresolved", "linewidth_unresolved", "momentum_bounds_unresolved",
           "weight_out_of_range", "response_out_of_range")
BASE = ["budget_not_met", "panel_floor_exceeds_budget"]


def h(x):
    return float(x).hex()


runs = []


def fam(family, cat, mats, e0s, moms, qs, lams, losses, *, settings="DEF", theta=0.0, entry="API",
        source=None, outdir="absent", outdir_name="ascii", grid=None, expect, routine, ref=False, t13a=False):
    """ref / routine may be bools or callables of (material, momentum, charge, length)."""
    if isinstance(losses, str):
        values = ARRAYS[losses][1]
        loss_spec = {"array": losses}
    else:
        values = np.atleast_1d(np.asarray(losses, float))
        loss_spec = {"values_hex": [h(t) for t in values], "values": [float(t) for t in values]}
    need(settings in SETTINGS, settings)
    for m, e, mo, q, la in product(mats, e0s, moms, qs, lams):
        lam = {"invalid": la} if isinstance(la, str) else {"hex": h(la), "value": float(la)}
        run = {
            "id": f"{family}#{len(runs)}", "family": family, "category": cat, "material": m,
            "electron_energy_eV": float(e), "epsilon_momentum": mo, "hole_charge": {"hex": h(q), "value": q},
            "weight_length_A": lam, "losses": loss_spec, "settings": settings, "theta_deg": theta,
            "entry": entry, "weight_length_source": source, "expect": expect,
            "routine": bool(routine(m, mo, q, la) if callable(routine) else routine),
        }
        if entry == "CLI":
            need(grid is not None, f"{family}: CLI run needs a config grid")
            run["output_dir_before"] = outdir
            run["output_dir_name"] = outdir_name
            run["cli_config"] = {"loss_grid_eV": grid,
                                 "integration": None if settings == "CLIDEFAULT" else settings}
        else:
            run["output_dir_before"] = None
        if ref if not callable(ref) else ref(m, mo, q, la):
            run["reference_keys"] = [ref_key(m, e, mo, la, t) for t in values]
        if t13a:
            run["term_13a_reference_keys"] = [ref_key(m, e, mo, 20.0, t) for t in values]
        run["_values"] = values
        runs.append(run)


def ref_key(m, e, mo, la, t):
    return f"{m}|{float(e)!r}|{mo}|{h(la)}|{h(t)}"


OK = {"status": "success"}
OK_NULLDIR = {"status": "success", "direct_sum_diagnostics": "null"}

# ---- MS must succeed: every point succeeds and all eleven outputs meet the run budget
fam("MS1", "MS", GAPLESS, (1000.,), BOTH, Q3, (5., 20., 80.), "GA", expect=OK,
    routine=lambda m, mo, q, la: la == 20. and q == 1.)
fam("MS2", "MS", ("Al",), (4000.,), BOTH, (1.,), (20.,), "GA", expect=OK, routine=True)
for arr in ("F", "F2", "FRD"):
    fam("MS3", "MS", ("Al", "Si"), (1000.,), ("k_total",), (1.,), (20.,), arr, expect=OK,
        routine=(lambda m, mo, q, la, a=arr: m == "Al" or a == "FRD"))
fam("MS4", "MS", ("Al", "Si", "Z"), (1000.,), ("k_total",), (1.,), (20.,), "GC", expect=OK, routine=True)
fam("MS4", "MS", ("Al",), (1000.,), ("k_total",), (1.,), (20.,), "GC", settings="APIDEFAULT", expect=OK, routine=True)
fam("MS4", "MS", ("Al",), (1000.,), ("k_total",), (1.,), (20.,), "GC", settings="CLIPREC", expect=OK, routine=True)
CLI_OK = {"status": "complete", "exit_code": 0, "output_dir_after": "created"}
fam("MS4", "MS", ("Al",), (1000.,), ("k_total",), (1.,), (20.,), "GC", entry="CLI", grid=GRID_GC,
    expect={**CLI_OK, "direct_sum_diagnostics": "null"}, routine=True)
fam("MS4", "MS", ("Al",), (1000.,), ("k_total",), (1.,), (20.,), "GC", entry="CLI", grid=GRID_GC,
    source=JA_SOURCE, outdir_name="ja", expect={**CLI_OK, "direct_sum_diagnostics": "null"}, routine=True)
fam("MS4", "MS", ("Al",), (1000.,), ("k_total",), (1.,), (20.,), "GC", entry="CLI", grid=GRID_GC,
    settings="CLIDEFAULT", expect=CLI_OK, routine=True)
fam("MS4", "MS", ("Si",), (1000.,), ("k_total",), (1.,), (20.,), "GC", entry="CLI", grid=GRID_GC,
    expect={**CLI_OK, "direct_sum_diagnostics": "nonnull"}, routine=True)
fam("MS4", "MS", ("Z",), (1000.,), ("k_total",), (1.,), (20.,), "GC", entry="CLI", grid=GRID_GC,
    expect={**CLI_OK, "low_loss_warning": False, "direct_sum_diagnostics": "null"}, routine=True)
for t in V7_GAPLESS_T:
    fam("MS5", "MS", GAPLESS, (1000.,), BOTH, Q3, (5., 20., 80.), [t], expect=OK, routine=True, ref=True)
for m in ("W4", "W6", "W7"):
    fam("MS6", "MS", (m,), (1000.,), BOTH, (0., 1.), (20.,), [SQ], expect=OK, routine=True, ref=True)
fam("MS7", "MS", ("Z",), (1000.,), BOTH, (0., 1.), (20.,), "GA", expect={"status": "success", "exact_zero": True},
    routine=True)
fam("MS7", "MS", ("Z",), (1000.,), ("k_total",), (1.,), (20.,), [5.0], settings="IMP",
    expect={"status": "success", "exact_zero": True}, routine=True)
for t in (1.0, 15.6, 30.0):
    fam("MS8", "MS", ("Al", "Y0"), (1000.,), BOTH, Q3, (1e-6, 1e-5), [t], expect=OK, routine=True, t13a=True)
    fam("MS9", "MS", ("Al", "Y0"), (1000.,), BOTH, Q3, (1e5, 1e6), [t], expect=OK, routine=True, t13a=True)
fam("MS10", "MS", ("Si", "S0"), (1000.,), BOTH, Q3, (5., 20., 80.), "GA", expect=OK,
    routine=lambda m, mo, q, la: la == 20. and (q == 1. or (m == "Si" and mo == "k_total" and q == 0.5)))
for t in V7_GAPPED_T:
    fam("MS11", "MS", ("Si", "S0"), (1000.,), BOTH, Q3, (5., 20., 80.), [t], expect=OK, ref=True,
        routine=(lambda m, mo, q, la, t=t: (la == 20. and q in (0., 1.))
                 or (m == "Si" and la == 5. and q == 1. and t == GA[0])))
fam("MS12", "MS", ("Z0",), (1000.,), BOTH, (0., 1.), (20.,), "GA10G",
    expect={"status": "success", "exact_zero": True, "gap_warnings": False}, routine=True)
fam("MS12", "MS", ("Z1",), (1000.,), BOTH, (1.,), (20.,), "GA",
    expect={"status": "success", "same_gap_acceptance_as": "Si"}, routine=lambda m, mo, q, la: mo == "k_total")
fam("MS13", "MS", ("N2",), (1000.,), BOTH, (0., 1.), (20.,), [4.0], expect=OK_NULLDIR, routine=True, ref=True)
for t in GA[:5]:
    fam("MS14", "MS", ("Si", "S0"), (1000.,), BOTH, (0., 0.5), (1e-6, 1e-5), [t], expect=OK,
        routine=lambda m, mo, q, la: la == 1e-5, ref=True, t13a=True)
    fam("MS14", "MS", ("Si", "S0"), (1000.,), BOTH, (1.,), (1e-6, 1e-5), [t], settings="LOOSE", expect=OK,
        routine=lambda m, mo, q, la: la == 1e-5, ref=True, t13a=True)
for t in (1.2, 5.0, 16.8):
    fam("MS14", "MS", ("Si", "S0"), (1000.,), BOTH, Q3, (1e-6, 1e-5), [t], expect=OK,
        routine=lambda m, mo, q, la: la == 1e-5, t13a=True)
for t in (1.0, 1.2, 5.0, 16.8):
    fam("MS15", "MS", ("Si",), (1000.,), BOTH, Q3, (1e5, 1e6), [t], expect=OK, routine=True, t13a=True,
        ref=(t == 1.0))
for t in V7_GAPLESS_T:
    fam("MS16", "MS", ("Al",), (1000.,), BOTH, (-0.5, 1.5), (20.,), [t], expect=OK, routine=True, ref=True)
for t in V7_GAPPED_T:
    fam("MS16", "MS", ("Si",), (1000.,), BOTH, (-0.5, 1.5), (20.,), [t], expect=OK, routine=True, ref=True)
for t in GA[:5]:
    fam("MS17", "MS", ("Si",), (1000.,), BOTH, (Q_NEXT,), (20.,), [t], settings="LOOSE", expect=OK_NULLDIR,
        routine=True, ref=True)

# ---- NC specific non-convergence
fam("NC1", "NC", ("Al", "Si"), (1000.,), ("k_total",), (1.,), (20.,), [5.0], settings="IMP",
    expect={"status": "not_converged", "exception": "DepthAverageConvergenceError",
            "reason": "panel_floor_exceeds_budget", "index": 0, "subdivisions": 0}, routine=True)
fam("NC2", "NC", ("Al",), (1000.,), ("k_total",), (1.,), (20.,), [5.0], settings="IMP", entry="CLI",
    grid={"min": 5.0, "max": 5.0, "step": 1.0},
    expect={"status": "not_converged", "exit_code": 2, "output_dir_after": "absent",
            "reason": "panel_floor_exceeds_budget"}, routine=True)
for t in GA[:5]:
    fam("NC3", "NC", ("Si", "S0"), (1000.,), BOTH, (1.,), (1e-6, 1e-5), [t],
        expect={"status": "not_converged", "exception": "DepthAverageConvergenceError",
                "reason": "direct_sum_floor_exceeds_budget", "index": 0, "column": "total", "subdivisions": 0},
        routine=(lambda m, mo, q, la, t=t: la == 1e-5 or t == GA[0]))


# ---- EI success matching the reference, or one of the named non-convergence reasons; one loss per run
def either(reasons):
    return {"status": "either", "exception": "DepthAverageConvergenceError", "allowed_reasons": reasons,
            "index": 0, "on_success": "finite nonnegative estimates within B_run; V7 three-valued comparison"}


for t in (0.2, 0.4):
    fam("EI1", "EI", ("Si",), (1000.,), ("k_total",), (1.,), (1.,), [t],
        expect=either(BASE + ["direct_sum_floor_exceeds_budget"]), routine=True, ref=True)
for t in (0.05, 0.2):
    fam("EI2", "EI", ("Al",), (1000.,), ("k_total",), (1.,), (20.,), [t], settings="TIGHT", expect=either(BASE),
        routine=True, ref=True)
SURF = BASE + ["surface_denominator_unresolved"]
fam("EI3", "EI", ("W8",), (1000.,), BOTH, (0., 1.), (20.,), [SQ], expect=either(SURF), routine=True, ref=True)
for t in (SQ - 1.5e-7, SQ + 1.5e-7):
    fam("EI3", "EI", ("W7",), (1000.,), ("k_parallel",), (1.,), (20.,), [t], settings="TIGHT", expect=either(SURF),
        routine=True, ref=True)
fam("EI3", "EI", ("W2",), (1000.,), ("k_parallel",), (1.,), (20.,), [12.0], settings="TIGHT", expect=either(SURF),
    routine=True, ref=True)
fam("EI3", "EI", ("W7a",), (1000.,), BOTH, (1.,), (20.,), [SQ], expect=either(SURF), routine=True, ref=True)
fam("EI4", "EI", ("N1",), (1000.,), BOTH, (0., 1.), (20.,), [4.0],
    expect=either(BASE + ["surface_resonance_unresolved", "d3_ill_conditioned"]), routine=True, ref=True)
for t in (0.01, 0.02, 0.1):
    fam("EI5", "EI", ("Al",), (1000.,), ("k_total",), (1.,), (20.,), [t], expect=either(BASE), routine=True, ref=True)


# ---- RJ specific rejection before any panel of the body
def rj(exception, reason=None, idx=None, **kw):
    return {"status": "rejected", "exception": exception, "reason": reason, "offending_indices": idx,
            "panel_evaluations": 0, **kw}


fam("RJ1", "RJ", ("Si",), (1000.,), BOTH, (1.,), (20.,), "GA10G", expect=rj("DielectricGapEdgeError", idx=[10]),
    routine=True)
fam("RJ1", "RJ", ("Si",), (1000.,), BOTH, (1.,), (20.,), [G + 1.12e-12 / 2],
    expect=rj("DielectricGapEdgeError", idx=[0]), routine=True)
for m, t in (("P1", 1.999), ("P1", 1.9), ("P0", 1.999)):
    fam("RJ2", "RJ", (m,), (1000.,), BOTH, (1.,), (20.,), [t],
        expect=rj("FixedDepthGapDomainError", "surface_pole_not_excluded", [0]), routine=True)
fam("RJ2", "RJ", ("M1",), (1000.,), BOTH, (1.,), (20.,), [2.5],
    expect=rj("FixedDepthGapDomainError", "monotonicity_condition_failed", [0]), routine=True)
fam("RJ2", "RJ", ("UE",), (3.,), ("k_parallel",), (1.,), (20.,), [1.5],
    expect=rj("FixedDepthGapDomainError", "surface_pole_not_excluded", [0], U_eV_present_and_null=True),
    routine=True)
fam("RJ2-SPedge", "RJ", ("SP",), (1000.,), BOTH, (1.,), (20.,), [1.0], expect=rj("DielectricGapEdgeError", idx=[0]),
    routine=True)
fam("RJ2-SPsign", "RJ", ("SP",), (1000.,), BOTH, (1.,), (20.,), [2.5],
    expect=rj("ValueError", label="static_sign", message_contains=["static sign tolerance"]), routine=True)
fam("RJ3", "RJ", ("P1",), (1000.,), ("k_total",), (1.,), (20.,), [1.999], entry="CLI",
    grid={"min": 1.999, "max": 1.999, "step": 1.0},
    expect={"status": "error", "exit_code": 2, "output_dir_after": "absent", "exception": "FixedDepthGapDomainError",
            "reason": "surface_pole_not_excluded", "offending_indices": [0], "panel_evaluations": 0},
    routine=True)
fam("RJ4", "RJ", ("Al",), (1000.,), ("k_total",), (1.,), BAD_LENGTHS, [5.0],
    expect=rj("ValueError", label="invalid_weight_length", message_contains=["weight_length_A"]), routine=True)
fam("RJ4", "RJ", ("Al",), (1000.,), ("k_total",), (1.,), (20.,), [5.0], theta=10.0,
    expect=rj("ValueError", label="unsupported_theta", message_contains=["theta_deg"]), routine=True)
fam("RJ4", "RJ", ("Al",), (1000.,), ("k_total",), (1.,), (20.,), [5.0], settings="LIM1",
    expect=rj("ValueError", label="limit_too_small",
              message_contains=["integration limit too small", "loss_eV=5, index=(0,)"]), routine=True)
fam("RJ5", "RJ", ("Al",), (1000.,), ("k_total",), (1.,), (20.,), "GC", entry="CLI", grid=GRID_GC, outdir="existing",
    outdir_name="ja",
    expect={"status": "error", "exit_code": 2, "output_dir_after": "unchanged", "exception": "FileExistsError",
            "offending_indices": None, "panel_evaluations": 0},
    routine=True)


# ---------------------------------------------------------------- checks
def point_key(r, t):
    return (r["material"], r["electron_energy_eV"], r["epsilon_momentum"], r["hole_charge"]["hex"],
            json.dumps(r["weight_length_A"], sort_keys=True), h(t), r["settings"], r["theta_deg"], r["entry"],
            r["weight_length_source"], r["output_dir_before"], r.get("output_dir_name"))


run_keys = {}
for r in runs:
    rk = point_key(r, 0.0)[:5] + point_key(r, 0.0)[6:] + (tuple(h(t) for t in r["_values"]),)
    need(rk not in run_keys, f"duplicate run key: {r['id']} and {run_keys.get(rk)}")
    run_keys[rk] = r["id"]

point_cat = {}
for r in runs:
    idx = r["expect"].get("offending_indices")
    n = len(r["_values"])
    if r["category"] == "RJ" and idx is not None:
        need(len(idx) > 0 and len(set(idx)) == len(idx) and all(0 <= i < n for i in idx),
             f"{r['id']}: bad offending indices {idx}")
    for i, t in enumerate(r["_values"]):
        if r["category"] == "RJ" and idx is not None and i not in idx:
            continue  # not evaluated: the run is rejected before any panel; indexless rejections count all
        key = point_key(r, t)
        old = point_cat.setdefault(key, (r["category"], r["id"]))
        need(old[0] == r["category"], f"point key in two categories: {key} {old} {r['id']}")

for r in runs:
    if r["category"] == "RJ":
        continue
    cats = {point_cat[point_key(r, t)][0] for t in r["_values"]}
    need(cats == {r["category"]}, f"{r['id']}: run mixes categories {cats}")
    if r["category"] == "EI":
        reasons = r["expect"]["allowed_reasons"]
        need(len(r["_values"]) == 1 and reasons and set(reasons) <= set(REASONS), f"{r['id']}: bad EI run")
    if r["category"] == "NC":
        need(r["expect"]["reason"] in REASONS, f"{r['id']}: unknown reason")

# ---------------------------------------------------------------- reference keys (V7, V11) and limit keys (V4, V5)
refs = {}
for r in runs:
    for k in r.get("reference_keys", []):
        refs.setdefault(k, set()).add(r["family"])
for r in runs:
    for k in r.get("term_13a_reference_keys", []):
        need(k in refs, f"{r['id']}: 13a cross-length key {k} has no reference")
need(all("reference_keys" in r for r in runs if r["category"] == "EI"), "EI run without a reference key")
reference_keys = []
for k in sorted(refs):
    m, e, mo, la, t = k.split("|")
    reference_keys.append({"key": k, "material": m, "electron_energy_eV": float(e), "epsilon_momentum": mo,
                           "weight_length_A": float.fromhex(la), "loss_eV": float.fromhex(t),
                           "used_by": sorted(refs[k])})
limits = {}
for r in runs:
    if r["family"] in ("MS8", "MS14", "NC3"):
        kind = "lambda_to_0"
    elif r["family"] in ("MS9", "MS15"):
        kind = "lambda_to_inf"
    else:
        continue
    for t in r["_values"]:
        k = f"{kind}|{r['material']}|{r['electron_energy_eV']!r}|{r['epsilon_momentum']}|{r['hole_charge']['hex']}|{h(t)}"
        limits.setdefault(k, set()).add(r["weight_length_A"]["value"])
limit_keys = []
for k in sorted(limits):
    kind, m, e, mo, q, t = k.split("|")
    limit_keys.append({"key": k, "kind": kind, "material": m, "electron_energy_eV": float(e), "epsilon_momentum": mo,
                       "hole_charge": float.fromhex(q), "loss_eV": float.fromhex(t),
                       "weight_lengths_A": sorted(limits[k])})

# ---------------------------------------------------------------- independent expectations (not derived from fam calls)
EI_REASONS = {
    "EI1": {"budget_not_met", "panel_floor_exceeds_budget", "direct_sum_floor_exceeds_budget"},
    "EI2": {"budget_not_met", "panel_floor_exceeds_budget"},
    "EI3": {"budget_not_met", "panel_floor_exceeds_budget", "surface_denominator_unresolved"},
    "EI4": {"budget_not_met", "panel_floor_exceeds_budget", "surface_resonance_unresolved", "d3_ill_conditioned"},
    "EI5": {"budget_not_met", "panel_floor_exceeds_budget"},
}
INDEX_REQUIRED = {"RJ1", "RJ2", "RJ2-SPedge", "RJ3"}
INDEX_ABSENT = {"RJ2-SPsign", "RJ4", "RJ5"}


def active(m):
    return any(line[1] != 0 for line in MATERIALS[m]["oscillators"])


def subgap(m, t):
    g = MATERIALS[m]["gap_eV"]
    return g > 0 and 0 < t < g


def expected_ref(r, t):
    """Physical condition for a V7 reference at (run, loss), written independently of the fam(ref=...) flags."""
    f, m, la = r["family"], r["material"], r["weight_length_A"].get("value")
    if f in ("MS5", "MS6", "MS11", "MS13", "MS16", "MS17") or f.startswith("EI"):
        return True
    if f == "MS14":
        return subgap(m, t)
    if f == "MS15":
        return t == 1.0
    return False


for r in runs:
    f = r["family"]
    if f in EI_REASONS:
        need(set(r["expect"]["allowed_reasons"]) == EI_REASONS[f], f"{r['id']}: EI reasons differ from the table")
    if f in INDEX_REQUIRED:
        need(r["expect"].get("offending_indices") is not None, f"{r['id']}: rejection indices required")
    if f in INDEX_ABSENT:
        need(r["expect"].get("offending_indices") is None, f"{r['id']}: rejection without indices")
    for t, has in zip(r["_values"], [True] * len(r.get("reference_keys", [])) or [False] * len(r["_values"])):
        need(expected_ref(r, t) == has, f"{r['id']}: reference key presence differs from the physical rule at {t!r}")

# ---------------------------------------------------------------- check plan: (run id, original index, key)
plan = {"V3d": [], "V4": [], "V5": [], "V7": [], "V11iii": [], "V11v": []}
for r in runs:
    if r["category"] not in ("MS", "EI"):
        continue
    m, q = r["material"], r["hole_charge"]["value"]
    la = r["weight_length_A"].get("value")
    for i, t in enumerate(r["_values"]):
        rk = ref_key(m, r["electron_energy_eV"], r["epsilon_momentum"], la, t)
        has = rk in refs
        if "reference_keys" in r:
            plan["V7"].append([r["id"], i, r["reference_keys"][i]])
        if "term_13a_reference_keys" in r:
            plan["V3d"].append([r["id"], i, r["term_13a_reference_keys"][i]])
        if r["family"] in ("MS8", "MS14"):
            plan["V4"].append([r["id"], i, f"lambda_to_0|{m}|{r['electron_energy_eV']!r}|{r['epsilon_momentum']}|{r['hole_charge']['hex']}|{h(t)}"])
        if r["family"] in ("MS9", "MS15"):
            plan["V5"].append([r["id"], i, f"lambda_to_inf|{m}|{r['electron_energy_eV']!r}|{r['epsilon_momentum']}|{r['hole_charge']['hex']}|{h(t)}"])
        if subgap(m, t):
            if q == 1.0 and active(m):
                plan["V11iii"].append([r["id"], i, rk if has else None])
            elif q != 1.0:
                need(has or not active(m), f"{r['id']}: sub-gap non-unit charge without a reference at {t!r}")
                plan["V11v"].append([r["id"], i, rk if active(m) else "analytic_zero"])
# Each linked key must carry the run's own physical fields at that original index (13a links use 20 A).
ref_fields = {k["key"]: k for k in reference_keys}
for r in runs:
    for field, length in (("reference_keys", None), ("term_13a_reference_keys", 20.0)):
        for i, k in enumerate(r.get(field, [])):
            f = ref_fields[k]
            la = r["weight_length_A"]["value"] if length is None else length
            need(f["material"] == r["material"] and f["electron_energy_eV"] == r["electron_energy_eV"]
                 and f["epsilon_momentum"] == r["epsilon_momentum"] and f["weight_length_A"] == la
                 and h(f["loss_eV"]) == h(r["_values"][i]), f"{r['id']}[{i}]: {field} fields differ from the run")
limit_names = {k["key"] for k in limit_keys}
need(all(e[2] in limit_names for e in plan["V4"] + plan["V5"]), "a V4/V5 plan entry has no limit key")
routine_ids = {r["id"] for r in runs if r["routine"]}
plan_counts = {k: {"full": len(v), "routine": sum(1 for e in v if e[0] in routine_ids)} for k, v in plan.items()}

for r in runs:
    del r["_values"]
counts = {c: sum(1 for v in point_cat.values() if v[0] == c) for c in ("MS", "NC", "EI", "RJ")}
doc = {
    "schema": "xpsloss-m3a-expected-states-v2",
    "prereg": "docs/M3A_PREREG.md",
    "materials": MATERIALS,
    "arrays": {k: {"expression": e, "length": int(len(a)),
                   "sha256_float64_le": hashlib.sha256(np.asarray(a, "<f8").tobytes()).hexdigest(),
                   "first_hex": h(a[0]), "last_hex": h(a[-1])} for k, (e, a) in ARRAYS.items()},
    "settings": {k: {"epsabs": v[0], "epsrel": v[1], "limit": v[2]} for k, v in SETTINGS.items()},
    "reason_codes": list(REASONS),
    "cli_strings": {"ja_source": JA_SOURCE, "ja_output_dir": JA_DIR},
    "counts": {"runs": len(runs), "routine_runs": sum(r["routine"] for r in runs), "point_keys": len(point_cat),
               "point_keys_by_category": counts, "reference_keys": len(reference_keys),
               "reference_point_keys": sum(len(r.get("reference_keys", [])) for r in runs),
               "routine_reference_keys": len({k for r in runs if r["routine"] for k in r.get("reference_keys", [])}),
               "limit_keys": len(limit_keys)},
    "runs": runs,
    "reference_keys": reference_keys,
    "limit_keys": limit_keys,
    "check_plan": plan,
    "check_plan_counts": plan_counts,
}
text = json.dumps(doc, ensure_ascii=False, indent=1, sort_keys=True) + "\n"
data = text.encode("utf-8")
if sys.argv[1:] == ["--check"]:
    same = OUT.read_bytes() == data
    print(json.dumps({"identical": same, "sha256": hashlib.sha256(data).hexdigest()}))
    sys.exit(0 if same else 1)
OUT.write_bytes(data)
print(json.dumps({"sha256": hashlib.sha256(data).hexdigest(), **doc["counts"]}, indent=1))
