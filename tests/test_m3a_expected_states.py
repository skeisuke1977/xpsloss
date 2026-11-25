"""M3a expected states (docs/M3A_PREREG.md sections 3, 4 and 6): V13, the MS judgment and the
R2 part of V15, written by R2 without the implementation.

The table benchmarks/numerics/m3a/expected_states_m3a.json (SHA-256 checked first) lists every
run.  The routine runs (flag "routine", 1329 runs) are executed through the API or the CLI and
judged by category:
  MS  success; every point and all eleven outputs have finite, nonnegative estimates with
      e_j <= B_run_j = max(epsabs, epsrel |K_j|); the extra conditions of "expect".  For every MS
      run the diagnostics error_budget_ratios, max_budget_ratio (tie rule of section 4) and
      quadrature_error are recomputed from the returned values (V15, R2 part).
  RJ  exception type, reason, offending indices, zero panel evaluations of the body.
  NC  DepthAverageConvergenceError with reason, loss_eV, index (column, subdivisions when given)
      and the failure-side observations of section 4 (recorded _panel_integrals calls).
  EI  success as MS, or an allowed reason at index 0 with the failure-side observations.
      The reference comparison of a successful EI run belongs to R1's tests.
  CLI state, exit code, output directory afterwards, stderr JSON by error type.
Before the runs the planned pair counts of check_plan are checked against check_plan_counts;
after the runs each planned pair is classified (value available, allowed EI non-convergence,
run failed, not executed) and the classes must add up to the plan.  A missing pair is not a pass.
Set XPSLOSS_M3A_FULL=1 to execute all 2041 runs.  benchmarks/numerics/m3a/run_expected_states_m3a.py
uses the same functions for the full record.
"""
from contextlib import ExitStack, redirect_stderr, redirect_stdout
import copy
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TABLE_FILE = ROOT / "benchmarks" / "numerics" / "m3a" / "expected_states_m3a.json"
TABLE_SHA256 = "d87dad8f02a2fd211679229bd7c5205d908dbb2b3b241b30a27e94469f414ed5"
INPUTS = {"Al": ROOT / "benchmarks" / "inputs" / "al_elfdb_v1.0.bulk.json",
          "Si": ROOT / "benchmarks" / "inputs" / "si_elfdb_v1.0.bulk.json"}
FULL = os.environ.get("XPSLOSS_M3A_FULL") == "1"
TERMS = ("term_13a", "term_13b", "term_14", "term_15_T1", "term_15_T2", "term_15_T3", "term_16_e", "term_16_h")
OUTPUTS = TERMS + ("no_hole", "hole", "total")
UR = 2.0 ** -53
INVALID_LENGTHS = {"0.0": 0.0, "-1.0": -1.0, "nan": float("nan"), "inf": float("inf"), "1+0j": 1 + 0j, "True": True}
CSV_COLUMNS = ("loss_eV", "total", "no_hole", "hole") + TERMS + ("quadrature_error",)
NODE_REASON_TEXT = {
    "surface_denominator_unresolved": ("surface denominator unresolved at quadrature node",),
    "oscillator_denominator_unresolved": ("oscillator denominator unresolved at quadrature node",),
    "linewidth_unresolved": ("linewidth unresolved on momentum panels", "surface oscillator sum unresolved",
                             "unresolved oscillator and surface lines"),
    "d3_ill_conditioned": ("closed-form D3",),
    "surface_resonance_unresolved": ("surface search", "surface resonance unresolved"),
}
FLOOR_REASONS = ("panel_floor_exceeds_budget", "direct_sum_floor_exceeds_budget")
GAP_WARNING = re.compile(r"\bgap\b", re.IGNORECASE)
LOW_LOSS_WARNING = re.compile(r"low[- ]?loss|lower (loss )?limit|loss window|window area", re.IGNORECASE)


# ----------------------------------------------------------------------------- table

def load_table():
    raw = TABLE_FILE.read_bytes()
    sha = hashlib.sha256(raw).hexdigest()
    if sha != TABLE_SHA256:
        raise AssertionError(f"expected-state table SHA-256 {sha} differs from {TABLE_SHA256}")
    return json.loads(raw.decode("utf-8"))


def _arrays(table):
    ga = 0.2 * np.arange(1, 500)
    f = 0.2 / 4 * np.arange(4, 1997)
    out = {"GA": ga, "GC": np.linspace(0.2, 99.8, 499), "F": f, "F2": f[::2], "FRD": np.r_[f[::-4], f[:40:4]],
           "GA10G": np.r_[ga[:10], 1.12]}
    for name, arr in out.items():
        meta = table["arrays"][name]
        sha = hashlib.sha256(np.asarray(arr, "<f8").tobytes()).hexdigest()
        if sha != meta["sha256_float64_le"] or len(arr) != meta["length"]:
            raise AssertionError(f"array {name} differs from the table")
    return out


_ARRAY_CACHE = {}


def run_losses(run, table):
    spec = run["losses"]
    if "array" in spec:
        if not _ARRAY_CACHE:
            _ARRAY_CACHE.update(_arrays(table))
        return np.array(_ARRAY_CACHE[spec["array"]], dtype=float)
    return np.array([float.fromhex(h) for h in spec["values_hex"]], dtype=float)


def material(name, table):
    from xpsloss import Material, Oscillator
    entry = table["materials"][name]
    return Material(name, float(entry["gap_eV"]), tuple(Oscillator(*map(float, line)) for line in entry["oscillators"]))


def effective_settings(run, table):
    s = table["settings"][run["settings"]]
    return {"epsabs": s["epsabs"], "epsrel": s["epsrel"], "limit": s["limit"]}


def weight_length(run):
    lam = run["weight_length_A"]
    if "invalid" in lam:
        return INVALID_LENGTHS[lam["invalid"]]
    return float.fromhex(lam["hex"])


def hole_charge(run):
    return float.fromhex(run["hole_charge"]["hex"])


def selected_runs(table, full=FULL):
    return [r for r in table["runs"] if full or r["routine"]]


def check_plan_counts(table):
    """Planned pairs per check: lengths and routine counts against check_plan_counts (before any run)."""
    problems = []
    ids = {r["id"]: r for r in table["runs"]}
    for check, entries in table["check_plan"].items():
        counts = table["check_plan_counts"][check]
        routine = sum(1 for e in entries if ids[e[0]]["routine"])
        if len(entries) != counts["full"]:
            problems.append(f"{check}: {len(entries)} pairs, table says {counts['full']}")
        if routine != counts["routine"]:
            problems.append(f"{check}: {routine} routine pairs, table says {counts['routine']}")
        for run_id, index, _ in entries:
            run = ids.get(run_id)
            if run is None:
                problems.append(f"{check}: unknown run {run_id}")
            elif not 0 <= index < len(run_losses(run, table)):
                problems.append(f"{check}: index {index} outside {run_id}")
        seen = set()
        for run_id, index, key in entries:
            if (run_id, index) in seen:
                problems.append(f"{check}: duplicate pair {run_id}[{index}]")
            seen.add((run_id, index))
    return problems


# ----------------------------------------------------------------------------- recording the body

class Recorder:
    """Counts and records the body's panel evaluations and the surface search (section 4)."""

    def __init__(self):
        self.calls = []
        self.surface_calls = 0
        self.raised_in_panels = False
        self.raised_in_surface = False

    def patches(self):
        import xpsloss.xps_depth_averaged as dav
        import xpsloss._gap_response as gap_response
        original = dav._panel_integrals
        recorder = self

        def panel_integrals(*args, **kwargs):
            try:
                out = original(*args, **kwargs)
            except BaseException:
                recorder.raised_in_panels = True
                recorder.calls.append({"args": args, "kwargs": kwargs, "out": None, "raised": True})
                raise
            recorder.calls.append({"args": args, "kwargs": kwargs, "out": out, "raised": False})
            return out

        surface_original = gap_response.resolve_surface

        def resolve_surface(*args, **kwargs):
            recorder.surface_calls += 1
            try:
                return surface_original(*args, **kwargs)
            except BaseException:
                recorder.raised_in_surface = True
                raise

        stack = [patch.object(dav, "_panel_integrals", panel_integrals),
                 patch.object(gap_response, "resolve_surface", resolve_surface)]
        for module_name in ("xpsloss.xps_depth_averaged", "xpsloss.xps_fixed_depth"):
            module = sys.modules.get(module_name)
            if module is not None and getattr(module, "resolve_surface", None) is surface_original:
                stack.append(patch.object(module, "resolve_surface", resolve_surface))
        return stack

    @property
    def count(self):
        return len(self.calls)

    @staticmethod
    def _values(out):
        if isinstance(out, tuple) and len(out) >= 2:
            values, estimates = np.asarray(out[0], dtype=float), np.asarray(out[1], dtype=float)
            if values.ndim == 2 and values.shape[1] == 11 and estimates.shape == values.shape:
                return values, estimates
        return None, None

    @staticmethod
    def _edges(call, n):
        arrays = [a for a in list(call["args"]) + list(call["kwargs"].values())
                  if isinstance(a, np.ndarray) and a.ndim == 1 and a.dtype.kind == "f"]
        pairs = [a for a in arrays if a.size == n]
        if len(pairs) >= 2:
            return list(zip(pairs[0].tolist(), pairs[1].tolist()))
        joined = [a for a in arrays if a.size == n + 1]
        if joined:
            e = joined[0].tolist()
            return list(zip(e[:-1], e[1:]))
        return None

    def panels(self):
        """Recorded panels: list of dicts per call with edges (or None), values and estimates."""
        out = []
        for call in self.calls:
            values, estimates = self._values(call["out"]) if not call["raised"] else (None, None)
            edges = None if values is None else self._edges(call, values.shape[0])
            out.append({"edges": edges, "values": values, "estimates": estimates, "raised": call["raised"]})
        return out

    def leaves(self):
        """Final mesh reconstructed from the calls: latest record per panel, without split parents.
        Returns None when the panel edges cannot be read from the call arguments."""
        latest = {}
        for rec in self.panels():
            if rec["values"] is None:
                continue
            if rec["edges"] is None:
                return None
            for p, edge in enumerate(rec["edges"]):
                latest[edge] = (rec["values"][p], rec["estimates"][p])
        keys = list(latest)
        leaves = [k for k in keys if not any(o != k and k[0] <= o[0] and o[1] <= k[1] for o in keys)]
        return {k: latest[k] for k in leaves}


# ----------------------------------------------------------------------------- execution

def _api(run, table, recorder):
    from xpsloss.xps_depth_averaged import effective_cross_section_depth_averaged
    kwargs = dict(weight_length_A=weight_length(run), epsilon_momentum=run["epsilon_momentum"],
                  theta_deg=run["theta_deg"], hole_charge=hole_charge(run))
    if run["settings"] != "APIDEFAULT":
        kwargs.update(effective_settings(run, table))
    with ExitStack() as stack:
        for p in recorder.patches():
            stack.enter_context(p)
        try:
            return {"result": effective_cross_section_depth_averaged(
                material(run["material"], table), run["electron_energy_eV"], run_losses(run, table), **kwargs)}
        except Exception as exc:  # noqa: BLE001 - every exception is judged against the table
            return {"exception": exc}


def cli_config(run, table, folder):
    """Input config: the material's input file (Al, Si) or the table material, with the run's
    loss_grid_eV and integration.  Returns the path and the problems found when reading it back."""
    from xpsloss.config import load_config
    name = run["material"]
    if name in INPUTS:
        base = json.loads(INPUTS[name].read_text("utf-8"))
        if base["material"]["name"] != name:
            raise AssertionError(f"input file for {name} names {base['material']['name']}")
    else:
        mat = table["materials"][name]
        base = {"model": "bulk_drude_lindhard_v0",
                "material": {"name": name, "gap_eV": mat["gap_eV"],
                             "oscillators": [dict(zip(("energy_eV", "strength_eV2", "width_eV", "dispersion"), line))
                                             for line in mat["oscillators"]]}}
    cfg = copy.deepcopy(base)
    cfg["electron_energy_eV"] = run["electron_energy_eV"]
    cfg["loss_grid_eV"] = dict(run["cli_config"]["loss_grid_eV"])
    cfg.pop("integration", None)
    if run["cli_config"]["integration"] is not None:
        cfg["integration"] = dict(table["settings"][run["cli_config"]["integration"]])
    path = Path(folder) / f"config_{re.sub(r'[^A-Za-z0-9]+', '_', run['id'])}.json"
    path.write_text(json.dumps(cfg, allow_nan=False), encoding="utf-8")
    problems = []
    loaded = load_config(path)
    grid = loaded.grid()
    want = run_losses(run, table)
    if grid.shape != want.shape or not np.array_equal(grid.view("<u8"), want.astype("<f8").view("<u8")):
        problems.append("config grid differs from the table losses")
    eff = effective_settings(run, table)
    if (loaded.epsabs, loaded.epsrel, loaded.limit) != (eff["epsabs"], eff["epsrel"], eff["limit"]):
        problems.append(f"config precision {(loaded.epsabs, loaded.epsrel, loaded.limit)} differs from {eff}")
    want_mat = material(name, table)
    if loaded.material.gap_eV != want_mat.gap_eV or [tuple(vars(o).values()) for o in loaded.material.oscillators] \
            != [tuple(vars(o).values()) for o in want_mat.oscillators]:
        problems.append("config material differs from the table")
    return path, problems


def _tree(folder):
    folder = Path(folder)
    if not folder.exists():
        return None
    return {str(p.relative_to(folder)): hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else "dir"
            for p in sorted(folder.rglob("*"))}


def _cli(run, table, recorder, workdir):
    from xpsloss import cli
    folder = Path(tempfile.mkdtemp(dir=workdir))
    config, problems = cli_config(run, table, folder)
    name = table["cli_strings"]["ja_output_dir"] if run.get("output_dir_name") == "ja" else "out_ascii"
    out = folder / name
    before = None
    if run["output_dir_before"] == "existing":
        out.mkdir()
        (out / "known.txt").write_bytes(b"known content for RJ5\n")
        before = _tree(out)
    argv = ["xps-depth-averaged", str(config), "--weight-length-A", repr(weight_length(run)),
            "--epsilon-momentum", run["epsilon_momentum"], "--theta-deg", repr(run["theta_deg"]),
            "--hole-charge", repr(hole_charge(run)), "--out", str(out)]
    if run["weight_length_source"] is not None:
        argv[4:4] = ["--weight-length-source", run["weight_length_source"]]
    stdout, stderr = io.StringIO(), io.StringIO()
    with ExitStack() as stack:
        for p in recorder.patches():
            stack.enter_context(p)
        with redirect_stdout(stdout), redirect_stderr(stderr):
            code = cli.main(argv)
    return {"cli": True, "exit_code": code, "stdout": stdout.getvalue(), "stderr": stderr.getvalue(),
            "out": out, "before": before, "after": _tree(out), "config_problems": problems}


# ----------------------------------------------------------------------------- judgments

def _rel_equal(a, b, rel=4 * UR):
    a, b = float(a), float(b)
    if a == b:
        return True
    return math.isfinite(a) and math.isfinite(b) and abs(a - b) <= rel * max(abs(a), abs(b))


def recompute_budget_fields(values, errors, settings):
    """error_budget_ratios, max_budget_ratio (tie rule) and quadrature_error from returned arrays.
    values, errors: dicts output name -> 1-D arrays."""
    eps_a, eps_r = settings["epsabs"], settings["epsrel"]
    ratios = {}
    for name in OUTPUTS:
        budget = np.maximum(eps_a, eps_r * np.abs(np.asarray(values[name], dtype=float)))
        ratios[name] = np.asarray(errors[name], dtype=float) / budget
    # The first maximum in output order, then in original index order (section 4).  With all
    # estimates zero this is {"value": 0.0, "column": "term_13a", "index": 0}.
    best = None
    for name in OUTPUTS:
        for i, r in enumerate(ratios[name]):
            if best is None or r > best["value"]:
                best = {"value": float(r), "column": name, "index": i}
    quad = np.max(np.array([np.asarray(errors[n], dtype=float) for n in OUTPUTS]), axis=0)
    return ratios, best, quad


def _check_budget_fields(values, errors, diagnostics, quadrature_error, settings, where):
    out = []
    ratios, best, quad = recompute_budget_fields(values, errors, settings)
    given = diagnostics.get("error_budget_ratios")
    if not isinstance(given, dict) or set(given) != set(OUTPUTS):
        return out + [f"{where}: error_budget_ratios missing or without the eleven outputs"]
    for name in OUTPUTS:
        g = np.asarray(given[name], dtype=float)
        if g.shape != ratios[name].shape or not all(_rel_equal(x, y) for x, y in zip(g, ratios[name])):
            out.append(f"{where}: error_budget_ratios[{name}] differs from e_j/B_run_j")
    mbr = diagnostics.get("max_budget_ratio")
    if not isinstance(mbr, dict) or set(mbr) != {"value", "column", "index"}:
        out.append(f"{where}: max_budget_ratio malformed {mbr!r}")
    elif not (_rel_equal(mbr["value"], best["value"]) and mbr["column"] == best["column"]
              and mbr["index"] == best["index"]):
        out.append(f"{where}: max_budget_ratio {mbr} differs from the recomputed {best}")
    q = np.asarray(quadrature_error, dtype=float)
    if q.shape != quad.shape or not np.array_equal(q, quad):
        out.append(f"{where}: quadrature_error is not the maximum of the eleven estimates")
    return out


def _check_ms_values(values, errors, settings, n, where):
    out = []
    eps_a, eps_r = settings["epsabs"], settings["epsrel"]
    for name in OUTPUTS:
        v = np.asarray(values.get(name), dtype=float) if values.get(name) is not None else None
        e = np.asarray(errors.get(name), dtype=float) if errors.get(name) is not None else None
        if v is None or e is None or v.shape != (n,) or e.shape != (n,):
            out.append(f"{where}: {name} values or component_errors missing or of the wrong length")
            continue
        if not np.all(np.isfinite(v)):
            out.append(f"{where}: {name} not finite at {np.flatnonzero(~np.isfinite(v)).tolist()[:5]}")
        if not np.all(np.isfinite(e) & (e >= 0)):
            out.append(f"{where}: {name} estimate not finite and nonnegative")
        bad = np.flatnonzero(~(e <= np.maximum(eps_a, eps_r * np.abs(v))))
        if bad.size:
            out.append(f"{where}: {name} e_j > B_run_j at {bad.tolist()[:5]}")
    return out


def _check_direct_sum(diagnostics, values, errors, expect, n, where):
    out = []
    ds = diagnostics.get("direct_sum")
    dmax = diagnostics.get("direct_sum_max", "missing")
    if not isinstance(ds, list) or len(ds) != n or dmax == "missing":
        return [f"{where}: direct_sum / direct_sum_max missing or of the wrong length"]
    used = [i for i, d in enumerate(ds) if d is not None]
    for i in used:
        d = ds[i]
        if not isinstance(d, dict) or set(d) != {"K_dir", "e_dir", "delta_rec", "e_total"}:
            out.append(f"{where}: direct_sum[{i}] malformed")
            continue
        nums = [d[k] for k in ("K_dir", "e_dir", "delta_rec", "e_total")]
        if not all(isinstance(x, (int, float)) and math.isfinite(x) for x in nums) or min(nums[1:]) < 0:
            out.append(f"{where}: direct_sum[{i}] not finite or a negative error")
            continue
        if d["e_total"] != float(errors["total"][i]):
            out.append(f"{where}: direct_sum[{i}].e_total differs from component_errors['total']")
        if d["delta_rec"] != abs(float(values["total"][i]) - d["K_dir"]):
            out.append(f"{where}: direct_sum[{i}].delta_rec differs from |total - K_dir|")
    if used:
        if not isinstance(dmax, dict):
            out.append(f"{where}: direct_sum_max is not a dict although direct sums were used")
    elif dmax is not None:
        out.append(f"{where}: direct_sum_max must be None without direct sums")
    want = expect.get("direct_sum_diagnostics")
    if want == "null" and used:
        out.append(f"{where}: direct sums used at {used[:5]} but expected none")
    if want == "nonnull" and not used:
        out.append(f"{where}: expected at least one direct sum")
    return out


def _result_maps(result):
    values = {name: np.asarray(getattr(result, name)) for name in OUTPUTS}
    errors = dict(getattr(result, "component_errors", {}) or {})
    return values, errors


def judge_ms_api(run, table, result, n, recorder):
    where = run["id"]
    out = []
    settings = effective_settings(run, table)
    losses = run_losses(run, table)
    got_loss = np.asarray(getattr(result, "loss_eV", None), dtype=float)
    if got_loss.shape != losses.shape or not np.array_equal(got_loss, losses):
        out.append(f"{where}: loss_eV differs from the input")
    for field in ("depth_A", "path_length_A"):
        if hasattr(result, field):
            out.append(f"{where}: result has {field}")
    values, errors = _result_maps(result)
    out += _check_ms_values(values, errors, settings, n, where)
    diagnostics = getattr(result, "diagnostics", None)
    if not isinstance(diagnostics, dict):
        return out + [f"{where}: diagnostics missing"]
    eff = diagnostics.get("effective_settings")
    if eff != settings:
        out.append(f"{where}: effective_settings {eff} differ from {settings}")
    if not out:
        out += _check_budget_fields(values, errors, diagnostics, result.quadrature_error, settings, where)
    out += _check_direct_sum(diagnostics, values, errors, run["expect"], n, where)
    expect = run["expect"]
    if expect.get("exact_zero"):
        nz = [name for name in OUTPUTS if np.any(values[name] != 0)]
        if nz:
            out.append(f"{where}: outputs not exactly zero: {nz}")
    return out


def _stderr_json(text):
    try:
        return json.loads(text)
    except ValueError:
        return None


def judge_ms_cli(run, table, res):
    where = run["id"]
    out = list(res["config_problems"])
    stdout = _stderr_json(res["stdout"])
    if res["exit_code"] != run["expect"].get("exit_code", 0) or not stdout or stdout.get("status") != "complete":
        return out + [f"{where}: CLI exit {res['exit_code']} stdout {res['stdout'][:300]!r} stderr {res['stderr'][:500]!r}"]
    if res["after"] is None or "run.json" not in res["after"] or "xps_depth_averaged.csv" not in res["after"]:
        return out + [f"{where}: output directory missing files: {res['after']}"]
    with open(res["out"] / "xps_depth_averaged.csv", newline="", encoding="ascii") as fh:
        rows = list(csv.reader(fh))
    if tuple(rows[0]) != CSV_COLUMNS:
        out.append(f"{where}: CSV header {rows[0]}")
        return out
    data = np.array([[float(x) for x in r] for r in rows[1:]])
    n = data.shape[0]
    losses = run_losses(run, table)
    if n != losses.size or not np.array_equal(data[:, 0], losses):
        out.append(f"{where}: CSV losses differ from the configured grid")
    manifest = json.loads((res["out"] / "run.json").read_text("utf-8"))
    diag = manifest.get("diagnostics", {})
    values = {name: data[:, CSV_COLUMNS.index(name)] for name in OUTPUTS}
    errs = diag.get("component_errors")
    if not isinstance(errs, dict) or set(errs) != set(OUTPUTS):
        return out + [f"{where}: run.json component_errors missing"]
    errors = {k: np.asarray(v, dtype=float) for k, v in errs.items()}
    settings = effective_settings(run, table)
    out += _check_ms_values(values, errors, settings, n, where)
    eff = diag.get("effective_settings")
    if eff != settings:
        out.append(f"{where}: run.json effective_settings {eff} differ from {settings}")
    out += _check_budget_fields(values, errors, diag, data[:, -1], settings, where)
    out += _check_direct_sum(diag, values, errors, run["expect"], n, where)
    weight = manifest.get("weight", {})
    if weight.get("length_source") != run["weight_length_source"]:
        out.append(f"{where}: weight.length_source {weight.get('length_source')!r}")
    return out


def _exc_name(exc):
    return type(exc).__name__


def _observe_failure(run, table, exc, recorder):
    """Failure-side observations of section 4 for DepthAverageConvergenceError."""
    where = run["id"]
    out = []
    reason = getattr(exc, "reason", None)
    settings = effective_settings(run, table)
    loss = run_losses(run, table)[0]
    if getattr(exc, "index", None) != 0 or getattr(exc, "loss_eV", None) != float(loss):
        out.append(f"{where}: exception loss/index {getattr(exc, 'loss_eV', None)!r}/{getattr(exc, 'index', None)!r}")
    for attr in ("reason", "loss_eV", "index", "column", "estimate", "budget", "subdivisions", "panels", "inner"):
        if not hasattr(exc, attr):
            out.append(f"{where}: exception has no attribute {attr}")
    if out:
        return out
    if not isinstance(exc.subdivisions, int) or not isinstance(exc.panels, int):
        out.append(f"{where}: subdivisions/panels not int")
    panels = recorder.panels()
    evaluated = [p for p in panels if p["values"] is not None]
    if reason in FLOOR_REASONS or reason == "budget_not_met":
        if not evaluated:
            return out + [f"{where}: {reason} without any recorded panel evaluation"]
        if exc.column not in OUTPUTS:
            out.append(f"{where}: column {exc.column!r}")
            return out
        c = OUTPUTS.index(exc.column)
        if not (isinstance(exc.estimate, float) and isinstance(exc.budget, float) and math.isfinite(exc.budget)
                and exc.estimate > exc.budget):
            out.append(f"{where}: estimate {exc.estimate!r} not above budget {exc.budget!r}")
    if reason in FLOOR_REASONS:
        last = evaluated[-1]
        total = math.fsum(last["values"][:, c].tolist())
        want = max(settings["epsabs"], settings["epsrel"] * abs(total))
        if not abs(exc.budget - want) <= 4 * UR * want:
            out.append(f"{where}: budget {exc.budget!r} differs from max(epsabs, epsrel|sum|) = {want!r}")
        if exc.subdivisions != 0:
            out.append(f"{where}: subdivisions {exc.subdivisions}")
        edges = [p["edges"] for p in evaluated]
        if all(e is not None for e in edges):
            if any(e != edges[0] for e in edges):
                out.append(f"{where}: panel sets other than the initial one were evaluated")
        elif len({p["values"].shape[0] for p in evaluated}) != 1:
            out.append(f"{where}: panel counts change between recorded calls")
        if exc.panels != evaluated[-1]["values"].shape[0]:
            out.append(f"{where}: panels {exc.panels} differ from the recorded initial set")
    elif reason == "budget_not_met":
        leaves = recorder.leaves()
        if leaves is None:
            out.append(f"{where}: panel edges not readable from the _panel_integrals arguments")
        else:
            limit = settings["limit"]
            if exc.panels != len(leaves) or exc.panels != limit:
                out.append(f"{where}: panels {exc.panels}, recorded final mesh {len(leaves)}, limit {limit}")
            total_est = math.fsum(float(v[1][c]) for v in leaves.values())
            if not exc.estimate >= total_est * (1 - 4 * UR):
                out.append(f"{where}: estimate {exc.estimate!r} below the recorded panel estimates {total_est!r}")
    elif reason in NODE_REASON_TEXT:
        if reason == "surface_resonance_unresolved":
            if not recorder.raised_in_surface:
                out.append(f"{where}: {reason} not raised inside resolve_surface")
        elif not recorder.raised_in_panels:
            out.append(f"{where}: {reason} not raised inside _panel_integrals")
        inner = exc.inner or ""
        if not any(text in inner for text in NODE_REASON_TEXT[reason]):
            out.append(f"{where}: inner {inner!r} lacks {NODE_REASON_TEXT[reason]}")
    else:
        out.append(f"{where}: reason {reason!r} has no registered failure-side observation")
    return out


def judge_api_failure(run, table, exc, recorder):
    where = run["id"]
    expect = run["expect"]
    out = []
    if _exc_name(exc) != expect["exception"]:
        return [f"{where}: {_exc_name(exc)}: {exc} (expected {expect['exception']})"]
    from xpsloss import NumericalConvergenceError
    if expect["exception"] == "DepthAverageConvergenceError":
        if not isinstance(exc, NumericalConvergenceError):
            out.append(f"{where}: not a NumericalConvergenceError")
        reason = getattr(exc, "reason", None)
        if run["category"] == "NC":
            if reason != expect["reason"]:
                out.append(f"{where}: reason {reason!r}, expected {expect['reason']!r}")
            for key in ("column", "subdivisions"):
                if key in expect and getattr(exc, key, None) != expect[key]:
                    out.append(f"{where}: {key} {getattr(exc, key, None)!r}, expected {expect[key]!r}")
        else:
            if reason not in expect["allowed_reasons"]:
                out.append(f"{where}: reason {reason!r} not in {expect['allowed_reasons']}")
        if getattr(exc, "index", None) != expect["index"]:
            out.append(f"{where}: index {getattr(exc, 'index', None)!r}")
        out += _observe_failure(run, table, exc, recorder)
        return out
    # RJ
    if recorder.count:
        out.append(f"{where}: {recorder.count} panel evaluations before the rejection")
    idx = expect.get("offending_indices")
    losses = run_losses(run, table)
    if expect["exception"] == "DielectricGapEdgeError":
        if list(getattr(exc, "offending_indices", [])) != idx:
            out.append(f"{where}: offending_indices {getattr(exc, 'offending_indices', None)} != {idx}")
    elif expect["exception"] == "FixedDepthGapDomainError":
        if getattr(exc, "reason", None) != expect["reason"]:
            out.append(f"{where}: reason {getattr(exc, 'reason', None)!r}")
        points = getattr(exc, "details", {}).get("points", [])
        if [p.get("index") for p in points] != idx:
            out.append(f"{where}: details points indices {[p.get('index') for p in points]} != {idx}")
        if list(getattr(exc, "offending_loss_eV", [])) != [float(losses[i]) for i in idx]:
            out.append(f"{where}: offending_loss_eV {getattr(exc, 'offending_loss_eV', None)}")
        if expect.get("U_eV_present_and_null") and not (points and "U_eV" in points[0] and points[0]["U_eV"] is None):
            out.append(f"{where}: details points[0] has no U_eV = None")
    elif expect["exception"] == "ValueError":
        for text in expect.get("message_contains", []):
            if text not in str(exc):
                out.append(f"{where}: message {str(exc)!r} lacks {text!r}")
    return out


def judge_cli_failure(run, table, res, recorder):
    where = run["id"]
    expect = run["expect"]
    out = list(res["config_problems"])
    if res["exit_code"] != expect["exit_code"]:
        out.append(f"{where}: exit code {res['exit_code']}")
    err = _stderr_json(res["stderr"])
    if err is None:
        return out + [f"{where}: stderr is not JSON: {res['stderr'][:300]!r}"]
    if not res["stderr"].isascii() or res["stdout"].strip():
        out.append(f"{where}: stderr not ASCII or unexpected stdout")
    if err.get("status") != expect["status"]:
        out.append(f"{where}: status {err.get('status')!r}, expected {expect['status']!r}")
    after = expect["output_dir_after"]
    if after == "absent" and res["after"] is not None:
        out.append(f"{where}: output directory was created")
    if after == "unchanged" and res["after"] != res["before"]:
        out.append(f"{where}: existing output directory changed")
    if recorder.count and run["category"] == "RJ":
        out.append(f"{where}: {recorder.count} panel evaluations before the rejection")
    if run["category"] == "NC":
        if err.get("error_type") != "DepthAverageConvergenceError" or err.get("reason") != expect["reason"]:
            out.append(f"{where}: stderr {err}")
        if err.get("index") != 0 or err.get("loss_eV") != float(run_losses(run, table)[0]):
            out.append(f"{where}: stderr loss/index {err.get('loss_eV')!r}/{err.get('index')!r}")
    elif expect.get("exception") == "FixedDepthGapDomainError":
        if err.get("error_type") != "FixedDepthGapDomainError" or err.get("reason") != expect["reason"]:
            out.append(f"{where}: stderr {err}")
        points = (err.get("details") or {}).get("points", [])
        if [p.get("index") for p in points] != expect["offending_indices"]:
            out.append(f"{where}: stderr indices {[p.get('index') for p in points]}")
        losses = run_losses(run, table)
        if err.get("offending_loss_eV") != [float(losses[i]) for i in expect["offending_indices"]]:
            out.append(f"{where}: stderr offending_loss_eV {err.get('offending_loss_eV')}")
    elif expect.get("exception") == "FileExistsError":
        if err.get("error_type") != "FileExistsError":
            out.append(f"{where}: stderr {err}")
        for key in ("reason", "offending_loss_eV", "offending_indices", "index", "loss_eV"):
            if key in err:
                out.append(f"{where}: FileExistsError stderr carries {key}")
    return out


def execute_and_judge(run, table, workdir):
    """Run one table entry and judge it.  Returns a picklable record."""
    recorder = Recorder()
    start = time.perf_counter()
    n = len(run_losses(run, table))
    record = {"id": run["id"], "family": run["family"], "category": run["category"], "entry": run["entry"],
              "routine": run["routine"], "failures": [], "outcome": None, "max_budget_ratio": None, "facts": {}}
    try:
        if run["entry"] == "CLI":
            res = _cli(run, table, recorder, workdir)
            if run["category"] == "MS":
                record["failures"] = judge_ms_cli(run, table, res)
                record["outcome"] = "success" if not record["failures"] else "failed"
                if res["after"] and "run.json" in res["after"]:
                    manifest = json.loads((res["out"] / "run.json").read_text("utf-8"))
                    record["facts"]["warnings"] = list(manifest.get("warnings", []))
                    record["max_budget_ratio"] = manifest.get("diagnostics", {}).get("max_budget_ratio")
            else:
                record["failures"] = judge_cli_failure(run, table, res, recorder)
                record["outcome"] = "expected_failure" if not record["failures"] else "failed"
        else:
            res = _api(run, table, recorder)
            if "result" in res:
                result = res["result"]
                if run["category"] in ("MS", "EI"):
                    record["failures"] = judge_ms_api(run, table, result, n, recorder)
                    record["outcome"] = "success" if not record["failures"] else "failed"
                    record["max_budget_ratio"] = (result.diagnostics or {}).get("max_budget_ratio")
                    record["facts"]["warnings"] = list(result.warnings)
                    record["facts"]["panel_evaluations"] = recorder.count
                    diag = result.diagnostics or {}
                    record["facts"]["gap"] = {k: diag.get(k, "missing") for k in
                                              ("subgap_proofs", "monotonicity", "static_sign")} \
                        if table["materials"][run["material"]]["gap_eV"] else None
                else:
                    record["failures"] = [f"{run['id']}: returned a result, expected {run['expect']['exception']}"]
                    record["outcome"] = "failed"
            else:
                exc = res["exception"]
                if run["category"] == "MS":
                    record["failures"] = [f"{run['id']}: {type(exc).__name__}: {str(exc)[:500]}"]
                    record["outcome"] = "failed"
                else:
                    record["failures"] = judge_api_failure(run, table, exc, recorder)
                    if record["failures"]:
                        record["outcome"] = "failed"
                    else:
                        record["outcome"] = "ei_not_converged" if run["category"] == "EI" else "expected_failure"
                    record["facts"]["reason"] = getattr(exc, "reason", None)
                    record["facts"]["exception"] = type(exc).__name__
    except Exception as exc:  # noqa: BLE001 - an error in the harness or an import error is a failure
        record["failures"].append(f"{run['id']}: {type(exc).__name__}: {str(exc)[:500]}")
        record["outcome"] = "failed"
    record["seconds"] = time.perf_counter() - start
    return record


def cross_run_checks(records, table):
    """Conditions that relate two runs: Z1 accepted like Si, Z0 without gap warnings, Z (CLI) without
    the low-loss warnings of the active-line Al CLI run, and a success run that calls the body."""
    out = []
    runs = {r["id"]: r for r in table["runs"]}
    by_key = {}
    for rec in records.values():
        run = runs[rec["id"]]
        key = (run["family"], run["material"], run["epsilon_momentum"], run["hole_charge"]["hex"],
               json.dumps(run["weight_length_A"], sort_keys=True), json.dumps(run["losses"], sort_keys=True),
               run["settings"], run["entry"], run["weight_length_source"])
        by_key[key] = rec
    for rec in records.values():
        run = runs[rec["id"]]
        expect = run["expect"]
        if rec["outcome"] != "success":
            continue
        if expect.get("same_gap_acceptance_as"):
            other = by_key.get(("MS10", expect["same_gap_acceptance_as"], run["epsilon_momentum"],
                                run["hole_charge"]["hex"], json.dumps(run["weight_length_A"], sort_keys=True),
                                json.dumps(run["losses"], sort_keys=True), run["settings"], "API", None))
            if other is None or other["outcome"] != "success":
                out.append(f"{rec['id']}: the Si run of the same key was not executed successfully")
            elif rec["facts"]["gap"] != other["facts"]["gap"] or "missing" in rec["facts"]["gap"].values():
                out.append(f"{rec['id']}: gap acceptance differs from {other['id']}")
        if expect.get("gap_warnings") is False:
            bad = [w for w in rec["facts"]["warnings"] if GAP_WARNING.search(w)]
            if bad:
                out.append(f"{rec['id']}: gap warnings {bad}")
        if expect.get("low_loss_warning") is False:
            al = [r for r in records.values() if runs[r["id"]]["family"] == "MS4" and runs[r["id"]]["entry"] == "CLI"
                  and runs[r["id"]]["material"] == "Al" and runs[r["id"]]["settings"] == "DEF"
                  and runs[r["id"]]["weight_length_source"] is None and r["outcome"] == "success"]
            if not al:
                out.append(f"{rec['id']}: no successful Al CLI control run")
            else:
                extra = set(al[0]["facts"]["warnings"]) - set(rec["facts"]["warnings"])
                if not extra:
                    out.append(f"{rec['id']}: the Al control has no warning that the zero-strength run lacks")
                bad = [w for w in rec["facts"]["warnings"] if LOW_LOSS_WARNING.search(w)]
                if bad:
                    out.append(f"{rec['id']}: low-loss warnings {bad}")
    api_ms = [r for r in records.values() if r["category"] == "MS" and r["entry"] == "API" and r["outcome"] == "success"]
    if api_ms and not any(r["facts"].get("panel_evaluations", 0) > 0 for r in api_ms):
        out.append("no successful API run evaluated xps_depth_averaged._panel_integrals (patched function not used)")
    return out


def plan_accounting(records, table, full=FULL):
    """Classify every planned pair of the executed subset; the classes must add up to the plan."""
    runs = {r["id"]: r for r in table["runs"]}
    summary, problems = {}, []
    for check, entries in table["check_plan"].items():
        cls = {"value": 0, "ei_not_converged": 0, "run_failed": 0, "not_executed": 0}
        planned = [e for e in entries if full or runs[e[0]]["routine"]]
        for run_id, index, _ in planned:
            rec = records.get(run_id)
            if rec is None:
                cls["not_executed"] += 1
            elif rec["outcome"] == "success":
                cls["value"] += 1
            elif rec["outcome"] == "ei_not_converged":
                cls["ei_not_converged"] += 1
            else:
                cls["run_failed"] += 1
        want = table["check_plan_counts"][check]["full" if full else "routine"]
        if len(planned) != want or sum(cls.values()) != want:
            problems.append(f"{check}: classes {cls} do not add up to the planned {want}")
        if cls["not_executed"] or cls["run_failed"]:
            problems.append(f"{check}: {cls['not_executed']} pairs not executed, {cls['run_failed']} in failed runs")
        summary[check] = {"planned": want, **cls}
    return summary, problems


def category_summary(records, table, full=FULL):
    out = {}
    for run in selected_runs(table, full):
        rec = records.get(run["id"])
        c = out.setdefault(run["category"], {"planned": 0, "passed": 0, "failed": 0, "not_executed": 0,
                                             "ei_success": 0, "ei_not_converged": 0})
        c["planned"] += 1
        if rec is None:
            c["not_executed"] += 1
        elif rec["failures"]:
            c["failed"] += 1
        else:
            c["passed"] += 1
            if run["category"] == "EI":
                c["ei_success" if rec["outcome"] == "success" else "ei_not_converged"] += 1
    return out


# ----------------------------------------------------------------------------- tests

class ExpectedStatesTable(unittest.TestCase):
    """Table identity and planned counts (no implementation needed)."""

    def test_sha_and_counts(self):
        table = load_table()
        self.assertEqual(len(table["runs"]), table["counts"]["runs"])
        self.assertEqual(sum(r["routine"] for r in table["runs"]), table["counts"]["routine_runs"])
        self.assertEqual(sum(1 for r in table["runs"] if r["category"] == "MS"), 1951)
        self.assertEqual(sum(1 for r in table["runs"] if r["category"] == "MS" and r["routine"]), 1255)
        self.assertEqual(check_plan_counts(table), [])

    def test_cli_configs_read_back(self):
        table = load_table()
        with tempfile.TemporaryDirectory() as folder:
            for run in table["runs"]:
                if run["entry"] != "CLI":
                    continue
                with self.subTest(run=run["id"]):
                    _, problems = cli_config(run, table, folder)
                    self.assertEqual(problems, [])

    def test_tie_rule_of_max_budget_ratio(self):
        values = {name: np.array([1.0, 1.0]) for name in OUTPUTS}
        errors = {name: np.zeros(2) for name in OUTPUTS}
        settings = {"epsabs": 1e-12, "epsrel": 1e-9}
        self.assertEqual(recompute_budget_fields(values, errors, settings)[1],
                         {"value": 0.0, "column": "term_13a", "index": 0})
        errors["hole"] = np.array([0.0, 1e-9])
        errors["total"] = np.array([1e-9, 1e-9])
        self.assertEqual(recompute_budget_fields(values, errors, settings)[1],
                         {"value": 1.0, "column": "hole", "index": 1})


class ExpectedStatesRuns(unittest.TestCase):
    """V13, MS and the R2 part of V15 over the routine runs (all runs with XPSLOSS_M3A_FULL=1)."""

    @classmethod
    def setUpClass(cls):
        cls.table = load_table()
        cls.records = {}
        cls.workdir = tempfile.mkdtemp(prefix="m3a_states_")
        for run in selected_runs(cls.table):
            cls.records[run["id"]] = execute_and_judge(run, cls.table, cls.workdir)

    def test_every_run_as_registered(self):
        failures = [f for rec in self.records.values() for f in rec["failures"]]
        self.assertEqual(failures[:200], [], f"{len(failures)} failures")

    def test_cross_run_conditions(self):
        self.assertEqual(cross_run_checks(self.records, self.table), [])

    def test_plan_accounting(self):
        summary, problems = plan_accounting(self.records, self.table)
        self.assertEqual(problems, [], json.dumps(summary))

    def test_categories_complete(self):
        summary = category_summary(self.records, self.table)
        want = {c: sum(1 for r in selected_runs(self.table) if r["category"] == c) for c in ("MS", "NC", "EI", "RJ")}
        for c, n in want.items():
            with self.subTest(category=c):
                self.assertEqual(summary[c]["planned"], n)
                self.assertEqual(summary[c]["passed"], n, json.dumps(summary[c]))

    def test_v15_recomputed_for_every_ms_run(self):
        ms = [r for r in self.records.values() if r["category"] == "MS"]
        self.assertEqual(len(ms), 1951 if FULL else 1255)
        self.assertTrue(all(r["outcome"] == "success" for r in ms))


if __name__ == "__main__":
    unittest.main()
