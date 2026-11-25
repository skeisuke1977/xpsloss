"""M3a checks owned by the reference author R1 (docs/M3A_PREREG.md, sections 5 and 6):
V1, V3, V4, V5, V6, V7 (with the successful EI runs), V11(iii)(1) and V11(v).

The expected values are frozen in benchmarks/numerics/m3a/references_m3a.json, written by
benchmarks/numerics/m3a/generate_references_m3a.py (mpmath, no package functions) and
identified by the SHA-256 in the appendix of docs/M3A_PREREG.md.  The runs and the check
plan come from benchmarks/numerics/m3a/expected_states_m3a.json (SHA-256 of section 3).
This file calls the public API (xpsloss.effective_cross_section_depth_averaged, the CLI
for the CLI runs of the plan) and the depth kernel
xpsloss._depth_kernel.depth_averaged_weights only, and does not import mpmath.

Arithmetic: the stored decimals are read as exact rationals; the returned binary64 values,
the hole charge, the weight length and the precision settings are used as exact rationals;
group sums and charge multiples are formed in rationals.

V7 is three-valued for every (check, run id, original index, output column):
    pass          d + u <= min(B_ref-, e)
    fail          max(0, d - u) > min(B_ref+, e)
    undetermined  otherwise (the reference has to be refined; never counted as a pass)
with d = |K - K_ref|, u the reference uncertainty, e the published error estimate
(component_errors) and B_ref-+ = max(epsabs, epsrel*(|K_ref| -+ u)).

Two scopes: the routine part uses the plan entries of the runs marked "routine" in the
table; the full part (all entries) runs only with XPSLOSS_M3A_FULL=1.  The numbers of
planned entries are compared with check_plan_counts before the runs, and the numbers of
pass / fail / undetermined / EI not converged / not run units are compared with the plan
afterwards.  A missing field or a missing reference is a failure, not a skip.
"""
from collections import Counter
from collections.abc import Mapping
import csv
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import re
import struct
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TABLE_FILE = ROOT / "benchmarks" / "numerics" / "m3a" / "expected_states_m3a.json"
TABLE_SHA256 = "d87dad8f02a2fd211679229bd7c5205d908dbb2b3b241b30a27e94469f414ed5"
REFERENCE_FILE = ROOT / "benchmarks" / "numerics" / "m3a" / "references_m3a.json"
PREREG = ROOT / "docs" / "M3A_PREREG.md"
FULL = os.environ.get("XPSLOSS_M3A_FULL") == "1"

TERMS = ("term_13a", "term_13b", "term_14", "term_15_T1", "term_15_T2",
         "term_15_T3", "term_16_e", "term_16_h")
NO_HOLE_TERMS = ("term_13a", "term_13b", "term_15_T1", "term_15_T3", "term_16_e")
HOLE_TERMS = ("term_14", "term_15_T2", "term_16_h")
OUTPUTS = TERMS + ("no_hole", "hole", "total")
KERNEL_COMPONENTS = ("W13", "N", "W14", "W15b", "W15c", "W16e", "W16h")
UR = Fraction(1, 2 ** 53)
CLI_INPUTS = {"Al": "al_elfdb_v1.0.bulk.json", "Si": "si_elfdb_v1.0.bulk.json"}
CSV_COLUMNS = ("loss_eV", "total", "no_hole", "hole") + TERMS + ("quadrature_error",)

_CACHE = {}


class MissingField(AssertionError):
    pass


def need(mapping, name, where):
    """mapping[name]; a missing field is a failure."""
    if not isinstance(mapping, Mapping) or name not in mapping:
        raise MissingField(f"missing field {name!r} in {where}")
    return mapping[name]


def table():
    if "table" not in _CACHE:
        raw = TABLE_FILE.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        if digest != TABLE_SHA256:
            raise AssertionError(f"expected_states_m3a.json has SHA-256 {digest}, not the pre-registered one")
        doc = json.loads(raw.decode("utf-8"))
        doc["_runs"] = {r["id"]: r for r in doc["runs"]}
        _CACHE["table"] = doc
    return _CACHE["table"]


def frozen_reference_sha256():
    """SHA-256 of references_m3a.json recorded in the appendix of the pre-registration."""
    text = PREREG.read_text(encoding="utf-8")
    appendix = text[text.rindex("## 付録"):]
    found = re.findall(r"references_m3a\.json` の SHA-256: `([0-9a-f]{64})`", appendix)
    if not found:
        raise AssertionError("the appendix of docs/M3A_PREREG.md does not record the SHA-256 of references_m3a.json")
    return found[0]


def references():
    if "refs" not in _CACHE:
        if not REFERENCE_FILE.exists():
            raise AssertionError("benchmarks/numerics/m3a/references_m3a.json is missing")
        raw = REFERENCE_FILE.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        frozen = frozen_reference_sha256()
        if digest != frozen:
            raise AssertionError(f"references_m3a.json has SHA-256 {digest}, the pre-registration records {frozen}")
        _CACHE["refs"] = json.loads(raw.decode("ascii"))
    return _CACHE["refs"]


def loss_array(name):
    key = ("array", name)
    if key not in _CACHE:
        GA = 0.2 * np.arange(1, 500)
        F = 0.2 / 4 * np.arange(4, 1997)
        arrays = {"GA": GA, "GC": np.linspace(0.2, 99.8, 499), "F": F, "F2": F[::2],
                  "FRD": np.r_[F[::-4], F[:40:4]], "GA10G": np.r_[GA[:10], 1.12]}
        arr = np.ascontiguousarray(arrays[name], dtype="<f8")
        want = table()["arrays"][name]["sha256_float64_le"]
        if hashlib.sha256(arr.tobytes()).hexdigest() != want:
            raise AssertionError(f"loss array {name} does not reproduce the pre-registered bytes")
        _CACHE[key] = arr.astype(float)
    return _CACHE[key]


def run_losses(run):
    spec = run["losses"]
    if "array" in spec:
        return loss_array(spec["array"])
    return np.array([float.fromhex(h) for h in spec["values_hex"]], dtype=float)


def api():
    import xpsloss
    fn = getattr(xpsloss, "effective_cross_section_depth_averaged", None)
    if fn is None:
        raise AssertionError("xpsloss.effective_cross_section_depth_averaged is not implemented")
    return xpsloss, fn


class Outputs:
    """The eleven outputs, their published errors and the direct-sum diagnostics of one run."""

    def __init__(self, losses, values, errors, direct_sum, where):
        n = losses.size
        self.losses = losses
        self.values, self.errors = {}, {}
        for name in OUTPUTS:
            v = np.atleast_1d(np.asarray(need(values, name, where), dtype=float))
            e = np.atleast_1d(np.asarray(need(errors, name, where + " component_errors"), dtype=float))
            if v.shape != (n,) or e.shape != (n,):
                raise MissingField(f"{where}: output {name} does not have one value and one error per loss")
            self.values[name], self.errors[name] = v, e
        if not isinstance(direct_sum, (list, tuple)) or len(direct_sum) != n:
            raise MissingField(f"{where}: diagnostics['direct_sum'] is not a list of length {n}")
        self.direct_sum = list(direct_sum)


def execute_api(run, losses):
    xpsloss, fn = api()
    tab = table()
    spec = tab["materials"][run["material"]]
    material = xpsloss.Material(run["material"], spec["gap_eV"],
                                tuple(xpsloss.Oscillator(*o) for o in spec["oscillators"]))
    kwargs = dict(weight_length_A=float.fromhex(run["weight_length_A"]["hex"]),
                  epsilon_momentum=run["epsilon_momentum"], theta_deg=run["theta_deg"],
                  hole_charge=float.fromhex(run["hole_charge"]["hex"]))
    if run["settings"] != "APIDEFAULT":
        kwargs.update(tab["settings"][run["settings"]])
    result = fn(material, run["electron_energy_eV"], losses, **kwargs)
    where = f"result of {run['id']}"
    got_losses = np.atleast_1d(np.asarray(result.loss_eV, dtype=float))
    if got_losses.shape != losses.shape or got_losses.tobytes() != losses.tobytes():
        raise AssertionError(f"{where}: loss_eV differs from the input")
    values = {name: getattr(result, name, None) for name in OUTPUTS}
    values = {name: v for name, v in values.items() if v is not None}
    errors = getattr(result, "component_errors", None)
    diagnostics = getattr(result, "diagnostics", None)
    return Outputs(losses, values, errors, need(diagnostics, "direct_sum", where + " diagnostics"), where)


def execute_cli(run, losses):
    tab = table()
    if run["material"] not in CLI_INPUTS:
        raise AssertionError(f"{run['id']}: no CLI input file is known to this test for {run['material']}")
    config = json.loads((ROOT / "benchmarks" / "inputs" / CLI_INPUTS[run["material"]]).read_text(encoding="utf-8"))
    spec = tab["materials"][run["material"]]
    got = [[o["energy_eV"], o["strength_eV2"], o["width_eV"], o["dispersion"]] for o in config["material"]["oscillators"]]
    if config["material"]["gap_eV"] != spec["gap_eV"] or got != spec["oscillators"]:
        raise AssertionError(f"{run['id']}: the CLI input file does not hold the material of the table")
    config["loss_grid_eV"] = run["cli_config"]["loss_grid_eV"]
    config.pop("integration", None)
    if run["cli_config"]["integration"] is not None:
        config["integration"] = dict(tab["settings"][run["cli_config"]["integration"]])
    with tempfile.TemporaryDirectory() as tmp:
        cfg = Path(tmp) / "config.json"
        cfg.write_text(json.dumps(config), encoding="utf-8")
        name = tab["cli_strings"]["ja_output_dir"] if run.get("output_dir_name") == "ja" else "out"
        out = Path(tmp) / name
        cmd = [sys.executable, "-m", "xpsloss", "xps-depth-averaged", str(cfg),
               "--weight-length-A", repr(float.fromhex(run["weight_length_A"]["hex"])),
               "--epsilon-momentum", run["epsilon_momentum"], "--theta-deg", repr(float(run["theta_deg"])),
               "--hole-charge", repr(float.fromhex(run["hole_charge"]["hex"])), "--out", str(out)]
        if run.get("weight_length_source") is not None:
            cmd += ["--weight-length-source", run["weight_length_source"]]
        env = dict(os.environ)
        env["PYTHONPATH"] = str(ROOT / "src") + os.pathsep + env.get("PYTHONPATH", "")
        done = subprocess.run(cmd, cwd=ROOT, env=env, capture_output=True, check=False)
        if done.returncode != 0:
            raise RuntimeError(f"CLI exit code {done.returncode}: {done.stderr[:400]!r}")
        with open(out / "xps_depth_averaged.csv", newline="", encoding="ascii") as handle:
            rows = list(csv.DictReader(handle))
        manifest = json.loads((out / "run.json").read_text(encoding="utf-8"))
    where = f"CLI output of {run['id']}"
    for column in CSV_COLUMNS:
        if not rows or column not in rows[0]:
            raise MissingField(f"{where}: CSV column {column} is missing")
    got_losses = np.array([float(r["loss_eV"]) for r in rows])
    if got_losses.shape != losses.shape or got_losses.tobytes() != losses.tobytes():
        raise AssertionError(f"{where}: the losses of the CSV differ from the pre-registered grid")
    values = {name: np.array([float(r[name]) for r in rows]) for name in OUTPUTS}
    diagnostics = need(manifest, "diagnostics", where + " run.json")
    return Outputs(losses, values, need(diagnostics, "component_errors", where + " run.json diagnostics"),
                   need(diagnostics, "direct_sum", where + " run.json diagnostics"), where)


def execute(run_id):
    """("ok", Outputs) or ("error", exception); each run is executed once."""
    key = ("run", run_id)
    if key not in _CACHE:
        run = table()["_runs"][run_id]
        losses = run_losses(run)
        try:
            out = execute_cli(run, losses) if run["entry"] == "CLI" else execute_api(run, losses)
            _CACHE[key] = ("ok", out)
        except MissingField:
            raise
        except Exception as error:          # judged by the caller (EI runs may fail in the allowed ways)
            _CACHE[key] = ("error", error)
    return _CACHE[key]


def exact(x):
    """Fraction of a finite binary64 number, or None."""
    x = float(x)
    return Fraction(x) if math.isfinite(x) else None


def judge(got, err, ref, u, epsabs, epsrel):
    """Three-valued V7 comparison in exact rationals."""
    g, e = exact(got), exact(err)
    if g is None or e is None or e < 0:
        return "fail"
    d = abs(g - ref)
    b_minus = max(epsabs, epsrel * max(abs(ref) - u, 0))
    b_plus = max(epsabs, epsrel * (abs(ref) + u))
    if d + u <= min(b_minus, e):
        return "pass"
    if max(Fraction(0), d - u) > min(b_plus, e):
        return "fail"
    return "undetermined"


def ref_output(entry, name):
    out = need(need(entry, "outputs", "reference entry"), name, "reference outputs")
    return Fraction(need(out, "value", name)), Fraction(need(out, "u", name))


def ref_total(entry, charge):
    """K_e + Q_h K_h with u_e + |Q_h| u_h; the analytic zero below the gap at unit charge."""
    if need(entry, "region", "reference entry") == "below_gap" and charge == 1:
        if "total_unit_charge" not in entry["analytic_zero"]:
            raise MissingField("analytic zero of the unit-charge total is not recorded")
        return Fraction(0), Fraction(0)
    (ve, ue), (vh, uh) = ref_output(entry, "no_hole"), ref_output(entry, "hole")
    return ve + charge * vh, ue + abs(charge) * uh


def settings_of(run):
    st = table()["settings"][run["settings"]]
    return Fraction(st["epsabs"]), Fraction(st["epsrel"])


def is_subgap(run, loss):
    gap = table()["materials"][run["material"]]["gap_eV"]
    return gap > 0 and loss < gap


class Tally:
    """Statuses per (run id, original index, column)."""

    def __init__(self, check):
        self.check = check
        self.units = {}
        self.notes = []

    def add(self, run_id, index, column, status, note=None):
        key = (run_id, index, column)
        if key in self.units:
            raise AssertionError(f"{self.check}: unit {key} judged twice")
        self.units[key] = status
        if status not in ("pass", "ei_not_converged") and len(self.notes) < 12:
            self.notes.append(f"{status} {run_id}[{index}] {column}" + (f": {str(note)[:160]}" if note else ""))

    def counts(self):
        return dict(Counter(self.units.values()))

    def finish(self, case, planned_units):
        counts = self.counts()
        print(f"M3A {case.scope} {self.check}: planned {planned_units}, " + json.dumps(counts, sort_keys=True))
        case.assertEqual(len(self.units), planned_units, f"{self.check}: judged units differ from the plan")
        bad = {k: v for k, v in counts.items() if k not in ("pass", "ei_not_converged")}
        case.assertFalse(bad, f"{self.check}: {counts}; first cases: {self.notes}")


class FrozenFiles(unittest.TestCase):
    def test_table_and_reference_hashes(self):
        tab = table()
        refs = references()
        self.assertEqual(need(refs, "table", "references")["sha256"], TABLE_SHA256)
        self.assertEqual(set(need(refs, "references", "references")), {r["key"] for r in tab["reference_keys"]})
        self.assertEqual(set(need(refs, "limits", "references")), {r["key"] for r in tab["limit_keys"]})
        self.assertEqual(len(refs["references"]), tab["counts"]["reference_keys"])
        self.assertEqual(len(refs["limits"]), tab["counts"]["limit_keys"])
        self.assertNotIn("partial", refs)

    def test_reference_uncertainty_is_small(self):
        """u_j <= 0.01 B_ref for the tightest use of every key, except the listed undetermined keys."""
        refs = references()
        listed = {item["key"] for item in need(need(refs, "summary", "references"), "undetermined", "summary")}
        for key, entry in refs["references"].items():
            if key not in listed:
                self.assertLessEqual(need(entry, "max_u_over_Bref", key), 0.01, key)
        self.assertFalse(listed, f"undetermined reference keys are listed: {sorted(listed)}")

    def test_plan_counts(self):
        tab = table()
        for check, counts in tab["check_plan_counts"].items():
            entries = tab["check_plan"][check]
            self.assertEqual(len(entries), counts["full"], check)
            self.assertEqual(sum(1 for e in entries if tab["_runs"][e[0]]["routine"]), counts["routine"], check)


class DepthKernelV1(unittest.TestCase):
    """V1: scaled weights p^2 W/(2b) against the exact rationals of R1."""

    @classmethod
    def setUpClass(cls):
        cls.v1 = need(references(), "v1", "references")

    def kernel(self):
        try:
            from xpsloss._depth_kernel import depth_averaged_weights
        except ImportError as error:
            self.fail(f"xpsloss._depth_kernel.depth_averaged_weights is not implemented: {error}")
        return depth_averaged_weights

    @staticmethod
    def inputs(point):
        return tuple(float.fromhex(point[name]) for name in ("k_hex", "b_hex", "p_hex"))

    @staticmethod
    def rational(item):
        return Fraction(int(item["numerator"]), int(item["denominator"]))

    def parts(self, weights, name):
        return (self.rational(weights[name]["re"]), self.rational(weights[name]["im"]))

    @staticmethod
    def widths(kf, bf, pf):
        k, b, p = Fraction(kf), Fraction(bf), Fraction(pf)
        L = p + k
        D, H = k * k + b * b, L * L + b * b
        scale = p * p / (2 * b)
        tol = Fraction(1, 10 ** 13)
        return {("W13", 0): tol * (k * L + b * b) / (p * D * H) * scale,
                ("W15c", 0): tol * (k * L + b * b) / (p * D * H) * scale,
                ("W16h", 1): tol * (k * L + b * b) / (D * H) * scale}

    def returned(self, out, name, where):
        value = need(out, name, where)
        self.assertEqual(len(value), 4, f"{where} {name}: (re, im, err_re, err_im) expected")
        re_, im_, ere, eim = (float(v) for v in value)
        for v in (re_, im_, ere, eim):
            self.assertTrue(math.isfinite(v), f"{where} {name}: non-finite value")
        self.assertGreaterEqual(ere, 0.0)
        self.assertGreaterEqual(eim, 0.0)
        return (Fraction(re_), Fraction(im_)), (Fraction(ere), Fraction(eim)), (re_, im_)

    def test_points_and_stored_decimals(self):
        points = need(self.v1, "points", "v1")
        self.assertEqual(len(points), 35)
        labels = [p["label"] for p in points]
        for must in ("k<<p", "k~p", "k>>p", "b>>p", "zero_crossing_before", "zero_crossing",
                     "zero_crossing_after", "underflow"):
            self.assertIn(must, labels)
        self.assertEqual(self.inputs(points[labels.index("zero_crossing_after")]),
                         (0.5, 1.0, math.nextafter(1.5, math.inf)))
        self.assertEqual(self.inputs(points[labels.index("underflow")]), (1.0, 1.0, 1e82))
        for point in points:
            for name in KERNEL_COMPONENTS:
                for part in ("re", "im"):
                    item = point["weights"][name][part]
                    value = self.rational(item)
                    self.assertLessEqual(abs(Fraction(item["decimal"]) - value), abs(value) * Fraction(1, 10 ** 39))
        for point in points[:8]:
            quad = need(point, "depth_quadrature", point["label"])
            for name in KERNEL_COMPONENTS:
                ref = self.parts(point["weights"], name)
                size = max(abs(ref[0]), abs(ref[1]))
                diff = Fraction(quad["components"][name]["abs_difference_from_exact"])
                self.assertLessEqual(diff, size * Fraction(1, 10 ** 20), f"{point['label']} {name}")

    def test_values_and_errors(self):
        fn = self.kernel()
        for point in need(self.v1, "points", "v1"):
            kf, bf, pf = self.inputs(point)
            where = f"V1 point {point['label']}"
            out = fn(kf, bf, pf)
            widths = self.widths(kf, bf, pf)
            for name in KERNEL_COMPONENTS:
                got, err, raw = self.returned(out, name, where)
                ref = self.parts(point["weights"], name)
                for part in (0, 1):
                    tag = f"{where} {name}.{'re' if part == 0 else 'im'}"
                    diff = abs(got[part] - ref[part])
                    if name == "N" and part == 1:
                        self.assertEqual(raw[1], 0.0, tag)
                        continue
                    width = widths.get((name, part), Fraction(1, 10 ** 13) * abs(ref[part]))
                    self.assertLessEqual(diff, width, tag + " value")
                    self.assertLessEqual(diff, err[part], tag + " published error")

    def test_underflow_and_out_of_range(self):
        fn = self.kernel()
        out = fn(1.0, 1.0, 1e82)
        got = float(need(out, "W14", "underflow point")[0])
        self.assertNotEqual(got, 0.0)
        target = Fraction("1.0000000000000002e-164")
        self.assertLessEqual(abs(Fraction(got) - target), target * Fraction(1, 10 ** 13))
        control = need(self.v1, "out_of_range_control", "v1")
        with self.assertRaises(FloatingPointError):
            fn(*(float.fromhex(control[name]) for name in ("k_hex", "b_hex", "p_hex")))

    def test_zero_crossing_error_floor(self):
        """Just after the zero crossing the published error of the three crossing parts
        is at least the exact value (only the e_n part is left there)."""
        fn = self.kernel()
        after = need(self.v1, "zero_crossing_after", "v1")
        self.assertTrue(after["e_n_covers_exact_n"])
        self.assertTrue(after["terms_proportional_to_n_hat_vanish"])
        self.assertTrue(after["e_n_part_covers_exact_value"])
        out = fn(0.5, 1.0, math.nextafter(1.5, math.inf))
        for name, part in (("W13", 0), ("W15c", 0), ("W16h", 1)):
            _, err, _ = self.returned(out, name, "zero crossing")
            size = self.rational(after["exact_magnitude"][f"{name}.{'re' if part == 0 else 'im'}"])
            self.assertGreater(size, 0)
            self.assertGreaterEqual(err[part], size, f"{name} part {part}")

    def test_input_error_box(self):
        fn = self.kernel()
        boxed = [p for p in need(self.v1, "points", "v1") if "box_corners" in p]
        self.assertEqual([p["label"] for p in boxed], ["k<<p", "k~p", "k>>p", "b>>p", "zero_crossing_after"])
        for point in boxed:
            kf, bf, pf = self.inputs(point)
            ulps = point["input_error_ulps"]
            self.assertEqual(ulps, 4)
            out = fn(kf, bf, pf, e_k=ulps * math.ulp(kf), e_b=ulps * math.ulp(bf), e_p=ulps * math.ulp(pf))
            self.assertEqual(len(point["box_corners"]), 8)
            for corner in point["box_corners"]:
                for name in KERNEL_COMPONENTS:
                    got, err, _ = self.returned(out, name, point["label"])
                    ref = self.parts(corner["weights"], name)
                    for part in (0, 1):
                        self.assertLessEqual(abs(got[part] - ref[part]), err[part],
                                             f"{point['label']} corner {corner['signs']} {name} part {part}")

    def test_recorded_mutations(self):
        mutations = need(self.v1, "mutations", "v1")
        self.assertEqual(set(mutations), {"K1", "K2", "K3", "K4"})
        for name, record in mutations.items():
            self.assertGreaterEqual(record["points_with_a_component_beyond_the_width"], 1, name)
            self.assertTrue(record["registered_detection_of_prereg_table_holds"], name)


class _Checks:
    """Checks against the multiprecision references for one scope ("routine" or "full")."""

    scope = None

    @classmethod
    def setUpClass(cls):
        cls.table = table()
        cls.refs = references()
        cls.runs = cls.table["_runs"]

    # ------------------------------------------------------------------ helpers
    def plan(self, check):
        entries = self.table["check_plan"][check]
        if self.scope == "routine":
            entries = [e for e in entries if self.runs[e[0]]["routine"]]
        self.assertEqual(len(entries), self.table["check_plan_counts"][check][self.scope],
                         f"{check}: planned entries differ from check_plan_counts")
        return entries

    def reference(self, key):
        return need(self.refs["references"], key, "references (reference not made)")

    def outputs(self, run_id, tally, index, columns):
        """Outputs of a successful run, or None after recording the units of a failed run."""
        api()
        state, value = execute(run_id)
        if state == "ok":
            return value
        run = self.runs[run_id]
        status, note = "not_run", f"{type(value).__name__}: {value}"
        if run["category"] == "EI":
            allowed = run["expect"]["allowed_reasons"]
            ok = (type(value).__name__ == run["expect"]["exception"] and getattr(value, "reason", None) in allowed
                  and getattr(value, "index", None) == run["expect"]["index"])
            status = "ei_not_converged" if ok else "fail"
        for column in columns:
            tally.add(run_id, index, column, status, note)
        return None

    def run_budget_ok(self, run, out, index):
        """Finite, non-negative estimates within B_run (success criterion of MS, used for EI)."""
        epsabs, epsrel = settings_of(run)
        for name in OUTPUTS:
            v, e = exact(out.values[name][index]), exact(out.errors[name][index])
            if v is None or e is None or e < 0 or e > max(epsabs, epsrel * abs(v)):
                return False
        return True

    def charge(self, run):
        return Fraction(float.fromhex(run["hole_charge"]["hex"]))

    def length(self, run):
        return Fraction(float.fromhex(run["weight_length_A"]["hex"]))

    # ------------------------------------------------------------------ V7
    def test_v7_reference(self):
        tally = Tally("V7")
        entries = self.plan("V7")
        for run_id, index, key in entries:
            run = self.runs[run_id]
            out = self.outputs(run_id, tally, index, OUTPUTS)
            if out is None:
                continue
            entry = self.reference(key)
            epsabs, epsrel = settings_of(run)
            ei_ok = run["category"] != "EI" or self.run_budget_ok(run, out, index)
            for name in OUTPUTS:
                ref, u = ref_total(entry, self.charge(run)) if name == "total" else ref_output(entry, name)
                status = judge(out.values[name][index], out.errors[name][index], ref, u, epsabs, epsrel)
                if status == "pass" and not ei_ok:
                    status = "fail"
                tally.add(run_id, index, name, status, f"got {out.values[name][index]!r} ref {float(ref)!r}")
        tally.finish(self, len(entries) * len(OUTPUTS))

    # ------------------------------------------------------------------ V6
    def test_v6_zero_response(self):
        api()
        runs = [r for r in self.table["runs"]
                if (r["family"] == "MS7" or (r["family"] == "MS12" and r["material"] == "Z0"))
                and (self.scope == "full" or r["routine"])]
        self.assertTrue(runs)
        for run in runs:
            state, out = execute(run["id"])
            if state != "ok":
                self.fail(f"{run['id']}: required success raised {type(out).__name__}: {out}")
            for name in OUTPUTS:
                self.assertEqual(out.values[name].astype("<f8").tobytes(), bytes(8 * out.losses.size),
                                 f"{run['id']} {name}: not bitwise zero")

    def test_v6_linearity(self):
        tally = Tally("V6")
        entries = [e for e in self.plan("V7")
                   if self.runs[e[0]]["family"] == "MS16"
                   or (self.runs[e[0]]["family"] in ("MS5", "MS11") and self.runs[e[0]]["weight_length_A"]["value"] == 20.0)]
        self.assertTrue(entries)
        planned = 0
        for run_id, index, key in entries:
            run = self.runs[run_id]
            charge = self.charge(run)
            columns = ("total",) + ((HOLE_TERMS + ("hole",)) if charge == 0 else ())
            planned += len(columns)
            out = self.outputs(run_id, tally, index, columns)
            if out is None:
                continue
            entry = self.reference(key)
            epsabs, epsrel = settings_of(run)
            for name in columns:
                ref, u = ref_total(entry, charge) if name == "total" else ref_output(entry, name)
                tally.add(run_id, index, name,
                          judge(out.values[name][index], out.errors[name][index], ref, u, epsabs, epsrel))
        tally.finish(self, planned)

    # ------------------------------------------------------------------ V3
    def test_v3_abc_term_13a(self):
        tally = Tally("V3abc")
        seen, planned = set(), 0
        kinds = {"Al0": "analytic_no_dispersion", "Y0": "analytic_no_dispersion", "S0": "analytic_no_dispersion",
                 "Al": "bulk_integral_total_momentum", "Y5": "bulk_integral_total_momentum",
                 "Si": "bulk_integral_total_momentum"}
        for check in ("V7", "V11iii", "V11v"):
            for run_id, index, key in self.plan(check):
                if (run_id, index) in seen:
                    continue
                seen.add((run_id, index))
                run = self.runs[run_id]
                loss = float(run_losses(run)[index])
                below = is_subgap(run, loss)
                independent = None
                if not below and check == "V7" and run["material"] in kinds:
                    dispersive = kinds[run["material"]] == "bulk_integral_total_momentum"
                    if not dispersive or run["epsilon_momentum"] == "k_total":
                        group = need(self.refs["groups"], self.reference(key)["group"], "groups")
                        independent = need(group, "term_13a_independent", f"group of {key}")
                        self.assertEqual(independent["kind"], kinds[run["material"]])
                if not below and independent is None:
                    continue
                planned += 1
                out = self.outputs(run_id, tally, index, ("term_13a",))
                if out is None:
                    continue
                got = out.values["term_13a"][index]
                if below:                               # (c): bitwise zero
                    status = "pass" if struct.pack("<d", float(got)) == bytes(8) else "fail"
                else:                                   # (a), (b): V7 condition
                    epsabs, epsrel = settings_of(run)
                    status = judge(got, out.errors["term_13a"][index], Fraction(independent["value"]),
                                   Fraction(independent["u"]), epsabs, epsrel)
                tally.add(run_id, index, "term_13a", status, f"got {got!r}")
        self.assertGreater(planned, 0)
        tally.finish(self, planned)

    def test_v3_d_term_13a_at_extreme_lengths(self):
        tally = Tally("V3d")
        entries = self.plan("V3d")
        for run_id, index, key in entries:
            run = self.runs[run_id]
            out = self.outputs(run_id, tally, index, ("term_13a",))
            if out is None:
                continue
            entry = self.reference(key)
            self.assertEqual(float.fromhex(entry["weight_length_A_hex"]), 20.0)
            ref, u = ref_output(entry, "term_13a")
            epsabs, epsrel = settings_of(run)
            tally.add(run_id, index, "term_13a",
                      judge(out.values["term_13a"][index], out.errors["term_13a"][index], ref, u, epsabs, epsrel))
        tally.finish(self, len(entries))

    # ------------------------------------------------------------------ V4
    def limit(self, key, kind, run):
        entry = need(self.refs["limits"], key, "limits (reference not made)")
        self.assertEqual(entry["kind"], kind)
        self.assertEqual(float.fromhex(entry["hole_charge_hex"]), float.fromhex(run["hole_charge"]["hex"]), key)
        self.assertIn(run["weight_length_A"]["hex"], [float.fromhex(h).hex() for h in entry["weight_lengths_A_hex"]])
        return entry

    @staticmethod
    def quantity(entry, field, name):
        item = need(need(entry, field, "limit entry"), name, field)
        return Fraction(item["value"]), Fraction(item["u"])

    def test_v4_short_length_limit(self):
        tally = Tally("V4")
        entries = self.plan("V4")
        for run_id, index, key in entries:
            run = self.runs[run_id]
            out = self.outputs(run_id, tally, index, OUTPUTS)
            if out is None:
                continue
            entry = self.limit(key, "lambda_to_0", run)
            lam = self.length(run)
            for name in OUTPUTS:
                k_hat, e = exact(out.values[name][index]), exact(out.errors[name][index])
                if k_hat is None or e is None or e < 0:
                    tally.add(run_id, index, name, "fail", "non-finite")
                    continue
                g0, u_lim = self.quantity(entry, "G0", name)
                env, u_env = self.quantity(entry, "envelope", name)
                kappa, u_kappa = self.quantity(entry, "kappa", name)
                bound = lam * (abs(kappa) + env) + lam * e + u_lim + lam * u_env + lam * u_kappa
                ok = abs(lam * k_hat - g0) <= bound
                tally.add(run_id, index, name, "pass" if ok else "fail",
                          f"|H-G0| = {float(abs(lam * k_hat - g0))!r}, bound {float(bound)!r}")
        tally.finish(self, len(entries) * len(OUTPUTS))

    # ------------------------------------------------------------------ V5
    @staticmethod
    def weights_of(name, charge):
        if name in TERMS:
            return {name: Fraction(1)}
        w = {}
        if name in ("no_hole", "total"):
            w.update({m: Fraction(1) for m in NO_HOLE_TERMS})
        if name == "hole":
            w.update({m: Fraction(1) for m in HOLE_TERMS})
        if name == "total":
            w.update({m: charge for m in HOLE_TERMS})
        return w

    def test_v5_long_length_limit(self):
        tally = Tally("V5")
        entries = self.plan("V5")
        deviations = {}
        for run_id, index, key in entries:
            run = self.runs[run_id]
            out = self.outputs(run_id, tally, index, OUTPUTS)
            if out is None:
                continue
            entry = self.limit(key, "lambda_to_inf", run)
            lam = self.length(run)
            charge = self.charge(run)
            epsabs, _ = settings_of(run)
            per_length = need(entry["lengths"], float.fromhex(run["weight_length_A"]["hex"]).hex(), "limit lengths")
            k_hat = {m: exact(out.values[m][index]) for m in OUTPUTS}
            e_hat = {m: exact(out.errors[m][index]) for m in OUTPUTS}
            if any(v is None for v in k_hat.values()) or any(v is None or v < 0 for v in e_hat.values()):
                for name in OUTPUTS:
                    tally.add(run_id, index, name, "fail", "non-finite")
                continue
            below_unit = is_subgap(run, float(out.losses[index])) and charge == 1
            for name in OUTPUTS:
                w = self.weights_of(name, charge)
                tr = {m: c for m, c in w.items() if m != "term_13a" and c != 0}
                h_tr = lam * sum((c * k_hat[m] for m, c in tr.items()), Fraction(0))
                g_inf, u_inf = self.quantity(entry, "G_inf", name)
                e_tr, u_tr = self.quantity(per_length, "E_tr", name)
                numeric = lam * sum((abs(c) * e_hat[m] for m, c in tr.items()), Fraction(0))
                rho_h = 4 * UR * len(tr) * lam * sum((abs(c * k_hat[m]) for m, c in tr.items()), Fraction(0))
                ok = abs(h_tr - g_inf) <= e_tr + u_tr + numeric + rho_h + u_inf
                note = f"|Htr-Ginf| = {float(abs(h_tr - g_inf))!r}"
                if name == "total" and below_unit:
                    ok = ok and abs(k_hat["total"]) <= min(e_hat["total"], epsabs)
                else:
                    full = {m: c for m, c in w.items() if c != 0}
                    rho_k = 4 * UR * len(full) * sum((abs(c * k_hat[m]) for m, c in full.items()), Fraction(0))
                    assembled = sum((c * k_hat[m] for m, c in full.items()), Fraction(0))
                    if abs(k_hat[name] - assembled) > rho_k:
                        ok, note = False, "returned column differs from the sum of its terms"
                tally.add(run_id, index, name, "pass" if ok else "fail", note)
                deviations[(key, index, name, run["weight_length_A"]["hex"])] = (
                    h_tr - g_inf, numeric + rho_h, lam, entry, per_length, bool(tr))
        tally.finish(self, len(entries) * len(OUTPUTS))

        # ratio of the deviations at the two lengths (imposed only where the first order dominates)
        imposed = skipped = 0
        failures = []
        pairs = {}
        for (key, index, name, lam_hex), value in deviations.items():
            pairs.setdefault((key, index, name), []).append(value)
        for (key, index, name), items in sorted(pairs.items()):
            if len(items) != 2 or not items[0][5]:
                skipped += 1
                continue
            short, long_ = sorted(items, key=lambda item: item[2])
            entry = short[3]
            a1, u_a1 = self.quantity(entry, "A1", name)
            _, u_inf = self.quantity(entry, "G_inf", name)
            eligible = abs(a1) > u_a1
            for dev, numeric, lam, _, per_length, _ in (short, long_):
                rem, u_rem = self.quantity(per_length, "remainder_bound_after_first_order", name)
                if rem + u_rem + numeric + u_inf + u_a1 / lam > abs(a1) / (10 * lam):
                    eligible = False
            if not eligible:
                skipped += 1
                continue
            imposed += 1
            d_short, d_long = short[0], long_[0]
            lo, hi = Fraction(90, 11), Fraction(110, 9)
            ok = (d_long > 0 and lo * d_long <= d_short <= hi * d_long) or \
                 (d_long < 0 and hi * d_long <= d_short <= lo * d_long)
            if not ok:
                failures.append(f"{key}[{index}] {name}: {float(d_short)!r} / {float(d_long)!r}")
        print(f"M3A {self.scope} V5 ratio: imposed {imposed}, not imposed {skipped}, failed {len(failures)}")
        self.assertFalse(failures, failures[:8])

    # ------------------------------------------------------------------ V11
    def test_v11_iii_1_subgap_terms_and_groups(self):
        tally = Tally("V11iii(1)")
        entries = self.plan("V11iii")
        columns = TERMS + ("no_hole", "hole")
        with_key = [e for e in entries if e[2] is not None]
        print(f"M3A {self.scope} V11iii: {len(entries)} planned entries, {len(with_key)} with a reference key")
        for run_id, index, key in with_key:
            run = self.runs[run_id]
            out = self.outputs(run_id, tally, index, columns)
            if out is None:
                continue
            entry = self.reference(key)
            self.assertEqual(entry["region"], "below_gap", key)
            self.assertEqual(float.fromhex(entry["loss_eV_hex"]).hex(), float(out.losses[index]).hex(), key)
            epsabs, epsrel = settings_of(run)
            ei_ok = run["category"] != "EI" or self.run_budget_ok(run, out, index)
            for name in columns:
                ref, u = ref_output(entry, name)
                status = judge(out.values[name][index], out.errors[name][index], ref, u, epsabs, epsrel)
                if status == "pass" and not ei_ok:
                    status = "fail"
                tally.add(run_id, index, name, status)
        tally.finish(self, len(with_key) * len(columns))

    def test_v11_v_subgap_total_other_charges(self):
        tally = Tally("V11v")
        entries = self.plan("V11v")
        for run_id, index, key in entries:
            run = self.runs[run_id]
            out = self.outputs(run_id, tally, index, ("total",))
            if out is None:
                continue
            charge = self.charge(run)
            self.assertNotEqual(charge, 1)
            if key == "analytic_zero":
                ref, u = Fraction(0), Fraction(0)
            else:
                entry = self.reference(key)
                self.assertEqual(entry["region"], "below_gap", key)
                ke, ue = ref_output(entry, "electron_group")
                ref, u = (1 - charge) * ke, abs(1 - charge) * ue
            epsabs, epsrel = settings_of(run)
            status = judge(out.values["total"][index], out.errors["total"][index], ref, u, epsabs, epsrel)
            if out.direct_sum[index] is not None:
                status = "fail"
            tally.add(run_id, index, "total", status, f"got {out.values['total'][index]!r} ref {float(ref)!r}")
        tally.finish(self, len(entries))


class ReferenceRoutine(_Checks, unittest.TestCase):
    scope = "routine"


@unittest.skipUnless(FULL, "full record: set XPSLOSS_M3A_FULL=1")
class ReferenceFull(_Checks, unittest.TestCase):
    scope = "full"


if __name__ == "__main__":
    unittest.main()
