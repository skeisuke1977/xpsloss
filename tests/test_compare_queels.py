"""Exercise the section 4.2 verdicts with synthetic depth tables only."""

from contextlib import redirect_stderr
from io import StringIO
import json
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

import numpy as np
from scipy.integrate import trapezoid

from benchmarks import compare_queels_fixed_depth as comparison
from xpsloss import effective_cross_section_fixed_depth
from xpsloss.config import load_config


ROOT = Path(__file__).resolve().parents[1]
COMPONENTS = {"total": "T", "no_hole": "e", "hole": "h"}
DEPTHS = np.arange(1, 66) * 1.25 + 0.1
COLUMN = 7  # 10.1 Å, the nearest depth to 10 Å in run 7


def write_run(folder, n, depths, curves):
    for name, letter in COMPONENTS.items():
        table = np.zeros((503, 65))
        table[0] = table[1] = depths  # normal emission: path equals depth
        table[3:-1, COLUMN] = curves[name]  # the final loss row is zero
        path = folder / f"K_XPS_{letter}_Al_1000_0_q1_n{n}.m"
        np.savetxt(path, table, fmt="%.7g", comments="",
                   header="synthetic\n@\n1000\n0\n20\n0.2\n500\n65")


class ComparisonTests(unittest.TestCase):
    def test_rejects_invalid_comparison_inputs(self):
        with TemporaryDirectory() as tmp:
            folder = Path(tmp)
            levels = {"total": 1., "no_hole": 0.75, "hole": 0.25}
            for points in (1000, 2000):
                write_run(folder, points, DEPTHS, levels)
            for case, points, component in (("infinite", 1000, "e"), ("path", 2000, "T"),
                                            ("depth", 1000, "h"), ("depth", 2000, "e"),
                                            ("zero_scale", 2000, "T"), ("negative_area", 2000, "T"),
                                            ("overflow", 1000, "e"), ("nonfinite_header", 2000, "e")):
                with self.subTest(case=case, points=points):
                    path = folder / f"K_XPS_{component}_Al_1000_0_q1_n{points}.m"
                    original = path.read_bytes()
                    header = original.decode("ascii").splitlines()[:8]
                    table = np.loadtxt(path, skiprows=8)
                    if case == "infinite":
                        table[17, COLUMN] = np.inf
                    elif case == "path":
                        table[1, COLUMN] += 0.1
                    elif case == "depth":
                        table[:2, COLUMN] += 0.1
                    elif case == "zero_scale":
                        table[3:, COLUMN] = 0
                    elif case == "negative_area":
                        table[20:-1, COLUMN] = -1
                    elif case == "overflow":
                        table[17, COLUMN] = 1e308
                    else:
                        header[4] = "nan"
                    np.savetxt(path, table, fmt="%.8g", comments="", header="\n".join(header))
                    out = folder / f"invalid_{case}_{points}.json"
                    errors = StringIO()
                    with redirect_stderr(errors):
                        code = comparison.main(["--reference-dir", str(folder), "--out", str(out)])
                    self.assertEqual(code, 2)
                    report = json.loads(out.read_text(encoding="utf-8"))
                    self.assertEqual(report["status"], "error")
                    self.assertEqual(report, json.loads(errors.getvalue()))
                    path.write_bytes(original)

    def test_synthetic_pass_and_perturbed_failure(self):
        losses = 0.2 * np.arange(1, 500)
        material = load_config(ROOT / "benchmarks/inputs/al_elfdb_v1.0.bulk.json").material
        model = effective_cross_section_fixed_depth(material, 1000, losses, DEPTHS[COLUMN],
                                                    epsilon_momentum="k_total")
        # The passing reference preserves total = no_hole + hole. The electron
        # bias needs the run-to-run uncertainty, while the hole uses the global floor.
        reference = {"no_hole": model.no_hole + 0.004 * model.total,
                     "hole": model.hole + 0.0008 * model.total}
        reference["total"] = reference["no_hole"] + reference["hole"]
        window = reference["total"][4:]  # 1.0 to 99.8 eV
        peak = window.max()
        mean = trapezoid(window, losses[4:]) / (losses[-1] - losses[4])

        # Change one reference component at a time to exercise the all-component
        # verdict independently of the diagnostic for total = no_hole + hole.
        for case in ("within", "hole_spike", "electron_area", "total_offset"):
            with self.subTest(case=case), TemporaryDirectory() as tmp:
                folder = Path(tmp)
                run7 = {name: curve.copy() for name, curve in reference.items()}
                if case == "hole_spike":
                    run7["hole"][99] = model.hole[99] + 1.25e-3 * peak
                if case == "electron_area":
                    run7["no_hole"] = model.no_hole + 1.25e-3 * mean
                if case == "total_offset":
                    run7["total"] = 1.009 * model.total

                run6 = dict(run7)
                run6["total"] = run7["total"] + 0.002 * model.total
                if case != "electron_area":
                    run6["no_hole"] = run7["no_hole"] + 0.002 * model.total
                write_run(folder, 2000, DEPTHS, run7)
                write_run(folder, 1000, DEPTHS + 0.0001, run6)

                out = folder / "verdict.json"
                run = subprocess.run(
                    [sys.executable, str(ROOT / "benchmarks/compare_queels_fixed_depth.py"),
                     "--reference-dir", str(folder), "--out", str(out)],
                    cwd=folder, env=dict(os.environ, PYTHONPATH=str(ROOT / "src")),
                    capture_output=True, text=True,
                )
                expected = {"k_parallel": False, "k_total": case == "within"}
                self.assertEqual(run.returncode, 0 if case == "within" else 1, run.stderr)
                self.assertEqual(json.loads(run.stdout), expected)
                report = json.loads(out.read_text(encoding="utf-8"))
                self.assertEqual(report["depth_A"], DEPTHS[COLUMN])
                self.assertEqual(report["depth_column"], COLUMN)
                self.assertAlmostEqual(report["depth_shift_run6_vs_run7_A"], 0.0001)
                self.assertEqual(report["loss_window_eV"], [1.0, 99.8])
                self.assertEqual(report["points"], 495)
                choices = report["choices"]
                self.assertEqual({key: value["pass"] for key, value in choices.items()}, expected)
                components = choices["k_total"]["components"]
                self.assertEqual(set(components), set(COMPONENTS))
                for name, metrics in components.items():
                    fails_point = (case, name) in {("hole_spike", "hole"), ("total_offset", "total")}
                    fails_area = (case, name) in {("electron_area", "no_hole"), ("total_offset", "total")}
                    self.assertEqual(metrics["pointwise_pass"], not fails_point, name)
                    self.assertEqual(metrics["area_pass"], not fails_area, name)
                if case == "within":
                    # All biases follow the total curve, so the peak and area
                    # ratios have the same analytical value. Allow 7-digit file rounding.
                    for name, ratio in (("no_hole", 0.004 / (3 * 0.002 + 0.001 * 1.0048)),
                                        ("hole", 0.0008 / (0.001 * 1.0048))):
                        for field in ("max_diff_over_tolerance", "area_diff_over_area_tolerance"):
                            self.assertAlmostEqual(components[name][field], ratio, delta=2e-4)
                if case == "hole_spike":
                    self.assertEqual(components["hole"]["points_failing"], 1)
