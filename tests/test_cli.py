import copy
import csv
import math
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from xpsloss.config import load_config, parse_config


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "synthetic_bulk.json"


class CLITests(unittest.TestCase):
    def call(self, *args):
        env = dict(os.environ)
        env["PYTHONPATH"] = str(ROOT / "src")
        return subprocess.run([sys.executable, "-m", "xpsloss", *map(str, args)],
                              cwd=ROOT, env=env, capture_output=True, text=True, check=False)

    def test_capabilities_and_validation(self):
        r = self.call("capabilities", "--json")
        self.assertEqual(r.returncode, 0, r.stderr)
        cap = json.loads(r.stdout)
        self.assertEqual(cap["package"], "xpsloss")
        self.assertFalse(cap["queels_xps_equivalence"])
        self.assertIn("core_hole", cap["bulk_unsupported_physics"])
        r = self.call("validate", EXAMPLE)
        self.assertEqual(r.returncode, 0, r.stderr)
        result = json.loads(r.stdout)
        self.assertEqual(result["input_sha256"], hashlib.sha256(EXAMPLE.read_bytes()).hexdigest())

    def test_bulk_end_to_end_and_overwrite_protection(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "bulk"
            input_arg = "./examples/synthetic_bulk.json"
            r = self.call("bulk", input_arg, "--out", output)
            self.assertEqual(r.returncode, 0, r.stderr)
            result = json.loads(r.stdout)
            self.assertEqual(result["point_count"], 160)
            manifest = json.loads((output / "run.json").read_text())
            self.assertEqual(manifest["status"], "complete")
            self.assertNotIn("applicability", manifest)
            self.assertEqual(manifest["input_file"], input_arg)
            self.assertEqual(manifest["input_sha256"], hashlib.sha256(EXAMPLE.read_bytes()).hexdigest())
            self.assertEqual(manifest["units"]["diimfp_per_A_eV"], "angstrom^-1 eV^-1")
            data = np.loadtxt(output / "bulk.csv", delimiter=",", skiprows=1)
            np.testing.assert_allclose(data[:, 0], np.arange(0.5, 80.01, 0.5), rtol=0, atol=1e-12)
            self.assertTrue(np.all(np.isfinite(data)))
            self.assertTrue(np.all(data[:, 1:] >= 0))
            before = (output / "bulk.csv").read_bytes()
            r = self.call("bulk", EXAMPLE, "--out", output)
            self.assertEqual(r.returncode, 2)
            self.assertEqual(json.loads(r.stderr)["error_type"], "FileExistsError")
            self.assertEqual(before, (output / "bulk.csv").read_bytes())

    def test_elf_end_to_end(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "elf"
            r = self.call("elf", EXAMPLE, "--out", output)
            self.assertEqual(r.returncode, 0, r.stderr)
            manifest = json.loads((output / "run.json").read_text())
            self.assertTrue(manifest["diagnostics"]["optical_limit"])
            data = np.loadtxt(output / "elf.csv", delimiter=",", skiprows=1)
            self.assertEqual(data.shape, (160, 2))

    def test_optical_line_csv_at_resonance(self):
        obj = json.loads(EXAMPLE.read_text(encoding="utf-8"))
        obj["material"]["oscillators"][0]["dispersion"] = 0
        obj["loss_grid_eV"] = {"min": 15, "max": 15, "step": 1}
        # At resonance L=A/(Gamma*T)=6. Without dispersion the q integral
        # factors into L times the logarithm of the kinematic bound ratio.
        log_ratio = 2 * math.log(math.sqrt(1000) + math.sqrt(985)) - math.log(15)
        bulk = 6 * log_ratio / (math.pi * 0.529177210903 * 1000)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "synthetic.json"
            path.write_text(json.dumps(obj), encoding="utf-8")
            for cmd, field, expected in (("elf", "elf_dimensionless", 6.),
                                         ("bulk", "diimfp_per_A_eV", bulk)):
                with self.subTest(cmd=cmd):
                    out = Path(tmp) / cmd
                    proc = self.call(cmd, path, "--out", out)
                    self.assertEqual(proc.returncode, 0, proc.stderr)
                    with (out / f"{cmd}.csv").open(encoding="utf-8") as stream:
                        row, = csv.DictReader(stream)
                    self.assertEqual(float(row["loss_eV"]), 15.)
                    self.assertAlmostEqual(float(row[field]), expected, delta=1e-13 * expected)
                    if cmd == "bulk":
                        self.assertGreater(float(row["quadrature_error_per_A_eV"]), 0)

    def test_static_sign_and_underflow_failure_status(self):
        cases = (
            ("xps-fixed-depth", (10, 110, 1, 0.5), 3., 0, {}, "error",
             ("--depth-A", "10", "--epsilon-momentum", "k_parallel")),
            ("elf", (1e-100, 1e-30, 1e-240, 0), 5e-101, 0, {}, "not_converged", ()),
            ("dielectric", (15, 1, 1e-6, 0), 1., 15.01,
             {"epsabs": 1e-6, "epsrel": 1e-5, "limit": 5}, "not_converged", ()),
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "synthetic.json"
            for cmd, line, loss, gap, integration, status, extra in cases:
                with self.subTest(cmd=cmd):
                    obj = json.loads(EXAMPLE.read_text(encoding="utf-8"))
                    obj["material"]["oscillators"] = [
                        dict(zip(("energy_eV", "strength_eV2", "width_eV", "dispersion"), line, strict=True))]
                    obj["material"]["gap_eV"] = gap
                    obj["loss_grid_eV"] = {"min": loss, "max": loss, "step": 1}
                    obj["integration"] = integration
                    path.write_text(json.dumps(obj), encoding="utf-8")
                    out = Path(tmp) / cmd
                    proc = self.call(cmd, path, *extra, "--out", out)
                    self.assertEqual(proc.returncode, 2, proc.stderr)
                    self.assertEqual(json.loads(proc.stderr)["status"], status)
                    self.assertFalse(out.exists())

    def test_dielectric_end_to_end(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "gapped.json"
            obj = json.loads(EXAMPLE.read_text())
            obj["material"]["gap_eV"] = 14.25
            path.write_text(json.dumps(obj), encoding="utf-8")
            output = root / "dielectric"
            r = self.call("dielectric", path, "--q-inv-A", 0.4, "--out", output)
            self.assertEqual(r.returncode, 0, r.stderr)
            manifest = json.loads((output / "run.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["model_id"], "causal_oscillator_elf_v1")
            self.assertEqual(manifest["input_config_model_id"], manifest["config"]["model"])
            self.assertEqual(manifest["input_config_model_id"], "bulk_drude_lindhard_v0")
            self.assertEqual(manifest["input_file"], str(path))
            self.assertEqual(manifest["input_sha256"], hashlib.sha256(path.read_bytes()).hexdigest())
            diagnostics = manifest["diagnostics"]
            self.assertNotIn("applicability", manifest)
            self.assertEqual(manifest["schema_version"], 1)
            self.assertEqual(diagnostics["static_sign"], "negative")
            self.assertNotIn("static_status", diagnostics)
            self.assertFalse(diagnostics["electron_energy_used"])
            self.assertEqual(diagnostics["q_inv_A"], 0.4)
            self.assertEqual(diagnostics["high_energy_tail"], "untruncated_oscillator_model")
            self.assertEqual(diagnostics["gap_edge_policy"], "reject_grid")
            self.assertEqual(diagnostics["gap_edge_exclusion_eV"], 14.25e-12)
            self.assertEqual(diagnostics["nearest_gap_edge_distance_eV"], 0.25)
            self.assertEqual(diagnostics["sub_gap_pole_status"], "located")
            self.assertTrue(0 < diagnostics["sub_gap_pole_eV"] < 14.25)
            self.assertTrue(manifest["warnings"])
            self.assertEqual(json.loads(r.stdout)["warning_count"], len(manifest["warnings"]))
            parse_config(manifest["config"])
            data = np.loadtxt(output / "dielectric.csv", delimiter=",", skiprows=1)
            self.assertEqual(data.shape, (160, 6))
            self.assertTrue(np.all(np.isfinite(data)))
            np.testing.assert_allclose((data[:, 1] + 1j * data[:, 2]) * (data[:, 3] + 1j * data[:, 4]), -1)
            np.testing.assert_array_equal(data[:, 0], np.arange(0.5, 80.01, 0.5))
            np.testing.assert_array_equal(data[data[:, 0] < 14.25, 2], np.zeros(28))
            self.assertGreater(data[-1, 2], 0)

    def test_dielectric_distinguishes_invalid_settings_from_quad_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "gapped.json"
            obj = json.loads(EXAMPLE.read_text())
            obj["material"]["gap_eV"] = 14.25
            cases = (
                ({"limit": 1}, "error", "too small for oscillator breakpoints"),
                ({"epsabs": 1e-30, "epsrel": 1e-15}, "not_converged", "failed to converge"),
            )
            for integration, status, message in cases:
                with self.subTest(integration=integration):
                    obj["integration"] = integration
                    path.write_text(json.dumps(obj), encoding="utf-8")
                    untouched = root / status
                    r = self.call("dielectric", path, "--out", untouched)
                    self.assertEqual(r.returncode, 2)
                    error = json.loads(r.stderr)
                    self.assertEqual(error["status"], status)
                    self.assertIn(message, error["message"])
                    self.assertFalse(untouched.exists())

    def test_dielectric_rejects_exact_and_rounded_gap_grid_with_remedy(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "gapped.json"
            obj = json.loads(EXAMPLE.read_text())
            grids = (
                (8.0, obj["loss_grid_eV"], 8.0),
                (3.3, {"min": 0.1, "max": 10, "step": 0.1}, 3.3000000000000003),
            )
            for gap, grid, offending in grids:
                with self.subTest(gap=gap):
                    obj["material"]["gap_eV"] = gap
                    obj["loss_grid_eV"] = grid
                    path.write_text(json.dumps(obj), encoding="utf-8")
                    output = root / "must-not-exist"
                    r = self.call("dielectric", path, "--out", output)
                    self.assertEqual(r.returncode, 2)
                    error = json.loads(r.stderr)
                    self.assertEqual(error["status"], "error")
                    self.assertEqual(error["error_type"], "DielectricGapEdgeError")
                    self.assertEqual(error["gap_eV"], gap)
                    self.assertEqual(error["offending_loss_eV"], [offending])
                    self.assertAlmostEqual(error["gap_edge_exclusion_eV"] / gap, 1e-12)
                    self.assertIn("Shift the loss grid", error["remedy"])
                    self.assertFalse(output.exists())

    def test_dielectric_warns_on_both_shifted_grids(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "gapped.json"
            obj = json.loads(EXAMPLE.read_text())
            obj["material"]["gap_eV"] = 3.3
            for offset, expected_real in ((1e-9, -12.2), (1e-4, 10.7)):
                with self.subTest(offset=offset):
                    obj["loss_grid_eV"] = {"min": 0.1 + offset, "max": 10 + offset, "step": 0.1}
                    path.write_text(json.dumps(obj), encoding="utf-8")
                    output = root / str(offset)
                    r = self.call("dielectric", path, "--out", output)
                    self.assertEqual(r.returncode, 0, r.stderr)
                    manifest = json.loads((output / "run.json").read_text(encoding="utf-8"))
                    point, = manifest["diagnostics"]["near_gap_edge_points"]
                    self.assertAlmostEqual(point["distance_eV"], offset, delta=1e-14)
                    self.assertGreater(point["edge_log_contribution"], 0)
                    warning, = [item for item in manifest["warnings"] if "grid offset" in item]
                    self.assertIn(f"loss_eV={point['loss_eV']:.17g}", warning)
                    self.assertIn(f"distance_eV={point['distance_eV']:.17g}", warning)
                    self.assertIn(f"edge_log_contribution={point['edge_log_contribution']:.9g}", warning)
                    data = np.loadtxt(output / "dielectric.csv", delimiter=",", skiprows=1)
                    edge_row = np.argmin(abs(data[:, 0] - 3.3))
                    self.assertAlmostEqual(data[edge_row, 3], expected_real, delta=0.1)

    def test_dielectric_gapless_omits_gap_edge_policy(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "gapless"
            r = self.call("dielectric", EXAMPLE, "--out", output)
            self.assertEqual(r.returncode, 0, r.stderr)
            manifest = json.loads((output / "run.json").read_text(encoding="utf-8"))
            self.assertNotIn("gap_edge_policy", manifest["diagnostics"])
            self.assertEqual(manifest["warnings"], [])
            self.assertEqual(json.loads(r.stdout)["warning_count"], 0)

    def test_unknown_and_invalid_config_rejected(self):
        base = json.loads(EXAMPLE.read_text())
        mutations = [
            lambda d: d.update(model="surface"),
            lambda d: d.update(surfaec=True),
            lambda d: d["material"].update(name=""),
            lambda d: d["material"].update(oscillators=[]),
            lambda d: d["material"]["oscillators"][0].update(strength_eV2=-1),
            lambda d: d["material"]["oscillators"][0].update(energy_eV=float("inf")),
            lambda d: d["material"]["oscillators"][0].update(width_eV=True),
            lambda d: d["loss_grid_eV"].update(min=0),
            lambda d: d["loss_grid_eV"].update(max=1000),
            lambda d: d["loss_grid_eV"].update(step=0.4),
            lambda d: d["integration"].update(epsrel=0),
            lambda d: d["integration"].update(limit=True),
        ]
        for i, mutate in enumerate(mutations):
            data = copy.deepcopy(base)
            mutate(data)
            with self.subTest(mutation=i):
                with self.assertRaises(ValueError):
                    parse_config(data)

    def test_lone_surrogates_rejected_before_creating_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "input.json"
            output = Path(tmp) / "output"
            obj = json.loads(EXAMPLE.read_text(encoding="utf-8"))
            obj["material"]["name"] = "lone\ud800"
            path.write_text(json.dumps(obj), encoding="ascii")
            for args in (("validate", path), ("bulk", path, "--out", output)):
                with self.subTest(command=args[0]):
                    result = self.call(*args)
                    self.assertEqual(result.returncode, 2)
                    self.assertEqual(result.stdout, "")
                    error = json.loads(result.stderr)
                    self.assertEqual((error["status"], error["error_type"]), ("error", "ValueError"))
                    self.assertIn("material.name", error["message"])
                    self.assertFalse(output.exists())

            obj["material"]["name"] = "paired \U0001f600"
            path.write_text(json.dumps(obj), encoding="ascii")
            self.assertEqual(load_config(path).material.name, obj["material"]["name"])
            with self.assertRaisesRegex(ValueError, "input path"):
                load_config(Path(tmp) / "bad\udcff.json")

    def test_nonfinite_and_duplicate_json_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.json"
            for text in ['{"model": "a", "model": "b"}', '{"x": NaN}', '{"x": Infinity}']:
                path.write_text(text)
                with self.assertRaises(ValueError):
                    load_config(path)
            r = self.call("validate", path)
            self.assertEqual(r.returncode, 2)
            self.assertEqual(json.loads(r.stderr)["status"], "error")


if __name__ == "__main__":
    unittest.main()
