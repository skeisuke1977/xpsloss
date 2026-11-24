import contextlib
from dataclasses import FrozenInstanceError, replace
import hashlib
import io
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from xpsloss import (
    Material, Oscillator, bulk_diimfp, effective_cross_section_fixed_depth as calculate,
    loss_function, negative_inverse_dielectric_function,
)
from xpsloss import cli, model, xps_fixed_depth as fixed
from xpsloss.config import parse_config


ROOT = Path(__file__).resolve().parents[1]
RECORD = ROOT / "benchmarks/inputs/al_elfdb_v1.0.json"
SOURCE = json.loads(RECORD.read_text(encoding="utf-8"))["material"]
AL = Material(SOURCE["name"], SOURCE["gap_eV"], tuple(Oscillator(**o) for o in SOURCE["oscillators"]))
SYNTHETIC_LOCAL = Material("synthetic dispersionless Al variant", 0,
                           tuple(replace(o, dispersion=0) for o in AL.oscillators))
GRID = 0.2 * np.arange(1, 500)
MOMENTA = ("k_parallel", "k_total")
COMPONENTS = ("total", "no_hole", "hole", *fixed.TERM_NAMES)


class FixedDepthPhysicsTests(unittest.TestCase):
    def test_i1_bulk_limit_and_i2_dispersive_bulk_identity(self):
        # The scalar bulk reference dominates test time. Sample at each eV
        # plus the fine low-loss and optical-resonance regions of the case grid.
        indices = np.unique(np.r_[np.arange(5, 501, 5), 5:16, 51:60, 73:84, 499])
        grid = 0.2 * indices
        bulk = np.array([bulk_diimfp(AL, 1000, t, epsabs=1e-13, epsrel=1e-11)[0] for t in grid])
        shallow, deep = (calculate(AL, 1000, grid, a, epsilon_momentum="k_total") for a in (1e4, 1e5))
        for result in (shallow, deep):
            np.testing.assert_allclose(result.term_13a, bulk, rtol=1e-8, atol=0)
        for name in ("total", "no_hole"):
            with self.subTest(component=name):
                x1 = shallow.path_length_A * (getattr(shallow, name) - bulk)
                x2 = deep.path_length_A * (getattr(deep, name) - bulk)
                np.testing.assert_allclose(x1, x2, rtol=1e-6, atol=0)

    def test_i2_dispersionless_analytic_bulk(self):
        grid = np.append(GRID, 100)
        log_ratio = 2 * np.log(math.sqrt(1000) + np.sqrt(1000-grid)) - np.log(grid)
        expected = loss_function(SYNTHETIC_LOCAL, grid) * log_ratio / (math.pi * model.BOHR_RADIUS_A * 1000)
        for choice in MOMENTA:
            result = calculate(SYNTHETIC_LOCAL, 1000, grid, 10, epsilon_momentum=choice)
            np.testing.assert_allclose(result.term_13a, expected, rtol=1e-10, atol=0)

    def test_i4_zero_response(self):
        zero = Material("synthetic zero strength", 0, (replace(AL.oscillators[0], strength_eV2=0),))
        for choice in MOMENTA:
            result = calculate(zero, 1000, [0.2, 11, 15.6, 30, 99.8], 2, epsilon_momentum=choice)
            for name in (*COMPONENTS, "quadrature_error"):
                np.testing.assert_array_equal(getattr(result, name), 0)

    def test_i6_total_nonnegative_on_full_grid(self):
        # Signed no_hole and hole ranges for these runs are recorded in the task report.
        for depth in (2, 10, 50, 200):
            for choice in MOMENTA:
                with self.subTest(depth_A=depth, epsilon_momentum=choice):
                    result = calculate(AL, 1000, GRID, depth, epsilon_momentum=choice)
                    self.assertTrue(np.all(result.total >= 0))

    def test_i7_double_log_panels_at_both_momentum_choices(self):
        original = fixed._log_breakpoints

        def bisect_every_panel(*args):
            edges = original(*args)
            return np.sort(np.r_[edges, (edges[:-1] + edges[1:]) / 2])

        for choice in MOMENTA:
            base = calculate(AL, 1000, GRID, 10, epsilon_momentum=choice)
            with patch.object(fixed, "_log_breakpoints", bisect_every_panel):
                refined = calculate(AL, 1000, GRID, 10, epsilon_momentum=choice)
            for name in COMPONENTS:
                diff = abs(getattr(base, name) - getattr(refined, name))
                self.assertLessEqual(diff.max(), 1e-8 * base.total.max(), (choice, name))

    def test_i8_energy_and_depth_scaling_and_unit_conversion(self):
        # With E,T,Γ scaled by c, q scales by sqrt(c). Scaling a by 1/sqrt(c)
        # preserves ka and Ωa. The response is invariant for A -> c² A.
        grid = np.array([1, 11, 15.6, 30, 99.8])
        for choice in MOMENTA:
            base = calculate(AL, 1000, grid, 10, epsilon_momentum=choice)
            for energy_factor, length_factor, convert_units in ((4, 0.5, False), (1000, 0.1, True)):
                oscillators = tuple(replace(o, energy_eV=o.energy_eV * energy_factor,
                                            width_eV=o.width_eV * energy_factor,
                                            strength_eV2=o.strength_eV2 * energy_factor**2)
                                    for o in AL.oscillators)
                material = Material("synthetic scaled oscillator", 0, oscillators)
                with contextlib.ExitStack() as context:
                    factor = energy_factor
                    if convert_units:
                        # Same physical system, using meV and nm throughout.
                        h = model.HBAR2_OVER_2ME_EV_A2 * energy_factor * length_factor**2
                        context.enter_context(patch.object(model, "HBAR2_OVER_2ME_EV_A2", h))
                        context.enter_context(patch.object(fixed, "HBAR2_OVER_2ME_EV_A2", h))
                        context.enter_context(patch.object(fixed, "BOHR_RADIUS_A", model.BOHR_RADIUS_A * length_factor))
                        factor *= length_factor
                    result = calculate(material, 1000 * energy_factor, grid * energy_factor,
                                       10 * length_factor, epsilon_momentum=choice)
                for name in COMPONENTS:
                    np.testing.assert_allclose(getattr(result, name) * factor, getattr(base, name),
                                               rtol=1e-10, atol=1e-17, err_msg=name)

    def test_e4b_agreement_with_dielectric_at_arbitrary_q(self):
        material = Material("synthetic two oscillators", 0, (*AL.oscillators, Oscillator(30, 10, 1.1, 0.3)))
        q = np.array([0, 1e-5, 0.1, 1.7, 5, 32])
        for loss in (0.2, 11, 15.6, 30, 90):
            expected = negative_inverse_dielectric_function(material, loss, q)
            # The scalar/refinement path shares its line evaluator with dielectric.
            # The batch evaluator forms its own bounded real-arithmetic denominator.
            for evaluate in (fixed._oscillator_sum, fixed._bounded_sum):
                response, _, _ = evaluate(material, loss, q)
                np.testing.assert_allclose(response - 1, expected, rtol=2e-14, atol=1e-15)

    def test_surface_line_for_sum_of_narrow_oscillators(self):
        # Partitioning a strength between coincident oscillators gives the
        # same dielectric function, including the location of the surface line.
        single = Material("synthetic narrow", 0, (Oscillator(15, 200, 0.001, 1),))
        split = Material("synthetic split", 0, tuple(Oscillator(15, s, 0.001, 1) for s in (80, 120)))
        for choice in MOMENTA:
            a, b = (calculate(m, 1000, [13, 14, 17, 20, 80], 10, epsilon_momentum=choice)
                    for m in (single, split))
            for name in COMPONENTS:
                np.testing.assert_allclose(getattr(a, name), getattr(b, name), rtol=1e-8, atol=1e-12)

    def test_small_cartesian_surface_coefficients(self):
        # Independent epsilon-polynomial references in generate_references_4f.py.
        # Re C nearly vanishes at S=1.5; Im R is small near a large bulk pole.
        for t, width, index, part, ref in (
                (math.sqrt(105.), 3e-7, 2, "real", 3.01615238848569e-14),
                (math.nextafter(15., 0.), 3e-18, 1, "imag", 5.000000000000005e-19)):
            with self.subTest(loss=t, width=width):
                mat = Material("synthetic", 0, (Oscillator(15, 180, width),))
                values, estimates = fixed._coefficients(mat, t, np.array([0.]), 1e-9)
                error = abs(getattr(values[index][0], part) - ref)
                self.assertLessEqual(error, 1e-9 * abs(ref))
                self.assertLessEqual(error, estimates[index][int(part == "imag")][0] + math.ulp(ref))

    def test_refining_adjacent_panels_preserves_their_shared_boundary(self):
        # Split close to the bulk peak, forcing adaptive refinements on both
        # sides of the same endpoint. Eq. (1) is an independent integral in q.
        def two_panels(material, loss, lower, upper, offset):
            return np.log([lower, 2, upper])

        with patch.object(fixed, "_log_breakpoints", two_panels):
            result = calculate(AL, 1000, 30, 10, epsilon_momentum="k_total")
        reference, _ = bulk_diimfp(AL, 1000, 30, epsabs=1e-13, epsrel=1e-11)
        np.testing.assert_allclose(result.term_13a, reference, rtol=1e-8, atol=0)

    def test_parallel_limits_at_small_loss_and_line_breakpoints(self):
        t = 1e-14
        low, high = fixed._parallel_limits(1000, t)
        # Leading small-loss limit, derived from k_minus = Ω(1 + T/(4E) + ...).
        omega = t / (2 * math.sqrt(model.HBAR2_OVER_2ME_EV_A2 * 1000))
        self.assertAlmostEqual(low / (omega * math.sqrt(t/2000)), 1, places=14)
        self.assertAlmostEqual(high / (2*math.sqrt(1000/model.HBAR2_OVER_2ME_EV_A2)), 1, places=14)
        loss = 30
        lo, hi = fixed._parallel_limits(1000, loss)
        for offset in (0, loss / (2 * math.sqrt(model.HBAR2_OVER_2ME_EV_A2 * 1000))):
            line = math.sqrt((loss - AL.oscillators[0].energy_eV) / model.HBAR2_OVER_2ME_EV_A2 - offset**2)
            edges = fixed._log_breakpoints(AL, loss, lo, hi, offset)
            self.assertLess(np.min(abs(edges - math.log(line))), 1e-14)

    def test_api_domain_and_required_choice(self):
        args = dict(material=AL, electron_energy_eV=1000, loss_eV=11, depth_A=10)
        with self.assertRaises(TypeError):
            calculate(**args)
        args["epsilon_momentum"] = "k_total"
        self.assertTrue(np.isfinite(calculate(**(args | {"material": replace(AL, gap_eV=1)})).total[0]))
        for key, value, message in (
            ("theta_deg", 1e-12, "Q3/Q15"),
            ("theta_deg", float("inf"), "finite"), ("epsilon_momentum", None, "epsilon_momentum"),
            ("epsilon_momentum", "optical", "epsilon_momentum"), ("depth_A", 0, "positive"),
            ("depth_A", -1, "positive"), ("depth_A", float("nan"), "finite"),
            ("loss_eV", [5, 0], "delta"), ("loss_eV", -0.2, "loss_eV"),
            ("loss_eV", [1, 1000], "loss_eV"), ("loss_eV", 1001, "loss_eV"),
            ("loss_eV", float("nan"), "loss_eV"), ("loss_eV", [], "nonempty"),
            ("loss_eV", [1j], "real"), ("loss_eV", [[1, 3]], "dimensional"),
            ("electron_energy_eV", float("inf"), "finite"), ("electron_energy_eV", -1, "positive"),
            ("hole_charge", float("inf"), "finite"), ("epsabs", 0, "positive"),
            ("hole_charge", np.complex128(0.5+1j), "real"),
            ("epsrel", float("nan"), "finite"), ("limit", True, "integer"),
        ):
            with self.subTest(key=key, value=value), self.assertRaisesRegex(ValueError, message):
                calculate(**(args | {key: value}))
        result = calculate(**args)
        self.assertEqual(result.total.shape, (1,))
        with self.assertRaises(FrozenInstanceError):
            result.depth_A = 0


class FixedDepthCLITests(unittest.TestCase):
    def call(self, *args):
        process = subprocess.run([sys.executable, "-m", "xpsloss", *map(str, args)],
                                 env=dict(os.environ, PYTHONPATH=str(ROOT / "src")),
                                 capture_output=True, text=True, check=False)
        if process.returncode:
            self.assertFalse(process.stdout)
            return process.returncode, json.loads(process.stderr)
        self.assertFalse(process.stderr)
        return 0, json.loads(process.stdout)

    def test_cli_manifest_csv_and_overwrite(self):
        with tempfile.TemporaryDirectory(dir=ROOT) as tmp:
            source = ROOT / "benchmarks/inputs" / "al_elfdb_v1.0.bulk.json"
            config = json.loads(source.read_text(encoding="utf-8"))
            config["loss_grid_eV"] = {"min": 5, "max": 30, "step": 5}
            config["integration"] = {"epsabs": 1e-11, "epsrel": 1e-8, "limit": 150}
            path, out = Path(tmp)/"config.json", Path(tmp)/"out"
            path.write_text(json.dumps(config), encoding="utf-8")
            args = ("xps-fixed-depth", path, "--depth-A", 10, "--theta-deg", 0,
                    "--hole-charge", 0.5, "--epsilon-momentum", "k_parallel", "--out", out)
            code, success = self.call(*args)
            self.assertEqual(code, 0, success)
            self.assertEqual(success["files"], ["xps_fixed_depth.csv", "run.json"])
            manifest = json.loads((out/"run.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["model_id"], "simonsen1997_fixed_depth_v0")
            self.assertEqual(manifest["input_config_model_id"], "bulk_drude_lindhard_v0")
            self.assertEqual(manifest["input_sha256"], hashlib.sha256(path.read_bytes()).hexdigest())
            self.assertFalse(manifest["queels_xps_equivalence"])
            self.assertEqual(manifest["warnings"], [])
            self.assertIn("https://doi.org/10.1103/PhysRevB.56.1612", manifest["references"])
            self.assertEqual(parse_config(manifest["config"]).material, AL)
            diag = manifest["diagnostics"]
            self.assertNotIn("velocity", diag)
            self.assertEqual(manifest["schema_version"], 1)
            for key, expected in (("depth_A", 10), ("path_length_A", 10), ("theta_deg", 0),
                                  ("hole_charge", 0.5), ("epsilon_momentum", "k_parallel"),
                                  ("kinematics", "nonrelativistic")):
                self.assertEqual(diag[key], expected)
            self.assertEqual(set(diag["component_definitions"]), {"total", "no_hole", "hole"})
            csv = out / "xps_fixed_depth.csv"
            before = csv.read_bytes()
            data = np.genfromtxt(csv, delimiter=",", names=True)
            self.assertEqual(data.dtype.names, ("loss_eV", *COMPONENTS, "quadrature_error"))
            expected = calculate(AL, 1000, np.arange(5, 31, 5), 10, hole_charge=0.5,
                                 epsilon_momentum="k_parallel", **config["integration"])
            for name in data.dtype.names:
                np.testing.assert_array_equal(data[name], getattr(expected, name))
                self.assertEqual(manifest["units"][name], "eV" if name == "loss_eV" else "angstrom^-1 eV^-1")
            np.testing.assert_array_equal(data["total"], data["no_hole"] + 0.5*data["hole"])
            code, error = self.call(*args)
            self.assertEqual((code, error["error_type"]), (2, "FileExistsError"))
            self.assertEqual(csv.read_bytes(), before)

    def test_cli_failures_and_model_scopes(self):
        code, cap = self.call("capabilities", "--json")
        self.assertEqual(code, 0)
        scope = cap["xps_fixed_depth"]
        for name in ("surface_losses", "core_hole", "vacuum_losses"):
            self.assertIn(name, scope["supported"])
        self.assertNotIn("gapped_materials", scope["unsupported_physics"])
        for name in ("oblique_emission", "depth_averaging", "elastic_scattering",
                     "sub_gap_surface_mode_delta_terms", "gap_broadening",
                     "multiple_inelastic_scattering", "primary_spectrum_convolution"):
            self.assertIn(name, scope["unsupported_physics"])
        with tempfile.TemporaryDirectory(dir=ROOT) as tmp:
            path, out = Path(tmp)/"config.json", Path(tmp)/"out"
            config = json.loads((ROOT / "examples/synthetic_bulk.json").read_text(encoding="utf-8"))
            config["loss_grid_eV"] = {"min": 30, "max": 30, "step": 1}
            args = ("xps-fixed-depth", path, "--depth-A", 10, "--out", out)
            for options, integration, status in (
                ((), {}, "error"), (("--epsilon-momentum", "bad"), {}, "error"),
                (("--epsilon-momentum", "k_total", "--theta-deg", "1"), {}, "error"),
                (("--epsilon-momentum", "k_total"), {"limit": 1}, "error"),
                (("--epsilon-momentum", "k_total"), {"epsabs": 1e-30, "epsrel": 1e-30}, "not_converged"),
            ):
                config["integration"] = integration
                path.write_text(json.dumps(config), encoding="utf-8")
                with contextlib.redirect_stderr(io.StringIO()) as stderr:
                    code = cli.main(list(map(str, (*args, *options))))
                error = json.loads(stderr.getvalue())
                self.assertEqual((code, error["status"]), (2, status), error)
                self.assertFalse(out.exists())
