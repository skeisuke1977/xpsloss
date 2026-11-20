import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from xpsloss.cli import capabilities


ROOT = Path(__file__).resolve().parents[1]


class FixedDepthApplicabilityTests(unittest.TestCase):
    def manifest(self, energy, first_loss, last_loss, momentum):
        config = json.loads((ROOT / "examples/synthetic_bulk.json").read_text(encoding="utf-8"))
        config["electron_energy_eV"] = energy
        config["loss_grid_eV"] = {"min": first_loss, "max": last_loss, "step": last_loss - first_loss}
        with tempfile.TemporaryDirectory() as tmp:
            path, out = Path(tmp) / "input.json", Path(tmp) / "out"
            path.write_text(json.dumps(config), encoding="utf-8")
            proc = subprocess.run(
                [sys.executable, "-m", "xpsloss", "xps-fixed-depth", str(path), "--depth-A", "10",
                 "--epsilon-momentum", momentum, "--out", str(out)],
                cwd=ROOT, env=dict(os.environ, PYTHONPATH=str(ROOT / "src")), capture_output=True, text=True,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            summary = json.loads(proc.stdout)
            self.assertEqual(summary["status"], "complete")
            manifest = json.loads((out / "run.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["status"], "complete")
            self.assertEqual(summary["warning_count"], len(manifest["warnings"]))
            return manifest

    def test_high_loss_completes_with_quantified_warning(self):
        manifest = self.manifest(1000, 30, 900, "k_total")
        record = manifest["applicability"]
        self.assertEqual(record["max_loss_over_E0"], 0.9)
        self.assertIs(record["within_checked_envelope"], False)
        warning, = manifest["warnings"]
        self.assertIn("The numerical checks in this repository do not cover this input, "
                      "and physical applicability is not verified.", warning)
        self.assertIn("max_loss_over_E0=0.9 exceeds 0.1 by 0.8", warning)

    def test_checked_bounds_and_assumptions(self):
        for energy, last_loss, momentum in ((1000, 30, "k_parallel"), (4000, 400, "k_total")):
            with self.subTest(energy=energy):
                manifest = self.manifest(energy, 5, last_loss, momentum)
                record = manifest["applicability"]
                self.assertIs(record["within_checked_envelope"], True)
                self.assertEqual(manifest["warnings"], [])
                self.assertEqual(record["max_loss_over_E0"], last_loss / energy)
                self.assertEqual(record["electron_energy_eV"], energy)
                self.assertEqual(record["epsilon_momentum"], momentum)
                self.assertEqual(record["assumptions"], {
                    "velocity": "constant_at_initial_kinetic_energy", "kinematics": "nonrelativistic",
                    "trajectory": "rectilinear", "core_hole": "static", "small_loss": "T << E0",
                    "pair_creation": "sudden_at_t0", "core_hole_lifetime": "infinite",
                    "surface": "specular_reflection_model",
                    "source": "Simonsen et al. (1997) p. 1613 and Eq. (3); docs/M1_EQUATION_MAP.md section 2",
                })
                self.assertIn("Only max(T)/E0 and E0 are compared with the bounds", record["note"])
                self.assertIn("neither accuracy for arbitrary inputs nor physical validity", record["note"])
                bounds = record["checked_envelope"]
                for bound in bounds.values():
                    self.assertIn("test_applicability_envelope.py::EnvelopeEndpointTests::"
                                  "test_terms_and_components_at_both_energies", bound["source"])
                    self.assertIn("independent route I3 checked all eight terms", bound["basis"])
                    self.assertIn("For dispersive k_total no independent check covers all finite-depth terms "
                                  "at any E0", bound["basis"])
                self.assertEqual(bounds["max_loss_over_E0"]["max"], 0.1)
                self.assertIn("test_i2_dispersionless_analytic_bulk", bounds["max_loss_over_E0"]["supplementary"])
                energies = bounds["electron_energy_eV"]
                self.assertEqual((energies["min"], energies["max"]), (1000, 4000))
                self.assertIn("test_i8_energy_and_depth_scaling_and_unit_conversion", energies["supplementary"])

    def test_energy_departures_and_simultaneous_loss_departure(self):
        cases = ((800, 40, "below 1000 eV by 200 eV", 1),
                 (4500, 400, "exceeds 4000 eV by 500 eV", 1),
                 (800, 720, "below 1000 eV by 200 eV", 2))
        for energy, last_loss, departure, count in cases:
            with self.subTest(energy=energy, last_loss=last_loss):
                manifest = self.manifest(energy, 20, last_loss, "k_total")
                self.assertIs(manifest["applicability"]["within_checked_envelope"], False)
                warnings = manifest["warnings"]
                self.assertEqual(len(warnings), count)
                self.assertTrue(all("physical applicability is not verified" in w for w in warnings))
                self.assertTrue(any("electron_energy_eV=" in w and departure in w for w in warnings))
                if count == 2:
                    self.assertTrue(any("max_loss_over_E0=0.9 exceeds 0.1 by 0.8" in w for w in warnings))

    def test_capabilities_identify_the_scope_of_each_check(self):
        verification = capabilities()["xps_fixed_depth"]["verification"]
        momentum = verification["epsilon_momentum"]
        self.assertIn("I3", momentum["k_parallel"])
        self.assertIn("dispersionless and dispersive", momentum["k_parallel"])
        for phrase in ("(13a)", "deep limit", "same closed forms", "dispersive inputs",
                       "no independent check", "charge and potential equations", "I3 covers both"):
            self.assertIn(phrase, momentum["k_total"])
        self.assertIn("applicability.checked_envelope", verification["checked_envelope"])
        figures = verification["published_figures"]
        self.assertIn("attempted to reproduce", figures)
        self.assertIn("outside this package's oscillator model", figures)
        self.assertIn("Of 84 comparison points, 36 are within reading uncertainty, "
                      "18 are outside and 30 depend on width assumptions", figures)
        self.assertIn("37/19/28 for k_parallel and 37/18/29 for k_total", figures)
        self.assertIn("benchmarks/paper_figures/RESULTS.md in the source distribution or repository", figures)
