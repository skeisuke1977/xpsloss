"""Preflight priority, structured CLI failures and the hard-gap manifest."""

from contextlib import redirect_stderr, ExitStack
from dataclasses import asdict, replace
from fractions import Fraction
import io
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from test_m22_physics import SI, MODES, NAMES, SETTINGS, GRID, ROOT, run_fixed
from xpsloss import cli, FixedDepthConvergenceError, FixedDepthGapDomainError, Material, Oscillator
from xpsloss import xps_fixed_depth as fixed
from xpsloss._gap_domain import _lower_q_squared, _monotonicity, _proof, _uniform_bounds
from xpsloss.dielectric import DielectricGapEdgeError, _gap_edge_exclusion_eV


def synthetic(label, gap, line):
    return Material("synthetic " + label, gap, (Oscillator(*line),))


P1 = synthetic("P1", 2, (3, 8.1, 0.5, 0.5))
P0 = synthetic("P0", 2, (3, 8.1, 0.5, 0))
M1 = synthetic("M1", 2, (3, 8.1, 1, 0))
UE = synthetic("UE", 2, (1.9, 1, 0.05, 20))
SP = synthetic("SP", 1, (5, 40, 1, 0))
N1 = synthetic("N1", 2, (3, 6, 0.02, 0.5))


def cli_args(folder, material, losses, mode="k_total", energy=1000.0, settings=SETTINGS):
    config = {
        "material": asdict(material),
        "electron_energy_eV": energy,
        "loss_grid_eV": {
            "min": float(losses[0]),
            "max": float(losses[-1]),
            "step": float(losses[1] - losses[0]) if len(losses) > 1 else 1.0,
        },
        "integration": settings,
    }
    p = Path(folder) / "input.json"
    p.write_text(json.dumps(config, allow_nan=False), encoding="utf8")
    return [
        "xps-fixed-depth",
        str(p),
        "--out",
        str(Path(folder) / "output"),
        "--depth-A",
        "10",
        "--epsilon-momentum",
        mode,
    ]


class DomainTests(unittest.TestCase):
    def test_global_q_minimum_for_monotonicity_and_static_sign(self):
        mono = synthetic("dispersive M1", 2, (3, 8.1, 1, 1))
        static = synthetic("dispersive SP", 1, (5, 40, 1, 2))
        for material in (mono, static):
            result = run_fixed(material, 60.0, "k_total")
            self.assertTrue(result.diagnostics["monotonicity"][0]["passed"])
            self.assertEqual(result.diagnostics["static_sign"], "negative")
        with self.assertRaises(FixedDepthGapDomainError) as failed:
            run_fixed(mono, [2.5, 60.0], "k_total")
        self.assertEqual(failed.exception.reason, "monotonicity_condition_failed")
        with self.assertRaisesRegex(ValueError, "exceeds the static sign tolerance"):
            run_fixed(static, [2.0, 60.0], "k_total")

    def test_d7_equality_and_adjacent_exterior_floats(self):
        # For Si the nominal g +/- delta rounds inside the boundary. In this
        # synthetic material delta is exactly 2^-40, representable on either side.
        exact_gap = math.ldexp(1e12, -40)
        exact = replace(SI, name="synthetic exact D7 boundary", gap_eV=exact_gap)
        for mat in (exact, SI):
            g, d = mat.gap_eV, _gap_edge_exclusion_eV(mat)
            inside = [g - d, g + d]
            exterior = [math.nextafter(inside[0], -math.inf), math.nextafter(inside[1], math.inf)]
            for t in inside:
                if mat is exact:
                    self.assertEqual(abs(t - g), d)
                else:
                    self.assertLess(abs(t - g), d)
            self.assertTrue(all(abs(t - g) > d for t in exterior))
            for mode in MODES:
                with self.subTest(material=mat.name, mode=mode):
                    with self.assertRaises(DielectricGapEdgeError) as caught:
                        run_fixed(mat, inside, mode)
                    self.assertEqual(caught.exception.offending_loss_eV, inside)
                    result = run_fixed(mat, exterior, mode)
                    self.assertTrue(np.all(np.isfinite(result.total)))
                    self.assertTrue(result.diagnostics["subgap_proofs"][0]["proved"])

    def test_proof_log_outward_at_acceptance_threshold(self):
        # ln(2) = 2 atanh(1/3); its positive rational partial sum is strictly
        # below ln(2), with tail < 1e-80. Put the true bound just ABOVE 2.
        # Checking only the final binary64 bound far from 2 hides a decimal ulp.
        ln2_lower = 2 * sum(Fraction(1, (2*k + 1) * 3**(2*k + 1)) for k in range(85))
        for power in (2, 3):
            loss = 1.0 - 2.0**(1 - power)
            # Scale the coefficient for ln(8) to keep the constant part positive.
            c = Fraction(1, 2)
            p = 2 - c * power * ln2_lower
            proof = _proof([(1.5, p, c)], 1.0, loss)
            self.assertFalse(proof["proved"])
            self.assertGreater(proof["S_upper"], 2.0)
            self.assertLess(proof["margin_lower"], 0.0)

    def test_uniform_edge_coefficient_bounds_pi_outward(self):
        # Truncated decimal expansion, strictly less than pi. The exact
        # rational M*g/pi_lower resolves the two adjacent binary64 pi values.
        pi_lower = Fraction("3.14159265358979323846264338327950288")
        line = P1.oscillators[0]
        for u, _, coefficient in _uniform_bounds(P1):
            moment = Fraction(line.strength_eV2) * Fraction(line.width_eV) / (
                Fraction(line.energy_eV)**2 - Fraction(u)**2
            )**2
            self.assertGreaterEqual(coefficient, moment * Fraction(P1.gap_eV) / pi_lower)

    def test_domain_details_keep_numeric_fields_separate_from_overflow_reason(self):
        extreme = synthetic("extreme line", 1, (1e-100, 1, 1, 0))
        with self.assertRaises(FixedDepthGapDomainError) as caught:
            run_fixed(extreme, 2.0, "k_total")
        row, = caught.exception.details["oscillators"]
        self.assertIsNone(row["lhs_lower"])
        self.assertEqual(row["lhs_lower_status"], "outside_binary64_range")
        self.assertIsInstance(row["rhs_upper"], float)
        json.dumps(caught.exception.details, allow_nan=False)

    def test_epsilon_pole_diagnostic_failure_preserves_result(self):
        normal = run_fixed(SI, [0.2, 4.0], "k_total")
        for failure in (fixed.NumericalConvergenceError("D3"), RuntimeError("iterations"), ValueError("bracket")):
            with self.subTest(failure=failure), patch("xpsloss._gap_domain._sub_gap_pole", side_effect=failure):
                result = run_fixed(SI, [0.2, 4.0], "k_total")
                self.assertEqual(result.diagnostics["sub_gap_pole_status"], "unresolved")
                self.assertIsNone(result.diagnostics["sub_gap_pole_eV"])
                self.assertEqual(result.diagnostics["ill_conditioned_D3_count"], 0)
                for name in (*NAMES, "quadrature_error"):
                    np.testing.assert_array_equal(getattr(result, name), getattr(normal, name))

    def test_v12_domain_and_cli_before_quadrature(self):
        g, delta = SI.gap_eV, _gap_edge_exclusion_eV(SI)
        for mat, t, kind, reason in [
            *((SI, t, DielectricGapEdgeError, None) for t in (g, g + delta / 2, g - delta / 2)),
            (P1, 1.999, FixedDepthGapDomainError, "surface_pole_not_excluded"),
            (P1, 1.9, FixedDepthGapDomainError, "surface_pole_not_excluded"),
            (P0, 1.999, FixedDepthGapDomainError, "surface_pole_not_excluded"),
            (M1, 3.0, FixedDepthGapDomainError, "monotonicity_condition_failed"),
            (UE, 1.5, FixedDepthGapDomainError, "surface_pole_not_excluded"),
            (SP, 1.5, ValueError, None),
        ]:
            modes = ("k_parallel",) if mat is UE else MODES
            for mode in modes:
                with self.subTest(mat=mat.name, t=t, mode=mode), ExitStack() as stack:
                    # Both batch and single-point fixed-depth paths use this evaluator.
                    stack.enter_context(
                        patch.object(fixed, "_panel_integrals", side_effect=AssertionError("early node"))
                    )
                    energy = 3.0 if mat is UE else 1000.0
                    with self.assertRaises(kind) as caught:
                        fixed.effective_cross_section_fixed_depth(
                            mat, energy, t, 10.0, epsilon_momentum=mode, **SETTINGS
                        )
                    if reason:
                        e = caught.exception
                        self.assertEqual(e.reason, reason)
                        self.assertEqual(e.offending_loss_eV, [t])
                        self.assertEqual(e.details["points"][0]["index"], 0)
                        self.assertEqual(len(e.details["points"][0]["q_interval_inv_A"]), 2)
                        if reason.startswith("monotonicity"):
                            self.assertFalse(e.details["oscillators"][0]["passed"])
                        elif mat is UE:
                            self.assertIsNone(e.details["points"][0]["U_eV"])
                        else:
                            self.assertGreaterEqual(e.details["points"][0]["S_upper"], 2.0)
                    elif mat is SP:
                        self.assertIn("exceeds the static sign tolerance", str(caught.exception))
                    folder = stack.enter_context(tempfile.TemporaryDirectory())
                    args = cli_args(folder, mat, [t], mode, energy)
                    message = io.StringIO()
                    with redirect_stderr(message):
                        self.assertEqual(cli.main(args), 2)
                    self.assertTrue(message.getvalue().isascii())
                    doc = json.loads(message.getvalue())
                    self.assertEqual(doc["status"], "error")
                    self.assertEqual(doc["error_type"], kind.__name__)
                    if reason:
                        self.assertEqual(doc["reason"], reason)
                        self.assertEqual(doc["details"], caught.exception.details)
                    self.assertFalse((Path(folder) / "output").exists())

    def test_v12_order_mixed_batches_and_duplicates(self):
        d = _gap_edge_exclusion_eV(P1)
        losses = [3.0] * 8 + [2 + d / 2, 1.999, 2.0, 2 - d / 2, 2.0]
        with patch.object(fixed, "_parallel_limits", side_effect=AssertionError("D7 must precede kinematics")):
            with self.assertRaises(DielectricGapEdgeError) as c:
                run_fixed(P1, losses, "k_parallel")
        self.assertEqual(c.exception.offending_indices, [8, 10, 11, 12])
        self.assertEqual(c.exception.offending_loss_eV, [losses[i] for i in [8, 10, 11, 12]])
        with patch.object(fixed, "_panel_integrals", side_effect=AssertionError("early node")):
            with self.assertRaises(FixedDepthGapDomainError) as c:
                run_fixed(P1, [4.0] * 8 + [1.999, 3.0, 1.9, 1.999], "k_total")
            self.assertEqual(c.exception.offending_loss_eV, [1.999, 1.9, 1.999])
            self.assertEqual([p["index"] for p in c.exception.details["points"]], [8, 10, 11])
            with self.assertRaises(FixedDepthConvergenceError) as c:
                run_fixed(M1, [3.0, math.ulp(0.0), 4.0, math.ulp(0.0)], "k_parallel")
            self.assertIn("index=(1,)", str(c.exception))
            self.assertIn("index=(3,)", str(c.exception))
            with self.assertRaisesRegex(ValueError, "static sign tolerance"):
                run_fixed(SP, [0.99999], "k_total")
            bad = Material("synthetic two invalid lines", 2, M1.oscillators * 2)
            with self.assertRaises(FixedDepthGapDomainError) as c:
                run_fixed(bad, [1.9, 3.0, 1.9], "k_total")
            self.assertEqual(c.exception.reason, "monotonicity_condition_failed")
            self.assertEqual(
                [o["oscillator_index"] for o in c.exception.details["oscillators"] if not o["passed"]], [0, 1]
            )
            self.assertEqual(c.exception.offending_loss_eV, [1.9, 3.0, 1.9])
        with tempfile.TemporaryDirectory() as folder:
            args = cli_args(folder, P1, [2.0])
            out = Path(folder) / "output"
            out.mkdir()
            (out / "existing").write_bytes(b"unchanged")
            buffer = io.StringIO()
            with redirect_stderr(buffer):
                self.assertEqual(cli.main(args), 2)
            self.assertEqual(json.loads(buffer.getvalue())["error_type"], "FileExistsError")
            self.assertEqual((out / "existing").read_bytes(), b"unchanged")

    def test_v12_ill_conditioned_and_precision_failures(self):
        tiny = dict(epsabs=1e-30, epsrel=1e-15, limit=200)
        # A dispersionless line needs no surface search. Exactly on its
        # resonance the D3 condition estimate fails after successful preflight.
        stiff = synthetic("D3 failure", 2.0, (4.0, 4.0, 1e-6, 0.0))
        for mat, t, settings, expected in ((SI, 5.0, tiny, "budget"), (stiff, 4.0, SETTINGS, "D3")):
            with self.subTest(mat=mat.name), self.assertRaisesRegex(FixedDepthConvergenceError, expected):
                run_fixed(mat, t, "k_total", **settings)
            with tempfile.TemporaryDirectory() as folder:
                args = cli_args(folder, mat, [t], settings=settings)
                proc = subprocess.run([sys.executable, "-m", "xpsloss", *args], capture_output=True, cwd=ROOT)
                self.assertEqual(proc.returncode, 2)
                message = json.loads(proc.stderr.decode("ascii"))
                self.assertEqual(message["status"], "not_converged")
                self.assertEqual(message["diagnostics"]["ill_conditioned_D3_count"], int(mat is stiff))
                self.assertFalse((Path(folder) / "output").exists())

    def test_monotonicity_true_endpoint_and_empty_U_above_gap(self):
        # sqrt bounds and both exact kinematic endpoints reconstructed at
        # higher precision, including a loss smaller than binary64 spacing at E.
        from decimal import Decimal, localcontext
        from xpsloss.model import HBAR2_OVER_2ME_EV_A2

        with localcontext() as ctx:
            ctx.prec = 90
            for energy, t in ((1000.0, 0.2), (3.0, 1.5), (1.0, 1e-20)):
                e, x, h = map(Decimal, (energy, t, HBAR2_OVER_2ME_EV_A2))
                total = (x / (e.sqrt() + (e - x).sqrt())) ** 2 / h
                for flag in (True, False):
                    exact = total if flag else total - x * x / (4 * e * h)
                    bound = _lower_q_squared(energy, t, flag)
                    self.assertLessEqual(bound, Fraction(exact))
                    self.assertAlmostEqual(float(bound) / float(exact), 1.0)
        self.assertEqual(_uniform_bounds(UE), [])
        self.assertTrue(_monotonicity(UE, _lower_q_squared(3.0, 2.5, False))[0]["passed"])
        r = fixed.effective_cross_section_fixed_depth(
            UE, 3.0, 2.5, 10.0, epsilon_momentum="k_parallel", **SETTINGS
        )
        self.assertEqual(r.diagnostics["subgap_proofs"], [])
        candidates = _uniform_bounds(SI)
        self.assertEqual(len(candidates), 8)
        self.assertTrue(all(SI.gap_eV < u < 10.0 for u, _, _ in candidates))
        proof = _proof(candidates, SI.gap_eV, 0.2)
        self.assertTrue(proof["proved"])
        self.assertGreater(proof["margin_lower"], 0.0)

    def test_above_gap_surface_search_refines_internal_crossings(self):
        from xpsloss._gap_response import resolve_surface, _surface_intervals

        # Both endpoints lie on the same side. The two interior crossings must
        # be resolved without relying on opposite endpoint signs.
        edges = np.log([0.1, 2.0])
        refined = resolve_surface(N1, 4.0, edges, 0.0, 200)
        self.assertGreater(len(refined), 3)
        self.assertFalse(np.any(_surface_intervals(N1, 4.0, refined[:-1], refined[1:], 0.0)))
        with self.assertRaisesRegex(fixed.NumericalConvergenceError, "surface resonance"):
            resolve_surface(N1, 4.0, edges, 0.0, 2)


class ManifestTests(unittest.TestCase):
    def test_required_cli_grid_and_budget_manifest(self):
        for mode in MODES:
            with tempfile.TemporaryDirectory() as folder:
                args = cli_args(folder, SI, GRID, mode)
                proc = subprocess.run([sys.executable, "-m", "xpsloss", *args], capture_output=True, cwd=ROOT)
                self.assertEqual(proc.returncode, 0, proc.stderr)
                self.assertEqual(json.loads(proc.stdout.decode("ascii"))["status"], "complete")
                out = Path(folder) / "output"
                doc = json.loads((out / "run.json").read_text(encoding="utf8"))
                self.assertEqual(doc["status"], "complete")
                self.assertTrue(all(w.isascii() for w in doc["warnings"]))
                self.assertGreaterEqual(len(doc["warnings"]), 2)
                diag = doc["diagnostics"]
                self.assertEqual(diag["broadening_eV"], 0)
                self.assertEqual(diag["real_part_method"], "closed_form_D3_retained_band")
                self.assertEqual(len(diag["subgap_proofs"]), 5)
                self.assertTrue(all(p["proved"] for p in diag["subgap_proofs"]))
                self.assertEqual(set(diag["component_errors"]), set(NAMES))
                self.assertEqual(set(diag["error_budget_ratios"]), set(NAMES))
                data = np.genfromtxt(out / "xps_fixed_depth.csv", delimiter=",", names=True)
                errors = np.array([diag["component_errors"][name] for name in NAMES])
                np.testing.assert_array_equal(data["quadrature_error"], errors.max(axis=0))
                for i, name in enumerate(NAMES):
                    ratios = errors[i] / np.maximum(1e-12, abs(data[name]) * 1e-9)
                    np.testing.assert_allclose(ratios, diag["error_budget_ratios"][name], rtol=1e-14, atol=0)
                    self.assertTrue(np.all(ratios <= 1))
                self.assertIn("gapped_checks", doc["applicability"])
                self.assertNotIn("checked_envelope", doc["applicability"])
                self.assertIs(doc["applicability"]["within_checked_envelope"], True)
                self.assertEqual(diag["ill_conditioned_D3_count"], 0)
                self.assertFalse(any("Outside the checked envelope" in w for w in doc["warnings"]))
        r = run_fixed(SI, [SI.gap_eV - 1e-6, SI.gap_eV + 1e-6], "k_total")
        self.assertEqual(len(r.diagnostics["near_gap_edge_points"]), 2)
        self.assertEqual(len(r.warnings), 4)

    def test_hard_gap_envelope_departures_in_cli_manifest(self):
        # A matching name alone does not establish the Si oscillator input.
        altered = replace(SI, oscillators=(replace(SI.oscillators[0], strength_eV2=0.1), *SI.oscillators[1:]))
        for mat, energy, depth, losses, departures in (
            (SI, 1000, 2, [0.2, 99.8], []),
            (SI, 1000, 200, [5.0], []),
            (SI, 120, 10, [0.5, 60.5], ["electron_energy_eV", "max_loss_over_E0"]),
            (SI, 4000, 1.9, [0.19], ["electron_energy_eV", "depth_A", "loss_eV"]),
            (SI, 1000, 201, [99.9], ["depth_A", "loss_eV", "max_loss_over_E0"]),
            (altered, 1000, 10, [5.0], ["material"]),
        ):
            with self.subTest(energy=energy, depth=depth, losses=losses), tempfile.TemporaryDirectory() as folder:
                args = cli_args(folder, mat, losses, energy=energy)
                args[args.index("--depth-A") + 1] = str(depth)
                proc = subprocess.run([sys.executable, "-m", "xpsloss", *args], capture_output=True, cwd=ROOT)
                self.assertEqual(proc.returncode, 0, proc.stderr)
                doc = json.loads((Path(folder) / "output/run.json").read_text(encoding="utf8"))
                app = doc["applicability"]
                self.assertIs(app["within_checked_envelope"], not departures)
                self.assertEqual(app["gapped_checks"]["theta_deg_from_normal"], 0.0)
                warnings = [w for w in doc["warnings"] if "Outside the checked envelope" in w]
                self.assertEqual(len(warnings), len(departures))
                for field in departures:
                    self.assertTrue(any(field in w for w in warnings), field)
                self.assertEqual(json.loads(proc.stdout)["warning_count"], len(doc["warnings"]))
