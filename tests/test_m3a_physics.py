"""Implementation checks V10 and V11. XPSLOSS_M3A_FULL=1 selects the full preregistered sets."""

from functools import lru_cache
import json
import math
import os
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

from xpsloss import Material, Oscillator, effective_cross_section_depth_averaged as average
from xpsloss import xps_depth_averaged as avg, xps_fixed_depth as fixed
from xpsloss.dielectric import DielectricGapEdgeError


ROOT = Path(__file__).resolve().parents[1]
FULL = os.environ.get("XPSLOSS_M3A_FULL") == "1"
GA = 0.2*np.arange(1, 500)
F = 0.2/4*np.arange(4, 1997)
GRIDS = dict(GA=GA, GC=np.linspace(.2, 99.8, 499), F=F, F2=F[::2],
             FRD=np.r_[F[::-4], F[:40:4]], GA10G=np.r_[GA[:10], 1.12])
DEF = dict(epsabs=1e-12, epsrel=1e-9, limit=200)


@lru_cache(maxsize=1)
def table():
    return json.loads((ROOT/"benchmarks/numerics/m3a/expected_states_m3a.json").read_bytes())


@lru_cache(maxsize=None)
def material(name):
    row = table()["materials"][name]
    return Material(name, row["gap_eV"], tuple(Oscillator(*o) for o in row["oscillators"]))


def losses(row):
    spec = row["losses"]
    return GRIDS[spec["array"]] if "array" in spec else np.array([float.fromhex(t) for t in spec["values_hex"]])


def options(row):
    return dict(weight_length_A=row["weight_length_A"]["value"], epsilon_momentum=row["epsilon_momentum"],
                hole_charge=row["hole_charge"]["value"], theta_deg=row["theta_deg"],
                **table()["settings"][row["settings"]])


def evaluate(row):
    return average(material(row["material"]), row["electron_energy_eV"], losses(row), **options(row))


def columns(r):
    return np.column_stack([getattr(r, c) for c in avg.COLUMNS])


class MeshChecks(unittest.TestCase):
    def test_v10_refine_all_final_panels(self):
        for name in ("Al", "Si"):
            for length in ((5., 20., 80.) if FULL else (20.,)):
                settings = dict(hole_charge=1., weight_length_A=length, epsilon_momentum="k_total", limit=200)
                for t in GA:
                    with self.subTest(material=name, length=length, loss=t):
                        def point(abs_tol, rel_tol, **mesh):
                            return avg._depth_averaged_point(material(name), 1000., t, epsabs=abs_tol,
                                                             epsrel=rel_tol, **settings, **mesh)
                        base, _, _ = point(1e-12, 1e-9)
                        budget = np.maximum(1e-12, 1e-9*abs(base))
                        ea, er = 5e-13, 5e-10
                        a, errors, edges = point(ea, er)
                        if np.any(errors > budget/2):
                            ea, er = ea/2, er/2
                            a, errors, edges = point(ea, er)
                        np.testing.assert_array_less(errors, np.nextafter(budget/2, np.inf))
                        split = np.sort(np.r_[edges, (edges[:-1]+edges[1:])/2])
                        b, errors, final = point(ea, er, panels=split, adapt=False)
                        np.testing.assert_array_equal(final, split)
                        self.assertTrue(np.all(errors <= budget/2))
                        self.assertTrue(np.all(abs(a-b) <= budget))

    def test_v10_grid_order_and_duplicates(self):
        for name in ("Al", "Si"):
            def run(grid):
                return average(material(name), 1000., grid, weight_length_A=20., epsilon_momentum="k_total")
            previous = {}
            def compare(result):
                values = columns(result)
                errors = np.column_stack([result.component_errors[c] for c in avg.COLUMNS])
                for j, t in enumerate(result.loss_eV):
                    occurrences = previous.setdefault(t.hex(), [])
                    for old, old_errors in occurrences:
                        self.assertTrue(np.all(abs(old-values[j]) <= old_errors+errors[j]), (name, t))
                    occurrences.append((values[j], errors[j]))
            compare(run(GA))
            for key in (("F", "F2", "FRD") if FULL else ("FRD",)):
                compare(run(GRIDS[key]))


class GapChecks(unittest.TestCase):
    def test_finite_floor_with_subnormal_budget(self):
        with self.assertRaises(avg.DepthAverageConvergenceError) as caught:
            average(material("Si"), 1000., .2, weight_length_A=1e-5, epsilon_momentum="k_total",
                    epsabs=math.ulp(0.), epsrel=math.ulp(0.))
        exc = caught.exception
        self.assertEqual((exc.reason, exc.column, exc.subdivisions), ("direct_sum_floor_exceeds_budget", "total", 0))
        self.assertTrue(math.isfinite(exc.estimate) and exc.estimate > exc.budget)

    def test_rejections_and_failure_stage_order(self):
        args = material("Al"), 1000., 5.
        with patch.object(avg, "_panel_integrals", wraps=avg._panel_integrals) as calls:
            for length in (0., -1., float("nan"), float("inf"), 1+0j, True):
                with self.assertRaisesRegex(ValueError, "weight_length_A"):
                    average(*args, weight_length_A=length, epsilon_momentum="k_total")
            with self.assertRaisesRegex(ValueError, "theta_deg"):
                average(*args, weight_length_A=20., epsilon_momentum="k_total", theta_deg=10.)
            with self.assertRaisesRegex(ValueError, r"loss_eV=5, index=\(0,\).*integration limit too small"):
                average(*args, weight_length_A=20., epsilon_momentum="k_total", limit=1)
            with self.assertRaises(avg.DepthAverageConvergenceError) as failure:
                average(material("SP"), 1000., [2.5, math.ulp(0.)], weight_length_A=20., epsilon_momentum="k_total")
            self.assertEqual((failure.exception.reason, failure.exception.index), ("momentum_bounds_unresolved", 1))
            calls.assert_not_called()
            # The later point's node failure precedes the earlier point's floor.
            with self.assertRaises(avg.DepthAverageConvergenceError) as failure:
                average(material("N1"), 1000., [1., 4.], weight_length_A=20., epsilon_momentum="k_total",
                        epsabs=1e-30, epsrel=1e-15)
            self.assertEqual((failure.exception.reason, failure.exception.index), ("d3_ill_conditioned", 1))
            self.assertIn("D3", failure.exception.inner)
            self.assertTrue(calls.called)

    def test_v11_preflight_matches_fixed_depth_before_panels(self):
        class Accepted(Exception):
            pass

        def stop_at_proof(*args, **kwargs):
            raise Accepted(*fixed_preflight(*args, **kwargs))

        fixed_preflight = fixed.gap_preflight
        seen = set()
        for row in table()["runs"]:
            mat = material(row["material"])
            if not mat.gap_eV or row["family"] == "RJ5":
                continue
            grid, mode, energy = losses(row), row["epsilon_momentum"], row["electron_energy_eV"]
            key = (mat, grid.tobytes(), mode, energy)
            if key in seen:
                continue
            seen.add(key)
            args = (mat, energy, grid)
            with self.subTest(run=row["id"]), patch.object(avg, "_panel_integrals") as panels:
                outcomes = []
                for module, function, kwargs in (
                    (fixed, fixed.effective_cross_section_fixed_depth, {"depth_A": 10.}),
                    (avg, average, {"weight_length_A": 20.}),
                ):
                    with patch.object(module, "gap_preflight", side_effect=stop_at_proof):
                        with self.assertRaises((Accepted, ValueError)) as caught:
                            function(*args, **kwargs, epsilon_momentum=mode)
                    exc = caught.exception
                    outcomes.append((type(exc), exc.args, vars(exc)))
                self.assertEqual(outcomes[0], outcomes[1])
                if row["category"] == "RJ":
                    self.assertEqual(type(exc).__name__, row["expect"]["exception"])
                    if row["expect"].get("reason"):
                        self.assertEqual(exc.reason, row["expect"]["reason"])
                    if isinstance(exc, fixed.FixedDepthGapDomainError):
                        self.assertEqual([p["index"] for p in exc.details["points"]], row["expect"]["offending_indices"])
                    if isinstance(exc, DielectricGapEdgeError):
                        self.assertEqual(exc.offending_indices, row["expect"]["offending_indices"])
                else:
                    self.assertIsInstance(exc, Accepted)
                panels.assert_not_called()

    def test_v11_zero_strength_controls_and_edge_values(self):
        for row in table()["runs"]:
            if row["family"] == "MS12":
                r = evaluate(row)
                if row["material"] == "Z0":
                    self.assertEqual(columns(r).tobytes(), np.zeros_like(columns(r)).tobytes())
                    self.assertFalse(any("gap" in w.lower() for w in r.warnings))
                else:
                    self.assertTrue(all(p["proved"] for p in r.diagnostics["subgap_proofs"]))
                    si = average(material("Si"), 1000., r.loss_eV, **options(row))
                    self.assertEqual(r.diagnostics["subgap_proofs"], si.diagnostics["subgap_proofs"])
        si = material("Si")
        with self.assertRaises(DielectricGapEdgeError):
            average(si, 1000., si.gap_eV, weight_length_A=20., epsilon_momentum="k_total")
        edge = average(si, 1000., [1., si.gap_eV+1e-6], weight_length_A=20., epsilon_momentum="k_total")
        self.assertLessEqual(abs(edge.total[0]), edge.component_errors["total"][0])
        self.assertGreater(edge.total[1], 0.)
        self.assertEqual(edge.diagnostics["near_gap_edge_points"][0]["loss_eV"], edge.loss_eV[1])

    def test_v11_direct_sum_plan(self):
        runs = {r["id"]: r for r in table()["runs"]}
        groups = {}
        for key, index, _ in table()["check_plan"]["V11iii"]:
            if FULL or runs[key]["routine"]:
                groups.setdefault(key, []).append(index)
        self.assertEqual(sum(map(len, groups.values())), 249 if FULL else 101)
        for key, indices in groups.items():
            row = runs[key]
            with self.subTest(run=key):
                try:
                    if row["entry"] == "CLI":
                        from test_m3a_cli import run_cli_capture
                        r, _ = run_cli_capture(row)
                    else:
                        r = evaluate(row)
                except avg.DepthAverageConvergenceError as exc:
                    self.assertEqual(row["category"], "EI")
                    self.assertIn(exc.reason, row["expect"]["allowed_reasons"])
                    continue
                direct = r.diagnostics["direct_sum"]
                for i, record in enumerate(direct):
                    self.assertEqual(record is not None, i in indices)
                    if record is None:
                        continue
                    self.assertLessEqual(abs(r.total[i]), min(r.component_errors["total"][i],
                                                             options(row)["epsabs"]))
                    self.assertEqual(record["delta_rec"], abs(r.total[i]-record["K_dir"]))
                    self.assertEqual(record["e_total"], r.component_errors["total"][i])
                    self.assertAlmostEqual(record["e_total"], record["e_dir"] + record["delta_rec"],
                                           delta=2*np.finfo(float).eps*record["e_total"])
                self.assertEqual(r.diagnostics["direct_sum_max"],
                                 {k: max(direct[i][k] for i in indices) for k in ("e_dir", "delta_rec")})
