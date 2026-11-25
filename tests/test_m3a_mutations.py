"""Controlled mutation witnesses from M3A_PREREG V12 (a), (a'), (c)-(h)."""

import math
import unittest
from unittest.mock import patch

import numpy as np

from xpsloss import _gap_response as gap
from xpsloss import xps_depth_averaged as avg
from test_m3a_physics import average, material


class ErrorPropagationTests(unittest.TestCase):
    def test_cartesian_product_cross_term(self):
        self.assertGreaterEqual(avg._project(0j, (2., 0.), (0., 0., 3., 0.)), 6.)

    def test_separate_13b_projections(self):
        k = np.ones((1, 1))
        weights = avg.depth_averaged_weights(k, 1., 1.)
        z = np.zeros_like(k, dtype=complex)
        for sr, si, expected in ((3., 0., weights["N"][0]*3), (0., 5., 10*weights["W13"][1])):
            coeff = (z, z, z), ((sr, si), (0., 0.), (0., 0.))
            bound = avg._term_errors(coeff, weights, k, k, 1., 0.)[:, 1]
            np.testing.assert_allclose(bound, expected, rtol=1e-13)

    def test_multiple_panel_aggregation_roundoff(self):
        values = np.zeros((3, 11))
        values[:, 0] = [2.**52, 1., -2.**52]
        details = avg._PanelDetails(direct_values=values[:, 0].copy())
        _, errors, direct, _ = avg._aggregate(values, np.zeros_like(values), details, 1.)
        floor = avg.EPS*sum(abs(values[:, 0]))
        self.assertGreaterEqual(errors[0], floor)
        self.assertGreaterEqual(direct["e_dir"], floor)

    def test_above_gap_three_coordinate_changes(self):
        # Isolate the changes in E4b, in the removed band and in the imaginary
        # E1 at one node. The optical and D3 roundoff terms cancel in differences.
        mat = material("N2")
        osc, = mat.oscillators
        t, q = 4., np.array([.8])
        shift = osc.dispersion*avg.HBAR2_OVER_2ME_EV_A2*q[0]**2
        center = osc.energy_eV+shift
        a, width, g = osc.strength_eV2, osc.width_eV, mat.gap_eV
        den = ((center-t)*(center+t))**2+(width*t)**2
        def changes(rz):
            dw = (rz+8*avg.EPS)*shift + avg.EPS*center
            dx = dw*(2*center+dw)*(1+8*avg.EPS)
            dist = max(width*t, math.sqrt(den)*(1-8*avg.EPS)-dx)
            m = math.nextafter((center-g)*(center+g)*(1-8*avg.EPS)-dx, -math.inf)
            return np.array([a*dx/dist**2, (4/math.pi)*a*width*dx/m**3, 2*a*width*t*dx/dist**3])
        rz = 1e-4
        bare, removed, imag = changes(rz) - changes(0.)
        def variation(moment):
            with patch.object(gap, "_removed_moment", return_value=moment):
                _, sr0, si0 = gap.gapped_sum(mat, t, q, 0.)
                _, sr, si = gap.gapped_sum(mat, t, q, rz)
                return (sr-sr0).item(), (si-si0).item()
        real0, imag0 = variation(0.)
        real1, _ = variation(1.)
        for actual, expected in ((real0, bare), (real1-real0, removed), (imag0, imag)):
            self.assertGreater(expected, 0)
            self.assertTrue(math.isclose(actual, expected, rel_tol=1e-12))

    def test_above_gap_node_errors_reach_the_averaged_terms(self):
        point = avg._Point(material("N2"), 1000., 4., 0, .04, .04, .05, 1., False, True)
        k = np.array([[.8]])
        rel = np.array([[1e-4]])
        args = (point, k, rel, np.array([-.3]), np.array([-.2]), None)
        real_coeff = avg.fixed._coefficients
        with patch.object(avg.fixed, "_coefficients", wraps=real_coeff) as called:
            _, base, _ = avg._densities(*args)
        np.testing.assert_array_equal(called.call_args.kwargs["relative_z"], rel)
        def perturb(*a, **kw):
            coeff = real_coeff(*a, **kw)
            coeff[1][0][1][:] += .01
            return coeff
        with patch.object(avg.fixed, "_coefficients", side_effect=perturb):
            _, new, _ = avg._densities(*args)
        measure = k*k/(math.pi*avg.BOHR_RADIUS_A*1000.)
        expected = .01*measure/(k*k+point.b**2)
        np.testing.assert_allclose(new[:, 0]-base[:, 0], expected, rtol=1e-13)


class DirectSumTests(unittest.TestCase):
    def test_reconstruction_difference_exhausts_budget(self):
        # The positive floor and e_dir alone fit. Only adding delta_rec fails.
        point = avg._Point(material("Si"), 1000., 1., 0, .01, .01, .05, 1., True, True)
        mesh = avg._Mesh(point, np.array([0., 1.]), 1)
        def panels(left, right, *, point, details, allowance):
            values = np.zeros((1, 11))
            values[0, [0, 8, 10]] = 2.
            errors = np.zeros_like(values)
            errors[0, -1] = .25
            details.floors = np.zeros_like(values)
            details.direct_values = np.array([0.])
            details.direct_floor = np.array([0.])
            return values, errors
        with patch.object(avg, "_panel_integrals", side_effect=panels):
            mesh.evaluate(1., 1e-9)
            mesh.check_floor()
            with self.assertRaises(avg.DepthAverageConvergenceError) as caught:
                mesh.adapt(1., 1e-9, 1)
        self.assertEqual(caught.exception.column, "total")
        self.assertGreater(caught.exception.estimate, caught.exception.budget)
        self.assertLess(mesh.direct["e_dir"], 1.)

    def test_independent_node_sum_in_mixed_initial_and_subdivided_panels(self):
        real_panels = avg._panel_integrals
        real_densities = avg._densities
        seen = {0: [], 1: []}
        def controlled(*args):
            terms, errors, floor = real_densities(*args)
            # Direct order gives 1. Group order gives zero after charge addition.
            terms.fill(0.)
            for i, v in enumerate((2.**54, 1., -2.**54, 1.)):
                terms[:, i, :] = v
            return terms, errors*0., floor

        def observe(left, right, *, point, details, allowance):
            values, errors = real_panels(left, right, point=point, details=details, allowance=allowance)
            old = seen[point.index]
            old.append((left.copy(), right.copy(), details.direct_values))
            if len(old) == 1:
                errors[0, 0] += 1e8  # force one subdivision for both losses
            return values, errors

        with patch.object(avg, "_densities", side_effect=controlled), patch.object(avg, "_panel_integrals", side_effect=observe):
            r = average(material("Si"), 1000., [1., 2.], weight_length_A=20., epsilon_momentum="k_total",
                        epsabs=1e6, epsrel=1e-14)
        for i, calls in seen.items():
            self.assertEqual(len(calls), 2)
            self.assertEqual(len(calls[-1][0]), len(calls[0][0])+1)
            for left, right, direct in calls:
                if i == 0:
                    np.testing.assert_allclose(direct, right-left, rtol=1e-14)
                else:
                    self.assertIsNone(direct)
        d = r.diagnostics["direct_sum"][0]
        self.assertAlmostEqual(d["K_dir"], sum(seen[0][-1][1]-seen[0][-1][0]), places=12)
        self.assertGreater(d["delta_rec"], 0.)
        self.assertIsNone(r.diagnostics["direct_sum"][1])

    def test_direct_branch_requires_exact_unit_charge(self):
        for q in (0., .5, math.nextafter(1., 2.)):
            r = average(material("Si"), 1000., [1., 2.], weight_length_A=20.,
                        epsilon_momentum="k_total", hole_charge=q, epsabs=1e-6)
            self.assertEqual(r.diagnostics["direct_sum"], [None, None])
            self.assertIsNone(r.diagnostics["direct_sum_max"])
