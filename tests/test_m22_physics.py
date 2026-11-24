"""Implementation-side checks specified before M2-2 (docs/M2_2_PREREG.md)."""

from dataclasses import replace
from functools import lru_cache
import math
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np
from scipy.integrate import quad

from xpsloss import Material, Oscillator, bulk_diimfp, effective_cross_section_fixed_depth
from xpsloss import xps_fixed_depth as fixed
from xpsloss._gap_response import _closed_band, gapped_sum
from xpsloss.config import load_config
from xpsloss.dielectric import _closed_kk_band
from xpsloss.model import BOHR_RADIUS_A, HBAR2_OVER_2ME_EV_A2


ROOT = Path(__file__).resolve().parents[1]
SI = load_config(ROOT / "benchmarks/inputs/si_elfdb_v1.0.bulk.json").material
S0 = Material("synthetic S0", SI.gap_eV, tuple(replace(o, dispersion=0.0) for o in SI.oscillators))
NARROW = Material("synthetic resolved surface line", 2.0, (Oscillator(3, 6, 0.025, 0.5),))
MODES = ("k_parallel", "k_total")
GRID = 0.2 * np.arange(1, 500)
REF_POINTS = [
    0.2,
    0.4,
    1.0,
    1.119,
    SI.gap_eV - 1e-6,
    SI.gap_eV + 1e-6,
    SI.gap_eV + 0.001,
    1.2,
    2.0,
    5.0,
    11.0,
    16.8,
    17.0,
    30.0,
    60.0,
    99.8,
]
NAMES = (*fixed.TERM_NAMES, "no_hole", "hole", "total")
SETTINGS = dict(epsabs=1e-12, epsrel=1e-9, limit=200)


def run_fixed(material, t, mode, depth=10.0, charge=1.0, **settings):
    return effective_cross_section_fixed_depth(
        material, 1000.0, t, depth, hole_charge=charge, epsilon_momentum=mode, **(SETTINGS | settings)
    )


def values(r):
    return np.array([getattr(r, name) for name in NAMES])


def budget(table):
    return np.maximum(1e-12, 1e-9 * abs(table))


@lru_cache(maxsize=2)
def si_10(mode):
    return run_fixed(SI, GRID, mode)


def e1(material, t, q):
    return sum(
        o.strength_eV2
        * o.width_eV
        * t
        / (
            ((o.energy_eV + o.dispersion * HBAR2_OVER_2ME_EV_A2 * q * q) ** 2 - t * t) ** 2
            + (o.width_eV * t) ** 2
        )
        for o in material.oscillators
    )


def retained(material, t, q, *, uncertainty=False):
    """Independent E4a. Subtract u*L(u) at t before quadrature; cancel u^2-t^2
    algebraically. Integrate in u/center so even weak, high-q bands are resolved.
    Returns S itself, preserving weak lines that disappear when added to -1.
    """
    g = material.gap_eV
    total, errors = [], []
    for o in material.oscillators:
        w = o.energy_eV + o.dispersion * HBAR2_OVER_2ME_EV_A2 * q**2
        width = o.width_eV

        def regular(v):
            u = w * v
            den = ((w - u) * (w + u)) ** 2 + (width * u) ** 2
            return w * (w**4 - (u * t) ** 2) / den

        # Integrate like-signed pieces on both sides of the numerator zero.
        edges = [g / w, *sorted(x for x in {1.0, w / t} if x > g / w), np.inf]
        integral, error = 0.0, 0.0
        for a, b in zip(edges, edges[1:]):
            value, estimate, *info = quad(regular, a, b, epsabs=1e-15, epsrel=1e-13, full_output=1)
            if len(info) != 1:
                raise AssertionError(info[-1])
            integral += value
            error += estimate
        pv = t / 2 * math.log1p(2 * min(t, g) / abs(g - t))
        dt = ((w - t) * (w + t)) ** 2 + (width * t) ** 2
        scale = (2 / math.pi) * o.strength_eV2 * width / dt
        total.append(scale * (integral + pv))
        errors.append(scale * (error + 8 * np.finfo(float).eps * (abs(integral) + pv)))
    return (math.fsum(total), math.fsum(errors)) if uncertainty else math.fsum(total)


class NodeResponseTests(unittest.TestCase):
    def test_d3_counts_only_ill_conditioned_nodes(self):
        # The dispersionless resonance fails twice, the off-resonance center passes.
        _closed_band(4.0, 2.0, np.array([8.0]), 4.0, 1e-6)
        with self.assertRaisesRegex(fixed.NumericalConvergenceError, "D3 is ill-conditioned") as caught:
            _closed_band(4.0, 2.0, np.array([4.0, 8.0, 4.0]), 4.0, 1e-6)
        self.assertEqual(caught.exception.count, 2)

    def test_v5_all_d3_branches_against_scalar_and_retained_pv(self):
        # Same array on either side of the divided-difference cut and critical damping.
        width = 4.0
        switch = math.sqrt((width / 2) ** 2 - (width / 8) ** 2)
        upper_switch = math.sqrt((width / 2) ** 2 + (width / 8) ** 2)
        centers = np.array(
            [
                0.2,
                1.5,
                switch - 1e-10,
                switch + 1e-10,
                2 - 1e-9,
                2,
                2 + 1e-9,
                upper_switch - 1e-10,
                upper_switch + 1e-10,
                3,
                9,
                30.0,
            ]
        )
        for g in (1.0, 3.0 + 1e-12, 15.0):
            for t in (0.6, g - 1e-6, g + 1e-6, 5.0):
                raw, error, _ = _closed_band(t, g, centers, 1.0, width)
                scalar = [_closed_kk_band(t, g, w, 1.0, width)[0] for w in centers]
                np.testing.assert_allclose(raw, scalar, rtol=0, atol=1e-14)
                for i, w in enumerate(centers):
                    ref, u = retained(
                        Material("synthetic", g, (Oscillator(w, 1.0, width),)), t, 0.0, uncertainty=True
                    )
                    self.assertLessEqual(abs(raw[i] - ref) + u, 1e-12, (t, g, w))
                    self.assertLessEqual(abs(raw[i] - ref), error[i], (t, g, w))

    def test_v5_si_mixed_nodes_weak_lines_and_regular_epsilon_pole(self):
        # Includes points on both sides of the F=0 pole in T, and cancellation
        # between the optical lines above the first resonance.
        t = np.array([0.2, 1.12 - 1e-6, 1.12 + 1e-6, 1.2, 11.0, 15.0, 17.0, 60.0, 5.0])
        q = np.array([0.0, 0.0, 0.1, 1.0, 0.3, 0.2, 0.5, 1.0, 4.0])
        for mat in (
            SI,
            replace(
                SI,
                name="synthetic weak",
                oscillators=tuple(replace(o, strength_eV2=1e-20 * o.strength_eV2) for o in SI.oscillators),
            ),
        ):
            s, er, ei = gapped_sum(mat, t, q, 0.0)
            ref, u = np.array([retained(mat, x, y, uncertainty=True) for x, y in zip(t, q)]).T
            self.assertTrue(np.all(abs(ref - s.real) + u <= 1e-12))
            im = np.array([e1(mat, x, y) if x > mat.gap_eV else 0.0 for x, y in zip(t, q)])
            scalar = [
                math.fsum(
                    _closed_kk_band(
                        x,
                        mat.gap_eV,
                        o.energy_eV + o.dispersion * HBAR2_OVER_2ME_EV_A2 * y * y,
                        o.strength_eV2,
                        o.width_eV,
                    )[0]
                    for o in mat.oscillators
                )
                for x, y in zip(t, q)
            ]
            np.testing.assert_allclose(s.real, scalar, rtol=0, atol=1e-14)
            self.assertTrue(np.all(abs(ref - s.real) <= np.minimum(1e-12, er)))
            self.assertTrue(np.all(abs(im - s.imag) <= np.minimum(1e-12, ei)))
            if mat is SI:
                self.assertLess(s.real[0], 1.0)
                self.assertGreater(s.real[1], 1.0)
            else:
                self.assertTrue(np.all(s.real - 1 == -1.0))
                self.assertTrue(np.all(s.real != 0))
        zero = np.zeros(1)
        with patch.object(fixed, "gapped_sum", return_value=(np.ones(1, complex), zero.copy(), zero.copy())):
            (s, r, c), _ = fixed._coefficients(SI, 0.2, zero)
        np.testing.assert_array_equal([s, r, c], [[1.0], [1.0], [1.0]])

    def test_coordinate_error_across_finite_intervals(self):
        q = np.array([0.3, 0.7, 1.2, 2.0, 6.0])
        for t in (1.0, SI.gap_eV + 1e-6, 5.0, 17.0, 60.0):
            relative_z = 1e-7
            center, sr, si = gapped_sum(SI, t, q, relative_z)
            for sign in (-1.0, 1.0):
                moved, _, _ = gapped_sum(SI, t, q * np.sqrt(1 + sign * relative_z), 0.0)
                self.assertTrue(np.all(abs(moved.real - center.real) <= sr))
                self.assertTrue(np.all(abs(moved.imag - center.imag) <= si))

    def test_v7_edge_on_both_sides(self):
        g = SI.gap_eV
        distances = np.array([1e-6, 1e-8, 1e-10])
        for q in (0.0, 0.1, 1.0):
            for sign in (-1, 1):
                t = g + sign * distances
                s, _, _ = gapped_sum(SI, t, q, 0.0)
                slope = np.diff(s.real) / np.diff(np.log(g / abs(g - t)))
                np.testing.assert_allclose(slope, e1(SI, g, q) / math.pi, rtol=1e-3, atol=0)


class FixedDepthHardGapTests(unittest.TestCase):
    def test_narrow_line_adaptive_refinement_and_surface_search(self):
        # The dense comparison starts on the oscillator mesh alone, so a
        # mistaken surface breakpoint is not shared by the two meshes.
        def dense(material, loss, edges, offset, limit):
            return np.r_[np.linspace(edges[:-1], edges[1:], 32, endpoint=False).T.ravel(), edges[-1]]

        for mode in MODES:
            with self.subTest(mode=mode):
                with patch.object(fixed, "_panel_integrals", wraps=fixed._panel_integrals) as panels:
                    searched = run_fixed(NARROW, 4.0, mode)
                if mode == "k_parallel":
                    # Single loss: the first call is the initial mesh, subsequent
                    # calls integrate the two children of an adaptive split.
                    self.assertGreater(panels.call_count, 1)
                    for call in panels.call_args_list[1:]:
                        self.assertEqual(len(call.args[0]), 2)
                with patch.object(fixed, "resolve_surface", side_effect=lambda m, t, edges, off, lim: edges):
                    unsearched = run_fixed(NARROW, 4.0, mode)
                with patch.object(fixed, "resolve_surface", side_effect=dense):
                    reference = values(run_fixed(NARROW, 4.0, mode, limit=4000))
                for result in (searched, unsearched):
                    diff = abs(values(result) - reference)
                    self.assertTrue(np.all(diff <= budget(reference)))
                    self.assertTrue(np.all(diff <= [result.component_errors[name] for name in NAMES]))

    def test_surface_search_resolves_a_crossing_created_by_the_removed_band(self):
        # On this descending flank the bare Lorentzian is already below 2.
        # E4a independently locates a crossing in the corrected Re S, and
        # the complex variation exceeds half the minimum Im S on the panel.
        t = 4.0
        q = np.array([0.844475, 0.849475])
        w = 3.0 + 0.5 * HBAR2_OVER_2ME_EV_A2 * q**2
        detuning = w**2 - t**2
        damping = 0.025 * t
        self.assertTrue(np.all(detuning > damping))
        self.assertTrue(np.all(6.0 * detuning / (detuning**2 + damping**2) < 2))
        s = np.array([retained(NARROW, t, qi) for qi in q]) + 1j * e1(NARROW, t, q)
        self.assertGreater(s[0].real, 2)
        self.assertLess(s[1].real, 2)
        self.assertGreater(abs(s[1] - s[0]), 0.5 * min(s.imag))
        edges = np.log(q)
        self.assertGreater(len(fixed.resolve_surface(NARROW, t, edges, 0.0, 200)), len(edges))

    def test_direct_total_error_includes_group_discrepancy(self):
        # One exact synthetic panel. The returned eight terms cancel; only
        # the separately integrated ninth column carries an offset.
        offset, uncertainty = 2.0**-42, 2.0**-44
        v = np.zeros((1, 9))
        v[0, [0, 2, 8]] = (1.0, -1.0, offset)
        err = np.zeros_like(v)
        err[0, 8] = uncertainty
        mesh = (np.array([0.0, 1.0]), (v, err, np.zeros_like(v)))
        terms, estimates = fixed._integrate_gap(
            SI, 0.2, 1000.0, 10.0, 1.0, 0.0, 1.0, math.e, 1.0, 1e-12, 1e-9, 1, mesh
        )
        self.assertEqual(terms[0] + terms[2], 0.0)
        self.assertGreaterEqual(estimates[-1], offset + uncertainty)
        self.assertLessEqual(estimates[-1], 1e-12)

    def test_near_unit_charge_subgap_arithmetic_floor(self):
        for mode in MODES:
            self.assertLessEqual(abs(run_fixed(SI, 0.2, mode).total[0]), SETTINGS["epsabs"])
            for charge in (0.999999, 1.0 + 1e-12):
                with self.subTest(mode=mode, charge=charge), self.assertRaisesRegex(
                    fixed.FixedDepthConvergenceError, "arithmetic error exceeds an output budget"
                ):
                    run_fixed(SI, 0.2, mode, charge=charge)

    def test_v1_bulk_and_dispersionless_analytic(self):
        for mat in (S0, SI):
            result = si_10("k_total") if mat is SI else run_fixed(mat, GRID, "k_total")
            bulk = np.array([bulk_diimfp(mat, 1000.0, t, **SETTINGS)[0] for t in GRID])
            np.testing.assert_allclose(result.term_13a, bulk, rtol=1e-8, atol=0)
            self.assertTrue(np.all(bulk[GRID < mat.gap_eV] == 0))
            if mat is S0:
                ratio = (math.sqrt(1000) + np.sqrt(1000 - GRID)) ** 2 / GRID
                analytic = np.where(
                    GRID > mat.gap_eV,
                    e1(mat, GRID, 0.0) * np.log(ratio) / (math.pi * BOHR_RADIUS_A * 1000),
                    0.0,
                )
                np.testing.assert_allclose(result.term_13a, analytic, rtol=1e-10, atol=0)

    def test_v2_deep_residuals(self):
        for mode in MODES:
            residuals = []
            for a in (1e4, 1e5):
                r = run_fixed(SI, GRID[4:], mode, a)
                residuals.append(
                    a
                    * np.array(
                        [r.no_hole - r.term_13a, r.term_13b + r.term_15_T1 + r.term_15_T3 + r.term_16_e, r.hole]
                    )
                )
            first, last = residuals
            self.assertTrue(np.all(abs(first - last) <= 1e-6 * np.maximum(abs(last), 1e-6)), mode)

    def test_v4_v10_required_grid_at_all_depths(self):
        for mode in MODES:
            for depth in (2.0, 10.0, 50.0, 200.0):
                r = si_10(mode) if depth == 10 else run_fixed(SI, GRID, mode, depth)
                output = values(r)
                errors = np.array([r.component_errors[j] for j in NAMES])
                self.assertTrue(np.all((errors >= 0) & (errors <= budget(output))))
                np.testing.assert_array_equal(r.quadrature_error, errors.max(axis=0))
                sub = GRID < SI.gap_eV
                mag = abs(output[:8, sub]).sum(axis=0)
                self.assertTrue(np.all(abs(r.total[sub]) <= np.maximum(1e-12, 1e-13 * mag)))
                self.assertTrue(np.all(r.total[~sub] > 0))
                np.testing.assert_allclose(
                    r.diagnostics["max_subgap_identity_residual"],
                    np.max(abs(r.no_hole[sub] + r.hole[sub]) / mag),
                    rtol=1e-14,
                    atol=0,
                )

    def test_v6_factorized_dispersionless_terms(self):
        t = 1.0
        s_half = replace(
            S0,
            name="synthetic half S0",
            oscillators=tuple(replace(o, strength_eV2=o.strength_eV2 / 2) for o in S0.oscillators),
        )
        scaled = []
        for mat in (S0, s_half):
            s = retained(mat, t, 0.0)
            refl = s / (2 - s)
            c = 2 * s - refl
            coeff = np.array([-s, -s, c, c, refl, -refl, -refl])
            self.assertTrue(np.all(abs(coeff) > 0.1))
            scaled.append(values(run_fixed(mat, [t], "k_parallel"))[1:8, 0] / coeff)
        np.testing.assert_allclose(*scaled, rtol=1e-9, atol=0)

    def test_v8_cubic_gap_scaling(self):
        t = np.array([5.0, 17.0])
        bare = values(run_fixed(replace(S0, gap_eV=0.0), t, "k_parallel"))
        corrections = []
        for g in (0.4, 0.2, 0.1):
            val = values(run_fixed(replace(S0, gap_eV=g), t, "k_parallel"))
            diff = val - bare
            self.assertTrue(np.all(abs(diff[0]) <= budget(bare)[0]))
            if corrections:
                mask = abs(diff) >= 100 * budget(bare)
                self.assertTrue(np.any(mask))
                np.testing.assert_allclose(corrections[-1][mask] / diff[mask], 8, rtol=0.05, atol=0)
            corrections.append(diff)
        d4, d2, d1 = corrections
        np.testing.assert_allclose((d4 - d2)[mask] / (d2 - d1)[mask], 8, rtol=0.05, atol=0)

    def test_v9_full_panel_bisection_and_loss_grids(self):
        initial = fixed._initial_edges

        def bisect(*args):
            edges = initial(*args)
            return np.sort(np.r_[edges, (edges[1:] + edges[:-1]) / 2])

        for mode in MODES:
            baseline = values(si_10(mode))
            with patch.object(fixed, "_initial_edges", side_effect=bisect):
                refined = values(run_fixed(SI, GRID, mode))
            self.assertTrue(np.all(abs(baseline - refined) <= budget(baseline)))
            previous = baseline
            for n in (2, 4):
                t = 0.2 / n * np.arange(n, 499 * n + 1)
                finer = values(run_fixed(SI, t, mode))
                common = finer[:, ::2]
                self.assertTrue(np.all(abs(common - previous) <= budget(common) + budget(previous)))
                self.assertTrue(
                    np.all(abs(finer[:, ::n] - baseline) <= budget(finer[:, ::n]) + budget(baseline))
                )
                previous = finer

    def test_v11_charge_linearity_with_frozen_unit_groups(self):
        for mode in MODES:
            groups = values(si_10(mode))
            indices = np.array([0, 3, 4, 5, 24, 83, 149, 498])
            for charge in (-0.5, 0.0, 0.5, 1.0, 1.5):
                result = run_fixed(SI, GRID[indices], mode, charge=charge)
                output = values(result)
                expected = groups[:, indices].copy()
                expected[-1] = expected[-3] + charge * expected[-2]
                self.assertTrue(np.all(abs(output - expected) <= budget(expected)))
                if charge in (-0.5, 1.5):
                    errors = result.component_errors
                    required = errors["no_hole"] + abs(charge) * errors["hole"]
                    # Allow for the different order of positive panel reductions.
                    self.assertTrue(np.all(errors["total"] >= (1 - 2e-14) * required))

    def test_preregistered_dispersionless_points_and_v16_zero_strength(self):
        zero = replace(
            SI, name="synthetic Z0", oscillators=tuple(replace(o, strength_eV2=0.0) for o in SI.oscillators)
        )
        extra = replace(SI, name="synthetic Z1", oscillators=SI.oscillators + (Oscillator(0.5, 0.0, 0.3, 0.5),))
        t = np.array(REF_POINTS)
        for mode in MODES:
            self.assertTrue(np.all(np.isfinite(values(run_fixed(S0, t, mode)))))
            r = run_fixed(zero, np.r_[0.2 * np.arange(1, 11), SI.gap_eV], mode)
            self.assertTrue(np.all(values(r) == 0.0))
            self.assertFalse(r.warnings)
            self.assertIsNone(r.diagnostics["max_subgap_identity_residual"])
            before = run_fixed(SI, t, mode)
            after = run_fixed(extra, t, mode)
            self.assertEqual(before.diagnostics["subgap_proofs"], after.diagnostics["subgap_proofs"])
            self.assertTrue(np.all(abs(values(after) - values(before)) <= 1e-3 * budget(values(before))))
