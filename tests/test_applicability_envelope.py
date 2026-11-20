"""Independent I3 comparison at both fixed-depth energy bounds on unscaled materials."""
from pathlib import Path
import sys
import unittest

import numpy as np

from xpsloss.xps_fixed_depth import effective_cross_section_fixed_depth

# Resolve the frozen I3 fixtures for both discovery and a named unittest run.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_xps_route_b import ABS, HOLE_CHARGES, MATERIALS, REL  # noqa: E402
from xps_route_b_reference import TERMS, reference_terms  # noqa: E402


ENDPOINTS = ((1000.0, (5.0, 15.6, 50.0, 100.0)),
             (4000.0, (5.0, 11.0, 15.6, 30.0, 100.0, 400.0)))


def endpoint_error_ratios(energy, losses):
    """Yield pointwise |package - reference| / I3 tolerance for each term and component."""
    depth = 10.0
    materials = list(MATERIALS.values())
    # The k_perp integrals are shared by both materials and both hole charges.
    references = reference_terms(materials, energy, losses, depth)
    for material, reference in zip(materials, references, strict=True):
        totals = {q: reference["no_hole"] + q * reference["hole"] for q in HOLE_CHARGES}
        scale = max(np.max(abs(total)) for total in totals.values())
        for charge in HOLE_CHARGES:
            result = effective_cross_section_fixed_depth(
                material, energy, losses, depth, theta_deg=0.0, hole_charge=charge, epsilon_momentum="k_parallel",
            )
            expected = {**reference, "total": totals[charge]}
            for name in (*TERMS, "no_hole", "hole", "total"):
                actual, target = getattr(result, name), expected[name]
                np.testing.assert_equal(actual.shape, target.shape)
                tolerance = REL * abs(target) + ABS * scale
                yield material.name, charge, name, abs(actual - target) / tolerance


class EnvelopeEndpointTests(unittest.TestCase):
    def test_terms_and_components_at_both_energies(self):
        for energy, losses in ENDPOINTS:
            for material, charge, name, ratios in endpoint_error_ratios(energy, losses):
                with self.subTest(energy=energy, material=material, hole_charge=charge, term=name):
                    worst = np.argmax(ratios)
                    self.assertLessEqual(ratios[worst], 1.0, f"loss = {losses[worst]} eV")
