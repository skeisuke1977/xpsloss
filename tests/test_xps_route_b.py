"""Test I3 of docs/M2_FIRST_CASE.md section 4.1: independent second route.

The package's closed forms (Simonsen, Yubero & Tougaard, PRB 56, 1612 (1997),
Eqs. (13)-(16), theta = 0) are compared term by term with a route that evaluates Eq. (2)
with Eqs. (3), (6), (8)-(12) and integrates k_perp numerically over (-inf, inf)
(helper module xps_route_b_reference.py, written without the package implementation).

The second route is valid only for eps independent of k_perp (constraint of Eq. (7)),
so the inputs are a dispersionless oscillator and a dispersion-1 oscillator with
epsilon_momentum="k_parallel".  Grid: theta = 0, E0 = 1000 eV, a = 2, 10, 50 A,
hbar*w = 5, 11, 15.6, 30 eV, hole charge 1 and 0.
Criterion for every term and component: |diff| <= 1e-6 |value| + 1e-6 max|total|, where
value is the second-route value and max|total| is the largest |total| of the second route
for the same material and depth over the four losses and both hole charges.
"""
from pathlib import Path
import sys
import unittest

import numpy as np

from xpsloss.model import Material, Oscillator
from xpsloss.xps_fixed_depth import effective_cross_section_fixed_depth

sys.path.insert(0, str(Path(__file__).resolve().parent))
import xps_route_b_reference as ref  # noqa: E402

E0 = 1000.0
DEPTHS = (2.0, 10.0, 50.0)
LOSSES = np.array([5.0, 11.0, 15.6, 30.0])
HOLE_CHARGES = (1.0, 0.0)
REL = 1e-6
ABS = 1e-6
MATERIALS = {
    "dispersion0": Material("synthetic_al_dispersion0", 0.0, (Oscillator(15.6, 243.3598, 1.5, 0.0),)),
    "dispersion1": Material("synthetic_al_dispersion1", 0.0, (Oscillator(15.6, 243.3598, 1.5, 1.0),)),
}


def _field(result, name):
    value = getattr(result, name) if hasattr(result, name) else result[name]
    return np.asarray(value, dtype=float)


class SecondRouteI3(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        keys = list(MATERIALS)
        mats = [MATERIALS[k] for k in keys]
        cls.ref = {}
        cls.scale = {}
        for a in DEPTHS:
            for key, r in zip(keys, ref.reference_terms(mats, E0, LOSSES, a)):
                cls.ref[key, a] = r
                cls.scale[key, a] = max(np.max(np.abs(r["no_hole"] + q * r["hole"])) for q in HOLE_CHARGES)
        cls.pkg = {}
        for key in keys:
            for a in DEPTHS:
                for q in HOLE_CHARGES:
                    cls.pkg[key, a, q] = effective_cross_section_fixed_depth(
                        MATERIALS[key], E0, LOSSES, a, theta_deg=0.0, hole_charge=q,
                        epsilon_momentum="k_parallel")

    def _compare(self, got, want, scale, label):
        got = np.asarray(got, dtype=float)
        want = np.asarray(want, dtype=float)
        self.assertEqual(got.shape, want.shape, label)
        tol = REL * np.abs(want) + ABS * scale
        diff = np.abs(got - want)
        bad = ~(diff <= tol)
        if np.any(bad):
            rows = "; ".join(
                f"{w:g} eV: package {g:.12e}, route 2 {v:.12e}, |diff| {d:.2e} > tol {t:.2e}"
                for w, g, v, d, t in zip(LOSSES[bad], got[bad], want[bad], diff[bad], tol[bad]))
            self.fail(f"{label}: {rows}")

    def _check_terms(self, key):
        for a in DEPTHS:
            r = self.ref[key, a]
            for q in HOLE_CHARGES:
                p = self.pkg[key, a, q]
                for name in ref.TERMS:   # hole terms are for unit hole charge
                    with self.subTest(material=key, depth=a, hole_charge=q, term=name):
                        self._compare(_field(p, name), r[name], self.scale[key, a],
                                      f"{key} a={a:g} q_hole={q:g} {name}")

    def test_terms_dispersion0(self):
        self._check_terms("dispersion0")

    def test_terms_dispersion1_k_parallel(self):
        self._check_terms("dispersion1")

    def test_components(self):
        for key in MATERIALS:
            for a in DEPTHS:
                r = self.ref[key, a]
                for q in HOLE_CHARGES:
                    p = self.pkg[key, a, q]
                    want = {"no_hole": r["no_hole"], "hole": r["hole"], "total": r["no_hole"] + q * r["hole"]}
                    for name, value in want.items():
                        with self.subTest(material=key, depth=a, hole_charge=q, component=name):
                            self._compare(_field(p, name), value, self.scale[key, a],
                                          f"{key} a={a:g} q_hole={q:g} {name}")

    def test_reference_route_reproduces_eqs_8_to_10(self):
        """Self-check of the helper: its numerical F integrals equal the printed right-hand
        sides of Eqs. (8)-(10) (theta = 0, pi*delta(w) = 0 for w > 0)."""
        kp = ref.KPerp(E0, 15.6, 10.0)
        e, V, Q, a, w = ref.ECH, kp.V, kp.Q, kp.a, kp.w
        for q in (0.003, 0.05, 0.3, 3.0, 30.0):
            I = kp.integrals(q)
            f_ee = 1j * e / V * (np.exp(1j * Q * a) - np.exp(-q * a)) / (Q - 1j * q)
            f_hh = e * (1j / w) * np.exp(-q * a)
            f_v = -1j * e / V * np.exp(1j * Q * a) / (Q + 1j * q)
            self.assertLess(abs(I["F_ee"] - f_ee), 1e-12 * abs(f_ee))
            self.assertLess(abs(I["F_V"] - f_v), 1e-12 * abs(f_v))
            self.assertLess(abs(I["F_hh"] - f_hh), 1e-12 * e / w)


if __name__ == "__main__":
    unittest.main()
