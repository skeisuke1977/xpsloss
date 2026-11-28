"""Input validation and exact coordinate arithmetic for the binding axis."""

from dataclasses import replace
from fractions import Fraction as Q
import math
import unittest

import numpy as np

import xpsloss as x
from xpsloss import _loss_orders
from xpsloss._binding_axis import convert


def axis(photon=32., offset=0., **changes):
    values = dict(photon_energy_eV=photon,
                  reference=x.EnergyAxisReference(kind="photon_and_offset", kinetic_energy_offset_eV=offset),
                  primary_reference_choice="synthetic line", primary_reference_binding_eV=0.,
                  source=dict(name="synthetic", monochromatic=True, satellite_treatment="none"),
                  calibration="not_recorded", charge_correction=dict(applied=False, method=None, shift_eV=None),
                  core_level="synthetic", emission_angle_deg=0., acceptance_half_angle_deg=1.)
    return x.BindingEnergyAxis(**(values | changes))


class AxisInput(unittest.TestCase):
    def test_single_rounding_and_reference(self):
        ax = axis(1486.6, 4.6, primary_reference_binding_eV=72.3)
        p = x.BindingDeltaPrimary([72.3], [1])
        edges = [70.1, 72.3, 74.6]
        record = convert(ax, p, edges, 1400.).info
        c = Q(1486.6) - Q(4.6)
        exact = [c - Q(b) for b in reversed(edges)]
        self.assertNotEqual(c, Q(float(c)))
        self.assertTrue(all(float(c) - b != float(c - Q(b)) for b in edges))
        self.assertEqual(record["conversion_constant_exact"], f"{c.numerator}/{c.denominator}")
        self.assertEqual(record["kinetic_output_edges_eV"].tolist(), list(map(float, exact)))
        self.assertEqual(record["electron_energy_minus_reference_eV"], float(1400 - (c - Q(72.3))))
        rounding = max(abs(Q(float(v)) - v) for v in exact)
        self.assertGreaterEqual(Q(record["rounding_max_eV"]), rounding)
        self.assertLess(Q(math.nextafter(record["rounding_max_eV"], -math.inf)), rounding)

    def test_anchor_is_independent_of_photon_and_corrections(self):
        ref = x.EnergyAxisReference(kind="paired_anchor", binding_anchor_eV=-2., kinetic_anchor_eV=18.,
                                    line="not_recorded", citation="not_recorded")
        ax = axis(reference=ref, primary_reference_binding_eV=30., primary_reference_choice="not_recorded",
                  charge_correction=dict(applied=True, method="synthetic", shift_eV=99.),
                  source=dict(name="synthetic", monochromatic="not_recorded", satellite_treatment="not_recorded"))
        p = x.BindingDeltaPrimary([0.], [1.])
        info = convert(ax, p, [0., 1.]).info
        self.assertEqual(info["conversion_constant_exact"], "16/1")
        self.assertEqual(info["primary_reference_kinetic_eV"], -14.)
        self.assertEqual(info["implied_offset_eV"], 16.)
        self.assertEqual(info["not_recorded"], ("reference.line", "reference.citation", "primary_reference_choice",
                         "source.monochromatic", "source.satellite_treatment", "calibration"))
        ax.source["name"] = "modified"
        ax.charge_correction["shift_eV"] = 0.
        self.assertEqual(info["record"]["source"]["name"], "synthetic")
        self.assertEqual(info["record"]["charge_correction"]["shift_eV"], 99.)

    def test_primary_copy_and_stable_permutation(self):
        edges = np.array([3., 1., -2.])
        p = x.BindingBinnedPrimary(edges, [.2, .8])
        edges[0] = 4.
        self.assertEqual(p.binding_edges_eV.tolist(), [3., 1., -2.])
        self.assertFalse(p.binding_edges_eV.flags.writeable)
        self.assertFalse(p.intensities.flags.writeable)
        p = x.BindingDeltaPrimary([1., 3., 1.], [2., 3., 4.], description="synthetic duplicate")
        info = convert(axis(), p, [-0., 9.]).info
        self.assertTrue(np.signbit(info["input_binding_edges_eV"][0]))
        self.assertEqual(info["kinetic_primary"].intensities.tolist(), [3., 2., 4.])
        self.assertIsNone(info["direction"]["primary"])
        self.assertIs(info["input_primary"], p)
        for pos, val in (([0., 2., 1.], [1, 2]), ([1, 1], [1]), ([True, 2], [1]),
                         ([1, 2], [2.**-901]), ([-1, 2], [-1])):
            with self.subTest(pos=pos, val=val), self.assertRaises(x.SingleLossInputError) as cm:
                x.BindingBinnedPrimary(pos, val)
            self.assertEqual(cm.exception.reason, "invalid_primary")

    def test_union_collision_before_kernel_and_window(self):
        p = x.BindingDeltaPrimary([2.+2.**-44, 1.+2.**-45], [1., 0.])
        with self.assertRaises(x.MultipleLossInputError) as cm:
            _loss_orders.integrate_orders(self.fail, (0, 1), features={}, weight_length_A=1.,
                                          primary=p, output_edges_eV=[2., 1., 0.], energy_axis=axis(8192.),
                                          model="unknown", max_order=-1)
        exc = cm.exception
        self.assertEqual(exc.reason, "axis_conversion_not_injective")
        self.assertEqual(exc.roles, ("output_edges_eV[1]", "primary.binding_energies_eV[1]"))
        self.assertEqual(exc.binding_eV, (1., 1.+2.**-45))
        self.assertEqual(exc.kinetic_eV, (8191., 8191.))
        self.assertNotIn("field", vars(exc))

    def test_axis_errors_precede_primary_and_edges(self):
        ref = axis().reference
        cases = [
            (dict(photon_energy_eV=True), "photon_energy_eV", "type"),
            (dict(photon_energy_eV="not_recorded"), "photon_energy_eV", "not_recorded_not_allowed"),
            (dict(reference=replace(ref, line="x")), "reference", "exclusive"),
            (dict(reference=replace(ref, kinetic_energy_offset_eV=None)), "reference.kinetic_energy_offset_eV", "missing"),
            (dict(core_level="   "), "core_level", "empty"),
            (dict(acceptance_half_angle_deg=-1), "acceptance_half_angle_deg", "out_of_range"),
            (dict(charge_correction=dict(applied="not_recorded", method="x", shift_eV=None)),
             "charge_correction.method", "exclusive"),
            (dict(source=dict(name="x", monochromatic=1, satellite_treatment="x")), "source.monochromatic", "type"),
        ]
        for changes, field, detail in cases:
            with self.subTest(field=field), self.assertRaises(x.MultipleLossInputError) as cm:
                convert(axis(**changes), x.DeltaPrimary([1], [1]), [])
            self.assertEqual((cm.exception.reason, cm.exception.field, cm.exception.detail),
                             ("invalid_energy_axis", field, detail))

    def test_unrepresentable_conversion_and_negative_ke(self):
        big = float.fromhex("0x1.fffffffffffffp+1023")
        p = x.BindingDeltaPrimary([0.], [1.])
        for ax, edges, field, detail in (
            (axis(big, -big), [0., 1.], "conversion_constant_eV", "not_finite"),
            (axis(big), [-big, 0.], "output_edges_eV", "not_finite"),
            (axis(), [31., 33.], "output_edges_eV", "negative_kinetic_energy"),
            (axis(big, primary_reference_binding_eV=-big), [0., 1.], "primary_reference_kinetic_eV", "not_finite"),
        ):
            with self.subTest(field=field), self.assertRaises(x.MultipleLossInputError) as cm:
                convert(ax, p, edges)
            self.assertEqual((cm.exception.field, cm.exception.detail), (field, detail))
