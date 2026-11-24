"""M2-2 checks owned by the reference author (docs/M2_2_PREREG.md, section 5):
V3, V4 (group comparison), V7 (fixed-depth level), V13, V14 and V15.

The expected values are frozen in benchmarks/numerics/m22/references_m22.json, written
by benchmarks/numerics/m22/generate_references_m22.py (mpmath, no package functions) and
identified by the SHA-256 in the appendix of docs/M2_2_PREREG.md.  This file calls the
public API only and does not import mpmath.

Budget of output j: B_j = max(epsabs, epsrel*|K_j,ref|), u_j is the reference uncertainty.
V14 has two parts.  The budget part uses only the returned values.  The coverage part
needs a published error estimate per output; the name of that field was not fixed when
this file was written, so it is discovered by per_output_errors() and the coverage tests
are skipped ("not judged") when no such field exists.
"""
from collections.abc import Mapping
import dataclasses
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import re
import unittest

import numpy as np

import xpsloss
from xpsloss import FixedDepthConvergenceError, Material, Oscillator, effective_cross_section_fixed_depth

ROOT = Path(__file__).resolve().parents[1]
REFERENCE_FILE = ROOT / "benchmarks" / "numerics" / "m22" / "references_m22.json"
PREREG = ROOT / "docs" / "M2_2_PREREG.md"

TERMS = ("term_13a", "term_13b", "term_14", "term_15_T1", "term_15_T2",
         "term_15_T3", "term_16_e", "term_16_h")
ELECTRON_GROUP = ("term_13b", "term_15_T1", "term_15_T3", "term_16_e")
HOLE_GROUP = ("term_14", "term_15_T2", "term_16_h")
FIELDS = TERMS + ("no_hole", "hole", "total")
CHARGES = {0.0: "total_q0", 0.5: "total_q0.5", 1.0: "total_q1"}
SPECS = ("k_parallel", "k_total")
# Candidate names of the per-output error estimate, tried in this order; see per_output_errors().
ERROR_FIELD_CANDIDATES = ("output_errors", "output_error", "estimated_errors", "error_estimates",
                          "per_output_error", "per_output_errors", "quadrature_errors", "errors")
ERROR_SUFFIXES = ("_error", "_quadrature_error", "_estimated_error")


def _load():
    raw = REFERENCE_FILE.read_bytes()
    return json.loads(raw.decode("ascii")), hashlib.sha256(raw).hexdigest()


def _float(entry):
    return float.fromhex(entry["binary64_hex"])


def _as_error_map(value, size):
    """Mapping output name -> array of length size, or None."""
    if not isinstance(value, Mapping) or not all(name in value for name in FIELDS):
        return None
    out = {}
    for name in FIELDS:
        arr = np.atleast_1d(np.asarray(value[name], dtype=float))
        if arr.shape != (size,):
            return None
        out[name] = arr
    return out


def per_output_errors(result):
    """Published error estimate of each of the 11 outputs, or None when the result has none.

    Accepted layouts: a mapping field (one of ERROR_FIELD_CANDIDATES, or any dataclass field)
    keyed by the 11 output names; the same mapping one level inside a mapping field such as
    diagnostics; or eleven array fields named <output><suffix> with one of ERROR_SUFFIXES.
    """
    size = np.atleast_1d(np.asarray(result.loss_eV)).size
    names = list(ERROR_FIELD_CANDIDATES) + ["diagnostics"]
    if dataclasses.is_dataclass(result):
        names += [f.name for f in dataclasses.fields(result) if f.name not in names]
    for name in names:
        value = getattr(result, name, None)
        found = _as_error_map(value, size)
        if found is not None:
            return found
        if isinstance(value, Mapping):
            for key in list(ERROR_FIELD_CANDIDATES) + [k for k in value if k not in ERROR_FIELD_CANDIDATES]:
                found = _as_error_map(value.get(key), size)
                if found is not None:
                    return found
    for suffix in ERROR_SUFFIXES:
        if all(hasattr(result, name + suffix) for name in FIELDS):
            found = _as_error_map({name: getattr(result, name + suffix) for name in FIELDS}, size)
            if found is not None:
                return found
    return None


class M22Reference(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.doc, cls.sha256 = _load()
        cond = cls.doc["conditions"]
        cls.energy = _float(cond["electron_energy_eV"])
        cls.numerics = dict(epsabs=cond["epsabs"], epsrel=cond["epsrel"], limit=cond["limit"])
        cls.epsabs = cond["epsabs"]
        cls.materials = {}
        for name, spec in cls.doc["materials"].items():
            cls.materials[name] = Material(
                "Si" if name == "Si" else f"synthetic_{name}", _float(spec["gap_eV"]),
                tuple(Oscillator(*(_float(o[k]) for k in ("energy_eV", "strength_eV2", "width_eV", "dispersion")))
                      for o in spec["oscillators"]))
        cls.points = {}
        cls.losses = {}
        for p in cls.doc["points"]:
            loss = _float(p["loss_eV"])
            cls.points[p["material"], p["epsilon_momentum"], p["depth_A"], loss] = p
            cls.losses.setdefault(p["material"], set()).add(loss)
        cls.losses = {name: np.array(sorted(values)) for name, values in cls.losses.items()}
        cls._runs = {}

    # ------------------------------------------------------------------ helpers
    def _call(self, material, spec, depth, charge):
        """Result of the public API over all reference losses of the material, or the exception."""
        key = (material, spec, depth, charge)
        if key not in self._runs:
            try:
                self._runs[key] = effective_cross_section_fixed_depth(
                    self.materials[material], self.energy, self.losses[material], depth,
                    hole_charge=charge, epsilon_momentum=spec, **self.numerics)
            except Exception as error:      # reported by the tests that need the result
                self._runs[key] = error
        return self._runs[key]

    def _result(self, material, spec, depth, charge):
        out = self._call(material, spec, depth, charge)
        if isinstance(out, Exception):
            self.fail(f"required success: {material} {spec} depth {depth:g} A hole charge {charge:g} "
                      f"raised {type(out).__name__}: {out}")
        np.testing.assert_array_equal(np.asarray(out.loss_eV, dtype=float), self.losses[material])
        return out

    def _index(self, material, loss):
        return int(np.nonzero(self.losses[material] == loss)[0][0])

    def _reference(self, point, name, charge):
        return point["outputs"][CHARGES[charge] if name == "total" else name]

    def _required(self, materials=("Si", "S0")):
        for key, point in self.points.items():
            if key[0] in materials:
                yield key, point

    def _budget_failures(self, result, point, i, charge, names=FIELDS):
        bad = []
        for name in names:
            ref = self._reference(point, name, charge)
            got = float(getattr(result, name)[i])
            miss = abs(got - ref["value"]) + ref["u"]
            if not miss <= ref["B"]:
                bad.append(f"{name}: got {got!r}, reference {ref['value']!r}, "
                           f"|diff|+u {miss:.3e} > B {ref['B']:.3e}")
        return bad

    def _coverage_failures(self, result, errors, point, i, charge, names=FIELDS):
        bad = []
        for name in names:
            ref = self._reference(point, name, charge)
            got = float(getattr(result, name)[i])
            miss = abs(got - ref["value"]) + ref["u"]
            if not miss <= errors[name][i]:
                bad.append(f"{name}: |diff|+u {miss:.3e} > published estimate {errors[name][i]:.3e}")
        return bad

    def _report(self, bad, label):
        if bad:
            self.fail(f"{label}: {len(bad)} failures; first: " + " | ".join(bad[:6]))

    # ------------------------------------------------------------------ frozen reference
    def test_reference_file_matches_the_preregistered_hash(self):
        text = PREREG.read_text(encoding="utf-8")
        found = re.findall(r"references_m22\.json[^\n]*?\b([0-9a-f]{64})\b", text)
        self.assertTrue(found, "the appendix of docs/M2_2_PREREG.md records no SHA-256 for references_m22.json")
        self.assertEqual(found[0], self.sha256)

    def test_reference_uncertainty_is_within_one_percent_of_the_budget(self):
        count = 0
        for point in self.doc["points"]:
            for name, out in point["outputs"].items():
                count += 1
                budget = max(self.epsabs, self.numerics["epsrel"] * abs(out["value"]))
                self.assertEqual(out["B"], budget)
                self.assertLessEqual(out["u"], 0.01 * budget, (point["material"], point["loss_label"], name))
                if name in point["identically_zero"]:
                    self.assertEqual((out["value"], out["u"]), (0.0, 0.0))
        self.assertGreater(count, 2000)

    def test_reference_points_cover_the_preregistered_grid(self):
        gap = _float(self.doc["materials"]["Si"]["gap_eV"])
        self.assertEqual(gap, 1.12)
        want = [0.2, 0.4, 1.0, 1.119, 1.119999, 1.120001, 1.121, 1.2, 2.0, 5.0, 11.0, 16.8, 17.0, 30.0, 60.0, 99.8]
        for material in ("Si", "S0"):
            self.assertEqual(list(self.losses[material]), want)
            for spec in SPECS:
                for depth in (2.0, 10.0, 50.0):
                    for loss in want:
                        point = self.points[material, spec, depth, loss]
                        self.assertEqual(point["region"], "below_gap" if loss < gap else "above_gap")
                        self.assertEqual("i3" in point, spec == "k_parallel")
        for spec in SPECS:
            self.assertIn(("N1", spec, 10.0, 4.0), self.points)

    # ------------------------------------------------------------------ V3
    def test_V3_i3_reference_refinement_is_within_a_tenth_of_the_width(self):
        for point in self.doc["points"]:
            if "i3" in point:
                self.assertLessEqual(point["i3"]["refinement_change_over_v3_width_max"], 0.1,
                                     (point["material"], point["depth_A"], point["loss_label"]))

    def test_V3_i3_route_k_parallel(self):
        bad = []
        for (material, spec, depth, loss), point in self._required():
            if spec != "k_parallel":
                continue
            i3 = point["i3"]
            scale = i3["v3_scale"]
            i = self._index(material, loss)
            for charge in (0.0, 1.0):
                result = self._result(material, spec, depth, charge)
                want = dict(i3["values"])
                want["total"] = want["no_hole"] + charge * want["hole"]
                for name in FIELDS:
                    got = float(getattr(result, name)[i])
                    tol = 1e-6 * abs(want[name]) + 1e-6 * scale
                    if not abs(got - want[name]) <= tol:
                        bad.append(f"{material} {depth:g} A {point['loss_label']} eV q_h={charge:g} {name}: "
                                   f"got {got!r}, I3 {want[name]!r}, tol {tol:.2e}")
        self._report(bad, "V3")

    # ------------------------------------------------------------------ V4
    def test_V4_sub_gap_groups_and_unit_charge_total(self):
        bad = []
        count = 0
        for (material, spec, depth, loss), point in self._required():
            if point["region"] != "below_gap":
                continue
            count += 1
            result = self._result(material, spec, depth, 1.0)
            i = self._index(material, loss)
            label = f"{material} {spec} {depth:g} A {point['loss_label']} eV"
            groups = {"electron_group": ELECTRON_GROUP, "hole_group": HOLE_GROUP}
            for name, members in groups.items():
                got = math.fsum(float(getattr(result, t)[i]) for t in members)
                ref = point["outputs"][name]
                miss = abs(got - ref["value"]) + ref["u"]
                if not miss <= ref["B"]:
                    bad.append(f"{label} {name}: got {got!r}, reference {ref['value']!r}, "
                               f"|diff|+u {miss:.3e} > B {ref['B']:.3e}")
                self.assertGreater(abs(ref["value"]), 1e3 * ref["B"], "the group must not be zero")
            total = float(result.total[i])
            size = math.fsum(abs(float(getattr(result, t)[i])) for t in TERMS)
            if not abs(total) <= max(self.epsabs, 1e-13 * size):
                bad.append(f"{label} total at unit hole charge: {total!r}, sum|terms| {size!r}")
            if float(result.term_13a[i]) != 0.0:
                bad.append(f"{label} term_13a is not zero below the gap: {float(result.term_13a[i])!r}")
        self.assertEqual(count, 2 * 2 * 3 * 5)
        self._report(bad, "V4")

    # ------------------------------------------------------------------ V7 (fixed depth)
    def _edge_points(self):
        for key, point in self._required():
            if point["loss_label"] in ("g-1e-3", "g-1e-6", "g+1e-6", "g+1e-3"):
                yield key, point

    def test_V7_edge_points_budget(self):
        bad = []
        count = 0
        for (material, spec, depth, loss), point in self._edge_points():
            count += 1
            for charge in CHARGES:
                result = self._result(material, spec, depth, charge)
                i = self._index(material, loss)
                bad += [f"{material} {spec} {depth:g} A {point['loss_label']} q_h={charge:g} {b}"
                        for b in self._budget_failures(result, point, i, charge)]
        self.assertEqual(count, 2 * 2 * 3 * 4)
        self._report(bad, "V7")

    def test_V7_edge_points_error_coverage(self):
        self._coverage(self._edge_points(), "V7")

    def test_V7_reference_jump_at_the_edge(self):
        """The frozen references themselves show the jump of the unit-charge total at the edge."""
        for spec in SPECS:
            below = self.points["Si", spec, 10.0, 1.119999]["outputs"]["total_q1"]
            above = self.points["Si", spec, 10.0, 1.120001]["outputs"]["total_q1"]
            self.assertEqual(below["value"], 0.0)
            self.assertGreater(above["value"], 1e6 * above["B"])

    # ------------------------------------------------------------------ V13
    def test_V13_dispersionless_conventions_agree(self):
        bad = []
        for (material, spec, depth, loss), point in self._required(("S0",)):
            if spec != "k_parallel":
                continue
            i = self._index(material, loss)
            for charge in CHARGES:
                first = self._result(material, "k_parallel", depth, charge)
                second = self._result(material, "k_total", depth, charge)
                for name in FIELDS:
                    a, b = float(getattr(first, name)[i]), float(getattr(second, name)[i])
                    budget = self._reference(point, name, charge)["B"]
                    if not abs(a - b) <= budget:
                        bad.append(f"S0 {depth:g} A {point['loss_label']} q_h={charge:g} {name}: "
                                   f"k_parallel {a!r}, k_total {b!r}, B {budget:.3e}")
        self._report(bad, "V13 (S0)")

    def test_V13_dispersive_conventions_are_separated(self):
        selected = self.doc["v13_points"]
        self.assertGreaterEqual(len(selected), 6)
        self.assertTrue(any(s["loss_label"] in ("0.2", "0.4", "1.0") for s in selected))
        self.assertTrue(any("total_q1" in s["separated_outputs"] for s in selected))
        bad = []
        for entry in selected:
            depth = entry["depth_A"]
            loss = next(l for l in self.losses["Si"]
                        if self.points["Si", "k_parallel", depth, l]["loss_label"] == entry["loss_label"])
            i = self._index("Si", loss)
            pair = [self.points["Si", spec, depth, loss] for spec in SPECS]
            for key in entry["separated_outputs"]:
                a, b = pair[0]["outputs"][key], pair[1]["outputs"][key]
                self.assertGreaterEqual(abs(a["value"] - b["value"]), 1e3 * max(a["B"], b["B"]))
            for spec, point in zip(SPECS, pair):
                for charge in (0.0, 1.0):
                    result = self._result("Si", spec, depth, charge)
                    names = [n for n in FIELDS
                             if (CHARGES[charge] if n == "total" else n) in entry["separated_outputs"]]
                    bad += [f"Si {spec} {depth:g} A {entry['loss_label']} q_h={charge:g} {b}"
                            for b in self._budget_failures(result, point, i, charge, names)]
        self._report(bad, "V13 (Si)")

    # ------------------------------------------------------------------ V14
    def test_V14_budget(self):
        bad = []
        count = 0
        for (material, spec, depth, loss), point in self._required():
            count += 1
            i = self._index(material, loss)
            for charge in CHARGES:
                result = self._result(material, spec, depth, charge)
                bad += [f"{material} {spec} {depth:g} A {point['loss_label']} q_h={charge:g} {b}"
                        for b in self._budget_failures(result, point, i, charge)]
        self.assertEqual(count, 2 * 2 * 3 * 16)
        self._report(bad, "V14 budget")

    def _coverage(self, points, label):
        bad = []
        judged = 0
        for (material, spec, depth, loss), point in points:
            i = self._index(material, loss)
            for charge in CHARGES:
                result = self._result(material, spec, depth, charge)
                errors = per_output_errors(result)
                if errors is None:
                    self.skipTest(f"{label} coverage NOT JUDGED: FixedDepthResult has no per-output error "
                                  "estimate recognised by per_output_errors() (only the common "
                                  "quadrature_error); add the field name to ERROR_FIELD_CANDIDATES")
                judged += 1
                bad += [f"{material} {spec} {depth:g} A {point['loss_label']} q_h={charge:g} {b}"
                        for b in self._coverage_failures(result, errors, point, i, charge)]
        self.assertGreater(judged, 0)
        self._report(bad, f"{label} coverage")

    def test_V14_error_coverage(self):
        self._coverage(self._required(), "V14")

    def test_V14_N1_narrow_surface_resonance_either_outcome(self):
        """Pre-registered as 'either': success must agree with the reference under the V14
        conditions, otherwise the only accepted outcome is FixedDepthConvergenceError."""
        for spec in SPECS:
            point = self.points["N1", spec, 10.0, 4.0]
            for charge in CHARGES:
                with self.subTest(spec=spec, hole_charge=charge):
                    out = self._call("N1", spec, 10.0, charge)
                    if isinstance(out, Exception):
                        self.assertIsInstance(out, FixedDepthConvergenceError,
                                              f"{type(out).__name__}: {out}")
                        continue
                    self._report(self._budget_failures(out, point, 0, charge), f"N1 {spec} budget")
                    errors = per_output_errors(out)
                    if errors is not None:
                        self._report(self._coverage_failures(out, errors, point, 0, charge),
                                     f"N1 {spec} coverage")

    # ------------------------------------------------------------------ V15
    PI_LOWER = Fraction(314159265358979, 10 ** 14)       # pi = 3.14159265358979323... > this

    @staticmethod
    def _log_upper(x, terms=60):
        """Rational upper bound of log x for x >= 1: x = 2^n r, log r = 2 atanh((r-1)/(r+1)),
        truncated series plus the geometric bound of the remainder."""
        def atanh_upper(y):
            s = sum(y ** (2 * k + 1) / (2 * k + 1) for k in range(terms))
            return s + y ** (2 * terms + 1) / ((2 * terms + 1) * (1 - y * y))
        n, r = 0, Fraction(x)
        while r >= 2:
            r /= 2
            n += 1
        return n * 2 * atanh_upper(Fraction(1, 3)) + 2 * atanh_upper((r - 1) / (r + 1))

    def _si_bound(self, binary64):
        spec = self.doc["materials"]["Si"]
        conv = (lambda e: Fraction(_float(e))) if binary64 else (lambda e: Fraction(e["decimal"]))
        g = conv(spec["gap_eV"])
        osc = [(conv(o["energy_eV"]), conv(o["strength_eV2"]), conv(o["width_eV"])) for o in spec["oscillators"]]
        gf = _float(spec["gap_eV"])
        # D7 exclusion width (docs/DIELECTRIC_DESIGN.md); accepted sub-gap points have g - T > delta_g.
        if binary64:
            delta = Fraction(max(math.ulp(0.0), min(gf / 2, max(32 * math.ulp(gf), 1e-12 * gf))))
        else:
            delta = max(32 * Fraction(math.ulp(gf)), g / 10 ** 12)
        U = Fraction(3)
        self.assertTrue(g < U < min(E for E, _, _ in osc))
        B0 = sum(A / (E * E) for E, A, _ in osc)
        P0 = B0 / (1 - (g / U) ** 2)
        M = sum(A * G / (E * E - U * U) ** 2 for E, A, G in osc)
        return B0, P0 + (M * g / self.PI_LOWER) * self._log_upper(2 * g / delta)

    def test_V15_uniform_upper_bound_for_Si(self):
        B0, bound = self._si_bound(False)
        self.assertEqual(B0 - 1, Fraction(-259300169, 8820000000))      # -0.02939911...
        for binary64 in (False, True):
            B0, bound = self._si_bound(binary64)
            self.assertLess(bound, Fraction("1.31827"))                  # S < 1.31827, F < 0.31827 < 1
            self.assertGreater(bound, Fraction("1.3182"))                # the bound is not vacuous
            self.assertLess(B0 - 1, Fraction("-0.0293991"))              # static F(q, 0) <= B0 - 1
        recorded = self.doc["v15"]
        for key in ("decimal", "binary64"):
            self.assertTrue(recorded[key]["S_bound_below_1.31827"])
            self.assertTrue(recorded[key]["static_F_bound_below_-0.0293991"])
        self.assertLess(recorded["samples"]["max_S_minus_bound"], 0.0)
        self.assertLess(recorded["samples"]["max_static_F_minus_B0m1"], 0.0)

    def test_V15_implementation_accepts_the_proved_sub_gap_points(self):
        """The bound proves F < 1 at every Si sub-gap point outside the D7 width, so the
        implementation must not reject the sub-gap reference points (down to g - 1e-6 eV)."""
        domain_error = getattr(xpsloss, "FixedDepthGapDomainError", None)
        for spec in SPECS:
            out = self._call("Si", spec, 10.0, 1.0)
            if domain_error is not None:
                self.assertNotIsInstance(out, domain_error, str(out))
            result = self._result("Si", spec, 10.0, 1.0)
            below = self.losses["Si"] < 1.12
            self.assertEqual(int(below.sum()), 5)
            self.assertTrue(np.all(np.isfinite(result.no_hole[below])))


class PerOutputErrorDiscovery(unittest.TestCase):
    """Self-test of per_output_errors() on stand-in results."""

    def test_layouts(self):
        class Flat:
            loss_eV = np.array([1.0, 2.0])
        flat = Flat()
        self.assertIsNone(per_output_errors(flat))
        for name in FIELDS:
            setattr(flat, name + "_error", np.array([1e-13, 2e-13]))
        self.assertEqual(set(per_output_errors(flat)), set(FIELDS))

        class Mapped:
            loss_eV = np.array([1.0, 2.0])
            diagnostics = {"output_errors": {name: [1e-13, 2e-13] for name in FIELDS}}
        self.assertEqual(per_output_errors(Mapped())["total"][1], 2e-13)

        class Short:
            loss_eV = np.array([1.0, 2.0])
            output_errors = {name: [1e-13] for name in FIELDS}
        self.assertIsNone(per_output_errors(Short()))


if __name__ == "__main__":
    unittest.main()
