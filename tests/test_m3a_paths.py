"""M3a checks owned by R2 (docs/M3A_PREREG.md sections 5 and 6): V2, V8 and V9.

V2 runs the depth quadrature of path 1 (tests/m3a_paths.py) on test functions with known
averages.  V8 and V9 compare the depth-averaged implementation with the frozen path values
in benchmarks/numerics/m3a/paths_m3a.json (written by
benchmarks/numerics/m3a/generate_paths_m3a.py; its SHA-256 is in the appendix of the
pre-registration).  Routine: V8 at the three required points, V9 at lambda = 20 A.
The full sets run when the environment variable XPSLOSS_M3A_FULL is 1.

V8: D_j = 1e-7 |R_j| + 1e-8 max_{T in P} |R_total| (frozen scale), path u_j <= D_j/3,
comparison |Delta_j| + u_j <= D_j.  V9: |Delta_j| <= 1e-6 |K_j,I3| + 1e-6 max_9 |K_total,I3|,
and every one-axis change of path 3 is at most 1/10 of that width.
"""
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import m3a_paths as mp3  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PATHS_FILE = ROOT / "benchmarks" / "numerics" / "m3a" / "paths_m3a.json"
TABLE_FILE = ROOT / "benchmarks" / "numerics" / "m3a" / "expected_states_m3a.json"
TABLE_SHA256 = "d87dad8f02a2fd211679229bd7c5205d908dbb2b3b241b30a27e94469f414ed5"
PREREG = ROOT / "docs" / "M3A_PREREG.md"
FULL = os.environ.get("XPSLOSS_M3A_FULL") == "1"
REQUIRED = (("p1_Al_1000_k_total_Q1_L80", 99.8), ("p1_Al_1000_k_total_Q1_L20", 1.0),
            ("p1_Al_1000_k_total_Q1_L20", 15.6))
DEF = dict(epsabs=1e-12, epsrel=1e-9, limit=200)


def _table():
    raw = TABLE_FILE.read_bytes()
    if hashlib.sha256(raw).hexdigest() != TABLE_SHA256:
        raise AssertionError("expected-state table SHA-256 differs from docs/M3A_PREREG.md section 3")
    return json.loads(raw.decode("utf-8"))


def _registered_paths_sha():
    text = PREREG.read_text("utf-8")
    found = re.findall(r"paths_m3a\.json`[^`\n]*`([0-9a-f]{64})`", text)
    return found[-1] if found else None


def _paths():
    raw = PATHS_FILE.read_bytes()
    sha = hashlib.sha256(raw).hexdigest()
    registered = _registered_paths_sha()
    if registered is None:
        raise AssertionError("docs/M3A_PREREG.md has no SHA-256 of paths_m3a.json in the appendix")
    if sha != registered:
        raise AssertionError(f"paths_m3a.json SHA-256 {sha} differs from the appendix {registered}")
    return json.loads(raw.decode("ascii"))


def _implementation():
    from xpsloss.xps_depth_averaged import effective_cross_section_depth_averaged
    return effective_cross_section_depth_averaged


def _set(paths, set_id):
    for s in paths["path1"]:
        if s["id"] == set_id:
            return s
    raise AssertionError(f"path-1 set {set_id} missing from paths_m3a.json")


def _point(entry, loss):
    for p in entry["points"]:
        if p["loss_hex"] == float(loss).hex():
            return p
    raise AssertionError(f"{entry['id']}: loss {loss!r} missing")


class V2DepthQuadrature(unittest.TestCase):
    """V2: the path-1 depth quadrature with 12 nodes per piece on test functions."""

    ORDER = 12
    TOL = 1e-12

    @classmethod
    def setUpClass(cls):
        table = _table()
        ga = 0.2 * np.arange(1, 500)
        cls.layouts = []
        ms5 = sorted({float.fromhex(h) for r in table["runs"] if r["family"] == "MS5" for h in r["losses"]["values_hex"]})
        ms11 = sorted({float.fromhex(h) for r in table["runs"] if r["family"] == "MS11" and r["material"] == "Si"
                       for h in r["losses"]["values_hex"]})
        cls.layouts = [("MS5 1000 eV", ms5, 1000.0, (5.0, 20.0, 80.0)), ("MS11 1000 eV", ms11, 1000.0, (5.0, 20.0, 80.0)),
                       ("MS2 4000 eV", [float(ga[i]) for i in (4, 77, 149, 498)], 4000.0, (5.0, 20.0, 80.0))]

    def test_v2_i_truncated_cumulative(self):
        r = 8.0
        for name, losses, e0, lams in self.layouts:
            for lam in lams:
                rule = mp3.depth_rule(losses, e0, lam, order=self.ORDER, r=r)
                with self.subTest(layout=name, lam=lam):
                    one = mp3.quadrature_average(rule, lambda x: np.ones_like(x))
                    inv = mp3.quadrature_average(rule, lambda x: 1.0 / x) * lam
                    self.assertLessEqual(abs(one - (1 - (1 + r) * math.exp(-r))), self.TOL)
                    self.assertLessEqual(abs(inv - (1 - math.exp(-r))), self.TOL)

    def test_v2_ii_shallow_part(self):
        for name, losses, e0, lams in self.layouts:
            for lam in lams:
                for h in (0.01, 0.1):
                    rule = mp3.depth_rule(losses, e0, lam, order=self.ORDER, extra_edges=(h,))
                    self.assertIn(h, rule.edges)
                    mask = rule.t < h
                    part = math.fsum((rule.alpha[mask] / (lam * rule.t[mask])).tolist()) * lam
                    with self.subTest(layout=name, lam=lam, h=h):
                        self.assertLessEqual(abs(part - (1 - math.exp(-h))), self.TOL)

    def test_v2_iii_full_averages(self):
        for name, losses, e0, lams in self.layouts:
            for lam in lams:
                rule = mp3.depth_rule(losses, e0, lam, order=self.ORDER)
                for fname, (func, exact, scale) in mp3.v2_functions(lam).items():
                    got = mp3.quadrature_average(rule, func)
                    with self.subTest(layout=name, lam=lam, function=fname):
                        self.assertLessEqual(abs((got - exact) * scale), self.TOL)

    def test_rule_layout(self):
        """The registered edges, the geometric edges and the phase rule of section 5."""
        for name, losses, e0, lams in self.layouts:
            for lam in lams:
                rule = mp3.depth_rule(losses, e0, lam)
                k_lo, k_hi = mp3.kpar_limits(e0, np.asarray(losses))
                b = np.asarray(losses) / mp3.hbar_v(e0)
                t_c = mp3.DECAY / (k_lo * lam)
                with self.subTest(layout=name, lam=lam):
                    self.assertTrue(set(mp3.REGISTERED_EDGES) <= set(rule.edges))
                    floor = 1 / (8 * k_hi.max() * lam)
                    geo = [math.ldexp(40.0, -m) for m in range(1, 80) if math.ldexp(40.0, -m) >= floor]
                    self.assertTrue(set(geo) <= set(rule.edges))
                    self.assertLess(min(e for e in rule.edges if e > 0) / 2, floor)
                    for lo, hi in rule.pieces:
                        active = b[lo < t_c]
                        beff = active.max() if active.size else 0.0
                        self.assertLessEqual(beff * lam * (hi - lo), math.pi / 4 * (1 + 1e-12))
                    self.assertEqual(rule.n, 10 * len(rule.pieces))
                    self.assertAlmostEqual(rule.pieces[-1][1], 40.0, places=12)


class V8PathOne(unittest.TestCase):
    """V8: path 1 against the implementation (routine: the three required points)."""

    @classmethod
    def setUpClass(cls):
        cls.paths = _paths()
        cls.table = _table()

    def _check_frozen(self, entry, point):
        self.assertEqual(point["status"], "ok", f"{entry['id']} {point['loss_eV']}: path 1 {point.get('reasons')}")
        scale = entry["scale_max_abs_R_total"]
        for name in mp3.OUTPUTS:
            col = point["columns"][name]
            d = mp3.V8_REL * abs(col["value"]) + mp3.V8_ABS * scale
            self.assertEqual(d, col["D"], f"{entry['id']} {name}: frozen D_j")
            self.assertLessEqual(col["u"], col["D"] / 3, f"{entry['id']} {point['loss_eV']} {name}: u_j > D_j/3")
            self.assertLessEqual(col["allocation_max_ratio"], 1.0,
                                 f"{entry['id']} {point['loss_eV']} {name}: alpha_i e_ij > f_i eta D_j")
            parts = max(abs(col["delta_A"]), abs(col["delta_B"]), abs(col["delta_C"])) + col["node_error"] \
                + col["rounding"] + col["tail"]
            self.assertLessEqual(parts, col["u"], f"{entry['id']} {name}: u_j smaller than its parts")

    def _compare(self, entry, losses, result, index_of):
        failures = []
        for t in losses:
            point = _point(entry, t)
            if point["status"] != "ok":
                continue
            i = index_of(t)
            for name in mp3.OUTPUTS:
                col = point["columns"][name]
                got = float(np.asarray(getattr(result, name))[i])
                delta = abs(got - col["value"])
                if not delta + col["u"] <= col["D"]:
                    failures.append(f"{entry['id']} T={t!r} {name}: |Delta| {delta:.3e} + u {col['u']:.3e} > D "
                                    f"{col['D']:.3e}")
        return failures

    def test_required_points_frozen_path(self):
        for set_id, loss in REQUIRED:
            entry = _set(self.paths, set_id)
            with self.subTest(set=set_id, loss=loss):
                self._check_frozen(entry, _point(entry, loss))

    def test_required_points_against_implementation(self):
        f = _implementation()
        failures = []
        for set_id, loss in REQUIRED:
            entry = _set(self.paths, set_id)
            mat = mp3.material_from_table(entry["material"], self.table)
            # The registered MS5 run key: one loss per call, DEF passed explicitly.
            result = f(mat, entry["electron_energy_eV"], [loss], weight_length_A=entry["weight_length_A"],
                       epsilon_momentum=entry["epsilon_momentum"], hole_charge=entry["hole_charge"], **DEF)
            failures += self._compare(entry, [loss], result, lambda t: 0)
        self.assertEqual(failures, [])

    @unittest.skipUnless(FULL, "full V8 set: set XPSLOSS_M3A_FULL=1")
    def test_full_set_against_implementation(self):
        """Every point of every set P.  Path failures are listed and are not passes; the required
        points must succeed.  Implementation values come from the registered runs (MS5, MS11 single
        losses; MS2 the GA array)."""
        f = _implementation()
        failures, path_not_ok, compared = [], [], 0
        ga = 0.2 * np.arange(1, 500)
        for entry in self.paths["path1"]:
            mat = mp3.material_from_table(entry["material"], self.table)
            losses = [float.fromhex(h) for h in entry["losses_hex"]]
            kw = dict(weight_length_A=entry["weight_length_A"], epsilon_momentum=entry["epsilon_momentum"],
                      hole_charge=entry["hole_charge"], **DEF)
            for p in entry["points"]:
                if p["status"] != "ok":
                    path_not_ok.append((entry["id"], p["loss_eV"], p["status"]))
            if entry["electron_energy_eV"] == 4000.0:
                result = f(mat, 4000.0, ga, **kw)
                failures += self._compare(entry, losses, result, lambda t: int(np.flatnonzero(ga == t)[0]))
            else:
                for t in losses:
                    if _point(entry, t)["status"] != "ok":
                        continue
                    result = f(mat, entry["electron_energy_eV"], [t], **kw)
                    failures += self._compare(entry, [t], result, lambda _: 0)
                    compared += 1
        for set_id, loss in REQUIRED:
            self.assertEqual(_point(_set(self.paths, set_id), loss)["status"], "ok")
        self.assertEqual(failures, [])
        # A node failure leaves V8 incomplete (design section 3.2): such points are not passes.
        self.assertEqual(path_not_ok, [], f"V8 incomplete: {len(path_not_ok)} points without a valid path value")


class V9PathThree(unittest.TestCase):
    """V9: the I3-type path (k_parallel) against the implementation (routine: lambda = 20 A)."""

    @classmethod
    def setUpClass(cls):
        cls.paths = _paths()
        cls.table = _table()

    def _run(self, lams):
        f = _implementation()
        failures, uncertified = [], []
        entries = [e for e in self.paths["path3"] if e["weight_length_A"] in lams]
        self.assertEqual(len(entries), 2 * 3 * len(lams), "path-3 entries (2 materials x 3 charges per length)")
        for entry in entries:
            if entry.get("status") == "not_computed":
                uncertified.append((entry["material"], entry["weight_length_A"], "not computed"))
                continue
            mat = mp3.material_from_table(entry["material"], self.table)
            for p in entry["points"]:
                if p["status"] != "ok":
                    uncertified.append((entry["material"], entry["weight_length_A"], entry["hole_charge"],
                                        p["loss_eV"], p["status"], p.get("missing_variants")))
                    continue
                result = f(mat, 1000.0, [p["loss_eV"]], weight_length_A=entry["weight_length_A"],
                           epsilon_momentum="k_parallel", hole_charge=entry["hole_charge"], **DEF)
                for name in mp3.OUTPUTS:
                    col = p["columns"][name]
                    got = float(np.asarray(getattr(result, name))[0])
                    if not abs(got - col["value"]) <= col["width"]:
                        failures.append(f"{entry['material']} L={entry['weight_length_A']:g} Q={entry['hole_charge']:g} "
                                        f"T={p['loss_eV']!r} {name}: |Delta| {abs(got - col['value']):.3e} > "
                                        f"{col['width']:.3e}")
        self.assertEqual(uncertified, [], "path 3 not certified (one-axis changes above width/10, or not computed)")
        self.assertEqual(failures, [])

    def test_frozen_path_criterion(self):
        for entry in self.paths["path3"]:
            if entry.get("status") == "not_computed" or (not FULL and entry["weight_length_A"] != 20.0):
                continue
            for p in entry["points"]:
                if p["status"] != "ok":
                    continue
                for name, col in p["columns"].items():
                    with self.subTest(material=entry["material"], lam=entry["weight_length_A"],
                                      q=entry["hole_charge"], loss=p["loss_eV"], column=name):
                        width = mp3.V9_REL * abs(col["value"]) + mp3.V9_ABS * entry["scale_max_abs_K_total"]
                        self.assertEqual(width, col["width"])
                        self.assertEqual(set(col["changes"]), {"A", "B", "C", "D", "E"})
                        self.assertLessEqual(max(abs(x) for x in col["changes"].values()), 0.1 * width)
                        self.assertLessEqual(col["tail"], 0.1 * width)

    def test_routine_lambda_20(self):
        self._run((20.0,))

    @unittest.skipUnless(FULL, "full V9 set: set XPSLOSS_M3A_FULL=1")
    def test_full(self):
        self._run((5.0, 20.0))


if __name__ == "__main__":
    unittest.main()
