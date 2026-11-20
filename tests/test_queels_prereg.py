"""Preregistered decisions on synthetic .m files and the frozen nominal predictions."""

from contextlib import redirect_stdout
import hashlib
from io import StringIO
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from benchmarks.queels_prereg import evaluate_prereg as prereg


class NominalPredictionTests(unittest.TestCase):
    def test_frozen_hash_and_recalculated_bytes(self):
        script = Path(prereg.__file__)
        frozen = script.with_name("predictions_a10.txt").read_bytes()
        self.assertEqual(hashlib.sha256(frozen).hexdigest(),
                         "715f4f455227beacce8880f5697a071a531b3f07fd7d93b4eacc8fd5163f2a3d")
        actual = subprocess.check_output([sys.executable, str(script), "predict"])
        self.assertEqual(actual, frozen)

    def test_binary_file_output_and_check_status(self):
        frozen = Path(prereg.__file__).with_name("predictions_a10.txt").read_bytes()
        with tempfile.TemporaryDirectory() as tmp, patch.object(prereg, "predictions", return_value=frozen):
            path = Path(tmp) / "result.txt"
            self.assertEqual(prereg.main(["predict", "--out", str(path)]), 0)
            self.assertEqual(path.read_bytes(), frozen)
        for sample, matches in ((frozen, True), (frozen.replace(b"\r\n", b"\n"), False)):
            with patch.object(prereg, "predictions", return_value=sample), redirect_stdout(StringIO()) as log:
                self.assertEqual(prereg.main(["predict", "--check"]), 0 if matches else 1)
                self.assertEqual(json.loads(log.getvalue()), {
                    "identical": matches, "sha256": hashlib.sha256(sample).hexdigest()})


class DecisionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temp.cleanup)
        cls.folder = Path(cls.temp.name)
        cls.expected_columns = {}
        for index, name in enumerate(("Alg075", "Ala025", "Al2osc")):
            depths = (1 + np.arange(65)) * (1 + 0.2 * index)
            column = int(np.argmin(abs(depths - 10)))
            cls.expected_columns[name] = column
            for letter, level in zip("Teh", [10000, 11000, -1000]):
                table = np.full((503, 65), level, dtype=float)
                table[:2] = depths
                table[2] = table[-1] = 0
                table[3:27] = 5e5  # Through 4.8 eV, outside the window.
                if letter == "h":
                    table[27, column] += 30  # Exactly at the criterion, at 5.0 eV.
                    table[501, 1] += 1000  # Diagnostic-only discrepancy at 99.8 eV.
                np.savetxt(cls.folder / f"K_XPS_{letter}_{name}_1000_0_q1_n2000.m",
                           table, fmt="%.10g", comments="",
                           header="synthetic\n@\n1000\n0\n20\n0.2\n500\n1.6 -3.6 ac and bc")

    def verdict(self, scenario):
        # Use deterministic synthetic curves for threshold tests. Independently
        # exercise both real solvers and all three oscillator sets in predict.
        def curves(oscs, losses, depth, panels=400):
            self.assertEqual(len(losses), 475)
            self.assertAlmostEqual(losses[0], 5.0)
            self.assertAlmostEqual(losses[-1], 99.8)
            self.assertIn(panels, (400, 800))
            values = np.tile([10000., 11000., -1000.], (len(losses), 1))
            if oscs == prereg.RUNS["Al2osc"]:
                if scenario == "candidate":
                    values[-1, 1] += 31
                if panels == 800:
                    values[-1, 2] += 2 if scenario == "selfcheck" else 1
            return values

        def published(oscs, losses, depth):
            values = curves(oscs, losses, depth)
            values[-1, 0] += 30 if scenario == "discrimination" and len(oscs) == 2 else 31
            return values

        stdout = StringIO()
        with patch.object(prereg, "candidate", curves), patch.object(prereg, "published", published):
            with redirect_stdout(stdout):
                code = prereg.main(["evaluate", str(self.folder)])
        result = json.loads(stdout.getvalue())
        self.assertEqual(code, 0 if scenario == "pass" else 1)
        self.assertIs(result["confirmed"], scenario == "pass")
        for name, rows in result["runs"].items():
            self.assertEqual(rows["criterion"]["column"], self.expected_columns[name])
            self.assertEqual(rows["diagnostic_shallow"]["column"], 1)
            self.assertEqual(rows["diagnostic_deep"]["column"], 64)
            self.assertAlmostEqual(rows["diagnostic_shallow"]["candidate"]["h"], 0.1)
        return result["runs"]["Al2osc"]["criterion"]

    def test_candidate_pass_at_inclusive_limits(self):
        row = self.verdict("pass")
        self.assertEqual(row["candidate"]["h"], 0.003)
        self.assertEqual(row["selfcheck"], 0.0001)
        self.assertTrue(row["passed"])

    def test_one_failed_candidate_component_at_last_loss(self):
        row = self.verdict("candidate")
        self.assertFalse(row["candidate_pass"])
        self.assertTrue(row["discriminating"] and row["selfcheck_pass"])

    def test_published_expression_must_strictly_exceed_limit(self):
        row = self.verdict("discrimination")
        self.assertFalse(row["discriminating"])
        self.assertTrue(row["candidate_pass"] and row["selfcheck_pass"])

    def test_selfcheck_can_fail_alone(self):
        row = self.verdict("selfcheck")
        self.assertFalse(row["selfcheck_pass"])
        self.assertTrue(row["candidate_pass"] and row["discriminating"])

    def test_missing_file_and_malformed_header_return_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "K_XPS_T_Alg075_1000_0_q1_n2000.m"
            for missing in (True, False):
                if not missing:
                    path.write_text("synthetic\n@\n", encoding="ascii")
                stdout = StringIO()
                with redirect_stdout(stdout):
                    code = prereg.main(["evaluate", tmp])
                result = json.loads(stdout.getvalue())
                self.assertEqual(code, 2)
                self.assertIsNone(result["confirmed"])
                self.assertIn(path.name, result["error"])

    def test_rejects_incompatible_depth_and_loss_grids(self):
        name = "Alg075"
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            for source in self.folder.glob(f"*_{name}_*.m"):
                (folder / source.name).write_bytes(source.read_bytes())
            target, = folder.glob("K_XPS_e*")
            lines = target.read_text(encoding="ascii").splitlines()
            for case, message in (("step", "loss step"), ("count", "loss count"),
                                  ("short", "shape"), ("nonfinite", "non-finite"),
                                  ("order", "increasing"), ("grid", "grids disagree")):
                with self.subTest(case=case):
                    changed = lines.copy()
                    if case == "step":
                        changed[5] = "0.3"
                    elif case == "count":
                        changed[6] = "499"
                    elif case == "short":
                        changed.pop()
                    elif case == "nonfinite":
                        changed[-2] = "nan " * 65
                    else:
                        grid = changed[8].split()
                        grid[0] = grid[1] if case == "order" else "0.75"
                        changed[8] = changed[9] = " ".join(grid)
                    target.write_text("\n".join(changed), encoding="ascii")
                    with self.assertRaisesRegex(ValueError, message):
                        prereg.read_run(folder, name)

    def test_header_path_lengths_and_last_row(self):
        source = next(self.folder.glob("K_XPS_T*"))
        original = source.read_text(encoding="ascii").splitlines()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / source.name
            for index, replacement, error in (
                (2, "999", "energy"), (3, "1", "angle"), (4, "nan", "IMFP"), (4, "0", "IMFP"),
                (7, "1.5 -3.6 ac and bc", "ac/bc"), (7, "1.6 -3.5 ac and bc", "ac/bc"),
                (7, "64", "depth count"), (9, "2 " * 65, "path lengths"), (-1, "1 " * 65, "final"),
            ):
                with self.subTest(index=index, replacement=replacement):
                    lines = original.copy()
                    lines[index] = replacement
                    path.write_text("\n".join(lines), encoding="ascii")
                    with self.assertRaisesRegex(ValueError, error):
                        prereg.read_m(path)
            original[7] = "65"  # Also support the earlier, count-only synthetic header.
            path.write_text("\n".join(original), encoding="ascii")
            depths, _ = prereg.read_m(path)
            self.assertEqual(len(depths), 65)

    def test_nonpositive_normalization_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            for p in self.folder.glob("*.m"):
                (folder / p.name).write_bytes(p.read_bytes())
            target, = folder.glob("K_XPS_T_Alg075*")
            lines = target.read_text(encoding="ascii").splitlines()
            for level in ("0", "-1"):
                lines[11:-1] = [" ".join([level] * 65)] * 499
                target.write_text("\n".join(lines), encoding="ascii")
                with redirect_stdout(StringIO()) as log:
                    self.assertEqual(prereg.main(["evaluate", tmp]), 2)
                    self.assertIn("maximum reference total must be positive", json.loads(log.getvalue())["error"])
