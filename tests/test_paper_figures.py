"""Fast numerical and reporting checks on the independently calculated figure table."""

from collections import Counter
from contextlib import redirect_stdout
import csv
from io import StringIO
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from benchmarks.paper_figures import generate


class PaperFigureTests(unittest.TestCase):
    def test_check_roundoff_and_differences_without_regenerating_figures(self):
        csv_text = "panel,k_parallel,k_total_grid,k_parallel_margin,status\na,0.0730029433454,,0.017997057,outside\n"
        md_text = "|3 single|1000|2|total bulk peak height|2.1|0.092 ± 0.001|0.07300294 ×|0.07258735 ×|誤差外|\n"

        def write_tables(folder, table, report):
            (folder / "comparison.csv").write_text(table, encoding="utf-8")
            (folder / "RESULTS.md").write_text(report, encoding="utf-8")

        with tempfile.TemporaryDirectory() as tmp, patch.object(generate, "HERE", Path(tmp)):
            write_tables(Path(tmp), csv_text, md_text)
            for table, report, diagnostics in (
                (csv_text.replace("3454", "3455").replace("057,", "058,"),
                 md_text.replace("294 ×", "295 ×"), ()),
                (csv_text.replace("0.073", "0.074").replace(",,", ",0,"),
                 md_text.replace("0.072", "0.071"), ("row 2, k_parallel:", "row 2, k_total_grid:", "line 1, k_total:")),
                (csv_text.replace("outside", "within"), md_text.replace("×", "○"), ("row 2, status:", "line 1:")),
                (csv_text.replace("0.0730029433454", "nan"), md_text, ("row 2, k_parallel:",)),
                (csv_text.splitlines(keepends=True)[0], md_text, ("row 2:",)),
            ):
                with self.subTest(table=table, report=report), redirect_stdout(StringIO()) as log:
                    with patch.object(generate, "generate", side_effect=lambda dest: write_tables(dest, table, report)):
                        self.assertEqual(generate.main(["--check"]), 1 if diagnostics else 0)
                    for detail in diagnostics:
                        self.assertIn(detail, log.getvalue())

    def test_200_ev_bulk_peak_at_policy_and_post_hoc_widths(self):
        def peak(row):
            return (row["figure"], row["panel"], row["E0_eV"], row["quantity"]) == (
                "5", "a", "200", "bulk_peak_height")

        with (generate.HERE / "readings.csv").open(encoding="utf-8") as src:
            point, = filter(peak, csv.DictReader(src))
        with (generate.HERE / "comparison.csv").open(encoding="utf-8") as src:
            saved = list(filter(peak, csv.DictReader(src)))
        self.assertEqual([float(r["width_eV"]) for r in saved], [3.4, 2.8, 3.28, 2.1])
        for row, expected in zip(generate.compare_point(point), saved):
            for choice in ("k_parallel", "k_total"):
                self.assertAlmostEqual(float(row[choice]), float(expected[choice]), delta=1e-10)
            self.assertNotEqual(row["k_parallel"], row["k_total"])
            self.assertEqual(row["status"], "outside")
            self.assertEqual(row["policy_width"], row["width_eV"] != 2.1)
            self.assertFalse(row["k_parallel_within"])
            self.assertEqual(row["k_total_within"], row["width_eV"] == 2.1)
            self.assertGreater(float(row["k_parallel"]), float(row["k_parallel_grid"]))

    def test_table_widths_and_interpolation(self):
        for e, width in ((175, 3.4), (300, 2.8), (500, 2.1), (1180, 1.5), (2000, 1.3)):
            self.assertEqual(generate.widths(e), [(width, "Table I")])
        for energy, low, high, expected in ((200, 3.4, 2.8, 3.28), (1000, 2.1, 1.5, 1.658823529411765),
                                             (1413.5, 1.5, 1.3, 1.443048780487805)):
            values = [w for w, _ in generate.widths(energy)]
            self.assertEqual(values[:2], [low, high])
            self.assertAlmostEqual(values[2], expected, places=12)
        self.assertEqual(generate.widths(5000), [(1.3, "assumed from 2000 eV")])

    def test_summary_counts_distinct_readings(self):
        with (generate.HERE / "readings.csv").open(encoding="utf-8") as src:
            readings = {tuple(row) for row in list(csv.reader(src))[1:]}
        with (generate.HERE / "comparison.csv").open(encoding="utf-8") as src:
            reader = csv.DictReader(src)
            keys = reader.fieldnames[:reader.fieldnames.index("width_eV")]
            outcomes = {}
            for row in reader:
                key = tuple(row[k] for k in keys)
                statuses = tuple(row[k] for k in ("status", "k_parallel_status", "k_total_status"))
                self.assertEqual(outcomes.setdefault(key, statuses), statuses)
        self.assertEqual(readings, set(outcomes))
        for statuses, expected in zip(zip(*outcomes.values()), ((36, 18, 30), (37, 19, 28), (37, 18, 29))):
            self.assertEqual(Counter(statuses), dict(zip((generate.WITHIN, generate.OUTSIDE, generate.DEPENDENT),
                                                        expected)))
