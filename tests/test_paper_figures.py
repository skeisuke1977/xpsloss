"""Fast numerical and reporting checks on the independently calculated figure table."""

from collections import Counter
import csv
import unittest

from benchmarks.paper_figures import generate


class PaperFigureTests(unittest.TestCase):
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
