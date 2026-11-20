"""Fast numerical and reporting checks on the independently calculated figure table."""

from collections import Counter
import csv
import unittest

from benchmarks.paper_figures import generate


class PaperFigureTests(unittest.TestCase):
    def test_200_ev_bulk_peak_at_both_widths_and_momenta(self):
        def peak(row):
            return (row["figure"], row["panel"], row["E0_eV"], row["quantity"]) == (
                "5", "a", "200", "bulk_peak_height")

        with (generate.HERE / "readings.csv").open(encoding="utf-8") as src:
            point, = filter(peak, csv.DictReader(src))
        with (generate.HERE / "comparison.csv").open(encoding="utf-8") as src:
            saved = list(filter(peak, csv.DictReader(src)))
        self.assertEqual([r["width_eV"] for r in saved], ["3.4", "2.1"])
        for row, expected in zip(generate.compare_point(point), saved):
            for choice in ("k_parallel", "k_total"):
                self.assertAlmostEqual(float(row[choice]), float(expected[choice]), delta=1e-10)
            self.assertNotEqual(row["k_parallel"], row["k_total"])
            self.assertEqual(row["status"], "depends on width assumption")
            self.assertFalse(row["k_parallel_within"])
            self.assertEqual(row["k_total_within"], row["width_eV"] == 2.1)

    def test_summary_counts_distinct_readings(self):
        with (generate.HERE / "readings.csv").open(encoding="utf-8") as src:
            readings = {tuple(row) for row in list(csv.reader(src))[1:]}
        with (generate.HERE / "comparison.csv").open(encoding="utf-8") as src:
            reader = csv.DictReader(src)
            keys = reader.fieldnames[:reader.fieldnames.index("width_eV")]
            outcomes = {}
            for row in reader:
                key = tuple(row[k] for k in keys)
                self.assertEqual(outcomes.setdefault(key, row["status"]), row["status"])
        self.assertEqual(readings, set(outcomes))
        self.assertEqual(Counter(outcomes.values()), {
            "within reading uncertainty": 37, "outside": 10, "depends on width assumption": 37})
