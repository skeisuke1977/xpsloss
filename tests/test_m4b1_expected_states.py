"""M4b-1: the 109 registered rows and their expected states (prereg 5; table expected_states_m4b1.json).

The routine (CI) run executes the 90 rows with "ci": true in the registered
execution order; XPSLOSS_M4B1_FULL=1 executes all 109 (the 19 full-only rows
run only on request).  This test owns the keys of the kinds state, identity,
exact_zero_or_null, bitwise, replay, diagnostic, recomputed_from_report and
inequality.  A key that cannot be read is a failure, never a skip; the planned
and the judged numbers are compared for every kind.  The keys that read an
external reference belong to test_m4b1_reference_ra.py and to the R-B test.

The rows fail until the implementation of M4b-1 exists.
"""
import sys
import time
import unittest

import m4b1_states as m

PLANNED_BY_KIND = {"state": 109, "direct_reference": 9608, "bound_reference": 38, "error_part": 6,
                   "interval_contains": 16, "identity": 5408, "exact_zero_or_null": 341, "bitwise": 37, "replay": 3,
                   "diagnostic": 79, "recomputed_from_report": 40, "inequality": 8}


class TableTests(unittest.TestCase):
    def test_table_hash_and_shape(self):
        table = m.table()
        self.assertEqual(m.table_sha256(), m.TABLE_SHA256)
        self.assertEqual(table["schema"], "xpsloss-m4b1-expected-states-v1-draft4")
        rows = table["rows"]
        self.assertEqual(len(rows), 109)
        self.assertEqual(sorted(table["execution_order"]), sorted(r["id"] for r in rows))
        self.assertEqual(sum(r["ci"] for r in rows), 90)
        self.assertEqual(tuple(table["columns"]), m.COLUMNS)
        self.assertEqual(m.POST_FREEZE_EXPECT, {})

    def test_ledger_counts(self):
        table = m.table()
        counts = {}
        for rid, records in m.ledger().items():
            self.assertEqual(len(set(records)), len(records), rid)
            for _, kind, _ in records:
                counts[kind] = counts.get(kind, 0) + 1
        self.assertEqual(counts, table["counts"]["planned_keys_by_kind"])
        self.assertEqual(counts, PLANNED_BY_KIND)
        own = m.planned_counts([r["id"] for r in table["rows"]])
        self.assertEqual(own, {k: PLANNED_BY_KIND[k] for k in m.OWN_KINDS})
        for row in table["rows"]:
            self.assertEqual(sum(p["count"] for p in row["check_plan"]), len(m.ledger()[row["id"]]), row["id"])

    def test_reuse_follows_the_execution_order(self):
        table = m.table()
        position = {rid: i for i, rid in enumerate(table["execution_order"])}
        for row in table["rows"]:
            for item in row["reuse"]:
                source, artefact = item.split(".")
                self.assertLess(position[source], position[row["id"]])
                self.assertIn(artefact, m.rows()[source]["provides"])
                if row["ci"]:
                    self.assertTrue(m.rows()[source]["ci"], (row["id"], source))

    def test_analytic_kernels_of_the_table(self):
        """The kernels handed to the entries are the registered ones (independent of the implementation)."""
        import numpy as np
        k = m.AnalyticKernel("K_A3", 2.0)
        values, errors = k(np.array([1.75, 2.0, 2.0, 5.5]), ["interior", "left", "right", "interior"])
        total = values[m.COLUMNS.index("total")].tolist()
        self.assertEqual(total, [0.03 + 2.0*(-0.05), 0.03 + 2.0*(-0.05), (0.04 + 0.03) + 2.0*(-0.05), 0.0 + 2.0*0.02])
        self.assertEqual(values[m.COLUMNS.index("term_13a")].tolist(), [0.0, 0.0, 0.04, 0.0])
        self.assertFalse(errors.any())
        self.assertEqual(k.calls, 1)
        with self.assertRaises(m.KernelContractError):
            k(np.array([2.0]), "interior")
        piece, = m.table()["kernels"]["K_A4_9"]["columns"]["term_13a"]
        self.assertEqual(m.legendre_value(piece, 5.0), 0.1)
        self.assertEqual(m.legendre_value(piece, 3.5), 0.05)


class ExpectedStateTests(unittest.TestCase):
    def test_rows_in_execution_order(self):
        ids = m.selected_ids()
        judgement = m.Judgement()
        start = time.perf_counter()
        seconds, skipped = {}, []
        for rid in ids:
            row = m.rows()[rid]
            with self.subTest(row=rid):
                out = m.judge_row(row, judgement)
                seconds[rid] = out.seconds
                bad = judgement.failures(rid)
                if any(v[0] == "skipped" for v in bad.values()):
                    skipped.append(rid)
                bad = {k: v for k, v in bad.items() if v[0] != "skipped"}
                if bad:
                    (_, key, kind), (status, detail) = next(iter(bad.items()))
                    self.fail(f"{rid}: {len(bad)} of {sum(k[0] == rid for k in judgement.items)} keys not passed; "
                              f"first: {key} [{kind}] {status}: {detail}")
        elapsed = time.perf_counter() - start
        profile = "full" if m.full_run() else "ci"
        base = m._CACHE.get("baseline", {})
        print(f"\nM4b-1 expected states: profile {profile}, rows {len(ids)}, elapsed {elapsed:.2f} s "
              f"(rows {sum(seconds.values()):.2f} s; baseline child {base.get('seconds')}); "
              f"within the 60 s CI limit: {elapsed <= 60.0 if profile == 'ci' else 'not applicable'}", file=sys.stderr)
        slow = sorted(seconds.items(), key=lambda item: -item[1])[:8]
        print("M4b-1 slowest rows: " + ", ".join(f"{rid} {t:.2f} s" for rid, t in slow), file=sys.stderr)
        print(f"M4b-1 judged keys by kind: {judgement.counts()}", file=sys.stderr)
        if skipped:
            print(f"M4b-1 rows not compared with the baseline tag ({base.get('unavailable')}): {skipped}", file=sys.stderr)
        # every planned key of the kinds of this test has a judgement; a missing one is never a pass
        planned = m.planned_counts(ids)
        judged = {kind: sum(c.values()) for kind, c in judgement.counts().items()}
        self.assertEqual(judged, planned)
        if not m.os.environ.get("XPSLOSS_M4B1_ROWS"):
            self.assertEqual(len(ids), 109 if m.full_run() else 90)
            self.assertEqual(planned["state"], len(ids))
            if m.full_run():
                self.assertEqual(planned, {k: PLANNED_BY_KIND[k] for k in m.OWN_KINDS})


if __name__ == "__main__":
    unittest.main()
