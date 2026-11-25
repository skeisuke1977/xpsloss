"""M4a: the 100 registered rows and their expected states (prereg 5; table expected_states_m4a.json).

The routine (CI) run executes the 72 rows with ``ci: true`` in the registered
execution order; ``XPSLOSS_M4A_FULL=1`` executes all 100.  This test owns the
keys of the kinds state, identity, exact_zero_or_null, replay, diagnostic and
record.  A key that cannot be read is a failure, never a skip; the planned and
the judged numbers are compared for every kind.
"""
import sys
import time
import unittest

import m4a_states as m


class TableTests(unittest.TestCase):
    def test_table_hash_and_shape(self):
        table = m.table()
        self.assertEqual(m.table_sha256(), m.TABLE_SHA256)
        self.assertEqual(table["schema"], "xpsloss-m4a-expected-states-v4")
        rows = table["rows"]
        self.assertEqual(len(rows), 100)
        self.assertEqual(sorted(table["execution_order"]), sorted(r["id"] for r in rows))
        self.assertEqual(sum(r["ci"] for r in rows), 72)
        self.assertEqual(tuple(table["columns"]), m.COLUMNS)

    def test_ledger_counts(self):
        table = m.table()
        counts = {}
        for row in table["rows"]:
            for key, kind, oracle in m.plan(row):
                counts[kind] = counts.get(kind, 0) + 1
        self.assertEqual(counts, table["counts"]["planned_keys_by_kind"])
        self.assertEqual(counts, dict(state=100, direct_reference=1986, inverse_interval=50, identity=132,
                                      exact_zero_or_null=169, replay=414, error_part=5, diagnostic=124, record=9))
        own = m.planned_counts([r["id"] for r in table["rows"]])
        self.assertEqual(sum(own.values()), 100 + 132 + 169 + 414 + 124 + 9)

    def test_reuse_follows_the_execution_order(self):
        table = m.table()
        position = {rid: i for i, rid in enumerate(table["execution_order"])}
        for row in table["rows"]:
            for item in row["reuse"]:
                source = item.split(".")[0]
                self.assertLess(position[source], position[row["id"]])
                if row["ci"]:
                    self.assertTrue(m.rows()[source]["ci"], (row["id"], source))

    def test_tag_of_the_v_rows(self):
        state = m.v_available()
        if state is not None and m.v_unavailable_verdict(state)[0] == "skipped":
            self.skipTest(state[1])
        self.assertIsNone(state)


class ExpectedStateTests(unittest.TestCase):
    def test_rows_in_execution_order(self):
        ids = m.selected_ids()
        judgement = m.Judgement()
        start = time.perf_counter()
        seconds, skipped = {}, []
        for rid in ids:
            row = m.rows()[rid]
            with self.subTest(row=rid):
                out = m.judge_row(row, m.run_row, judgement)
                seconds[rid] = out.seconds
                bad = judgement.failures(rid)
                gone = {k: v for k, v in bad.items() if v[0] == "skipped"}
                if gone:
                    skipped.append(rid)
                bad = {k: v for k, v in bad.items() if v[0] != "skipped"}
                shown = out.extras.get("matched_display")
                if shown is not None and shown[0] != "pass":
                    bad[(rid, "comparison display against the completed table (full run)", "in_row")] = shown
                if bad:
                    (_, key, kind), (status, detail) = next(iter(bad.items()))
                    self.fail(f"{rid}: {len(bad)} of {sum(k[0] == rid for k in judgement.items)} keys not passed; "
                              f"first: {key} [{kind}] {status}: {detail}")
        elapsed = time.perf_counter() - start
        v = m._CACHE.get("v", {})
        profile = "full" if m.full_run() else "ci"
        print(f"\nM4a expected states: profile {profile}, rows {len(ids)}, elapsed {elapsed:.2f} s "
              f"(rows {sum(seconds.values()):.2f} s; V children {v.get('seconds')}); "
              f"within the 60 s CI limit: {elapsed <= 60.0 if profile == 'ci' else 'not applicable'}; "
              f"V profile {v.get('profile')}, V skipped: {skipped}", file=sys.stderr)
        slow = sorted(seconds.items(), key=lambda item: -item[1])[:8]
        print("M4a slowest rows: " + ", ".join(f"{rid} {t:.2f} s" for rid, t in slow), file=sys.stderr)
        print(f"M4a judged keys by kind: {judgement.counts()}", file=sys.stderr)
        print(f"M4a auxiliary controls of V07 (not rows): {m.auxiliary_controls() if any(r.startswith('V') for r in ids) else 'not run'}",
              file=sys.stderr)
        # every planned key of the kinds of this test has a judgement; a missing one is never a pass
        planned = m.planned_counts(ids)
        judged = {kind: sum(c.values()) for kind, c in judgement.counts().items()}
        self.assertEqual(judged, planned)
        if not m.os.environ.get("XPSLOSS_M4A_ROWS"):
            self.assertEqual(len(ids), 100 if m.full_run() else 72)
            self.assertEqual(planned["state"], len(ids))
        if skipped:
            print(f"M4a V rows not compared with the tag (reason: {v.get('unavailable')}): {skipped}", file=sys.stderr)


if __name__ == "__main__":
    unittest.main()
