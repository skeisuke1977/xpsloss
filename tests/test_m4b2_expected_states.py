"""M4b-2: the 82 registered rows and their expected states (prereg 5, 7; table expected_states_m4b2.json).

The routine (CI) run executes the 70 rows with "ci": true in the registered execution order;
XPSLOSS_M4B2_FULL=1 executes all 82 (the 12 full-only rows run only on request).  Every ledger
record (row, key, kind, oracle) of the frozen table gets one judgement; a key that cannot be
read is a failure, never a skip.  Counted apart and never as a pass: ``not_run`` (for example
``sdist_without_git_history`` for the baseline comparisons of an unpacked sdist, or a physical
control that did not converge) and ``no_reference``.  An unused or twice-used reference key
fails its row.

The BE rows fail until the implementation of M4b-2 exists.
"""
import os
import sys
import time
import unittest

import m4b2_states as m

PLANNED_BY_KIND = {"bin_membership": 347, "bit_identity": 23243, "bound_bracket": 9, "diagnostic": 144,
                   "exact_conversion": 599, "inequality": 1, "record_identity": 469, "state": 268, "truncation": 868}


class TableTests(unittest.TestCase):
    def test_table_hash_and_shape(self):
        table = m.table()
        self.assertEqual(m.table_sha256(), m.TABLE_SHA256)
        self.assertEqual(table["schema"], "xpsloss-m4b2-expected-states-v1-draft4")
        rows = table["rows"]
        self.assertEqual(len(rows), 82)
        self.assertEqual(sorted(table["execution_order"]), sorted(r["id"] for r in rows))
        self.assertEqual(sum(r["ci"] for r in rows), 70)
        self.assertEqual(tuple(table["key_names"]["columns"]), m.COLUMNS)
        self.assertEqual(table["baseline"]["commit"], m.BASELINE_COMMIT)
        self.assertEqual(tuple(table["vocabulary"]["warning_codes_in_order"]), m.BE_CODES)
        self.assertEqual(m.POST_FREEZE_EXPECT, {})

    def test_ledger_counts(self):
        table = m.table()
        counts = {}
        for rid, records in m.ledger().items():
            self.assertEqual(len(set(records)), len(records), rid)
            for _, kind, _ in records:
                counts[kind] = counts.get(kind, 0) + 1
        self.assertEqual(counts, table["counts"]["ledger_records_by_kind"])
        self.assertEqual(counts, PLANNED_BY_KIND)
        self.assertEqual(sum(v for k, v in counts.items() if k != "record_identity"),
                         table["counts"]["ledger_fixed_schema_values"])
        for row in table["rows"]:
            plan = {}
            for key, kind, oracle in m.ledger()[row["id"]]:
                plan[(kind, oracle)] = plan.get((kind, oracle), 0) + 1
            registered = {}
            for p in row["check_plan"]:
                registered[(p["kind"], p["oracle"])] = registered.get((p["kind"], p["oracle"]), 0) + p["count"]
            self.assertEqual(plan, registered, row["id"])

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
            for _key, _kind, oracle in m.ledger()[row["id"]]:
                if oracle.endswith(".result"):
                    self.assertIn(oracle, row["reuse"], row["id"])

    def test_analytic_kernels_of_the_table(self):
        """The kernels handed to the entries are the registered ones (independent of the implementation)."""
        import numpy as np
        k = m.AnalyticKernel("K_A3", 2.0)
        values, errors = k(np.array([1.75, 2.0, 2.0, 5.5]), ["interior", "left", "right", "interior"])
        total = values[m.COLUMNS.index("total")].tolist()
        self.assertEqual(total, [0.03 + 2.0*(-0.05), 0.03 + 2.0*(-0.05), (0.04 + 0.03) + 2.0*(-0.05), 0.0 + 2.0*0.02])
        self.assertFalse(errors.any())
        with self.assertRaises(m.KernelContractError):
            k(np.array([2.0]), "interior")
        e = m.AnalyticKernel("K_T1e")
        lo = 1.0 + 2.0**-45
        self.assertEqual(m.kernel_window("K_T1e"), (lo, 2.5))
        values, _ = e(np.array([1.0, lo, lo, 2.0]), ["interior", "left", "right", "interior"])
        self.assertEqual(values[0].tolist(), [0.0, 0.0, 0.2, 0.2])

    def test_inputs_of_the_rows_can_be_built(self):
        """Plain inputs (no object of the package): axes, defects, CLI files, the not_recorded lists."""
        table = m.table()
        for row in table["rows"]:
            i = row["inputs"]
            if row["entry"] == "cli":
                data = m.spectrum_json(row)
                self.assertTrue(m.json.dumps(data).isascii(), row["id"])
            elif isinstance(i.get("energy_axis"), dict):
                axis = m.axis_plain(i["energy_axis"])
                m.apply_defects(axis, i["energy_axis"].get("defects", ()))
                self.assertEqual(set(axis) - set(table["record_fields"]),
                                 {"photon_energy_eV", "reference", "primary_reference_binding_eV",
                                  "primary_reference_choice"}, row["id"])
        nr = m.axis_plain(m.rows()["BW1"]["inputs"]["energy_axis"])
        self.assertEqual(list(m.not_recorded_names(nr)), m.rows()["BW1"]["expect"]["not_recorded"])
        nest = m.axis_plain(m.rows()["BW4"]["inputs"]["energy_axis"])
        self.assertEqual(list(m.not_recorded_names(nest)), m.rows()["BW4"]["expect"]["not_recorded"])
        self.assertEqual(m.not_recorded_names(m.axis_plain(m.rows()["BE1"]["inputs"]["energy_axis"])), ())


class ExpectedStateTests(unittest.TestCase):
    def test_rows_in_execution_order(self):
        ids = m.selected_ids()
        judgement = m.Judgement()
        start = time.perf_counter()
        not_run = {}
        for rid in ids:
            row = m.rows()[rid]
            with self.subTest(row=rid):
                m.judge_row(row, judgement)
                bad = judgement.not_passed(rid)
                for (_, key, kind, _oracle), (status, detail) in bad.items():
                    if status == "not_run" and detail == m.NO_HISTORY:
                        not_run.setdefault(rid, []).append(key)
                bad = {k: v for k, v in bad.items() if not (v[0] == "not_run" and v[1] == m.NO_HISTORY)}
                faults = judgement.faults.get(rid)
                if bad:
                    (_, key, kind, oracle), (status, detail) = next(iter(bad.items()))
                    total = sum(k[0] == rid for k in judgement.items)
                    self.fail(f"{rid}: {len(bad)} of {total} keys not passed; first: {key} [{kind}; {oracle}] "
                              f"{status}: {detail}")
                if faults:
                    self.fail(f"{rid}: reference keys {faults}")
        elapsed = time.perf_counter() - start
        profile = "full" if m.full_run() else "ci"
        base = m._CACHE.get("baseline", {})
        print(f"\nM4b-2 expected states: profile {profile}, rows {len(ids)}, elapsed {elapsed:.2f} s "
              f"(baseline child {base.get('seconds')}); "
              f"within the 60 s CI limit: {elapsed <= 60.0 if profile == 'ci' else 'not applicable'}", file=sys.stderr)
        slow = sorted(judgement.seconds.items(), key=lambda item: -item[1])[:8]
        print("M4b-2 slowest rows (run, controls and judgement): "
              + ", ".join(f"{rid} {t:.2f} s" for rid, t in slow), file=sys.stderr)
        print(f"M4b-2 judged keys by kind: {judgement.counts()}", file=sys.stderr)
        print(f"M4b-2 judged keys by status: {judgement.by_status()}; not passed by cause: {judgement.categories()}",
              file=sys.stderr)
        if not_run:
            print(f"M4b-2 keys not run and not counted as passed ({m.NO_HISTORY}): "
                  f"{ {rid: len(keys) for rid, keys in not_run.items()} }", file=sys.stderr)
        # every planned record has a judgement; a missing one is never a pass
        planned = m.planned_counts(ids)
        judged = {kind: sum(c.values()) for kind, c in judgement.counts().items()}
        self.assertEqual(judged, planned)
        if not os.environ.get("XPSLOSS_M4B2_ROWS"):
            self.assertEqual(len(ids), 82 if m.full_run() else 70)
            self.assertEqual(planned["state"] >= len(ids), True)
            if m.full_run():
                self.assertEqual(planned, PLANNED_BY_KIND)


if __name__ == "__main__":
    unittest.main()
