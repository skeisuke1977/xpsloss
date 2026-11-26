"""CI by discovery. Set XPSLOSS_M4A2_FULL=1 for all normal rows, and
XPSLOSS_M4A2_MUTATIONS=1 to include the 26 control/mutant pairs. Optional
XPSLOSS_M4A2_ROWS selects comma-separated IDs; XPSLOSS_M4A2_LEDGER saves JSONL.
XPSLOSS_M4A2_NO_TIMING=1 omits timing measurements (the registered deadline still applies).
"""

from collections import Counter
from copy import deepcopy
from fractions import Fraction
import json
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from unittest.mock import patch
import time
import unittest

import m4a2_states as m


class ComparatorTests(unittest.TestCase):
    def test_covered_partial_underestimated_and_precision_failure(self):
        ref = dict(value="0", u_r="0.125", d_abs="1", d_abs_exact="1", d_rel=None)
        self.assertEqual(m.judge_numeric(0.25, 0.5, ref)["status"], "pass")
        self.assertEqual(m.judge_numeric(0.25, 0.25, ref)["status"], "unresolved")
        self.assertEqual(m.judge_numeric(0.25, 0.0, ref)["coverage"], "underestimate")
        self.assertEqual(m.judge_numeric(1.0, 2.0, ref)["status"], "fail")
        self.assertEqual(m.judge_numeric(0.25, 0.5, ref, budget=0.125)["budget"], False)
        self.assertEqual(m.judge_numeric(0.0, 0.0, dict(ref, u_r=".26"))["status"], "no_reference")
        self.assertEqual(m.judge_numeric(0.0, -1.0, ref)["status"], "fail")

    def test_fixed_terms_rounding_exact_floor_and_eligibility(self):
        row = dict(id="synthetic")
        rec = dict(status="parametric", counted_as="derived_density",
                   parameter=dict(symbol="m_hat", reference_length_A="2"),
                   terms=dict(constant=dict(value="-3", u_r="0"), proportional=dict(value="2", u_r="0"),
                              divide_by="3", d_abs=dict(constant_exact="1", per_m_hat_exact="1/4"),
                              epsrel_exact="1/1000000"))
        pk = (row["id"], "fixed_m", "density")
        table = m.ReferenceTable("synthetic", {pk: rec})
        with patch.object(m, "references", return_value=table):
            log = m.Ledger(row)
            m.mapped(log, *pk[1:], 1/3, 1e-15, multiplier=4.0)
            checked = log.items[pk[1:]]
            self.assertEqual((checked["status"], checked["multiplier_hex"]), ("pass", (4.).hex()))
            evaluated = checked["evaluated"]
            self.assertEqual(Fraction(evaluated["value"]), Fraction(int("3"*40), 10**40))
            self.assertEqual(Fraction(evaluated["u_r"]), Fraction(int("3"*39 + "4"), 10**80))
            self.assertEqual(Fraction(evaluated["d_abs"]), Fraction(int("6"*40), 10**40))
            self.assertEqual(Fraction(evaluated["d_rel"]), Fraction("1e-6"))
            self.assertEqual(Fraction(checked["d_abs_exact"]), Fraction(2, 3))
            with self.assertRaisesRegex(ValueError, "repeated comparison"):
                m.mapped(log, *pk[1:], 1/3, 1e-15, multiplier=4.0)

            # 8u/|r| equals the limit exactly. Decimal rounding widens the stored interval.
            rec["terms"]["constant"]["u_r"] = "1/80000"
            log = m.Ledger(row)
            m.mapped(log, *pk[1:], 1/3, 1e-15, multiplier=4.0)
            self.assertEqual(log.status("density"), "no_reference")
            self.assertEqual(Fraction(log.items[pk[1:]]["evaluated"]["d_rel"]), Fraction("1e-4"))
            rec["terms"]["constant"] = dict(value="-4", u_r="1")
            result = m.fixed_reference(rec, 4.0)
            self.assertEqual((result["status"], result["value"], result["d_rel"]), ("no_reference", "0", None))

    def test_exact_interval_and_empty_intersection(self):
        pk = ("synthetic", "inverse_interval", "common_interval")
        rec = dict(status="value", counted_as="inverse_interval", interval=["0.25", "0.5"],
                   interval_exact=["1/4", "1/2"])
        with patch.object(m, "references", return_value=m.ReferenceTable("synthetic", {pk: rec})):
            def compare(bounds):
                log = m.Ledger(dict(id=pk[0]))
                m.mapped(log, *pk[1:], None, interval=bounds)
                return log.items[pk[1:]]

            for bounds, expected in [((0.25, 0.5), "pass"), ((0.5, 1), "unresolved"),
                                     ((1, 2), "fail"), ((0.4, 0.3), "fail")]:
                self.assertEqual(compare(bounds)["status"], expected)
            # A displayed endpoint touches P, while the exact interval is disjoint.
            rec["interval_exact"][0] = str(Fraction(1, 4) + Fraction(1, 10**42))
            self.assertEqual(compare((0, 0.25))["status"], "fail")
            rec["interval"][0] = "0.3"
            self.assertEqual(compare((0, 1))["reason"], "mapped_interval_rounding")

    def test_absent_reference_file(self):
        with TemporaryDirectory() as tmp, patch.object(m, "DIRECTORY", Path(tmp)), \
                patch.object(m, "references", m.references.__wrapped__):
            log = m.Ledger(dict(id="synthetic"))
            m.mapped(log, "fixed_m", "density", 1.0, 0.0, multiplier=4.0)
            self.assertEqual(log.items["fixed_m", "density"], dict(status="not_run", reason="references_absent"))

    def test_frozen_hashes_and_reference_inventory(self):
        rid, kind, key = pk = ("V4-RB1", "ratio", "single_loss_weight")
        inputs = {"benchmarks/numerics/m4a/" + name: sha for name, sha in m.TABLE["frozen_inputs"].items()}
        inputs["benchmarks/numerics/m4a2/expected_states_m4a2.json"] = m.STATES_SHA256
        plan = {kind: {"direct_reference": 1}}
        doc = dict(metadata=dict(inputs=inputs, check_kinds={kind: "4.3"},
                                 planned_keys=dict(by_row_and_kind={rid: plan}, by_kind=plan)),
                   records=[dict(row_id=rid, check_kind=kind, key=key, counted_as="identity_4.4(5)",
                                 source_row="RB1", source_sha256=m.TABLE["frozen_inputs"]["references_rb_m4a.json"])])
        table = m.ReferenceTable("synthetic", m.index_references(doc))
        table.audit_usage({rid}, Counter({pk: 1}))
        for used in [Counter(), Counter({pk: 2}), Counter({pk: 1, (rid, kind, "extra"): 1})]:
            with self.assertRaisesRegex(ValueError, "reference usage"):
                table.audit_usage({rid}, used)
        for field, value, message in [("source_sha256", "bad", "source hash"),
                                      ("key", "surprise", "unregistered"), ("counted_as", "wrong", "identity")]:
            changed = deepcopy(doc)
            changed["records"][0][field] = value
            with self.assertRaisesRegex(ValueError, message):
                m.index_references(changed)
        changed = deepcopy(doc)
        changed["metadata"]["planned_keys"]["by_kind"][kind]["direct_reference"] = 2
        with self.assertRaisesRegex(ValueError, "planned_keys"):
            m.index_references(changed)
        with patch.object(Path, "read_bytes", return_value=b"changed"), self.assertRaisesRegex(ValueError, "input hash"):
            m.index_references(doc)
        del inputs["benchmarks/numerics/m4a2/expected_states_m4a2.json"]
        with self.assertRaisesRegex(ValueError, "input hashes"):
            m.index_references(doc)


def profile_failures(row, items, *, strict, references_present):
    bad = [v for v in items if v["status"] == "fail"]
    if strict:
        bad.extend(v for v in items if v["status"] == "not_run"
                   and (v.get("reason") == "unmapped_registered_key"
                        or (references_present and v.get("reason") != "time_limit")))
        if row["rule"] == "mutation":
            for key in row["keys"]:
                detected = any(v["kind"] == "mutation" and v["key"] == key and v["status"] == "pass"
                               and v.get("control") == "pass" and v.get("mutant") == "fail" for v in items)
                if not detected:
                    bad.append(dict(key=key, reason="mutation_not_detected"))
    return bad


class RegistrationTests(unittest.TestCase):
    def test_missing_evidence_cannot_be_satisfied_by_an_execution_marker(self):
        log = m.Ledger(dict(id="synthetic", rule="analytic_map", keys=["H", "Hpm", "unseen"]))
        log.check("cutoff_sensitivity.total", True)
        log.record("execution", "unseen", "pass")
        mapping = m.registered_checks(log)
        self.assertEqual(mapping, {"H": ["identity/cutoff_sensitivity.total"], "Hpm": [], "unseen": []})
        self.assertEqual(log.items["registration", "Hpm"], dict(status="not_run", reason="unmapped_registered_key"))

    def test_profile_failures_preserve_three_way_results_and_require_mutation_detection(self):
        row = dict(rule="external_map", keys=["k"])
        def failures(status, reason=None, strict=True, present=True, **witness):
            item = dict(key="k", kind=row["rule"], status=status, reason=reason, **witness)
            return profile_failures(row, [item], strict=strict, references_present=present)

        for reason in ("references_absent", "execution_failed", "unmapped_registered_key"):
            self.assertTrue(failures("not_run", reason))
            self.assertFalse(failures("not_run", reason, strict=False))
        self.assertTrue(failures("not_run", "unmapped_registered_key", present=False))
        self.assertFalse(failures("not_run", "references_absent", present=False))
        self.assertFalse(failures("not_run", "time_limit"))
        for status in ("no_reference", "unresolved"):
            self.assertFalse(failures(status))
        row["rule"] = "mutation"
        for status in ("not_run", "unresolved", "no_reference", "fail"):
            self.assertTrue(failures(status, control=status, mutant="fail"))
        self.assertTrue(failures("pass", control="pass", mutant="pass"))
        self.assertFalse(failures("pass", control="pass", mutant="fail"))


class RegisteredRows(unittest.TestCase):
    def test_registered(self):
        rows = m.profile()
        self.assertTrue(rows)
        self.assertEqual(len(m.ROWS), 465)
        self.assertEqual(sum(r["ci"] for r in m.ROWS.values()), 236)
        table = m.references()  # Verify the frozen inputs before any comparison.
        mutations = os.environ.get("XPSLOSS_M4A2_MUTATIONS") == "1"
        strict = mutations or os.environ.get("XPSLOSS_M4A2_FULL") == "1"
        if strict:
            # The table is committed, so the full and mutation profiles require it.
            self.assertIsNotNone(table, "references_m4a2.json is required in the full and mutation profiles")
        if mutations and not os.environ.get("XPSLOSS_M4A2_ROWS"):
            self.assertCountEqual([r["id"] for r in rows if r["rule"] == "mutation"],
                                  [f"V11-M{i:02d}" for i in range(1, 27)])
        undecided = []
        used = Counter()
        seed_statuses = {}
        timing = os.environ.get("XPSLOSS_M4A2_NO_TIMING") != "1"
        start = group_start = time.monotonic()
        in_mutations = False
        counts = Counter()
        seconds = Counter()
        ledger_path = os.environ.get("XPSLOSS_M4A2_LEDGER")
        # The inherited V child selects a denser grid in the full profile.
        env = {"XPSLOSS_M4A_FULL": "1"} if os.environ.get("XPSLOSS_M4A2_FULL") == "1" else {}
        with patch.dict(os.environ, env), open(ledger_path or os.devnull, "w", encoding="utf-8") as output:
            for row in rows:
                if row["rule"] == "mutation" and not in_mutations:
                    in_mutations = True
                    group_start = time.monotonic()
                then = time.monotonic()
                if then - group_start < 7200:
                    log = m.execute(row, seed_statuses)
                else:
                    log = m.Ledger(row)
                    for key in row["keys"]:
                        log.record("unavailable", key, "not_run", reason="time_limit")
                mapping = m.registered_checks(log) if strict else None
                if row["rule"] == "seed_reference":
                    seed_statuses[row["id"]] = log.status()
                group = row["id"].split("-")[0]
                items = [dict(kind=k, key=name, **v) for (k, name), v in log.items.items()]
                used.update((row["id"], v["kind"], v["key"]) for v in items if "counted_as" in v)
                counts.update((group, v["status"]) for v in items)
                record = dict(row=row["id"], checks=items)
                if strict:
                    record["registered_checks"] = mapping
                    undecided.extend((row["id"], v["kind"], v["key"], v["status"]) for v in items
                                     if v["status"] in ("unresolved", "no_reference"))
                if timing:
                    record["seconds"] = time.monotonic() - then
                    seconds["inherited" if row["rule"] == "m4a_ledger" else "additional"] += record["seconds"]
                if hasattr(log, "reference_sha256"):
                    record["reference_sha256"] = log.reference_sha256
                output.write(json.dumps(record, default=str) + "\n")
                output.flush()
                bad = profile_failures(row, items, strict=strict, references_present=table is not None)
                with self.subTest(row=row["id"]):
                    self.assertFalse(bad, ascii(bad[:2]))
        if strict:
            print("M4a-2 unresolved/no_reference (row, kind, key, status): " + json.dumps(undecided), file=sys.stderr)
        if table is not None:
            table.audit_usage({r["id"] for r in rows}, used)
        timing_info = f" in {time.monotonic() - start:.2f} s; row times {dict(seconds)}" if timing else ""
        print(f"\nM4a-2: {len(rows)} rows{timing_info}; checks {dict(counts)}", file=sys.stderr)
