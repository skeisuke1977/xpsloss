"""M4b-1, R-A: the analytic rows against the frozen analytic references (prereg 6.1, 6.3, 7.2, 7.3).

The rows run through the surface of prereg 4 only.  The accuracy (7.2) and the
error coverage (7.3) read the same reference record and are counted separately.
An accuracy verdict is pass, fail, unresolved, no_reference or not_run; only
pass counts as agreement.  A coverage verdict is covered, underestimate,
unresolved or not_run.  A key that cannot be read is a failure, never a skip.

The routine run compares the CI rows; XPSLOSS_M4B1_FULL=1 adds the full-only
row A5.  The rows fail until the implementation of M4b-1 exists.
"""
from fractions import Fraction as Fr
import re
import sys
import unittest

import m4b1_states as m

RA_SHA256 = "036701007a3d3634b0561a911c1730d54c3ffa81c0d4e3b03f03f09c96ae107a"
GENERATOR = m.DIRECTORY/"generate_references_ra_m4b1.py"
RA_KEYS = 3572
RA_KEYS_BY_KIND = {"direct_reference": 3521, "bound_reference": 38, "error_part": 6, "interval_contains": 7}
WRONG_ANSWERS = ("X01", "X02", "X04", "X05", "X06", "X07", "X08", "X09", "X10", "X12", "X14", "X15", "X16",
                 "X18", "X19", "X20", "X23", "X33")
INV_NOTE = "reference: reaching-side bin outside the width without the tail term"


def planned():
    out = {}
    for rid, records in m.ledger().items():
        for key, kind, oracle in records:
            if kind in m.EXTERNAL and oracle.startswith(("R-A", "R-INV")):
                out.setdefault(rid, []).append((key, kind, oracle))
    return out


def reference_rows():
    ra = m.ra_references()["rows"]
    return [rid for rid in m.selected_ids() if rid in ra and ra[rid]["keys"]]


class ReferenceFileTests(unittest.TestCase):
    def test_hashes(self):
        self.assertEqual(m.table_sha256(), m.TABLE_SHA256)
        self.assertEqual(m.ra_sha256(), RA_SHA256)
        ra = m.ra_references()
        self.assertEqual(ra["registration_table_sha256"], m.TABLE_SHA256)
        self.assertEqual(ra["registration_table_script_sha256"], m.TABLE_PY_SHA256)
        raw = m.RA_PATH.read_bytes()
        self.assertTrue(raw.isascii() and b"\r" not in raw and raw.endswith(b"\n"))

    def test_generator_does_not_import_the_package(self):
        source = GENERATOR.read_text(encoding="ascii")
        self.assertIsNone(re.search(r"^\s*(import|from)\s+xpsloss", source, re.M))
        self.assertNotIn("__import__", source)
        self.assertNotIn("import_module", source)
        self.assertIs(m.ra_references()["imports_package_under_test"], False)

    def test_keys_match_the_ledger(self):
        ra, plan = m.ra_references(), planned()
        self.assertEqual(set(plan), {rid for rid, row in ra["rows"].items() if row["keys"]})
        by_kind, total = {}, 0
        for rid, records in plan.items():
            have = sorted((k, v["kind"], v["oracle"]) for k, v in ra["rows"][rid]["keys"].items())
            self.assertEqual(have, sorted(records), rid)
            self.assertEqual(ra["counts"]["by_row"][rid]["total"], len(records), rid)
            for _, kind, _ in records:
                by_kind[kind] = by_kind.get(kind, 0) + 1
            total += len(records)
        self.assertEqual(total, RA_KEYS)
        self.assertEqual(by_kind, RA_KEYS_BY_KIND)
        counts = m.table()["counts"]
        rb = sum(1 for records in m.ledger().values() for _, kind, oracle in records
                 if kind in m.EXTERNAL and oracle.startswith("R-B"))
        self.assertEqual(total + rb, counts["planned_external_reference_keys"])
        for kind, n in RA_KEYS_BY_KIND.items():
            others = sum(1 for records in m.ledger().values() for _, k, oracle in records
                         if k == kind and oracle.startswith("R-B"))
            self.assertEqual(n + others, counts["planned_keys_by_kind"][kind], kind)

    def test_uncertainties_and_no_reference_rule(self):
        """Prereg 6.3 applied before any implementation value: u_r <= width/4, or no_reference."""
        ra = m.ra_references()
        listed = {(x["row"], x["key"]) for x in ra["no_reference"]}
        for rid, row in ra["rows"].items():
            for key, ref in row["keys"].items():
                r, u = Fr(ref["value"]), Fr(ref["u_r"])
                self.assertEqual(ref["zero_reference"], r == 0, (rid, key))
                if r == 0:
                    self.assertEqual(u, 0, (rid, key))
                else:
                    digits = ref["value"].split("e")[0].replace(".", "").lstrip("-")
                    self.assertGreaterEqual(len(digits), 30, (rid, key))
                self.assertEqual(ref["status"] == "no_reference", (rid, key) in listed, (rid, key))
                if ref["width"] is not None:
                    width = Fr(ref["width"]["d_abs"]) + Fr(ref["width"]["d_rel"])*abs(r)
                    self.assertEqual(ref["status"] == "reference", 4*u <= width, (rid, key))
                else:
                    self.assertEqual(ref["status"], "reference", (rid, key))
                    self.assertIn("rule", ref, (rid, key))
        self.assertEqual(ra["rows_short_of_references"], [])

    def test_frozen_checks_of_the_reference(self):
        ra = m.ra_references()
        self.assertTrue(ra["self_checks"] and all(ra["self_checks"].values()))
        mutations = {w["mutation"] for w in ra["known_wrong_answers"]}
        self.assertTrue(set(WRONG_ANSWERS) <= mutations, sorted(set(WRONG_ANSWERS) - mutations))
        registered = {x["id"]: x for x in m.table()["mutations"]}
        hit = {}
        for wrong in ra["known_wrong_answers"]:
            self.assertTrue(wrong["outside_width"], wrong)
            targets = {(d["row"], r["key"]) for d in registered[wrong["mutation"]]["detected_by"]
                       for r in d.get("records", ())}
            # a wrong answer is on a registered key, or on a further bin of the same registered row (X29)
            self.assertIn(wrong["row"], {row for row, _ in targets}, wrong)
            hit[wrong["mutation"]] = hit.get(wrong["mutation"], False) or (wrong["row"], wrong["key"]) in targets
        self.assertTrue(all(hit.values()), hit)
        extra = [w["key"] for w in ra["known_wrong_answers"] if (w["row"], w["key"]) not in
                 {(d["row"], r["key"]) for d in registered[w["mutation"]]["detected_by"] for r in d.get("records", ())}]
        self.assertEqual(sorted(extra), [f"truncation.bound.bin[{i}].grid_minimum_with_error_terms" for i in (1, 2, 3)])
        for rid in ("INV1", "INV2"):
            note = ra["rows"][rid]["notes"][INV_NOTE]
            self.assertTrue(note["passes"] and note["outside_width_without_tail"] and note["normal_controls_exactly_equal"])
        self.assertIn("bin[9]", ra["rows"]["INV1"]["notes"][INV_NOTE]["reaching_side"] +
                      ra["rows"]["INV1"]["notes"][INV_NOTE]["normal_controls_exactly_equal"])


class ReferenceComparisonTests(unittest.TestCase):
    def test_rows_against_the_references(self):
        accuracy, coverage, judged, short = {}, {}, 0, []
        for rid in reference_rows():
            row = m.rows()[rid]
            keys = m.ra_references()["rows"][rid]["keys"]
            with self.subTest(row=rid):
                out = m.run_row(rid)
                verdicts, covers = m.judge_reference_row(row, out)
                self.assertEqual(set(verdicts), set(keys))
                by_status = {}
                for status, _ in verdicts.values():
                    accuracy[status] = accuracy.get(status, 0) + 1
                    by_status[status] = by_status.get(status, 0) + 1
                    judged += 1
                for status in covers.values():
                    coverage[status] = coverage.get(status, 0) + 1
                missing = by_status.get("no_reference", 0)
                if 10*missing > len(keys) or (missing and not any(v["kind"] == "direct_reference" for v in keys.values())):
                    short.append(rid)
                budgets = m.judge_budgets(row, out)
                bad = {k: v for k, v in verdicts.items() if v[0] != "pass"}
                uncovered = {k: v for k, v in covers.items() if v != "covered"}
                if bad or uncovered or (budgets and budgets[0] != "pass"):
                    first = next(iter(bad.items()), None) or next(iter(uncovered.items()), None) or budgets
                    self.fail(f"{rid}: accuracy {by_status} of {len(keys)} keys; coverage not covered: {len(uncovered)}; "
                              f"local budgets: {budgets}; first: {first}")
        print(f"\nM4b-1 R-A comparison: judged {judged} keys; accuracy {accuracy}; coverage {coverage}; "
              f"rows short of references: {short}", file=sys.stderr)
        if not m.os.environ.get("XPSLOSS_M4B1_ROWS"):
            a5 = len(m.ra_references()["rows"]["A5"]["keys"])
            self.assertEqual(judged, RA_KEYS if m.full_run() else RA_KEYS - a5)
            self.assertEqual(accuracy.get("pass", 0), judged, f"only pass counts as agreement: {accuracy}")
            self.assertEqual(short, [])


if __name__ == "__main__":
    unittest.main()
