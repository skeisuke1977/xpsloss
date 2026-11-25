"""M4a, R-A: the analytic rows against the frozen analytic references (prereg 6.1, 7.2, 7.3).

The rows run through the test ports of prereg 4.7 only.  A key is pass, fail,
unresolved, no_reference or not_run; only pass counts as agreement.
"""
import hashlib
import json
import re
import unittest

import m4a_states as m

RA_SHA256 = "d4b9939bc5a623ed794462a58584b4dd2b7e324452e9486f1d6f10f4e12d4c96"
GENERATOR = m.ROOT/"benchmarks/numerics/m4a/generate_references_ra_m4a.py"


def reference_rows():
    return [rid for rid in m.table()["execution_order"] if rid in m.ra_references()["rows"]]


class ReferenceFileTests(unittest.TestCase):
    def test_hashes(self):
        self.assertEqual(m.table_sha256(), m.TABLE_SHA256)
        self.assertEqual(m.ra_sha256(), RA_SHA256)
        self.assertEqual(m.ra_references()["registration_table_sha256"], m.TABLE_SHA256)
        raw = m.RA_PATH.read_bytes()
        self.assertTrue(raw.isascii() and b"\r" not in raw)

    def test_generator_does_not_import_the_package(self):
        source = GENERATOR.read_text(encoding="ascii")
        self.assertIsNone(re.search(r"^\s*(import|from)\s+xpsloss", source, re.M))
        self.assertNotIn("__import__", source)

    def test_keys_match_the_ledger(self):
        ra, table = m.ra_references(), m.table()
        planned = {}
        for row in table["rows"]:
            for key, kind, oracle in m.plan(row):
                if kind in ("direct_reference", "inverse_interval", "error_part") and oracle.startswith("R-A"):
                    planned.setdefault(row["id"], []).append((key, kind))
        self.assertEqual(set(planned), set(ra["rows"]))
        numeric = 0
        for rid, keys in planned.items():
            self.assertEqual([(k, v["kind"]) for k, v in sorted(ra["rows"][rid]["keys"].items())], sorted(keys), rid)
            count = sum(kind != "error_part" for _, kind in keys)
            self.assertEqual(count, table["counts"]["planned_numeric_keys_by_row"].get(rid, 0), rid)
            numeric += count
        self.assertEqual(numeric, table["counts"]["planned_numeric_keys_by_family"]["RA"])
        self.assertEqual(numeric, 922)
        self.assertEqual(sum(kind == "error_part" for keys in planned.values() for _, kind in keys), 4)

    def test_uncertainties_and_no_reference_rule(self):
        """Prereg 6.3 applied before any implementation value: u_r <= width/4, or no_reference."""
        from fractions import Fraction as Fr
        ra = m.ra_references()
        listed = {(x["row"], x["key"]) for x in ra["no_reference"]}
        for rid, row in ra["rows"].items():
            for key, ref in row["keys"].items():
                if ref["kind"] == "inverse_interval":
                    self.assertEqual(ref["reference_key"], "window_area.total")
                    continue
                r, u = Fr(ref["value"]), Fr(ref["u_r"])
                width = Fr(ref["d_abs"]) + Fr(ref["d_rel"])*abs(r) if r else Fr(ref["d_abs"])
                self.assertEqual(ref["zero_reference"], r == 0, (rid, key))
                self.assertGreaterEqual(len(ref["value"].split("e")[0].replace(".", "").lstrip("-")), 30 if r else 1)
                self.assertEqual(ref["status"] == "reference", 4*u <= width, (rid, key))
                self.assertEqual(ref["status"] == "no_reference", (rid, key) in listed)
                if r == 0:
                    self.assertEqual(u, 0, (rid, key))

    def test_frozen_checks_of_the_reference(self):
        ra = m.ra_references()
        check = ra["ra13_prefreeze_check"]
        self.assertEqual(check["initial_panels"], 6)
        self.assertTrue(check["normal"]["passes"])
        self.assertFalse(check["truncated_after_initial_panels"]["passes"])
        mutations = {w["mutation"] for w in ra["known_wrong_answers"]}
        self.assertTrue({"M01", "M02a", "M03", "M05b", "M05c", "M07", "M10", "M11a", "M12", "M13"} <= mutations)
        for wrong in ra["known_wrong_answers"]:
            self.assertTrue(wrong["outside_width"], wrong)
        registered = {x["id"]: x for x in m.table()["mutations"]}
        for wrong in ra["known_wrong_answers"]:
            targets = {(d["row"], k) for d in registered[wrong["mutation"]]["detected_by"] for k in d["keys"]}
            key = wrong["key"].replace(".kernel_point_error", "") if wrong["row"] == "RA9" else wrong["key"]
            hit = (wrong["row"], key) in targets or (wrong["row"], wrong["key"]) in targets
            self.assertTrue(hit or wrong["row"] in {r for r, _ in targets}, wrong)


class ReferenceComparisonTests(unittest.TestCase):
    def test_rows_against_the_references(self):
        counts, problems, judged = {}, {}, 0
        selected = set(m.selected_ids()) if m.os.environ.get("XPSLOSS_M4A_ROWS") else None
        for rid in reference_rows():
            if selected is not None and rid not in selected:
                continue
            row = m.rows()[rid]
            with self.subTest(row=rid):
                out = m.run_row(rid)
                verdicts = m.judge_reference_row(row, out)
                keys = m.ra_references()["rows"][rid]["keys"]
                self.assertEqual({k for k in verdicts if not k.startswith("@")}, set(keys))
                bad = {}
                for key, (status, detail) in verdicts.items():
                    if not key.startswith("@"):
                        counts[status] = counts.get(status, 0) + 1
                        judged += 1
                    if status != "pass":
                        bad[key] = (status, detail)
                if bad:
                    problems[rid] = bad
                    by_status = {}
                    for status, _ in bad.values():
                        by_status[status] = by_status.get(status, 0) + 1
                    first = next(iter(bad.items()))
                    self.fail(f"{rid}: {by_status} of {len(keys)} keys; first: {first[0]} -> {first[1]}")
        if selected is None or "N27" in selected:
            with self.subTest(row="N27"):
                (key, (status, detail)), = m.judge_n27(m.run_row("N27")).items()
                counts[status] = counts.get(status, 0) + 1
                judged += 1
                self.assertEqual(status, "pass", f"N27 {key}: {detail}")
        print(f"\nM4a R-A comparison: judged {judged} keys, {counts}", file=m.sys.stderr)
        if selected is None:
            # 922 numeric keys (direct_reference, inverse_interval), the 4 error_part keys of R-A
            # and the analytic error_part key of N27
            self.assertEqual(judged, 927)
            self.assertEqual(counts.get("pass", 0), 927, f"only pass counts as agreement: {counts}")


if __name__ == "__main__":
    unittest.main()
