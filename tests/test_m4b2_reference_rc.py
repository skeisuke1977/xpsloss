"""M4b-2: the reference R-C (prereg 6) is reproducible, self-consistent and agrees with the frozen table.

Nothing here calls the implementation: the tests read the generator, the reference file and the
frozen table only.  They check the byte-for-byte regeneration (``--check``), the independence of
the generator (standard library only, no ``xpsloss``), the self checks, that every known wrong
answer lies outside the tolerance of prereg 7.2, that the reference holds exactly the ledger keys
whose oracle is R-C, and that the reference and the hand calculations inside the table generator
(two authors, two codes) give the same geometry.
"""
import ast
from fractions import Fraction as Fr
import hashlib
import subprocess
import sys
import unittest

import m4b2_states as m

RC_SHA256 = "6f285cdbc942934b0535ec119d1ff8df308362780eea6c7f7212c2549c2a4a7d"
RC_KEY_TOTAL = 1771
ALLOWED_IMPORTS = {"fractions", "hashlib", "json", "math", "pathlib", "struct", "sys"}
GEOMETRY_ROWS = 38
WRONG_IDS = {"X25(1)", "X25(2)", "X25(3)", "X26(1)", "X26(2)", "X27", "X32", "X34", "X36", "X41", "X43(1)", "X43(2)",
             "BR1-binding-overlap"}


def value(obj):
    if isinstance(obj, dict) and "hex" in obj:
        return float.fromhex(obj["hex"]).hex()
    if isinstance(obj, dict):
        return {k: value(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [value(v) for v in obj]
    return obj


class ReferenceFileTests(unittest.TestCase):
    def test_check_reproduces_the_bytes(self):
        done = subprocess.run([sys.executable, str(m.RC_PY), "--check"], capture_output=True, text=True, cwd=m.ROOT)
        self.assertEqual(done.returncode, 0, done.stderr[-500:])
        self.assertIn(RC_SHA256, done.stdout)
        raw = m.RC_PATH.read_bytes()
        self.assertEqual(hashlib.sha256(raw).hexdigest(), RC_SHA256)
        self.assertTrue(raw.isascii() and b"\r" not in raw and raw.endswith(b"\n"))

    def test_generator_is_independent_of_the_package(self):
        source = m.RC_PY.read_bytes()
        self.assertTrue(source.isascii() and b"\r" not in source)
        tree = ast.parse(source.decode("ascii"))
        imported = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported |= {a.name.split(".")[0] for a in node.names}
            elif isinstance(node, ast.ImportFrom):
                imported.add((node.module or "").split(".")[0])
            elif isinstance(node, ast.Call) and getattr(node.func, "id", "") in ("__import__", "exec", "eval"):
                self.fail("dynamic import or execution in the reference generator")
        self.assertLessEqual(imported, ALLOWED_IMPORTS)
        self.assertNotIn("xpsloss", imported)
        self.assertNotIn("numpy", imported)
        self.assertIn(m.TABLE_SHA256.encode("ascii"), source)

    def test_header_and_counts(self):
        ref = m.references()
        self.assertEqual(ref["schema"], "xpsloss-m4b2-references-rc-v1")
        self.assertEqual(ref["table"]["sha256"], m.TABLE_SHA256)
        self.assertEqual(ref["no_reference"], [])
        self.assertEqual(ref["key_count_total"], RC_KEY_TOTAL)
        self.assertEqual(ref["key_counts"], {rid: len(rec["keys"]) for rid, rec in ref["rows"].items()})
        groups = ref["row_groups"]
        self.assertEqual((len(groups["geometry"]), groups["collision"], groups["support_orders"]),
                         (GEOMETRY_ROWS, ["BN01", "BN23", "BN18", "BC06"], ["BT2b", "BS1x", "BS1dl"]))
        self.assertEqual(sorted(ref["rows"]), sorted(sum(groups.values(), [])))
        self.assertGreaterEqual(len(ref["self_checks"]), 20)

    def test_keys_are_the_ledger_keys_with_the_oracle_rc(self):
        """No registered R-C key lacks a reference and no reference key lacks a consumer."""
        ref = m.references()
        total = 0
        for rid, records in m.ledger().items():
            wanted = [key for key, _kind, oracle in records if oracle.startswith("R-C")]
            self.assertEqual(len(wanted), len(set(wanted)), rid)
            if not wanted and rid not in ref["rows"]:
                continue
            self.assertIn(rid, ref["rows"], rid)
            self.assertEqual(set(ref["rows"][rid]["keys"]), set(wanted), rid)
            total += len(wanted)
        self.assertEqual(total, RC_KEY_TOTAL)

    def test_known_wrong_answers_lie_outside_the_tolerance(self):
        ref = m.references()
        wrong = ref["known_wrong_answers"]
        self.assertEqual({w["id"] for w in wrong}, WRONG_IDS)
        for w in wrong:
            label = (w["id"], w["row"], w["key"])
            self.assertTrue(w["outside_tolerance"], label)
            if w["kind"] == "exact_conversion":
                self.assertNotEqual(value(w["wrong"]), value(w["right"]), label)
            elif w["kind"] == "bin_membership":
                y, r = Fr(w["wrong"]), Fr(w["right"])
                self.assertEqual(Fr(float(y)), y, label)
                self.assertFalse(m.within_membership(float(y), r), label)
            else:
                self.assertNotEqual(w["wrong"], w["right"], label)
            # the right value is the one of the reference (where the key is a reference key)
            keys = ref["rows"][w["row"]]["keys"]
            if w["key"] in keys and w["kind"] != "state":
                self.assertEqual(value(keys[w["key"]]), value(w["right"]), label)
        be8 = [w for w in wrong if w["id"] == "X34" and "kinetic_output_edges" in w["key"]]
        self.assertEqual(len(be8), 6)

    def test_registered_examples(self):
        rows = m.references()["rows"]
        self.assertEqual({k: rows["BR1"]["keys"][f"zero_loss.bin[{k}]"] for k in (0, 1)}, {0: "1/2", 1: "1/2"})
        diag = rows["BR1"]["zero_loss_exact_binding_overlap_diagnostic_only"]
        self.assertTrue(diag["differs_from_kinetic_widths"])
        self.assertEqual((diag["values"]["bin[0]"], diag["values"]["bin[1]"]), ("3/8", "5/8"))
        self.assertEqual([rows["BE4"]["keys"][f"zero_loss.bin[{k}]"] for k in range(3)], ["1/4", "1/2", "1/4"])
        self.assertEqual([rows["BE9"]["keys"][f"zero_loss.bin[{k}]"] for k in range(3)], ["1/2", "1/2", "0/1"])
        s = rows["BS1"]["support"]
        self.assertEqual((s["support_order_binding"], s["support_order_kinetic"], s["support_order"], s["computed"]),
                         (2, 1, 2, 1))
        self.assertEqual(rows["BS1d"]["control"]["max_order"], 2)
        self.assertEqual(rows["BZ0"]["control"]["max_order"], "display_support")
        self.assertEqual(rows["BDS"]["control"]["max_order"], 3)
        self.assertEqual(rows["BE2"]["keys"]["energy_axis.conversion_constant_exact"], "5929/4")
        self.assertEqual(rows["BC02"]["keys"]["run.json: energy_axis.implied_offset_eV"], {"hex": (4.25).hex()})
        self.assertEqual(rows["BE7"]["keys"]["energy_axis.electron_energy_minus_reference_eV"], {"hex": (1.0).hex()})
        self.assertEqual(rows["BE1"]["map"]["kinetic_bin_of_binding_bin"], [6, 5, 4, 3, 2, 1, 0])
        self.assertEqual(rows["BE1r"]["map"]["kinetic_bin_of_binding_bin"], [0, 1, 2, 3, 4, 5, 6])


class TableAgreementTests(unittest.TestCase):
    """The reference against the hand calculations of the table generator (prereg 6, independence)."""

    def test_geometry_of_every_row(self):
        ref, compared = m.references(), 0
        for rid in ref["row_groups"]["geometry"]:
            row, rec = m.rows()[rid], ref["rows"][rid]
            g, s, keys = row["geometry"], rec["support"], rec["keys"]
            cli = row["entry"] == "cli"
            pre, zero = ("run.json: ", "CSV ") if cli else ("", "")
            for name in ("support_order_binding", "support_order_kinetic", "support_order", "computed", "bound_class",
                         "region_loss_upper_eV", "width_equals_kinetic", "width_excess_over_kinetic_eV",
                         "kinetic_side_reaches_next_order"):
                self.assertEqual(value(g[name]), value(s[name]), (rid, name))
                compared += 1
            self.assertEqual(g["state_without_tolerance"], s["state_without_tolerance"], rid)
            self.assertEqual(g["support_disagreement_warning"], s["warning"], rid)
            sc = g["support_comparison"]
            for name in ("computed_relation_mismatch_count", "next_relation_mismatch_count",
                         "next_reachability_mismatch_bins"):
                self.assertEqual(sc[name], s[name], (rid, name))
            first = s["first_mismatches"] if cli else keys["energy_axis.support_comparison.first_mismatches"]
            self.assertEqual(value(sc["first_mismatches"]), value(first), rid)
            self.assertEqual(g["conversion_constant_exact"], keys[pre + "energy_axis.conversion_constant_exact"], rid)
            self.assertEqual(g["conversion_exact"], rec["conversion"]["exact"], rid)
            self.assertEqual(g["implied_offset_is_none"], keys[pre + "energy_axis.implied_offset_eV"] is None, rid)
            for region, text in g["zero_loss_kinetic_widths"].items():
                self.assertEqual(text, keys[zero + "zero_loss." + region], (rid, region))
                compared += 1
            if "zero_loss_exact_binding_overlap_diagnostic_only" in g:
                self.assertEqual(g["zero_loss_exact_binding_overlap_diagnostic_only"],
                                 rec["zero_loss_exact_binding_overlap_diagnostic_only"]["values"], rid)
            if "electron_energy_minus_reference_exact" in g:
                got = keys[pre + "energy_axis.electron_energy_minus_reference_eV"]
                self.assertEqual(Fr(float.fromhex(got["hex"])), Fr(g["electron_energy_minus_reference_exact"]), rid)
            e = row["expect"]
            self.assertEqual(e["support_order"], s["support_order"], rid)
            if "computed" in e:
                self.assertEqual(e["computed"], s["computed"], rid)
            if "truncation.state" in keys:
                self.assertEqual(e["truncation_state"], keys["truncation.state"], rid)
            compared += 12
        self.assertGreater(compared, 1000)

    def test_collisions_and_support_orders(self):
        ref = m.references()
        for rid in ref["row_groups"]["collision"]:
            e, c = m.rows()[rid]["expect"], ref["rows"][rid]["collision"]
            self.assertEqual([float.fromhex(h).hex() for h in e["binding_eV"]["hex"]], value(c["binding_eV"]), rid)
            self.assertEqual([float.fromhex(h).hex() for h in e["kinetic_eV"]["hex"]], value(c["kinetic_eV"]), rid)
            self.assertEqual(e["roles"], c["roles"], rid)
        rows = ref["rows"]
        self.assertEqual(rows["BT2b"]["support"], {"support_order_binding": 50, "support_order_kinetic": 50,
                                                   "support_order": 50})
        self.assertEqual(rows["BS1x"]["support"]["support_order"], m.rows()["BS1x"]["expect"]["required_order_estimate"])
        self.assertEqual(max(rows["BS1dl"]["support"]["support_order_binding"],
                             rows["BS1dl"]["support"]["support_order_kinetic"]),
                         m.rows()["BS1dl"]["expect"]["requested_order"])


if __name__ == "__main__":
    unittest.main()
