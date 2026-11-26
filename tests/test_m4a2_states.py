"""CI by discovery. Set XPSLOSS_M4A2_FULL=1 for all normal rows, and
XPSLOSS_M4A2_MUTATIONS=1 to include the 26 control/mutant pairs. Optional
XPSLOSS_M4A2_ROWS selects comma-separated IDs; XPSLOSS_M4A2_LEDGER saves JSONL.
"""

from collections import Counter
import json
import os
import sys
from unittest.mock import patch
import time
import unittest

import m4a2_states as m


class ComparatorTests(unittest.TestCase):
    def test_covered_partial_underestimated_and_precision_failure(self):
        ref = dict(value="0", u_r="0.125", d_abs="1", d_rel=None)
        self.assertEqual(m.judge_numeric(0.25, 0.5, ref)["status"], "pass")
        self.assertEqual(m.judge_numeric(0.25, 0.25, ref)["status"], "unresolved")
        self.assertEqual(m.judge_numeric(0.25, 0.0, ref)["coverage"], "underestimate")
        self.assertEqual(m.judge_numeric(1.0, 2.0, ref)["status"], "fail")
        self.assertEqual(m.judge_numeric(0.25, 0.5, ref, budget=0.125)["budget"], False)
        self.assertEqual(m.judge_numeric(0.0, 0.0, dict(ref, u_r=".26"))["status"], "no_reference")
        self.assertEqual(m.judge_numeric(0.0, -1.0, ref)["status"], "fail")
        row = dict(id="synthetic_comparator", rule="seed_reference", inputs=dict(c=0))
        log = m.Ledger(row)
        keys = {
            (row["id"], "fixed_m", "satellite"): dict(
                status="parametric",
                parameter=dict(
                    L_ref=2,
                    constant=dict(value=1, u_r=0),
                    proportional=dict(value=-3, u_r=0),
                    d_abs_constant=0,
                    d_abs_proportional=0,
                    epsrel="1e-6",
                ),
            ),
            (row["id"], "inverse_interval", "length"): dict(
                status="value", source_row="synthetic", interval=["0.25", "0.5"]
            ),
        }
        source = dict(synthetic=dict(keys={"window_area.term_13a": dict(r="3", u_r="1")}))
        with (
            patch.object(m, "references", return_value=("synthetic", keys)),
            patch.object(m, "rb_sources", return_value=source),
        ):
            m.mapped(log, "fixed_m", "satellite", -5.0, 0.0, multiplier=4.0)
            record = log.items["fixed_m", "satellite"]
            self.assertEqual(
                (record["status"], record["evaluated"]["value"], record["multiplier_hex"]),
                ("pass", "-5", float(4).hex()),
            )
            for interval, expected in [
                ((0.25, 0.5), "pass"),
                ((0.5, 1.0), "unresolved"),
                ((1.0, 2.0), "fail"),
            ]:
                m.mapped(log, "inverse_interval", "length", None, interval=interval)
                self.assertEqual(log.status("length"), expected)


class RegisteredRows(unittest.TestCase):
    def test_registered(self):
        rows = m.profile()
        self.assertTrue(rows)
        self.assertEqual(len(m.ROWS), 465)
        self.assertEqual(sum(r["ci"] for r in m.ROWS.values()), 236)
        start = time.monotonic()
        group_start = start
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
                    log = m.execute(row)
                else:
                    log = m.Ledger(row)
                    for key in row["keys"]:
                        log.record("unavailable", key, "not_run", reason="time_limit")
                elapsed = time.monotonic() - then
                group = row["id"].split("-")[0]
                seconds["inherited" if row["rule"] == "m4a_ledger" else "additional"] += elapsed
                items = [dict(kind=k, key=name, **v) for (k, name), v in log.items.items()]
                counts.update((group, v["status"]) for v in items)
                record = dict(row=row["id"], seconds=elapsed, checks=items)
                if hasattr(log, "reference_sha256"):
                    record["reference_sha256"] = log.reference_sha256
                output.write(json.dumps(record, default=str) + "\n")
                output.flush()
                bad = [v for v in items if v["status"] == "fail"]
                with self.subTest(row=row["id"]):
                    self.assertFalse(bad, ascii(bad[:2]))
        print(
            f"\nM4a-2: {len(rows)} rows in {time.monotonic() - start:.2f} s; "
            f"inherited/additional row times {dict(seconds)}; checks {dict(counts)}",
            file=sys.stderr,
        )
