"""Run the registered M4b-1 rows once in the execution order and write a JSON record.

    python benchmarks/numerics/m4b1/run_expected_states_m4b1.py --out RECORD.json [--ci] [--rows A1,A1s]

Without --ci all 109 rows run (the 19 full-only rows included).  The record has,
per row, the observed state, the time, every judged key of the state test
(state, identity, exact_zero_or_null, bitwise, replay, diagnostic,
recomputed_from_report, inequality) and, for the rows with R-A references, the
accuracy verdicts (prereg 7.2) and the coverage verdicts (prereg 7.3) in
separate fields.  The keys of R-B are not judged here.  PYTHONPATH must point
to the src directory of the version under test.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/"tests"))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", required=True)
    parser.add_argument("--ci", action="store_true", help="only the 90 rows of the routine run")
    parser.add_argument("--rows", help="comma separated row ids (debugging; the counts are then partial)")
    args = parser.parse_args(argv)
    if args.ci:
        os.environ.pop("XPSLOSS_M4B1_FULL", None)
    else:
        os.environ["XPSLOSS_M4B1_FULL"] = "1"
    if args.rows:
        os.environ["XPSLOSS_M4B1_ROWS"] = args.rows
    import numpy
    import scipy
    import m4b1_states as m
    try:
        import xpsloss
        package = dict(version=xpsloss.__version__, file=str(Path(xpsloss.__file__).resolve()))
    except Exception as exc:
        package = dict(error=f"{type(exc).__name__}: {exc}")
    ids = m.selected_ids()
    references = m.ra_references()["rows"]
    judgement = m.Judgement()
    start = time.perf_counter()
    record_rows, accuracy_counts, coverage_counts, short = [], {}, {}, []
    for rid in ids:
        row = m.rows()[rid]
        out = m.judge_row(row, judgement)
        keys = {f"{key}|{kind}": dict(status=status, detail=detail)
                for (r, key, kind), (status, detail) in judgement.items.items() if r == rid}
        item = dict(id=rid, group=row["group"], entry=row["entry"], ci=row["ci"], denominator=row["denominator"],
                    expected=row["expect"].get("status"), observed=out.summary(), seconds=out.seconds, keys=keys,
                    not_passed=sorted(k for k, v in keys.items() if v["status"] != "pass"))
        reference_ok = True
        if rid in references and references[rid]["keys"]:
            accuracy, coverage = m.judge_reference_row(row, out)
            item["reference_ra_accuracy"] = {k: dict(status=s, detail=d) for k, (s, d) in accuracy.items()}
            item["reference_ra_coverage"] = coverage
            item["reference_ra_local_budgets"] = m.judge_budgets(row, out)
            by_status = {}
            for status, _ in accuracy.values():
                accuracy_counts[status] = accuracy_counts.get(status, 0) + 1
                by_status[status] = by_status.get(status, 0) + 1
            for status in coverage.values():
                coverage_counts[status] = coverage_counts.get(status, 0) + 1
            if 10*by_status.get("no_reference", 0) > len(accuracy):
                short.append(rid)
            reference_ok = (by_status.get("pass", 0) == len(accuracy) and all(v == "covered" for v in coverage.values())
                            and (item["reference_ra_local_budgets"] or ("pass",))[0] == "pass")
        item["row_passed"] = not item["not_passed"] and reference_ok
        record_rows.append(item)
        print(f"{rid:12s} {str(out.status):14s} {out.seconds:8.2f} s  not passed: {len(item['not_passed'])}", flush=True)
    by_denominator = {}
    for item in record_rows:
        box = by_denominator.setdefault(item["denominator"], dict(rows=0, state_passed=0, rows_passed=0))
        box["rows"] += 1
        box["state_passed"] += item["keys"].get("terminal_state|state", {}).get("status") == "pass"
        box["rows_passed"] += item["row_passed"]
    base = m._CACHE.get("baseline", {})
    record = dict(
        schema="xpsloss-m4b1-expected-states-run-v1", registration_table_sha256=m.table_sha256(),
        references_ra_sha256=m.ra_sha256(), profile="ci" if args.ci else "full", partial=bool(args.rows),
        environment=dict(python=sys.version, platform=platform.platform(), numpy=numpy.__version__,
                         scipy=scipy.__version__, package=package),
        elapsed_seconds=time.perf_counter() - start, rows_run=len(ids),
        planned_keys_by_kind=m.planned_counts(ids), judged_keys_by_kind=judgement.counts(),
        reference_ra_accuracy=accuracy_counts, reference_ra_coverage=coverage_counts,
        rows_short_of_references=short, rows_by_denominator=by_denominator,
        baseline=dict(tag=m.BASELINE_TAG, commit=m.BASELINE_COMMIT, seconds=base.get("seconds"),
                      unavailable=base.get("unavailable")),
        rows=record_rows)
    data = (json.dumps(record, indent=1, ensure_ascii=True, default=str) + "\n").encode("ascii")
    Path(args.out).write_bytes(data)
    print("wrote", args.out, "sha256", hashlib.sha256(data).hexdigest(),
          "rows passed", sum(r["row_passed"] for r in record_rows), "of", len(record_rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
