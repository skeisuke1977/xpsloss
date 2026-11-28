"""Run the registered M4b-2 rows once in the execution order and write a JSON record.

    python benchmarks/numerics/m4b2/run_expected_states_m4b2.py --out RECORD.json [--ci] [--rows BE1,BS1]

Without --ci all 82 rows run (the 12 full-only rows included).  The record has, per row, the
observed state, the time, and every judged ledger record (row, key, kind, oracle) with its
status: pass, fail, not_run (with the reason) or no_reference.  The counts keep these apart
(prereg 7.3); a row passes only when every record of it passed and no reference key of it was
left unused or used twice.  PYTHONPATH must point to the src directory of the version under test.
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
    parser.add_argument("--ci", action="store_true", help="only the 70 rows of the routine run")
    parser.add_argument("--rows", help="comma separated row ids (debugging; the counts are then partial)")
    args = parser.parse_args(argv)
    if args.ci:
        os.environ.pop("XPSLOSS_M4B2_FULL", None)
    else:
        os.environ["XPSLOSS_M4B2_FULL"] = "1"
    if args.rows:
        os.environ["XPSLOSS_M4B2_ROWS"] = args.rows
    import numpy
    import scipy
    import m4b2_states as m
    try:
        import xpsloss
        package = dict(version=xpsloss.__version__, file=str(Path(xpsloss.__file__).resolve()))
    except Exception as exc:
        package = dict(error=f"{type(exc).__name__}: {exc}")
    ids = m.selected_ids()
    judgement = m.Judgement()
    start = time.perf_counter()
    record_rows = []
    for rid in ids:
        row = m.rows()[rid]
        out = m.judge_row(row, judgement)
        keys = [dict(key=key, kind=kind, oracle=oracle, status=status, detail=detail)
                for (r, key, kind, oracle), (status, detail) in judgement.items.items() if r == rid]
        statuses = {}
        for k in keys:
            statuses[k["status"]] = statuses.get(k["status"], 0) + 1
        observed = dict(status=out.status, seconds=round(out.seconds, 4))
        if out.exc is not None:
            observed.update(exception=type(out.exc).__name__, reason=getattr(out.exc, "reason", None),
                            message=str(out.exc)[:400])
        if out.harness_error:
            observed["harness_error"] = out.harness_error[-800:]
        faults = judgement.faults.get(rid, {})
        item = dict(id=rid, group=row["group"], entry=row["entry"], ci=row["ci"], denominator=row["denominator"],
                    expected=row["expect"].get("status"), observed=observed, seconds=judgement.seconds[rid],
                    records=len(keys), by_status=statuses, reference_key_faults=faults,
                    not_passed=[k for k in keys if k["status"] != "pass"],
                    row_passed=statuses.get("pass", 0) == len(keys) and not faults)
        record_rows.append(item)
        print(f"{rid:16s} {str(out.status):14s} {judgement.seconds[rid]:8.2f} s  "
              f"not passed: {len(item['not_passed'])} of {len(keys)}", flush=True)
    base = m._CACHE.get("baseline", {})
    record = dict(
        schema="xpsloss-m4b2-expected-states-run-v1", registration_table_sha256=m.table_sha256(),
        references_rc_sha256=m.rc_sha256(), profile="ci" if args.ci else "full", partial=bool(args.rows),
        environment=dict(python=sys.version, platform=platform.platform(), numpy=numpy.__version__,
                         scipy=scipy.__version__, package=package),
        elapsed_seconds=time.perf_counter() - start, rows_run=len(ids),
        planned_records_by_kind=m.planned_counts(ids), judged_records_by_kind=judgement.counts(),
        judged_records_by_status=judgement.by_status(), not_passed_by_cause=judgement.categories(),
        rows_passed=sum(r["row_passed"] for r in record_rows),
        baseline=dict(commit=m.BASELINE_COMMIT, seconds=base.get("seconds"), unavailable=base.get("unavailable")),
        rows=record_rows)
    data = (json.dumps(record, indent=1, ensure_ascii=True, default=str) + "\n").encode("ascii")
    Path(args.out).write_bytes(data)
    print("wrote", args.out, "sha256", hashlib.sha256(data).hexdigest(),
          "rows passed", record["rows_passed"], "of", len(record_rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
