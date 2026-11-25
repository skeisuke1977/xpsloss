"""Full record of the M3a expected states (docs/M3A_PREREG.md section 3): all 2041 runs.

Runs every entry of benchmarks/numerics/m3a/expected_states_m3a.json through the API or the CLI
and judges it with the functions of tests/test_m3a_expected_states.py (the routine test uses the
same functions on the routine subset).  Writes one JSON with, per run, the state, the time, the
failure reasons and, for MS runs, the maximum budget ratio, followed by the cross-run checks,
the accounting of check_plan and the per-category summary.

Usage (PYTHONPATH pointing at src):
    python benchmarks/numerics/m3a/run_expected_states_m3a.py --out RECORD.json [--processes 8] [--routine]
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import platform
import sys
import tempfile
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "tests"))
import test_m3a_expected_states as states  # noqa: E402


def _one(args):
    run_id, workdir = args
    table = states.load_table()
    run = next(r for r in table["runs"] if r["id"] == run_id)
    return states.execute_and_judge(run, table, workdir)


def main(argv=None):
    parser = argparse.ArgumentParser(description="M3a expected states: full record")
    parser.add_argument("--out", required=True)
    parser.add_argument("--processes", type=int, default=8)
    parser.add_argument("--routine", action="store_true", help="only the routine runs")
    args = parser.parse_args(argv)
    table = states.load_table()
    full = not args.routine
    problems = states.check_plan_counts(table)
    if problems:
        raise SystemExit("check_plan differs from check_plan_counts: " + "; ".join(problems))
    runs = states.selected_runs(table, full)
    start = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix="m3a_states_full_") as workdir:
        with ProcessPoolExecutor(max_workers=args.processes) as pool:
            records = {rec["id"]: rec for rec in pool.map(_one, [(r["id"], workdir) for r in runs], chunksize=4)}
    cross = states.cross_run_checks(records, table)
    plan, plan_problems = states.plan_accounting(records, table, full)
    categories = states.category_summary(records, table, full)
    import numpy
    import scipy
    import xpsloss
    doc = {
        "schema": "xpsloss-m3a-expected-states-record-v1",
        "table_sha256": states.TABLE_SHA256,
        "subset": "all" if full else "routine",
        "environment": {"python": platform.python_version(), "numpy": numpy.__version__, "scipy": scipy.__version__,
                        "xpsloss": xpsloss.__version__, "platform": platform.platform()},
        "seconds": time.perf_counter() - start,
        "categories": categories,
        "check_plan": plan,
        "check_plan_problems": plan_problems,
        "cross_run_failures": cross,
        "runs": [{"id": r["id"], "family": r["family"], "category": r["category"], "entry": r["entry"],
                  "routine": r["routine"], "outcome": r["outcome"], "seconds": r["seconds"],
                  "failures": r["failures"], "reason": r["facts"].get("reason"),
                  "max_budget_ratio": r["max_budget_ratio"]} for r in (records[x["id"]] for x in runs)],
    }
    data = (json.dumps(doc, indent=1, ensure_ascii=True, allow_nan=False, default=str) + "\n").encode("ascii")
    Path(args.out).write_bytes(data)
    failed = sum(1 for r in records.values() if r["failures"])
    print(json.dumps({"runs": len(records), "failed": failed, "cross_run_failures": len(cross),
                      "plan_problems": len(plan_problems), "sha256": hashlib.sha256(data).hexdigest(),
                      "seconds": doc["seconds"]}))
    return 0 if not (failed or cross or plan_problems) else 1


if __name__ == "__main__":
    raise SystemExit(main())
