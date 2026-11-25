"""Run the registered M4a rows once in the execution order and write a JSON record.

    python benchmarks/numerics/m4a/run_expected_states_m4a.py --out RECORD.json [--ci] [--rows RA1,RA2]

Without --ci all 100 rows run (the V rows with the 499-point profile).  The
record has, per row, the observed state, the time, every judged key of the
state test (state, identity, exact_zero_or_null, replay, diagnostic, record)
and of the R-A comparison, and the reasons of the failures.  PYTHONPATH must
point to the src directory of the version under test.
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
    parser.add_argument("--ci", action="store_true", help="only the 72 rows of the routine run")
    parser.add_argument("--rows", help="comma separated row ids (debugging; the counts are then partial)")
    args = parser.parse_args(argv)
    if args.ci:
        os.environ.pop("XPSLOSS_M4A_FULL", None)
    else:
        os.environ["XPSLOSS_M4A_FULL"] = "1"
    if args.rows:
        os.environ["XPSLOSS_M4A_ROWS"] = args.rows
    import numpy
    import scipy
    import m4a_states as m
    try:
        import xpsloss
        package = dict(version=xpsloss.__version__, file=str(Path(xpsloss.__file__).resolve()))
    except Exception as exc:
        package = dict(error=f"{type(exc).__name__}: {exc}")
    ids = m.selected_ids()
    references = m.ra_references()["rows"]
    judgement = m.Judgement()
    start = time.perf_counter()
    record_rows, reference_counts = [], {}
    for rid in ids:
        row = m.rows()[rid]
        out = m.judge_row(row, m.run_row, judgement)
        keys = {f"{key}|{kind}": dict(status=status, detail=m._canonical(detail) if not isinstance(detail, str) else detail)
                for (r, key, kind), (status, detail) in judgement.items.items() if r == rid}
        item = dict(id=rid, family=row["family"], entry=row["entry"], ci=row["ci"], denominator=row["denominator"],
                    expected=row["expect"].get("status"), observed=out.summary(), seconds=out.seconds, keys=keys,
                    not_passed=sorted(k for k, v in keys.items() if v["status"] not in ("pass", "recorded")))
        if rid in references or rid == "N27":
            verdicts = m.judge_n27(out) if rid == "N27" else m.judge_reference_row(row, out)
            item["reference_ra"] = {k: dict(status=s, detail=m._canonical(d) if not isinstance(d, str) else d)
                                    for k, (s, d) in verdicts.items()}
            for key, (status, _) in verdicts.items():
                if not key.startswith("@"):
                    reference_counts[status] = reference_counts.get(status, 0) + 1
        if "matched_display" in out.extras:
            item["in_row_matched_display"] = m._canonical(list(out.extras["matched_display"]))
        item["row_passed"] = (not item["not_passed"] and out.extras.get("matched_display", ("pass",))[0] == "pass"
                              and all(v["status"] == "pass" for v in item.get("reference_ra", {}).values()))
        record_rows.append(item)
        print(f"{rid:6s} {str(out.status):14s} {out.seconds:8.2f} s  not passed: {len(item['not_passed'])}", flush=True)
    by_denominator = {}
    for item in record_rows:
        box = by_denominator.setdefault(item["denominator"], dict(rows=0, state_passed=0))
        box["rows"] += 1
        box["state_passed"] += item["keys"].get("terminal_state|state", {}).get("status") in ("pass", "recorded")
    v = m._CACHE.get("v", {})
    record = dict(
        schema="xpsloss-m4a-expected-states-run-v1", registration_table_sha256=m.table_sha256(),
        references_ra_sha256=m.ra_sha256(), profile="ci" if args.ci else "full", partial=bool(args.rows),
        environment=dict(python=sys.version, platform=platform.platform(), numpy=numpy.__version__,
                         scipy=scipy.__version__, package=package),
        elapsed_seconds=time.perf_counter() - start, rows_run=len(ids),
        planned_keys_by_kind=m.planned_counts(ids), judged_keys_by_kind=judgement.counts(),
        reference_ra_keys=reference_counts, rows_by_denominator=by_denominator,
        v_rows=dict(profile=v.get("profile"), seconds=v.get("seconds"), unavailable=v.get("unavailable"),
                    points={k: v[k].get("points") for k in ("tag", "current") if isinstance(v.get(k), dict)},
                    auxiliary_controls=m.auxiliary_controls() if any(r.startswith("V") for r in ids) else None),
        rows=record_rows)
    data = (json.dumps(record, indent=1, ensure_ascii=True, default=str) + "\n").encode("ascii")
    Path(args.out).write_bytes(data)
    print("wrote", args.out, "sha256", hashlib.sha256(data).hexdigest(),
          "rows passed", sum(r["row_passed"] for r in record_rows), "of", len(record_rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
