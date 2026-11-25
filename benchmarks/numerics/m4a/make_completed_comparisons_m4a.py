"""Write the completed-comparison table of M4a (prereg 4.6) from the frozen R-B references.

Each record is one public entry that was compared key by key with the
independent R-B references (RB1-RB8 for single_loss_spectrum, RB1a-RB4a
for loss_window_area). The conditions are the normalized dictionary the
implementation itself builds for the display, captured while the row runs,
so the display matches only these exact conditions. CE rows are not
recorded: their success is not promoted to "matched" for the general area
entry (prereg 4.6). R-A rows use the analytic layer, not a public entry.

Usage (from the repository root, PYTHONPATH=src):
    python benchmarks/numerics/m4a/make_completed_comparisons_m4a.py          (write)
    python benchmarks/numerics/m4a/make_completed_comparisons_m4a.py --check  (compare)
"""
import hashlib
import json
import pathlib
import sys
from unittest import mock

ROOT = pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests"))
OUT = ROOT / "src" / "xpsloss" / "data" / "m4a_completed_comparisons.json"
TABLE = ROOT / "benchmarks" / "numerics" / "m4a" / "expected_states_m4a.json"
REFS = {"r_a": ROOT / "benchmarks" / "numerics" / "m4a" / "references_ra_m4a.json",
        "r_b": ROOT / "benchmarks" / "numerics" / "m4a" / "references_rb_m4a.json"}
ROWS = ("RB1", "RB2", "RB3", "RB4", "RB5", "RB6", "RB7", "RB8", "RB1a", "RB2a", "RB3a", "RB4a")


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build():
    import test_m4a_reference_rb as rb
    from xpsloss import xps_single_loss as sl
    records = []
    for row_id in ROWS:
        captured = []
        original = sl._comparison

        def spy(conditions, edges, basis):
            captured.append(json.loads(json.dumps(conditions)))
            return original(conditions, edges, basis)

        with mock.patch.object(sl, "_comparison", spy):
            spec, result, _ = rb.run_row(row_id)
        assert len(captured) == 1, (row_id, len(captured))
        statuses = rb.compare_result(row_id, spec, result, rb.Tally())
        groups = result.diagnostics["comparison"]["key_status"]
        key_status = {}
        for group, keys in groups.items():
            key_status[group] = {}
            for key in keys:
                assert key in statuses, (row_id, group, key)
                s = statuses[key]
                key_status[group][key] = {"status": s, "reason": None if s == "pass" else f"R-B comparison: {s}"}
        assert set(statuses) == {k for g in key_status.values() for k in g}, row_id
        entry = "single_loss_spectrum" if "spectrum" in groups else "loss_window_area"
        records.append({"entry": entry, "row": row_id, "conditions": captured[0], "key_status": key_status})
    return {"version": "m4a-1", "registration_sha256": sha256(TABLE),
            "reference_sha256": {k: sha256(p) for k, p in REFS.items()}, "records": records}


def main(argv):
    data = (json.dumps(build(), indent=1, sort_keys=True, ensure_ascii=True) + "\n").encode("ascii")
    if "--check" in argv:
        same = OUT.read_bytes() == data
        print(json.dumps({"identical": same, "sha256": hashlib.sha256(data).hexdigest()}))
        return 0 if same else 1
    OUT.write_bytes(data)
    print(json.dumps({"written": str(OUT.relative_to(ROOT)), "sha256": hashlib.sha256(data).hexdigest()}))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
