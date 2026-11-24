"""Compare the M2-1 fixed-depth model with QUEELS-XPS outputs held locally.

The acceptance rule is the one fixed in docs/M2_FIRST_CASE.md section 4.2 before any
comparison. This script contains no reference data: it reads the user's own QUEELS-XPS
output files from a directory given on the command line (they are not distributed,
see NOTICES.md). The JSON report contains comparison metrics, verdicts, run conditions
and reference-file hashes. Exit 0 if either momentum choice passes section 4.2,
1 if both fail, and 2 for invalid inputs. It omits the reference value tables.

Usage (repository root, PYTHONPATH=src):
    python benchmarks/compare_queels_fixed_depth.py --reference-dir DIR --out RESULT.json
"""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import trapezoid

from xpsloss import NumericalConvergenceError
from xpsloss.config import load_config
from xpsloss.xps_fixed_depth import effective_cross_section_fixed_depth

# Both comparisons consume the same normal-emission QUEELS table layout.
if __package__:
    from .queels_prereg.evaluate_prereg import read_m
else:
    from queels_prereg.evaluate_prereg import read_m


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "benchmarks" / "inputs" / "al_elfdb_v1.0.bulk.json"

# Section 4.2 of docs/M2_FIRST_CASE.md (fixed 2026-11-19, before any comparison).
REFERENCE_STEM = "Al_1000_0_q1_n2000"      # run 7
UNCERTAINTY_STEM = "Al_1000_0_q1_n1000"    # run 6
TARGET_DEPTH_A = 10.0                      # nearest column of the run-7 grid
LOSS_MIN_EV, LOSS_MAX_EV, LOSS_STEP_EV = 1.0, 99.8, 0.2
K_UNCERTAINTY = 3.0
GLOBAL_FRACTION = 1e-3
N_LOSS = 500
COMPONENTS = {"total": "T", "no_hole": "e", "hole": "h"}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def loss_axis(n):
    return LOSS_STEP_EV * np.arange(1, n + 1)


def compare(directory):
    files, ref, unc = {}, {}, {}
    for name, letter in COMPONENTS.items():
        for stem, store in ((REFERENCE_STEM, ref), (UNCERTAINTY_STEM, unc)):
            path = directory / f"K_XPS_{letter}_{stem}.m"
            files[path.name] = sha256(path)
            depths, curves = read_m(path)
            if store and not np.array_equal(depths, store["total"][0]):
                raise ValueError(f"{path.name}: component depth grids disagree")
            store[name] = depths, curves

    depth_grid = ref["total"][0]
    column = int(np.argmin(np.abs(depth_grid - TARGET_DEPTH_A)))
    depth = float(depth_grid[column])
    depth_shift = float(abs(unc["total"][0][column] - depth))

    loss = loss_axis(N_LOSS)
    window = (loss >= LOSS_MIN_EV - 1e-9) & (loss <= LOSS_MAX_EV + 1e-9)
    loss_w = loss[window]
    k_ref = {name: ref[name][1][window, column] for name in COMPONENTS}
    k_unc = {name: unc[name][1][window, column] for name in COMPONENTS}

    # Reference self-consistency (section 3): T = e + h within run 7.
    t_scale = float(np.max(k_ref["total"]))
    total_area = float(trapezoid(k_ref["total"], loss_w))
    if t_scale <= 0 or total_area <= 0:
        raise ValueError("maximum and integrated reference total must both be positive")
    self_check = float(np.max(np.abs(k_ref["total"] - k_ref["no_hole"] - k_ref["hole"])) / t_scale)

    material = load_config(CONFIG).material
    result = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "rule": "docs/M2_FIRST_CASE.md section 4.2",
        "reference_files_sha256": files,
        "depth_A": depth, "depth_column": column, "depth_shift_run6_vs_run7_A": depth_shift,
        "loss_window_eV": [LOSS_MIN_EV, LOSS_MAX_EV], "points": int(window.sum()),
        "reference_self_check_max_abs_T_minus_e_minus_h_over_maxT": self_check,
        "choices": {},
    }
    for choice in ("k_parallel", "k_total"):
        ours = effective_cross_section_fixed_depth(
            material, 1000.0, loss_w, depth, theta_deg=0.0, hole_charge=1.0, epsilon_momentum=choice)
        verdicts = {}
        for name in COMPONENTS:
            mine = np.asarray(getattr(ours, name), dtype=float)
            u = np.abs(k_ref[name] - k_unc[name])
            tol = K_UNCERTAINTY * u + GLOBAL_FRACTION * t_scale
            diff = np.abs(mine - k_ref[name])
            area_ref = float(trapezoid(k_ref[name], loss_w))
            area_unc = float(trapezoid(k_unc[name], loss_w))
            area_mine = float(trapezoid(mine, loss_w))
            area_tol = K_UNCERTAINTY * abs(area_ref - area_unc) + GLOBAL_FRACTION * total_area
            if not np.isfinite(tol).all() or not np.isfinite(area_tol) or min(tol.min(), area_tol) <= 0:
                raise ValueError("comparison tolerances must be positive and finite")
            worst = int(np.argmax(diff / tol))
            verdicts[name] = {
                "pointwise_pass": bool(np.all(diff <= tol)),
                "points_failing": int(np.sum(diff > tol)),
                "max_diff_over_tolerance": float(np.max(diff / tol)),
                "worst_loss_eV": float(loss_w[worst]),
                "max_abs_diff_over_maxT": float(np.max(diff) / t_scale),
                "area_pass": bool(abs(area_mine - area_ref) <= area_tol),
                "area_diff_over_area_tolerance": float(abs(area_mine - area_ref) / area_tol),
                "area_relative_diff": abs(area_mine - area_ref) / abs(area_ref) if area_ref else None,
            }
        passed = all(v["pointwise_pass"] and v["area_pass"] for v in verdicts.values())
        result["choices"][choice] = {"pass": passed, "components": verdicts}
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open("x", encoding="utf-8") as output:
            try:
                with np.errstate(over="raise", invalid="raise", divide="raise"):
                    result = compare(args.reference_dir)
                report = json.dumps(result, indent=2, allow_nan=False)
                summary = {c: v["pass"] for c, v in result["choices"].items()}
                code = 0 if any(summary.values()) else 1
            except (OSError, ValueError, ArithmeticError, NumericalConvergenceError) as exc:
                summary = {"status": "error", "message": str(exc)}
                report = json.dumps(summary, indent=2)
                code = 2
            output.write(report + "\n")
    except OSError as exc:
        message = str(exc)
        if isinstance(exc, FileExistsError) and exc.filename == str(args.out):
            message = f"output exists: {args.out}"
        summary = {"status": "error", "message": message}
        code = 2
    print(json.dumps(summary), file=sys.stderr if code == 2 else sys.stdout)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
