"""Predict or evaluate the three synthetic inputs of the frozen dispersion test."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np

from xpsloss.model import Material, NumericalConvergenceError, Oscillator
from xpsloss.xps_fixed_depth import effective_cross_section_fixed_depth

if __package__:
    from . import numerator_rule as rule
else:
    import numerator_rule as rule


RUNS = {
    "Alg075": [(15.6, 243.3598, 0.75, 1.0)],
    "Ala025": [(15.6, 243.3598, 1.5, 0.25)],
    "Al2osc": [(15.6, 121.68, 1.5, 1.0), (25.0, 312.5, 2.0, 0.0)],
}
E0 = 1000.0
LIMIT = 0.003
SELF_LIMIT = 1e-4
COMPONENTS = ("T", "e", "h")
STEP, N_LOSS, N_DEPTH = 0.2, 500, 65


def candidate(oscs, losses, depth, panels=400):
    return np.array([rule.combine(rule.terms(t, depth, "opt", mom="kt", panels=panels,
                                             oscillators=oscs)) for t in losses])


def published(oscs, losses, depth):
    material = Material("synthetic", 0.0, tuple(Oscillator(*o) for o in oscs))
    result = effective_cross_section_fixed_depth(material, E0, losses, depth, epsilon_momentum="k_total")
    return np.column_stack((result.total, result.no_hole, result.hole))


def predictions():
    """Render the nominal predictions as the registered ASCII/CRLF byte stream."""
    lines = []
    losses = np.array([11.0, 16.0, 30.0, 60.0])
    for name, oscs in RUNS.items():
        lines.append(f"{name}: a = 10.0 A, columns T e h, eV^-1 A^-1")
        for t, c, p in zip(losses, candidate(oscs, losses, 10.0), published(oscs, losses, 10.0)):
            lines.append(f"  {t:5.1f} eV  candidate {c[0]: .4e} {c[1]: .4e} {c[2]: .4e}"
                         f"   published {p[0]: .4e} {p[1]: .4e} {p[2]: .4e}")
    return ("\r\n".join(lines) + "\r\n").encode("ascii")


def read_m(path):
    """Read the normal-emission depth table with its implicit, uniformly spaced loss axis."""
    try:
        lines = path.read_text(encoding="ascii").splitlines()
        start = lines.index("@") + 1
        energy, angle, imfp, step, count = (float(lines[start + i].split()[0]) for i in range(5))
        if (energy, angle, step, count) != (E0, 0.0, STEP, N_LOSS):
            raise ValueError("unexpected header energy, angle, loss step or loss count")
        if not np.isfinite(imfp) or imfp <= 0:
            raise ValueError("invalid header IMFP")
        # The export records ac/bc. A minimal numeric header can declare the
        # column count instead (the synthetic table format of test_compare_queels).
        extra = lines[start + 5].split()
        if len(extra) == 1:
            if int(extra[0]) != N_DEPTH:
                raise ValueError("unexpected header depth count")
        elif tuple(map(float, extra[:2])) != (1.6, -3.6):
            raise ValueError("unexpected header ac/bc")
        table = np.loadtxt(lines[start + 6:], comments=None, ndmin=2)
        if table.shape != (3 + N_LOSS, N_DEPTH) or not np.isfinite(table).all():
            raise ValueError("invalid depth/loss table shape or non-finite entry")
        depths = table[0]
        if not (depths[0] > 0 and np.all(np.diff(depths) > 0)):
            raise ValueError("expected positive increasing depths")
        if not np.array_equal(depths, table[1]):
            raise ValueError("path lengths disagree with depths at normal emission")
        if np.any(table[-1] != 0):
            raise ValueError("expected final uncomputed loss row to be zero")
    except (ValueError, IndexError) as exc:
        raise ValueError(f"{path.name}: invalid header or table ({exc})") from exc
    return depths, table[3:]


def read_run(directory, name):
    curves = []
    for component in COMPONENTS:
        path = directory / f"K_XPS_{component}_{name}_1000_0_q1_n2000.m"
        depths, table = read_m(path)
        if curves and not np.array_equal(depths, grid):
            raise ValueError(f"{path.name}: component depth grids disagree")
        grid = depths
        curves.append(table)
    return grid, np.stack(curves, axis=-1)


def metrics(model, reference, scale):
    errors = np.max(np.abs(model - reference), axis=0) / scale
    if not np.isfinite(errors).all():
        raise ValueError("non-finite normalized difference")
    return dict(zip(COMPONENTS, map(float, errors)))


def evaluate(directory):
    inputs = {name: read_run(directory, name) for name in RUNS}
    losses = STEP * np.arange(1, N_LOSS + 1)
    window = (losses >= 5.0 - 1e-9) & (losses <= 99.8 + 1e-9)
    selected = losses[window]
    runs = {}
    for name, (depths, table) in inputs.items():
        nearest = int(np.argmin(np.abs(depths - 10.0)))
        rows = {}
        # Preserve the actual diagnostic columns of the frozen evaluator.
        for role, col in (("criterion", nearest), ("diagnostic_shallow", 1), ("diagnostic_deep", 64)):
            depth = float(depths[col])
            ref = table[window, col]
            scale = float(np.max(ref[:, 0]))
            if scale <= 0:
                raise ValueError(f"{name} {role}: maximum reference total must be positive")
            coarse = candidate(RUNS[name], selected, depth)
            mc = metrics(coarse, ref, scale)
            mp = metrics(published(RUNS[name], selected, depth), ref, scale)
            row = dict(column=col, depth_A=depth, candidate=mc, published=mp)
            if role == "criterion":
                fine = candidate(RUNS[name], selected, depth, panels=800)
                selfcheck = max(metrics(fine, coarse, scale).values())
                checks = dict(candidate_pass=max(mc.values()) <= LIMIT,
                              discriminating=max(mp.values()) > LIMIT,
                              selfcheck_pass=selfcheck <= SELF_LIMIT)
                row.update(checks, selfcheck=selfcheck, passed=all(checks.values()))
            rows[role] = row
        runs[name] = rows
    return dict(confirmed=all(rows["criterion"]["passed"] for rows in runs.values()),
                loss_window_eV=[5.0, 99.8], points=len(selected), runs=runs)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("predict")
    commands.add_parser("evaluate").add_argument("reference_dir", type=Path)
    args = parser.parse_args(argv)
    if args.command == "predict":
        sys.stdout.buffer.write(predictions())
        return 0
    try:
        result = evaluate(args.reference_dir)
    except (OSError, ValueError, NumericalConvergenceError) as exc:
        print(json.dumps({"confirmed": None, "error": str(exc)}))
        return 2
    print(json.dumps(result, allow_nan=False))
    return 0 if result["confirmed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
