"""Regenerate the figure comparisons from the recorded readings and independent equations."""

import argparse
from collections import Counter
import csv
from functools import lru_cache
from itertools import zip_longest
import math
from pathlib import Path
import re
import tempfile

import numpy as np
from scipy.optimize import minimize_scalar

if __package__:
    from . import model
else:
    import model

HERE = Path(__file__).resolve().parent
WITHIN = "within reading uncertainty"
OUTSIDE = "outside"
DEPENDENT = "depends on width assumption"

# compare_point writes 12 significant digits for values and 8 for margins.
CSV_REL_TOL = 1e-9
MARGIN_REL_TOL = 1e-7
ABS_TOL = 1e-12
NUMERIC_COLUMNS = {
    "figure", "E0_eV", "depth_A", "loss_eV", "reading", "uncertainty", "width_eV",
    "k_parallel", "k_parallel_grid", "k_parallel_margin", "k_total", "k_total_grid", "k_total_margin",
}
# write_results prints each calculated value to 7 significant digits next to its verdict.
RESULT_REL_TOL = 1e-6
RESULT_VALUE = re.compile(r"(?<=\|)(-?\d+(?:\.\d+)?(?:e[+-]?\d+)?)(?= [○×]\|)")


TABLE_I = {175: 3.4, 300: 2.8, 500: 2.1, 1180: 1.5, 2000: 1.3}
POST_HOC = {(3, 1000): (1.6, 1.7), (5, 200): (2.1,)}


def widths(energy):
    """Classification widths for every figure, from Table I on p. 1615."""
    if energy in TABLE_I:
        return [(TABLE_I[energy], "Table I")]
    if energy > max(TABLE_I):
        return [(TABLE_I[2000], "assumed from 2000 eV")]
    for low, high in zip(sorted(TABLE_I), sorted(TABLE_I)[1:]):
        if low < energy < high:
            interpolated = TABLE_I[low] + (TABLE_I[high] - TABLE_I[low]) * (energy - low) / (high - low)
            return [(TABLE_I[low], f"neighbour at {low} eV"),
                    (TABLE_I[high], f"neighbour at {high} eV"),
                    (interpolated, f"linear interpolation ({low}/{high} eV)")]
    raise ValueError(f"energy below Table I: {energy} eV")


@lru_cache(maxsize=None)
def components(figure, energy, depth, width, momentum, loss):
    if figure == 8:
        electron, hole = model.depth_average(loss, width, momentum)
    elif depth is None:
        return {"bulk": model.bulk(energy, loss, width)}
    else:
        panels = 48 if figure == 3 else 96
        electron, hole = model.split(energy, loss, depth, width, momentum, panels=panels)
    return {"electron": electron, "hole": hole, "total": electron + hole}


@lru_cache(maxsize=None)
def peak_values(figure, energy, depth, width, momentum, component, surface):
    """Return grid and refined (position, height) pairs for the same local maximum."""
    def at(loss):
        return components(figure, energy, depth, width, momentum, loss)[component]

    if figure == 3:
        low, high, step = (8 if component == "bulk" else 13), 25, 0.05
    elif surface:
        low, high, step = 10, 12, 0.1
    else:
        low = 14
        high = 17 if (figure == 5 and component == "electron") else 17.5
        step = 0.25 if figure == 8 else 0.1
    # The reviewers used grid maxima. Keep the original open peak-search
    # windows for Figs. 5/8 and include the upper endpoint for Fig. 3.
    first = round(low/step) + int(component != "bulk")
    grid = np.arange(first, round(high/step) + int(figure == 3)) * step
    grid = grid.round(10)
    values = [at(loss) for loss in grid]
    peak = np.argmax(values)
    if not 0 < peak < len(grid) - 1:
        raise ValueError("peak is not bracketed inside its search window")
    bounds = (grid[peak - 1], grid[peak + 1])
    optimum = minimize_scalar(lambda x: -at(x), bounds=bounds, method="bounded",
                              options={"xatol": 1e-8})
    if not optimum.success:
        raise ValueError(f"peak refinement failed: {optimum.message}")
    # Check the optimizer against a three-point vertex on a decimal grid.
    # Re-evaluate the curve itself at the reported micro-eV position.
    spacing = 0.0001
    center = round(float(optimum.x), 4)
    left, middle, right = at(center - spacing), at(center), at(center + spacing)
    curvature = left - 2 * middle + right
    if curvature >= 0 or middle < max(left, right):
        raise ValueError("local grid does not bracket a concave maximum")
    position = round(center + spacing * (left - right) / (2 * curvature), 6)
    height = at(position)
    if abs(position - optimum.x) > 1e-5 or abs(height + optimum.fun) > 1e-10:
        raise ValueError("peak refinements disagree beyond the position/height tolerances")
    return (grid[peak], values[peak]), (position, height)


def computed(point, width, momentum):
    figure = int(point["figure"])
    energy = float(point["E0_eV"])
    depth = float(point["depth_A"]) if point["depth_A"] else None
    component, quantity = point["component"], point["quantity"]
    if quantity == "value":
        value = components(figure, energy, depth, width, momentum, float(point["loss_eV"]))[component]
        return None, value
    grid, refined = peak_values(figure, energy, depth, width, None if component == "bulk" else momentum,
                                component, quantity.startswith("surface"))
    index = 0 if quantity.endswith("position") else 1
    return grid[index], refined[index]


def classify(policy_rows, choices, extrapolated):
    outcomes = [{row[ch + "_within"] for row in policy_rows} for ch in choices]
    if any(len(flags) > 1 for flags in outcomes):
        return DEPENDENT
    if any(False in flags for flags in outcomes):
        return OUTSIDE
    return DEPENDENT if extrapolated else WITHIN


def compare_point(point):
    observed, error = float(point["reading"]), float(point["uncertainty"])
    rows = []
    energy = float(point["E0_eV"])
    widths_to_check = widths(energy) + [(w, "sensitivity (post hoc)")
                                        for w in POST_HOC.get((int(point["figure"]), energy), ())]
    for width, source in widths_to_check:
        row = dict(point, width_eV=width, width_source=source, policy_width=source != "sensitivity (post hoc)")
        for choice in model.MOMENTA:
            grid, value = computed(point, width, choice)
            row[choice] = f"{value:.12g}"
            row[choice + "_grid"] = "" if grid is None else f"{grid:.12g}"
            margin = abs(value - observed) - error
            row[choice + "_margin"] = f"{margin:.8g}"
            row[choice + "_within"] = bool(margin <= 0)
        rows.append(row)
    policy_rows = [row for row in rows if row["policy_width"]]
    extrapolated = energy > max(TABLE_I)
    statuses = {ch + "_status": classify(policy_rows, (ch,), extrapolated) for ch in model.MOMENTA}
    statuses["status"] = classify(policy_rows, model.MOMENTA, extrapolated)
    for row in rows:
        row.update(statuses)
    return rows


def write_results(groups, destination):
    lines = [
        "# 原著の図3・5・8との比較",
        "",
        "原図は A. Cohen Simonsen, F. Yubero & S. Tougaard, Phys. Rev. B 56, 1612 (1997) の"
        "図3（p. 1616）、図5(a,b)（p. 1617）、図8（p. 1618）。[原著](https://doi.org/10.1103/PhysRevB.56.1612)",
        "",
        "独立実装で式(20)を使い、法線出射の比較点ごとに紙面の読取りと照合した。"
        "同じ点の線幅違いは一度だけ数える。",
        "",
        "| 図 | 判定に使う運動量 | 比較点 | 読取り誤差内 | 誤差外 | 幅の仮定に依存 |",
        "|---|---|---:|---:|---:|---:|",
    ]
    points = [rows[0] for rows in groups]
    for figure in ("全図", "3", "5", "8"):
        selected = [p for p in points if figure == "全図" or p["figure"] == figure]
        for label, field in (("両方", "status"), ("k_parallel", "k_parallel_status"), ("k_total", "k_total_status")):
            counts = Counter(p[field] for p in selected)
            lines.append(f"| {figure} | {label} | {len(selected)} | {counts[WITHIN]} |"
                         f" {counts[OUTSIDE]} | {counts[DEPENDENT]} |")
    labels = {WITHIN: "誤差内", OUTSIDE: "誤差外", DEPENDENT: "幅依存"}
    lines.extend(["", "片方の運動量だけが全判定幅で誤差外となった点を示す。", ""])
    for p in points:
        failed = [ch for ch in model.MOMENTA if p[ch + "_status"] == OUTSIDE]
        if len(failed) == 1:
            name = p["quantity"].replace("_", " ") if not p["loss_eV"] else f'{p["loss_eV"]} eV の値'
            lines.append(f'- 図{p["figure"]}({p["panel"]})、E₀ = {p["E0_eV"]} eV、'
                         f'{p["component"]} {name}。誤差外は {failed[0]}。総合判定は{labels[p["status"]]}。')
    lines.extend([
        "",
        "図5(b)の2000 eV付近では上位の正孔曲線が重なり、読み取ったピークの帰属に不確かさがある。",
        "",
        "各行に一つの幅を示す。○はその幅で読取り誤差内、×は誤差外。「事後」の行は集計から除いた。"
        "全体判定の規則、読取り方、単位、再計算コマンドは [README](README.md) に記載した。"
        "数値と各運動量の判定は [比較表](comparison.csv) から取得できる。",
        "",
        "|図・パネル|E₀ (eV)|深さ (Å)|成分・量|ħγ₀ (eV)|読取り ± 誤差|k_parallel|k_total|点の判定|",
        "|---|---:|---:|---|---:|---|---|---|---|",
    ])
    for rows in groups:
        for r in rows:
            name = r["component"] + " " + r["quantity"].replace("_", " ")
            if r["loss_eV"]:
                name += f' at {r["loss_eV"]} eV'
            values = [f'{float(r[ch]):.7g} {"○" if r[ch + "_within"] else "×"}'
                      for ch in model.MOMENTA]
            lines.append(f'|{r["figure"]} {r["panel"]}|{r["E0_eV"]}|{r["depth_A"] or "—"}|{name}|'
                         f'{r["width_eV"]:.7g}{" (事後)" if not r["policy_width"] else ""}|'
                         f'{r["reading"]} ± {r["uncertainty"]}|'
                         f'{values[0]}|{values[1]}|{labels[r["status"]]}|')
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def generate(output):
    with (HERE / "readings.csv").open(encoding="utf-8", newline="") as src:
        groups = [compare_point(point) for point in csv.DictReader(src)]
    output.mkdir(parents=True, exist_ok=True)
    with (output / "comparison.csv").open("w", encoding="utf-8", newline="") as dst:
        writer = csv.DictWriter(dst, fieldnames=list(groups[0][0]), lineterminator="\n")
        writer.writeheader()
        for rows in groups:
            writer.writerows(rows)
    write_results(groups, output / "RESULTS.md")


def close_numbers(expected, actual, rel_tol):
    try:
        a, b = float(expected), float(actual)
    except ValueError:
        return False
    return math.isfinite(a) and math.isfinite(b) and math.isclose(a, b, rel_tol=rel_tol, abs_tol=ABS_TOL)


def csv_differences(reference, generated):
    with reference.open(encoding="utf-8", newline="") as ref, generated.open(encoding="utf-8", newline="") as new:
        expected_rows, actual_rows = csv.reader(ref), csv.reader(new)
        header, actual_header = next(expected_rows, None), next(actual_rows, None)
        if not header or header != actual_header:
            yield f"header: expected {header!a}, got {actual_header!a}"
            return
        for line, (expected, actual) in enumerate(zip_longest(expected_rows, actual_rows), start=2):
            if expected is None or actual is None or len(expected) != len(header) or len(actual) != len(header):
                yield f"row {line}: expected {expected!a}, got {actual!a}"
                continue
            for column, a, b in zip(header, expected, actual):
                if column in NUMERIC_COLUMNS and (a or b):
                    tolerance = MARGIN_REL_TOL if column.endswith("_margin") else CSV_REL_TOL
                    same = close_numbers(a, b, tolerance)
                else:
                    same = a == b
                if not same:
                    yield f"row {line}, {column}: expected {a!a}, got {b!a}"


def results_differences(reference, generated):
    with reference.open(encoding="utf-8") as ref, generated.open(encoding="utf-8") as new:
        for line, (expected, actual) in enumerate(zip_longest(ref, new), start=1):
            if expected == actual:
                continue
            # Even fields contain all prose, delimiters, counts and verdict marks.
            left, right = RESULT_VALUE.split(expected or ""), RESULT_VALUE.split(actual or "")
            if left[::2] != right[::2]:
                yield f"line {line}: expected {expected!a}, got {actual!a}"
                continue
            for column, a, b in zip(model.MOMENTA, left[1::2], right[1::2]):
                if not close_numbers(a, b, RESULT_REL_TOL):
                    yield f"line {line}, {column}: expected {a!a}, got {b!a}"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--out", type=Path, default=HERE)
    mode.add_argument("--check", action="store_true", help="regenerate and compare both distributed files")
    args = parser.parse_args(argv)
    if args.check:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)
            generate(output)
            ok = True
            for name, compare in (("comparison.csv", csv_differences), ("RESULTS.md", results_differences)):
                differences = list(compare(HERE / name, output / name))
                print(f"{name}: {'differs' if differences else 'matches within tolerance'}")
                for detail in differences:
                    print(detail)
                if differences:
                    ok = False
            return 0 if ok else 1
    generate(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
