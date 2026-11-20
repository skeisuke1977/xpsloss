"""Regenerate the figure comparisons from the recorded readings and independent equations."""

import argparse
from collections import Counter
import csv
from functools import lru_cache
from pathlib import Path

import numpy as np

if __package__:
    from . import model
else:
    import model

HERE = Path(__file__).resolve().parent
WITHIN = "within reading uncertainty"
OUTSIDE = "outside"
DEPENDENT = "depends on width assumption"


def widths(figure, energy):
    """Table I is on p. 1615. It has 1180 eV, not 1000 eV."""
    if figure == 3:
        return [(1.5, "assumed from 1180 eV"), (1.6, "assumed sensitivity"),
                (1.7, "assumed sensitivity")]
    if figure == 8:
        return [(1.44, "assumed interpolation (1180/2000 eV)"), (1.5, "assumed from 1180 eV")]
    return {
        200: [(3.4, "assumed from 175 eV"), (2.1, "assumed from 500 eV")],
        500: [(2.1, "Table I")],
        1000: [(1.5, "assumed from 1180 eV")],
        2000: [(1.3, "Table I")],
        5000: [(1.3, "assumed from 2000 eV")],
        10000: [(1.3, "assumed from 2000 eV")],
    }[energy]


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


def computed(point, width, momentum):
    figure = int(point["figure"])
    energy = float(point["E0_eV"])
    depth = float(point["depth_A"]) if point["depth_A"] else None
    component = point["component"]
    quantity = point["quantity"]
    if component == "bulk":
        momentum = None  # Eq. (1) has a single integration variable k.

    def at(loss):
        return components(figure, energy, depth, width, momentum, loss)[component]

    if quantity == "value":
        return at(float(point["loss_eV"]))
    if figure == 3:
        low, high, step = (8 if component == "bulk" else 13), 25, 0.05
    elif quantity.startswith("surface"):
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
    return grid[peak] if quantity.endswith("position") else values[peak]


def compare_point(point):
    observed, error = float(point["reading"]), float(point["uncertainty"])
    rows = []
    for width, source in widths(int(point["figure"]), float(point["E0_eV"])):
        row = dict(point, width_eV=width, width_source=source)
        for choice in model.MOMENTA:
            value = computed(point, width, choice)
            row[choice] = f"{value:.12g}"
            # Decimal grid positions exactly one reading interval away are inside.
            row[choice + "_within"] = bool(abs(value - observed) <= error + 1e-14)
        rows.append(row)
    outcomes = [{row[ch + "_within"] for row in rows} for ch in model.MOMENTA]
    if any(len(flags) > 1 for flags in outcomes):
        status = DEPENDENT
    elif any(False in flags for flags in outcomes):
        status = OUTSIDE
    elif len(rows) == 1 and rows[0]["width_source"] != "Table I":
        status = DEPENDENT
    else:
        status = WITHIN
    for row in rows:
        row["status"] = status
    return rows


def write_results(groups, destination):
    counts = Counter(rows[0]["status"] for rows in groups)
    lines = [
        "# 原著の図3・5・8との比較",
        "",
        "原図は A. Cohen Simonsen, F. Yubero & S. Tougaard, Phys. Rev. B 56, 1612 (1997) の"
        "図3（p. 1616）、図5(a,b)（p. 1617）、図8（p. 1618）。[原著](https://doi.org/10.1103/PhysRevB.56.1612)",
        "",
        "独立実装で式(20)を使い、法線出射の比較点ごとに紙面の読取りと照合した。"
        "同じ点の線幅違いは一度だけ数える。",
        "",
        "| 比較点 | 読取り誤差内 | 誤差外 | 幅の仮定に依存 |",
        "|---:|---:|---:|---:|",
        f"| {len(groups)} | {counts[WITHIN]} | {counts[OUTSIDE]} | {counts[DEPENDENT]} |",
        "",
        "各行に一つの幅を示す。○はその幅で読取り誤差内、×は誤差外。"
        "全体判定の規則、読取り方、単位、再計算コマンドは [README](README.md) に記載した。"
        "数値と各運動量の判定は [比較表](comparison.csv) から取得できる。",
        "",
        "|図・パネル|E₀ (eV)|深さ (Å)|成分・量|ħγ₀ (eV)|読取り ± 誤差|k_parallel|k_total|点の判定|",
        "|---|---:|---:|---|---:|---|---|---|---|",
    ]
    labels = {WITHIN: "誤差内", OUTSIDE: "誤差外", DEPENDENT: "幅依存"}
    for rows in groups:
        for r in rows:
            name = r["component"] + " " + r["quantity"].replace("_", " ")
            if r["loss_eV"]:
                name += f' at {r["loss_eV"]} eV'
            values = [f'{float(r[ch]):.7g} {"○" if r[ch + "_within"] else "×"}'
                      for ch in model.MOMENTA]
            lines.append(f'|{r["figure"]} {r["panel"]}|{r["E0_eV"]}|{r["depth_A"] or "—"}|{name}|'
                         f'{r["width_eV"]}|{r["reading"]} ± {r["uncertainty"]}|'
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=HERE)
    generate(parser.parse_args().out)
