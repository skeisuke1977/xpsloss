"""Mapping between the single-loss record and dimensionless order outputs."""

from fractions import Fraction as Q
import math

import numpy as np

from ._fold import _up
from ._order_bounds import corrected, multiplier_rho, sensitivity
from ._order_quadrature import area_terms, edges_of
from ._order_types import OrderAreas, PARTS, SpectrumOrder
from ._single_loss_types import COLUMNS, COMPONENTS, EPS


def seed_errors(single):
    m = Q(single.multiplier_A)
    return {
        c: _up(m * Q(single.window_area_error[c]) + Q(EPS) * abs(Q(float(m * Q(single.window_area[c])))))
        for c in COLUMNS
    }


def unpack_regions(geometry, rows):
    result = dict(outside={"low": rows[0], "high": rows[-1]})
    for c in COMPONENTS:
        for key in (c, c + "_error"):
            result[key] = np.array([row[key] for row in rows[1:-1]])
        key = c + "_error_parts"
        result[key] = {p: np.array([row[key][p] for row in rows[1:-1]]) for p in (*PARTS, "multiplier_interval")}
    return result


def first_order(single, panels, geometry, features, seeds, zero):
    m, rho = Q(single.multiplier_A), Q(single.diagnostics.get("multiplier", {}).get("rho", 0.0))
    area, errors = {}, {}
    parts = {key: dict.fromkeys(COLUMNS, 0.0) for key in (*PARTS, "multiplier_interval", "convention")}
    for c in COLUMNS:
        area[c] = float(m * Q(single.window_area[c]))
        base = seeds[c]
        errors[c], parts["multiplier_interval"][c] = corrected(area[c], base, rho)
        parts["inherited"][c] = base
    area["total"] = single.single_loss_weight
    errors["total"] = single.single_loss_weight_error
    if single.multiplier_convention == "inverse_window_area":
        parts["convention"]["total"] = errors["total"]
        parts["inherited"]["total"] = parts["multiplier_interval"]["total"] = 0.0
    elif single.multiplier_interval_A is None:
        parts["inherited"]["total"] = errors["total"]
    else:
        parts["inherited"]["total"] = single.multiplier_A * single.window_area_error["total"] + EPS * abs(
            area["total"]
        )
        parts["multiplier_interval"]["total"] = _up(Q(errors["total"]) - Q(parts["inherited"]["total"]))
    edges = edges_of(panels)
    diag = single.diagnostics
    knots = sorted(set((*single.loss_window_eV, *features.get("jumps_eV", ()), *features.get("kinks_eV", ()))))
    mesh = dict(
        knots_eV=[x for x in knots if edges[0] <= x <= edges[-1]],
        initial_panel_edges_eV=diag.get("gradings", edges),
        panel_edges_eV=edges,
        panels=len(panels),
        builds=1,
        stages=diag["stages"],
        history=[],
        local_area_terms=[area_terms(p).tolist() for p in panels],
        formal_folds=0 if not geometry else 1,
        quadrature_points=0,
        cells=diag.get("cells", 0),
    )
    kw = dict(
        order=1,
        area=area,
        area_error=errors,
        area_error_parts=parts,
        area_zero_kind=zero_kinds(area, zero),
        cutoff_sensitivity=single.cutoff_sensitivity.copy(),
        cutoff_sensitivity_error=single.cutoff_sensitivity_error.copy(),
        mesh=mesh,
    )
    if geometry is None:
        return OrderAreas(**kw)
    rows = []
    for i, name in enumerate(geometry.names):
        side = name.split(".")[-1]
        row = {}
        for c in COMPONENTS:
            key = "single_loss_" + c

            def read(field):
                return single.outside[side][field] if name.startswith("outside") else getattr(single, field)[i - 1]

            row[c], row[c + "_error"] = read(key), read(key + "_error")
            base = diag["error_budget"]["outputs"][name + "." + key]["estimate"]
            row[c + "_error_parts"] = dict(
                inherited=base,
                propagated=0.0,
                local=0.0,
                multiplier_interval=_up(Q(row[c + "_error"]) - Q(base)) if rho else 0.0,
            )
        rows.append(row)
    return SpectrumOrder(**kw, **unpack_regions(geometry, rows))


def zero_kinds(areas, structural):
    return {c: ("exact_zero" if structural else "rounded_zero") if x == 0 else None for c, x in areas.items()}


def balances(orders, single, geometry, ratio):
    """Residuals of the area product, charge decomposition and detector mass identities."""
    report = {}
    for o in orders:
        if o.order == 1:
            continue
        power = Q(ratio["value"]) ** (o.order - 1)
        row = dict(area_product={c: float(Q(o.area[c]) - Q(orders[0].area[c]) * power) for c in COLUMNS})
        if geometry:
            row["region_sum"] = {}
            for c in COMPONENTS:
                terms = (o.outside["low"][c], *getattr(o, c), o.outside["high"][c])
                row["region_sum"][c] = float(sum(map(Q, terms)) - Q(single.primary_integral) * Q(o.area[c]))
            row["components"] = {}
            for i, name in enumerate(geometry.names):
                if name.startswith("outside"):
                    v = o.outside[name.split(".")[-1]]
                    total, electron, hole = (v[c] for c in COMPONENTS)
                else:
                    total, electron, hole = (getattr(o, c)[i - 1] for c in COMPONENTS)
                row["components"][name] = float(Q(total) - Q(electron) - Q(single.hole_charge) * Q(hole))
        report[o.order] = row
    return report


def higher_order(n, single, first, panels1, geometry, values, parts, folded, mesh, seeds, zero):
    rho = multiplier_rho(single, n)
    areas = dict(zip(COLUMNS, values))
    area_parts = {p: dict(zip(COLUMNS, v)) for p, v in zip(PARTS, parts)}
    area_parts["convention"] = dict.fromkeys(COLUMNS, 0.0)
    area_parts["multiplier_interval"] = {}
    errors = {}
    for c, y in areas.items():
        fixed = _up(sum(Q(area_parts[p][c]) for p in PARTS))
        errors[c], area_parts["multiplier_interval"][c] = corrected(y, fixed, rho)
    h, he = sensitivity(first, panels1, seeds, single, n, rho)
    kw = dict(
        order=n,
        area=areas,
        area_error=errors,
        area_error_parts=area_parts,
        area_zero_kind=zero_kinds(areas, zero),
        cutoff_sensitivity=h,
        cutoff_sensitivity_error=he,
        mesh=mesh,
    )
    if geometry is None:
        return OrderAreas(**kw)
    rows = []
    for name in geometry.names:
        row = {}
        for c in COMPONENTS:
            data = folded[name + "." + c]
            y = data["value"]
            fixed = _up(sum(Q(data[p]) for p in PARTS))
            err, extra = corrected(y, fixed, rho)
            row[c], row[c + "_error"] = y, err
            row[c + "_error_parts"] = {**{p: data[p] for p in PARTS}, "multiplier_interval": extra}
        rows.append(row)
    return SpectrumOrder(**kw, **unpack_regions(geometry, rows))


def spectrum_sum(single, orders, structural_zero):
    def add(y, error, terms):
        terms = list(terms)
        if not terms:
            return y, error
        exact = Q(y) + sum(Q(term) for term, _ in terms)
        value = float(exact)
        uncertainty = Q(error) + sum(Q(e) for _, e in terms) + abs(exact - Q(value))
        return value, _up(uncertainty)

    if len(orders) == 1:
        return dict(
            spectrum=single.spectrum.copy(),
            spectrum_error=single.spectrum_error.copy(),
            outside={s: {k: v[k] for k in ("spectrum", "spectrum_error")} for s, v in single.outside.items()},
        )
    active = () if structural_zero else orders
    y, e = np.array(
        [
            add(v, err, ((o.total[i], o.total_error[i]) for o in active))
            for i, (v, err) in enumerate(zip(single.zero_loss, single.zero_loss_error))
        ]
    ).T
    outside = {}
    for s, row in single.outside.items():
        v, err = add(
            row["zero_loss"],
            row["zero_loss_error"],
            ((o.outside[s]["total"], o.outside[s]["total_error"]) for o in active),
        )
        outside[s] = dict(spectrum=v, spectrum_error=err)
    if (
        not np.isfinite(y).all()
        or not np.isfinite(e).all()
        or not all(math.isfinite(v) for row in outside.values() for v in row.values())
    ):
        raise OverflowError("nonfinite sum of loss orders")
    return dict(spectrum=y, spectrum_error=e, outside=outside)
