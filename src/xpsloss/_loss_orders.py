"""Finite-window first-event/13a chains on continuous loss-energy panels."""

from fractions import Fraction as Q
from itertools import chain, pairwise
import math

import numpy as np

from . import _fold, _loss_panels as lp
from ._order_bounds import ratio_record, truncation
from ._order_quadrature import Polynomials, area_terms, build_panel, first_envelopes, mandatory_knots, minimum
from ._order_results import balances, first_order, higher_order, seed_errors, spectrum_sum
from ._order_types import (
    MULTIPLE_LOSS_MODEL_ID,
    PARTS,
    MultipleLossConvergenceError,
    MultipleLossOrderAreasResult,
    MultipleLossSpectrumResult,
    preflight,
)
from ._single_loss_types import COLUMNS, COMPONENTS, EPS


def _order_cell(panel, parts, left, right, g0, g1, multiplier_power):
    weights, products, norm, mass_error, coord_error = _fold._cell_weights(panel, left, right, g0, g1)
    width = right - left
    scale = width * norm
    values = panel.values[_fold.COMP_INDEX]
    folded = _scaled(values @ weights, multiplier_power * scale)
    # Use the same numerical weights as the single-loss fold. Keep all error
    # products in panel units, even when an individual product is subnormal.
    abs_weights = list(map(Q, abs(weights)))
    abs_products = list(map(Q, products))
    envelope, slope = panel.component_bounds
    errors = []
    for j, i in enumerate(_fold.COMP_INDEX):
        inherited, propagated, local = (
            sum(Q(v) * w for v, w in zip(parts[part][i], abs_weights)) for part in PARTS
        )
        uncertain = (abs(Q(v)) + Q(e) for v, e in zip(values[j], parts["local"][i]))
        rounding = 48 * Q(EPS) * sum(v * w for v, w in zip(uncertain, abs_products))
        interpolation = min(width * Q(panel.coeff_diff[i]), Q(panel.interpolation[i]))
        local_error = (
            scale * (local + Q(panel.verification[i]) + rounding + coord_error * Q(slope[j]))
            + norm * interpolation + width * mass_error * Q(envelope[j])
        )
        errors.append((scale * inherited, scale * propagated, local_error))
    return folded, errors


def _fold_order(panels, geometry, *, multiplier_power, error_parts):
    """Direct satellite integrals; nodal error sources share the same cell weights."""
    rows = {}
    for region, name in enumerate(geometry.names):
        values = np.zeros((len(panels), len(COMPONENTS)))
        errors = np.zeros((len(panels), len(COMPONENTS), len(PARTS)), dtype=object)
        for i, (p, parts) in enumerate(zip(panels, error_parts)):
            a, b = Q(p.left_eV), Q(p.right_eV)
            for left, right in pairwise(chain((a,), geometry.breaks(region, a, b), (b,))):
                g0 = geometry.mass(region, left, "right")
                g1 = geometry.mass(region, right, "left")
                if not (g0 or g1):
                    continue
                cell_values, cell_errors = _order_cell(p, parts, left, right, g0, g1, multiplier_power)
                values[i] += cell_values
                errors[i] += cell_errors
        for j, c in enumerate(COMPONENTS):
            totals = _scaled(errors[:, j, :].sum(axis=0), multiplier_power, up=True)
            rows[name + "." + c] = dict(
                value=math.fsum(values[:, j]), **dict(zip(PARTS, totals)),
                local_by_panel=_scaled(errors[:, j, 2], multiplier_power, up=True),
            )
    return rows


def _scaled(values, factor, *, up=False):
    rounding = _fold._up if up else float
    return np.array([rounding(factor * Q(v)) for v in values])


def _judge_inherited(rows, n, phi, s, mass, m, diagnostics):
    exceeded = []
    for key, c, value, error, is_area in rows:
        absolute = phi[COLUMNS.index(c)] * (
            m * s["area_epsabs"] if is_area else mass * s["spectrum_epsabs_fraction"]
        )
        relative = n * s["epsrel"] * abs(value)
        if error > max(absolute, relative):
            side = "relative" if relative >= absolute else "absolute"
            suggestion = (
                "epsrel" if side == "relative" else "area_epsabs" if is_area else "spectrum_epsabs_fraction"
            )
            exceeded.append(
                dict(
                    key=key,
                    inherited=error,
                    budget_relative=relative,
                    budget_absolute=absolute,
                    side=side,
                    suggestion=suggestion,
                )
            )
    if exceeded:

        def severity(row):
            budget = max(row["budget_absolute"], row["budget_relative"])
            return row["inherited"] / budget if budget else math.inf

        worst = max(exceeded, key=severity)
        raise MultipleLossConvergenceError(
            "inherited_error_exceeds_order_budget",
            order=n,
            stage="order_fold"
            if "." in worst["key"] and not worst["key"].startswith("area.")
            else "order_refinement",
            exceeded=exceeded,
            key=worst["key"],
            column=worst["key"].split(".")[-1],
            side=worst["side"],
            suggestion=worst["suggestion"],
            diagnostics=diagnostics,
        )


def _choose_split(contributions, total_error, value, absolute, relative, split):
    budget = max(absolute, relative * abs(value))
    if total_error > budget:
        split |= contributions > 0.5 * budget / len(split)


def _preliminary_bins(panels, geometry, mass, factor, settings, split):
    """Use a binary64 primary CDF only to nominate panels for bisection."""
    if not mass:
        return
    points = np.array([float(x) for x, _ in geometry.lines])
    perm = np.argsort(points)
    positions = points[perm]
    cumulative = np.r_[0.0, np.cumsum([float(geometry.lines[i][1]) for i in perm])]
    bins = [(float(a), float(b), float((b - a) * density)) for a, b, density in geometry.bins]
    edges = np.asarray(geometry.edges)
    try:
        scale = float(factor)
    except OverflowError:
        return  # The exact-rational formal fold still applies.
    # Bound working memory independently of the number of detector edges.
    for start in range(0, len(geometry.names), 128):
        stop = min(start + 128, len(geometry.names))
        boundaries = np.r_[-np.inf, edges, np.inf][start : stop + 1]
        contributions, values = [], []
        for p in panels:
            t = lp._nodes(p.left_eV, p.right_eV)[0]
            at = boundaries[:, None] + t
            cdf = cumulative[np.searchsorted(positions, at, side="left")]
            for a, b, intensity in bins:
                cdf += intensity * np.clip((at - a) / (b - a), 0.0, 1.0)
            weights = np.diff(cdf, axis=0) * (lp.CC * (p.right_eV - p.left_eV) / 2)
            v = p.values[_fold.COMP_INDEX]
            density_error = p.errors[_fold.COMP_INDEX] + (p.coeff_diff + p.verification)[_fold.COMP_INDEX, None]
            contributions.append((weights @ density_error.T) * scale)
            values.append((weights @ v.T) * scale)
        terms = np.array(contributions)
        budgets = np.maximum(
            settings["order_epsabs_fraction"] * mass, settings["order_epsrel"] * abs(np.sum(values, axis=0))
        )
        failing = terms.sum(axis=0) > budgets
        if failing.any():
            split |= (terms[:, failing] > 0.5 * budgets[failing] / len(panels)).any(axis=1)


def _adapt(
    left,
    right,
    window,
    features,
    n,
    s,
    geometry,
    factor,
    diagnostics,
    order_inject,
    order_limit,
    min_width,
    inherited_judge=None,
    mass=0.0,
    transform=None,
):
    knots = mandatory_knots(features, window, n, order_limit)
    bounds = list(pairwise(knots))
    cache = {}
    resources = diagnostics["resources"]
    start = resources["quadrature_points_total"]
    history = []
    snapshot = dict(order=n, history=history, panels=len(bounds), stages=0)
    diagnostics["current_mesh"] = snapshot
    formal = 0

    def charge(points):
        resources["quadrature_points_total"] += points
        measured = resources["quadrature_points_total"]
        if measured > s["max_convolution_points"]:
            raise MultipleLossConvergenceError(
                "convolution_work_limit_reached",
                stage="order_refinement",
                order=n,
                resource="max_convolution_points",
                measured=measured,
                limit=s["max_convolution_points"],
                diagnostics=diagnostics,
            )
        return s["max_convolution_points"] - measured

    while True:
        edges = [bounds[0][0], *(b for a, b in bounds)]
        snapshot.update(panels=len(bounds), stages=len(history), panel_edges_eV=edges)
        cells = geometry.cell_count(edges, s["max_cells"]) if geometry else 0
        resources["cells_max"] = max(resources["cells_max"], cells)
        if geometry and cells > s["max_cells"]:
            raise MultipleLossConvergenceError(
                "convolution_work_limit_reached",
                order=n,
                stage="order_fold",
                resource="max_cells",
                measured=cells,
                limit=s["max_cells"],
                diagnostics=diagnostics,
            )
        items = [
            cache[pair] if pair in cache else build_panel(left, right, pair, n, charge, order_inject, transform)
            for pair in bounds
        ]
        panels, envelopes, parts = zip(*items)
        local = np.array([area_terms(p) for p in panels])
        areas = _scaled(
            np.sum([p.values @ (lp.CC * ((p.right_eV - p.left_eV) / 2)) for p in panels], axis=0), factor
        )
        split = np.zeros(len(panels), dtype=bool)
        local_error = _scaled([sum(map(Q, col)) for col in local.T], factor, up=True)
        physical_local = np.array([_scaled(terms, factor, up=True) for terms in local])
        for c, v in enumerate(areas):
            _choose_split(physical_local[:, c], local_error[c], v, s["order_epsabs_fraction"], s["order_epsrel"], split)
        if geometry and not split.any():
            _preliminary_bins(panels, geometry, mass, factor, s, split)
            # A cheap prediction cannot reject a run at a resource boundary.
            if len(panels) + split.sum() > order_limit or any(
                refine
                and (
                    not lp._nodes(a, lp._midpoint(a, b))[2]
                    or not lp._nodes(lp._midpoint(a, b), b)[2]
                    or (min_width is not None and (b - a) / 2 < min_width)
                )
                for refine, (a, b) in zip(split, bounds)
            ):
                split.fill(False)
        folded = None
        if not split.any():
            if geometry:
                formal += 1
                try:
                    folded = _fold_order(list(panels), geometry, multiplier_power=factor, error_parts=list(parts))
                except ArithmeticError as exc:
                    raise MultipleLossConvergenceError(
                        "nonfinite_output",
                        order=n,
                        stage="order_fold",
                        inner=dict(type=type(exc).__name__, reason="nonfinite_output", stage="order_fold"),
                        diagnostics=diagnostics,
                    ) from exc
                for row in folded.values():
                    _choose_split(
                        row["local_by_panel"],
                        row["local"],
                        row["value"],
                        mass * s["order_epsabs_fraction"],
                        s["order_epsrel"],
                        split,
                    )
            if inherited_judge is not None:
                inherited_judge(areas, folded)
        history.append(
            dict(
                panels=len(panels),
                formal=folded is not None,
                passed=not split.any(),
                split_panels=int(split.sum()),
            )
        )
        if not split.any():
            mesh = dict(
                knots_eV=knots,
                initial_panel_edges_eV=knots,
                panel_edges_eV=edges,
                panels=len(panels),
                builds=1,
                stages=len(history) - 1,
                history=history,
                local_area_terms=local.tolist(),
                formal_folds=formal,
                cells=cells,
                quadrature_points=resources["quadrature_points_total"] - start,
            )
            return (
                list(panels),
                np.array(envelopes),
                areas,
                local_error,
                folded,
                mesh,
            )
        if len(panels) + split.sum() > order_limit:
            raise MultipleLossConvergenceError(
                "order_budget_not_met",
                order=n,
                stage="order_refinement",
                panel_limit=order_limit,
                panel_eV=bounds[int(np.argmax(split))],
                diagnostics=diagnostics,
            )
        cache = {pair: item for pair, item, refine in zip(bounds, items, split) if not refine}
        next_bounds = []
        for (a, b), refine in zip(bounds, split):
            if refine:
                mid = lp._midpoint(a, b)
                if not a < mid < b or (min_width is not None and min(mid - a, b - mid) < min_width):
                    raise MultipleLossConvergenceError(
                        "order_mesh_unresolved",
                        order=n,
                        stage="order_refinement",
                        panel_eV=(a, b),
                        diagnostics=diagnostics,
                    )
                next_bounds.extend(((a, mid), (mid, b)))
            else:
                next_bounds.append((a, b))
        bounds = next_bounds


def _empty_mesh():
    return dict(
        knots_eV=[],
        initial_panel_edges_eV=[],
        panel_edges_eV=[],
        panels=0,
        builds=0,
        stages=0,
        history=[],
        local_area_terms=[],
        formal_folds=0,
        cells=0,
        quadrature_points=0,
    )


def _orders(
    single, panels1, nmax, features, settings, geometry, diagnostics, zero, order_inject, panel_limit, min_width
):
    seeds = seed_errors(single)
    diagnostics["first_order"]["fixed_multiplier_area_error"] = seeds
    first = first_order(single, panels1, geometry, features, seeds, zero)
    orders = [first] if nmax else []
    panel_map = {1: panels1}
    m = Q(single.multiplier_A)
    mass = single.primary_abs_integral if geometry else 0.0
    env = first_envelopes(panels1)
    left = Polynomials(panels1, env, slice(None))
    right = Polynomials(panels1, env, slice(0, 1))
    rnorm = m * right.norms()[0]
    first_norm = [m * v for v in left.norms()]
    prev_errors = [list(map(Q, seeds.values())), [0] * 11, [0] * 11]
    er = Q(seeds["term_13a"])
    for n in range(2, nmax + 1):
        try:
            factor = m**n
            if zero:
                values = np.zeros(11)
                parts = np.zeros((3, 11))
                folded = (
                    {
                        key + "." + c: dict(value=0.0, inherited=0.0, propagated=0.0, local=0.0)
                        for key in geometry.names
                        for c in COMPONENTS
                    }
                    if geometry
                    else None
                )
                mesh = _empty_mesh()
                panel_map[n] = []
            else:
                prev_norm = [m ** (n - 1) * v for v in left.norms()]
                inherited, propagated = [], []
                for norm, (ei, ep, el) in zip(prev_norm, zip(*prev_errors)):
                    inherited.append(_fold._up(rnorm * ei + norm * er + (ei + ep + el) * er))
                    propagated.append(_fold._up(rnorm * (ep + el)))
                phi = [float(rnorm ** (n - 1) + (n - 1) * v * rnorm ** (n - 2)) for v in first_norm]

                def judge(areas, folded):
                    rows = [("area." + c, c, v, e, True) for c, v, e in zip(COLUMNS, areas, inherited)]
                    if folded:
                        rows.extend(
                            (key, key.split(".")[-1], r["value"], r["inherited"], False)
                            for key, r in folded.items()
                        )
                    _judge_inherited(rows, n, phi, settings, mass, float(m), diagnostics)

                panels, env, values, local, folded, mesh = _adapt(
                    left,
                    right,
                    single.loss_window_eV,
                    features,
                    n,
                    settings,
                    geometry,
                    factor,
                    diagnostics,
                    order_inject,
                    panel_limit,
                    min_width,
                    judge,
                    mass,
                )
                panel_map[n] = panels
                left = Polynomials(panels, env, slice(None))
                parts = [inherited, propagated, local]
                prev_errors = [list(map(Q, part)) for part in parts]
            orders.append(
                higher_order(n, single, first, panels1, geometry, values, parts, folded, mesh, seeds, zero)
            )
        except MultipleLossConvergenceError as exc:
            if exc.diagnostics is None:
                exc.diagnostics = diagnostics
            raise
        except (ArithmeticError, FloatingPointError) as exc:
            raise MultipleLossConvergenceError(
                "nonfinite_output", order=n, stage="order_refinement", diagnostics=diagnostics
            ) from exc
    return tuple(orders), panel_map


_WARNINGS = {
    "model_is_approximation": "The finite-window loss-order expansion is approximate.",
    "excess_enters_once": "Surface, hole and interference contributions enter only in the first event.",
    "depth_correlation_dropped": "The expansion omits trajectory-dependent depth correlations.",
    "no_elastic_scattering": "The expansion omits elastic scattering.",
    "order_areas_not_probabilities": "Order areas are dimensionless weights.",
    "truncated": "Higher orders may reach the displayed kinetic energies.",
    "order_area_ratio_not_below_one": "The ratio interval reaches or exceeds one.",
    "multiplier_differs_from_weight_length": "The multiplier differs from the depth weight length.",
    "cumulative_loss_over_E0_gt_0p1": "The largest cumulative loss exceeds 0.1 times E0.",
    "negative_total_kernel": "The total kernel takes negative values.",
    "negative_13a_kernel": "The repeated term 13a takes negative values.",
    "dispersive_k_parallel_13a": "For dispersive k_parallel, term 13a differs from the infinite-medium bulk kernel.",
    "E0_above_5keV": "The kernel energy exceeds 5 keV.",
}


def finish(
    single,
    panels,
    plan,
    requested,
    settings,
    features,
    *,
    deliveries=1,
    zero=False,
    dispersive=False,
    order_inject=None,
    order_panel_limit=None,
    order_min_panel_width_eV=None,
    trace=None,
):
    n, support, geometry, lower = plan
    estimate = 0 if zero else 280.5 * sum(k * (len(panels) + k - 1) for k in range(2, n + 1))
    resources = dict(
        work_lower_bound_points=lower,
        work_estimate_after_first_order=estimate,
        quadrature_points_total=0,
        cells_max=single.diagnostics.get("cells", 0),
    )
    diag = dict(
        first_order=dict(
            beta_deliveries=deliveries,
            final_panels=len(panels),
            area_constraint="inverse_window_area"
            if single.multiplier_convention == "inverse_window_area"
            else None,
            kernel_calls=single.diagnostics["kernel_calls"],
        ),
        resources=resources,
        effective_settings=settings,
    )
    if estimate > settings["max_convolution_points"]:
        raise MultipleLossConvergenceError(
            "convolution_work_limit_exceeded",
            stage="after_first_order",
            resource="max_convolution_points",
            estimate=estimate,
            limit=settings["max_convolution_points"],
            diagnostics=diag,
        )
    try:
        orders, panel_map = _orders(
            single,
            panels,
            n,
            features,
            settings,
            geometry,
            diag,
            zero,
            order_inject,
            order_panel_limit or settings["panel_limit"],
            order_min_panel_width_eV,
        )
        ratio = ratio_record(single)
        diag["balance"] = balances(orders, single, geometry, ratio)
    except ArithmeticError as exc:
        raise MultipleLossConvergenceError(
            "nonfinite_output", order=n, stage="after_first_order", diagnostics=diag
        ) from exc
    diag.pop("current_mesh", None)
    if orders and trace is not None:
        orders[0].mesh.update(
            initial_panel_edges_eV=trace["initial_panel_edges_eV"],
            history=trace["history"],
            formal_folds=sum(row["formal"] for row in trace["history"]),
        )
    convention = "not_determined_by_convention"
    if single.weight_length_kind == "inverse_13a_window_area" and single.multiplier_convention == "weight_length":
        convention = "below_one" if single.extra_inverse_length_per_A else "exactly_one"
    minima = {}
    for c in ("total", "term_13a"):
        i = COLUMNS.index(c)
        node_min = float(trace["sampled_min"][i]) if trace else min(float(p.values[i].min()) for p in panels)
        cont_min = min(minimum(p, i) for p in panels)
        minima[c] = dict(sampled_min=node_min, continuous_min=cont_min)
    diag["negative_values"] = minima
    negative = {k: min(v.values()) < 0 for k, v in minima.items()}
    energy = single.electron_energy_eV
    cumulative = float(n * Q(single.loss_window_eV[1]) / Q(energy)) if energy else None
    diag["cumulative_loss_over_E0"] = cumulative
    kw = dict(
        single_loss=single,
        orders=orders,
        order_area_ratio=ratio,
        convention_ratio=convention,
        model_id=MULTIPLE_LOSS_MODEL_ID,
        depends_on_model_id=single.model_id,
        diagnostics=diag,
    )
    state = None
    if geometry:
        try:
            kw["truncation"] = truncation(
                single,
                panels,
                geometry,
                requested,
                n,
                support,
                ratio,
                settings["neglect_tolerance_fraction"],
                negative["term_13a"],
                diag,
            )
            state = kw["truncation"]["state"]
            kw.update(
                spectrum_sum(single, orders, zero),
                zero_loss=single.zero_loss.copy(),
                zero_loss_error=single.zero_loss_error.copy(),
            )
        except ArithmeticError as exc:
            raise MultipleLossConvergenceError("nonfinite_output", stage="truncation", diagnostics=diag) from exc
    codes = list(_WARNINGS)[:5]
    codes.extend(
        k
        for k, active in (
            ("truncated", state == "truncated"),
            ("order_area_ratio_not_below_one", ratio["interval"][1] >= 1),
            ("multiplier_differs_from_weight_length", single.multiplier_A != single.weight_length_A),
            ("cumulative_loss_over_E0_gt_0p1", cumulative is not None and cumulative > 0.1),
            ("negative_total_kernel", negative["total"]),
            ("negative_13a_kernel", negative["term_13a"]),
            ("dispersive_k_parallel_13a", dispersive),
            ("E0_above_5keV", energy is not None and energy > 5000),
        )
        if active
    )
    diag["warnings"] = [
        *single.diagnostics["warnings"],
        *(dict(source="xps_multiple_loss", code=c, message=_WARNINGS[c]) for c in codes),
    ]
    kw["warnings"] = tuple(w["message"] for w in diag["warnings"])
    return (MultipleLossSpectrumResult(**kw) if geometry else MultipleLossOrderAreasResult(**kw)), panel_map


def integrate_orders(
    kernel,
    loss_window_eV,
    *,
    features,
    weight_length_A,
    hole_charge=1.0,
    primary=None,
    output_edges_eV=None,
    model,
    max_order,
    neglect_tolerance_fraction=None,
    order_limit=64,
    max_convolution_points=200_000_000,
    epsrel=1e-6,
    area_epsabs=1e-10,
    spectrum_epsabs_fraction=1e-9,
    panel_limit=400,
    max_cells=2_000_000,
    order_epsrel=1e-6,
    order_epsabs_fraction=1e-9,
    multiplier="weight_length",
    multiplier_A=None,
    initial_edges=None,
    inject=None,
    order_inject=None,
    order_panel_limit=None,
    order_min_panel_width_eV=None,
):
    window = lp.checked_window(loss_window_eV)
    first = dict(
        epsrel=epsrel,
        area_epsabs=area_epsabs,
        spectrum_epsabs_fraction=spectrum_epsabs_fraction,
        panel_limit=panel_limit,
        max_cells=max_cells,
    )
    settings = dict(
        **first,
        order_limit=order_limit,
        max_convolution_points=max_convolution_points,
        order_epsrel=order_epsrel,
        order_epsabs_fraction=order_epsabs_fraction,
        neglect_tolerance_fraction=neglect_tolerance_fraction,
    )
    plan = preflight(model, max_order, window, primary, output_edges_eV, settings)
    with lp.observe_mesh() as trace:
        single, panels = lp.integrate_window(
            kernel,
            window,
            features=features,
            weight_length_A=weight_length_A,
            hole_charge=hole_charge,
            primary=primary,
            output_edges_eV=output_edges_eV,
            initial_edges=initial_edges,
            inject=inject,
            multiplier=multiplier,
            multiplier_A=multiplier_A,
            **first,
        )
    return finish(
        single,
        panels,
        plan,
        max_order,
        settings,
        features,
        order_inject=order_inject,
        order_panel_limit=order_panel_limit,
        order_min_panel_width_eV=order_min_panel_width_eV,
        trace=trace,
    )


def convolve_panels(
    first_panels,
    loss_window_eV,
    *,
    model,
    max_order,
    features,
    order_epsrel=1e-6,
    order_epsabs_fraction=1e-9,
    order_panel_limit=400,
):
    from ._order_types import ORDER_DEFAULTS

    window = lp.checked_window(loss_window_eV)
    settings = dict(ORDER_DEFAULTS, order_epsrel=order_epsrel, order_epsabs_fraction=order_epsabs_fraction)
    preflight(model, max_order, window, None, None, settings)
    env = first_envelopes(first_panels)
    right = Polynomials(first_panels, env, slice(0, 1))
    panels = first_panels
    result = {1: panels}
    diag = dict(resources=dict(quadrature_points_total=0, cells_max=0))
    for n in range(2, max_order + 1):
        panels, env, *_ = _adapt(
            Polynomials(panels, env, slice(None)),
            right,
            window,
            features,
            n,
            settings,
            None,
            Q(1),
            diag,
            None,
            order_panel_limit,
            None,
        )
        result[n] = panels
    return result
