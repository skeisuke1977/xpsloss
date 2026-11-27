"""Private comparison models from M4b design 4.5, for unit hole charge only."""

from fractions import Fraction as Q
import math
from types import SimpleNamespace

import numpy as np

from . import _loss_panels as lp
from ._loss_orders import _adapt
from ._order_bounds import corrected, log_laplace, multiplier_rho, truncation
from ._order_quadrature import Polynomials, first_envelopes
from ._order_results import first_order, seed_errors
from ._order_types import MULTIPLE_LOSS_MODEL_ID, MultipleLossConvergenceError, preflight
from ._single_loss_types import COLUMNS, COMPONENTS, EPS

GEOMETRIC = "full_kernel_geometric_ref"
FACTORIAL = "full_kernel_factorial_ref"
EXCESS = "excess_poisson_ref"
TOTAL, ELECTRON, HOLE = (COLUMNS.index(c) for c in COMPONENTS)
# Two auxiliary columns carry powers of g for total and no_hole respectively.
AUX = [1, 2]
COMP = [TOTAL, ELECTRON]


def _pack(panels, model):
    left, right = [], []
    env = first_envelopes(panels)
    for p, e in zip(panels, env[:, 0]):
        v, errors = np.zeros_like(p.values), np.zeros_like(p.errors)
        v[COMP], errors[COMP] = p.values[COMP], e[COMP, None]
        repeat, uncertainty = v.copy(), errors.copy()
        if model == EXCESS:
            v[AUX] = p.values[COMP] - p.values[0]
            errors[AUX] = e[COMP, None] + e[0] + EPS * abs(v[AUX])
            repeat[COMP], uncertainty[COMP] = p.values[0], e[0]
            repeat[AUX], uncertainty[AUX] = v[AUX], errors[AUX]
        left.append(lp.Panel.from_lobatto(p.left_eV, p.right_eV, v, errors))
        right.append(lp.Panel.from_lobatto(p.left_eV, p.right_eV, repeat, uncertainty))
    return left, right


def _transform(block, n, model):
    if model == FACTORIAL:
        block /= n
    elif model == EXCESS:
        # C_n = C_(n-1)*R + g^(*n)/n!, applied to each full kernel.
        block[:, AUX] /= n
        block[:, COMP] += block[:, AUX]
    block[0, HOLE] = block[0, TOTAL] - block[0, ELECTRON]
    block[1:, HOLE] = block[1:, TOTAL] + block[1:, ELECTRON]
    block[3, HOLE] += EPS * (abs(block[0, TOTAL]) + abs(block[0, ELECTRON]))


def _log_tail(model, panels, m):
    def evaluate(taus, n, logk, logr):
        if model == GEOMETRIC:
            ans = np.full(len(logk), np.inf)
            ok = logk < 0
            ans[ok] = (n + 1) * logk[ok] - np.log(-np.expm1(logk[ok]))
            return ans
        if model == FACTORIAL:
            # exp(k) bounds the ratio of the remaining series to its first term.
            return (n + 1) * logk - math.lgamma(n + 2) + np.exp(logk)
        ok = logr < 0
        ans = np.full(len(logk), np.inf)
        logg, logr = log_laplace(panels, AUX[0], taus, m)[ok], logr[ok]
        terms = [(n + 1) * logr]
        for j in range(1, n + 1):
            terms.append(j * logg - math.lgamma(j + 1) + (n + 1 - j) * logr)
        terms.append((n + 1) * logg - math.lgamma(n + 2) + np.exp(logg))
        ans[ok] = np.logaddexp.reduce(terms) - np.log(-np.expm1(logr))
        return ans

    return evaluate


def integrate_orders_reference(
    kernel,
    loss_window_eV,
    *,
    model,
    features,
    weight_length_A,
    max_order,
    hole_charge=1.0,
    primary=None,
    output_edges_eV=None,
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
    if model not in (GEOMETRIC, FACTORIAL, EXCESS) or hole_charge != 1:
        raise ValueError("a private reference model and unit hole charge are required")
    first = dict(
        epsrel=epsrel,
        area_epsabs=area_epsabs,
        spectrum_epsabs_fraction=spectrum_epsabs_fraction,
        panel_limit=panel_limit,
        max_cells=max_cells,
    )
    s = dict(
        **first,
        order_limit=order_limit,
        max_convolution_points=max_convolution_points,
        order_epsrel=order_epsrel,
        order_epsabs_fraction=order_epsabs_fraction,
        neglect_tolerance_fraction=neglect_tolerance_fraction,
    )
    window = lp.checked_window(loss_window_eV)
    count, support, geometry, lower = preflight(
        MULTIPLE_LOSS_MODEL_ID, max_order, window, primary, output_edges_eV, s
    )
    single, panels = lp.integrate_window(
        kernel,
        window,
        features=features,
        weight_length_A=weight_length_A,
        primary=primary,
        output_edges_eV=output_edges_eV,
        initial_edges=initial_edges,
        inject=inject,
        multiplier=multiplier,
        multiplier_A=multiplier_A,
        **first,
    )
    m = Q(single.multiplier_A)
    estimate = 280.5 * sum(n * (len(panels) + n - 1) for n in range(2, count + 1))
    resources = dict(
        work_lower_bound_points=lower,
        work_estimate_after_first_order=estimate,
        quadrature_points_total=0,
        cells_max=single.diagnostics.get("cells", 0),
    )
    diag = dict(resources=resources)
    if estimate > max_convolution_points:
        raise MultipleLossConvergenceError(
            "convolution_work_limit_exceeded",
            stage="after_first_order",
            resource="max_convolution_points",
            estimate=estimate,
            limit=max_convolution_points,
            diagnostics=diag,
        )
    order1 = first_order(single, panels, geometry, features, seed_errors(single), False)
    first_record = dict(
        order=1,
        area={c: order1.area[c] for c in COMPONENTS},
        area_error={c: order1.area_error[c] for c in COMPONENTS},
    )
    if geometry:
        for c in COMPONENTS:
            for key in (c, c + "_error"):
                first_record[key] = getattr(order1, key)
        first_record["outside"] = {
            s: {k: v for k, v in row.items() if not k.endswith("_parts")} for s, row in order1.outside.items()
        }
    orders = [SimpleNamespace(**first_record)] if count else []
    left_panels, right_panels = _pack(panels, model)
    env = first_envelopes(left_panels)
    right = Polynomials(right_panels, first_envelopes(right_panels), slice(None))
    for n in range(2, count + 1):
        left = Polynomials(left_panels, env, slice(None))
        left_panels, env, area, local, folded, mesh = _adapt(
            left,
            right,
            window,
            features,
            n,
            s,
            geometry,
            m**n,
            diag,
            order_inject,
            order_panel_limit or panel_limit,
            order_min_panel_width_eV,
            mass=single.primary_integral if geometry else 0.0,
            transform=lambda block, n: _transform(block, n, model),
        )
        rho = multiplier_rho(single, n)
        record = dict(order=n, area={}, area_error={}, mesh=mesh)
        for c in COMPONENTS:
            i = COLUMNS.index(c)
            record["area"][c] = area[i]
            point_error = sum(
                (Q(p.right_eV) - Q(p.left_eV)) * Q(float(e[:2, i].sum())) for p, e in zip(left_panels, env)
            )
            record["area_error"][c] = corrected(area[i], float(point_error * m**n) + local[i], rho)[0]
        if geometry:
            record["outside"] = {"low": {}, "high": {}}
            for c in COMPONENTS:
                rows = [folded[name + "." + c] for name in geometry.names]
                errors = [
                    corrected(r["value"], math.fsum(r[k] for k in ("inherited", "propagated", "local")), rho)[0]
                    for r in rows
                ]
                for key, items in ((c, [r["value"] for r in rows]), (c + "_error", errors)):
                    record[key] = np.array(items[1:-1])
                    record["outside"]["low"][key], record["outside"]["high"][key] = items[0], items[-1]
        orders.append(SimpleNamespace(**record))
    result = SimpleNamespace(reference_only=True, model_id=model, orders=tuple(orders), diagnostics=diag)
    if geometry:
        upper = (single.multiplier_interval_A or (float(m),) * 2)[1]
        tail_model = _log_tail(model, right_panels, upper)
        result.truncation = truncation(
            single,
            panels,
            geometry,
            max_order,
            count,
            support,
            None,
            neglect_tolerance_fraction,
            False,
            diag,
            tail_model=tail_model,
        )
        # At tau=0 this is the model's absolute all-axis remainder estimate.
        t = np.array([0.0])
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            logs = tail_model(t, count, log_laplace(panels, TOTAL, t, upper), log_laplace(panels, 0, t, upper))
            bound = single.primary_abs_integral * np.exp(logs[0])
        result.truncation["neglected_total_area_bound"] = float(bound) if math.isfinite(bound) else None
    return result
