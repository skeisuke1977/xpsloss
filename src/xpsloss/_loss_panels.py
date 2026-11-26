"""Loss-energy quadrature and shared error-budget driver (M4a design section 3)."""

from dataclasses import dataclass
from fractions import Fraction as Q
from functools import cached_property
import math
import operator

import numpy as np
from numpy.polynomial.chebyshev import chebvander
from numpy.polynomial.legendre import leggauss

from . import _fold
from ._loss_roots import distance

# Clearance rule (design 3.2 (4)): a panel needs dist(z, P) >= |P|/2 for every
# candidate. The width cap Gamma_min is the boundary case of this rule for a
# pole at distance Gamma_min/2, which binary64 roots reach only up to their
# last bits. Roots within a relative 2**-48 inside the boundary count as on it,
# so the panel count does not depend on those bits (prereg 11, CE5 and CE6).
CLEARANCE_FACTOR = 1. - 2.**-48
from ._single_loss_types import (
    BinnedPrimary, COLUMNS, COMPONENTS, DeltaPrimary, EPS, LossWindowAreaResult,
    PARTS, QUANTITIES, SingleLossConvergenceError, SingleLossInputError,
    SingleLossSpectrumResult, inverse_area, vector,
)
from .model import finite_number

VERIFY = -np.cos((np.arange(16)+.5)*np.pi/16)
GAUSS, GAUSS_W = leggauss(40)
GAUSS_POLY = chebvander(GAUSS, 16).T
VERIFY_POLY = chebvander(VERIFY, 16).T
# Clenshaw-Curtis weights from the exact integrals of the Chebyshev basis.
CC = np.array([0. if k % 2 else 2/(1-k*k) for k in range(17)]) @ _fold.CARDINAL
# Map values at the nine even nodes to Chebyshev coefficients.
CARDINAL9 = np.cos(np.pi*np.outer(np.arange(9), np.arange(8, -1, -1))/8)/4
CARDINAL9[[0, -1], :] *= .5
CARDINAL9[:, [0, -1]] *= .5


def _midpoint(a, b):
    return a + (b-a)/2 if b > np.finfo(float).max-a else (a+b)/2


def _nodes(a, b):
    mid, half = _midpoint(a, b), (b-a)/2
    nodes = np.empty(17)
    nodes[0], nodes[-1] = a, b
    nodes[1:-1] = mid + half*_fold.LOBATTO[1:-1]
    check = mid + half*VERIFY
    separated = np.all((nodes[:-1] < check) & (check < nodes[1:]))
    return nodes, check, separated


def checked_window(window, energy=None):
    details = {}
    try:
        a, b = (finite_number(x, "loss_window_eV") for x in window)
        details["loss_window_eV"] = (a, b)
        valid = 0 < a < b and (energy is None or b < energy) and _nodes(a, b)[2]
    except (TypeError, ValueError, OverflowError) as exc:
        details["detail"] = str(exc)
        valid = False
    if not valid:
        raise SingleLossInputError("invalid_window", loss_window_eV=window, message_details=details)
    return a, b


def validate_settings(settings):
    for key, value in settings.items():
        if key.endswith("limit") or key == "max_cells":
            if isinstance(value, (bool, np.bool_)):
                raise SingleLossInputError("invalid_settings", setting=key)
            try:
                settings[key] = operator.index(value)
            except TypeError as exc:
                raise SingleLossInputError("invalid_settings", setting=key) from exc
            ok = settings[key] >= 1
        else:
            try:
                ok = finite_number(value, key) > 0
            except ValueError:
                ok = False
        if not ok or (key == "epsrel" and value < 64*EPS):
            raise SingleLossInputError("invalid_settings", setting=key)


def _width_failure(a, b, stage, diag):
    exc = SingleLossConvergenceError if stage else SingleLossInputError
    raise exc("panel_width_unresolved", panel_eV=(a, b), panel_width_ulp=(b-a)/math.ulp(a),
              **({"diagnostics": diag} if stage else {}))


def _initial_mesh(window, features, limit, diag, *, energy=None, gap=0., initial_edges=None):
    a, b = map(Q, window)
    knots = {a, b}
    if initial_edges is not None:
        supplied = list(map(Q, initial_edges))
        if not supplied or supplied[0] != a or supplied[-1] != b or supplied != sorted(set(supplied)):
            raise SingleLossInputError("invalid_settings", setting="initial_edges")
        knots.update(supplied)
        width, poles = b-a, []
    else:
        for key in ("jumps_eV", "kinks_eV"):
            knots.update(Q(v) for v in features.get(key, ()) if a < v < b)
        g = Q(gap)
        step = a-g
        while g+step < b:
            knots.add(g+step)
            step *= 2
        if energy is not None:
            e = Q(energy)
            step = e-b
            while e-step > a:
                knots.add(e-step)
                step *= 2
        width = b-a
        for key in ("min_half_width_eV", "min_variation_scale_eV"):
            v = features.get(key)
            if v is not None:
                if not math.isfinite(v) or v <= 0:
                    raise SingleLossInputError("invalid_settings", setting=key)
                width = min(width, 2*Q(v))
        poles = [complex(x, y) for x, y in features.get("singularities", ())]
    knots = sorted(knots)
    diag["gradings"] = list(map(float, knots))
    intervals = [(l, r, math.ceil((r-l)/width)) for l, r in zip(knots, knots[1:])]
    count = sum(n for _, _, n in intervals)
    mesh = []
    for left, right, n in intervals:
        step = (right-left)/n
        if count <= limit:
            indices = range(n)
        else:
            # The base count needs no large allocation. Only panels near a
            # candidate can add further bisections to the exact required count.
            nearby = set()
            for pole in poles:
                i = math.floor((Q(pole.real)-left)/step)
                nearby.update(range(max(0, i-2), min(n, i+3)))
            indices = sorted(nearby)
        for i in indices:
            stack = [(float(left+i*step), float(left+(i+1)*step))]
            while stack:
                l, r = stack.pop()
                if any(distance(z, l, r) < (r-l)/2*CLEARANCE_FACTOR for z in poles):
                    mid = _midpoint(l, r)
                    if not l < mid < r:
                        _width_failure(l, r, 0, diag)
                    diag["tolerance_bisected_panels"].append((l, r))
                    count += 1
                    stack.extend(((mid, r), (l, mid)))
                elif count <= limit:
                    mesh.append((l, r))
    diag["required_panels"] = count
    if count > limit:
        raise SingleLossInputError("resolution_exceeds_panel_limit", required_panels=count,
                                   message_details=dict(panel_limit=limit))
    return mesh


@dataclass
class Panel:
    left_eV: float
    right_eV: float
    values: np.ndarray
    errors: np.ndarray
    coeff_diff: np.ndarray
    interpolation: np.ndarray
    verification: np.ndarray
    rho: np.ndarray
    epsilon: np.ndarray

    @classmethod
    def from_lobatto(cls, left_eV, right_eV, values, errors):
        v, e = _checked_arrays((values, errors), 17)
        return cls(left_eV, right_eV, v.copy(), e.copy(), *(np.zeros(11) for _ in range(5)))

    @cached_property
    def coefficients(self):
        return self.values @ _fold.CARDINAL.T

    @cached_property
    def component_bounds(self):
        coeff = abs(self.coefficients[_fold.COMP_INDEX])
        errors = self.errors[_fold.COMP_INDEX].max(axis=1)
        # Lobatto Lebesgue constant < 4. Include coefficient roundoff
        # and nodal uncertainty in both the amplitude and derivative.
        uncertainty = 4*errors + 64*EPS*abs(self.values[_fold.COMP_INDEX]).max(axis=1)
        return coeff.sum(axis=1)+uncertainty, coeff @ (np.arange(17)**2)+256*uncertainty

    def check(self, values, errors):
        coeff = self.coefficients
        diff = coeff.copy()
        diff[:, :9] -= self.values[:, ::2] @ CARDINAL9.T
        self.coeff_diff = abs(diff).sum(axis=1)
        self.interpolation = (abs(diff @ GAUSS_POLY) @ GAUSS_W)*((self.right_eV-self.left_eV)/2)
        residual = abs(values-coeff @ VERIFY_POLY)
        self.rho = residual.max(axis=1)
        self.epsilon = errors.max(axis=1)
        self.verification = (residual+errors).max(axis=1)


def _checked_arrays(pair, n):
    v, e = (np.asarray(x, dtype=float) for x in pair)
    if (v.shape != (11, n) or e.shape != v.shape or not np.isfinite(v).all()
            or not np.isfinite(e).all() or np.any(e < 0)):
        raise ArithmeticError("nonfinite or invalid kernel outputs")
    return v, e


def _areas(panels):
    contributions, estimates = [], []
    for p in panels:
        width = p.right_eV-p.left_eV
        weights = CC*(width/2)
        contributions.append(p.values @ weights)
        estimates.append(np.stack((p.interpolation, width*p.verification,
                                   p.errors @ weights, np.zeros(11),
                                   12*EPS*(abs(p.values) @ weights)), axis=1))
    errors = np.stack(estimates)
    outputs, local = {}, {}
    for j, name in enumerate(COLUMNS):
        key = "window_area."+name
        outputs[key] = _fold._entry(math.fsum(c[j] for c in contributions), errors[:, j].sum(axis=0))
        local[key] = errors[:, j].sum(axis=1)
    return outputs, local


def _diagnostics(settings):
    return dict(resolution_basis=None, resolution_basis_reason="analytic_layer", candidates=[], remote_root_checks=[],
                required_panels=0, final_panels=0, stages=0, gradings=[], tolerance_bisected_panels=[],
                verification=dict(rho_max=0., epsilon_max=0.), error_budget=dict(outputs={}, max_ratio=None),
                balance={}, kernel_calls=dict(calls=0, points_sent_per_stage=[], points_sent_total=0,
                                             distinct_points=0, batch_sizes=[]),
                effective_settings=settings, comparison={}, warnings=[], LbK_b={})


def _panel_keys(nodes):
    return [(float(t), "right" if i == 0 else "left" if i == 16 else "interior")
            for i, t in enumerate(nodes)]


class _Sampler:
    def __init__(self, kernel, diag, physical):
        self.kernel, self.diag, self.physical = kernel, diag, physical
        self.cache = {}

    def key(self, pair):
        return (pair[0], "interior") if self.physical else pair

    def fill(self, pairs, stage):
        keys = sorted({self.key(p) for p in pairs}.difference(self.cache))
        if not keys:
            return
        points = np.array([t for t, _ in keys])
        if self.physical:
            pair = self.kernel(points, stage)
        else:
            calls = self.diag["kernel_calls"]
            calls["calls"] += 1
            calls["points_sent_per_stage"].append(len(keys))
            calls["points_sent_total"] += len(keys)
            calls["batch_sizes"].append(len(keys))
            pair = self.kernel(points, np.array([s for _, s in keys]))
        v, e = _checked_arrays(pair, len(keys))
        self.cache.update((k, (v[:, i], e[:, i])) for i, k in enumerate(keys))
        self.diag["kernel_calls"]["distinct_points"] = len({t for t, _ in self.cache})

    def read(self, keys, role, inject):
        v, e = (np.stack(x, axis=1) for x in zip(*(self.cache[self.key(k)] for k in keys)))
        if inject is not None:
            t = np.array([p for p, _ in keys])
            return _checked_arrays(inject(role, t, v, e), len(keys))
        return v, e


def _budget(outputs, local, diag, widths, support, epsrel, area_epsabs, spec_abs):
    records, worst, split = {}, None, np.zeros(len(widths), dtype=bool)
    for key, row in outputs.items():
        budget = max(area_epsabs if key.startswith("window_area.") else spec_abs, epsrel*abs(row["value"]))
        estimate = row["estimate"]
        if not all(math.isfinite(v) for v in (*row.values(), budget)):
            raise ArithmeticError("nonfinite loss-window result")
        ratio = estimate/budget if budget else (math.inf if estimate else None)
        records[key] = {p: row[p] for p in (*PARTS, "estimate")}
        records[key].update(budget=budget, ratio=ratio,
                            flag="zero_budget_zero_error" if ratio is None else None)
        if ratio is not None and (worst is None or ratio > worst["value"]):
            worst = dict(key=key, value=ratio)
        if estimate > budget:
            spans = support.get(key, widths)
            apportioned = .5*budget*(spans/spans.sum())
            contributions = local[key] if key in local else estimate*(widths/widths.sum())
            split |= contributions > apportioned
    diag["error_budget"] = dict(outputs=records, max_ratio=worst)
    return split


def _finish(window, length, charge, panels, areas, folded, geometry, primary, edges, diag, energy, mode):
    budgets = diag["error_budget"]["outputs"]
    area = {k: areas["window_area."+k]["value"] for k in COLUMNS}
    errors = {k: budgets["window_area."+k]["estimate"] for k in COLUMNS}
    a, b = window
    v, e = panels[0].values[:, 0], panels[0].errors[:, 0]
    def scaled(scale, vals):
        return {k: float(scale*Q(x)) for k, x in zip(COLUMNS, vals)}
    h = scaled(Q(length)*Q(a), v)
    he = scaled(Q(length)*Q(a), e)
    he = {k: he[k]+2*EPS*abs(h[k]) for k in COLUMNS}
    diag["LbK_b"] = scaled(Q(length)*Q(b), panels[-1].values[:, -1])
    w = length*area["total"]
    inv, interval, reason = inverse_area(area["total"], errors["total"])
    params = dict(electron_energy_eV=energy, loss_window_eV=window, weight_length_A=length,
                  epsilon_momentum=mode, theta_deg=0., hole_charge=charge, window_area=area,
                  window_area_error=errors, window_area_budget={k: budgets["window_area."+k]["budget"]
                                                               for k in COLUMNS},
                  single_loss_weight=w, single_loss_weight_error=length*errors["total"]+EPS*abs(w),
                  inverse_window_area_A=inv, inverse_window_area_interval_A=interval, inverse_reason=reason,
                  cutoff_sensitivity=h, cutoff_sensitivity_error=he, diagnostics=diag,
                  warnings=tuple(row["message"] for row in diag["warnings"]))
    area_rounding = {k: areas["window_area."+k]["rounding"] for k in COLUMNS}
    _balance(diag, "window_area.components", [area["no_hole"], charge*area["hole"]], area["total"],
             area_rounding["total"]+area_rounding["no_hole"]+abs(charge)*area_rounding["hole"])
    if primary is None:
        return LossWindowAreaResult(**params)
    mass = math.fsum(primary.intensities)
    outside = {"low": {}, "high": {}}
    outputs = folded.outputs
    params.update(output_edges_eV=edges, primary=primary, primary_integral=mass, primary_abs_integral=mass,
                  outside=outside, cutoff_sensitivity_bins={}, cutoff_sensitivity_bins_error={})
    for name in QUANTITIES:
        rows = [outputs[r+"."+name] for r in geometry.names]
        for side, row in zip(outside, (rows[0], rows[-1])):
            outside[side][name] = row["value"]
            outside[side][name+"_error"] = row["estimate"]
        bins = rows[1:-1]
        vals = np.array([r["value"] for r in bins])
        err = np.array([r["estimate"] for r in bins])
        params[name], params[name+"_error"] = vals, err
        widths = [Q(r)-Q(l) for l, r in zip(edges, edges[1:])]
        params[name+"_per_eV"] = np.array([float(Q(v)/d) for v, d in zip(vals, widths)])
        params[name+"_per_eV_error"] = np.array([float((Q(e)+Q(EPS)*abs(Q(v)))/d)
                                                 for v, e, d in zip(vals, err, widths)])
        rounding = math.fsum(r["rounding"]+r["edge"] for r in rows)
        if name == "zero_loss":
            expected = mass
        elif name == "spectrum":
            expected = mass*(1+w)
            rounding += float(Q(mass)*Q(length)*Q(area_rounding["total"]))
        else:
            component = name.removeprefix("single_loss_")
            expected = float(Q(mass)*Q(length)*Q(area[component]))
            rounding += float(Q(mass)*Q(length)*Q(area_rounding[component]))
        _balance(diag, name, [r["value"] for r in rows], expected, rounding)
    for region in geometry.names:
        total, electron, hole = [outputs[region+".single_loss_"+c] for c in COMPONENTS]
        rounding = math.fsum(factor*(r["edge"]+r["rounding"])
                             for factor, r in ((1., total), (1., electron), (abs(charge), hole)))
        _balance(diag, region+".components", [electron["value"], charge*hole["value"]], total["value"], rounding)
    for sign in ("plus", "minus"):
        for suffix, attr in (("", "value"), ("_error", "estimate")):
            outname = "cutoff_sensitivity"+suffix
            values = {}
            for c in COMPONENTS:
                rows = [outputs[f"{region}.cutoff_sensitivity.{sign}.{c}"][attr] for region in geometry.names]
                values[c] = np.array(rows[1:-1])
                for side, val in zip(outside, (rows[0], rows[-1])):
                    outside[side].setdefault(outname, {}).setdefault(sign, {})[c] = val
            params["cutoff_sensitivity_bins"+suffix][sign] = values
    return SingleLossSpectrumResult(**params)


def _balance(diag, key, terms, target, rounding):
    # Include the roundoff accumulated inside each signed integral.
    diag["balance"][key] = dict(difference=math.fsum([*terms, -target]),
                                rounding_bound=rounding+8*EPS*(math.fsum(map(abs, terms))+abs(target)))


def integrate_window(kernel, loss_window_eV, *, features, weight_length_A, hole_charge=1.,
                     primary=None, output_edges_eV=None, initial_edges=None, inject=None,
                     epsrel, area_epsabs, spectrum_epsabs_fraction, panel_limit, max_cells,
                     multiplier="weight_length", multiplier_A=None):
    settings = dict(epsrel=epsrel, area_epsabs=area_epsabs, spectrum_epsabs_fraction=spectrum_epsabs_fraction,
                    panel_limit=panel_limit, max_cells=max_cells)
    return _integrate(kernel, checked_window(loss_window_eV), features, weight_length_A, hole_charge,
                      primary, output_edges_eV, settings, _diagnostics(settings),
                      kernel_errstate=np.geterr(), initial_edges=initial_edges, inject=inject,
                      multiplier=multiplier, multiplier_A=multiplier_A)


def _prepare(window, features, primary, edges, settings, diag, energy, gap, initial_edges):
    geometry = None
    if primary is not None:
        if not isinstance(primary, (DeltaPrimary, BinnedPrimary)):
            raise SingleLossInputError("invalid_primary")
        edges = vector(edges, "invalid_output_edges", increasing=True)
        geometry = _fold.Geometry(primary, edges)
        support = geometry.points
        lo, hi = (support[0], support[-1]) if support else (None, None)
        diag["primary_support"] = dict(min_eV=lo, max_eV=hi,
                                       min_minus_E0_eV=None if energy is None or lo is None else lo-energy,
                                       max_minus_E0_eV=None if energy is None or hi is None else hi-energy)
    mesh = _initial_mesh(window, features, settings["panel_limit"], diag, energy=energy, gap=gap,
                         initial_edges=initial_edges)
    return geometry, edges, mesh


def preflight(prepared, settings, diag):
    geometry, _, mesh = prepared
    if geometry:
        bounds = [mesh[0][0], *(r for _, r in mesh)]
        count = geometry.cell_count(bounds, settings["max_cells"])
        if count > settings["max_cells"]:
            raise SingleLossInputError("cell_limit_exceeded", message_details=dict(
                cells_at_least=count, max_cells=settings["max_cells"]))
    for l, r in mesh:
        if not _nodes(l, r)[2]:
            _width_failure(l, r, 0, diag)


def _denominator(row, extra, epsrel, name, stage, diag):
    """Decide on the exact binary64 inputs, before rounding the reciprocal."""
    area, error = Q(row["value"]), Q(row["estimate"])
    d = area + Q(extra)
    threshold = Q(epsrel)*abs(d)
    resolved = d-error > 0 and error <= threshold
    def display(value):
        # A sum or a threshold may exceed binary64 even for a finite reciprocal.
        try:
            return float(value)
        except OverflowError:
            return str(value)
    diag["denominator"] = dict(value=display(d), error=float(error), threshold=display(threshold),
                               resolved=resolved, stage=stage)
    argument = "weight_length_kind" if name == "weight_length" else "multiplier"
    attr = dict(stage=stage, diagnostics=diag, message_details=dict(argument=argument))
    if d+error < 0 and error <= threshold:
        raise SingleLossInputError(name+"_undefined", **attr)
    if not resolved:
        return None, None
    center, interval, reason = inverse_area(d, error)
    if reason:
        raise SingleLossConvergenceError(name+"_not_representable", **attr)
    return center, interval


@np.errstate(over="raise", divide="raise", invalid="raise")
def _integrate(kernel, window, features, length, charge, primary, edges, settings, diag, *,
               kernel_errstate, energy=None, mode=None, gap=0., zero=False, initial_edges=None, inject=None,
               prepared=None, seed=False, extra=0., multiplier="weight_length", multiplier_A=None,
               weight_interval=None):
    validate_settings(settings)
    geometry, edges, mesh = prepared or _prepare(window, features, primary, edges, settings, diag, energy,
                                                 gap, window if zero else initial_edges)
    m = multiplier_A if multiplier == "external" else length
    interval = weight_interval if multiplier == "weight_length" else None
    denominator_key = "window_area.term_13a" if seed else "window_area.total"
    name = "weight_length" if seed else "multiplier"
    stage = 0
    sampler = _Sampler(kernel, diag, energy is not None)
    previous = {}
    try:
        while True:
            bounds = [mesh[0][0], *(p[1] for p in mesh)]
            if geometry:
                count = geometry.cell_count(bounds, settings["max_cells"])
                if count > settings["max_cells"]:
                    details = dict(cells_at_least=count, max_cells=settings["max_cells"])
                    if stage:
                        raise SingleLossConvergenceError("cell_limit_reached", diagnostics=diag, message_details=details)
                    raise SingleLossInputError("cell_limit_exceeded", message_details=details)
                if not stage:
                    diag["initial_cells"] = count
                diag["cells"] = count
            nodes = {}
            requests = []
            for l, r in mesh:
                if (l, r) in previous:
                    continue
                inter, ver, separated = _nodes(l, r)
                if not separated:
                    _width_failure(l, r, stage, diag)
                inter_keys = _panel_keys(inter)
                ver_keys = [(float(t), "interior") for t in ver]
                nodes[l, r] = inter_keys, ver_keys
                requests.extend(inter_keys+ver_keys)
            diag.update(final_panels=len(mesh), stages=stage)
            lower_key = (window[0], "left")
            if geometry:
                requests.append(lower_key)
            if not zero:
                with np.errstate(**kernel_errstate):
                    sampler.fill(requests, stage)
            current = {}
            for bounds_pair in mesh:
                if bounds_pair in previous:
                    current[bounds_pair] = previous[bounds_pair]
                    continue
                inter, ver = nodes[bounds_pair]
                pair = (np.zeros((11, 17)),)*2 if zero else sampler.read(inter, "interpolation", inject)
                p = Panel.from_lobatto(*bounds_pair, *pair)
                if not zero:
                    p.check(*sampler.read(ver, "verification", inject))
                current[bounds_pair] = p
            panels = list(current.values())
            diag["verification"] = dict(rho_max=max(float(p.rho.max()) for p in panels),
                                        epsilon_max=max(float(p.epsilon.max()) for p in panels))
            outputs, local = _areas(panels)
            pending = False
            if seed or multiplier == "inverse_window_area":
                if seed:
                    outputs = {denominator_key: outputs[denominator_key]}
                m, interval = _denominator(outputs[denominator_key], extra if seed else 0.,
                                           settings["epsrel"], name, stage, diag)
                pending = m is None
            folded, spec_abs, support = None, 0., {}
            if geometry and not pending:
                if not zero and inject is not None:
                    lower_minus = tuple(x[:, 0] for x in sampler.read([lower_key], "interpolation", inject))
                else:
                    lower_minus = None if zero else sampler.cache[sampler.key(lower_key)]
                folded = _fold._fold(panels, geometry, m, lower_minus)
                support = folded.support_widths
                for region in geometry.names:
                    for quantity in QUANTITIES[:-1]:
                        key = region+"."+quantity
                        outputs[key] = folded.outputs[key]
                local.update(folded.panel_errors)
                spec_abs = settings["spectrum_epsabs_fraction"]*math.fsum(primary.intensities)
            split = _budget(outputs, local, diag, np.diff(bounds), support, settings["epsrel"],
                            settings["area_epsabs"], spec_abs)
            if pending:
                index = max(range(len(mesh)), key=lambda i: (local[denominator_key][i],
                            Q(mesh[i][1])-Q(mesh[i][0]), -i))
                split.fill(False)
                split[index] = True
            if not np.any(split):
                if seed:
                    row = outputs[denominator_key]
                    return dict(area13=row["value"], error13=row["estimate"], length_A=m,
                                interval_A=interval), panels
                result = _finish(window, m, charge, panels, outputs, folded, geometry,
                                 primary, edges, diag, energy, mode)
                result.weight_length_A = length
                result.multiplier_A = m
                result.multiplier_convention = multiplier
                result.multiplier_interval_A = interval
                if interval is not None:
                    diag["multiplier"] = dict(convention=multiplier, value_A=m, interval_A=interval)
                    _fold.propagate_multiplier(result)
                _check_finite_result(result)
                return result, panels
            if len(mesh)+int(split.sum()) > settings["panel_limit"]:
                budget = diag["error_budget"]
                if pending:
                    argument = "weight_length_kind" if seed else "multiplier"
                    raise SingleLossConvergenceError(name+"_unresolved", stage=stage, diagnostics=diag,
                                                     message_details=dict(argument=argument, worst=budget["max_ratio"]))
                key = budget["max_ratio"]["key"]
                worst = budget["outputs"][key]
                details = dict(output=key, estimate=worst["estimate"], budget=worst["budget"],
                               ratio=worst["ratio"], panel_limit=settings["panel_limit"])
                if worst["kernel_point_error"] > worst["budget"]:
                    raise SingleLossConvergenceError("kernel_error_exceeds_budget", diagnostics=diag,
                                                     message_details=details)
                raise SingleLossConvergenceError("loss_budget_not_met", diagnostics=diag,
                                                 dominant_part=max(PARTS, key=worst.get), message_details=details)
            mesh = []
            for (l, r), refine in zip(current, split):
                if refine:
                    mid = _midpoint(l, r)
                    mesh.extend(((l, mid), (mid, r)))
                else:
                    mesh.append((l, r))
            previous = current
            stage += 1
    except (ArithmeticError, FloatingPointError) as exc:
        raise SingleLossConvergenceError("nonfinite_output", diagnostics=diag) from exc


def _check_finite_result(result):
    def check(obj):
        if isinstance(obj, dict):
            for value in obj.values():
                check(value)
        elif isinstance(obj, (float, np.ndarray)) and not np.isfinite(obj).all():
            raise ArithmeticError("derived quantity outside binary64")
    for name, value in vars(result).items():
        if name != "diagnostics":
            check(value)
    check(result.diagnostics["balance"])
    check(result.diagnostics["LbK_b"])

