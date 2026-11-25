"""Normal-emission single-scattering cross section, Eq. (18) of Simonsen et al. (1997)."""

from dataclasses import dataclass
import math

import numpy as np

from . import _gap_response, xps_fixed_depth as fixed
from ._depth_kernel import assemble_terms, depth_averaged_weights
from ._gap_domain import check_edges, gap_preflight
from .dielectric import _checked_static_response
from .model import (
    BOHR_RADIUS_A, DECIMAL_DIGITS, FLOAT_EPSILON as EPS, HBAR2_OVER_2ME_EV_A2,
    Material, NumericalConvergenceError, _PanelLimitError, _accuracy, finite_number,
)

DEPTH_AVERAGED_MODEL_ID = "simonsen1997_depth_average_external_v0"
COLUMNS = (*fixed.TERM_NAMES, "no_hole", "hole", "total")
_BASE_WARNINGS = (
    "This cross section describes the first single scattering, including intrinsic losses.",
    "Simonsen et al. (1997) caution that the REELS path weighting requires examination for XPS.",
    "The weight length is supplied externally and is constant over all paths.",
)
_LOW_LOSS = ("The low-loss behavior is unverified. A loss-window area may depend strongly on "
             "its lower cutoff. Only point values are returned.")


@dataclass(frozen=True)
class DepthAveragedResult:
    """All cross sections in eV^-1 angstrom^-1 per solid path, including vacuum losses.

    The path weight is x*exp(-x/L)/L^2, with a positive external length L.
    """
    electron_energy_eV: float
    weight_length_A: float
    theta_deg: float
    hole_charge: float
    epsilon_momentum: str
    loss_eV: np.ndarray
    term_13a: np.ndarray
    term_13b: np.ndarray
    term_14: np.ndarray
    term_15_T1: np.ndarray
    term_15_T2: np.ndarray
    term_15_T3: np.ndarray
    term_16_e: np.ndarray
    term_16_h: np.ndarray
    no_hole: np.ndarray
    hole: np.ndarray
    total: np.ndarray
    quadrature_error: np.ndarray
    component_errors: dict
    diagnostics: dict
    warnings: tuple


class DepthAverageConvergenceError(NumericalConvergenceError):
    def __init__(self, reason, loss_eV, index, *, column=None, estimate=None, budget=None,
                 subdivisions=0, panels=0, inner=None):
        self.reason, self.loss_eV, self.index = reason, float(loss_eV), int(index)
        self.column, self.estimate, self.budget = column, estimate, budget
        self.subdivisions, self.panels, self.inner = int(subdivisions), int(panels), inner
        message = f"{reason} at loss_eV={loss_eV:.17g}, index={index}"
        if column is not None:
            message += f", {column} estimate={estimate:.9g}, budget={budget:.9g}"
        super().__init__(message if inner is None else f"{message}: {inner}")


class _WeightRangeError(ArithmeticError):
    pass


def _reason(exc):
    if isinstance(exc, _WeightRangeError):
        return "weight_out_of_range"
    if isinstance(exc, _gap_response._D3ConditionError):
        return "d3_ill_conditioned"
    message = str(exc)
    for words, reason in (
        (("surface denominator", "surface pole at", "surface coefficient unresolved"), "surface_denominator_unresolved"),
        (("oscillator denominator",), "oscillator_denominator_unresolved"),
        (("linewidth unresolved", "surface oscillator sum unresolved", "unresolved oscillator and surface lines"),
         "linewidth_unresolved"),
        (("surface resonance unresolved", "surface search"), "surface_resonance_unresolved"),
        (("retained band", "removed-band", "damping unresolved"), "coefficient_unresolved"),
    ):
        if any(word in message for word in words):
            return reason
    return "response_out_of_range"


@dataclass(frozen=True)
class _Point:
    material: Material
    energy: float
    loss: float
    index: int
    b: float
    offset: float
    p: float
    charge: float
    direct: bool
    dispersive: bool


def _product_error(z, ez, w, ew):
    """Finite perturbation of one Cartesian product, including the mixed term."""
    return abs(w)*ez + (abs(z) + ez)*ew + EPS*abs(z*w)


def _project(z, ez, w):
    x, y, ex, ey = w
    return (_product_error(z.real, ez[0], x, ex) + _product_error(z.imag, ez[1], y, ey)
            + EPS*(abs(z.real*x) + abs(z.imag*y)))


def _term_errors(coeff, weights, measure, bulk_weight, b, rel_measure):
    (s, r, c), (ds, dr, dc) = coeff
    n, w13 = weights["N"], weights["W13"]
    # 13b has two independent Cartesian projections, not a complex dot product.
    error13 = (2*_product_error(s.imag, ds[1], w13[1], w13[3])
               + _product_error(s.real, ds[0], n[0], n[2]))
    errors = (
        bulk_weight*(ds[1] + (2*rel_measure + 10*EPS)*abs(s.imag)),
        measure*error13, measure*_project(s, ds, weights["W14"])/b,
        measure*_product_error(c.real, dc[0], n[0], n[2]),
        measure*_project(c, dc, weights["W15b"])/b, measure*_project(r, dr, weights["W15c"]),
        measure*_project(r, dr, weights["W16e"]), measure*_project(r, dr, weights["W16h"])/b,
    )
    return np.stack(errors, axis=1)


def _densities(point, k, rel_z, left, right, allowance):
    m, b = point.material, point.b
    q = np.hypot(k, point.offset) if point.offset else k
    q_eval = q if point.dispersive else np.zeros(1)
    coeff = fixed._coefficients(m, point.loss, q_eval, relative_z=rel_z if m.gap_eV and point.dispersive else 0.)
    changes = None
    if point.dispersive and not m.gap_eV:
        changes = fixed._include_node_roundoff(m, point.loss, q, rel_z, coeff)
    try:
        weights = depth_averaged_weights(k, b, point.p, e_k=rel_z*k/2, e_b=3*EPS*b, e_p=EPS*point.p)
    except FloatingPointError as exc:
        raise _WeightRangeError(str(exc)) from exc
    measure = k**2/(math.pi*BOHR_RADIUS_A*point.energy)
    bulk_weight = measure/(k*k + b*b)
    rel_measure = rel_z + 4*EPS

    def project_errors():
        return _term_errors(coeff, weights, measure, bulk_weight, b, rel_measure)

    errors = project_errors()
    if allowance is not None:
        mask = (errors > allowance[:, None]).any(axis=1)
        if mask.any():
            if point.dispersive:
                squares = fixed._node_squares(left, right, mask, point.loss, point.energy, bool(point.offset))
                fixed._exact_points(m, point.loss, q, mask, coeff, node_squares=squares)
                ds, dr = changes
                residual = 4*EPS*(right-left) + 10.0**(4-DECIMAL_DIGITS)
                ratio = (residual/rel_z[:, 0])[np.nonzero(mask)[0]]
                for (er, ei), change in zip(coeff[1], (ds, dr, 2*ds + dr)):
                    er[mask] += ratio*change[mask]
                    ei[mask] += ratio*change[mask]
            else:
                fixed._exact_points(m, point.loss, q_eval, np.ones(1, bool), coeff)
            errors = project_errors()
    def complex_weight(name):
        x, y = weights[name][:2]
        return x + 1j*y

    terms = assemble_terms(*coeff[0], measure, b, bulk_weight, 2*weights["W13"][1], weights["N"][0],
                           *(complex_weight(name) for name in ("W14", "W15b", "W15c", "W16e", "W16h")))
    errors += (rel_measure[:, None, :] + 8*EPS)*abs(terms)
    direct_floor = None
    if point.direct:
        # For any common k,b,p and any common real S,R the exact sum is zero.
        # Bound the operations AFTER choosing those common inputs, including C=2S-R.
        # Replace |kL-b^2| by kL+b^2 before cancellation (M3A_DESIGN 3.1).
        # kL+b^2 <= sqrt(D H) by Cauchy-Schwarz. Its scaled bound is
        # |W15c|, already available without reconstructing an unscaled power.
        unsigned_w15c = np.hypot(*weights["W15c"][:2])
        sr, rr = abs(coeff[0][0].real), abs(coeff[0][1].real)
        norm = weights["N"][0]
        scale = measure*(sr*(norm + weights["W14"][0]/b)
                         + (2*sr + rr)*(norm + abs(weights["W15b"][0])/b)
                         + rr*(unsigned_w15c + abs(weights["W16e"][0]) + weights["W16h"][0]/b))
        # The longest product/quotient path plus the projection and seven additions.
        direct_floor = 24*EPS*scale
    return terms, errors, direct_floor


@dataclass
class _PanelDetails:
    floors: np.ndarray | None = None
    direct_values: np.ndarray | None = None
    direct_floor: np.ndarray | None = None


def _panel_integrals(left, right, *, point, details, allowance=None):
    """Evaluate one loss on the supplied log panels; return two (n_panels, 11) arrays."""
    half = (right-left)/2
    mid = (left+right)/2
    k = np.exp(mid[:, None] + half[:, None]*fixed._NODES)
    rel_z = fixed.MOMENTUM_ULPS*EPS*(2 + np.maximum(abs(left), abs(right))[:, None])
    terms, arithmetic, direct_floor = _densities(point, k, rel_z, left, right, allowance)

    def integrate(array, high):
        rule = fixed._HIGH_RULE if high else fixed._LOW_RULE
        split = len(fixed._LOW_RULE[0])
        nodes = array[..., split:] if high else array[..., :split]
        return half[:, None]*np.einsum("ptn,n->pt", nodes, rule[1])

    def estimates(v, e):
        hi = integrate(v, True)
        floor = fixed.PANEL_ROUNDOFF*integrate(abs(v), True)
        arith = 2*integrate(e, True) + integrate(e, False)
        return hi, np.maximum(abs(hi-integrate(v, False)), floor) + arith, floor + arith

    values, errors, floors = estimates(terms, arithmetic)
    charge = point.charge
    def combine(e):
        groups = fixed._component_sums(e, abs(charge))
        groups += 8*EPS*fixed._component_sums(abs(values), abs(charge))
        return np.column_stack((e, groups))

    all_values = np.column_stack((values, fixed._component_sums(values, charge)))
    all_errors = combine(errors)
    details.floors = combine(floors)
    if point.direct:
        # Independent reduction at EACH node, before integration or grouping.
        node_sum = terms.sum(axis=1, keepdims=True)
        direct, err, floor = estimates(node_sum, direct_floor[:, None, :])
        details.direct_values = direct[:, 0]
        details.direct_floor = floor[:, 0]
        all_errors[:, -1] = err[:, 0]
        details.floors[:, -1] = details.direct_floor
    return all_values, all_errors


def _aggregate(values, errors, details, charge):
    summed = np.array([math.fsum(v) for v in values[:, :8].T])
    output = np.r_[summed, fixed._component_sums(summed, charge)]
    # fsum is rounded once. Charge products and group additions have their own floor.
    rounding = EPS*abs(values)
    rounding[:, 8:] += 8*EPS*fixed._component_sums(abs(values[:, :8]), abs(charge))
    direct = None
    if details.direct_values is not None:
        rounding[:, -1] = EPS*abs(details.direct_values)
    contributions = errors + rounding
    estimates = contributions.sum(axis=0)*(1 + len(values)*EPS)
    if details.direct_values is not None:
        k_dir = math.fsum(details.direct_values)
        delta = abs(output[-1] - k_dir)
        e_dir = float(estimates[-1])
        estimates[-1] = e_dir + delta
        direct = dict(K_dir=k_dir, e_dir=e_dir, delta_rec=float(delta), e_total=float(estimates[-1]))
    return output, estimates, direct, contributions


@dataclass
class _Mesh:
    point: _Point
    edges: np.ndarray
    initial_count: int
    allowance: np.ndarray | None = None

    def evaluate(self, epsabs, epsrel):
        details = _PanelDetails()
        with np.errstate(over="raise", divide="raise", invalid="raise", under="ignore"):
            try:
                values, errors = _panel_integrals(self.edges[:-1], self.edges[1:], point=self.point,
                                                 details=details, allowance=self.allowance)
                self.output, self.errors, self.direct, self.contributions = _aggregate(
                    values, errors, details, self.point.charge)
                self.budgets = _accuracy(self.output, epsabs, epsrel)
                self.floors = details.floors.sum(axis=0)
                if not np.all(np.isfinite([self.output, self.errors, self.floors, self.budgets])):
                    raise ArithmeticError("depth-averaged response or error estimate outside numerical range")
                # An arbitrarily small positive budget can overflow a *ranking*
                # ratio even when the response, error and floor are all finite.
                with np.errstate(over="ignore"):
                    self.priority = (self.contributions/self.budgets).max(axis=1)
                    if self.direct is not None:
                        local = abs(values[:, -1] - details.direct_values)
                        remaining = max(self.budgets[-1] - self.direct["delta_rec"],
                                        EPS*self.budgets[-1], math.ulp(0.))
                        self.priority = np.maximum(self.priority, (self.contributions[:, -1] + local)/remaining)
            except (ArithmeticError, NumericalConvergenceError) as exc:
                raise self.failure(_reason(exc), inner=str(exc)) from exc

    def failure(self, reason, column=None, estimates=None, inner=None):
        return DepthAverageConvergenceError(
            reason, self.point.loss, self.point.index, inner=inner, panels=len(self.edges)-1,
            subdivisions=len(self.edges)-1-self.initial_count,
            column=COLUMNS[column] if column is not None else None,
            estimate=float(estimates[column]) if column is not None else None,
            budget=float(self.budgets[column]) if column is not None else None,
        )

    def refine_coefficients(self, epsabs, epsrel):
        # Replace only gapless nodes whose projected uncertainty could use the
        # remaining budget. No unrefined, reducible bound may trigger a floor failure.
        if self.point.material.gap_eV or np.all(self.errors <= self.budgets):
            return
        share = self.budgets[:8].copy()
        charge = abs(self.point.charge)
        for group, indices, scale in ((8, fixed._ELECTRON, 1.), (9, fixed._HOLE, 1.),
                                      (10, fixed._ELECTRON, 1.), (10, fixed._HOLE, charge)):
            if scale:
                share[indices] = np.minimum(share[indices], self.budgets[group]/(8*scale))
        self.allowance = share/(8*(self.edges[-1]-self.edges[0]))
        self.evaluate(epsabs, epsrel)
        if np.any(self.floors > self.budgets):
            # Rule out any remaining reducible part before labeling a floor.
            self.allowance.fill(0.)
            self.evaluate(epsabs, epsrel)

    def check_floor(self):
        if self.direct is not None and self.floors[-1] > self.budgets[-1]:
            raise self.failure("direct_sum_floor_exceeds_budget", -1, self.floors)
        bad = np.flatnonzero(self.floors > self.budgets)
        if bad.size:
            raise self.failure("panel_floor_exceeds_budget", int(bad[0]), self.floors)

    def adapt(self, epsabs, epsrel, limit):
        while np.any(self.errors > self.budgets):
            if len(self.edges)-1 >= limit:
                j = int(np.argmax(self.errors/self.budgets))
                raise self.failure("budget_not_met", j, self.errors)
            i = int(np.argmax(self.priority))
            midpoint = (self.edges[i] + self.edges[i+1])/2
            self.edges = np.insert(self.edges, i+1, midpoint)
            self.evaluate(epsabs, epsrel)



def _solve(material, electron_energy_eV, loss_eV, length, mode, theta, charge,
           epsabs, epsrel, limit, *, panels=None, adapt=True):
    if not isinstance(material, Material):
        raise ValueError("material must be a Material")
    if not isinstance(mode, str) or mode not in fixed.EPSILON_MOMENTA:
        raise ValueError("epsilon_momentum must be 'k_parallel' or 'k_total'")
    if finite_number(theta, "theta_deg") != 0:
        raise ValueError("only normal emission (theta_deg=0) is supported")
    length = finite_number(length, "weight_length_A")
    energy = finite_number(electron_energy_eV, "electron_energy_eV")
    charge = finite_number(charge, "hole_charge")
    if length <= 0 or energy <= 0:
        raise ValueError("weight_length_A and electron_energy_eV must be positive")
    grid = np.atleast_1d(loss_eV)
    if grid.ndim != 1 or not grid.size or grid.dtype.kind not in "fiu":
        raise ValueError("loss_eV must be a real scalar or a nonempty one-dimensional grid")
    grid = np.array(grid, dtype=float)
    if not np.all((grid > 0) & (grid < energy)):
        raise ValueError("require 0 < loss_eV < electron_energy_eV")
    epsabs, epsrel = finite_number(epsabs, "epsabs"), finite_number(epsrel, "epsrel")
    if epsabs <= 0 or epsrel <= 0:
        raise ValueError("epsabs and epsrel must be positive")
    if isinstance(limit, bool) or not isinstance(limit, int) or limit < 1:
        raise ValueError("limit must be a positive integer")

    diagnostics, warnings = {}, ()
    if material.gap_eV:
        check_edges(material, grid)
    with np.errstate(all="ignore"):
        lower, upper = fixed._parallel_limits(energy, grid)
        b = grid/(2*math.sqrt(HBAR2_OVER_2ME_EV_A2*energy))
        offsets = b if mode == "k_total" else np.zeros_like(b)
        qlo, qhi = np.hypot(lower, offsets), np.hypot(upper, offsets)
        valid = (lower > 0) & (lower < upper) & np.isfinite(upper) & (b > 0) & np.isfinite(b)
        if material.gap_eV:
            valid &= (qlo < qhi) & np.isfinite(qhi)
    if not valid.all():
        i = int(np.flatnonzero(~valid)[0])
        raise DepthAverageConvergenceError("momentum_bounds_unresolved", grid[i], i,
                                          inner="parallel momentum bounds unresolved")
    try:
        if material.gap_eV:
            diagnostics, warnings = gap_preflight(material, energy, grid, qlo, qhi, mode == "k_total")
        else:
            _checked_static_response(material, qlo.min())
    except (NumericalConvergenceError, ArithmeticError) as exc:
        i = int(qlo.argmin())
        raise DepthAverageConvergenceError(_reason(exc), grid[i], i, inner=str(exc)) from exc

    active = any(o.strength_eV2 for o in material.oscillators)
    dispersive = any(o.strength_eV2 and o.dispersion for o in material.oscillators)
    meshes = []
    for i, t in enumerate(grid):
        point = _Point(material, energy, float(t), i, b[i], offsets[i], 1/length, charge,
                       bool(active and t < material.gap_eV and charge == 1), dispersive)
        try:
            if panels is None:
                edges = fixed._log_breakpoints(material, t, lower[i], upper[i], offsets[i])
                if len(edges)-1 > limit:
                    raise _PanelLimitError(f"integration limit too small for {len(edges)-1} initial oscillator panels")
                if material.gap_eV:
                    edges = _gap_response.resolve_surface(material, t, edges, offsets[i], limit)
            else:
                edges = np.asarray(panels, dtype=float)
                if (len(grid) != 1 or edges.ndim != 1 or len(edges) < 2 or not np.all(np.diff(edges) > 0)
                        or edges[0] != math.log(lower[i]) or edges[-1] != math.log(upper[i])):
                    raise ValueError("panels must partition the full log momentum interval")
        except _PanelLimitError as exc:
            raise ValueError(f"depth average at loss_eV={t:.17g}, index=({i},): {exc}") from exc
        except (NumericalConvergenceError, ArithmeticError) as exc:
            raise DepthAverageConvergenceError(_reason(exc), t, i, inner=str(exc)) from exc
        meshes.append(_Mesh(point, edges, len(edges)-1))
    # Stage barriers apply to the original grid, even when it crosses the gap.
    # Each observable panel call belongs to exactly one original loss index.
    for mesh in meshes:
        mesh.evaluate(epsabs, epsrel)
        mesh.refine_coefficients(epsabs, epsrel)
    for mesh in meshes:
        mesh.check_floor()
    if adapt:
        for mesh in meshes:
            mesh.adapt(epsabs, epsrel, limit)
    return grid, meshes, diagnostics, _BASE_WARNINGS + warnings + ((_LOW_LOSS,) if active else ())


def _depth_averaged_point(material, electron_energy_eV, loss_eV, *, weight_length_A,
                          epsilon_momentum, hole_charge, epsabs, epsrel, limit, panels=None, adapt=True):
    """Return values, estimates (both in COLUMNS order), and final log-panel endpoints."""
    _, (mesh,), _, _ = _solve(material, electron_energy_eV, loss_eV, weight_length_A, epsilon_momentum,
                              0., hole_charge, epsabs, epsrel, limit, panels=panels, adapt=adapt)
    return mesh.output, mesh.errors, mesh.edges


def effective_cross_section_depth_averaged(
    material, electron_energy_eV, loss_eV, *, weight_length_A, epsilon_momentum,
    theta_deg=0.0, hole_charge=1.0, epsabs=fixed.FIXED_ABS, epsrel=fixed.FIXED_REL, limit=200,
) -> DepthAveragedResult:
    """Eq. (18) with weight x*exp(-x/L)/L^2 and external constant L in angstrom.

    Returns eV^-1 angstrom^-1 per unit path in the solid, including vacuum losses.
    The angle is measured from the surface normal. Hole terms use unit positive charge.
    """
    grid, meshes, diag, warnings = _solve(material, electron_energy_eV, loss_eV, weight_length_A,
                                          epsilon_momentum, theta_deg, hole_charge, epsabs, epsrel, limit)
    outputs = np.stack([m.output for m in meshes], axis=1)
    errors = np.stack([m.errors for m in meshes], axis=1)
    ratios = errors/_accuracy(outputs, epsabs, epsrel)
    j, i = np.unravel_index(ratios.argmax(), ratios.shape)
    direct = [m.direct for m in meshes]
    used = [d for d in direct if d is not None]
    diag.update(effective_settings=dict(epsabs=float(epsabs), epsrel=float(epsrel), limit=limit),
                direct_sum=direct, direct_sum_max={k: max(d[k] for d in used)
                                                  for k in ("e_dir", "delta_rec")} if used else None,
                error_budget_ratios=dict(zip(COLUMNS, ratios)),
                max_budget_ratio=dict(value=float(ratios[j, i]), column=COLUMNS[j], index=int(i)))
    if material.gap_eV:
        magnitude = abs(outputs[:8]).sum(axis=0)
        mask = (grid < material.gap_eV) & (magnitude > 0)
        residual = abs(outputs[8, mask] + outputs[9, mask])/magnitude[mask]
        point_residuals = [None]*len(grid)
        for index, value in zip(np.flatnonzero(mask), residual):
            point_residuals[index] = float(value)
        diag.update(ill_conditioned_D3_count=0,
                    subgap_identity_residual=point_residuals,
                    max_subgap_identity_residual=float(residual.max()) if mask.any() else None)
    return DepthAveragedResult(float(electron_energy_eV), float(weight_length_A), float(theta_deg),
                               float(hole_charge), epsilon_momentum, grid, *outputs, errors.max(axis=0),
                               dict(zip(COLUMNS, errors)), diag, warnings)
