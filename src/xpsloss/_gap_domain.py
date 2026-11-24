"""Pre-quadrature checks for the retarded hard-gap fixed-depth model."""

from dataclasses import replace
from decimal import Context, Decimal, ROUND_CEILING
from fractions import Fraction
import math

import numpy as np

from .dielectric import (
    DielectricGapEdgeError,
    GAP_EDGE_WARNING_FRACTION,
    _checked_static_response,
    _gap_edge_exclusion_eV,
    _sub_gap_pole,
)
from .model import DEFAULT_ABS, DEFAULT_REL, HBAR2_OVER_2ME_EV_A2, loss_function


class FixedDepthGapDomainError(ValueError):
    """The sufficient hard-gap domain conditions could not be established."""

    def __init__(self, reason, points, **details):
        super().__init__(reason)
        self.reason = reason
        self.offending_loss_eV = [p["loss_eV"] for p in points]
        self.details = {"points": points, **details}


def check_edges(material, grid):
    width = _gap_edge_exclusion_eV(material)
    indices = np.flatnonzero(abs(grid - material.gap_eV) <= width) if width else []
    if len(indices):
        exc = DielectricGapEdgeError(material.gap_eV, width, grid[indices])
        exc.offending_indices = indices.tolist()
        raise exc


# Rational operations on the exact binary64 inputs. Decimal.sqrt and ln are
# correctly rounded to nearest, independent of context rounding. Enclose each
# transcendental result by its adjacent decimal before re-entering rationals.
_CTX = Context(prec=60, rounding=ROUND_CEILING)
_PI_LO = Fraction(math.nextafter(math.pi, 0.0))


def _sqrt_interval(x):
    ctx = _CTX.copy()
    upper_x = ctx.divide(Decimal(x.numerator), Decimal(x.denominator))
    lower_x = ctx.next_minus(upper_x)
    return (Fraction(ctx.next_minus(ctx.sqrt(lower_x))), Fraction(ctx.next_plus(ctx.sqrt(upper_x))))


def _lower_q_squared(energy, loss, total_momentum):
    e, t, h = map(Fraction, (energy, loss, HBAR2_OVER_2ME_EV_A2))
    s_lo, s_hi = _sqrt_interval(e)
    r_lo, r_hi = _sqrt_interval(e - t)
    # Both lower endpoints increase with T. Factor the small difference
    # exactly as in Eq. (17), using separate enclosures of numerator and denominator.
    if total_momentum:
        return t**2 / (h * (s_hi + r_hi) ** 2)
    return (t / (s_hi + r_hi)) ** 3 * (3 * s_lo + r_lo) / (4 * e * h)


def _directed_float(x, upward):
    try:
        rounded = float(x)
    except OverflowError:
        return math.inf if x > 0 else -math.inf
    if rounded < x if upward else rounded > x:
        rounded = math.nextafter(rounded, math.inf if upward else -math.inf)
    return rounded


def _record_number(x, upward):
    rounded = _directed_float(x, upward)
    return rounded if math.isfinite(rounded) else str(_CTX.divide(Decimal(x.numerator), Decimal(x.denominator)))


def _monotonicity(material, q2_lower):
    gap = Fraction(material.gap_eV)
    rows = []
    for i, osc in enumerate(material.oscillators):
        if not osc.strength_eV2:
            continue
        center = Fraction(osc.energy_eV) + Fraction(osc.dispersion) * Fraction(HBAR2_OVER_2ME_EV_A2) * q2_lower
        lhs = (1 - (gap / center) ** 2) ** 3
        rhs = 4 * Fraction(osc.width_eV) * gap / (_PI_LO * center**2)
        rows.append(
            {
                "oscillator_index": i,
                "lhs_lower": _record_number(lhs, False),
                "rhs_upper": _record_number(rhs, True),
                "passed": center > gap and lhs > rhs,
            }
        )
    return rows


def _uniform_bounds(material):
    """Pairs (U, independent constant part, edge-log coefficient) in ascending U."""
    lines = [
        (Fraction(o.energy_eV), Fraction(o.strength_eV2), Fraction(o.width_eV))
        for o in material.oscillators
        if o.strength_eV2
    ]
    if not lines:
        return [(None, Fraction(0), Fraction(0))]
    g = Fraction(material.gap_eV)
    emin = min(e for e, _, _ in lines)
    if g >= emin:
        return []
    a, b = math.log(g), math.log(emin)
    lo, hi = math.nextafter(float(g), math.inf), math.nextafter(float(emin), 0.0)
    if lo > hi:
        return []
    candidates = sorted({min(hi, max(lo, math.exp(a + (b - a) * i / 9))) for i in range(1, 9)})
    b0 = sum(strength / e**2 for e, strength, _ in lines)
    result = []
    for u in map(Fraction, candidates):
        m = sum(strength * width / (e**2 - u**2) ** 2 for e, strength, width in lines)
        result.append((float(u), b0 / (1 - (g / u) ** 2), m * g / _PI_LO))
    return result


def _proof(candidates, gap, loss):
    if not candidates:
        return {"U_eV": None, "S_upper": None, "margin_lower": None, "proved": False}
    g, t = Fraction(gap), Fraction(loss)
    ratio = 2 * g / (g - t)
    ctx = _CTX.copy()
    argument = ctx.divide(Decimal(ratio.numerator), Decimal(ratio.denominator))
    log_upper = Fraction(ctx.next_plus(ctx.ln(argument)))
    bound, u = min((p + c * log_upper, u) for u, p, c in candidates)
    return {
        "U_eV": u,
        "S_upper": _record_number(bound, True),
        "margin_lower": _record_number(2 - bound, False),
        "proved": bound < 2,
    }


def gap_preflight(material, energy, grid, qlo, qhi, total_momentum):
    gap = material.gap_eV
    qmin = float(min(qlo))
    q2_lower = _lower_q_squared(energy, float(min(grid)), total_momentum)
    rows = _monotonicity(material, q2_lower)

    def point(i):
        return {"index": int(i), "loss_eV": float(grid[i]), "q_interval_inv_A": [float(qlo[i]), float(qhi[i])]}

    if any(not r["passed"] for r in rows):
        raise FixedDepthGapDomainError(
            "monotonicity_condition_failed",
            [point(i) for i in range(len(grid))],
            oscillators=rows,
            scope="global_minimum_q",
            q_min_inv_A=qmin,
            remedy="Choose oscillators satisfying the hard-gap monotonicity condition.",
        )
    static, tol = _checked_static_response(material, np.asarray(qmin))
    sub = np.flatnonzero(grid < gap)
    proofs = []
    if sub.size:
        candidates = _uniform_bounds(material)
        proofs = [{**point(i), **_proof(candidates, gap, float(grid[i]))} for i in sub]
        failures = [p for p in proofs if not p["proved"]]
        if failures:
            raise FixedDepthGapDomainError(
                "surface_pole_not_excluded",
                failures,
                remedy="Exclude unproved sub-gap losses or choose a supported material.",
            )

    exclusion = _gap_edge_exclusion_eV(material)
    sign = "marginal" if abs(static) <= tol else "negative"
    pole, status = None, "absent"
    if sign == "marginal":
        status = "marginal_static"
    elif exclusion:
        pole, status = _sub_gap_pole(material, np.asarray(qmin), exclusion, (DEFAULT_ABS, DEFAULT_REL, 200))
    distances = abs(grid - gap)
    near = []
    warnings = []
    if rows:
        warnings.append(
            "The hard gap gives a logarithmic divergence of Re F and logarithmic components near the edge; the unit-hole total jumps from zero."
        )
        if sub.size:
            warnings.append(
                "Below the gap the medium is lossless and the unit-hole total is zero. The separate components describe reactive contributions, not absorption."
            )
        for i in np.flatnonzero(distances <= GAP_EDGE_WARNING_FRACTION * gap):
            c = loss_function(replace(material, gap_eV=0.0), gap, qlo[i]) / math.pi
            term = c * (math.log(gap) - math.log(distances[i]))
            near.append({**point(i), "distance_eV": float(distances[i]), "edge_log_contribution": float(term)})
            warnings.append(
                f"Near the hard gap edge at loss_eV={grid[i]:.17g}, the fixed-q logarithmic term at q_minus is {term:.9g}."
            )
    return {
        "gap_eV": gap,
        "gap_edge_exclusion_eV": exclusion,
        "nearest_gap_edge_distance_eV": float(min(distances)),
        "near_gap_edge_points": near,
        "q_min_inv_A": qmin,
        "static_sign": sign,
        "static_negative_inverse_epsilon": float(static),
        "static_sign_tolerance": float(tol),
        "q_min_squared_lower_inv_A2": _record_number(q2_lower, False),
        "monotonicity": rows,
        "subgap_proofs": proofs,
        "sub_gap_pole_eV": pole,
        "sub_gap_pole_status": status,
        "epsilon_pole_q_inv_A": qmin,
        "subgap_policy": "retarded_point_values_proved_surface_delta_absent",
        "broadening_eV": 0.0,
        "real_part_method": "closed_form_D3_retained_band",
        "dielectric_model_id": "causal_oscillator_elf_v1",
        "ill_conditioned_D3_count": 0,
    }, tuple(warnings)
