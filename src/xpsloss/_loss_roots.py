"""Resolution candidates at fixed and moving momentum endpoints (M4a design 3.2)."""

import cmath
from fractions import Fraction
import math

import numpy as np
from numpy.polynomial import Polynomial as Poly


def distance(z, left, right):
    return math.hypot(max(left-z.real, 0., z.real-right), z.imag)


def _family_polynomial(lines, energy, mode, end, family, index):
    # t = 1-sqrt(1-z/E0) for both moving endpoints. At the upper
    # endpoint the momentum variable is 2-t.
    t = Poly([0, 1])
    loss = t if end == "both" else energy*t*(2-t)
    momentum = t if end == "lower" else 2-t
    shift = energy*momentum**2
    if mode == "k_parallel":
        shift -= loss**2/(4*energy)
    product = Poly([1])
    numerator = Poly([0])
    for osc in ((lines[index],) if family == "bulk" else lines):
        w = osc.energy_eV + osc.dispersion*shift
        denominator = w*w - loss*loss - 1j*osc.width_eV*loss
        if family == "bulk":
            return denominator, loss
        numerator = numerator*denominator + osc.strength_eV2*product
        product = product*denominator
    return index*product - numerator, loss


def _response(z, energy, mode, end, lines, family, index, *, factored=False):
    """Original rational equation and its z derivative on the physical branch."""
    shift = derivative = 0j
    if end != "both":
        s = cmath.sqrt(1-z/energy)
        if s == 0:
            raise ZeroDivisionError
        t = z/energy/(1+s)
        if end == "lower":
            if mode == "k_parallel":
                shift = energy*t**3*(1-t/4)
                derivative = t*t*(3-t)/(2*s)
            else:
                shift, derivative = energy*t*t, t/s
        else:
            shift, derivative = energy*(1+s)**2, -(1+s)/s
            if mode == "k_parallel":
                shift -= z*z/(4*energy)
                derivative -= z/(2*energy)
    terms, slopes = [], []
    cancelled = False
    for osc in ((lines[index],) if family == "bulk" else lines):
        w = osc.energy_eV + osc.dispersion*shift
        real_part = (w-z)*(w+z) if factored else w*w-z*z
        d = real_part - 1j*osc.width_eV*z
        dp = 2*w*osc.dispersion*derivative - 2*z - 1j*osc.width_eV
        scale = abs(w)**2 + abs(z)**2 + osc.width_eV*abs(z)
        if family == "bulk":
            return d, dp, abs(d)/scale, False
        cancelled |= abs(d) <= 1e-10*scale
        if d == 0:
            return 0j, 0j, math.inf, True
        terms.append(osc.strength_eV2/d)
        slopes.append(terms[-1]*dp/d)
    f = index - sum(terms)
    return f, sum(slopes), abs(f)/(1+sum(map(abs, terms))), cancelled


def _polish(root, neighbors, equation):
    # Keep each iterate in the original root's Voronoi cell. A Newton
    # trajectory attracted to a different simple root cannot replace this one.
    z = root
    for iteration in range(12):
        value, slope, residual, cancelled = equation(z)
        if cancelled or not slope or (iteration and residual < 8*np.finfo(float).eps):
            break
        trial = z - value/slope
        if (abs(trial-root) > (1+abs(root))/2 or
                any(abs(trial-other) < abs(trial-root) for other in neighbors)):
            z = root
            break
        if trial == z:
            break
        z = trial
    _, _, residual, cancelled = equation(z)
    return z, residual, cancelled


def _remote_single_root(t0, energy, mode, end, osc, family, index):
    """Confirm a remote, ill-conditioned one-oscillator root in local t.

    A correction can be smaller than an ulp of z while changing A/D by
    order one. Carry the real and imaginary parts of t as exact rationals.
    The 2-by-2 multiplication matrices implement complex arithmetic.
    """
    def matrix(x, y=0):
        x, y = Fraction(x), Fraction(y)
        return np.array([[x, -y], [y, x]], dtype=object)

    def scalar(m):
        return complex(float(m[0, 0]), float(m[1, 0]))

    t = matrix(t0.real, t0.imag)
    identity = matrix(1)
    e = Fraction(energy)
    target = 0 if family == "bulk" else Fraction(osc.strength_eV2)/index
    for _ in range(8):
        z = e*(t @ (2*identity-t))
        zp = 2*e*(identity-t)
        u = t if end == "lower" else 2*identity-t
        q = e*(u @ u)
        qp = (2 if end == "lower" else -2)*e*u
        if mode == "k_parallel":
            q -= (z @ z)/(4*e)
            qp -= (z @ zp)/(2*e)
        w = matrix(osc.energy_eV)+Fraction(osc.dispersion)*q
        wp = Fraction(osc.dispersion)*qp
        damping = matrix(0, osc.width_eV)
        d = (w @ w)-(z @ z)-(damping @ z)
        dp = 2*(w @ wp)-2*(z @ zp)-(damping @ zp)
        v, slope = scalar(d-matrix(target)), scalar(dp)
        scale = abs(scalar(w))**2+abs(scalar(z))**2+osc.width_eV*abs(scalar(z))
        if family == "bulk":
            residual, cancelled = abs(v)/scale, False
        else:
            s = osc.strength_eV2/scalar(d)
            residual, cancelled = abs(index-s)/(1+abs(s)), abs(scalar(d)) <= 1e-10*scale
        if residual <= 1e-8 or not slope:
            break
        step = v/slope
        t -= matrix(step.real, step.imag)
    return scalar(z), residual, cancelled


def _family_candidates(lines, energy, mode, end, family, index, window):
    left, right = window
    polynomial, loss = _family_polynomial(lines, energy, mode, end, family, index)
    if not np.isfinite(polynomial.coef).all():
        raise ArithmeticError("nonfinite polynomial coefficients")
    variables = polynomial.trim().roots()
    if not np.isfinite(variables).all():
        raise ArithmeticError("nonfinite polynomial roots")
    variables = [complex(v) for v in variables if end == "both" or v.real < 1]
    roots = [complex(loss(v)) for v in variables]
    if not all(map(cmath.isfinite, roots)):
        raise ArithmeticError("nonfinite loss roots")

    def equation(z):
        return _response(z, energy, mode, end, lines, family, index)

    for i, root in enumerate(roots):
        neighbors = roots[:i]+roots[i+1:]
        z, residual, cancelled = _polish(root, neighbors, equation)
        method = "enumerated"
        if not cancelled:
            # Check both forms where nearly coincident poles cost digits in D.
            residual = max(residual, _response(z, energy, mode, end, lines, family, index,
                                               factored=True)[2])
        if (end != "both" and len(lines) == 1 and residual > 1e-8
                and distance(root, left, right) >= right-left):
            trial, err, cancel = _remote_single_root(variables[i], energy, mode, end,
                                                     lines[0], family, index)
            if not cmath.isfinite(trial):
                raise ArithmeticError("nonfinite rationally checked root")
            if all(abs(trial-root) <= abs(trial-other) for other in neighbors):
                z, residual, cancelled = trial, err, cancel
                method = "rational_remote_check"
        if not cmath.isfinite(z):
            raise ArithmeticError("nonfinite polished root")
        if cancelled:
            continue
        if not math.isfinite(residual):
            raise ArithmeticError("nonfinite root residual")
        yield dict(z_re=z.real, z_im=z.imag, family=family, end=end, residual=residual,
                   confirmed=residual <= 1e-8, near_origin=abs(z)<left/2, method=method,
                   distance_to_window_eV=distance(z, left, right))


def candidates(material, energy, mode, window):
    lines = [o for o in material.oscillators if o.strength_eV2]
    fixed = all(o.dispersion == 0 for o in lines)
    families = [("bulk", i) for i in range(len(lines))] + [("surface", 2), ("epsilon", 1)]
    records, remote_checks, warnings = [], [], []
    all_confirmed = True
    for end in (("both",) if fixed else ("lower", "upper")):
        for family, index in families:
            try:
                with np.errstate(over="raise", divide="raise", invalid="raise"):
                    for row in _family_candidates(lines, energy, mode, end, family, index, window):
                        all_confirmed &= row["confirmed"]
                        if row["distance_to_window_eV"] >= window[1]-window[0]:
                            if row["method"] == "rational_remote_check":
                                remote_checks.append(row)
                            continue
                        z = complex(row["z_re"], row["z_im"])
                        if any(r["family"] == family and r["end"] == end
                               and abs(z-complex(r["z_re"], r["z_im"])) <= 1e-9*(1+abs(z))
                               for r in records):
                            continue
                        records.append(row)
            except (ArithmeticError, np.linalg.LinAlgError, ValueError, TypeError) as exc:
                label = f"bulk[{index}]" if family == "bulk" else family
                warning = dict(source="xps_single_loss", code="enumeration_failed", family=family,
                               end=end, reason=type(exc).__name__,
                               message=f"Root enumeration failed for {label} at {end}: {type(exc).__name__}: {exc}")
                if family == "bulk":
                    warning["active_oscillator_index"] = index
                warnings.append(warning)
    reason = None
    if warnings:
        reason = "enumeration_failed"
    elif material.gap_eV:
        reason = "hard_gap"
    elif not all_confirmed:
        reason = "unconfirmed_root"
    basis = "fixed_poles_enumerated" if fixed else "endpoint_roots_enumerated"
    return dict(candidates=records, remote_root_checks=remote_checks, warnings=warnings,
                resolution_basis="computed_candidates_only" if reason else basis, resolution_basis_reason=reason)

