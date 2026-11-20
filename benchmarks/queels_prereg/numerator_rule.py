"""Eight terms of the preregistered numerator-rule comparison (normal emission)."""

import math

import numpy as np
from scipy.special import roots_legendre

from xpsloss.model import BOHR_RADIUS_A as a0, HBAR2_OVER_2ME_EV_A2 as H

W, A, GAMMA, E = 15.6, 243.3598, 1.5, 1000.


# Preserve the preregistration's Eq. (17) quadrature limits and operation order
# independently of later changes to the package's private helper.
def _parallel_limits(energy, loss):
    """Eq. (17), factoring k±² − Ω0² before evaluating either square root."""
    s, r = math.sqrt(energy), np.sqrt(energy - loss)
    plus = s + r
    minus = loss / plus
    scale = 2 * s * math.sqrt(H)
    return (minus * np.sqrt(minus * (3*s + r)) / scale,
            plus * np.sqrt(plus * (3*s - r)) / scale)


def terms(t, a, kind="opt", alpha=1., mom="kt", panels=400, energy=E, gamma=GAMMA, amp=A,
          oscillators=None):
    o = t/(2*np.sqrt(H*energy))
    lo, hi = _parallel_limits(energy, np.asarray(t))
    edges = np.linspace(np.log(lo), np.log(hi), panels+1)
    unodes, uw = roots_legendre(16)
    half = (edges[1:] - edges[:-1])/2
    k = np.exp((edges[1:]+edges[:-1])[:,None]/2 + half[:,None]*unodes).ravel()
    weights = (half[:,None]*uw).ravel()
    q2 = k*k+(o*o if mom=="kt" else 0.)
    s = np.zeros_like(k, dtype=complex)
    s_opt = np.zeros_like(s)
    for w0, strength, width, slope in oscillators or [(W, amp, gamma, alpha)]:
        center = w0 + slope*H*q2
        det = center**2-t*t
        invdet = strength/(det**2+(width*t)**2)
        s += invdet*(det+1j*width*t)
        s_opt += invdet*((w0**2-t*t)+1j*width*t)
    bulk_s = s_opt if kind in ("opt","bulk") else s
    refl_s = s_opt if kind in ("opt","refl") else s
    b = -bulk_s
    r = refl_s/(2-refl_s)
    surf = 2*bulk_s-r
    D = k*k+o*o
    p = (o+1j*k)/D
    m = p.conj()
    u = np.exp(-a*k-1j*o*a)
    d = 1-u
    mix = u*d.conj()
    norm = abs(d)**2/D
    meas = k*k/(np.pi*a0*energy)
    pre = meas/(2*a*o)
    vals = [meas*s.imag/D, pre*(2*b.imag*(p**2*d).imag+b.real*norm),
            pre*(b*(p*d+m*mix)).real/o, pre*surf.real*norm,
            pre*(surf*m*mix).real/o, pre*(r*m**2*d.conj()).real,
            -pre*(r*p*(2j*k/D-p*u)).real, -pre*(r*p*u).real/o]
    return np.array([np.dot(weights, v) for v in vals])

def combine(v):
    e = np.sum(v[[0,1,3,5,6]])
    h = np.sum(v[[2,4,7]])
    return np.array([e+h,e,h])
