"""Normal-emission figure calculations, Simonsen, Yubero & Tougaard (1997).

Adapted from the reviewers' independent transcription of Eqs. (1), (13)–(18)
and (20). Energies are eV, lengths Å, k is in inverse Å. K is per eV per Å.
"""
from functools import lru_cache

import numpy as np
from numpy.polynomial.legendre import leggauss

H = 3.80998212  # ħ²/(2m), eV Å²
A0 = 2 * H / 14.399645  # ħ²/(m e²), with e² in eV Å
MOMENTA = ("k_parallel", "k_total")


def eq20(loss, q, width):
    recoil = H * q**2
    return 1.0 - 15.6**2 / (loss**2 - recoil**2 + 1j * width * loss)


@lru_cache(maxsize=4)
def gauss_rule(order):
    return leggauss(order)


def log_nodes(lo, hi, panels, order):
    nodes, weights = gauss_rule(order)
    edges = np.linspace(np.log(lo), np.log(hi), panels + 1)
    half = np.diff(edges)[:, None] / 2
    center = (edges[1:] + edges[:-1])[:, None] / 2
    return np.exp(center + half * nodes).ravel(), (half * weights).ravel()


def bulk(E0, W, width):
    """Infinite medium, Eq. (1). It has only one momentum argument."""
    s, r = np.sqrt(E0), np.sqrt(E0 - W)
    q, weights = log_nodes((s - r) / np.sqrt(H), (s + r) / np.sqrt(H), 1000, 40)
    return np.sum(weights * (-1.0 / eq20(W, q, width)).imag) / (E0 * np.pi * A0)


def split(E0, W, depth, width, momentum, panels=96):
    """Eqs. (13)–(17), positive losses, unit hole charge, 0° from the normal.

    Return (electron, hole), each scalar or an array over creation depths.
    The last axis inside the integral is always the parallel wave number.
    """
    s, r_final = np.sqrt(E0), np.sqrt(E0 - W)
    lo = np.sqrt((s - r_final)**3 * (3*s + r_final) / (4*E0*H))
    hi = np.sqrt((s + r_final)**3 * (3*s - r_final) / (4*E0*H))
    q, weights = log_nodes(lo, hi, panels, 16)
    weights *= q  # dq = q d(log q)
    a = np.asarray(depth)[..., None]
    wperp = Om = W / (2 * np.sqrt(H * E0))
    if momentum == "k_parallel":
        eps = eq20(W, q, width)
    elif momentum == "k_total":
        eps = eq20(W, np.sqrt(q**2 + Om**2), width)
    else:
        raise ValueError(f"unknown momentum: {momentum}")
    inv1 = 1.0 / eps - 1.0
    N = Om**2 + q**2
    Dp, Dm = Om + 1j*q, Om - 1j*q
    hf = 1j / wperp
    P = 1.0 / (2 * np.pi * E0 * A0 * a * wperp)
    r = 0.0  # v_parallel / v_perp at normal emission
    ex = np.exp
    out = {}
    out["13a"] = (1 / (np.pi * A0 * E0), 1j * q * inv1 / N)
    brace13 = ((1j + r) * Dp ** 2 * (ex(-1j * Om * a - q * a) - 1)
               + (-1j + r) * Dm ** 2 * (ex(1j * Om * a - q * a) - 1)
               + (1j - r) * (ex(-1j * Om * a) - ex(-q * a)) * (ex(1j * Om * a) - ex(-q * a)) * N)
    out["13b"] = (-P, 1j * q * inv1 / N ** 2 * brace13)
    brace14 = ((1j + r) * Dp * (1 - ex(-q * a - 1j * Om * a))
               + (-1j + r) * Dm * (ex(-2 * q * a) - ex(-q * a - 1j * Om * a)))
    out["14"] = (-P, hf * q * inv1 / N * brace14)
    C1 = (eps - 1) * (eps + 2) / (eps * (eps + 1))
    pre15 = 1j * q * (-1j + r) * Dm * (ex(-q * a) - ex(-1j * Om * a)) / N
    out["15e"] = (-P, pre15 * (C1 * (ex(1j * Om * a) - ex(-q * a)) * Dp / N
                              + (eps - 1) / (eps + 1) * Dm / N * ex(1j * Om * a)))
    out["15h"] = (-P, pre15 * C1 * (-1j * hf * ex(-q * a)))
    pre16 = q * (1 - eps) / (1 + eps) * (1 - 1j * r) * Dp * ex(-1j * Om * a) / N
    out["16e"] = (P, pre16 * (2j * q * ex(1j * Om * a) - Dp * ex(-q * a)) / N)
    out["16h"] = (P, pre16 * (-1j * hf * ex(-q * a)))
    def integrate(label):
        prefactor, integrand = out[label]
        return np.sum(weights * (prefactor * integrand).real, axis=-1)

    electron = sum(integrate(term) for term in ("13a", "13b", "15e", "16e"))
    hole = sum(integrate(term) for term in ("14", "15h", "16h"))
    return electron, hole


def depth_average(W, width, momentum):
    """Fig. 8, Eq. (18), K_sc = integral t exp(-t) K_eff(λt) dt."""
    edges = np.array([0, .01, .04, .12, .3, .7, 1.4, 2.5, 4.5, 8, 14, 25, 40])
    nodes, weights = gauss_rule(6)
    half = np.diff(edges)[:, None] / 2
    center = (edges[1:] + edges[:-1])[:, None] / 2
    t = (center + half * nodes).ravel()
    wt = (half * weights).ravel() * t * np.exp(-t)
    electron, hole = split(1413.5, W, 25.4*t, width, momentum, panels=48)
    return np.sum(wt * electron), np.sum(wt * hole)
