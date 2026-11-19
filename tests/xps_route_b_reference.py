"""Independent second route for the fixed-depth XPS cross section (test I3 of
docs/M2_FIRST_CASE.md).

Source: A. Cohen Simonsen, F. Yubero, S. Tougaard, Phys. Rev. B 56, 1612 (1997).
This module does NOT use the closed forms Eqs. (13)-(16).  It evaluates Eq. (2)

    K_eff = 2/((2 pi)^4 x hbar^2 w) Int dt Int d^3r rho_e(r,t)
            Re[ i Int d^3k (k.v) Phi_ind(k,w) exp(i k.r - i w t) ]

with the charges of Eq. (3), the Fourier transforms of Eq. (6) (time integrals done
analytically), the F functions of Eqs. (8)-(10) taken from their k_perp-integral
DEFINITIONS (integrated numerically), the induced potentials of Eqs. (11)-(12), and
d^3k = 2 pi k_par dk_par dk_perp with k_perp integrated NUMERICALLY over (-inf, inf)
(the extension stated after Eq. (12)).  The k_par limits are those of Eq. (17).

Scope: theta = 0, hbar w > 0 (the pi delta(w) parts vanish), gapless oscillator
materials, eps evaluated at k_par (the constraint eps(k, w) = eps(k_par, w) implied by
Eq. (7); the k_perp integral is only defined for an eps that does not depend on k_perp).

Conventions: energies eV, lengths Angstrom, hbar = 1 (times in hbar/eV), Gaussian units.
hbar v = 2 (H E0)^(1/2) (non-relativistic E0 = m v^2/2), H = hbar^2/2m.
The only constants are xpsloss.model.HBAR2_OVER_2ME_EV_A2 (H) and BOHR_RADIUS_A (a0);
e^2 = 2 H / a0 so that this route and the closed forms share one constant source.

Time-integral convergence: the charge transforms use w -> w + i0 (causal, as in Eq. (9):
Int_0^inf exp(i w t) dt = pi delta(w) + i/w); the electron's probe factor for t > t_s uses
exp(-eta t), eta -> 0+, i.e. 1/(v(k_perp - Omega0) + i0).

Terms returned (unit hole charge for the hole terms):
  term_13a  bulk part of the electron/image term: the a -> inf limit of the k_perp
            integral, where (k v) S_M rho_e -> 2 pi e^2 a Omega0/v delta(k_perp - Omega0)
  term_13b  (electron+image term from Eq. (2)) - term_13a
  term_14   (1/eps - 1) rho_hh'         with the medium probe factor S_M
  term_15_T1  C1 F^M_ee'                with S_M,  C1 = (eps-1)(eps+2)/(eps(eps+1))
  term_15_T2  C1 F^M_hh'                with S_M
  term_15_T3  C2 F^V_ee'                with S_M,  C2 = (1-eps)/(1+eps)
  term_16_e   C2 (F^M_ee' + F^V_ee')    with the vacuum probe factor S_V
  term_16_h   C2 F^M_hh'                with S_V
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.laguerre import laggauss

from xpsloss.model import BOHR_RADIUS_A, HBAR2_OVER_2ME_EV_A2

H = HBAR2_OVER_2ME_EV_A2
A0 = BOHR_RADIUS_A
E2 = 2.0 * H / A0            # e^2 [eV A]
ECH = np.sqrt(E2)

TERMS = ("term_13a", "term_13b", "term_14", "term_15_T1", "term_15_T2",
         "term_15_T3", "term_16_e", "term_16_h")
NO_HOLE_TERMS = ("term_13a", "term_13b", "term_15_T1", "term_15_T3", "term_16_e")
HOLE_TERMS = ("term_14", "term_15_T2", "term_16_h")

_XG, _WG = leggauss(20)
_XL, _WL = laggauss(80)
_XS, _WS = leggauss(32)
_XS, _WS = 0.5 * (_XS + 1.0), 0.5 * _WS          # nodes on (0, 1)
_XO, _WO = leggauss(16)


# ----------------------------------------------------------------------------- input

def inverse_epsilon_kpar(material, loss_eV, q):
    """1/eps(q, w) for a gapless oscillator material, e^{-i w t} convention:
    -1/eps = -1 + sum_i A_i / (Omega_i(q)^2 - w^2 - i G_i w),
    Omega_i(q) = E_i + alpha_i H q^2."""
    if material.gap_eV != 0.0:
        raise ValueError("reference route implemented for gapless materials only")
    q = np.asarray(q, dtype=float)
    s = np.zeros(q.shape, dtype=complex)
    w = float(loss_eV)
    for o in material.oscillators:
        om = o.energy_eV + o.dispersion * H * q ** 2
        s += o.strength_eV2 / (om ** 2 - w ** 2 - 1j * o.width_eV * w)
    return 1.0 - s


def kpar_limits_theta0(E0, w):
    """Eq. (17) at theta = 0, written without cancellation:
    k_par-^2 = (s-r)^3 (3s+r)/(4 s^2 H), k_par+^2 = (s+r)^3 (3s-r)/(4 s^2 H)."""
    s, r = np.sqrt(E0), np.sqrt(E0 - w)
    d = w / (s + r)                                   # s - r
    lo = np.sqrt(d ** 3 * (3 * s + r) / (4 * s * s * H))
    hi = np.sqrt((s + r) ** 3 * (3 * s - r) / (4 * s * s * H))
    return lo, hi


# ----------------------------------------------------------------------------- k_perp

def _phi(z):
    """(exp(z) - 1)/z, stable near z = 0."""
    out = np.empty_like(z)
    small = np.abs(z) < 1e-3
    zs = z[small]
    out[small] = 1 + zs / 2 + zs ** 2 / 6 + zs ** 3 / 24 + zs ** 4 / 120
    zl = z[~small]
    out[~small] = np.expm1(zl) / zl
    return out


def _mul(d1, d2):
    out = {}
    for n1, c1 in d1.items():
        for n2, c2 in d2.items():
            out[n1 + n2] = out.get(n1 + n2, 0) + c1 * c2
    return out


def _add(d1, d2):
    out = dict(d1)
    for n, c in d2.items():
        out[n] = out.get(n, 0) + c
    return out


class KPerp:
    """k_perp integrals at theta = 0 for one (E0, w, a).  All quantities do not depend
    on the material, because eps depends on k_par only.

    Charge transforms (Eq. (6) with Eq. (3); a = (-a, 0) so k.a = -k_perp a; t_s = a/v):
      rho^M_ee'(k) = -e[exp(+i k a) t_s phi(i(Q-k)a) + exp(-i k a) t_s phi(i(Q+k)a)]
      rho^M_hh'(k) = e (exp(i k a) + exp(-i k a)) * (i/w)          (unit hole charge)
      rho^V_ee'(k) = -i e exp(iQa) [1/(v(Q-k)+i0) + 1/(v(Q+k)+i0)]
    Probe factors of Eq. (2), Int dt Int d^3r rho_e exp(i k.r - i w t):
      S_M(k) = -e exp(-i k a) t_s phi(i(k-Q)a)                       (0 < t < t_s)
      S_V(k) = -i e exp(-iQa) / (v(k-Q) + i0)                         (t > t_s)
    with k = k_perp, Q = Omega0 = w/v, phi(z) = (e^z - 1)/z.
    """

    def __init__(self, E0, w, a):
        self.V = 2.0 * np.sqrt(H * E0)          # hbar v [eV A]
        self.w = float(w)
        self.a = float(a)
        self.Q = self.w / self.V
        self.ts = self.a / self.V

    # -- stable (core) forms
    def _SM(self, k):
        a, Q = self.a, self.Q
        return -ECH * np.exp(-1j * k * a) * self.ts * _phi(1j * (k - Q) * a)

    def _rho_ee(self, k):
        a, Q = self.a, self.Q
        return (-ECH * np.exp(1j * k * a) * self.ts * _phi(1j * (Q - k) * a)
                - ECH * np.exp(-1j * k * a) * self.ts * _phi(1j * (Q + k) * a))

    def _rho_hh(self, k):
        return ECH * 2.0 * np.cos(k * self.a) * (1j / self.w)

    # -- tail forms: {n: c_n(k)} meaning sum_n c_n(k) exp(i n a k); valid for |Re k| > Q
    def _SM_t(self, k):
        e, V, Q, a = ECH, self.V, self.Q, self.a
        return {0: -e * np.exp(-1j * Q * a) / (1j * V * (k - Q)), -1: e / (1j * V * (k - Q))}

    def _rho_ee_t(self, k):
        e, V, Q, a = ECH, self.V, self.Q, self.a
        d1 = {0: -e * np.exp(1j * Q * a) / (1j * V * (Q - k)), 1: e / (1j * V * (Q - k))}
        d2 = {0: -e * np.exp(1j * Q * a) / (1j * V * (Q + k)), -1: e / (1j * V * (Q + k))}
        return _add(d1, d2)

    def _rho_hh_t(self, k):
        c = ECH * 1j / self.w * np.ones_like(k)
        return {1: c, -1: c}

    # -- quadrature pieces
    def _core(self, q, K, extra=()):
        """Composite Gauss-Legendre nodes on [-K, K]; panels graded towards the poles
        +-i q and no wider than pi/a (resolves exp(+-2 i k a))."""
        hmax = np.pi / self.a
        e = [0.0]
        x = 0.0
        while x < K:
            x = min(x + min(0.25 * np.hypot(x, q), hmax), K)
            e.append(x)
        e = np.asarray(e)
        e = np.concatenate([-e[::-1], e[1:]])
        if extra:
            e = np.unique(np.concatenate([e, np.asarray(extra, dtype=float)]))
        lo, hi = e[:-1], e[1:]
        k = (0.5 * (hi - lo)[:, None] * _XG + 0.5 * (hi + lo)[:, None]).ravel()
        wts = (0.5 * (hi - lo)[:, None] * _WG).ravel()
        return k, wts

    def _tail0(self, f, K, q):
        """Int_{|k|>K} f(k) dk for a non-oscillatory f decaying at least like 1/k^2:
        log-spaced panels (ratio 2) up to Kb = 8 max(q, K), then k = Kb/s."""
        Kb = 8.0 * max(q, K)
        n = max(1, int(np.ceil(np.log2(Kb / K))))
        e = K * (Kb / K) ** (np.arange(n + 1) / n)
        lo, hi = e[:-1], e[1:]
        k = (0.5 * (hi - lo)[:, None] * _XG + 0.5 * (hi + lo)[:, None]).ravel()
        wts = (0.5 * (hi - lo)[:, None] * _WG).ravel()
        tot = 0.0 + 0.0j
        for side in (1.0, -1.0):
            tot += np.sum(wts * f(side * k))
            ks = side * Kb / _XS
            tot += np.sum(_WS * f(ks) * Kb / _XS ** 2)
        return tot

    def _tail_points(self, K, q, keys):
        """Nodes, weights and bookkeeping for _tail (all evaluated in one call)."""
        Kb = 8.0 * max(q, K)
        n = max(1, int(np.ceil(np.log2(Kb / K))))
        e = K * (Kb / K) ** (np.arange(n + 1) / n)
        lo, hi = e[:-1], e[1:]
        kl = (0.5 * (hi - lo)[:, None] * _XG + 0.5 * (hi + lo)[:, None]).ravel()
        wl = (0.5 * (hi - lo)[:, None] * _WG).ravel()
        ks = Kb / _XS
        wsub = _WS * Kb / _XS ** 2
        k0 = np.concatenate([kl, ks, -kl, -ks]).astype(complex)
        w0 = np.concatenate([wl, wsub, wl, wsub]).astype(complex)
        pts = [k0]
        blocks = [(0, slice(0, k0.size), w0)]
        start = k0.size
        a = self.a
        for nn in keys:
            if nn == 0:
                continue
            b = nn * a
            y = _XL / abs(b)
            wy = _WL / abs(b)
            sgn = 1.0 if nn > 0 else -1.0
            for side in (1.0, -1.0):
                kk = side * K + 1j * sgn * y
                # right tail: +i*sgn ; left tail: -i*sgn ; exp(i b k) = exp(i b side K) exp(-|b| y)
                fac = (1j * sgn if side > 0 else -1j * sgn) * np.exp(1j * b * side * K)
                pts.append(kk)
                blocks.append((nn, slice(start, start + kk.size), wy * fac))
                start += kk.size
        return np.concatenate(pts), blocks

    def _tail(self, dict_fun, K, q):
        """Int_{|k|>K} sum_n c_n(k) exp(i n a k) dk.  n = 0: log-spaced panels up to
        8 max(q, K) and k = Kb/s beyond; n != 0: the path is turned onto Re k = +-K (into
        the half plane where exp(i n a k) decays; the only singularities, k = +-i q and
        k = +-Q, have |Re k| < K) and integrated with Gauss-Laguerre in |n| a Im k."""
        keys = list(dict_fun(np.array([K + 0j])).keys())
        kk, blocks = self._tail_points(K, q, keys)
        d = dict_fun(kk)
        tot = 0.0 + 0.0j
        for nn, sl, wts in blocks:
            if nn in d:
                tot += np.sum(wts * d[nn][sl])
        return tot

    def _pv(self, g, c, K, q, k, wts):
        """P Int_R g(k)/(k - c) dk with |c| < K: subtraction on [-K, K] (nodes k, wts with
        c as a panel edge) plus tails."""
        core = np.sum(wts * (g(k) - g(c)) / (k - c)) + g(c) * np.log((K - c) / (K + c))
        return core + self._tail0(lambda kk: g(kk) / (kk - c), K, q)

    def integrals(self, q):
        """k_perp integrals at k_par = q (complex numbers):
        I_ee = Int (k v) S_M 4pi/k^2 rho^M_ee'   I_hh = Int (k v) S_M 4pi/k^2 rho^M_hh'
        I_S  = Int (k v) S_M 4pi/k^2             I_V  = Int (k v) S_V 4pi/k^2
        F_ee, F_hh, F_V: (q/2pi) Int rho/(k^2+q^2) for rho^M_ee', rho^M_hh', rho^V_ee'."""
        Q, V, a = self.Q, self.V, self.a
        # Split point K of the k_perp axis (> Q, so the tail forms are regular).  For
        # q a >= 50 the turned paths pass the poles +-i q only where exp(-|n| a Im k) <
        # e^-50, so K need not exceed q; otherwise K is kept well beyond q.
        if q * a >= 50.0:
            K = max(2.0 * Q, 10.0 / a)
        else:
            K = 2.0 * max(q, Q) + 10.0 / a
        k, wts = self._core(q, K, extra=(Q, -Q))
        k2 = k * k + q * q
        base = k * V * self._SM(k) * 4 * np.pi / k2
        r_ee = self._rho_ee(k)
        r_hh = self._rho_hh(k)
        out = dict(
            I_ee=np.sum(wts * base * r_ee), I_hh=np.sum(wts * base * r_hh), I_S=np.sum(wts * base),
            F_ee=q / (2 * np.pi) * np.sum(wts * r_ee / k2),
            F_hh=q / (2 * np.pi) * np.sum(wts * r_hh / k2))

        def base_t(kk):
            return _mul({0: kk * V * 4 * np.pi / (kk * kk + q * q)}, self._SM_t(kk))
        inv = lambda kk: {0: 1.0 / (kk * kk + q * q)}
        out["I_ee"] += self._tail(lambda kk: _mul(base_t(kk), self._rho_ee_t(kk)), K, q)
        out["I_hh"] += self._tail(lambda kk: _mul(base_t(kk), self._rho_hh_t(kk)), K, q)
        out["I_S"] += self._tail(base_t, K, q)
        out["F_ee"] += q / (2 * np.pi) * self._tail(lambda kk: _mul(inv(kk), self._rho_ee_t(kk)), K, q)
        out["F_hh"] += q / (2 * np.pi) * self._tail(lambda kk: _mul(inv(kk), self._rho_hh_t(kk)), K, q)

        # vacuum: 1/(x + i0) = P(1/x) - i pi delta(x)
        gV = lambda kk: kk * V * 4 * np.pi / (kk * kk + q * q)
        J = self._pv(gV, Q, K, q, k, wts) - 1j * np.pi * gV(Q)
        out["I_V"] = -1j * ECH * np.exp(-1j * Q * a) / V * J
        h = lambda kk: 1.0 / (kk * kk + q * q)
        JF = -self._pv(h, Q, K, q, k, wts) + self._pv(h, -Q, K, q, k, wts) - 1j * np.pi * (h(Q) + h(-Q))
        out["F_V"] = q / (2 * np.pi) * (-1j * ECH * np.exp(1j * Q * a) / V) * JF
        return out


# ----------------------------------------------------------------------------- k_par

def _resonances(materials, w):
    """k_par positions of narrow lines for dispersive oscillators: bulk Omega(q) = w and
    the single-oscillator surface condition Omega(q)^2 = w^2 + A/2; half-widths G/(4 alpha H q)."""
    pts = []
    for m in materials:
        for o in m.oscillators:
            if o.dispersion <= 0:
                continue
            for target in (w, np.sqrt(w * w + o.strength_eV2 / 2.0)):
                if target > o.energy_eV:
                    qr = np.sqrt((target - o.energy_eV) / (o.dispersion * H))
                    pts.append((qr, o.width_eV / (4 * o.dispersion * H * qr)))
    return pts


def kpar_nodes(E0, w, materials, n_base=12):
    """Composite 16-point Gauss-Legendre in u = ln k_par between the Eq. (17) limits,
    with extra panel edges clustered around the narrow lines of dispersive inputs."""
    lo, hi = kpar_limits_theta0(E0, w)
    ue = list(np.linspace(np.log(lo), np.log(hi), n_base + 1))
    for qr, hw in _resonances(materials, w):
        for f in (0, 0.25, 0.5, 1, 2, 4, 8, 16):
            for sgn in (1, -1):
                qq = qr + sgn * f * hw
                if lo < qq < hi:
                    ue.append(np.log(qq))
    ue = np.unique(np.asarray(ue))
    u0, u1 = ue[:-1], ue[1:]
    u = (0.5 * (u1 - u0)[:, None] * _XO + 0.5 * (u1 + u0)[:, None]).ravel()
    wu = (0.5 * (u1 - u0)[:, None] * _WO).ravel()
    q = np.exp(u)
    return q, wu * q


def reference_terms(materials, E0, losses_eV, a):
    """Route-2 terms for several gapless materials at theta = 0, depth a [A].
    Returns a list (one dict per material) of arrays over losses_eV, keys TERMS plus
    'no_hole' and 'hole' (unit hole charge).  Units eV^-1 A^-1."""
    losses = np.atleast_1d(np.asarray(losses_eV, dtype=float))
    res = [{t: np.zeros(losses.shape) for t in TERMS} for _ in materials]
    x = a                                              # theta = 0: x = a / cos(theta)
    for iw, w in enumerate(losses):
        kp = KPerp(E0, w, a)
        qs, ws = kpar_nodes(E0, w, materials)
        acc = [{t: 0.0 + 0.0j for t in TERMS + ("M1",)} for _ in materials]
        for q, wq in zip(qs, ws):
            I = kp.integrals(q)
            for im, mat in enumerate(materials):
                inv = inverse_epsilon_kpar(mat, w, q)
                eps = 1.0 / inv
                inv1 = inv - 1.0
                C1 = (eps - 1) * (eps + 2) / (eps * (eps + 1))
                C2 = (1 - eps) / (1 + eps)
                d = acc[im]
                f = wq * 2 * np.pi * q                 # d^3k = 2 pi k_par dk_par dk_perp
                d["M1"] += f * inv1 * I["I_ee"]
                d["term_14"] += f * inv1 * I["I_hh"]
                d["term_15_T1"] += f * I["I_S"] * C1 * I["F_ee"]
                d["term_15_T2"] += f * I["I_S"] * C1 * I["F_hh"]
                d["term_15_T3"] += f * I["I_S"] * C2 * I["F_V"]
                d["term_16_e"] += f * I["I_V"] * C2 * (I["F_ee"] + I["F_V"])
                d["term_16_h"] += f * I["I_V"] * C2 * I["F_hh"]
                d["term_13a"] += wq * q * inv1 / (kp.Q ** 2 + q ** 2)
        pref = 2.0 / ((2 * np.pi) ** 4 * x * w)       # 2/((2pi)^4 x hbar^2 w), hbar = 1
        # term_13a (bulk part).  S_M rho_e = (4 e^2/V^2) sin^2(u a/2)/u^2, u = k_perp - Q,
        # -> (2 pi e^2 a / V^2) delta(u) for a -> inf, so (k V) S_M rho_e -> (2 pi e^2 a Q/V) delta(u).
        # Inserted in Eq. (2) with x = a and Q/w = 1/V:
        # pref * 2pi * 4pi * 2pi e^2 a Q/V * Re[i Int q dq (1/eps-1)/(Q^2+q^2)]
        #   = (2 e^2/(pi V^2)) Re[...]   (= 1/(pi a0 E0) Re[...] since V^2 = 4 H E0, e^2 = 2H/a0).
        pre13a = 2.0 * E2 / (np.pi * kp.V ** 2)
        for im in range(len(materials)):
            d = acc[im]
            m1 = pref * np.real(1j * d["M1"])
            t13a = pre13a * np.real(1j * d["term_13a"])
            res[im]["term_13a"][iw] = t13a
            res[im]["term_13b"][iw] = m1 - t13a
            for t in ("term_14", "term_15_T1", "term_15_T2", "term_15_T3", "term_16_e", "term_16_h"):
                res[im][t][iw] = pref * np.real(1j * d[t])
    for r in res:
        r["no_hole"] = sum(r[t] for t in NO_HOLE_TERMS)
        r["hole"] = sum(r[t] for t in HOLE_TERMS)
    return res
