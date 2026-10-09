"""Buchstab's function omega(u) and its integral Omega(u) = int_1^u omega(t) dt.

omega(u) = 0 for u < 1, omega(u) = 1/u on [1, 2], and u omega(u) = 1 + Omega(u - 1) for u >= 2, i.e.
(u omega(u))' = omega(u - 1).  This is the omega of HL's (17), p. 123 (rstar_claim.md Section 2), which
turns every discarded sum into an integral of omega((1 - a - b)/b)/(a b^2) or of its four-dimensional analogue.

Construction.  On each unit interval [k, k + 1] the functions are Chebyshev series of degree DEG, built in
order of k: omega_k(u) = (1 + Omega_{k-1}(u - 1))/u and Omega_k(u) = Omega(k) + int_k^u omega_k.  On a closed
unit interval omega is analytic (its singularities lie at the integers to the left), so the series converge
geometrically.  For vectorized evaluation, the values and one-sided derivatives are tabulated from the series
on the grid k + j/NT of each unit interval, and evaluated by cubic Hermite interpolation: Omega with
derivative omega, and omega with derivative omega'(u) = (omega(u - 1) - omega(u))/u (from the delay
equation).  For u >= UMAX, omega(u) = e^{-gamma} to double precision and Omega is continued linearly.

Closed forms used by `check`: omega(u) = 1/u and Omega(u) = log u on [1, 2]; on [2, 3],
omega(u) = (1 + log(u - 1))/u and Omega(u) = log u + log(u - 1) log u + Li2(1 - u) + pi^2/12, where
Li2(1 - u) = scipy.special.spence(u).  The limit is omega(u) -> e^{-gamma} as u -> infinity.
"""
import numpy as np
from numpy.polynomial import Chebyshev
from scipy.special import spence

DEG = 40
UMAX = 40
NT = 2000
EXP_MINUS_GAMMA = np.exp(-np.euler_gamma)


def _series():
    """Chebyshev series (omega_k, Omega_k) on [k, k + 1] for k = 1, ..., UMAX - 1."""
    om = {1: Chebyshev.interpolate(lambda u: 1.0 / u, DEG, domain=[1, 2])}
    Om = {1: om[1].integ(lbnd=1)}
    for k in range(2, UMAX):
        prev = Om[k - 1]
        om[k] = Chebyshev.interpolate(lambda u, prev=prev: (1.0 + prev(u - 1.0)) / u, DEG, domain=[k, k + 1])
        Om[k] = om[k].integ(lbnd=k) + Om[k - 1](float(k))
    return om, Om


OMEGA_SERIES, BIG_OMEGA_SERIES = _series()


def _tables():
    """Per unit interval k (row k - 1): values of omega, omega', Omega at k + j/NT, j = 0, ..., NT."""
    t = np.arange(NT + 1) / NT
    w, dw, W = (np.empty((UMAX - 1, NT + 1)) for _ in range(3))
    for k in range(1, UMAX):
        u = k + t
        w[k - 1] = OMEGA_SERIES[k](u)
        dw[k - 1] = OMEGA_SERIES[k].deriv()(u)
        W[k - 1] = BIG_OMEGA_SERIES[k](u)
    return w, dw, W


_W, _DW, _OM = _tables()
_OMEGA_END = float(BIG_OMEGA_SERIES[UMAX - 1](float(UMAX)))


def _hermite(u, val, der):
    """Cubic Hermite interpolation of a tabulated function on [1, UMAX) (per-interval one-sided derivatives)."""
    k = np.floor(u).astype(np.int64)
    x = (u - k) * NT
    j = np.minimum(x.astype(np.int64), NT - 1)
    s = x - j
    r = k - 1
    y0, y1 = val[r, j], val[r, j + 1]
    d0, d1 = der[r, j] / NT, der[r, j + 1] / NT
    s2 = s * s
    s3 = s2 * s
    return ((2 * s3 - 3 * s2 + 1) * y0 + (s3 - 2 * s2 + s) * d0
            + (3 * s2 - 2 * s3) * y1 + (s3 - s2) * d1)


def omega(u):
    """Buchstab's omega(u), vectorized; 0 for u < 1, e^{-gamma} for u >= UMAX."""
    u = np.asarray(u, dtype=float)
    uc = np.clip(u, 1.0, UMAX - 1e-9)
    val = _hermite(uc, _W, _DW)
    return np.where(u < 1.0, 0.0, np.where(u >= UMAX, EXP_MINUS_GAMMA, val))


def Omega(u):
    """Omega(u) = int_1^u omega(t) dt, vectorized; 0 for u <= 1, linear with slope e^{-gamma} beyond UMAX."""
    u = np.asarray(u, dtype=float)
    uc = np.clip(u, 1.0, UMAX - 1e-9)
    val = _hermite(uc, _OM, _W)
    return np.where(u <= 1.0, 0.0, np.where(u >= UMAX, _OMEGA_END + EXP_MINUS_GAMMA * (u - UMAX), val))


def check():
    """Accuracy of omega and Omega against the closed forms on [1, 3], the limit e^{-gamma}, and the
    delay equation u omega(u) = 1 + Omega(u - 1) on [2, UMAX)."""
    u12 = np.linspace(1.0, 2.0, 20001)
    u23 = np.linspace(2.0, 3.0, 20001)
    Om23 = np.log(u23) + np.log(u23 - 1.0) * np.log(u23) + spence(u23) + np.pi ** 2 / 12
    rng = np.random.default_rng(0)
    ur = rng.uniform(2.0, UMAX - 1e-6, 200000)
    return {
        "max |omega - 1/u| on [1,2]": float(np.max(np.abs(omega(u12) - 1.0 / u12))),
        "max |Omega - log u| on [1,2]": float(np.max(np.abs(Omega(u12) - np.log(u12)))),
        "max |omega - (1+log(u-1))/u| on [2,3]": float(np.max(np.abs(omega(u23) - (1.0 + np.log(u23 - 1.0)) / u23))),
        "max |Omega - closed form| on [2,3]": float(np.max(np.abs(Omega(u23) - Om23))),
        "max |u omega(u) - 1 - Omega(u-1)| on [2,40)": float(np.max(np.abs(ur * omega(ur) - 1.0 - Omega(ur - 1.0)))),
        "omega(10) - e^-gamma": float(omega(10.0) - EXP_MINUS_GAMMA),
        "omega(39.9) - e^-gamma": float(omega(39.9) - EXP_MINUS_GAMMA),
        "omega(2.5)": float(omega(2.5)),
        "omega(3)": float(omega(3.0)),
    }


if __name__ == "__main__":
    for name, value in check().items():
        print(f"{name:46s} {value: .3e}" if abs(value) < 1e-3 else f"{name:46s} {value:.15f}", flush=True)
    print(f"{'e^-gamma':46s} {EXP_MINUS_GAMMA:.15f}", flush=True)
