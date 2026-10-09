"""Parameters, loss regions and rescue sets of the Harman--Lewis sieve and its refinements.

Sources.  "HL" is G. Harman and P. Lewis, Gaussian primes in narrow sectors, Mathematika 48 (2001) 119-135
(journal page numbers).  "R*" is routes/R1/rstar_claim.md, "R*2" is routes/R1/rstar2_claim.md; section
numbers below refer to these files.  All exponents are exponents of norms: |p|^2 = X^a, |q|^2 = X^b,
|r|^2 = X^g, |s|^2 = X^d, with a > b > g > d.

A two-dimensional region is {a0 <= a <= a1, blo(a) <= b <= bhi(a)}; its loss is the integral of
omega((1 - a - b)/b)/(a b^2).  A four-dimensional region adds glo(a, b) <= g <= ghi(a, b) and
dlo(a, b, g) <= d <= dhi(a, b, g); its loss is the integral of omega((1 - a - b - g - d)/d)/(a b g d^2), where
omega(u) = 0 for u < 1 caps d at (1 - a - b - g)/2.  `kinks` lists the values of a at which a limit of the
region has a kink or the region opens or closes; the primary integrator splits there.

A rescue set is a set of 4-tuples whose discarded sum gets an asymptotic formula instead (R*2 Sections 3, 4).
It is given in three forms: pointwise with the original definitions of R*2 (used by the QMC cross-check),
pointwise with the reduced definitions stated in the paper, and as a union of d-intervals of the reduced
definitions for fixed (a, b, g) (used by the primary integrator).  The three are compared pointwise in
losses.classification_mismatches.

Everything here accepts floats or numpy arrays; `params` also accepts fractions.Fraction.
"""
import itertools
from dataclasses import dataclass
from typing import Callable

import numpy as np

mn, mx = np.minimum, np.maximum


@dataclass(frozen=True)
class Params:
    """HL's parameters (HL p. 120; R* Section 1): sigma = 1/2 - theta, tau = 1/2 + theta,
    rho = (1 + 22 theta)/12, zeta = (1 + 4 theta)/3, nu = (10 theta + 7)/12, chi = (34 theta + 7)/24.
    HL print chi(0.119) = 0.46205; the formula gives 0.46025 and is used (R* Section 1)."""
    theta: object
    sigma: object
    tau: object
    rho: object
    zeta: object
    nu: object
    chi: object


def params(theta):
    return Params(theta=theta, sigma=(1 - 2 * theta) / 2, tau=(1 + 2 * theta) / 2, rho=(1 + 22 * theta) / 12,
                  zeta=(1 + 4 * theta) / 3, nu=(10 * theta + 7) / 12, chi=(34 * theta + 7) / 24)


@dataclass(frozen=True)
class Region2:
    name: str
    a0: float
    a1: float
    blo: Callable
    bhi: Callable
    kinks: tuple = ()


@dataclass(frozen=True)
class Region4:
    name: str
    a0: float
    a1: float
    blo: Callable
    bhi: Callable
    glo: Callable
    ghi: Callable
    dlo: Callable
    dhi: Callable
    kinks: tuple = ()


PHI2_READINGS = ("printed", "D2", "D2-R*")


def regions_2d(theta, phi2):
    """The two-dimensional losses phi1, phi2, phi4, phi6, phi8.

    phi2 = "printed": HL's printed region {1/3 <= a < sigma, zeta - a <= b < (1 - a)/2} (R* Section 4);
    phi2 = "D2": the exact discarded set D2, lower limit max(zeta - a, (tau - a)/2) (R* Section 4);
    phi2 = "D2-R*": D2 minus the regrouped set R* = {b <= rho/2} (R* Section 5.1), reported as "D2-R*"."""
    p = params(theta)
    s, t, r, z, n, c = p.sigma, p.tau, p.rho, p.zeta, p.nu, p.chi
    d2lo = lambda a: mx(z - a, (t - a) / 2)
    phi2 = {
        "printed": Region2("phi2", 1 / 3, s, lambda a: z - a, lambda a: (1 - a) / 2),
        "D2": Region2("phi2", 1 / 3, s, d2lo, lambda a: (1 - a) / 2, kinks=(2 * z - t,)),
        "D2-R*": Region2("D2-R*", 1 / 3, s, lambda a: mx(d2lo(a), r / 2), lambda a: (1 - a) / 2,
                         kinks=(2 * z - t, z - r / 2, t - r)),
    }[phi2]
    return [
        # phi1, HL Section 1 (p. 132), three pieces as printed; R* Section 4 (last paragraph) and
        # R1/notes.md (Phase 1 transcription).
        Region2("phi1", s, c, lambda a: n - a, lambda a: (1 - a) / 2),
        Region2("phi1", 0.25 + 1.5 * theta, c, lambda a: t - a, lambda a: a - 2 * theta),
        Region2("phi1", c, 0.5, lambda a: t - a, lambda a: (1 - a) / 2),
        phi2,
        # phi4, HL Section 3 (p. 132); R* Section 4 (last paragraph), Section 5.6.
        Region2("phi4", r, 1 / 3, lambda a: z - a, lambda a: a),
        # phi6, HL Section 4 (p. 132); R* Section 4 (last paragraph), Section 5.6.
        Region2("phi6", 0.5 - 2 * theta, r, lambda a: (n - a) / 2, lambda a: a),
        # phi8, HL Section 5 (p. 133); R* Section 4 (last paragraph), Section 5.6.
        Region2("phi8", n / 3, s - theta, lambda a: (n - a) / 2, lambda a: a),
    ]


def regions_4d(theta, rstar, cut=True):
    """The four-dimensional losses phi3, phi5, phi7, phi9, phi10, phi11 (R*2 Section 2, table), and T3 if
    rstar.  phi10 and phi11 are integrated over their full a-ranges for every theta, which includes the strip
    [tau/4, sigma - 2 theta) below theta = 3/26 (R* Section 6).  cut=False drops HL's cut
    g <= a + b - 2 theta in phi3 and phi9 (see the comment below)."""
    p = params(theta)
    s, t, r, z, n = p.sigma, p.tau, p.rho, p.zeta, p.nu
    # HL's cut g <= a + b - 2 theta in phi3 and phi9 (Lemma 11 applied to pq and r) is implied by the
    # d-limits tau - a - b <= d <= (1 - a - b - g)/2: these give tau - a - b <= (1 - a - b - g)/2, i.e.
    # g <= 1 - 2 tau + a + b = a + b - 2 theta.  Where the cut binds, the d-range is empty, so the cut changes
    # no loss (checked numerically by run_all.py, step "cut").
    regs = [
        # phi3: Lemma 13, (35), HL Section 2 (p. 132); R*2 Section 2.  The b-range closes at a = 2 zeta - tau.
        Region4("phi3", 1 / 3, s, lambda a: (t - a) / 2, lambda a: z - a,
                lambda a, b: t - a - b, (lambda a, b: mn(b, a + b - 2 * theta)) if cut else (lambda a, b: b),
                lambda a, b, g: t - a - b, lambda a, b, g: g, kinks=(2 * z - t,)),
        # phi5: Lemma 13, (35), HL Section 3 (p. 132); R*2 Section 2.
        Region4("phi5", r, 1 / 3, lambda a: (t - a) / 2, lambda a: z - a,
                lambda a, b: t - a - b, lambda a, b: b,
                lambda a, b, g: t - a - b, lambda a, b, g: g),
        # phi7: Lemma 13, (35), HL Section 4 (p. 133); R*2 Section 2.  HL's exclusion "d > 1 - a - g - sigma
        # if a + g >= sigma" (R*2 Section 1, p. 133) is kept in the region; it lies inside G4 (x = s,
        # S = {p, r}), so it changes nothing when G4 is used, and it reproduces R* alone otherwise.
        Region4("phi7", 0.5 - 2 * theta, r, lambda a: (t - a) / 2, lambda a: (n - a) / 2,
                lambda a, b: t - a - b, lambda a, b: b,
                lambda a, b, g: np.where(a + g >= s, mx(t - a - b, 1 - a - g - s), t - a - b),
                lambda a, b, g: g),
        # phi9: Lemma 13, (35), HL Section 5 with pq >= W (p. 133); R*2 Section 2.
        Region4("phi9", t / 3, s - theta, lambda a: (t - a) / 2, lambda a: mn((n - a) / 2, a),
                lambda a, b: t - a - b, (lambda a, b: mn(b, b + a - 2 * theta)) if cut else (lambda a, b: b),
                lambda a, b, g: t - a - b, lambda a, b, g: g, kinks=(n / 3, 2 * theta)),
        # phi10: Lemma 14, (36), HL Section 5 with pq < W (p. 133), with the corrected upper b-limit
        # min(a, sigma - a) (R* Section 4); R*2 Section 2.
        Region4("phi10", t / 4, s - theta, lambda a: (t - a) / 3, lambda a: mn(a, s - a),
                lambda a, b: (t - a - b) / 2, lambda a, b: b,
                lambda a, b, g: t - a - b - g, lambda a, b, g: g, kinks=(s / 2,)),
        # phi11: HL Section 6 (p. 134), the j = 1 analogue of (36); R* Section 6 item 3; R*2 Section 2.
        Region4("phi11", t / 5, s - 2 * theta, lambda a: (t - a) / 4, lambda a: mn(a, s - theta - a),
                lambda a, b: (t - a - b) / 3, lambda a, b: b,
                lambda a, b, g: (t - a - b - g) / 2, lambda a, b, g: g, kinks=((s - theta) / 2,)),
    ]
    if rstar:
        # T3: third term of (35) on R* = {1/3 <= a < sigma, max(zeta - a, (tau - a)/2) <= b <= rho/2}
        # (R* Sections 5.1, 5.2; R*2 Section 2).
        regs.append(Region4("T3", 1 / 3, s, lambda a: mx(z - a, (t - a) / 2), lambda a: r / 2 + 0 * a,
                            lambda a, b: t - a - b, lambda a, b: b,
                            lambda a, b, g: t - a - b, lambda a, b, g: g,
                            kinks=(2 * z - t, z - r / 2, t - r)))
    return regs


# ----------------------------------------------------------------------------------------------------------
# Rescue sets of R*2.  v = [a, b, g, d] (arrays of one shape).  Modes:
#   "none"            no rescue;
#   "G4"              G4;          "G4+P"   G4 union P;          "G4+L"   G4 union L;
#   "G4+L[qrs|p]"     G4 union L, with the Type I condition I restricted to Lemma 6 and the single
#                     Lemma 7 split U = qrs, V = p (reduced form only; it has no counterpart in the claim files).

MODES = ("none", "G4", "G4+P", "G4+L", "G4+L[qrs|p]")
_I4 = (0, 1, 2, 3)
ALL_SPLITS = None                    # every split U | V of {p, q, r, s} (U and V non-empty)
SPLIT_QRS_P = ((1, 2, 3),)           # U = {q, r, s}, V = {p}
_MODE_SPLITS = {"G4+L": ALL_SPLITS, "G4+L[qrs|p]": SPLIT_QRS_P}


def _subsets(items):
    return [S for k in range(1, len(items) + 1) for S in itertools.combinations(items, k)]


# Original definitions, as stated in the claim file R*2.

def g4_original(v, p):
    """G4 of R*2 Section 3: some prime x of p, q, r, s and a non-empty set S of the other three with
    sigma <= sum_S <= tau - x (28 conditions)."""
    ok = np.zeros(np.shape(v[0]), dtype=bool)
    for x in _I4:
        for S in _subsets([i for i in _I4 if i != x]):
            sS = sum(v[i] for i in S)
            ok |= (sS >= p.sigma) & (sS <= p.tau - v[x])
    return ok


def typeI_original(v, p):
    """(I) of R*2 Section 4: Lemma 6 (a + b + g + d <= tau), or, for some split A | B of {p, q, r, s},
    max(u, sigma) + max(v, sigma/2) <= nu with u = sum_A, v = sum_B (the form of R* Section 3)."""
    s4 = v[0] + v[1] + v[2] + v[3]
    ok = s4 <= p.tau
    for A in _subsets(_I4)[:-1]:
        u = sum(v[i] for i in A)
        ok |= mx(u, p.sigma) + mx(s4 - u, p.sigma / 2) <= p.nu
    return ok


def cover_original(v, p):
    """(C) of R*2 Section 4: every e in [0, d] lies in one of the 75 e-intervals given by a prime x among
    p, q, r, s, v (v of exponent e) and a non-empty set S of the other four with sigma <= sum_S <= tau - x.
    The cover of [0, d] is decided by growing the reach from 0 until it is stable."""
    v5 = list(v) + [None]
    lo, hi = [], []
    shape = np.shape(v[0])
    for x in range(5):
        for S in _subsets([i for i in range(5) if i != x]):
            if x == 4:                                    # e <= tau - sum_S, if sum_S >= sigma
                sS = sum(v5[i] for i in S)
                lo.append(np.full(shape, -np.inf))
                hi.append(np.where(sS >= p.sigma, p.tau - sS, -np.inf))
            elif 4 in S:                                  # sigma - S' <= e <= tau - x - S'
                sp = sum((v5[i] for i in S if i != 4), np.zeros(shape))
                lo.append(p.sigma - sp)
                hi.append(p.tau - v5[x] - sp)
            else:                                         # independent of e
                sS = sum(v5[i] for i in S)
                ok = (sS >= p.sigma) & (sS <= p.tau - v5[x])
                lo.append(np.where(ok, -np.inf, np.inf))
                hi.append(np.where(ok, np.inf, -np.inf))
    lo, hi = np.stack(lo, -1), np.stack(hi, -1)
    reach = np.zeros(shape)
    while True:
        new = mx(reach, np.max(np.where(lo <= reach[..., None], hi, -np.inf), axis=-1))
        if np.array_equal(new, reach):
            return reach >= v[3]
        reach = new


def p_set(v, p):
    """P of R*2 Section 4 (the same in both forms): a + b <= sigma, g + d <= rho, a + b + d >= sigma,
    a + b + 2d <= tau."""
    a, b, g, d = v
    return (a + b <= p.sigma) & (g + d <= p.rho) & (a + b + d >= p.sigma) & (a + b + 2 * d <= p.tau)


def rescued_original(mode, v, p):
    """Indicator of the rescued 4-tuples with the definitions of R*2 (modes none, G4, G4+P, G4+L)."""
    if mode not in ("none", "G4", "G4+P", "G4+L"):
        raise ValueError(f"mode {mode} has no original definition in R*2")
    if mode == "none":
        return np.zeros(np.shape(v[0]), dtype=bool)
    ok = g4_original(v, p)
    if mode == "G4+P":
        ok |= p_set(v, p)
    elif mode == "G4+L":
        ok |= typeI_original(v, p) & cover_original(v, p)
    return ok


# Reduced definitions (the form stated in the paper).
#   G4 = union over non-empty S of {p, q, r} of {sigma <= sum_S <= tau - d}: (30) with the smallest prime s
#        as the prime of (30).  Equal to the original G4 by the swap lemma (R*2 Section 3.1).
#   I  = {a + b + g + d < tau} (Lemma 6) union, over the splits U | V of {p, q, r, s}, of
#        {sum_U <= zeta, sum_V <= rho, sum_U + sum_V <= nu} (Lemma 7 in HL's form, HL p. 127).  Equal to the
#        original (I) for theta >= 1/14 (R* Section 3).
#   L  = I intersected with the union over non-empty S of {p, q, r, s} of {sigma <= sum_S <= tau - d}.  This
#        is (C): by the swap lemma the prime of (30) may be the fifth prime v, and the condition on its
#        exponent e is monotone in e, so it is needed at e = d only (R*2 Section 4, "Computation of (C)").
#   P  as above.

def g4_reduced(v, p):
    """G4 in the reduced form (the prime of (30) is s)."""
    ok = np.zeros(np.shape(v[0]), dtype=bool)
    for S in _subsets((0, 1, 2)):
        sS = sum(v[i] for i in S)
        ok |= (sS >= p.sigma) & (sS <= p.tau - v[3])
    return ok


def typeI_reduced(v, p, splits=ALL_SPLITS):
    """I in the reduced form; `splits` lists the U-sets of the allowed splits (ALL_SPLITS: all 14)."""
    s4 = v[0] + v[1] + v[2] + v[3]
    ok = s4 < p.tau
    for U in (_subsets(_I4)[:-1] if splits is None else splits):
        u = sum(v[i] for i in U)
        ok |= (u <= p.zeta) & (s4 - u <= p.rho) & (s4 <= p.nu)
    return ok


def l_reduced(v, p, splits=ALL_SPLITS):
    """L in the reduced form: I (with the given splits) and the cover condition at e = d."""
    cover = np.zeros(np.shape(v[0]), dtype=bool)
    for S in _subsets(_I4):
        sS = sum(v[i] for i in S)
        cover |= (sS >= p.sigma) & (sS <= p.tau - v[3])
    return typeI_reduced(v, p, splits) & cover


def rescued_reduced(mode, v, p):
    """Indicator of the rescued 4-tuples with the reduced definitions (every mode of MODES)."""
    if mode not in MODES:
        raise ValueError(f"unknown mode {mode}")
    if mode == "none":
        return np.zeros(np.shape(v[0]), dtype=bool)
    ok = g4_reduced(v, p)
    if mode == "G4+P":
        ok |= p_set(v, p)
    elif mode in _MODE_SPLITS:
        ok |= l_reduced(v, p, _MODE_SPLITS[mode])
    return ok


# The reduced sets as unions of d-intervals for fixed (a, b, g), for the primary integrator.
#   G4: d in [0, tau - sum_S] for each non-empty S of {p, q, r} with sum_S >= sigma.
#   I:  a down-set [0, dI]: Lemma 6 gives d < tau - (a + b + g); a split with s on the U-side gives, if
#       sum_V <= rho, d <= min(zeta - u', nu - (a + b + g)); with s on the V-side, if sum_U <= zeta,
#       d <= min(rho - v', nu - (a + b + g)); u', v' are the sums of the other primes on the side of s.
#   L:  I intersected with the G4 intervals (S without s) and with [sigma - S', (tau - S')/2] for S = S' + {s},
#       S' a subset of {p, q, r} (possibly empty).  So G4 union L = G4 union (I intersected with the latter).
#   P:  d in [sigma - a - b, min(rho - g, (tau - a - b)/2)] if a + b <= sigma.
# The d-sets jump where a subset sum of (a, b, g) crosses sigma (G4, P, phi7's exclusion), rho or zeta
# (the two sides of Lemma 7); these are JUMP_LEVELS, and the primary integrator splits every axis there.

JUMP_LEVELS = ("sigma", "rho", "zeta")


def _sums3(a, b, g):
    return {S: sum((a, b, g)[i] for i in S) for S in [()] + _subsets((0, 1, 2))}


def typeI_threshold(a, b, g, p, splits=ALL_SPLITS):
    """dI with I <=> d <= dI (up to the boundary), for fixed (a, b, g)."""
    s3 = a + b + g
    dI = p.tau - s3                                                        # Lemma 6
    for U in (_subsets(_I4)[:-1] if splits is None else splits):
        u = sum((a, b, g)[i] for i in U if i != 3)                        # U without s
        v = s3 - u                                                          # V without s
        if 3 in U:
            dI = mx(dI, np.where(v <= p.rho, mn(p.zeta - u, p.nu - s3), -np.inf))
        else:
            dI = mx(dI, np.where(u <= p.zeta, mn(p.rho - v, p.nu - s3), -np.inf))
    return dI


def rescued_intervals(mode, a, b, g, p):
    """List of (lo, hi) arrays: the d-intervals rescued in mode (reduced definitions), for fixed (a, b, g)."""
    if mode == "none":
        return []
    sums = _sums3(a, b, g)
    out = [(np.zeros_like(a), np.where(sS >= p.sigma, p.tau - sS, -np.inf))
           for S, sS in sums.items() if S]                                         # G4
    if mode == "G4+P":
        ok = a + b <= p.sigma
        out.append((np.where(ok, p.sigma - a - b, np.inf), mn(p.rho - g, (p.tau - a - b) / 2)))
    elif mode in _MODE_SPLITS:
        dI = typeI_threshold(a, b, g, p, _MODE_SPLITS[mode])
        out += [(p.sigma - sS + 0 * a, mn((p.tau - sS) / 2, dI)) for sS in sums.values()]
    return out


def jump_points(level_values, fixed):
    """Values x of the next integration variable with x + sum_S' = c, for every level c and every subset S'
    of the already fixed variables: the jump planes of the rescued d-sets met along this axis."""
    out = []
    for c in level_values:
        for k in range(len(fixed) + 1):
            for S in itertools.combinations(fixed, k):
                out.append(c - sum(S, 0 * fixed[0] if fixed else 0.0))
    return out
