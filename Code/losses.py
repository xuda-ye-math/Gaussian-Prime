"""Loss integrals of the decompositions, two independent integrators, and theta_min.

Decompositions (each is a list of regions from regions.py and a rescue mode):
    "HL"          HL's Section 5 with the exact phi2 region D2 and the corrected phi10 (R* Section 4);
    "HL printed"  the same with HL's printed phi2 region;
    "R*"          HL with D2 replaced by D2 - R* and the new loss T3 (R* Section 5);
    "R*+G4", "R*+G4+P", "R*+G4+L"   R* with the four-dimensional losses reduced by the rescue sets of
                  R*2 Sections 3 and 4.
The lower-bound constant is 1 - (total loss); theta_min is the root of total(theta) = 1.

Primary integrator (`losses_gl`), deterministic:
  * 2-dim: the b-integral is exact, int_{b0}^{b1} omega((1-a-b)/b)/b^2 db = Omega((1-a)/b0 - 1) - Omega((1-a)/b1 - 1)
    divided by (1 - a); composite Gauss-Legendre in a, split at the region's kinks.
  * 4-dim: the d-integral is exact over [dlo, dhi] minus the union of the rescued d-intervals, through
    int_{d0}^{d1} omega((1-s-d)/d)/d^2 dd = (Omega((1-s)/d0 - 1) - Omega((1-s)/d1 - 1))/(1 - s), s = a + b + g;
    composite Gauss-Legendre in a, b, g.  The a-axis is split at the region's kinks, and every axis at the
    planes where the rescued d-set jumps (regions.JUMP_LEVELS), so each piece has a continuous integrand.
    The g-axis is also split where dlo or dhi meets the cap (1 - a - b - g)/2 (omega jumps from 0 to 1
    there) and where dlo meets dhi (the d-range opens): the g-integrand has a kink at these points.  dlo and
    dhi are affine in g between the jump planes, so the points are found exactly.
  Resolution (k, panels): each axis is cut into `panels` equal panels, the panels are cut further at the
  points above, and each resulting piece gets k Gauss-Legendre nodes.  The 2-dim losses use k nodes per piece
  and 16 * panels panels on the a-axis.  At the production setting k = 10, panels = 16 (run_all.py) this is
  256 panels on the a-axis for the 2-dim losses and 16 panels per axis for the 4-dim losses.
  With a process pool the a-nodes are distributed over the workers.

Cross-check integrator (`losses_check`), independent of the Omega substitution and of the d-interval
reduction:
  * 2-dim: nested adaptive quadrature (scipy quad) of omega((1-a-b)/b)/(a b^2) itself, with breakpoints at
    the kinks in a and at the b where (1-a-b)/b is an integer;
  * 4-dim: scrambled Sobol quasi-Monte Carlo on the unit cube mapped onto each region, omega evaluated
    directly, rescued tuples removed with the original definitions of R*2 (regions.rescued_original);
    standard error from independent scramblings.
"""
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.stats import qmc

import regions as R
from buchstab import Omega, omega

DECOMPOSITIONS = {
    "HL": ("D2", False, "none"),
    "HL printed": ("printed", False, "none"),
    "R*": ("D2-R*", True, "none"),
    "R*+G4": ("D2-R*", True, "G4"),
    "R*+G4+P": ("D2-R*", True, "G4+P"),
    "R*+G4+L": ("D2-R*", True, "G4+L"),
    "R*+G4+L[qrs|p]": ("D2-R*", True, "G4+L[qrs|p]"),
}
ORDER = ("phi1", "phi2", "D2-R*", "phi4", "phi6", "phi8", "phi3", "phi5", "phi7", "phi9", "phi10", "phi11", "T3")
mx = np.maximum


def _breakpoints(lo, hi, cuts, panels=1):
    """lo, the interior ends of `panels` equal panels, the cuts (those inside (lo, hi) for some entry,
    clipped to [lo, hi]) and hi, sorted along a new last axis (lo, hi arrays)."""
    lo = np.asarray(lo, dtype=float)
    hi = mx(np.asarray(hi, dtype=float), lo)
    even = [lo + (hi - lo) * j / panels for j in range(1, panels)]
    inside = [np.clip(c, lo, hi) for c in cuts if np.any((c > lo) & (c < hi))]
    return np.sort(np.stack(np.broadcast_arrays(lo, *even, *inside, hi), -1), -1)


def _composite(pts, k):
    """k Gauss-Legendre nodes and weights on each piece between consecutive breakpoints;
    shape pts.shape[:-1] + (m,)."""
    x, w = np.polynomial.legendre.leggauss(k)
    left, length = pts[..., :-1, None], np.diff(pts, axis=-1)[..., None]
    shape = pts.shape[:-1] + (-1,)
    return (left + length * (x + 1) / 2).reshape(shape), (length * w / 2).reshape(shape)


def _zeros(pts, h):
    """First and last zero of a function h that is affine on each piece between consecutive breakpoints
    (it may jump at a breakpoint, so it is sampled just inside each piece); +-inf where there is none."""
    x0, x1 = pts[..., :-1], pts[..., 1:]
    y0, y1 = x0 + 1e-9 * (x1 - x0), x1 - 1e-9 * (x1 - x0)
    h0, h1 = h(y0), h(y1)
    cross = h0 * h1 < 0
    z = y0 - h0 * (y1 - y0) / np.where(cross, h1 - h0, 1.0)
    return [np.min(np.where(cross, z, np.inf), -1), np.max(np.where(cross, z, -np.inf), -1)]


def _F(s, d):
    """Omega((1 - s)/d - 1): the d-antiderivative, up to the factor -1/(1 - s), of omega((1-s-d)/d)/d^2."""
    return Omega((1 - s) / d - 1)


def _union_integral(s, d0, d1, intervals):
    """Sum of F(lo) - F(hi) over the union of the intervals, clipped to [d0, d1] (sort, then running max)."""
    if not intervals:
        return np.zeros_like(s)
    lo = np.clip(np.stack([np.broadcast_to(l, s.shape) for l, _ in intervals], -1), d0[..., None], d1[..., None])
    hi = np.clip(np.stack([np.broadcast_to(h, s.shape) for _, h in intervals], -1), d0[..., None], d1[..., None])
    hi = mx(hi, lo)
    order = np.argsort(lo, axis=-1)
    lo, hi = np.take_along_axis(lo, order, -1), np.take_along_axis(hi, order, -1)
    before = np.concatenate([d0[..., None], np.maximum.accumulate(hi, axis=-1)[..., :-1]], -1)
    new_lo, new_hi = mx(lo, before), mx(hi, before)
    return np.sum(_F(s[..., None], new_lo) - _F(s[..., None], new_hi), axis=-1)


def loss2_gl(reg, k, panels):
    """2-dim loss of reg: exact b-integral through Omega, composite Gauss-Legendre in a."""
    a, wa = _composite(_breakpoints(reg.a0, reg.a1, reg.kinks, panels), k)
    b0, b1 = reg.blo(a), reg.bhi(a)
    ok = b1 > b0
    b0, b1 = np.where(ok, b0, 0.25), np.where(ok, b1, 0.25)
    val = (Omega((1 - a) / b0 - 1) - Omega((1 - a) / b1 - 1)) / (a * (1 - a))
    return float(np.sum(np.where(ok, wa * val, 0.0)))


def _alpha_nodes(reg, p, k, panels):
    """a-nodes and weights of a 4-dim region: split at its kinks and at the jump levels."""
    levels = [getattr(p, name) for name in R.JUMP_LEVELS]
    return _composite(_breakpoints(reg.a0, reg.a1, list(reg.kinks) + R.jump_points(levels, []), panels), k)


def _loss4_part(reg, mode, p, k, panels, A, WA, chunk=4):
    """Contribution of the a-nodes A (weights WA) to the 4-dim loss of reg."""
    levels = [getattr(p, name) for name in R.JUMP_LEVELS]
    total = 0.0
    for i in range(0, A.size, chunk):
        a, wa = A[i:i + chunk, None], WA[i:i + chunk, None]
        b, wb = _composite(_breakpoints(reg.blo(a)[:, 0], reg.bhi(a)[:, 0], R.jump_points(levels, [a[:, 0]]),
                                        panels), k)
        a2 = np.broadcast_to(a, b.shape)
        glo, ghi, jumps = reg.glo(a2, b), reg.ghi(a2, b), R.jump_points(levels, [a2, b])
        A3, B3 = a2[..., None], b[..., None]
        pts = _breakpoints(glo, ghi, jumps)
        kinks = [z for h in (lambda g: reg.dlo(A3, B3, g) - (1 - A3 - B3 - g) / 2,
                             lambda g: reg.dhi(A3, B3, g) - (1 - A3 - B3 - g) / 2,
                             lambda g: reg.dlo(A3, B3, g) - reg.dhi(A3, B3, g)) for z in _zeros(pts, h)]
        g, wg = _composite(_breakpoints(glo, ghi, jumps + kinks, panels), k)
        a3, b3 = np.broadcast_to(a2[..., None], g.shape), np.broadcast_to(b[..., None], g.shape)
        w = wa[..., None] * wb[..., None] * wg
        live = w > 0
        a3, b3, g = np.where(live, a3, 0.3), np.where(live, b3, 0.2), np.where(live, g, 0.1)
        s = a3 + b3 + g
        d0, d1 = reg.dlo(a3, b3, g), reg.dhi(a3, b3, g)
        live &= d1 > d0
        d0, d1 = np.where(live, d0, 0.05), np.where(live, d1, 0.05)
        inner = (_F(s, d0) - _F(s, d1)
                 - _union_integral(s, d0, d1, R.rescued_intervals(mode, a3, b3, g, p)))
        total += float(np.sum(np.where(live, w * inner / ((1 - s) * a3 * b3 * g), 0.0)))
    return total


def _loss4_job(job):
    """Worker entry: rebuild the region from theta (regions hold lambdas, which do not pickle)."""
    theta, rstar, cut, index, mode, k, panels, A, WA = job
    return _loss4_part(R.regions_4d(theta, rstar, cut)[index], mode, R.params(theta), k, panels, A, WA)


def losses_gl(theta, decomposition, k=10, panels=16, pool=None, cut=True):
    """Per-region losses (dict name -> value) with the primary integrator.  cut=False drops HL's cut
    g <= a + b - 2 theta in phi3 and phi9 (regions.regions_4d)."""
    phi2, rstar, mode = DECOMPOSITIONS[decomposition]
    p = R.params(theta)
    out = {}
    for reg in R.regions_2d(theta, phi2):
        out[reg.name] = out.get(reg.name, 0.0) + loss2_gl(reg, k, 16 * panels)
    jobs, names = [], []
    for index, reg in enumerate(R.regions_4d(theta, rstar, cut)):
        A, WA = _alpha_nodes(reg, p, k, panels)
        for part in np.array_split(np.arange(A.size), max(1, A.size // 4)):
            jobs.append((theta, rstar, cut, index, mode, k, panels, A[part], WA[part]))
            names.append(reg.name)
    for name, value in zip(names, (pool.map if pool else map)(_loss4_job, jobs)):
        out[name] = out.get(name, 0.0) + value
    return {name: out[name] for name in ORDER if name in out}


def total_gl(theta, decomposition, k=10, panels=16, pool=None):
    return sum(losses_gl(theta, decomposition, k, panels, pool).values())


# ----------------------------------------------------------------------------------------------------------
# Cross-check integrator.

def loss2_quad(reg):
    """Nested adaptive quadrature of omega((1-a-b)/b)/(a b^2) over the region (omega evaluated directly)."""
    def inner(a):
        b0, b1 = float(reg.blo(a)), float(reg.bhi(a))
        if b1 <= b0:
            return 0.0
        pts = [(1 - a) / (j + 1) for j in range(1, 8) if b0 < (1 - a) / (j + 1) < b1]
        f = lambda b: float(omega((1 - a - b) / b)) / (a * b * b)
        return quad(f, b0, b1, points=pts or None, epsabs=1e-13, epsrel=1e-12, limit=200)[0]
    pts = [x for x in reg.kinks if reg.a0 < x < reg.a1]
    return quad(inner, reg.a0, reg.a1, points=pts or None, epsabs=1e-12, epsrel=1e-11, limit=200)[0]


def _qmc_estimate(reg, mode, p, m, seed, chunk=2 ** 16):
    """One scrambled Sobol QMC estimate of the 4-dim loss of reg, with 2^m points."""
    u = qmc.Sobol(4, scramble=True, seed=seed).random_base2(m)
    acc = 0.0
    for i in range(0, len(u), chunk):
        uu = u[i:i + chunk]
        a = reg.a0 + (reg.a1 - reg.a0) * uu[:, 0]
        b0, b1 = reg.blo(a), reg.bhi(a)
        lb = mx(b1 - b0, 0.0)
        b = b0 + lb * uu[:, 1]
        g0, g1 = reg.glo(a, b), reg.ghi(a, b)
        lg = mx(g1 - g0, 0.0)
        g = g0 + lg * uu[:, 2]
        s = a + b + g
        d0 = reg.dlo(a, b, g)
        ld = mx(np.minimum(reg.dhi(a, b, g), (1 - s) / 2) - d0, 0.0)
        d = d0 + ld * uu[:, 3]
        jac = (reg.a1 - reg.a0) * lb * lg * ld
        live = jac > 0
        a, b, g, d = (np.where(live, x, y) for x, y in ((a, 0.3), (b, 0.2), (g, 0.1), (d, 0.05)))
        f = omega((1 - a - b - g - d) / d) / (a * b * g * d * d)
        f = np.where(live & ~R.rescued_original(mode, [a, b, g, d], p), jac * f, 0.0)
        acc += float(np.sum(f))
    return acc / len(u)


def _qmc_job(job):
    theta, rstar, index, mode, m, seed = job
    return _qmc_estimate(R.regions_4d(theta, rstar)[index], mode, R.params(theta), m, seed)


def losses_check(theta, decomposition, m=20, reps=8, seed=1, pool=None):
    """Per-region (value, standard error) with the cross-check integrator: nested adaptive quadrature for the
    2-dim losses (s.e. reported as 0), and the mean of `reps` scrambled Sobol estimates with 2^m points each
    for the 4-dim losses (s.e. = sample standard deviation / sqrt(reps))."""
    phi2, rstar, mode = DECOMPOSITIONS[decomposition]
    out = {}
    for reg in R.regions_2d(theta, phi2):
        out[reg.name] = (out.get(reg.name, (0.0, 0.0))[0] + loss2_quad(reg), 0.0)
    regs = R.regions_4d(theta, rstar)
    jobs = [(theta, rstar, index, mode, m, seed + r) for index in range(len(regs)) for r in range(reps)]
    est = np.array(list((pool.map if pool else map)(_qmc_job, jobs))).reshape(len(regs), reps)
    for reg, e in zip(regs, est):
        out[reg.name] = (float(e.mean()), float(e.std(ddof=1) / np.sqrt(reps)))
    return {name: out[name] for name in ORDER if name in out}


CLASSIFICATION_TESTS = ("G4", "I", "L", "G4+P", "G4+L", "intervals")


def _classification_job(job):
    """Mismatch counts {test: (comparisons, mismatches)} at n random points of one 4-dim region."""
    theta, index, n, seed, chunk = job
    p = R.params(theta)
    reg = R.regions_4d(theta, True)[index]
    rng = np.random.default_rng(seed)
    count = {test: 0 for test in CLASSIFICATION_TESTS}
    points = 0
    for _ in range(0, n, chunk):
        u = rng.random((chunk, 4))
        a = reg.a0 + (reg.a1 - reg.a0) * u[:, 0]
        b = reg.blo(a) + (reg.bhi(a) - reg.blo(a)) * u[:, 1]
        g = reg.glo(a, b) + (reg.ghi(a, b) - reg.glo(a, b)) * u[:, 2]
        d0, d1 = reg.dlo(a, b, g), np.minimum(reg.dhi(a, b, g), (1 - a - b - g) / 2)
        d = d0 + (d1 - d0) * u[:, 3]
        keep = (reg.bhi(a) > reg.blo(a)) & (reg.ghi(a, b) > reg.glo(a, b)) & (d1 > d0)
        v = [x[keep] for x in (a, b, g, d)]
        points += v[0].size
        typeI, cover = R.typeI_original(v, p), R.cover_original(v, p)
        count["G4"] += np.sum(R.g4_original(v, p) != R.g4_reduced(v, p))
        count["I"] += np.sum(typeI != R.typeI_reduced(v, p))
        count["L"] += np.sum((typeI & cover) != R.l_reduced(v, p))
        for mode in ("G4+P", "G4+L"):
            count[mode] += np.sum(R.rescued_original(mode, v, p) != R.rescued_reduced(mode, v, p))
        for mode in R.MODES:
            inside = np.zeros(v[0].shape, dtype=bool)
            for lo, hi in R.rescued_intervals(mode, *v[:3], p):
                inside |= (v[3] >= lo) & (v[3] <= hi)
            count["intervals"] += np.sum(inside != R.rescued_reduced(mode, v, p))
    return {test: (points * (len(R.MODES) if test == "intervals" else 1), int(c)) for test, c in count.items()}


def classification_mismatches(theta, n=10 ** 6, seed=2, chunk=2 ** 16, pool=None):
    """Brute-force comparison at n random points of each 4-dim region (inside the cap d <= (1-a-b-g)/2):
    "G4", "I", "L": original definition of R*2 against the reduced definition;
    "G4+P", "G4+L": the rescued sets of these modes, original against reduced;
    "intervals": the reduced definition against its d-interval form, in every mode of regions.MODES.
    Returns {region: {test: (comparisons, mismatches)}}."""
    regs = R.regions_4d(theta, True)
    jobs = [(theta, index, n, seed + index, chunk) for index in range(len(regs))]
    return dict(zip((reg.name for reg in regs), (pool.map if pool else map)(_classification_job, jobs)))


def theta_min(decomposition, bracket, k=10, panels=16, pool=None, log=print):
    """Root of total(theta) = 1 (brentq, xtol 1e-10), logging every evaluation."""
    def f(theta):
        t0 = time.time()
        value = total_gl(theta, decomposition, k, panels, pool)
        log(f"    {decomposition:10s} theta = {theta:.10f}  total = {value:.9f}  ({time.time() - t0:.1f}s)", flush=True)
        return value - 1.0
    return brentq(f, *bracket, xtol=1e-10)


if __name__ == "__main__":
    theta, decomposition = float(sys.argv[1]), sys.argv[2]
    with ProcessPoolExecutor() as pool:
        out = losses_gl(theta, decomposition, pool=pool)
    for name, value in out.items():
        print(f"{name:6s} {value:.7f}", flush=True)
    print(f"TOTAL  {sum(out.values()):.7f}   1 - TOTAL {1 - sum(out.values()):.7f}", flush=True)
