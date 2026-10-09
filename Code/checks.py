"""Structural inequalities and lemma margins, in exact rational arithmetic.

Each check is a margin m(theta), required to be > 0, >= 0 or = 0, or reported for information.  Every margin
is evaluated exactly (fractions.Fraction) at theta = 0.1129, 0.1133, 0.1142, 0.119, and written as an exact
linear form c1 theta + c0, obtained from the values at 0.1129 and 0.119 and confirmed exactly at 0.1133 and
0.1142.  Margins that are infima over a region (rows "(U1)", "(U2)", "Q >= X^eps", "Q <= 2P") are computed
exactly as linear programs over the region's polygon in (a, b), by enumerating its vertices.

Sources: "gaps 2.6", "gaps 3.3" are the tables of Sections 2.6 and 3.3 of routes/R1-gaps/gaps_v3.md, whose
values (7 decimals where printed) are compared with the exact values; "R* 3", "R* 6", "R* 7" are sections of
routes/R1/rstar_claim.md.
"""
import itertools
from dataclasses import dataclass, field
from fractions import Fraction as F
from typing import Callable

from regions import params

THETAS = (F("0.1129"), F("0.1133"), F("0.1142"), F("0.119"))


@dataclass(frozen=True)
class Check:
    source: str
    name: str
    kind: str                    # ">", ">=", "=" or "info"
    margin: Callable             # Params -> Fraction
    gaps: dict = field(default_factory=dict)   # {theta string: printed value string} from gaps_v3.md


def polygon_min(halfplanes, objective):
    """Exact minimum of an affine objective (c0, ca, cb) over {(a, b): c0 + ca a + cb b >= 0 for each
    half-plane}, by enumerating the vertices (intersections of two boundary lines that satisfy all)."""
    best = None
    for (p0, pa, pb), (q0, qa, qb) in itertools.combinations(halfplanes, 2):
        det = pa * qb - pb * qa
        if det == 0:
            continue
        a = (-p0 * qb + pb * q0) / det
        b = (-pa * q0 + p0 * qa) / det
        if all(c0 + ca * a + cb * b >= 0 for c0, ca, cb in halfplanes):
            value = objective[0] + objective[1] * a + objective[2] * b
            best = value if best is None or value < best else best
    return best


def phi10_polygon(p):
    """(a, b) of phi10: tau/4 <= a <= sigma - theta, (tau - a)/3 <= b <= min(a, sigma - a) (regions.py)."""
    return [(-p.tau / 4, F(1), F(0)), (p.sigma - p.theta, F(-1), F(0)), (-p.tau / 3, F(1, 3), F(1)),
            (F(0), F(1), F(-1)), (p.sigma, F(-1), F(-1))]


def section6_polygon(p):
    """(a, b) of HL's Section 6 (phi11): tau/5 <= a <= sigma - 2 theta, (tau - a)/4 <= b <= min(a, sigma - theta - a)."""
    return [(-p.tau / 5, F(1), F(0)), (p.sigma - 2 * p.theta, F(-1), F(0)), (-p.tau / 4, F(1, 4), F(1)),
            (F(0), F(1), F(-1)), (p.sigma - p.theta, F(-1), F(-1))]


def _g(a, b, c):
    return {"0.1129": a, "0.1142": b, "0.119": c}


CHECKS = [
    # Ranges of theta.
    Check("gaps 3.3", "theta > 1/10 (Lemma 14'; tau > 3 sigma/2)", ">", lambda p: p.theta - F(1, 10),
          _g("0.0129000", "0.0142000", "0.0190000")),
    Check("gaps 3.3", "theta < 1/8 (Lemma 11'; nu + rho < 1)", ">", lambda p: F(1, 8) - p.theta,
          _g("0.0121000", "0.0108000", "0.0060000")),
    Check("gaps 3.3", "theta < 1/6 (Lemma A3)", ">", lambda p: F(1, 6) - p.theta,
          _g("0.0537667", "0.0524667", "0.0476667")),
    Check("gaps 3.3", "theta > 1/14 (Lemma 7)", ">", lambda p: p.theta - F(1, 14),
          _g("0.0414714", "0.0427714", "0.0475714")),
    Check("R* 7", "theta > 1/9", ">", lambda p: p.theta - F(1, 9)),
    Check("R* 7", "theta >= 3/26 (HL's own range; replaced below it by the strip argument, R* 6)", "info",
          lambda p: p.theta - F(3, 26)),
    # gaps_v3 Section 2.6: R*'s use of Lemma 13', and HL's uses of Lemma 13.
    Check("gaps 2.6", "R*: PQ >= W, (1/3 + tau)/2 - sigma = 3 theta/2 - 1/12", ">",
          lambda p: (F(1, 3) + p.tau) / 2 - p.sigma, _g("0.0860167", "0.0879667", "0.0951667")),
    Check("gaps 2.6", "R*: PQ <= X^{tau-eps}, w >= X^eps: tau - sigma - rho/2 = (26 theta - 1)/24", ">",
          lambda p: p.tau - p.sigma - p.rho / 2, _g("0.0806417", "0.0820500", "0.0872500")),
    Check("gaps 2.6", "R*: P <= X^{zeta-eps}, zeta - sigma = (14 theta - 1)/6", ">",
          lambda p: p.zeta - p.sigma, _g("0.0967667", "0.0998000", "0.1110000")),
    Check("gaps 2.6", "R* non-empty: rho/2 - (zeta - sigma)", ">",
          lambda p: p.rho / 2 - (p.zeta - p.sigma), _g("0.0483917", "0.0465500", "0.0397500")),
    Check("gaps 2.6", "R*: PQ^2 <= X^{nu-eps} from 2b <= rho: sigma + rho - nu", "=",
          lambda p: p.sigma + p.rho - p.nu),
    Check("gaps 2.6", "R*: Q >= X^eps, (tau - sigma)/2 - theta", "=", lambda p: (p.tau - p.sigma) / 2 - p.theta),
    Check("gaps 2.6", "R*: Q <= 2P, b <= rho/2 < 1/3 <= a: 1/3 - rho/2", ">", lambda p: F(1, 3) - p.rho / 2),
    Check("gaps 2.6", "Lemma 13: Q <= X^{rho-eps}, rho - zeta/2 = (14 theta - 1)/12", ">",
          lambda p: p.rho - p.zeta / 2, _g("0.0484", "0.0499", "0.0555")),
    Check("gaps 2.6", "Lemma 13: PQ <= X^{tau-eps}, tau - zeta = (1 - 2 theta)/6", ">",
          lambda p: p.tau - p.zeta, _g("0.1290", "0.1286", "0.1270")),
    # gaps_v3 Section 3.3: Lemma 14' (U1: j = 0, phi10; U2: j = 1, Section 6).
    Check("gaps 3.3", "j = 0: products of (i) below X^tau, tau - 3 sigma/2 = 5 theta/2 - 1/4", ">",
          lambda p: p.tau - 3 * p.sigma / 2, _g("0.0322500", "0.0355000", "0.0475000")),
    Check("gaps 3.3", "j = 0: |pqr|^2 >= W (ii), (sigma - theta + tau)/2 - sigma = theta/2", ">",
          lambda p: (p.sigma - p.theta + p.tau) / 2 - p.sigma, _g("0.0564500", "0.0571000", "0.0595000")),
    Check("gaps 3.3", "j = 1: products of (i) below X^tau, tau - 3(sigma - theta)/2 = 4 theta - 1/4", ">",
          lambda p: p.tau - 3 * (p.sigma - p.theta) / 2, _g("0.2016000", "0.2068000", "0.2260000")),
    Check("gaps 3.3", "j = 1: |pqr|^2 >= W X^-theta (ii), (2(sigma - 2 theta) + tau)/3 - (sigma - theta) = theta/3",
          ">", lambda p: (2 * (p.sigma - 2 * p.theta) + p.tau) / 3 - (p.sigma - p.theta),
          _g("0.0376333", "0.0380667", "0.0396667")),
    Check("gaps 3.3", "(U2): |pq|^2 >= W X^-2theta on Section 6, 2 tau/5 - (sigma - 2 theta)", ">",
          lambda p: 2 * p.tau / 5 - (p.sigma - 2 * p.theta), _g("0.0838600", "0.0882800", "0.1046000")),
    Check("gaps 3.3", "(U2): the same as min over Section 6 of a + b - (sigma - 2 theta) (LP)", ">",
          lambda p: polygon_min(section6_polygon(p), (-(p.sigma - 2 * p.theta), F(1), F(1))),
          _g("0.0838600", "0.0882800", "0.1046000")),
    Check("gaps 3.3", "(U2): j(p) = 2 on Section 6, tau/5 - (sigma - 3 theta)", ">",
          lambda p: p.tau / 5 - (p.sigma - 3 * p.theta), _g("0.0741800", "0.0796400", "0.0998000")),
    Check("gaps 3.3", "(U1): min over phi10 of a + b - (sigma - theta) (LP) = 5 theta/2 - 1/4", ">",
          lambda p: polygon_min(phi10_polygon(p), (-(p.sigma - p.theta), F(1), F(1))),
          _g("0.0322500", "0.0355000", "0.0475000")),
    Check("gaps 3.3", "(U1): Q >= X^eps, min over phi10 of b - theta (LP)", ">=",
          lambda p: polygon_min(phi10_polygon(p), (-p.theta, F(0), F(1)))),
    Check("gaps 3.3", "(U2): Q >= X^eps, min over Section 6 of b - theta (LP)", ">=",
          lambda p: polygon_min(section6_polygon(p), (-p.theta, F(0), F(1)))),
    Check("gaps 3.3", "(U1): Q <= 2P, min over phi10 of a - b (LP)", ">=",
          lambda p: polygon_min(phi10_polygon(p), (F(0), F(1), F(-1)))),
    Check("gaps 3.3", "(U2): Q <= 2P, min over Section 6 of a - b (LP)", ">=",
          lambda p: polygon_min(section6_polygon(p), (F(0), F(1), F(-1)))),
    # rstar_claim.md Section 7: HL's Section 5 decomposition, Lemmas 7, 11, 13, 14.
    Check("R* 7", "Lemma 7: zeta - (1 - 4 sigma/3)", "=", lambda p: p.zeta - (1 - 4 * p.sigma / 3)),
    Check("R* 7", "Lemma 7: rho - (1 - 11 sigma/6)", "=", lambda p: p.rho - (1 - 11 * p.sigma / 6)),
    Check("R* 7", "Lemma 7: nu - (1 - 5 sigma/6)", "=", lambda p: p.nu - (1 - 5 * p.sigma / 6)),
    Check("R* 3", "Lemma 7, max form = HL form: rho - sigma/2", ">=", lambda p: p.rho - p.sigma / 2),
    Check("R* 3", "Lemma 7, max form = HL form: nu - 3 sigma/2", ">=", lambda p: p.nu - 3 * p.sigma / 2),
    Check("R* 7", "Lemma 11: chi - (nu + 2 theta)/2", "=", lambda p: p.chi - (p.nu + 2 * p.theta) / 2),
    Check("R* 7", "Lemma 11: chi - sigma", ">=", lambda p: p.chi - p.sigma),
    Check("R* 7", "Lemma 11: zeta - chi", ">", lambda p: p.zeta - p.chi),
    Check("R* 7", "Lemma 11: 1/2 - chi", ">", lambda p: F(1, 2) - p.chi),
    Check("R* 7", "Lemma 11: (1 - sigma)/2 - (nu - sigma)", ">=", lambda p: (1 - p.sigma) / 2 - (p.nu - p.sigma)),
    Check("R* 7", "phi1 piece 2: 1/4 + 3 theta/2 - sigma", ">=", lambda p: F(1, 4) + 3 * p.theta / 2 - p.sigma),
    Check("R* 7", "phi1 piece 2: chi - (1/4 + 3 theta/2)", ">=", lambda p: p.chi - (F(1, 4) + 3 * p.theta / 2)),
    Check("R* 7", "Section 2, j = 0: 1/3 - (sigma - theta)", ">=", lambda p: F(1, 3) - (p.sigma - p.theta)),
    Check("R* 7", "Section 2, Lemma 13: nu - (2 zeta - 1/3)", ">=", lambda p: p.nu - (2 * p.zeta - F(1, 3))),
    Check("R* 7", "D2: (1 + 10 theta)/6 - 1/3 (also phi4's lower limit >= z)", ">=",
          lambda p: (1 + 10 * p.theta) / 6 - F(1, 3)),
    Check("R* 7", "D2: sigma - (1 + 10 theta)/6", ">=", lambda p: p.sigma - (1 + 10 * p.theta) / 6),
    Check("R* 7", "1/2 - zeta", ">", lambda p: F(1, 2) - p.zeta),
    Check("R* 7", "Section 3, j = 0: rho - (sigma - theta)", ">=", lambda p: p.rho - (p.sigma - p.theta)),
    Check("R* 7", "Section 3: 1/3 - rho", ">", lambda p: F(1, 3) - p.rho),
    Check("R* 7", "(rho + tau)/2 - chi", "=", lambda p: (p.rho + p.tau) / 2 - p.chi),
    Check("R* 7", "2 zeta - nu - rho", "=", lambda p: 2 * p.zeta - p.nu - p.rho),
    Check("R* 7", "Section 4: (sigma - theta) - (1/2 - 3 theta)", ">=", lambda p: p.sigma - p.theta - (F(1, 2) - 3 * p.theta)),
    Check("R* 7", "Section 4: (sigma - theta) - nu/3", ">", lambda p: p.sigma - p.theta - p.nu / 3),
    Check("R* 7", "phi8: nu/3 - (2 sigma - nu)", ">=", lambda p: p.nu / 3 - (2 * p.sigma - p.nu)),
    Check("R* 7", "phi9: tau/3 - (1/2 - 3 theta)", ">=", lambda p: p.tau / 3 - (F(1, 2) - 3 * p.theta)),
    Check("R* 7", "Lemma 14 at the top of phi10: (tau - (sigma - theta))/3 - theta", "=",
          lambda p: (p.tau - (p.sigma - p.theta)) / 3 - p.theta),
    Check("R* 7", "phi10 case: tau/4 - (1 - 7 theta)/2", ">=", lambda p: p.tau / 4 - (1 - 7 * p.theta) / 2),
    Check("R* 6", "Section 6: j(pqr) = 0, 3 tau/5 - (sigma - theta)", ">=", lambda p: 3 * p.tau / 5 - (p.sigma - p.theta)),
    Check("R* 7", "strip: sigma/2 - (sigma - 2 theta)", ">=", lambda p: p.sigma / 2 - (p.sigma - 2 * p.theta)),
    Check("R* 7", "strip: nu/3 - (sigma - 2 theta)", ">=", lambda p: p.nu / 3 - (p.sigma - 2 * p.theta)),
    Check("R* 7", "strip: tau/3 - (sigma - 2 theta)", ">=", lambda p: p.tau / 3 - (p.sigma - 2 * p.theta)),
    Check("R* 7", "R*: (1 - sigma)/2 - rho/2", ">=", lambda p: (1 - p.sigma) / 2 - p.rho / 2),
]

STATUS = {">": lambda m: m > 0, ">=": lambda m: m >= 0, "=": lambda m: m == 0, "info": lambda m: True}


def linear_form(check):
    """(c1, c0) with margin = c1 theta + c0 at 0.1129 and 0.119, and whether this holds at 0.1133 and 0.1142."""
    t1, t2 = F("0.1129"), F("0.119")
    m1, m2 = check.margin(params(t1)), check.margin(params(t2))
    c1 = (m2 - m1) / (t2 - t1)
    c0 = m1 - c1 * t1
    return c1, c0, all(check.margin(params(t)) == c1 * t + c0 for t in THETAS)


def evaluate():
    """One row per check: values at THETAS, status, linear form, agreement with the printed gaps_v3 value."""
    rows = []
    for c in CHECKS:
        values = {t: c.margin(params(t)) for t in THETAS}
        c1, c0, affine = linear_form(c)
        printed = {}
        for t_str, s in c.gaps.items():
            half_unit = F(1, 2) * F(10) ** -(len(s.split(".")[1]))
            printed[t_str] = (s, abs(values[F(t_str)] - F(s)) <= half_unit)
        rows.append(dict(check=c, values=values, ok={t: STATUS[c.kind](m) for t, m in values.items()},
                         c1=c1, c0=c0, affine=affine, printed=printed))
    return rows


def _form(c1, c0):
    if c1 == 0:
        return f"{c0}"
    slope = {1: "theta", -1: "-theta"}.get(c1, f"{c1} theta")
    return slope if c0 == 0 else f"{slope} {'+' if c0 > 0 else '-'} {abs(c0)}"


def report(rows):
    """Plain-text table (all four theta values) and a list of failures."""
    head = f"{'source':9s} {'kind':4s} {'linear form':24s} " + " ".join(f"{str(float(t)):>11s}" for t in THETAS)
    lines = [head, "-" * len(head)]
    failures = []
    for r in rows:
        c = r["check"]
        vals = " ".join(f"{float(r['values'][t]):+11.7f}" for t in THETAS)
        lines.append(f"{c.source:9s} {c.kind:4s} {_form(r['c1'], r['c0']):24s} {vals}   {c.name}")
        extra = []
        if not r["affine"]:
            extra.append("not affine on the four values")
        for t_str, (s, match) in r["printed"].items():
            extra.append(f"gaps_v3 at {t_str}: {s} {'agrees' if match else 'DISAGREES'}")
            if not match:
                failures.append(f"{c.name}: gaps_v3 value {s} at {t_str}")
        if extra:
            lines.append(" " * 40 + "; ".join(extra))
        for t, ok in r["ok"].items():
            if not ok:
                failures.append(f"{c.name}: fails at theta = {t} (margin {float(r['values'][t]):+.7f})")
        if not r["affine"]:
            failures.append(f"{c.name}: not affine")
    return "\n".join(lines), failures


def exact_table(rows):
    """The table requested for the paper: linear form and exact values at 0.1129 and 0.1133."""
    t1, t2 = F("0.1129"), F("0.1133")
    lines = [f"{'source':9s} {'kind':4s} {'linear form':24s} {'at 0.1129 (exact = decimal)':36s} "
             f"{'at 0.1133 (exact = decimal)':36s} check"]
    for r in rows:
        c = r["check"]
        cells = [f"{str(r['values'][t]):>22s} = {float(r['values'][t]):+.7f}" for t in (t1, t2)]
        lines.append(f"{c.source:9s} {c.kind:4s} {_form(r['c1'], r['c0']):24s} {cells[0]:36s} {cells[1]:36s} {c.name}")
    return "\n".join(lines)


if __name__ == "__main__":
    rows = evaluate()
    text, failures = report(rows)
    print(text, flush=True)
    print(f"\n{len(rows)} checks; failures: {failures if failures else 'none'}", flush=True)
