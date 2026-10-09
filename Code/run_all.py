"""Driver: every computation of the code companion, written to Codes/results/.

    python -u run_all.py

Steps and outputs (all plain text, plus results/summary.json; progress in results/run_all.log):
  buchstab        buchstab.log                      accuracy of omega and Omega
  checks          checks.log, table_structural.txt   structural inequalities and lemma margins (exact)
  classification  classification.log                 original against reduced rescue sets, and the
                                                     reduced sets against their d-interval form
  hl_0119         table_hl_0119.txt                  HL's decomposition at theta = 0.119 against HL's bounds
  roots           roots.log, table_thetamin.txt      theta_min of every decomposition
  losses          table_losses_0.1129_G4L.txt, table_losses_0.1133_G4P.txt, table_losses_0.1142_Rstar.txt
                                                     per-region losses, both integrators
  variants        variants.log                       single Lemma 7 split; HL's cut in phi3, phi9
  convergence     table_convergence_0.1129.txt       primary integrator under refinement
Every reference value of the companion (REFERENCES) is compared with the computed value in summary.json.
"""
import json
import platform
import threading
import time
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from pathlib import Path

import numpy as np
import scipy

import buchstab
import checks
import losses as L

OUT = Path(__file__).resolve().parent / "results"
RES = dict(k=10, panels=16)                 # production resolution of the primary integrator
QMC = dict(m=20, reps=8)                    # cross-check: 8 scramblings of 2^20 Sobol points
LADDER = [(6, 2), (6, 4), (8, 4), (8, 8), (10, 8), (8, 16), (10, 16), (12, 16), (12, 32)]
TOL = 1e-6
BRACKETS = {"HL": (0.116, 0.119), "HL printed": (0.116, 0.119), "R*": (0.1115, 0.1175),
            "R*+G4": (0.1115, 0.1175), "R*+G4+P": (0.1115, 0.1175), "R*+G4+L": (0.1115, 0.1175)}

# Reference values: (value, source).
RV2 = "routes/R1-review2/rv2_totals.log"
REFERENCES = {
    "total HL at 0.119": (0.942558, "routes/R1/final_0119.log; rstar_claim.md Section 4"),
    "theta_min HL": (0.1167824, "routes/R1/thetamin.log; rstar_claim.md Section 8"),
    "theta_min HL printed": (0.1172286, "routes/R1/thetamin.log"),
    "theta_min R*": (0.1141073, "rstar_claim.md Section 8"),
    "total R* at 0.1142": (0.9977385, "routes/R1-review/rv_01142_p80.log"),
    "theta_min R*+G4": (0.1139313, "rstar2_claim.md Section 0"),
    "theta_min R*+G4+P": (0.1132050, "rstar2_claim.md Section 0"),
    "total R*+G4+P at 0.1133": (0.9978243, RV2),
    "1 - total R*+G4+P at 0.1133": (0.0021757, RV2),
    "theta_min R*+G4+L": (0.1128555, "rstar2_claim.md Section 0"),
    "total R*+G4+L at 0.1129": (0.9990032, RV2),
    "1 - total R*+G4+L at 0.1129": (0.0009968, RV2),
    "total R*+G4+L[qrs|p] at 0.1129": (0.9990035, "routes/R1-review2/rv2_sens2.log (resolution 16)"),
}
REGIONS_0129 = {   # R* + G4 + L at 0.1129: RV2 (rv2_totals.log), and rstar2_claim.md Section 6, integrators A and B
    "RV2": dict(phi1=0.3201605, **{"D2-R*": 0.3437028}, phi4=0.2132922, phi6=0.0480512, phi8=0.0684437,
                phi3=0.0000098, phi5=0.0006199, phi7=0.0004451, phi9=0.0002727, phi10=0.0028858,
                phi11=0.0000988, T3=0.0010206),
    "claim A": dict(phi1=0.3201605, **{"D2-R*": 0.3437028}, phi4=0.2132922, phi6=0.0480512, phi8=0.0684437,
                    phi3=0.0000098, phi5=0.0006207, phi7=0.0004465, phi9=0.0002733, phi10=0.0028859,
                    phi11=0.0000983, T3=0.0010216),
    "claim B": dict(phi3=0.0000098, phi5=0.0006200, phi7=0.0004451, phi9=0.0002727, phi10=0.0028851,
                    phi11=0.0000988, T3=0.0010206),
}
REGIONS_0133 = {   # R* + G4 + P at 0.1133: RV2 (rv2_totals.log, resolution 80)
    "RV2": dict(phi1=0.3171641, **{"D2-R*": 0.3398803}, phi4=0.2094263, phi6=0.0524427, phi8=0.0659758,
                phi3=0.0000111, phi5=0.0006172, phi7=0.0004751, phi9=0.0002408, phi10=0.0104998,
                phi11=0.0000862, T3=0.0010050),
}
REGIONS_0142 = {   # R* at 0.1142: first referee (rv_01142_p80.log)
    "RV1": dict(phi1=0.31043458, **{"D2-R*": 0.33136148}, phi4=0.20076836, phi6=0.06220925, phi8=0.06058557,
                phi3=0.00002911, phi5=0.00128465, phi7=0.00059226, phi9=0.00203658, phi10=0.02578030,
                phi11=0.00029967, T3=0.00235665),
}
HL_PRINTED = dict(phi1=0.275, phi2=0.341, phi3=0.001, phi4=0.156, phi5=0.001, phi6=0.112, phi7=0.001,
                  phi8=0.036, phi9=0.001, phi10=0.017, phi11=0.001)   # HL pp. 132-134; total < 0.943
R1_0119 = dict(phi1=0.274784, phi2=0.345646, phi3=0.000092, phi4=0.155456, phi5=0.001186, phi6=0.111625,
               phi7=0.000885, phi8=0.035639, phi9=0.001075, phi10=0.016142, phi11=0.000028)  # final_0119.log

LOG = OUT / "run_all.log"


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as f:
        f.write(line + "\n")


class Heartbeat:
    """Prints a progress line every 20 s while a step runs."""
    def __init__(self, label):
        self.label, self.stop = label, threading.Event()

    def __enter__(self):
        self.t0 = time.time()
        def beat():
            while not self.stop.wait(20):
                log(f"  ... {self.label}: running, {time.time() - self.t0:.0f}s")
        threading.Thread(target=beat, daemon=True).start()
        return self

    def __exit__(self, *exc):
        self.stop.set()


def write(name, text):
    (OUT / name).write_text(text + "\n")
    log(f"  wrote results/{name}")


def compare(summary, name, value):
    ref, source = REFERENCES[name]
    diff = value - ref
    summary["references"].append(dict(name=name, value=value, reference=ref, source=source, diff=diff,
                                      tolerance=TOL, agrees=abs(diff) <= TOL))
    log(f"  {name}: {value:.7f}  reference {ref:.7f}  diff {diff:+.1e}  {'agrees' if abs(diff) <= TOL else 'DIFFERS'}")


def loss_table(title, gl, check, refs):
    """Per-region table: primary, cross-check with s.e., their difference, and reference columns."""
    names = list(gl)
    head = f"{'region':7s} {'primary':>11s} {'cross-check':>11s} {'s.e.':>8s} {'diff':>9s}" + "".join(
        f" {r:>11s} {'diff':>9s}" for r in refs)
    lines = [title, "", head, "-" * len(head)]
    for n in names:
        v, (c, e) = gl[n], check[n]
        row = f"{n:7s} {v:11.7f} {c:11.7f} {e:8.1e} {v - c:+9.1e}"
        for r in refs.values():
            row += f" {r[n]:11.7f} {v - r[n]:+9.1e}" if n in r else f" {'':11s} {'':9s}"
        lines.append(row)
    tot, ctot = sum(gl.values()), sum(c for c, _ in check.values())
    se = float(np.sqrt(sum(e * e for _, e in check.values())))
    lines.append(f"{'TOTAL':7s} {tot:11.7f} {ctot:11.7f} {se:8.1e} {tot - ctot:+9.1e}")
    lines.append(f"{'1-TOTAL':7s} {1 - tot:11.7f} {1 - ctot:11.7f}")
    return "\n".join(lines)


def variants(pool, full):
    """The single Lemma 7 split U = qrs, V = p at 0.1129 (against `full`, all splits), and HL's cut
    g <= a + b - 2 theta in phi3 and phi9 (regions.regions_4d, cut=False)."""
    single = L.losses_gl(0.1129, "R*+G4+L[qrs|p]", **RES, pool=pool)
    lines = ["Single Lemma 7 split: I restricted to Lemma 6 and the split U = qrs, V = p, at theta = 0.1129.",
             f"  total with all splits      {sum(full.values()):.9f}",
             f"  total with U = qrs, V = p  {sum(single.values()):.9f}",
             "  regions that differ: " + (", ".join(n for n in full if abs(full[n] - single[n]) > 1e-12) or "none"),
             "", "HL's cut g <= a + b - 2 theta in phi3 and phi9: losses with and without it."]
    cut = {}
    for theta, dec in [(0.119, "HL"), (0.1129, "R*"), (0.1129, "R*+G4+L")]:
        on = L.losses_gl(theta, dec, **RES, pool=pool)
        off = L.losses_gl(theta, dec, **RES, pool=pool, cut=False)
        cut[f"{dec} at {theta}"] = {n: (on[n], off[n]) for n in ("phi3", "phi9")}
        for n in ("phi3", "phi9"):
            lines.append(f"  {dec:8s} theta = {theta:<7} {n}: with cut {on[n]:.12f}  without {off[n]:.12f}  "
                         f"difference {off[n] - on[n]:+.1e}")
    return lines, single, cut


def main():
    OUT.mkdir(exist_ok=True)
    LOG.write_text("")
    summary = dict(generated=time.strftime("%Y-%m-%d %H:%M:%S"), python=platform.python_version(),
                   numpy=np.__version__, scipy=scipy.__version__, resolution=RES, qmc=QMC, tolerance=TOL,
                   references=[], runtimes={})
    t_all = time.time()
    with ProcessPoolExecutor() as pool:
        # 1. Buchstab's function.
        t0 = time.time()
        log("step buchstab")
        b = buchstab.check()
        write("buchstab.log", "\n".join(f"{k:46s} {v:.3e}" if abs(v) < 1e-3 else f"{k:46s} {v:.15f}"
                                        for k, v in b.items()))
        summary["buchstab"] = b
        summary["runtimes"]["buchstab"] = time.time() - t0

        # 2. Structural inequalities and margins.
        t0 = time.time()
        log("step checks")
        rows = checks.evaluate()
        text, failures = checks.report(rows)
        write("checks.log", text + f"\n\n{len(rows)} checks; failures: {failures if failures else 'none'}")
        write("table_structural.txt", "Structural inequalities and lemma margins: exact linear form in theta, "
              "exact value and decimal at 0.1129 and 0.1133.\nkind: '>' strict, '>=' non-strict, '=' identity, "
              "'info' not a condition.\n\n" + checks.exact_table(rows))
        root_checks = {}
        summary["checks"] = dict(count=len(rows), failures=failures, rows=[
            dict(source=r["check"].source, name=r["check"].name, kind=r["check"].kind,
                 linear_form=[str(r["c1"]), str(r["c0"])],
                 values={str(t): float(v) for t, v in r["values"].items()}) for r in rows])
        log(f"  {len(rows)} checks, failures: {len(failures)}")
        summary["runtimes"]["checks"] = time.time() - t0

        # 3. Brute-force comparison of the rescue-set definitions.
        t0 = time.time()
        log("step classification (10^6 random points per region and theta)")
        lines, cls = [], {}
        with Heartbeat("classification"):
            for theta in (0.1129, 0.1133):
                res = L.classification_mismatches(theta, pool=pool)
                cls[str(theta)] = res
                lines.append(f"theta = {theta}   (comparisons, mismatches)")
                for reg, tests in res.items():
                    lines.append(f"  {reg:6s} " + "  ".join(f"{t} {n} {m}" for t, (n, m) in tests.items()))
        total_mis = sum(m for r in cls.values() for tests in r.values() for _, m in tests.values())
        lines.append(f"\ntotal mismatches: {total_mis}")
        write("classification.log", "\n".join(lines))
        summary["classification"] = dict(results=cls, total_mismatches=total_mis)
        log(f"  total mismatches: {total_mis}")
        summary["runtimes"]["classification"] = time.time() - t0

        # 4. HL's decomposition at 0.119.
        t0 = time.time()
        log("step hl_0119")
        with Heartbeat("hl_0119"):
            d2 = L.losses_gl(0.119, "HL", **RES, pool=pool)
            pr = L.losses_gl(0.119, "HL printed", **RES, pool=pool)
            chk = L.losses_check(0.119, "HL", **QMC, pool=pool)
        lines = ["HL's decomposition at theta = 0.119: exact phi2 region D2 and HL's printed phi2 region,",
                 "against HL's printed bounds (pp. 132-134) and R1's values (routes/R1/final_0119.log).", "",
                 f"{'region':7s} {'D2 primary':>11s} {'cross-check':>11s} {'s.e.':>8s} {'printed phi2':>12s} "
                 f"{'HL bound':>9s} {'R1':>9s}  note"]
        for n in d2:
            note = "exceeds HL's bound" if d2[n] >= HL_PRINTED[n] else ""
            lines.append(f"{n:7s} {d2[n]:11.7f} {chk[n][0]:11.7f} {chk[n][1]:8.1e} {pr[n]:12.7f} "
                         f"{HL_PRINTED[n]:9.3f} {R1_0119[n]:9.6f}  {note}")
        lines.append(f"{'TOTAL':7s} {sum(d2.values()):11.7f} {sum(c for c, _ in chk.values()):11.7f} {'':8s} "
                     f"{sum(pr.values()):12.7f} {'< 0.943':>9s} {sum(R1_0119.values()):9.6f}")
        write("table_hl_0119.txt", "\n".join(lines))
        compare(summary, "total HL at 0.119", sum(d2.values()))
        summary["hl_0119"] = dict(D2=d2, printed=pr, check=chk)
        summary["runtimes"]["hl_0119"] = time.time() - t0

        # 5. theta_min of every decomposition.
        t0 = time.time()
        log("step roots")
        roots, trace = {}, []
        def tlog(msg, **kw):
            trace.append(msg)
            log(msg)
        for dec, bracket in BRACKETS.items():
            t1 = time.time()
            with Heartbeat(f"theta_min {dec}"):
                roots[dec] = L.theta_min(dec, bracket, **RES, pool=pool, log=tlog)
            tlog(f"  theta_min {dec}: {roots[dec]:.7f}  ({time.time() - t1:.0f}s)")
            failing = [r["check"].name for r in rows if r["check"].kind != "info"
                       and not checks.STATUS[r["check"].kind](r["check"].margin(checks.params(Fraction(roots[dec]))))]
            root_checks[dec] = failing
            tlog(f"  structural checks at the root: {'all hold' if not failing else failing}")
        write("roots.log", "\n".join(trace))
        lines = ["theta_min (total loss = 1; brentq, xtol 1e-10; primary integrator at "
                 f"k = {RES['k']}, panels = {RES['panels']}). alpha = 1/2 - 2 theta_min.", "",
                 f"{'decomposition':24s} {'theta_min':>10s} {'alpha':>10s} {'gain in alpha':>14s} "
                 f"{'reference':>10s} {'diff':>9s}"]
        labels = {"HL": "HL, exact D2", "HL printed": "HL, printed phi2", "R*": "+ R*", "R*+G4": "+ R* + G4",
                  "R*+G4+P": "+ R* + G4 + P", "R*+G4+L": "+ R* + G4 + L"}
        for dec, root in roots.items():
            ref = REFERENCES[f"theta_min {dec}"][0]
            gain = "" if dec.startswith("HL") else f"{2 * (roots['HL'] - root):14.7f}"
            lines.append(f"{labels[dec]:24s} {root:10.7f} {0.5 - 2 * root:10.7f} {gain:>14s} {ref:10.7f} "
                         f"{root - ref:+9.1e}")
            compare(summary, f"theta_min {dec}", root)
        lines.append("\ngain in alpha: against HL with the exact D2.  Structural checks (checks.py) at each root: "
                     + ("all hold." if not any(root_checks.values()) else str(root_checks)))
        write("table_thetamin.txt", "\n".join(lines))
        summary["roots"] = dict(values=roots, structural_failures=root_checks)
        summary["runtimes"]["roots"] = time.time() - t0

        # 6. Per-region losses at the three stated values of theta, both integrators.
        t0 = time.time()
        log("step losses")
        summary["losses"] = {}
        for theta, dec, refs, fname in [(0.1129, "R*+G4+L", REGIONS_0129, "table_losses_0.1129_G4L.txt"),
                                        (0.1133, "R*+G4+P", REGIONS_0133, "table_losses_0.1133_G4P.txt"),
                                        (0.1142, "R*", REGIONS_0142, "table_losses_0.1142_Rstar.txt")]:
            with Heartbeat(f"losses {dec} at {theta}"):
                gl = L.losses_gl(theta, dec, **RES, pool=pool)
                chk = L.losses_check(theta, dec, **QMC, pool=pool)
            title = (f"{dec} at theta = {theta}.  Primary: Gauss-Legendre, k = {RES['k']}, panels = {RES['panels']} "
                     f"(reduced definitions).  Cross-check: nested adaptive quadrature (2-dim) and scrambled Sobol "
                     f"QMC, {QMC['reps']} x 2^{QMC['m']} points (4-dim, original definitions of R*2).  "
                     f"Reference columns: " + ", ".join(refs))
            write(fname, loss_table(title, gl, chk, refs))
            tot = sum(gl.values())
            summary["losses"][f"{dec} at {theta}"] = dict(primary=gl, check=chk, total=tot,
                                                         check_total=sum(c for c, _ in chk.values()))
            if dec == "R*":
                compare(summary, "total R* at 0.1142", tot)
            else:
                compare(summary, f"total {dec} at {theta}", tot)
                compare(summary, f"1 - total {dec} at {theta}", 1 - tot)
        summary["runtimes"]["losses"] = time.time() - t0

        # 7. Variants: the single Lemma 7 split U = qrs, V = p; HL's cut in phi3 and phi9.
        t0 = time.time()
        log("step variants")
        with Heartbeat("variants"):
            lines, single, cut = variants(pool, summary["losses"]["R*+G4+L at 0.1129"]["primary"])
        compare(summary, "total R*+G4+L[qrs|p] at 0.1129", sum(single.values()))
        write("variants.log", "\n".join(lines))
        summary["variants"] = dict(single_split=dict(total=sum(single.values()), regions=single), cut=cut)
        summary["runtimes"]["variants"] = time.time() - t0

        # 8. Convergence of the primary integrator at 0.1129 (R* + G4 + L).
        t0 = time.time()
        log("step convergence")
        ladder = {}
        with Heartbeat("convergence"):
            for k, panels in LADDER:
                t1 = time.time()
                ladder[(k, panels)] = L.losses_gl(0.1129, "R*+G4+L", k, panels, pool=pool)
                log(f"  k = {k}, panels = {panels}: total {sum(ladder[(k, panels)].values()):.10f} "
                    f"({time.time() - t1:.0f}s)")
        finest = ladder[LADDER[-1]]
        names = [n for n in finest if n in ("phi3", "phi5", "phi7", "phi9", "phi10", "phi11", "T3")]
        lines = [f"Primary integrator at theta = 0.1129, R* + G4 + L: values, and differences from the finest rule "
                 f"(k = {LADDER[-1][0]}, panels = {LADDER[-1][1]}).  2-dim losses use k nodes and 16 * panels panels on the "
                 f"a-axis.", "",
                 f"{'k':>3s} {'panels':>6s} " + " ".join(f"{n:>12s}" for n in names) + f" {'TOTAL':>14s} {'diff':>9s}"]
        for (k, panels), out in ladder.items():
            tot = sum(out.values())
            lines.append(f"{k:3d} {panels:6d} " + " ".join(f"{out[n]:12.9f}" for n in names)
                         + f" {tot:14.10f} {tot - sum(finest.values()):+9.1e}")
        lines.append("")
        lines.append("differences from the finest rule:")
        for (k, panels), out in ladder.items():
            lines.append(f"{k:3d} {panels:6d} " + " ".join(f"{out[n] - finest[n]:+12.1e}" for n in names))
        prod = sum(ladder[(RES["k"], RES["panels"])].values())
        lines.append(f"\nproduction rule (k = {RES['k']}, panels = {RES['panels']}): total differs from the finest "
                     f"by {prod - sum(finest.values()):+.1e}")
        write("table_convergence_0.1129.txt", "\n".join(lines))
        summary["convergence"] = {f"k={k},panels={p}": out for (k, p), out in ladder.items()}
        summary["runtimes"]["convergence"] = time.time() - t0

    summary["runtimes"]["total"] = time.time() - t_all
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, default=str) + "\n")
    bad = [r["name"] for r in summary["references"] if not r["agrees"]]
    log(f"wrote results/summary.json; references outside {TOL:g}: {bad if bad else 'none'}; "
        f"total runtime {summary['runtimes']['total']:.0f}s")


if __name__ == "__main__":
    main()
