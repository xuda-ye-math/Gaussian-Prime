# Infinitely Many Primes $p = m^2 + n^2$ with $n < p^{0.1129}$

Xuda Ye

Harman and Lewis proved that infinitely many primes $p = m^2 + n^2$ have $n < p^{0.119}$. This project improves the exponent to $0.1129$, so that infinitely many primes also satisfy $\{\sqrt{p}\} < p^{-0.2742}$. Both results follow from a lower bound for the number of Gaussian primes of norm between $X/2$ and $X$ in a sector of angle $X^{\theta - 1/2}$, proved with Harman's sieve at $\theta = 0.1129$.

## Contents

| Folder | Contents |
|---|---|
| [`Paper/`](Paper) | the manuscript `main.tex`, the bibliography `references.bib` and the build `main.pdf` |
| [`Code/`](Code) | Python evaluation of the sieve losses: Buchstab's function, the loss regions, the loss integrals and $\theta_{\min}$, and the structural inequalities in exact rational arithmetic; outputs in `Code/results/` |
| [`Lean/`](Lean) | Lean 4 formalization: the deduction of Theorem 2 from Theorem 1 and the parameter table of Appendix D are proved; the plan for the remaining parts is in `Lean/README.md` |

`Code/README.md` and `Lean/README.md` give the details of each part.

## Reproduce

Code, with Python 3, `numpy` and `scipy` (about 19 min on 32 cores):

```
cd Code
pip install numpy scipy
python -u run_all.py
```

Lean, with Lean 4 v4.33.1 and Mathlib v4.33.1:

```
cd Lean
lake build
```

## Status

The manuscript is under revision.
