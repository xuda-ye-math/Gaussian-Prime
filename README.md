# Infinitely Many Primes $p = m^2 + n^2$ with $n < p^{0.1129}$

Xuda Ye

Harman and Lewis proved that infinitely many primes $p = m^2 + n^2$ have $n < p^{0.119}$. This project improves the exponent to $0.1129$, so that infinitely many primes also satisfy $\{\sqrt{p}\} < p^{-0.2742}$. Both results follow from a lower bound for the number of Gaussian primes of norm between $X/2$ and $X$ in a sector of angle $X^{\theta - 1/2}$, proved with Harman's sieve at $\theta = 0.1129$.

## Contents

| Folder | Contents |
|---|---|
| [`Paper/`](Paper) | the manuscript `main.tex`, the bibliography `references.bib` and the build `main.pdf` |
| [`Code/`](Code) | Julia computation of every number of Section 10 and Appendix D of the paper: the losses by both methods of Section 10, Tables 2–4, and a comparison of each number with the paper; outputs in `Code/results/` |
| [`Lean/`](Lean) | Lean 4 formalization: the deduction of Theorem 1 from Theorem 2 and the parameter table of Appendix D are proved; the plan for the remaining parts is in `Lean/README.md` |

`Code/README.md` and `Lean/README.md` give the details of each part.

## Reproduce

Code, with Julia 1.13 (37 s with 32 threads):

```
cd Code
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -t 32,1 run_all.jl
```

Lean, with Lean 4 v4.33.1 and Mathlib v4.33.1:

```
cd Lean
lake build
```

## Status

The manuscript is under revision.
