# Code: the numbers of Section 10 and Appendix D

This folder computes every number of Section 10 and Appendix D of the paper (`../Paper/main.tex`), in Julia. Labels follow the paper; table, display and lemma numbers are those of its current build.

## Files

- `buchstab.jl`: Buchstab's function ω and Ω(u) = ∫₁ᵘ ω. On each unit interval, Chebyshev series of degree 40 are built from the relation that defines ω; ω, ω′ and Ω are tabulated on a grid of step 10⁻⁵ and evaluated by cubic Hermite interpolation. Includes the accuracy check of Section 10.
- `regions.jl`: the parameters (4); the regions D1, D2, D4, D6, D8 and E3, ER, E5, E7, E9, E10, E11 of Section 7 (Table 1, display (42)); the region R (15), the sets G, I, L (16) and the polytope P (17), pointwise and as unions of δ-intervals. In the decomposition of Harman and Lewis, where R is not grouped, the region of Range 2 is D2 ∪ R, labelled `D2+R`.
- `losses.jl`: the losses (18) by the two methods of Section 10, and the θ at which the total loss is 1 (Table 3).
  - First method: the innermost integral exact through Ω; composite Gauss–Legendre in the outer variables, split at the kinks of the region limits and at the planes where G ∪ L or G ∪ P changes its formula; k = 10 nodes per piece, 256 panels on the α-axis for the regions D_k, 16 panels per axis for the regions E_k.
  - Second method: nested adaptive Gauss–Kronrod quadrature (21-point rule, at most 200 subintervals per piece) for the regions D_k; for the regions E_k, the mean of 8 scrambled Sobol estimates with 2²⁰ points each (Joe–Kuo direction numbers, Matoušek's linear matrix scramble and a digital shift, scrambling r from the generator `Xoshiro(r)`), with the sets tested pointwise.
- `inequalities.jl`: Table 4 (Appendix D) in exact rational arithmetic.
- `run_all.jl`: runs everything, writes `results/`, and compares every number of Section 10 with the number printed in the paper.
- `figure_data.jl`: the data of the figures of the paper, from the same definitions and the first method: `figures/domain.json` (the domain D at θ = 0.1129, its Ranges 1–6, the regions D_k and R) and `figures/total_loss.json` (the total loss of each decomposition of Table 3 for 0.112 ≤ θ ≤ 0.118).
- `figures/plot_figures.py`: draws `figures/domain.pdf` and `figures/total_loss.pdf` (Figures 1 and 2 of the paper) from these data with matplotlib.
- `Project.toml`, `Manifest.toml`: the Julia environment.

Every computation is deterministic: the first method sums the contributions in a fixed order, independent of the number of threads, and the scramblings of the second method are fixed by their seeds.

## Results

| output | paper | result |
|---|---|---|
| `buchstab.log` | Section 10 | ω and Ω agree with their closed forms on [1, 3] to 5.9 × 10⁻¹⁵, and satisfy the relation defining ω on [2, 40) to 1.5 × 10⁻¹³ |
| `table_inequalities.txt` | Table 4 | all 22 rows: the margin of each inequality of a row is a positive multiple of its quantity, the quantity is positive at θ = 0.1129 and 0.1133, and its value at 0.1129 is the printed value |
| `table_hl_0119.txt` | Section 10, Table 6 | the decomposition of Harman and Lewis at θ = 0.119, region by region against their bounds: total 0.9425585, below their bound 0.943 |
| `table_roots.txt`, `roots.log` | Table 3 | θ = 0.1167825, 0.1141073, 0.1139312, 0.1132049, 0.1128554; the quantities of Table 4 are positive at each |
| `table_losses_0.1129.txt` | Table 2 | first method: every region, and Φ(0.1129) = 0.9990033; second method: every region with its error estimate, and Φ(0.1129) = 0.9990044, error estimate 1.0 × 10⁻⁶ |
| `table_losses_P_0.1133.txt` | Section 10 | Φ_P(0.1133) = 0.9978243 (first method) and 0.9978250 (second method) |
| `table_convergence_0.1129.txt` | Section 10 | Φ(0.1129) by the rule of Table 2 differs from the finest rule (k = 12, 32 panels per axis, 512 panels on the α-axis for the D_k) by 2.6 × 10⁻¹⁰, and each loss by at most 2.5 × 10⁻⁹; over the 8 coarser rules the difference of Φ(0.1129) decreases from 8.7 × 10⁻⁶ to 6.8 × 10⁻¹¹ |
| `paper_numbers.txt` | Section 10, Tables 2, 3 and 6 | the 87 numbers of Section 10 and of Tables 2, 3 and 6, as printed and as computed: all agree |
| `summary.json` | | every value above, the runtimes and the versions |

The two methods differ by 1.2 × 10⁻⁶ in Φ(0.1129), 1.1 times the error estimate of the second method; the margin 1 − Φ(0.1129) = 9.97 × 10⁻⁴ is 855 times the largest of the change under refinement, the error estimate and this difference.

## Requirements

Julia 1.13.0 with the packages of `Project.toml` (FastGaussQuadrature, JSON, PolyLog, QuadGK, Roots); `Manifest.toml` fixes their versions. The figures need Python 3 with matplotlib (drawn with Python 3.14.7 and matplotlib 3.11.2).

```
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Reproduce

From this folder:

```
julia --project=. -t 32,1 run_all.jl                  # everything; 37 s with 32 threads (summary.json)
julia --project=. buchstab.jl                         # the accuracy check
julia --project=. inequalities.jl                     # Table 4
julia --project=. -t 32 losses.jl 0.1129 HL+R+G+L     # the first method, region by region
julia --project=. -t 32 figure_data.jl                # the data of the figures; about 1 min
python figures/plot_figures.py                        # the figures
```

The run used 32 threads, with a peak memory of 1.8 GB. Step times from `summary.json`: the decomposition of Harman and Lewis at 0.119, 25 s, most of it compilation; Table 3, 7 s; the convergence, 4 s; every other step under 1 s.
