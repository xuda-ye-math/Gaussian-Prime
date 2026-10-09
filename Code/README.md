# Code companion: the Harman–Lewis sieve at θ = 0.1129

## Scripts

- `buchstab.py`: Buchstab's function ω(u) and Ω(u) = ∫₁ᵘ ω, with an accuracy check.
- `regions.py`: HL's parameters σ, τ, ρ, ζ, ν, χ; the loss regions φ1–φ11 (exact D2, corrected φ10), R* and T3; the rescue sets G4, L and P in three forms: as in the claim files, in the reduced form of the paper, and as δ-intervals.
- `losses.py`: the loss integrals and θ_min. Primary integrator: Gauss–Legendre, with the δ-integral exact through Ω, split at the kinks of the region limits and at the planes where the rescue sets jump. At the production setting it uses k = 10 nodes per piece, with 16 panels per axis for the 4-dim losses and 16 · 16 = 256 panels on the α-axis for the 2-dim losses. Cross-check: nested adaptive quadrature for the 2-dim losses, and scrambled Sobol QMC (8 × 2²⁰ points, original definitions) for the 4-dim losses.
- `checks.py`: 59 structural inequalities and lemma margins in exact rational arithmetic.
- `run_all.py`: runs everything and writes `results/`.

## Tests and outcomes

The tolerance is 1e-6. Every reference value listed in `summary.json` agrees within it.

| test | file | result | reference |
|---|---|---|---|
| Buchstab accuracy | `buchstab.log` | closed forms on [1,3] to 1.7e-14; ω(39.9) − e^−γ = 2.2e-16 | – |
| structural inequalities | `checks.log`, `table_structural.txt` | all 59 hold at θ = 0.1129, 0.1133, 0.1142, 0.119; each is an exact linear form in θ | the 54 printed values of gaps_v3.md §2.6, §3.3 all agree; (U1) = 5θ/2 − 1/4 = 0.0322500 at 0.1129 |
| reduced vs original sets | `classification.log` | 0 mismatches; 7 regions, 0.42–1.05 × 10⁶ points each, θ = 0.1129 and 0.1133 | – |
| HL at θ = 0.119 | `table_hl_0119.txt` | total 0.9425585 | 0.942558 (R1) |
| θ_min | `table_thetamin.txt` | HL 0.1167825; HL printed φ2 0.1172286; + R* 0.1141073; + G4 0.1139312; + G4 + P 0.1132049; + G4 + L 0.1128554 | 0.1167824, 0.1172286, 0.1141073, 0.1139313, 0.1132050, 0.1128555 |
| R* at 0.1142 | `table_losses_0.1142_Rstar.txt` | total 0.9977384 | 0.9977385 (RV1) |
| R* + G4 + P at 0.1133 | `table_losses_0.1133_G4P.txt` | total 0.9978243; 1 − total 0.0021757 | 0.9978243 (RV2) |
| R* + G4 + L at 0.1129 | `table_losses_0.1129_G4L.txt` | total 0.9990033; 1 − total 0.0009967; every region within 4.1e-8 of RV2 | 0.9990032 (RV2) |
| integrator cross-check | the three `table_losses_*` | totals differ by 6.7e-7 (0.7 s.e.), 2.1e-7 and 8.9e-8 at 0.1129, 0.1133, 0.1142 | – |
| convergence | `table_convergence_0.1129.txt` | the production rule differs from the finest rule (k = 12, 32 panels) by 2.6e-10 in the total | – |
| single Lemma 7 split U = qrs, V = p | `variants.log` | total 0.999003254, equal to the total with all splits | 0.9990035 (RV2, at its resolution 16) |
| HL's cut γ ≤ α + β − 2θ in φ3, φ9 | `variants.log` | dropping it changes φ3 and φ9 by 0 to 12 digits, at 0.119 (HL) and 0.1129 (R*, R* + G4 + L); proof in `regions.py` | – |

One discrepancy concerns the claim itself. The integrator-A column of rstar2_claim.md §6 at 0.1129 differs from these values by up to 1.4e-6 (φ7: 0.0004465 against 0.0004451), and its total 0.9990065 is 3.2e-6 too high. Its integrator-B column agrees.

## Provenance

The tables and logs in `results/` come from a run that was interrupted at about 22:53, during the convergence step. `table_convergence_0.1129.txt` and `summary.json` come from the complete run just before it. Every number in the interrupted run's tables was checked and is identical to `summary.json`. `run_all.log` is the interrupted run's log. The header of `table_convergence_0.1129.txt` still gives the old wording "(k, 16 panels)" for the 2-dim resolution; the docstring of `losses.py` is corrected.

## Requirements

Python 3 with the packages `numpy` and `scipy`; every other import is from the standard library or from this folder. The results were produced with Python 3.14.7, numpy 2.5.3 and scipy 1.18.1.

```
pip install numpy scipy
```

## Reproduce

From this folder:

```
python -u run_all.py                 # everything; about 19 min (1125 s, summary.json)
python -u buchstab.py                # under 1 s
python -u checks.py                  # under 1 s
python -u losses.py 0.1129 'R*+G4+L'  # per-region losses; about 30 s (run_all.log)
```

The run used 32 cores. Step times from `summary.json`: roots 372 s, convergence 577 s, variants 91 s, losses 73 s, all other steps under 15 s.
