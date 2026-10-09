# Lean formalization of the paper (`../Paper/main.tex`)

Lean 4 (v4.33.1) with Mathlib v4.33.1. Build with `lake build`. There is no `sorry`, no `admit`
and no `native_decide`; computations are checked by the kernel. The sources are in
`GaussianPrime/`. Section, table, display and lemma numbers are those of the paper's current build.

## What can be verified

Lean checks that every proof follows from its stated assumptions. It cannot check that a Lean
definition says what the paper means: the definitions `sectorPrimes` and `SectorTheorem`
(`Statements.lean`), and later the regions and `Φ`, are to be compared with the paper by a reader.

| Part of the paper | What Lean verifies | Status |
|---|---|---|
| Theorem 1 from Theorem 2 (Section 1) | the full deduction: a Gaussian prime in the sector at `ψ = 0` gives a prime `p = m² + n²` with `n < p^θ`, hence `{√p} < p^(2θ − 1/2)` | proved: `theorem1` |
| Parameters (3) and Table 4 (Section 2, Appendix D) | every relation for `1/9 < θ < 1/8`, every quantity, every printed value at `θ = 0.1129` | proved |
| Buchstab's function (Sections 3 and 10) | the definition by `uω(u) = 1 + ∫₁^{u-1} ω`, the closed forms on `[1, 3]`, the bounds used | verifiable; layer 2 |
| Decomposition and loss regions (Sections 3 and 7, Table 1) | the regions as sets cut out by linear inequalities, and the facts the paper checks on them: the hypotheses of each range, `E₁₀` inside Range 5, `P ⊆ L`, `R ≠ ∅` | verifiable by linear arithmetic; layer 3 |
| `Φ(0.1129) < 1` (Section 10, Table 2) | a rigorous upper bound by exact rational arithmetic checked by the kernel, in place of the floating-point computation of `../Code` | verifiable; layer 4 |
| Sieve combinatorics (Buchstab's identity (5), the reduction to short intervals of Section 5, the fundamental lemma, Lemma 3, proved by Brun's pure sieve in Appendix B) | elementary but long proofs | verifiable in principle; not yet planned |
| Asymptotic formulas and the lower bound (Sections 4–9) | they rest on the arithmetic information below and on the prime number theorem in progressions modulo `4` | assumed: one axiom `lower_bound`, layer 5 |
| Arithmetic information (Lemmas 20–23 of Appendix A: the sector inequality, the mean-value estimate and the prime number theorem for Hecke characters, the linear estimates of Harman and Lewis) | none: these results are not in Mathlib, and the paper quotes them without proof | not verifiable in this project |

When layer 5 is done, `#print axioms GaussianPrime.theorem1` will list `lower_bound` as the only
non-standard axiom. Theorem 1 is then verified relative to the analytic part of Sections 4–9, with
the numerical claim and every elementary step checked by Lean.

## Proved so far (layer 1)

| Lean | Statement | File |
|---|---|---|
| `SectorTheorem θ` | Theorem 2 at exponent `θ` (a definition, not yet proved) | `Statements.lean` |
| `theorem1` | `SectorTheorem θ₀ → ` Theorem 1: infinitely many primes `p = m² + n²` with `0 < m`, `0 < n`, `n < p^0.1129`, and infinitely many primes with `{√p} < p^(-0.2742)` | `Deduction.lean` |
| `infinite_small_n`, `infinite_fract_sqrt` | the same deduction for every `0 < θ < 1/2`, resp. `0 < θ < 1/4` | `Deduction.lean` |
| `norm_natAbs_prime` | a Gaussian prime `m + in` with `m, n ≠ 0` has a rational prime norm `m² + n²` | `GaussianPrimeNorm.lean` |
| `appD_relations`, `appD_quantities`, `appD_quantities_pos`, `appD_values` | Table 4 (Appendix D): every relation for `1/9 < θ < 1/8`, every first-column quantity, and every printed value at `θ = 0.1129` to `5 · 10⁻⁸` | `Params.lean` |

`#print axioms GaussianPrime.theorem1` gives `propext`, `Classical.choice`, `Quot.sound`.

`SectorTheorem θ` follows the paper's convention: a Gaussian prime is represented by its associate
with positive real part and non-negative imaginary part, and the argument condition is taken
modulo `π/2` (`toIcoMod`). The threshold for `X` may depend on `ψ`, which is the weaker reading of
Theorem 2 and all that the deduction uses.

## Plan for the next layers

The proof of Theorem 2 has two parts of very different nature.

1. **The analytic part (Sections 4–9).** From the arithmetic information quoted from Harman and
   Lewis, it proves the lower bound: if `Φ(θ) < 1`, then Theorem 2 holds at `θ`. Its inputs
   (Hecke `L`-functions of `ℚ(i)`, the sector inequality, the mean-value estimate, the prime number
   theorem with Hecke characters) are not in Mathlib. This part enters first as one explicitly
   stated axiom in `Axioms.lean`, under the standing assumption `1/9 < θ < 1/8` of Sections 3 and 7:
   `lower_bound : ∀ θ, 1/9 < θ → θ < 1/8 → Φ θ < 1 → SectorTheorem θ`, so that
   `#print axioms` shows exactly what is assumed. Later layers can replace it, lemma by lemma,
   by proofs from smaller axioms (the arithmetic information alone).
2. **The numerical part (Section 10): `Φ(0.1129) < 1`.** This is finite and is to be proved in
   full, by the kernel, at the single value `θ = 0.1129` (no root finding).

Layers, each building on a product that already builds:

- **Layer 2 — Buchstab's function.** Define `ω` by `uω(u) = 1 + ∫₁^{u-1} ω` (`ω = 0` on
  `u < 1`) and `Ω(u) = ∫₁^u ω`; prove `ω = 1/u` on `[1, 2]`, `ω = (1 + log(u - 1))/u` on `[2, 3]`,
  `ω(u) = (1 + Ω(u - 1))/u` on `[3, 4]`, `0 ≤ ω ≤ 1`, and a bound `ω ≤ M` on `[2, ∞)` for the
  four-dimensional regions.
- **Layer 3 — the regions and the loss.** Define `D₁, D₂, D₄, D₆, D₈`, `R` (two-dimensional) and
  `E₃, E_R, E₅, E₇, E₉, E₁₀, E₁₁`, `G`, `I`, `L` (four-dimensional) exactly as in Sections 3
  and 7, and `Φ θ` as display (18), with the integrals as lower Lebesgue integrals of non-negative
  integrands (no integrability side conditions).
- **Layer 4 — the certified bound `Φ θ₀ < 1`.** The margin is `1 − 0.9990033 ≈ 9.97 · 10⁻⁴`.
  - Two-dimensional losses (`99.4 %` of `Φ`): the substitution `u = (1 − α − β)/β` gives
    `dβ/β² = du/(1 − α)`, so each loss is the one-dimensional integral
    `∫ (Ω(u_hi(α)) − Ω(u_lo(α))) / (α(1 − α)) dα`. On these regions `u ∈ [1, 3.56]`, where `Ω` is
    built from `log` and `∫ log(t − 1)/t`. Rational enclosures of `log` come from Mathlib's series
    bounds; the outer integral is bounded by monotone upper sums.
  - Four-dimensional losses (`0.0054` in total, so about `15 %` relative accuracy suffices): the
    `δ`-integral is exact through `Ω` on a finite union of intervals; the outer `(α, β, γ)`-integral
    is bounded cell by cell. Outside `G ∪ L`, random sampling (`8 · 10⁶` points per region) finds `u` up
    to about `4.43` on `E₁₁` and `3.43` on `E₁₀`, and at most `2.94` on the other regions.
  - Each finite computation is a function on `ℚ` with a soundness theorem, evaluated by the kernel.
- **Layer 5 — assembly.** `theorem1` from `lower_bound` and layer 4, with `#print axioms` listing
  `lower_bound` as the only non-standard axiom.
