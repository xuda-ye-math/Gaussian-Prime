import Mathlib

/-!
# The statement of Theorem 2

A non-zero Gaussian integer is represented by its associate with argument in `[0, π/2)`, that is,
with positive real part and non-negative imaginary part, and arguments are taken modulo `π/2`.
`sectorPrimes θ ψ X` is the set of Gaussian primes counted in Theorem 2: norm in `[X/2, X)` and
argument in `[ψ, ψ + X^(θ - 1/2))` modulo `π/2`.

`SectorTheorem θ` is Theorem 2 at exponent `θ`, in a weaker form: the paper's threshold for `X`
does not depend on `ψ`, while here it may. The weaker form follows from the paper's statement,
and it is all that the proof of Theorem 1 uses.
-/

local notation "ℤ[i]" => GaussianInt

noncomputable section

namespace GaussianPrime

open Real Filter

def sectorPrimes (θ ψ X : ℝ) : Set ℤ[i] :=
  {z | Prime z ∧ 0 < z.re ∧ 0 ≤ z.im ∧ X / 2 ≤ (z.norm : ℝ) ∧ (z.norm : ℝ) < X ∧
    toIcoMod (div_pos pi_pos two_pos) 0 (Complex.arg (z : ℂ) - ψ) < X ^ (θ - 1 / 2)}

def SectorTheorem (θ : ℝ) : Prop :=
  ∃ c > 0, ∀ ψ : ℝ, ∀ᶠ X in atTop,
    c * X ^ (1 / 2 + θ) / Real.log X ≤ ((sectorPrimes θ ψ X).ncard : ℝ)

end GaussianPrime
