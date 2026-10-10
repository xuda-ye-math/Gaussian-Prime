import Mathlib
import GaussianPrime.Params
import GaussianPrime.Statements
import GaussianPrime.GaussianPrimeNorm

/-!
# Theorem 1 from Theorem 2

The proof of Theorem 1 in Section 1: apply Theorem 2 with `ψ = 0`. Fewer than `√X + 1` of the
counted Gaussian primes have imaginary part `0`, so for large `X` the sector contains a Gaussian
prime `m + in` with `m, n > 0`; its norm `p = m² + n²` is a rational prime with `n < p^θ`. Then
`0 < √p - m = n² / (√p + m) < p^(2θ - 1/2)` gives the statement on `{√p}`.
-/

local notation "ℤ[i]" => GaussianInt

noncomputable section

namespace GaussianPrime

open Real Filter Asymptotics

/-- If every counted Gaussian prime has imaginary part `0`, there are at most `√X + 1` of them. -/
lemma ncard_le_of_im_eq_zero {θ X : ℝ}
    (h : ∀ z ∈ sectorPrimes θ 0 X, z.im = 0) :
    ((sectorPrimes θ 0 X).ncard : ℝ) ≤ √X + 1 := by
  classical
  set F : Finset ℤ[i] := (Finset.range (⌊√X⌋₊ + 1)).image (fun a : ℕ => (⟨(a : ℤ), 0⟩ : ℤ[i]))
  have hsub : sectorPrimes θ 0 X ⊆ (F : Set ℤ[i]) := by
    intro z hz
    have him := h z hz
    obtain ⟨-, hre, -, -, hlt, -⟩ := hz
    have hnorm : (z.norm : ℝ) = (z.re : ℝ) ^ 2 := by
      rw [Zsqrtd.norm_def, him]; push_cast; ring
    simp only [F, Finset.coe_image, Finset.coe_range, Set.mem_image, Set.mem_Iio]
    refine ⟨z.re.toNat, ?_, ?_⟩
    · have hre' : ((z.re.toNat : ℕ) : ℝ) = (z.re : ℝ) := by
        exact_mod_cast Int.toNat_of_nonneg hre.le
      have hlt' : (z.re : ℝ) < √X := by
        rw [Real.lt_sqrt (by exact_mod_cast hre.le), ← hnorm]; exact hlt
      have : (z.re.toNat : ℝ) ≤ √X := by rw [hre']; exact hlt'.le
      exact Nat.lt_succ_of_le (Nat.le_floor this)
    · ext
      · simp [Int.toNat_of_nonneg hre.le]
      · simp [him]
  calc ((sectorPrimes θ 0 X).ncard : ℝ) ≤ (F.card : ℝ) := by
        exact_mod_cast (Set.ncard_le_ncard hsub (Finset.finite_toSet F)).trans
          (Set.ncard_coe_finset F).le
    _ ≤ ((Finset.range (⌊√X⌋₊ + 1)).card : ℝ) := by exact_mod_cast Finset.card_image_le
    _ = (⌊√X⌋₊ : ℝ) + 1 := by simp
    _ ≤ √X + 1 := by gcongr; exact Nat.floor_le (Real.sqrt_nonneg X)

/-- The lower bound of Theorem 2 eventually exceeds `√X + 1`. -/
lemma eventually_sqrt_lt {θ c : ℝ} (hθ : 0 < θ) (hc : 0 < c) :
    ∀ᶠ X in atTop, √X + 1 < c * X ^ (1 / 2 + θ) / Real.log X := by
  have h1 := (isLittleO_log_rpow_atTop hθ).def (by positivity : (0 : ℝ) < c / 4)
  filter_upwards [h1, eventually_gt_atTop (1 : ℝ)] with X hX hX1
  have hlog : 0 < Real.log X := Real.log_pos hX1
  have hX0 : 0 < X := by linarith
  rw [Real.norm_of_nonneg hlog.le, Real.norm_of_nonneg (by positivity)] at hX
  have hsq : X ^ (1 / 2 + θ) = √X * X ^ θ := by
    rw [Real.rpow_add hX0, Real.sqrt_eq_rpow]
  have hs : 1 < √X := by
    rw [show (1 : ℝ) = √1 by simp]; exact Real.sqrt_lt_sqrt zero_le_one hX1
  have hP : 0 < √X * X ^ θ := by positivity
  rw [hsq, lt_div_iff₀ hlog]
  calc (√X + 1) * Real.log X < (2 * √X) * Real.log X := by
        exact mul_lt_mul_of_pos_right (by linarith) hlog
    _ ≤ (2 * √X) * (c / 4 * X ^ θ) := by gcongr
    _ < c * (√X * X ^ θ) := by nlinarith [mul_pos hc hP]

/-- For every bound `B`, Theorem 2 with `ψ = 0` gives an `X > 2B` and a counted Gaussian prime
with non-zero imaginary part. -/
lemma exists_far {θ : ℝ} (hθ : 0 < θ) (h : SectorTheorem θ) (B : ℝ) :
    ∃ X : ℝ, B < X / 2 ∧ ∃ z ∈ sectorPrimes θ 0 X, z.im ≠ 0 := by
  obtain ⟨c, hc, hX⟩ := h
  obtain ⟨X, hX1, hX2, hX3⟩ :=
    ((hX 0).and ((eventually_sqrt_lt hθ hc).and (eventually_gt_atTop (max (2 * B) 1)))).exists
  have hB := le_max_left (2 * B) 1
  have h1 := le_max_right (2 * B) 1
  refine ⟨X, by linarith, ?_⟩
  by_contra hcon
  push Not at hcon
  have := ncard_le_of_im_eq_zero hcon
  linarith

/-- A counted Gaussian prime `m + in` with `n ≠ 0` gives a rational prime `p = m² + n²` with
`X/2 ≤ p < X` and `n < p^θ`. -/
lemma prime_of_mem {θ X : ℝ} (hθ : θ < 1 / 2) {z : ℤ[i]} (hz : z ∈ sectorPrimes θ 0 X)
    (him : z.im ≠ 0) :
    ∃ m n : ℕ, 0 < m ∧ 0 < n ∧ (m ^ 2 + n ^ 2).Prime ∧ X / 2 ≤ ((m ^ 2 + n ^ 2 : ℕ) : ℝ) ∧
      (n : ℝ) < ((m ^ 2 + n ^ 2 : ℕ) : ℝ) ^ θ := by
  obtain ⟨hp, hre, him0, hlo, hhi, harg⟩ := hz
  have him' : 0 < z.im := lt_of_le_of_ne him0 (Ne.symm him)
  set m := z.re.toNat with hm
  set n := z.im.toNat with hn
  have hmz : (m : ℤ) = z.re := Int.toNat_of_nonneg hre.le
  have hnz : (n : ℤ) = z.im := Int.toNat_of_nonneg him0
  have hnormZ : z.norm = (m : ℤ) ^ 2 + (n : ℤ) ^ 2 := by
    rw [Zsqrtd.norm_def, hmz, hnz]; ring
  have hnatAbs : z.norm.natAbs = m ^ 2 + n ^ 2 := by
    rw [hnormZ]; exact_mod_cast Int.natAbs_natCast (m ^ 2 + n ^ 2)
  have hnormR : (z.norm : ℝ) = ((m ^ 2 + n ^ 2 : ℕ) : ℝ) := by
    rw [hnormZ]; push_cast; ring
  refine ⟨m, n, by omega, by omega, ?_, by rw [← hnormR]; exact hlo, ?_⟩
  · rw [← hnatAbs]; exact norm_natAbs_prime hp hre.ne' him
  -- the analytic part: `n = ‖z‖ sin(arg z) ≤ ‖z‖ arg z < √p · X^(θ - 1/2) < p^θ`
  set w : ℂ := (z : ℂ)
  set P : ℝ := ((m ^ 2 + n ^ 2 : ℕ) : ℝ) with hP
  have hwre : w.re = (z.re : ℝ) := (GaussianInt.intCast_re z).symm
  have hwim : w.im = (z.im : ℝ) := (GaussianInt.intCast_im z).symm
  have hPpos : 0 < P := by rw [← hnormR]; exact_mod_cast (GaussianInt.norm_pos.2 hp.ne_zero)
  have hnormw : ‖w‖ = √P := by
    rw [Complex.norm_def, ← GaussianInt.intCast_real_norm, hnormR]
  have harg0 : 0 ≤ Complex.arg w := Complex.arg_nonneg_iff.2 (by rw [hwim]; exact_mod_cast him0)
  have harg1 : Complex.arg w < π / 2 :=
    Complex.arg_lt_pi_div_two_iff.2 (Or.inl (by rw [hwre]; exact_mod_cast hre))
  have hargX : Complex.arg w < X ^ (θ - 1 / 2) := by
    have h := harg
    rw [sub_zero, (toIcoMod_eq_self _).2 ⟨harg0, by simpa using harg1⟩] at h
    exact h
  have hw0 : ‖w‖ ≠ 0 := by rw [hnormw]; exact (Real.sqrt_pos.2 hPpos).ne'
  have hsin : (n : ℝ) = ‖w‖ * Real.sin (Complex.arg w) := by
    rw [Complex.sin_arg, mul_div_cancel₀ _ hw0, hwim, ← hnz]; norm_cast
  have hPX : P < X := by rw [← hnormR]; exact hhi
  have hrpow : X ^ (θ - 1 / 2) < P ^ (θ - 1 / 2) :=
    Real.rpow_lt_rpow_of_neg hPpos hPX (by linarith)
  have hsqrtP : 0 < √P := Real.sqrt_pos.2 hPpos
  calc (n : ℝ) = ‖w‖ * Real.sin (Complex.arg w) := hsin
    _ ≤ ‖w‖ * Complex.arg w := by gcongr; exact Real.sin_le harg0
    _ < √P * X ^ (θ - 1 / 2) := by rw [hnormw]; gcongr
    _ < √P * P ^ (θ - 1 / 2) := by gcongr
    _ = P ^ θ := by
        rw [Real.sqrt_eq_rpow, ← Real.rpow_add hPpos]; ring_nf

/-- Theorem 1, first statement, from Theorem 2 at any `0 < θ < 1/2`. -/
theorem infinite_small_n {θ : ℝ} (hθ0 : 0 < θ) (hθ : θ < 1 / 2) (h : SectorTheorem θ) :
    {p : ℕ | p.Prime ∧ ∃ m n : ℕ, 0 < m ∧ 0 < n ∧ p = m ^ 2 + n ^ 2 ∧
      (n : ℝ) < (p : ℝ) ^ θ}.Infinite := by
  refine Set.infinite_of_forall_exists_gt fun a => ?_
  obtain ⟨X, hX, z, hz, him⟩ := exists_far hθ0 h a
  obtain ⟨m, n, hm, hn, hp, hlo, hlt⟩ := prime_of_mem hθ hz him
  refine ⟨m ^ 2 + n ^ 2, ⟨hp, m, n, hm, hn, rfl, hlt⟩, ?_⟩
  have : (a : ℝ) < ((m ^ 2 + n ^ 2 : ℕ) : ℝ) := by linarith
  exact_mod_cast this

/-- The step from `n < p^θ` to `{√p} < p^(2θ - 1/2)`, displays (1) and (2) of the paper. -/
lemma fract_sqrt_lt {θ : ℝ} (hθ : θ < 1 / 4) {p m n : ℕ} (hp : p.Prime) (hn : 0 < n)
    (hpmn : p = m ^ 2 + n ^ 2) (hlt : (n : ℝ) < (p : ℝ) ^ θ) :
    Int.fract (√(p : ℝ)) < (p : ℝ) ^ (2 * θ - 1 / 2) := by
  have hp1 : (1 : ℝ) < p := by exact_mod_cast hp.one_lt
  have hp0 : (0 : ℝ) < p := by linarith
  set s := √(p : ℝ) with hs
  have hs0 : 0 < s := Real.sqrt_pos.2 hp0
  have hss : s ^ 2 = (p : ℝ) := Real.sq_sqrt hp0.le
  have hpR : (p : ℝ) = (m : ℝ) ^ 2 + (n : ℝ) ^ 2 := by rw [hpmn]; push_cast; ring
  have hn0 : (0 : ℝ) < n := by exact_mod_cast hn
  have hm0 : (0 : ℝ) ≤ m := Nat.cast_nonneg m
  have hms : (m : ℝ) < s := by
    by_contra hcon
    push Not at hcon
    have : s ^ 2 ≤ (m : ℝ) ^ 2 := by gcongr
    nlinarith
  -- `s - m = n² / (s + m) ≤ n² / s < p^(2θ) / √p = p^(2θ - 1/2)`
  have hdiff : s - m = (n : ℝ) ^ 2 / (s + m) := by
    rw [eq_div_iff (by positivity)]; nlinarith
  have hn2 : (n : ℝ) ^ 2 < (p : ℝ) ^ (2 * θ) := by
    rw [show 2 * θ = θ * 2 by ring, Real.rpow_mul hp0.le, Real.rpow_two]
    exact pow_lt_pow_left₀ hlt (by positivity) two_ne_zero
  have hbound : s - m < (p : ℝ) ^ (2 * θ - 1 / 2) := by
    rw [hdiff, Real.rpow_sub hp0, ← Real.sqrt_eq_rpow]
    calc (n : ℝ) ^ 2 / (s + m) ≤ (n : ℝ) ^ 2 / s := by
          gcongr; linarith
      _ < (p : ℝ) ^ (2 * θ) / s := by gcongr
  have hlt1 : (p : ℝ) ^ (2 * θ - 1 / 2) < 1 :=
    Real.rpow_lt_one_of_one_lt_of_neg hp1 (by linarith)
  have hfloor : ⌊s⌋ = (m : ℤ) := by
    rw [Int.floor_eq_iff]; push_cast; constructor <;> linarith
  rw [Int.fract, hfloor]; push_cast; exact hbound

/-- Theorem 1, second statement, from Theorem 2 at any `0 < θ < 1/4`. -/
theorem infinite_fract_sqrt {θ : ℝ} (hθ0 : 0 < θ) (hθ : θ < 1 / 4) (h : SectorTheorem θ) :
    {p : ℕ | p.Prime ∧ Int.fract (√(p : ℝ)) < (p : ℝ) ^ (2 * θ - 1 / 2)}.Infinite := by
  refine (infinite_small_n hθ0 (by linarith) h).mono ?_
  rintro p ⟨hp, m, n, -, hn, hpmn, hlt⟩
  exact ⟨hp, fract_sqrt_lt hθ hp hn hpmn hlt⟩

/-- **Theorem 1** of the paper, from Theorem 2 at `θ = 0.1129`. -/
theorem theorem1 (h : SectorTheorem θ₀) :
    {p : ℕ | p.Prime ∧ ∃ m n : ℕ, 0 < m ∧ 0 < n ∧ p = m ^ 2 + n ^ 2 ∧
      (n : ℝ) < (p : ℝ) ^ (0.1129 : ℝ)}.Infinite ∧
    {p : ℕ | p.Prime ∧ Int.fract (√(p : ℝ)) < (p : ℝ) ^ (-0.2742 : ℝ)}.Infinite := by
  have e1 : (0.1129 : ℝ) = θ₀ := by norm_num [θ₀]
  have e2 : (-0.2742 : ℝ) = 2 * θ₀ - 1 / 2 := by norm_num [θ₀]
  have h0 : 0 < θ₀ := by norm_num [θ₀]
  rw [e1, e2]
  exact ⟨infinite_small_n h0 (by norm_num [θ₀]) h, infinite_fract_sqrt h0 (by norm_num [θ₀]) h⟩

end GaussianPrime
