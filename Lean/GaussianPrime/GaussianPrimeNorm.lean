import Mathlib

/-!
# Norms of Gaussian primes

A Gaussian prime `m + in` with `m ≠ 0` and `n ≠ 0` has a rational prime `m² + n²` as its norm.
This is the step "Its norm `p = m² + n²` is a rational prime" in the proof of Theorem 1 from
Theorem 2.
-/

local notation "ℤ[i]" => GaussianInt

namespace GaussianPrime

theorem norm_natAbs_prime {z : ℤ[i]} (hz : Prime z) (hre : z.re ≠ 0) (him : z.im ≠ 0) :
    z.norm.natAbs.Prime := by
  set N := z.norm.natAbs with hNdef
  have hz0 : z ≠ 0 := hz.ne_zero
  have hNz : (N : ℤ) = z.norm := GaussianInt.abs_natCast_norm z
  have hN0 : N ≠ 0 := by
    intro h
    apply hz0
    rw [← GaussianInt.norm_eq_zero, ← hNz, h, Nat.cast_zero]
  have hN1 : N ≠ 1 := fun h => hz.not_isUnit (Zsqrtd.norm_eq_one_iff.1 h)
  -- `z` divides `N = z * conj z`
  have hzN : z ∣ (N : ℤ[i]) := by
    refine ⟨star z, ?_⟩
    rw [hNdef, GaussianInt.natCast_natAbs_norm, Zsqrtd.norm_eq_mul_conj]
  -- `z` divides a prime factor `q` of `N`
  have hprod : ((N.primeFactorsList.map (fun q : ℕ => (q : ℤ[i]))).prod) = (N : ℤ[i]) := by
    rw [← Nat.cast_list_prod, Nat.prod_primeFactorsList hN0]
  obtain ⟨a, ha, hza⟩ := (hz.dvd_prod_iff).1 (hprod ▸ hzN)
  obtain ⟨q, hq, rfl⟩ := List.mem_map.1 ha
  have hqp : q.Prime := Nat.prime_of_mem_primeFactorsList hq
  obtain ⟨w, hw⟩ := hza
  have hnorm : (q : ℤ) * q = z.norm * w.norm := by
    have := congrArg Zsqrtd.norm hw
    rwa [Zsqrtd.norm_natCast, Zsqrtd.norm_mul] at this
  -- `N` divides `q²`
  have hdvd : N ∣ q ^ 2 := by
    refine ⟨w.norm.natAbs, ?_⟩
    have := congrArg Int.natAbs hnorm
    rw [Int.natAbs_mul, Int.natAbs_mul, Int.natAbs_natCast] at this
    rw [sq, this]
  obtain ⟨k, hk, hNk⟩ := (Nat.dvd_prime_pow hqp).1 hdvd
  interval_cases k
  · exact absurd (by simpa using hNk) hN1
  · simpa [hNk] using hqp
  · -- `N = q²` forces `w` to be a unit, so `z` is an associate of `q`
    exfalso
    have hq0 : (q : ℤ) ≠ 0 := by exact_mod_cast hqp.ne_zero
    have hw1 : w.norm = 1 := by
      have h2 : z.norm = (q : ℤ) * q := by rw [← hNz, hNk]; push_cast; ring
      rw [h2] at hnorm
      have := mul_self_ne_zero.2 hq0
      exact (mul_right_eq_self₀.1 hnorm.symm).resolve_right this
    have hz' : z = (q : ℤ[i]) * star w := by
      have h := congrArg (· * star w) hw
      rw [mul_assoc, ← Zsqrtd.norm_eq_mul_conj, hw1, Int.cast_one, mul_one] at h
      exact h.symm
    have hre' : z.re = q * w.re := by rw [hz']; simp
    have him' : z.im = -(q * w.im) := by rw [hz']; simp
    have hsq : w.re * w.re + w.im * w.im = 1 := by
      have := hw1; rw [Zsqrtd.norm_def] at this; linarith
    rcases eq_or_ne w.re 0 with h0 | h0
    · exact hre (by rw [hre', h0, mul_zero])
    rcases eq_or_ne w.im 0 with h1 | h1
    · exact him (by rw [him', h1, mul_zero, neg_zero])
    have a1 : 0 < w.re * w.re := mul_self_pos.2 h0
    have a2 : 0 < w.im * w.im := mul_self_pos.2 h1
    omega

end GaussianPrime
