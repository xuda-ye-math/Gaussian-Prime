# Parameters, discarded regions and the sets G, I, L, P of the paper (Paper/main.tex).
#
# Labels follow the paper: the two-dimensional regions D1, D2, D4, D6, D8 and the four-dimensional regions
# E3, ER, E5, E7, E9, E10, E11 of Section 7 (Table 1 and display (37)), the region R (13), the sets G, I, L (14)
# and the polytope P (15).  "HL" is G. Harman and P. Lewis, Gaussian primes in narrow sectors, Mathematika 48
# (2001) 119-135 (journal page numbers); the index k of D_k and E_k is that of their loss phi_k.  All exponents
# are exponents of norms: |p|^2 = X^a, |q|^2 = X^b, |r|^2 = X^g, |s|^2 = X^d, with a > b > g > d.
#
# A two-dimensional region is {a0 <= a <= a1, blo(a) <= b <= bhi(a)}; its loss is the integral of
# omega((1 - a - b)/b)/(a b^2).  A four-dimensional region adds glo(a, b) <= g <= ghi(a, b) and
# dlo(a, b, g) <= d <= dhi(a, b, g); its loss is the integral of omega((1 - a - b - g - d)/d)/(a b g d^2), where
# omega(u) = 0 for u < 1 caps d at (1 - a - b - g)/2.  `kinks` lists the values of a at which a limit of the
# region has a kink or the region opens or closes; the first method splits there.
#
# The sets G, L and P are sets of 4-tuples whose discarded sum gets an asymptotic formula.  They are given in two
# forms: pointwise, as defined in the paper (the second method of Section 10), and as a union of d-intervals for
# fixed (a, b, g) (the first method).

"""The parameters (3) of the paper, and chi = (nu + 2 theta)/2 of Section 7, in the arithmetic of theta."""
function params(theta)
    return (theta = theta, sigma = (1 - 2theta) / 2, tau = (1 + 2theta) / 2, rho = (1 + 22theta) / 12,
            zeta = (1 + 4theta) / 3, nu = (10theta + 7) / 12, chi = (34theta + 7) / 24)
end

struct Region2{F1, F2}
    name::String
    a0::Float64
    a1::Float64
    blo::F1
    bhi::F2
    kinks::Vector{Float64}
end

struct Region4{F1, F2, F3, F4, F5, F6}
    name::String
    a0::Float64
    a1::Float64
    blo::F1
    bhi::F2
    glo::F3
    ghi::F4
    dlo::F5
    dhi::F6
    kinks::Vector{Float64}
end

Region2(name, a0, a1, blo, bhi) = Region2(name, Float64(a0), Float64(a1), blo, bhi, Float64[])
Region4(name, a0, a1, blo, bhi, glo, ghi, dlo, dhi) =
    Region4(name, Float64(a0), Float64(a1), blo, bhi, glo, ghi, dlo, dhi, Float64[])

"""The two-dimensional regions D1, D2, D4, D6, D8 of Section 7.  With rstar, D2 excludes R, as in the paper;
without it (the decomposition of HL), the whole sum is discarded on R too, and the region is "D2+R"."""
function regions_2d(theta::Float64, rstar::Bool)
    p = params(theta)
    s, t, r, z, n, c = p.sigma, p.tau, p.rho, p.zeta, p.nu, p.chi
    d2lo(a) = max(z - a, (t - a) / 2)
    d2 = rstar ? Region2("D2", 1 / 3, s, a -> max(d2lo(a), r / 2), a -> (1 - a) / 2, [2z - t, z - r / 2, t - r]) :
                 Region2("D2+R", 1 / 3, s, d2lo, a -> (1 - a) / 2, [2z - t])
    return Any[
        # D1 (phi1 of HL, p. 132): the pieces above and below the set where Lemma 9 applies, and a >= chi.
        Region2("D1", s, c, a -> n - a, a -> (1 - a) / 2),
        Region2("D1", 0.25 + 1.5theta, c, a -> t - a, a -> a - 2theta),
        Region2("D1", c, 0.5, a -> t - a, a -> (1 - a) / 2),
        d2,
        # D4 (phi4 of HL, p. 132).
        Region2("D4", r, 1 / 3, a -> z - a, a -> a),
        # D6 (phi6 of HL, p. 132).
        Region2("D6", 0.5 - 2theta, r, a -> (n - a) / 2, a -> a),
        # D8 (phi8 of HL, p. 133).
        Region2("D8", n / 3, s - theta, a -> (n - a) / 2, a -> a),
    ]
end

"""The four-dimensional regions E3, E5, E7, E9, E10, E11 of display (37), and ER if rstar.  E10 and E11 are
integrated over their full a-ranges for every theta, which for theta < 3/26 includes the interval
[tau/4, sigma - 2 theta)."""
function regions_4d(theta::Float64, rstar::Bool)
    p = params(theta)
    s, t, r, z, n = p.sigma, p.tau, p.rho, p.zeta, p.nu
    # E3 and E9 carry HL's cut g <= a + b - 2 theta (Lemma 9 applied to pq and r).  It is implied by the
    # d-limits tau - a - b <= d <= (1 - a - b - g)/2: these give tau - a - b <= (1 - a - b - g)/2, i.e.
    # g <= 1 - 2 tau + a + b = a + b - 2 theta.  Where the cut binds, the d-range is empty, so the cut changes
    # no loss (Appendix E).
    regs = Any[
        # E3: Lemma 11, identity (35) of HL; phi3 of HL (p. 132).  The b-range closes at a = 2 zeta - tau.
        Region4("E3", 1 / 3, s, a -> (t - a) / 2, a -> z - a, (a, b) -> t - a - b, (a, b) -> min(b, a + b - 2theta),
                (a, b, g) -> t - a - b, (a, b, g) -> g, [2z - t]),
        # E5: Lemma 11, identity (35) of HL; phi5 of HL (p. 132).
        Region4("E5", r, 1 / 3, a -> (t - a) / 2, a -> z - a, (a, b) -> t - a - b, (a, b) -> b,
                (a, b, g) -> t - a - b, (a, b, g) -> g),
        # E7: Lemma 11, identity (35) of HL; phi7 of HL (p. 133).  HL's exclusion "d > 1 - a - g - sigma
        # if a + g >= sigma" (HL p. 133) is kept in the region; it lies inside G (S = {p, r}), so it changes
        # nothing when G is used, and it reproduces HL's region otherwise (Appendix E).
        Region4("E7", 0.5 - 2theta, r, a -> (t - a) / 2, a -> (n - a) / 2, (a, b) -> t - a - b, (a, b) -> b,
                (a, b, g) -> a + g >= s ? max(t - a - b, 1 - a - g - s) : t - a - b, (a, b, g) -> g),
        # E9: Lemma 11, identity (35) of HL, with pq >= W; phi9 of HL (p. 133).
        Region4("E9", t / 3, s - theta, a -> (t - a) / 2, a -> min((n - a) / 2, a), (a, b) -> t - a - b,
                (a, b) -> min(b, b + a - 2theta), (a, b, g) -> t - a - b, (a, b, g) -> g, [n / 3, 2theta]),
        # E10: Lemma 12 with j = 0, identity (36) of HL, with pq < W; phi10 of HL (p. 133), with the upper
        # b-limit min(a, sigma - a) (Remark 27).
        Region4("E10", t / 4, s - theta, a -> (t - a) / 3, a -> min(a, s - a), (a, b) -> (t - a - b) / 2,
                (a, b) -> b, (a, b, g) -> t - a - b - g, (a, b, g) -> g, [s / 2]),
        # E11: Lemma 12 with j = 1; phi11 of HL (p. 134).
        Region4("E11", t / 5, s - 2theta, a -> (t - a) / 4, a -> min(a, s - theta - a), (a, b) -> (t - a - b) / 3,
                (a, b) -> b, (a, b, g) -> (t - a - b - g) / 2, (a, b, g) -> g, [(s - theta) / 2]),
    ]
    if rstar
        # ER: the third sum of identity (35) of HL on R = {1/3 <= a < sigma, max(zeta - a, (tau - a)/2) <= b <= rho/2}.
        push!(regs, Region4("ER", 1 / 3, s, a -> max(z - a, (t - a) / 2), a -> r / 2, (a, b) -> t - a - b,
                            (a, b) -> b, (a, b, g) -> t - a - b, (a, b, g) -> g, [2z - t, z - r / 2, t - r]))
    end
    return regs
end

# ----------------------------------------------------------------------------------------------------------
# The sets G, I, L, P.  Modes of the rescued set: "none"; "G"; "G+P" (G union P); "G+L" (G union L).
# A subset of {p, q, r, s} is a bit mask: 1 = p, 2 = q, 4 = r, 8 = s.

const MODES = ("none", "G", "G+P", "G+L")

@inline masked_sum(mask, a, b, g, d) =
    (mask & 1 != 0 ? a : 0.0) + (mask & 2 != 0 ? b : 0.0) + (mask & 4 != 0 ? g : 0.0) + (mask & 8 != 0 ? d : 0.0)

"""G of (14): the union over non-empty S of {p, q, r} of {sigma <= sum_S <= tau - d}."""
@inline function g_set(a, b, g, d, p)
    for S in 1:7
        sS = masked_sum(S, a, b, g, d)
        (sS >= p.sigma) & (sS <= p.tau - d) && return true
    end
    return false
end

"""I of (14): {a + b + g + d <= tau}, union, over the splits U | V of {p, q, r, s}, of
{sum_U <= zeta, sum_V <= rho, a + b + g + d <= nu}."""
@inline function i_set(a, b, g, d, p)
    s4 = a + b + g + d
    s4 < p.tau && return true
    for U in 1:14
        u = masked_sum(U, a, b, g, d)
        (u <= p.zeta) & (s4 - u <= p.rho) & (s4 <= p.nu) && return true
    end
    return false
end

"""L of (14): I intersected with the union over non-empty S of {p, q, r, s} of {sigma <= sum_S <= tau - d}."""
@inline function l_set(a, b, g, d, p)
    i_set(a, b, g, d, p) || return false
    for S in 1:15
        sS = masked_sum(S, a, b, g, d)
        (sS >= p.sigma) & (sS <= p.tau - d) && return true
    end
    return false
end

"""P of (15): a + b <= sigma, g + d <= rho, a + b + d >= sigma, a + b + 2d <= tau."""
@inline p_set(a, b, g, d, p) =
    (a + b <= p.sigma) & (g + d <= p.rho) & (a + b + d >= p.sigma) & (a + b + 2d <= p.tau)

"""Whether (a, b, g, d) is rescued in mode (pointwise; the second method).  Mode codes: 1 none, 2 G, 3 G+P,
4 G+L (the index in MODES)."""
@inline function rescued(mode::Int, a, b, g, d, p)
    mode == 1 && return false
    g_set(a, b, g, d, p) && return true
    mode == 3 && return p_set(a, b, g, d, p)
    mode == 4 && return l_set(a, b, g, d, p)
    return false
end

# The sets as unions of d-intervals for fixed (a, b, g), for the first method.
#   G:  d in [0, tau - sum_S] for each non-empty S of {p, q, r} with sum_S >= sigma.
#   I:  a down-set [0, dI]: a + b + g + d <= tau gives d < tau - (a + b + g); a split with s on the U-side
#       gives, if sum_V <= rho, d <= min(zeta - u', nu - (a + b + g)); with s on the V-side, if sum_U <= zeta,
#       d <= min(rho - v', nu - (a + b + g)); u', v' are the sums of the other primes on the side of s.
#   L:  I intersected with the G intervals (S without s) and with [sigma - S', (tau - S')/2] for S = S' + {s},
#       S' a subset of {p, q, r} (possibly empty).  So G union L = G union (I intersected with the latter).
#   P:  d in [sigma - a - b, min(rho - g, (tau - a - b)/2)] if a + b <= sigma.
# The d-sets jump where a subset sum of (a, b, g) crosses sigma (G, P, E7's exclusion), rho or zeta (the two
# sides of the splits of I); these are the jump levels, and the first method splits every axis there.

jump_levels(p) = (p.sigma, p.rho, p.zeta)

"""dI with I <=> d <= dI (up to the boundary), for fixed (a, b, g)."""
@inline function i_threshold(a, b, g, p)
    s3 = a + b + g
    dI = p.tau - s3
    for U in 1:14
        u = masked_sum(U & 7, a, b, g, 0.0)
        v = s3 - u
        if U & 8 != 0
            dI = max(dI, v <= p.rho ? min(p.zeta - u, p.nu - s3) : -Inf)
        else
            dI = max(dI, u <= p.zeta ? min(p.rho - v, p.nu - s3) : -Inf)
        end
    end
    return dI
end

"""The d-intervals rescued in mode for fixed (a, b, g), written into lo, hi; returns their number."""
@inline function rescued_intervals!(lo, hi, mode::Int, a, b, g, p)
    mode == 1 && return 0
    n = 0
    for S in 1:7
        sS = masked_sum(S, a, b, g, 0.0)
        n += 1
        @inbounds lo[n], hi[n] = 0.0, sS >= p.sigma ? p.tau - sS : -Inf
    end
    if mode == 3
        n += 1
        @inbounds lo[n], hi[n] = a + b <= p.sigma ? p.sigma - a - b : Inf, min(p.rho - g, (p.tau - a - b) / 2)
    elseif mode == 4
        dI = i_threshold(a, b, g, p)
        for S in 0:7
            sS = masked_sum(S, a, b, g, 0.0)
            n += 1
            @inbounds lo[n], hi[n] = p.sigma - sS, min((p.tau - sS) / 2, dI)
        end
    end
    return n
end
