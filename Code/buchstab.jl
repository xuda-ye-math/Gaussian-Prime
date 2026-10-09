# Buchstab's function omega and its integral Omega(u) = int_1^u omega(t) dt (Section 10 of the paper).
#
# omega(u) = 0 for u < 1, omega(u) = 1/u on [1, 2], and u omega(u) = 1 + Omega(u - 1) for u >= 2.  On each unit
# interval [k, k + 1] the functions are Chebyshev series of degree DEG, built in order of k from this relation:
# omega_k(u) = (1 + Omega_{k-1}(u - 1))/u and Omega_k(u) = Omega(k) + int_k^u omega_k.  On a closed unit interval
# omega is analytic, so the series converge geometrically.  The values of omega, omega' and Omega are tabulated
# on the grid of step 1/NT = 1e-5 and evaluated by cubic Hermite interpolation (Omega with derivative omega, omega
# with derivative omega'(u) = (omega(u - 1) - omega(u))/u).  For u >= UMAX, omega(u) = e^{-gamma} to double
# precision and Omega is continued linearly.

using Printf
using Random
using PolyLog: li2

const DEG = 40
const UMAX = 40
const NT = 100_000
const EXP_MINUS_GAMMA = exp(-Base.MathConstants.eulergamma)

"""A Chebyshev series sum_m c[m + 1] T_m(x) on [lo, hi], x the affine image of u in [-1, 1]."""
struct Cheb
    c::Vector{Float64}
    lo::Float64
    hi::Float64
end

function (s::Cheb)(u)
    x = (2u - s.lo - s.hi) / (s.hi - s.lo)
    b1, b2 = 0.0, 0.0
    for j in length(s.c):-1:2
        b1, b2 = 2x * b1 - b2 + s.c[j], b1
    end
    return x * b1 - b2 + s.c[1]
end

"""Interpolation of f at the DEG + 1 Chebyshev points of the first kind on [lo, hi]."""
function chebinterp(f, lo, hi)
    n = DEG + 1
    t = [pi * (j + 0.5) / n for j in 0:DEG]
    fx = [f(lo + (hi - lo) * (cos(tj) + 1) / 2) for tj in t]
    c = [2 / n * sum(fx[j] * cos(m * t[j]) for j in 1:n) for m in 0:DEG]
    c[1] /= 2
    return Cheb(c, lo, hi)
end

"""The antiderivative of s that vanishes at s.lo."""
function chebinteg(s::Cheb)
    n = length(s.c)
    c = [s.c; 0.0; 0.0]
    C = zeros(n + 1)
    for m in 1:n
        C[m + 1] = ((m == 1 ? 2c[1] : c[m]) - c[m + 2]) / (2m) * (s.hi - s.lo) / 2
    end
    C[1] = -sum(C[m + 1] * (-1)^m for m in 1:n)
    return Cheb(C, s.lo, s.hi)
end

"""The series omega_k and Omega_k on [k, k + 1], k = 1, ..., UMAX - 1."""
function buchstab_series()
    om = Vector{Cheb}(undef, UMAX - 1)
    Om = Vector{Cheb}(undef, UMAX - 1)
    om[1] = chebinterp(u -> 1 / u, 1.0, 2.0)
    Om[1] = chebinteg(om[1])
    for k in 2:UMAX - 1
        prev = Om[k - 1]
        om[k] = chebinterp(u -> (1 + prev(u - 1)) / u, Float64(k), Float64(k + 1))
        I = chebinteg(om[k])
        Om[k] = Cheb([I.c[1] + Om[k - 1](Float64(k)); I.c[2:end]], I.lo, I.hi)
    end
    return om, Om
end

"""Columns k = 1, ..., UMAX - 1: omega, omega' and Omega at k + j/NT, j = 0, ..., NT."""
function buchstab_tables()
    om, Om = buchstab_series()
    W, DW, OM = (Matrix{Float64}(undef, NT + 1, UMAX - 1) for _ in 1:3)
    Threads.@threads for k in 1:UMAX - 1
        for j in 0:NT
            u = k + j / NT
            w = om[k](u)
            W[j + 1, k] = w
            DW[j + 1, k] = k == 1 ? -1 / u^2 : (om[k - 1](u - 1) - w) / u
            OM[j + 1, k] = Om[k](u)
        end
    end
    return W, DW, OM, Om[UMAX - 1](Float64(UMAX))
end

const W_TABLE, DW_TABLE, OM_TABLE, OMEGA_END = buchstab_tables()

@inline function hermite(val, der, u)
    k = unsafe_trunc(Int, u)
    x = (u - k) * NT
    j = min(unsafe_trunc(Int, x), NT - 1)
    s = x - j
    @inbounds y0, y1 = val[j + 1, k], val[j + 2, k]
    @inbounds d0, d1 = der[j + 1, k] / NT, der[j + 2, k] / NT
    s2 = s * s
    s3 = s2 * s
    return (2s3 - 3s2 + 1) * y0 + (s3 - 2s2 + s) * d0 + (3s2 - 2s3) * y1 + (s3 - s2) * d1
end

"""Buchstab's omega(u): 0 for u < 1, e^{-gamma} for u >= UMAX."""
@inline omega(u::Float64) = u < 1.0 ? 0.0 : u >= UMAX ? EXP_MINUS_GAMMA :
                             hermite(W_TABLE, DW_TABLE, min(u, UMAX - 1e-9))

"""Omega(u) = int_1^u omega: 0 for u <= 1, linear with slope e^{-gamma} beyond UMAX."""
@inline Omega(u::Float64) = u <= 1.0 ? 0.0 : u >= UMAX ? OMEGA_END + EXP_MINUS_GAMMA * (u - UMAX) :
                             hermite(OM_TABLE, W_TABLE, min(u, UMAX - 1e-9))

"""Accuracy of omega and Omega against the closed forms on [1, 3], the limit e^{-gamma}, and the relation
u omega(u) = 1 + Omega(u - 1) on [2, UMAX)."""
function buchstab_check()
    u12 = range(1.0, 2.0, length = 20001)
    u23 = range(2.0, 3.0, length = 20001)
    Om23(u) = log(u) + log(u - 1) * log(u) + li2(1 - u) + pi^2 / 12
    ur = 2.0 .+ (UMAX - 1e-6 - 2.0) .* rand(Xoshiro(0), 200_000)
    return [
        "max |omega - 1/u| on [1,2]" => maximum(abs(omega(u) - 1 / u) for u in u12),
        "max |Omega - log u| on [1,2]" => maximum(abs(Omega(u) - log(u)) for u in u12),
        "max |omega - (1+log(u-1))/u| on [2,3]" => maximum(abs(omega(u) - (1 + log(u - 1)) / u) for u in u23),
        "max |Omega - closed form| on [2,3]" => maximum(abs(Omega(u) - Om23(u)) for u in u23),
        "max |u omega(u) - 1 - Omega(u-1)| on [2,40)" => maximum(abs(u * omega(u) - 1 - Omega(u - 1)) for u in ur),
        "omega(10) - e^-gamma" => omega(10.0) - EXP_MINUS_GAMMA,
        "omega(39.9) - e^-gamma" => omega(39.9) - EXP_MINUS_GAMMA,
        "omega(2.5)" => omega(2.5),
        "omega(3)" => omega(3.0),
    ]
end

format_check(name, v) = abs(v) < 1e-3 ? @sprintf("%-46s % .3e", name, v) : @sprintf("%-46s %.15f", name, v)

if abspath(PROGRAM_FILE) == @__FILE__
    for (name, v) in buchstab_check()
        println(format_check(name, v))
    end
end
