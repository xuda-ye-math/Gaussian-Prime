# The losses (16) of the paper, by the two methods of Section 10, and the theta at which the total loss is 1.
#
# Decompositions (the rows of Table 3; each is a list of regions from regions.jl and a rescued set):
#     "HL"        the decomposition of Harman and Lewis, with the regions derived from the lemmas (Remark 27):
#                 D1, D2+R, D4, D6, D8, E3, E5, E7, E9, E10, E11;
#     "HL+R"      with R: D2 in place of D2+R, and the region ER;
#     "HL+R+G", "HL+R+G+P", "HL+R+G+L"   with R, and with G, G union P or G union L removed from the regions E_k.
# The total loss of "HL+R+G+L" is Phi(theta) of (18), and that of "HL+R+G+P" is Phi_P(theta).
#
# First method (`losses_first`), deterministic:
#   * D_k: the b-integral is exact, int_{b0}^{b1} omega((1-a-b)/b)/b^2 db = (Omega((1-a)/b0 - 1) - Omega((1-a)/b1 - 1))
#     /(1 - a); composite Gauss-Legendre in a, split at the region's kinks.
#   * E_k: the d-integral is exact over [dlo, dhi] minus the union of the rescued d-intervals, through
#     int_{d0}^{d1} omega((1-s-d)/d)/d^2 dd = (Omega((1-s)/d0 - 1) - Omega((1-s)/d1 - 1))/(1 - s), s = a + b + g;
#     composite Gauss-Legendre in a, b, g.  The a-axis is split at the region's kinks, and every axis at the
#     planes where the rescued d-set jumps (jump_levels), so each piece has a continuous integrand.  The g-axis
#     is also split where dlo or dhi meets the cap (1 - a - b - g)/2 (omega jumps from 0 to 1 there) and where
#     dlo meets dhi (the d-range opens): the g-integrand has a kink at these points.  dlo and dhi are affine in g
#     between the jump planes, so the points are found exactly.
#   Resolution (k, panels): each axis is cut into `panels` equal panels, the panels are cut further at the points
#   above, and each resulting piece gets k Gauss-Legendre nodes.  The D_k use k nodes per piece and 16 * panels
#   panels on the a-axis.  The rule of Table 2 is k = 10, panels = 16: 256 panels on the a-axis for the D_k and
#   16 panels per axis for the E_k.  The a-nodes are distributed over the threads; each a-node's contribution is
#   stored and the contributions are summed in a fixed order, so the result does not depend on the threads.
#
# Second method (`losses_second`), independent of the Omega substitution and of the d-intervals:
#   * D_k: nested adaptive Gauss-Kronrod quadrature (QuadGK; on each piece between breakpoints the 21-point rule
#     and at most 200 subintervals, as QUADPACK's QAG with limit 200) of omega((1-a-b)/b)/(a b^2) itself, with
#     breakpoints at the kinks in a and at the b where (1-a-b)/b is an integer; its error estimate is that of the
#     outer integral;
#   * E_k: scrambled Sobol quasi-Monte Carlo on the unit cube mapped onto each region, omega evaluated directly,
#     the sets G, L, P tested pointwise (`rescued`).  The points are the first 2^m points of the Sobol sequence
#     in dimension 4 (Joe-Kuo direction numbers) with Matousek's linear matrix scramble and a random digital
#     shift, generated in Gray-code order; scrambling r uses the random number generator Xoshiro(seed + r).
#     The error estimate is the standard deviation of the scramblings divided by the square root of their number.

using FastGaussQuadrature: gausslegendre
using QuadGK: quadgk
using Roots: find_zero, Brent
using Random

@isdefined(Omega) || include("buchstab.jl")
@isdefined(params) || include("regions.jl")

const DECOMPOSITIONS = Dict(                # name => (with R, rescued set)
    "HL" => (false, "none"), "HL+R" => (true, "none"), "HL+R+G" => (true, "G"),
    "HL+R+G+P" => (true, "G+P"), "HL+R+G+L" => (true, "G+L"))
const ORDER = ("D1", "D2+R", "D2", "D4", "D6", "D8", "E3", "ER", "E5", "E7", "E9", "E10", "E11")   # Table 2

mode_code(mode) = findfirst(==(mode), MODES)

# ----------------------------------------------------------------------------------------------------------
# First method.

"""lo, the interior ends of `panels` equal panels, the cuts strictly inside (lo, hi), and hi, sorted."""
function breakpoints!(pts, lo, hi, cuts, panels)
    hi = max(hi, lo)
    empty!(pts)
    push!(pts, lo)
    for j in 1:panels - 1
        push!(pts, lo + (hi - lo) * j / panels)
    end
    for c in cuts
        (c > lo) & (c < hi) && push!(pts, c)
    end
    push!(pts, hi)
    return sort!(pts)
end

"""k Gauss-Legendre nodes and weights on each piece between consecutive breakpoints."""
function composite!(x, w, pts, glx, glw)
    empty!(x)
    empty!(w)
    for i in 1:length(pts) - 1
        len = pts[i + 1] - pts[i]
        len > 0 || continue
        for j in eachindex(glx)
            push!(x, pts[i] + len * (glx[j] + 1) / 2)
            push!(w, len * glw[j] / 2)
        end
    end
    return x, w
end

"""First and last zero of h, affine on each piece between consecutive breakpoints (sampled just inside each
piece, since it may jump at a breakpoint); +-Inf where there is none."""
function zero_range(pts, h)
    zmin, zmax = Inf, -Inf
    for i in 1:length(pts) - 1
        x0, x1 = pts[i], pts[i + 1]
        y0, y1 = x0 + 1e-9 * (x1 - x0), x1 - 1e-9 * (x1 - x0)
        h0, h1 = h(y0), h(y1)
        if h0 * h1 < 0
            z = y0 - h0 * (y1 - y0) / (h1 - h0)
            zmin, zmax = min(zmin, z), max(zmax, z)
        end
    end
    return zmin, zmax
end

"""Omega((1 - s)/d - 1): the d-antiderivative, up to the factor -1/(1 - s), of omega((1-s-d)/d)/d^2."""
@inline Fd(s, d) = Omega((1 - s) / d - 1)

"""Sum of Fd(lo) - Fd(hi) over the union of the n intervals (lo, hi), clipped to [d0, d1]."""
@inline function union_integral(s, d0, d1, lo, hi, n)
    @inbounds for i in 1:n
        l = clamp(lo[i], d0, d1)
        lo[i], hi[i] = l, max(clamp(hi[i], d0, d1), l)
    end
    @inbounds for i in 2:n                           # insertion sort by lo
        l, h, j = lo[i], hi[i], i - 1
        while j >= 1 && lo[j] > l
            lo[j + 1], hi[j + 1] = lo[j], hi[j]
            j -= 1
        end
        lo[j + 1], hi[j + 1] = l, h
    end
    total, before = 0.0, d0
    @inbounds for i in 1:n
        nl, nh = max(lo[i], before), max(hi[i], before)
        nh > nl && (total += Fd(s, nl) - Fd(s, nh))
        before = max(before, hi[i])
    end
    return total
end

function loss2_first(reg, k, panels, glx, glw)
    pts, x, w = Float64[], Float64[], Float64[]
    composite!(x, w, breakpoints!(pts, reg.a0, reg.a1, reg.kinks, panels), glx, glw)
    total = 0.0
    for i in eachindex(x)
        a = x[i]
        b0, b1 = reg.blo(a), reg.bhi(a)
        b1 > b0 && (total += w[i] * (Omega((1 - a) / b0 - 1) - Omega((1 - a) / b1 - 1)) / (a * (1 - a)))
    end
    return total
end

"""The a-nodes and weights of a four-dimensional region: split at its kinks and at the jump levels."""
function alpha_nodes(reg, p, panels, glx, glw)
    pts, x, w = Float64[], Float64[], Float64[]
    composite!(x, w, breakpoints!(pts, reg.a0, reg.a1, [reg.kinks; collect(jump_levels(p))], panels), glx, glw)
    return x, w
end

"""The contribution of one a-node (weight wa) to the loss of a four-dimensional region."""
function loss4_anode(reg, mode::Int, p, panels, glx, glw, a, wa)
    levels = jump_levels(p)
    pts, bx, bw, gx, gw, cuts = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]
    lo, hi = zeros(16), zeros(16)
    composite!(bx, bw, breakpoints!(pts, reg.blo(a), reg.bhi(a), [c - S for c in levels for S in (0.0, a)], panels),
               glx, glw)
    total = 0.0
    for (b, wb) in zip(bx, bw)
        glo, ghi = reg.glo(a, b), reg.ghi(a, b)
        jumps = [c - S for c in levels for S in (0.0, a, b, a + b)]
        breakpoints!(pts, glo, ghi, jumps, 1)
        empty!(cuts)
        append!(cuts, jumps)
        for h in (g -> reg.dlo(a, b, g) - (1 - a - b - g) / 2, g -> reg.dhi(a, b, g) - (1 - a - b - g) / 2,
                  g -> reg.dlo(a, b, g) - reg.dhi(a, b, g))
            append!(cuts, zero_range(pts, h))
        end
        composite!(gx, gw, breakpoints!(pts, glo, ghi, cuts, panels), glx, glw)
        for (g, wg) in zip(gx, gw)
            weight = wa * wb * wg
            weight > 0 || continue
            s = a + b + g
            d0, d1 = reg.dlo(a, b, g), reg.dhi(a, b, g)
            d1 > d0 || continue
            n = rescued_intervals!(lo, hi, mode, a, b, g, p)
            inner = Fd(s, d0) - Fd(s, d1) - union_integral(s, d0, d1, lo, hi, n)
            total += weight * inner / ((1 - s) * a * b * g)
        end
    end
    return total
end

"""Per-region losses (name => value, in the order of Table 2) with the first method."""
function losses_first(theta, decomposition; k = 10, panels = 16)
    rstar, mode = DECOMPOSITIONS[decomposition]
    code = mode_code(mode)
    p = params(theta)
    glx, glw = gausslegendre(k)
    out = Dict{String, Float64}()
    regs2 = regions_2d(theta, rstar)
    vals2 = zeros(length(regs2))
    Threads.@threads :dynamic for i in eachindex(regs2)
        vals2[i] = loss2_first(regs2[i], k, 16panels, glx, glw)
    end
    for (reg, v) in zip(regs2, vals2)
        out[reg.name] = get(out, reg.name, 0.0) + v
    end
    regs4 = regions_4d(theta, rstar)
    tasks = [(i, a, wa) for (i, reg) in enumerate(regs4) for (a, wa) in zip(alpha_nodes(reg, p, panels, glx, glw)...)]
    vals4 = zeros(length(tasks))
    Threads.@threads :dynamic for t in eachindex(tasks)
        i, a, wa = tasks[t]
        vals4[t] = loss4_anode(regs4[i], code, p, panels, glx, glw, a, wa)
    end
    for (t, (i, _, _)) in enumerate(tasks)
        out[regs4[i].name] = get(out, regs4[i].name, 0.0) + vals4[t]
    end
    return [n => out[n] for n in ORDER if haskey(out, n)]
end

total_first(theta, decomposition; kw...) = sum(last, losses_first(theta, decomposition; kw...))

# ----------------------------------------------------------------------------------------------------------
# Second method.

"""Adaptive quadrature of f over [x[1], x[end]], piece by piece between consecutive breakpoints: (value, sum of
the error estimates).  Each piece uses the 21-point Gauss-Kronrod rule and at most 200 subintervals, as QUADPACK's
QAG with limit 200."""
function quad_pieces(f, x; atol, rtol)
    value, error = 0.0, 0.0
    for i in 1:length(x) - 1
        x[i + 1] > x[i] || continue
        I, E = quadgk(f, x[i], x[i + 1]; atol, rtol, order = 10, maxevals = 200 * 21)
        value += I
        error += E
    end
    return value, error
end

"""The integrand b -> omega((1-a-b)/b)/(a b^2) at fixed a."""
struct InnerIntegrand
    a::Float64
end
(f::InnerIntegrand)(b) = omega((1 - f.a - b) / b) / (f.a * b * b)

"""a -> the integral over b in [blo(a), bhi(a)], split at the b where (1-a-b)/b is an integer."""
struct OuterIntegrand
    blo::Function
    bhi::Function
end
function (f::OuterIntegrand)(a)
    b0, b1 = f.blo(a), f.bhi(a)
    b1 <= b0 && return 0.0
    x = [b0; sort!([(1 - a) / (j + 1) for j in 1:7 if b0 < (1 - a) / (j + 1) < b1]); b1]
    return quad_pieces(InnerIntegrand(a), x; atol = 1e-13, rtol = 1e-12)[1]
end

"""Nested adaptive quadrature of omega((1-a-b)/b)/(a b^2) over a two-dimensional region (omega evaluated
directly), split at the kinks in a: (value, error estimate of the outer integral)."""
loss2_second(reg) = quad_pieces(OuterIntegrand(reg.blo, reg.bhi),
                                [reg.a0; sort!([x for x in reg.kinks if reg.a0 < x < reg.a1]); reg.a1];
                                atol = 1e-12, rtol = 1e-11)

# Sobol direction numbers in dimension 4 (Joe and Kuo, new-joe-kuo-6.21201): dimension 1 is the van der Corput
# sequence; the others are (s, a, m_1, ..., m_s).
const SOBOL_JOE_KUO = ((1, 0, (1,)), (2, 1, (1, 3)), (3, 1, (1, 3, 1)))

"""The 32 direction numbers (as 32-bit integers, most significant bit first) of each of the 4 dimensions."""
function sobol_directions()
    V = zeros(UInt32, 32, 4)
    for k in 1:32
        V[k, 1] = UInt32(1) << (32 - k)
    end
    for (j, (s, a, m0)) in enumerate(SOBOL_JOE_KUO)
        m = zeros(UInt64, 32)
        m[1:s] .= m0
        for k in s + 1:32
            mk = m[k - s] ⊻ (m[k - s] << s)
            for i in 1:s - 1
                (a >> (s - 1 - i)) & 1 == 1 && (mk ⊻= m[k - i] << i)
            end
            m[k] = mk
        end
        for k in 1:32
            V[k, j + 1] = UInt32(m[k] << (32 - k))
        end
    end
    return V
end

const SOBOL_V = sobol_directions()

"""Matousek's linear matrix scramble and a digital shift of the direction numbers, from rng: each output bit r of
a dimension is the parity of a random lower-triangular row (diagonal 1) applied to the input bits."""
function scrambled_directions(rng)
    V = similar(SOBOL_V)
    shift = zeros(UInt32, 4)
    for j in 1:4
        rows = [(rand(rng, UInt32) & ~(typemax(UInt32) >> (r - 1))) | (UInt32(1) << (32 - r)) for r in 1:32]
        for k in 1:32
            v = UInt32(0)
            for r in 1:32
                isodd(count_ones(rows[r] & SOBOL_V[k, j])) && (v |= UInt32(1) << (32 - r))
            end
            V[k, j] = v
        end
        shift[j] = rand(rng, UInt32)
    end
    return V, shift
end

"""One scrambled Sobol estimate of the loss of a four-dimensional region, with the first 2^m points."""
function qmc_estimate(reg, mode::Int, p, m, seed)
    V, x = scrambled_directions(Xoshiro(seed))
    acc = 0.0
    N = 1 << m
    for i in 0:N - 1
        if i > 0
            c = trailing_zeros(i) + 1                       # Gray code: flip the direction of the lowest set bit
            for j in 1:4
                x[j] ⊻= V[c, j]
            end
        end
        u1, u2, u3, u4 = (x[j] / 2.0^32 for j in 1:4)
        a = reg.a0 + (reg.a1 - reg.a0) * u1
        b0 = reg.blo(a)
        lb = max(reg.bhi(a) - b0, 0.0)
        b = b0 + lb * u2
        g0 = reg.glo(a, b)
        lg = max(reg.ghi(a, b) - g0, 0.0)
        g = g0 + lg * u3
        s = a + b + g
        d0 = reg.dlo(a, b, g)
        ld = max(min(reg.dhi(a, b, g), (1 - s) / 2) - d0, 0.0)
        d = d0 + ld * u4
        jac = (reg.a1 - reg.a0) * lb * lg * ld
        jac > 0 || continue
        rescued(mode, a, b, g, d, p) && continue
        acc += jac * omega((1 - s - d) / d) / (a * b * g * d * d)
    end
    return acc / N
end

"""Per-region (value, error estimate) with the second method, in the order of Table 2."""
function losses_second(theta, decomposition; m = 20, reps = 8, seed = 1)
    rstar, mode = DECOMPOSITIONS[decomposition]
    code = mode_code(mode)
    p = params(theta)
    out = Dict{String, Tuple{Float64, Float64}}()
    regs2 = regions_2d(theta, rstar)
    vals2 = Vector{Tuple{Float64, Float64}}(undef, length(regs2))
    Threads.@threads :dynamic for i in eachindex(regs2)
        vals2[i] = loss2_second(regs2[i])
    end
    for (reg, (v, e)) in zip(regs2, vals2)
        v0, e0 = get(out, reg.name, (0.0, 0.0))
        out[reg.name] = (v0 + v, sqrt(e0^2 + e^2))
    end
    regs4 = regions_4d(theta, rstar)
    est = zeros(length(regs4), reps)
    Threads.@threads :dynamic for t in 1:length(regs4) * reps
        i, r = divrem(t - 1, reps) .+ 1
        est[i, r] = qmc_estimate(regs4[i], code, p, m, seed + r - 1)
    end
    for (i, reg) in enumerate(regs4)
        e = est[i, :]
        mu = sum(e) / reps
        out[reg.name] = (mu, sqrt(sum((e .- mu) .^ 2) / (reps - 1)) / sqrt(reps))
    end
    return [n => out[n] for n in ORDER if haskey(out, n)]
end

# ----------------------------------------------------------------------------------------------------------

"""The theta at which the total loss of the first method is 1 (Brent's method, xtol 1e-10), logging every
evaluation."""
function theta_root(decomposition, bracket; k = 10, panels = 16, log = println)
    function f(theta)
        t0 = time()
        v = total_first(theta, decomposition; k, panels)
        log(@sprintf("    %-10s theta = %.10f  total = %.9f  (%.1fs)", decomposition, theta, v, time() - t0))
        return v - 1
    end
    return find_zero(f, bracket, Brent(); xatol = 1e-10, xrtol = 0.0)
end

if abspath(PROGRAM_FILE) == @__FILE__
    theta, decomposition = parse(Float64, ARGS[1]), ARGS[2]
    out = losses_first(theta, decomposition)
    for (n, v) in out
        @printf("%-6s %.7f\n", n, v)
    end
    total = sum(last, out)
    @printf("total  %.7f   1 - total %.7f\n", total, 1 - total)
end
