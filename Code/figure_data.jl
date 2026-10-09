# The data of the figures of the paper, from the definitions in regions.jl and the first method of losses.jl.
# The figures themselves are drawn by figures/plot_figures.py.
#
#     julia --project=. -t 32 figure_data.jl
#
#   figures/domain.json       the domain D of (12) at theta = 0.1129, its Ranges 1-6 (Table 1), the discarded regions
#                             D1, D2, D4, D6, D8 and the region R (13), as sampled boundaries
#   figures/total_loss.json   the total loss of each decomposition of Table 3 for 0.112 <= theta <= 0.118, and the
#                             theta of Table 3 at which it equals 1

using JSON

include("losses.jl")

const N_SAMPLES = 600
const ROWS = ["HL" => "Harman and Lewis", "HL+R" => "with R", "HL+R+G" => "with R and G",
              "HL+R+G+P" => "with R, G and P", "HL+R+G+L" => "with R, G and L"]   # Table 3

"""beta_0(alpha) of Section 3: tau - alpha for alpha >= sigma, (tau - alpha)/(j + 2) below, j = floor((sigma - alpha)/theta)."""
beta0(a, p) = a >= p.sigma ? p.tau - a : (p.tau - a) / (floor((p.sigma - a) / p.theta) + 2)

"""{a0 <= a <= a1, lo(a) <= b <= hi(a)} sampled at N_SAMPLES values of a, keeping the a with hi > lo."""
function sampled(name, a0, a1, lo, hi)
    as = collect(range(a0, a1, length = N_SAMPLES))
    keep = [hi(a) > lo(a) for a in as]
    return JSON.Object{String, Any}(["name" => name, "alpha" => as[keep], "lower" => lo.(as[keep]),
                                     "upper" => hi.(as[keep])])
end

function domain_data()
    theta = 0.1129
    p = params(theta)
    s, t, r, z = p.sigma, p.tau, p.rho, p.zeta
    upper(a) = min(a, (1 - a) / 2)
    # D, one piece per value of j, since beta_0 jumps where j changes; D starts where beta_0(alpha) = alpha
    pieces = [sampled("D", a0 + 1e-12, a1 - 1e-12, a -> beta0(a, p), upper)
              for (a0, a1) in ((t / 5, s - 2theta), (s - 2theta, s - theta), (s - theta, s), (s, 0.5))]
    discarded = [sampled(reg.name, reg.a0, reg.a1, reg.blo, reg.bhi) for reg in regions_2d(theta, true)]
    R = sampled("R", 1 / 3, s, a -> max(z - a, (t - a) / 2), a -> r / 2)
    # Range boundaries inside D: alpha = sigma, 1/3, rho, sigma - theta, and alpha + beta = sigma - theta
    verticals = [[a, beta0(a + 1e-9, p), upper(a)] for a in (s, 1 / 3, r, s - theta)]
    as = range(t / 5, s - theta, length = N_SAMPLES)
    keep = [beta0(a + 1e-12, p) <= s - theta - a <= upper(a) for a in as]
    diagonal = JSON.Object{String, Any}(["alpha" => collect(as[keep]), "beta" => s - theta .- collect(as[keep])])
    return JSON.Object{String, Any}(["theta" => theta, "sigma" => s, "tau" => t, "rho" => r, "zeta" => z,
                                     "D" => pieces, "discarded" => discarded, "R" => R,
                                     "range_verticals" => verticals, "range_diagonal" => diagonal])
end

function total_loss_data()
    thetas = collect(0.112:0.0001:0.118)
    roots = JSON.parsefile(joinpath(@__DIR__, "results", "summary.json"))["roots"]["values"]
    curves = [JSON.Object{String, Any}(["decomposition" => dec, "label" => lab, "root" => roots[dec],
                                        "total" => [total_first(th, dec) for th in thetas]]) for (dec, lab) in ROWS]
    return JSON.Object{String, Any}(["method" => "first method, k = 10, panels = 16", "theta" => thetas,
                                     "curves" => curves])
end

if abspath(PROGRAM_FILE) == @__FILE__
    for (file, data) in (("domain.json", domain_data()), ("total_loss.json", total_loss_data()))
        write(joinpath(@__DIR__, "figures", file), JSON.json(data; pretty = true) * "\n")
        println("wrote figures/", file)
    end
end
