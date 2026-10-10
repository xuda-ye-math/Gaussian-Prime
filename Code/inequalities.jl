# Table 4 of the paper (Appendix D): the inequalities between the parameters, in exact rational arithmetic.
#
# Each row of the table is a quantity, linear in theta, that must be positive, the inequality between the
# parameters (params in regions.jl) that it expresses, and its use; a row may express two inequalities.  For each
# row, `inequalities_evaluate` computes the margin of each inequality (right side minus left side) from the
# parameters, and checks exactly that
#   * each margin is a positive multiple of the quantity, as linear forms in theta;
#   * the quantity is positive at theta = 0.1129 and at theta = 0.1133 (Section 10 uses both);
#   * the value of the quantity at theta = 0.1129, rounded to 7 decimals, is the value printed in the table.
# Lemma numbers are those printed in the paper.

using Printf

@isdefined(params) || include("regions.jl")

const Q = Rational{BigInt}
const THETA, THETA_P = Q(1129, 10000), Q(1133, 10000)

struct Row
    quantity::String        # the quantity, as printed in the table
    form::Tuple{Q, Q}       # (c1, c0): the quantity is c1 theta + c0
    inequality::String      # the inequality and its use, as printed in the table
    margins::Vector{Function}  # parameters -> right side minus left side, one per inequality of the row
    printed::String         # the value at theta = 0.1129, as printed in the table
end

const ROWS = [
    Row("theta - 1/10", (Q(1), Q(-1, 10)), "3 sigma/2 < tau, sigma - 3 theta < theta; Lemmas 12 and 13",
        [p -> p.tau - 3 * p.sigma / 2, p -> p.theta - (p.sigma - 3 * p.theta)], "0.0129000"),
    Row("1/8 - theta", (Q(-1), Q(1, 8)), "tau/5 < sigma - 2 theta; Lemma 13",
        [p -> p.sigma - 2 * p.theta - p.tau / 5], "0.0121000"),
    Row("(1 - 8 theta)/3", (Q(-8, 3), Q(1, 3)), "nu + rho < 1, 2 rho < tau; Lemma 9, case (b) fails in Ranges 4 and 5",
        [p -> 1 - p.nu - p.rho, p -> p.tau - 2 * p.rho], "0.0322667"),
    Row("1/6 - theta", (Q(-1), Q(1, 6)), "2 theta < sigma, tau < 2 sigma, 1/3 < sigma; Lemmas 2 and 9 to 13, Range 1",
        [p -> p.sigma - 2 * p.theta, p -> 2 * p.sigma - p.tau, p -> p.sigma - Q(1, 3)], "0.0537667"),
    Row("2 theta - 1/6", (Q(2), Q(-1, 6)), "sigma - theta < 1/3; j(X^alpha) = 0 in Range 2",
        [p -> Q(1, 3) - (p.sigma - p.theta)], "0.0591333"),
    Row("3 theta/2 - 1/12", (Q(3, 2), Q(-1, 12)), "sigma < (1/3 + tau)/2; PQ >= W in Range 2 and on R",
        [p -> (Q(1, 3) + p.tau) / 2 - p.sigma], "0.0860167"),
    Row("(3 - 22 theta)/12", (Q(-22, 12), Q(3, 12)), "2 zeta - 1/3 < nu, rho < 1/3; PQ^2 <= X^nu in Range 2, beta <= alpha on E3",
        [p -> p.nu - (2 * p.zeta - Q(1, 3)), p -> Q(1, 3) - p.rho], "0.0430167"),
    Row("(1 + 6 theta)/12", (Q(6, 12), Q(1, 12)), "zeta - 1/3 < rho; Q <= X^rho in Range 2",
        [p -> p.rho - (p.zeta - Q(1, 3))], "0.1397833"),
    Row("(1 - 2 theta)/6", (Q(-2, 6), Q(1, 6)), "zeta < tau; PQ <= X^(tau - eps) where alpha + beta <= zeta",
        [p -> p.tau - p.zeta], "0.1290333"),
    Row("(14 theta - 1)/6", (Q(14, 6), Q(-1, 6)), "sigma < zeta, zeta < 2 rho; P <= X^zeta on R, Q <= X^rho in Range 3, beta <= alpha on E3 and E5",
        [p -> p.zeta - p.sigma, p -> 2 * p.rho - p.zeta], "0.0967667"),
    Row("(5 - 34 theta)/24", (Q(-34, 24), Q(5, 24)), "zeta - sigma < rho/2, rho < sigma; R non-empty, Lemma 9",
        [p -> p.rho / 2 - (p.zeta - p.sigma), p -> p.sigma - p.rho], "0.0483917"),
    Row("(26 theta - 1)/24", (Q(26, 24), Q(-1, 24)), "sigma + rho/2 < tau; PQ <= X^(tau - eps) on R",
        [p -> p.tau - p.sigma - p.rho / 2], "0.0806417"),
    Row("(7 - 22 theta)/24", (Q(-22, 24), Q(7, 24)), "rho/2 < 1/3; beta < alpha on E_R",
        [p -> Q(1, 3) - p.rho / 2], "0.1881750"),
    Row("(1 - 7 theta)/3", (Q(-7, 3), Q(1, 3)), "tau/3 < sigma - theta; Lemma 13",
        [p -> p.sigma - p.theta - p.tau / 3], "0.0699000"),
    Row("(46 theta - 5)/12", (Q(46, 12), Q(-5, 12)), "sigma - theta < rho; j(X^alpha) = 0 in Range 3, Q <= X^rho in Range 5",
        [p -> p.rho - (p.sigma - p.theta)], "0.0161167"),
    Row("(58 theta - 5)/24", (Q(58, 24), Q(-5, 24)), "sigma < (rho + tau)/2; PQ >= W in Range 3",
        [p -> (p.rho + p.tau) / 2 - p.sigma], "0.0645083"),
    Row("(10 theta - 1)/12", (Q(10, 12), Q(-1, 12)), "rho < 8 theta/3; case (b) fails in Range 3",
        [p -> 8 * p.theta / 3 - p.rho], "0.0107500"),
    Row("(10 theta - 1)/6", (Q(10, 6), Q(-1, 6)), "1/3 < 2 zeta - tau; beta_0(alpha) < zeta - alpha on D4",
        [p -> 2 * p.zeta - p.tau - Q(1, 3)], "0.0215000"),
    Row("(1 - 2 theta)/24", (Q(-2, 24), Q(1, 24)), "chi < zeta; Lemma 9 in Range 1",
        [p -> p.zeta - p.chi], "0.0322583"),
    Row("(11 - 82 theta)/36", (Q(-82, 36), Q(11, 36)), "nu/3 < sigma - theta; D8 non-empty, beta <= alpha on E7",
        [p -> p.sigma - p.theta - p.nu / 3], "0.0483944"),
    Row("(1 - 2 theta)/12", (Q(-2, 12), Q(1, 12)), "tau < nu; beta_0(alpha) < (nu - alpha)/2 on D6 and D8",
        [p -> p.nu - p.tau], "0.0645167"),
    Row("(14 theta - 1)/9", (Q(14, 9), Q(-1, 9)), "sigma < 2 nu/3; alpha + beta > sigma on D8",
        [p -> 2 * p.nu / 3 - p.sigma], "0.0645111"),
    Row("(46 theta - 5)/24", (Q(46, 24), Q(-5, 24)), "(nu + sigma - theta)/2 < zeta; PQ <= X^zeta in Range 5",
        [p -> p.zeta - (p.nu + p.sigma - p.theta) / 2], "0.0080583"),
    Row("5 theta/3 - 1/6", (Q(5, 3), Q(-1, 6)), "sigma < 2 tau/3; PQ >= W on E9",
        [p -> 2 * p.tau / 3 - p.sigma], "0.0215000"),
    Row("5 theta/2 - 1/4", (Q(5, 2), Q(-1, 4)), "sigma - theta < tau/2; E10 lies in Range 5",
        [p -> p.tau / 2 - (p.sigma - p.theta)], "0.0322500"),
    Row("5 theta - 1/2", (Q(5), Q(-1, 2)), "sigma < 4 theta; alpha + beta < sigma for alpha < sigma - 2 theta",
        [p -> 4 * p.theta - p.sigma], "0.0645000"),
    Row("17 theta/5 - 3/10", (Q(17, 5), Q(-3, 10)), "sigma - 2 theta < 2 tau/5; |pq|^2 >= W X^(-2 theta) in Range 6",
        [p -> 2 * p.tau / 5 - (p.sigma - 2 * p.theta)], "0.0838600"),
    Row("21 theta/5 - 2/5", (Q(21, 5), Q(-2, 5)), "sigma - 3 theta < tau/5; j(X^alpha) <= 2 in Range 6",
        [p -> p.tau / 5 - (p.sigma - 3 * p.theta)], "0.0741800"),
]

quantity(row, theta) = row.form[1] * theta + row.form[2]

"""c with margin = c * quantity for every theta, or nothing if the margin is not a multiple of the quantity."""
function multiple(row, margin)
    t1, t2 = Q(1, 9), Q(1, 8)
    m1, m2 = margin(params(t1)), margin(params(t2))
    c = (m2 - m1) / (quantity(row, t2) - quantity(row, t1))
    all(margin(params(t)) == c * quantity(row, t) for t in (t1, t2, THETA)) || return nothing
    return c
end

"""The multiples of the quantity, one per inequality of the row."""
multiples(row) = [multiple(row, m) for m in row.margins]

"""The quantities that are not positive at theta."""
failing(theta) = [r.quantity for r in ROWS if quantity(r, Q(theta)) <= 0]

ratstring(x::Q) = denominator(x) == 1 ? string(numerator(x)) : "$(numerator(x))/$(denominator(x))"

"""One named tuple per row: exact values at THETA and THETA_P, the multiple, and the agreement with the table."""
function inequalities_evaluate()
    return [begin
                q, qp, cs = quantity(r, THETA), quantity(r, THETA_P), multiples(r)
                (quantity = r.quantity, inequality = r.inequality, value = q, value_p = qp, multiples = cs,
                 printed = r.printed, agrees = @sprintf("%.7f", Float64(q)) == r.printed,
                 ok = all(c -> c !== nothing && c > 0, cs) && q > 0 && qp > 0)
            end for r in ROWS]
end

"""Plain-text table of the rows, and the list of failures."""
function inequalities_table(rows)
    head = @sprintf("%-20s %10s %18s %10s %8s %10s  inequality and use", "quantity", "at 0.1129", "exact at 0.1129",
                    "at 0.1133", "multiple", "printed")
    lines = [head, "-"^length(head)]
    failures = String[]
    for r in rows
        c = join([m === nothing ? "none" : ratstring(m) for m in r.multiples], ",")
        push!(lines, @sprintf("%-20s %10.7f %18s %10.7f %8s %10s  %s", r.quantity, Float64(r.value),
                              ratstring(r.value), Float64(r.value_p), c, r.printed, r.inequality))
        r.ok || push!(failures, "$(r.quantity): not a positive multiple of its margin, or not positive")
        r.agrees || push!(failures, @sprintf("%s: %.7f differs from the printed %s", r.quantity, Float64(r.value),
                                             r.printed))
    end
    return join(lines, "\n"), failures
end

if abspath(PROGRAM_FILE) == @__FILE__
    text, failures = inequalities_table(inequalities_evaluate())
    println(text)
    println("\n$(length(ROWS)) rows; failures: ", isempty(failures) ? "none" : failures)
end
