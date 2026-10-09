# Driver: every number of Section 10 and Appendix D of the paper (Paper/main.tex), written to Code/results/.
#
#     julia --project=. -t 32,1 run_all.jl
#
# Steps and outputs (plain text, plus results/summary.json; progress in results/run_all.log):
#   buchstab      buchstab.log                     accuracy of omega and Omega
#   inequalities  table_inequalities.txt           Table 4 (Appendix D), exact
#   hl_0119       table_hl_0119.txt                the decomposition of Harman and Lewis at theta = 0.119
#   roots         table_roots.txt, roots.log       Table 3: the theta at which the total loss is 1
#   losses        table_losses_0.1129.txt          Table 2: the losses and Phi(0.1129), both methods
#                 table_losses_P_0.1133.txt        the losses and Phi_P(0.1133), both methods
#   convergence   table_convergence_0.1129.txt     the first method under refinement
#   paper         paper_numbers.txt                every number of Section 10 against the number printed there

using Dates
using JSON
using Printf

include("losses.jl")
include("inequalities.jl")

const OUT = joinpath(@__DIR__, "results")
const RES = (k = 10, panels = 16)               # the rule of Section 10 for the first method
const QMC = (m = 20, reps = 8)                  # second method: 8 scramblings of 2^20 Sobol points
const LADDER = [(6, 2), (6, 4), (8, 4), (8, 8), (10, 8), (8, 16), (10, 16), (12, 16), (12, 32)]
const BRACKETS = ["HL" => (0.116, 0.119), "HL+R" => (0.1115, 0.1175), "HL+R+G" => (0.1115, 0.1175),
                  "HL+R+G+P" => (0.1115, 0.1175), "HL+R+G+L" => (0.1115, 0.1175)]
const ROWS_TABLE3 = Dict("HL" => "Harman and Lewis", "HL+R" => "with R", "HL+R+G" => "with R and G",
                         "HL+R+G+P" => "with R, G and P", "HL+R+G+L" => "with R, G and L")
const HL_BOUNDS = Dict("D1" => 0.275, "D2+R" => 0.341, "D4" => 0.156, "D6" => 0.112, "D8" => 0.036, "E3" => 0.001,
                       "E5" => 0.001, "E7" => 0.001, "E9" => 0.001, "E10" => 0.017, "E11" => 0.001)
                                                # HL's bounds for phi_k, pp. 132-134; total < 0.943

# Every number of Section 10, as printed: (location, name, printed).  A printed "10^-k" is a bound on the order of
# magnitude, a printed ">x" a lower bound; every other number must equal the computed value rounded to the
# printed digits.
const TABLE2 = ["D1" => ("0.3201605", "0.3201605"), "D2" => ("0.3437028", "0.3437028"), "D4" => ("0.2132922", "0.2132922"),
                "D6" => ("0.0480512", "0.0480512"), "D8" => ("0.0684437", "0.0684437"), "E3" => ("0.0000098", "0.0000098"),
                "ER" => ("0.0010206", "0.0010207"), "E5" => ("0.0006199", "0.0006199"), "E7" => ("0.0004451", "0.0004451"),
                "E9" => ("0.0002727", "0.0002726"), "E10" => ("0.0028858", "0.0028870"), "E11" => ("0.0000988", "0.0000988")]
const TABLE3 = ["HL" => ("0.1167825", "0.2664351"), "HL+R" => ("0.1141073", "0.2717854"),
                "HL+R+G" => ("0.1139312", "0.2721376"), "HL+R+G+P" => ("0.1132049", "0.2735902"),
                "HL+R+G+L" => ("0.1128554", "0.2742892")]
const PAPER = vcat(
    [("Section 10", "omega and Omega against their closed forms on [1, 3]", "10^-14")],
    [("Table 2", "$n at 0.1129, $m method", v[i]) for (n, v) in TABLE2 for (i, m) in enumerate(("first", "second"))],
    [("Table 2", "Phi(0.1129), first method", "0.9990033"), ("Table 2", "Phi(0.1129), second method", "0.9990044"),
     ("Section 10", "total of HL at 0.119, first method", "0.9425585")],
    [("Table 3", "$(ROWS_TABLE3[dec]): $c", v[i]) for (dec, v) in TABLE3 for (i, c) in enumerate(("theta", "1/2 - 2 theta"))],
    [("Section 10", "change of Phi(0.1129) from the rule of Table 2 to the finest rule", "2.6e-10"),
     ("Section 10", "error estimate of the second method over the E_k at 0.1129", "1.0e-6"),
     ("Section 10", "difference of the two methods in Phi(0.1129)", "1.2e-6"),
     ("Section 10", "difference of the two methods over the error estimate", "1.1"),
     ("Section 10", "1 - Phi(0.1129)", "9.97e-4"),
     ("Section 10", "1 - Phi(0.1129) over the largest of the three quantities above", ">800"),
     ("Section 10", "Phi_P(0.1133), first method", "0.9978243"),
     ("Section 10", "Phi_P(0.1133), second method", "0.9978250"),
     ("Section 10", "1/2 - 2 * 0.1133", "0.2734")])

const LOGFILE = joinpath(OUT, "run_all.log")
const LOGLOCK = ReentrantLock()

function logline(msg)
    line = "[$(Dates.format(now(), "HH:MM:SS"))] $msg"
    lock(LOGLOCK) do
        println(line)
        flush(stdout)
        open(io -> println(io, line), LOGFILE, "a")
    end
end

"""Run f, logging a progress line every 20 s from the interactive thread while it runs."""
function with_heartbeat(f, label)
    done = Threads.Atomic{Bool}(false)
    t0 = time()
    if Threads.nthreads(:interactive) > 0
        Threads.@spawn :interactive while true
            sleep(20)
            done[] && break
            logline("  ... $label: running, $(round(Int, time() - t0))s")
        end
    end
    try
        return f()
    finally
        done[] = true
    end
end

function write_result(name, text)
    write(joinpath(OUT, name), text * "\n")
    logline("  wrote results/$name")
end

"""The value written with the digits of the printed number."""
function formatted(value, printed)
    (startswith(printed, "10^") || startswith(printed, ">")) && return @sprintf("%.3g", value)
    if occursin('e', printed)
        mantissa = split(printed, 'e')[1]
        digits = occursin('.', mantissa) ? length(split(mantissa, '.')[2]) : 0
        s = Printf.format(Printf.Format("%.$(digits)e"), value)
        m = match(r"^(.*)e([+-])0*(\d+)$", s)
        return m[1] * "e" * (m[2] == "-" ? "-" : "") * m[3]
    end
    return Printf.format(Printf.Format("%.$(length(split(printed, '.')[2]))f"), value)
end

"""Whether the computed value is the printed number (see PAPER)."""
function agrees(value, printed)
    startswith(printed, "10^") && return floor(Int, log10(abs(value))) <= parse(Int, printed[4:end])
    startswith(printed, ">") && return value > parse(Float64, printed[2:end])
    return formatted(value, printed) == printed
end

"""Per-region table: first method, second method with its error estimate, and their difference."""
function loss_table(title, fm, sm, total)
    head = @sprintf("%-12s %11s %11s %10s %10s", "region", "first", "second", "error est.", "difference")
    lines = [title, "", head, "-"^length(head)]
    s = Dict(sm)
    for (n, v) in fm
        c, e = s[n]
        push!(lines, @sprintf("%-12s %11.7f %11.7f %10.1e %+10.1e", n, v, c, e, v - c))
    end
    tot, ctot = sum(last, fm), sum(first, values(s))
    se = sqrt(sum(last(x)^2 for x in values(s)))
    push!(lines, @sprintf("%-12s %11.7f %11.7f %10.1e %+10.1e", total, tot, ctot, se, tot - ctot))
    push!(lines, @sprintf("%-12s %11.7f %11.7f", "1 - " * split(total, '(')[1], 1 - tot, 1 - ctot))
    return join(lines, "\n")
end

obj(pairs) = JSON.Object{String, Any}(pairs)

function main()
    mkpath(OUT)
    write(LOGFILE, "")
    summary = JSON.Object{String, Any}()
    summary["generated"] = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    summary["julia"] = string(VERSION)
    summary["threads"] = Threads.nthreads()
    summary["resolution"] = obj(["k" => RES.k, "panels" => RES.panels])
    summary["qmc"] = obj(["m" => QMC.m, "reps" => QMC.reps])
    runtimes = JSON.Object{String, Any}()
    summary["runtimes"] = runtimes
    computed = Dict{String, Float64}()          # name in PAPER => computed value
    t_all = time()

    # 1. Buchstab's function.
    t0 = time()
    logline("step buchstab")
    b = buchstab_check()
    write_result("buchstab.log", join([format_check(n, v) for (n, v) in b], "\n"))
    summary["buchstab"] = obj(b)
    computed["omega and Omega against their closed forms on [1, 3]"] =
        maximum(v for (n, v) in b if occursin("on [1,2]", n) || occursin("on [2,3]", n))
    runtimes["buchstab"] = time() - t0

    # 2. Table 4 (Appendix D).
    t0 = time()
    logline("step inequalities")
    rows = inequalities_evaluate()
    text, failures = inequalities_table(rows)
    write_result("table_inequalities.txt",
        "Table 4 (Appendix D): the inequalities between the parameters, written as quantities that must be positive.\n" *
        "For each row, in exact arithmetic: the margin of the inequality (right side minus left side) is the\n" *
        "quantity times 'multiple'; the quantity at theta = 0.1129 and 0.1133; and its value at 0.1129 against\n" *
        "the value printed in the table.\n\n" * text * "\n\n$(length(rows)) rows; failures: " *
        (isempty(failures) ? "none" : string(failures)))
    summary["inequalities"] = obj(["failures" => failures, "rows" => [obj([
        "quantity" => r.quantity, "inequality" => r.inequality, "value_0_1129" => ratstring(r.value),
        "value_0_1133" => ratstring(r.value_p), "multiple" => r.multiple === nothing ? nothing : ratstring(r.multiple),
        "printed" => r.printed, "agrees" => r.agrees, "ok" => r.ok]) for r in rows]])
    logline("  $(length(rows)) rows, failures: $(length(failures))")
    runtimes["inequalities"] = time() - t0

    # 3. The decomposition of Harman and Lewis at 0.119.
    t0 = time()
    logline("step hl_0119")
    fm, sm = with_heartbeat("hl_0119") do
        losses_first(0.119, "HL"; RES...), losses_second(0.119, "HL"; QMC...)
    end
    s = Dict(sm)
    lines = ["The decomposition of Harman and Lewis at theta = 0.119, with the regions derived from the lemmas (Remark 27),",
             "by both methods, against the bounds of Harman and Lewis (pp. 132-134).", "",
             @sprintf("%-7s %11s %11s %10s %9s  note", "region", "first", "second", "error est.", "HL bound")]
    for (n, v) in fm
        push!(lines, @sprintf("%-7s %11.7f %11.7f %10.1e %9.3f  %s", n, v, s[n][1], s[n][2], HL_BOUNDS[n],
                              v >= HL_BOUNDS[n] ? "exceeds the bound" : ""))
    end
    push!(lines, @sprintf("%-7s %11.7f %11.7f %10s %9s", "total", sum(last, fm), sum(first, values(s)), "", "< 0.943"))
    write_result("table_hl_0119.txt", join(lines, "\n"))
    computed["total of HL at 0.119, first method"] = sum(last, fm)
    summary["hl_0119"] = obj(["first" => obj(fm), "second" => obj([n => collect(v) for (n, v) in sm])])
    runtimes["hl_0119"] = time() - t0

    # 4. Table 3: the theta at which the total loss is 1.
    t0 = time()
    logline("step roots")
    roots, root_checks, trace = JSON.Object{String, Any}(), JSON.Object{String, Any}(), String[]
    tlog(msg) = (push!(trace, msg); logline(msg))
    for (dec, bracket) in BRACKETS
        t1 = time()
        roots[dec] = with_heartbeat("root $dec") do
            theta_root(dec, bracket; RES..., log = tlog)
        end
        tlog(@sprintf("  %s: total loss 1 at theta = %.7f  (%.0fs)", dec, roots[dec], time() - t1))
        root_checks[dec] = failing(roots[dec])
        tlog("  quantities of Table 4 at this theta: " * (isempty(root_checks[dec]) ? "all positive" : string(root_checks[dec])))
    end
    write_result("roots.log", join(trace, "\n"))
    lines = ["Table 3: the theta at which the total loss is 1 (Brent's method, xtol 1e-10; first method, " *
             "k = $(RES.k), panels = $(RES.panels)),", "and the exponent 1/2 - 2 theta.", "",
             @sprintf("%-24s %10s %14s", "decomposition", "theta", "1/2 - 2 theta")]
    for (dec, root) in roots
        push!(lines, @sprintf("%-24s %10.7f %14.7f", ROWS_TABLE3[dec], root, 0.5 - 2root))
        computed["$(ROWS_TABLE3[dec]): theta"] = root
        computed["$(ROWS_TABLE3[dec]): 1/2 - 2 theta"] = 0.5 - 2root
    end
    push!(lines, "\nThe quantities of Table 4 at each of these theta: " *
                 (all(isempty, values(root_checks)) ? "all positive." : string(root_checks)))
    write_result("table_roots.txt", join(lines, "\n"))
    summary["roots"] = obj(["values" => roots, "failing_quantities" => root_checks])
    runtimes["roots"] = time() - t0

    # 5. Table 2, Phi(0.1129), and Phi_P(0.1133), by both methods.
    t0 = time()
    logline("step losses")
    losses = JSON.Object{String, Any}()
    for (theta, dec, fname, total) in ((0.1129, "HL+R+G+L", "table_losses_0.1129.txt", "Phi(0.1129)"),
                                       (0.1133, "HL+R+G+P", "table_losses_P_0.1133.txt", "Phi_P(0.1133)"))
        fm, sm = with_heartbeat("losses $dec at $theta") do
            losses_first(theta, dec; RES...), losses_second(theta, dec; QMC...)
        end
        removed = endswith(dec, "L") ? "G and L" : "G and P"
        title = "The losses (13) at theta = $theta, with $removed removed from the regions E_k, and their sum $total.\n" *
                "First method: Gauss-Legendre, k = $(RES.k), panels = $(RES.panels).  Second method: nested adaptive " *
                "quadrature (D_k) and scrambled Sobol QMC,\n$(QMC.reps) x 2^$(QMC.m) points (E_k), with its error estimate."
        write_result(fname, loss_table(title, fm, sm, total))
        s = Dict(sm)
        tot, ctot = sum(last, fm), sum(first, values(s))
        losses["$dec at $theta"] = obj(["first" => obj(fm), "second" => obj([n => collect(v) for (n, v) in sm]),
                                        "first_total" => tot, "second_total" => ctot])
        if dec == "HL+R+G+L"
            for (n, v) in fm
                computed["$n at 0.1129, first method"] = v
                computed["$n at 0.1129, second method"] = s[n][1]
            end
            computed["Phi(0.1129), first method"] = tot
            computed["Phi(0.1129), second method"] = ctot
            computed["error estimate of the second method over the E_k at 0.1129"] =
                sqrt(sum(e^2 for (n, (_, e)) in sm if startswith(n, "E")))
            computed["difference of the two methods in Phi(0.1129)"] = abs(tot - ctot)
            computed["1 - Phi(0.1129)"] = 1 - tot
        else
            computed["Phi_P(0.1133), first method"] = tot
            computed["Phi_P(0.1133), second method"] = ctot
            computed["1/2 - 2 * 0.1133"] = 0.5 - 2 * 0.1133
        end
    end
    summary["losses"] = losses
    runtimes["losses"] = time() - t0

    # 6. Convergence of the first method at 0.1129.
    t0 = time()
    logline("step convergence")
    ladder = Pair{Tuple{Int, Int}, Vector{Pair{String, Float64}}}[]
    with_heartbeat("convergence") do
        for (k, panels) in LADDER
            t1 = time()
            out = losses_first(0.1129, "HL+R+G+L"; k, panels)
            push!(ladder, (k, panels) => out)
            logline(@sprintf("  k = %d, panels = %d: Phi(0.1129) = %.10f (%.0fs)", k, panels, sum(last, out), time() - t1))
        end
    end
    finest = Dict(last(ladder[end]))
    names = [n for (n, _) in last(ladder[end]) if startswith(n, "E")]
    ftot = sum(values(finest))
    lines = ["The first method at theta = 0.1129: the losses of the regions E_k (with G and L removed) and Phi(0.1129),",
             "and their differences from the finest rule (k = $(LADDER[end][1]), panels = $(LADDER[end][2])).  " *
             "The regions D_k use k nodes and 16 * panels panels on the alpha-axis.", "",
             @sprintf("%3s %6s ", "k", "panels") * join([@sprintf("%12s", n) for n in names], " ") *
             @sprintf(" %14s %10s", "Phi(0.1129)", "difference")]
    for ((k, panels), out) in ladder
        d, tot = Dict(out), sum(last, out)
        push!(lines, @sprintf("%3d %6d ", k, panels) * join([@sprintf("%12.9f", d[n]) for n in names], " ") *
                     @sprintf(" %14.10f %+10.1e", tot, tot - ftot))
    end
    push!(lines, "", "differences from the finest rule:")
    for ((k, panels), out) in ladder
        d = Dict(out)
        push!(lines, @sprintf("%3d %6d ", k, panels) * join([@sprintf("%+12.1e", d[n] - finest[n]) for n in names], " "))
    end
    change = sum(last, last(ladder[findfirst(x -> first(x) == (RES.k, RES.panels), ladder)])) - ftot
    push!(lines, @sprintf("\nrule of Table 2 (k = %d, panels = %d): Phi(0.1129) differs from the finest rule by %+.1e",
                          RES.k, RES.panels, change))
    write_result("table_convergence_0.1129.txt", join(lines, "\n"))
    computed["change of Phi(0.1129) from the rule of Table 2 to the finest rule"] = abs(change)
    summary["convergence"] = obj(["k=$k,panels=$panels" => obj(out) for ((k, panels), out) in ladder])
    runtimes["convergence"] = time() - t0

    # 7. Every number of Section 10 against the number printed in the paper.
    largest = max(computed["change of Phi(0.1129) from the rule of Table 2 to the finest rule"],
                  computed["error estimate of the second method over the E_k at 0.1129"],
                  computed["difference of the two methods in Phi(0.1129)"])
    computed["1 - Phi(0.1129) over the largest of the three quantities above"] = computed["1 - Phi(0.1129)"] / largest
    computed["difference of the two methods over the error estimate"] =
        computed["difference of the two methods in Phi(0.1129)"] / computed["error estimate of the second method over the E_k at 0.1129"]
    paper = [obj(["location" => loc, "name" => name, "printed" => printed, "computed" => computed[name],
                  "agrees" => agrees(computed[name], printed)]) for (loc, name, printed) in PAPER]
    head = @sprintf("%-11s %10s %10s  %-7s  number", "location", "printed", "computed", "")
    lines = ["Every number of Section 10 of the paper, as printed and as computed by this run (rounded to the printed digits).",
             "", head, "-"^length(head)]
    for r in paper
        push!(lines, @sprintf("%-11s %10s %10s  %-7s  %s", r["location"], r["printed"], formatted(r["computed"], r["printed"]),
                              r["agrees"] ? "agrees" : "DIFFERS", r["name"]))
    end
    differ = [r["name"] for r in paper if !r["agrees"]]
    push!(lines, "\n$(length(paper)) numbers; differ: " * (isempty(differ) ? "none" : string(differ)))
    write_result("paper_numbers.txt", join(lines, "\n"))
    summary["paper"] = paper

    runtimes["total"] = time() - t_all
    write(joinpath(OUT, "summary.json"), JSON.json(summary; pretty = true) * "\n")
    logline("wrote results/summary.json; numbers of Section 10 that differ from the paper: " *
            (isempty(differ) ? "none" : string(differ)) * @sprintf("; total runtime %.0fs; max RSS %d MB",
                                                                 runtimes["total"], round(Int, Sys.maxrss() / 2^20)))
end

main()
