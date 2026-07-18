# ─────────────────────────────────────────────────────────────
# coefficients_to_moments.jl — turn the published coefficient files into
# moments: decile cut points and averages, levels, shares, Ginis, and the
# copula density at the joint median.
#
# Per measure (consum, income, wealth) and date:
#   * decile cut points  Ξ⁻¹(u), u = 0.1 … 0.9      (relative to per-HH mean)
#   * decile averages    mean of Ξ⁻¹ within each decile (32-node quadrature,
#                        mirroring the quadgk decile integral in the
#                        replication package's CreateTimeSeries.jl)
#   * levels             group average × <measure>_per_hh from
#                        data/aggregate_anchors.csv (dollars)
#   * shares             group share of the total, groups = bottom 50%
#                        (deciles 1–5), next 40% (6–9), top 10% (10) — the
#                        package's define_sequences grouping
#   * Gini               generalized (robust) Gini of Raffinetti, Siletti &
#                        Vernizzi (2015) on the 10 decile averages with equal
#                        weights — well-defined with negative wealth; ported
#                        from Distributional_Counterfactuals/Support.jl
#
# Usage:
#   julia --project=. julia/coefficients_to_moments.jl
#   # options: DATASET=PSID|SCF|CEX  TREND=normal|average  DATE=2008-Q3|all
#
# Output: julia/coefficients_to_moments_output.csv, long format
#   date, measure, stat, index, value
# (stat ∈ decile_cut / decile_avg / level_per_hh / share / gini /
#  copula_density_median).
#
# Only deps: CSV, DataFrames (wraps julia/reconstruct.jl).
# Python twin: coefficients_to_moments.py.
# ─────────────────────────────────────────────────────────────

const HERE = @__DIR__
const REPO = abspath(joinpath(HERE, ".."))

include(joinpath(HERE, "reconstruct.jl"))
using .DistributionalReconstruction
using CSV, DataFrames, Statistics, Dates

const MEASURES = [:consum, :income, :wealth]
const GRID = 10
const NODES_PER_DECILE = 32
const GROUPS = [("bot50", 1:5), ("mid40", 6:9), ("top10", 10:10)]

has_measure(r, date, m) = try
    quantile_at(r, date, m, [0.5]); true
catch
    false
end

"Mean of the quantile function within each decile (package convention)."
function decile_averages(r, date, m)
    avgs = zeros(GRID)
    for d in 1:GRID
        us = (d - 1) / GRID .+ (collect(1:NODES_PER_DECILE) .- 0.5) ./ (GRID * NODES_PER_DECILE)
        avgs[d] = mean(quantile_at(r, date, m, us))
    end
    return avgs
end

"""
Generalized (robust) Gini of Raffinetti, Siletti & Vernizzi (2015) — handles
negative values, stays in [0,1]. Ported from Distributional_Counterfactuals.
    G = Σᵢⱼ wᵢwⱼ|xᵢ−xⱼ| / [ 2(Σw − min w)(T⁺ + T⁻) ]
"""
function gini_robust(x::AbstractVector{<:Real}, w::AbstractVector{<:Real} = ones(length(x)))
    N = length(x)
    S = sum(w)
    Tplus = sum(w[i] * max(x[i], 0) for i in 1:N)
    Tminus = sum(w[i] * (-min(x[i], 0)) for i in 1:N)
    num = sum(w[i] * w[j] * abs(x[i] - x[j]) for i in 1:N, j in 1:N)
    denom = 2 * (S - minimum(w)) * (Tplus + Tminus)
    return num / denom
end

"per-HH anchor lookup: \"1962-Q3\" ⇒ Dict(measure ⇒ per_hh)."
function load_anchors()
    df = CSV.read(joinpath(REPO, "data", "aggregate_anchors.csv"), DataFrame)
    out = Dict{String,Dict{Symbol,Float64}}()
    for row in eachrow(df)
        dt = Date(row.time)
        key = string(year(dt), "-Q", quarterofyear(dt))
        out[key] = Dict(m => float(row[Symbol(string(m, "_per_hh"))]) for m in MEASURES)
    end
    return out
end

function moments_for_date!(rows, r, anchors, date; verbose = true)
    measures = filter(m -> has_measure(r, date, m), MEASURES)
    anchor = get(anchors, date, nothing)
    us = collect(0.1:0.1:0.9)

    for m in measures
        cuts = quantile_at(r, date, m, us)
        avgs = decile_averages(r, date, m)
        tot = sum(avgs)
        g = gini_robust(avgs, fill(0.1, GRID))

        for (u, v) in zip(us, cuts)
            push!(rows, (date, string(m), "decile_cut", string(u), v))
        end
        for d in 1:GRID
            push!(rows, (date, string(m), "decile_avg", string(d), avgs[d]))
        end
        for (name, ids) in GROUPS
            push!(rows, (date, string(m), "share", name, sum(avgs[ids]) / tot))
            if anchor !== nothing
                push!(rows, (date, string(m), "level_per_hh", name, mean(avgs[ids]) * anchor[m]))
            end
        end
        push!(rows, (date, string(m), "gini", "", g))

        if verbose
            println(rpad(string(m), 8),
                "deciles: ", join([string(round(v, digits = 3)) for v in avgs], "  "))
            println(rpad("", 8), "shares (bot50/mid40/top10): ",
                join([string(round(sum(avgs[ids]) / tot, digits = 3)) for (_, ids) in GROUPS], " / "),
                "   gini: ", round(g, digits = 3))
        end
    end

    if length(measures) == 3
        cmed = copula_density_at(r, date, 0.5, 0.5, 0.5)
        push!(rows, (date, "joint", "copula_density_median", "", cmed))
        verbose && println("\ncopula density at joint median: ", round(cmed, digits = 4))
    end
end

function main()
    dataset = get(ENV, "DATASET", "PSID")
    trend   = get(ENV, "TREND", "normal")
    date    = get(ENV, "DATE", "2008-Q3")

    csv = joinpath(REPO, "data", "$(dataset)_coefficients_$(trend).csv")
    isfile(csv) || error("Coefficient file not found: $csv")
    r = Reconstruction(csv)
    anchors = load_anchors()

    rows = DataFrame(date = String[], measure = String[], stat = String[],
                     index = String[], value = Float64[])

    if lowercase(date) == "all"
        for (n, d) in enumerate(available_dates(r))
            moments_for_date!(rows, r, anchors, d; verbose = false)
            n % 50 == 0 && println("  $n quarters done")
        end
    else
        println("Dataset $dataset ($trend trend), date $date\n")
        moments_for_date!(rows, r, anchors, date)
    end

    out = joinpath(HERE, "coefficients_to_moments_output.csv")
    CSV.write(out, rows)
    println("\nWrote $(nrow(rows)) rows → $out")
end

main()
