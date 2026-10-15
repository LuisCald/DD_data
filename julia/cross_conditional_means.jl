# ─────────────────────────────────────────────────────────────
# cross_conditional_means.jl — cross-conditional group means from the
# published coefficient files: E[outcome | conditioning variable in group],
# e.g. average wealth of the bottom 50% / next 40% / top 10% of INCOME.
#
# Default method (METHOD=maxent): fit the positive maximum-entropy density
#     f(u, v) = exp( Σ_kl λ_kl Q_k(u) Q_l(v) ) / Z
# whose moments match the pairwise copula coefficients m_kl = E[Q_k(u)Q_l(v)],
# and integrate the EXACT group indicator over it. The truncated series
# Σ m_kl Q_k Q_l (METHOD=series) is equivalent to integrating an order-12
# approximation of the indicator, which rings and lets every group pick up a
# few percent of the other variable's top tail — biased for heavy-tailed
# outcomes (SCF 2019 microdata: wealth by middle-40% income −9% with the
# series, within ±1% with max-entropy). Damped Newton on the convex dual with
# a small ridge ε, so it converges even if the smoothed coefficients are not
# exactly the moments of a positive density. See docs/maxent_note.pdf; the
# python/cross_conditional_means.py is the twin.
#
#   julia --project=. julia/cross_conditional_means.jl
#   DATASET=SCF DATE=all julia --project=. julia/cross_conditional_means.jl
#   GROUPS=40,40,20 METHOD=series julia --project=. julia/cross_conditional_means.jl
#
# Output: julia/cross_conditional_means_output.csv, long format
#   date, outcome, conditioning, group, rel_mean, level, moment_error
# (rel_mean relative to the outcome's per-HH mean; level × <outcome>_per_hh).
#
# Only deps: CSV, DataFrames, LinearAlgebra (wraps julia/reconstruct.jl).
# ─────────────────────────────────────────────────────────────

const HERE = @__DIR__
const REPO = abspath(joinpath(HERE, ".."))

include(joinpath(HERE, "reconstruct.jl"))
using .DistributionalReconstruction
using CSV, DataFrames, Dates, LinearAlgebra
const DR = DistributionalReconstruction

const MEASURES = [:consum, :income, :wealth]       # axis order of the copula tensor
const K = 12
const DENSE_TOP = [0.0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999, 1.0]
const GL_X = [-0.9815606342467192, -0.9041172563704749, -0.7699026741943047, -0.5873179542866175,
              -0.3678314989981802, -0.1252334085114689, 0.1252334085114689, 0.3678314989981802,
              0.5873179542866175, 0.7699026741943047, 0.9041172563704749, 0.9815606342467192]
const GL_W = [0.0471753363865118, 0.1069393259953184, 0.1600783285433462, 0.2031674267230659,
              0.2334925365383548, 0.2491470458134028, 0.2491470458134028, 0.2334925365383548,
              0.2031674267230659, 0.1600783285433462, 0.1069393259953184, 0.0471753363865118]

Qb(k, u) = DR.Q_m(k, u)

"Gauss–Legendre grid, dense at the top, aligned with every group edge."
function make_grid(group_edges)
    edges = sort(unique(round.(vcat(DENSE_TOP, group_edges), digits = 12)))
    x = Float64[]; w = Float64[]
    for s in 1:length(edges)-1
        a, b = edges[s], edges[s+1]
        append!(x, (a + b) / 2 .+ (b - a) / 2 .* GL_X); append!(w, (b - a) / 2 .* GL_W)
    end
    Qx = vcat([Qb(k, x)' for k in 0:K-1]...)                                # K × n
    P = reshape([Qx[k, i] * Qx[k2, i] for k in 1:K, k2 in 1:K, i in 1:length(x)], K * K, length(x))
    return (x = x, w = w, Qx = Qx, P = P, W2 = w * w')
end

"""
    maxent_density_2d(M, g; λ0, ε, tol, maxit) -> (f, λ, err)

Masses f (n×n) of exp(Σ λ_kl Q_k(u)Q_l(v)) matching M (K×K, M[1,1] = 1), from the
regularised dual log Z(λ) − λ·M + ε/2‖λ‖² by damped Newton. Separable algebra:
log f = Qxᵀ Λ Qx, moments = Qx f Qxᵀ, Hessian = reshape(P f Pᵀ).
"""
function maxent_density_2d(M, g; λ0 = nothing, ε = 1e-6, tol = 1e-7, maxit = 200)
    free = trues(K, K); free[1, 1] = false
    idx = findall(vec(free))
    Λ = λ0 === nothing ? zeros(K, K) : copy(λ0)
    function state(Lv)                     # local names on purpose (no capture of the outer f/μ/D)
        ev = g.Qx' * Lv * g.Qx
        mx = maximum(ev)
        ew = exp.(ev .- mx) .* g.W2
        Z = sum(ew)
        fv = ew ./ Z
        μv = g.Qx * fv * g.Qx'
        Dv = log(Z) + mx - sum(Lv[idx] .* M[idx]) + ε / 2 * sum(abs2, Lv[idx])
        return fv, μv, Dv
    end
    f, μ, D = state(Λ)
    gr = (μ .- M)[idx] .+ ε .* Λ[idx]
    for _ in 1:maxit
        maximum(abs.(gr)) < tol && break
        S4 = reshape(g.P * f * g.P', K, K, K, K)
        H = reshape(permutedims(S4, (1, 3, 2, 4)), K * K, K * K) .- vec(μ) * vec(μ)'
        step = Symmetric(H[idx, idx] + ε * I) \ gr
        t = 1.0; accepted = false
        while t > 1e-8
            Ln = copy(Λ); Ln[idx] .-= t .* step
            fn, μn, Dn = state(Ln)
            if Dn <= D - 1e-4 * t * dot(gr, step)
                Λ, f, μ, D = Ln, fn, μn, Dn; accepted = true; break
            end
            t /= 2
        end
        accepted || break
        gr = (μ .- M)[idx] .+ ε .* Λ[idx]
    end
    E = abs.(μ .- M); E[1, 1] = 0.0
    return f, Λ, maximum(E)
end

"K×K moments, rows = conditioning degree, cols = outcome degree."
function pair_moments(κ, cond, out)
    ic, io = findfirst(==(cond), MEASURES), findfirst(==(out), MEASURES)
    ix = Any[1, 1, 1]; ix[ic] = Colon(); ix[io] = Colon()
    S = κ[ix...]
    return ic < io ? S : permutedims(S)
end

function exact_mean(r, date, m)
    s = 0.0
    for s_ in 1:length(DENSE_TOP)-1
        a, b = DENSE_TOP[s_], DENSE_TOP[s_+1]
        x = (a + b) / 2 .+ (b - a) / 2 .* GL_X
        s += sum((b - a) / 2 .* GL_W .* quantile_at(r, date, m, x))
    end
    return s
end

function group_labels(edges)
    n = length(edges) - 1
    [(i == 1 ? "bot" : i == n ? "top" : "mid") * string(round(Int, 100 * (edges[i+1] - edges[i]))) for i in 1:n]
end

function group_means(r, date, cond, out, edges, g, method, warm)
    M = pair_moments(DR._kappa(r, date), cond, out)
    labels = group_labels(edges)
    if method == "series"
        A = zeros(K)
        for s in 1:length(DENSE_TOP)-1
            a, b = DENSE_TOP[s], DENSE_TOP[s+1]
            x = (a + b) / 2 .+ (b - a) / 2 .* GL_X
            v = quantile_at(r, date, out, x)
            A .+= [sum((b - a) / 2 .* GL_W .* Qb(l, x) .* v) for l in 0:K-1]
        end
        return [(labels[i], dot([DR.I_m(k, edges[i+1]) - DR.I_m(k, edges[i]) for k in 0:K-1], M * A) / (edges[i+1] - edges[i]), NaN)
                for i in 1:length(labels)]
    end
    f, λ, err = maxent_density_2d(M, g; λ0 = get(warm, (cond, out), nothing))
    warm[(cond, out)] = λ
    EV = f * quantile_at(r, date, out, g.x)
    scale = exact_mean(r, date, out) / sum(EV)          # groups add up to the exact mean
    return [(labels[i], scale * sum(EV[(g.x .> edges[i]) .& (g.x .<= edges[i+1])]) / (edges[i+1] - edges[i]), err)
            for i in 1:length(labels)]
end

has_measure(r, date, m) = try
    all(isfinite, quantile_at(r, date, m, [0.5]))
catch
    false
end

function load_anchors()
    p = joinpath(REPO, "data", "aggregate_anchors.csv")
    isfile(p) || return Dict{String,Dict{Symbol,Float64}}()
    df = CSV.read(p, DataFrame)
    out = Dict{String,Dict{Symbol,Float64}}()
    for row in eachrow(df)
        dt = Date(string(row.time))
        out[string(year(dt), "-Q", quarterofyear(dt))] = Dict(m => float(row[Symbol(string(m, "_per_hh"))]) for m in MEASURES)
    end
    return out
end

function main()
    dataset = get(ENV, "DATASET", "PSID")
    trend   = get(ENV, "TREND", "normal")
    date    = get(ENV, "DATE", "2019-Q4")
    method  = get(ENV, "METHOD", "maxent")
    sizes   = parse.(Float64, split(get(ENV, "GROUPS", "50,40,10"), ",")) ./ 100
    abs(sum(sizes) - 1) < 1e-9 || error("GROUPS must add up to 100")
    edges = vcat(0.0, cumsum(sizes)); edges[end] = 1.0

    r = Reconstruction(joinpath(REPO, "data", "$(dataset)_coefficients_$(trend).csv"))
    dates = lowercase(date) == "all" ? available_dates(r) : [date]
    anchors = load_anchors()
    g = make_grid(edges)
    warm = Dict{Tuple{Symbol,Symbol},Matrix{Float64}}()
    rows = DataFrame(date = String[], outcome = String[], conditioning = String[], group = String[],
                     rel_mean = Float64[], level = Float64[], moment_error = Float64[])
    for d in dates
        present = [m for m in MEASURES if has_measure(r, d, m)]
        for cond in present, out in present
            cond == out && continue
            for (lab, rel, err) in group_means(r, d, cond, out, edges, g, method, warm)
                a = get(anchors, d, nothing)
                lvl = a === nothing ? NaN : rel * a[out]
                push!(rows, (d, string(out), string(cond), lab, rel, lvl, err))
            end
        end
    end
    outp = joinpath(HERE, "cross_conditional_means_output.csv")
    CSV.write(outp, rows)
    length(dates) == 1 && show(stdout, MIME("text/plain"), unstack(rows, [:outcome, :conditioning], :group, :rel_mean))
    println("\nWrote $(nrow(rows)) rows → $outp (method = $method)")
    method == "maxent" && nrow(rows) > 0 && println("max moment error of the max-entropy fits: ", maximum(rows.moment_error))
end

main()
