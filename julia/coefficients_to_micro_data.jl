# ─────────────────────────────────────────────────────────────
# coefficients_to_micro_data.jl — synthetic micro data from the
# published coefficient files.
# ─────────────────────────────────────────────────────────────
#
# WHAT THIS PRODUCES
# One weighted synthetic cross-section per quarter, in the same layout the
# model pipeline writes (see `construct_micro_dataset` in Correlations.jl and
# `create_micro_df`, called from CreateTimeSeries.jl): one row per cell of the
# 10×10×10 decile copula grid,
#
#   cop_share, grid_point, consum, income, wealth,
#   consumgrid, incomegrid, wealthgrid, time
#
#   cop_share  — probability mass of the cell (the household weight share;
#                model-implied, so isolated cells can be slightly negative)
#   grid_point — concatenated cell index, e.g. 372 = (consum decile 3,
#                income decile 7, wealth decile 2); deciles are 1 (bottom)
#                to 10 (top)
#   consum/income/wealth — the MEAN OF THE MEASURE AMONG THE HOUSEHOLDS IN
#                THE CELL (copula-weighted; VALUES=cell, default since Oct 2026),
#                RELATIVE TO THE PER-HOUSEHOLD MEAN of that quarter.
#                VALUES=decile gives the earlier layout: each cell carries the
#                unconditional decile average of each measure
#   *grid      — the cell indices as separate columns
#
# HOW IT FOLLOWS THE MODEL CODE
# * Weights: `copula_pmf_grid` reproduces the published `ciw_*` probability
#   masses (the copula block of the coefficient files).
# * Values (VALUES=cell): E[m | cell] = (1/w) ∫_cell sinh(Ξ⁻¹_m(u_m)) c(u) du.
#   The copula density is a Legendre polynomial, so this factors into
#   one-dimensional integrals per decile — the cell version of
#   RecessionPanels.box_conditional_means, which the paper's recession figure
#   uses. Gauss–Legendre (64 nodes per decile) on the outcome side, closed form
#   on the conditioning side. Cells with mass below MIN_CELL_MASS fall back to the
#   decile average; values a truncated copula pushes outside their own decile
#   are clipped at its bounds. A quarter without the consumption marginal
#   conditions on (income, wealth) only.
# * Values (VALUES=decile): the decile average of each measure (32 nodes).
# * Row layout: identical to `construct_micro_dataset` (Correlations.jl),
#   including the concatenated `grid_point` index.
#
# UNITS — SCALING TO LEVELS
# Values are relative to the per-household mean (the same convention as
# `quantile_at`; e.g. 0.73 = 73% of mean consumption that quarter). The model
# pipeline multiplies by a per-household aggregate series, shipped as
# `data/aggregate_anchors.csv` (time, consum_per_hh, income_per_hh,
# wealth_per_hh, tot_hhs — the same series the estimation uses). To get
# dollar levels, join on the quarter and compute
#
#   level_it = value_it × <measure>_per_hh_t
#
# (The published data/PSID_synthetic_microdata.csv and dd_microdata are already
# scaled this way, i.e. in 2019Q4 USD (CPI-U) per household.)
#
# To expand into unit-record data (N synthetic households), sample rows with
# probability `cop_share` (clip negatives to 0 and renormalize) and assign
# each household its row's values.
#
# USAGE
#   julia --project=. coefficients_to_micro_data.jl
#   # options: DATASET=PSID TREND=normal DATE=2008-Q3 (omit DATE for all quarters)
#   #          VALUES=cell (default) | decile
#   #          OUT=path.csv to write elsewhere
#   # output:  julia/coefficients_to_micro_data_output.csv (git-ignored; never
#   #          overwrites the published data/ files)
#
# Only deps: CSV, DataFrames (wraps the canonical helper julia/reconstruct.jl).
# Python twin: coefficients_to_micro_data.py.
# ─────────────────────────────────────────────────────────────

const HERE = @__DIR__
const REPO = abspath(joinpath(HERE, ".."))

include(joinpath(HERE, "reconstruct.jl"))
using .DistributionalReconstruction
using CSV, DataFrames, Statistics, LinearAlgebra
const DR = DistributionalReconstruction

const MEASURES = [:consum, :income, :wealth]   # alphabetical — model convention
const GRID = 10
const NODES_PER_DECILE = 32

"Decile-average quantiles: mean of the quantile function within each decile."
function decile_averages(r, date, measure)
    avgs = zeros(GRID)
    for d in 1:GRID
        us = (d - 1) / GRID .+ (collect(1:NODES_PER_DECILE) .- 0.5) ./ (GRID * NODES_PER_DECILE)
        avgs[d] = mean(quantile_at(r, date, measure, us))
    end
    return avgs
end

# Gauss–Legendre nodes/weights on [-1, 1] (Golub–Welsch)
function _gauss_legendre(n)
    β = [k / sqrt(4k^2 - 1) for k in 1:n-1]
    E = eigen(SymTridiagonal(zeros(n), β))
    return E.values, 2 .* E.vectors[1, :] .^ 2
end
const GL_X, GL_W = _gauss_legendre(64)
const EDGES = collect(0:GRID) ./ GRID
# below this copula mass a cell's conditional mean is numerically unstable (ratio of
# near-zero integrals); it carries the decile average instead (≈12% of cells, ≈0.1% of weight)
const MIN_CELL_MASS = 5e-5

"Mean of each measure among the households in each 10×10×10 cell (see header)."
function cell_means(r, date)
    κ = DR._kappa(r, date)
    K = size(κ, 1)
    U = [(EDGES[d] + EDGES[d+1]) / 2 + x / (2GRID) for d in 1:GRID, x in GL_X]      # (10, n)
    Wq = repeat((GL_W ./ (2GRID))', GRID)
    Iint = [I_m(o, [EDGES[d+1]])[1] - I_m(o, [EDGES[d]])[1] for o in 0:K-1, d in 1:GRID]   # (12, 10)
    has = Dict(m => !all(isnan, DR._xi(r, date, m)) for m in MEASURES)
    if !has[:consum]
        κ = copy(κ); κ[2:end, :, :] .= 0.0     # only the (income, wealth) slice is defined
    end
    Ic = has[:consum] ? Iint : repeat([o == 0 ? 1.0 : 0.0 for o in 0:K-1], 1, GRID)  # integrate consumption out
    side = Dict(:consum => Ic, :income => Iint, :wealth => Iint)
    contract(A, B, C) = [sum(κ[a, b, c] * A[a, i] * B[b, j] * C[c, k] for a in 1:K, b in 1:K, c in 1:K)
                         for i in 1:GRID, j in 1:GRID, k in 1:GRID]
    mass = contract(side[:consum], side[:income], side[:wealth])
    out = Dict{Symbol,Any}()
    for (mi, m) in enumerate(MEASURES)
        if !has[m]
            out[m] = nothing
            continue
        end
        f = reshape(quantile_at(r, date, m, vec(U)), size(U))
        dec = vec(sum(f .* Wq; dims = 2)) .* GRID
        F = [sum(Q_m(a - 1, U[d, :]) .* f[d, :] .* Wq[d, :]) for a in 1:K, d in 1:GRID]
        mats = [side[:consum], side[:income], side[:wealth]]
        mats[mi] = F
        num = contract(mats...)
        # range the quantile function spans within each decile (the fitted function
        # need not be monotone near an atom, e.g. zero wealth)
        lo, hi = vec(minimum(f; dims = 2)), vec(maximum(f; dims = 2))
        v = similar(num)
        for I in CartesianIndices(num)
            d = I[mi]
            v[I] = mass[I] > MIN_CELL_MASS ? clamp(num[I] / mass[I], lo[d], hi[d]) : dec[d]
        end
        out[m] = v
    end
    return out
end

"Decile-average layout (VALUES=decile): each cell carries the decile averages."
function _decile_values(r, date)
    q = Dict(m => decile_averages(r, date, m) for m in MEASURES)
    return Dict(m => [q[m][(i, j, k)[mi]] for i in 1:GRID, j in 1:GRID, k in 1:GRID]
                for (mi, m) in enumerate(MEASURES))
end

"10×10 (income, wealth) cell masses from the copula's consumption-free slice κ[1, :, :]."
function income_wealth_pmf(r, date)
    k2 = replace(DR._kappa(r, date)[1, :, :], NaN => 0.0)
    K = size(k2, 1)
    Imat = [I_m(o, [EDGES[d+1]])[1] for o in 0:K-1, d in 1:GRID]
    pad = zeros(GRID + 1, GRID + 1)
    pad[2:end, 2:end] = Imat' * k2 * Imat
    pmf = pad[2:end, 2:end] .- pad[1:end-1, 2:end] .- pad[2:end, 1:end-1] .+ pad[1:end-1, 1:end-1]
    pmf[pmf .< 0] .= 0.0
    return pmf ./ sum(pmf)
end

"One quarter's synthetic cross-section (1000 rows) in the model's layout."
function micro_rows(r, date; digits = 6, values = "cell")
    if all(isnan, DR._xi(r, date, :consum))
        # no consumption marginal: the (income, wealth) cells sit in the consumgrid == 1
        # rows; the other nine copies carry weight 0
        pmf = zeros(GRID, GRID, GRID)
        pmf[1, :, :] = income_wealth_pmf(r, date)
    else
        pmf = copula_pmf_grid(r, date)                   # 10×10×10, matches ciw_*
    end
    vals = values == "cell" ? cell_means(r, date) : _decile_values(r, date)
    val(m, i, j, k) = vals[m] === nothing ? NaN : round(vals[m][i, j, k]; digits)
    df = DataFrame(cop_share = Float64[], grid_point = Int[],
                   consum = Float64[], income = Float64[], wealth = Float64[],
                   consumgrid = Int[], incomegrid = Int[], wealthgrid = Int[])
    for i in 1:GRID, j in 1:GRID, k in 1:GRID            # same order as construct_micro_dataset
        push!(df, (round(pmf[i, j, k]; digits),
                   parse(Int, string(i) * string(j) * string(k)),
                   val(:consum, i, j, k), val(:income, i, j, k), val(:wealth, i, j, k), i, j, k))
    end
    df[!, :time] .= date
    return df
end

function main()
    dataset = get(ENV, "DATASET", "PSID")
    trend   = get(ENV, "TREND", "normal")
    date    = get(ENV, "DATE", "")            # empty → all quarters
    dest    = get(ENV, "OUT", joinpath(HERE, "coefficients_to_micro_data_output.csv"))
    values  = get(ENV, "VALUES", "cell")      # cell (default) | decile

    csv = joinpath(REPO, "data", "$(dataset)_coefficients_$(trend).csv")
    isfile(csv) || error("Coefficient file not found: $csv")
    r = Reconstruction(csv)

    dates = isempty(date) ? available_dates(r) : [date]
    out = DataFrame()
    for (n, d) in enumerate(dates)
        append!(out, micro_rows(r, d; values))
        n % 50 == 0 && println("  $n / $(length(dates)) quarters")
    end

    CSV.write(dest, out)
    println("Wrote $(nrow(out)) rows ($(length(dates)) quarters × $(GRID^3) cells) → $dest")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
