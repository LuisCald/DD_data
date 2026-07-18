# ─────────────────────────────────────────────────────────────
# factors_to_coefficients.jl — reconstruct coefficient rows from the
# smoothed factors (the FactorMap bridge).
#
# WHY THIS WORKS (and its one hard limitation)
# The estimation model's own factors → coefficients chain is
#     coef = Gⱼ · x_smoothed, rescaled by per-block stds, plus means and trend
# (Reconstruction.jl in the DD_replication repo). Gⱼ, stds, means, and the
# HP trend are model exports that are NOT shipped here — so the model-based
# reconstruction is not reproducible from this repo. `FactorMap` sidesteps
# that: it RE-ESTIMATES the composite affine map by OLS from the two published
# files themselves (coef_t = α + Λ̂·F_4q). Under the `_average` trend
# convention the trend is a constant per coefficient, so the whole chain is
# one affine map and OLS recovers it exactly — median R² = 1.000 across all
# 1730 coefficients.
#
# LIMITATION: this only holds for the `_average` files. The `_normal` files
# carry the date-anchored HP trend, which is not a function of the factors —
# fitting on them fails (R² ≈ 0.4). Do not use this bridge for `_normal`.
#
# This is the middle stage of the pipeline
#
#     factors  →  coefficients  →  moments / micro data
#
# and the hook for counterfactuals: perturb a factor before predicting.
# Only the FIRST 8 state columns (x1..x8) are the distributional factors —
# x9..x32 are their lags (used via the 4-quarter average), x33..x43 the
# aggregate block. Example:
#
#     include("factors_to_coefficients.jl")          # in a session
#     F = factors_at(fm, "2008-Q3"); F[1] += 1.0     # shock factor 1
#     row = predict(fm, F)                           # counterfactual coefficients
#     quantile_from_row(row, :wealth, [0.1, 0.5, 0.9])
#
# Usage:
#   julia --project=. factors_to_coefficients.jl
#   # options: DATASET=PSID (fit target; also SCF, CEX)
#   # output:  <DATASET>_coefficients_reconstructed.csv (next to this script)
#
# Only deps: CSV, DataFrames.  Python twin: factors_to_coefficients.py.
# ─────────────────────────────────────────────────────────────

const HERE = @__DIR__
const REPO = abspath(joinpath(HERE, ".."))

include(joinpath(HERE, "reconstruct.jl"))
using .DistributionalReconstruction
using CSV, DataFrames

dataset = get(ENV, "DATASET", "PSID")
fm = FactorMap(joinpath(REPO, "data", "$(dataset)_coefficients_average.csv"),
               joinpath(REPO, "data", "smoothed_factors.csv"))
println(DistributionalReconstruction.summary(fm))

if abspath(PROGRAM_FILE) == @__FILE__
    rows = [predict(fm, fm.factors_4q[i, :]) for i in eachindex(fm.dates_used)]
    out = DataFrame(reduce(vcat, permutedims.(rows)), ["x$i" for i in 1:fm.n_coefs])
    out[!, "time"] = fm.dates_used
    dest = joinpath(HERE, "$(dataset)_coefficients_reconstructed.csv")
    CSV.write(dest, select(out, "time", :))
    println("Wrote $(nrow(out)) reconstructed coefficient rows → $dest")
end
