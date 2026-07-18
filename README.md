# Distributional Dynamics — Data

**Quarterly estimates of the joint distribution of consumption, income, and
wealth for U.S. households, 1962-Q3 to 2024-Q1**, with posterior uncertainty —
plus small scripts that turn them into whatever moment or micro dataset you
need.

This is the *data product* of:

> Bayer, Christian, Luis Calderon, and Moritz Kuhn. "Distributional Dynamics."
> CEPR Discussion Paper 19829, 2026.

Please cite the paper when using the data ([`CITATION.cff`](CITATION.cff)).
Methodology, estimation code, and the full replication package live in
[**DD_replication**](https://github.com/LuisCald/DD_replication).

---

## The data (`data/`)

| File | What it is |
|---|---|
| `smoothed_factors.csv` | **Latent states** $\hat F_t$ (43 columns × 248 quarters, point estimate / posterior mode). **The first 8 columns (`x1`–`x8`) are the distributional factors — the main objects.** `x9`–`x32` are their lags (the model conditions annual surveys on 4-quarter averages), `x33`–`x43` the aggregate block. |
| `smoothed_factor_draws.csv` | **400 posterior draws** of the states (parameter + state uncertainty), same column layout — again, `x1`–`x8` are the factors that matter. The basis for credible bands on any moment. |
| `smoothed_factors_bands.csv` | Ready-made per-quarter 5/50/95 percentiles of each factor. |
| `{PSID,SCF,CEX}_coefficients_{normal,average}.csv` | **Model-implied Legendre coefficients** of the joint distribution, dense across all quarters. `_normal` = HP trend re-added (use in-sample); `_average` = time-averaged trend (use for extrapolation). |
| `{PSID,SCF}_functional_data{,_detrended}.csv` | The raw survey-based coefficient estimates that *enter* the model (NaN where no survey). |
| `PSID_synthetic_microdata.csv` | **Synthetic micro data**: one weighted cross-section per quarter (1000 rows = 10³ decile-copula cells) — cell weight, decile-average consumption/income/wealth (relative to per-HH mean), cell indices. |
| `aggregate_anchors.csv` | Per-household aggregates (`consum/income/wealth_per_hh`, `tot_hhs`) — multiply relative values by these for dollar levels. |

**Units:** distributional values are relative to the quarter's per-household
mean (e.g. 0.80 = 80% of mean consumption). Levels = value ×
`<measure>_per_hh` from `aggregate_anchors.csv`.

## The scripts (`julia/`, `python/`)

Named for the pipeline stage they implement — same API in both languages,
agreeing to machine precision:

```
factors  →  coefficients  →  moments / micro data
```

| Script | In → Out |
|---|---|
| `factors_to_coefficients` | factors → coefficient rows (`FactorMap`; perturb a factor for counterfactuals). **`_average` trend convention only** — the map is re-estimated from the published files (exact, R² = 1.0) since the model's own loadings/means/trend are not shipped; the `_normal` HP trend is not a function of the factors. |
| `coefficients_to_moments` | coefficients → decile cut points + copula density |
| `coefficients_to_micro_data` | coefficients → weighted synthetic cross-sections |
| `posterior_bands` | factor draws → credible band on any moment |
| `plot_factor_bands` | factor draws → factor paths with bands |
| `reconstruct` | the core library both languages share (`Reconstruction`, `quantile_at`, `copula_density_at`, `copula_pmf_grid`, `FactorMap`, …) |

## Quick start

**Julia** (deps pinned in `Project.toml`):

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # once
julia --project=. julia/coefficients_to_moments.jl    # deciles + copula density
DATE=2008-Q3 julia --project=. julia/coefficients_to_micro_data.jl
julia --project=. julia/plot_factor_bands.jl          # factors with 5–95% bands
```

**Python** (needs `numpy`, `pandas`, `matplotlib`):

```bash
python python/coefficients_to_moments.py
python python/posterior_bands.py --measure wealth --date 2008-Q3
```

**As a library:**

```python
import sys; sys.path.insert(0, "python")
from reconstruct import Reconstruction, FactorMap

r = Reconstruction("data/PSID_coefficients_normal.csv")
r.quantile_at("2008-Q3", "consum", [0.1, 0.5, 0.9])   # → [0.325, 0.802, 1.874]
r.copula_density_at("2008-Q3", 0.5, 0.5, 0.5)          # → 1.754
r.copula_pmf_grid("2008-Q3")                           # 10³ cell masses
```

```julia
include("julia/reconstruct.jl"); using .DistributionalReconstruction
r = Reconstruction("data/PSID_coefficients_normal.csv")
quantile_at(r, "2008-Q3", :consum, [0.1, 0.5, 0.9])
```

## Posterior uncertainty

`smoothed_factors.csv` is the posterior mode. For bands, push each row-group of
`smoothed_factor_draws.csv` through `FactorMap` — `posterior_bands` does this
for you. Each draw combines parameter uncertainty (posterior θ) with the
Kalman smoother's state uncertainty; draws are marginal per quarter (exact for
pointwise bands and per-date moments; joint cross-date statistics would need a
joint simulation smoother).

## Technical notes

- Marginal coefficients are stored on an asinh scale; the helpers apply the
  sinh back-transform, so returned values are already in natural (relative)
  units. The copula block is stored as-is.
- Basis: orthonormal Legendre on [0,1], $Q_o(u) = \sqrt{2o+1}\,P_o(2u-1)$;
  12 orders per margin, 12³ copula tensor (34 immutable entries reconstructed
  by the helpers).
- The synthetic micro data reproduces the replication pipeline's
  `construct_micro_dataset` layout; the degree-11 basis truncates the extreme
  wealth tail (weighted mean ≈ 0.93 of the true mean for wealth — the
  Forbes-type top is not representable in this basis).

## License

See [`LICENSE`](LICENSE). Data are model-generated estimates; underlying
surveys (PSID, SCF, CEX, CPS, SIPP) are public — see DD_replication for
sources and the full data availability statement.
