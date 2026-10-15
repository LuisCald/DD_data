# Distributional Dynamics — Data

Quarterly estimates of the joint distribution of consumption, income, and wealth for
U.S. households, 1962-Q3 to 2024-Q1 (247 quarters), with posterior uncertainty.

**Vintage:** October 2026 · 8 distributional factors · Legendre order 12 · 200 posterior draws ·
[distributionaldynamics.com](https://www.distributionaldynamics.com)

## Load the data in one line

```python
import pandas as pd
df    = pd.read_csv("https://www.distributionaldynamics.com/data/dd_series.csv")
items = pd.read_csv("https://www.distributionaldynamics.com/data/items/dd_items_SCF.csv")
```

```julia
using CSV, DataFrames, Downloads
df = CSV.read(Downloads.download("https://www.distributionaldynamics.com/data/dd_series.csv"), DataFrame)
```

```stata
use "https://www.distributionaldynamics.com/data/dd_series.dta", clear
use "https://www.distributionaldynamics.com/data/items/dd_items_SCF.dta", clear
```

## Files (`data/`)

| File | Contents |
|---|---|
| `dd_series.{csv,dta}` | **Main file.** 247 quarters × 406 columns: decile-group means, totals and shares of consumption, income and wealth in the PSID and SCF views (plus `_dt` detrended versions), and the 43 smoothed states as `factor_1`–`factor_43`. |
| `dd_master.{csv,dta}` | Synthetic micro data, PSID view, 1984-Q2 to 2024-Q1: 1,000 cells (10×10×10 deciles) per quarter. Cell weight `cop_share`, decile-group means `consum/income/wealth` in 2019 USD per household, and cell indices. Copy of `PSID_synthetic_microdata.csv` with `year`/`quarter` added. |
| `items/dd_items_SCF.{csv,dta}` | Balance-sheet items: real estate, mortgages, consumer credit, pensions (pensions from 1983-Q3). There is also one CSV per item. |
| `smoothed_factors.csv` | Posterior mode of the states. `x1`–`x8` are the distributional factors, `x9`–`x32` their lags 1–3, and `x33`–`x43` the 11 aggregate factors. |
| `smoothed_factor_draws.csv`, `smoothed_factors_bands.csv` | 200 posterior draws (`draw, time, x1…x43`) and per-quarter 5/50/95% bands. |
| `{PSID,SCF,CEX}_coefficients_{normal,average}.csv` | Legendre coefficients of the joint distribution, 1,730 per quarter. `_normal` adds the HP trend back (use in sample); `_average` uses the time-averaged trend (use for extrapolation and with `FactorMap`). `{PSID,SCF}_functional_data{,_detrended}.csv` are the same coefficients evaluated: copula cell masses `ciw_ijk` and decile profiles. |
| `aggregate_anchors.csv` | Per-household consumption, income and wealth (2019 USD) and the number of households (`tot_hhs`). Dates are quarter starts (`1962-07-01` = 1962-Q3). The published range starts at 1962-Q3; rows from 1947 to 1962-Q2 are pre-sample, and the file runs to 2024-Q2. |

## Column names

`dd_series` columns follow `<var>_<stat>_<pctl>_<survey>[_dt]`:

- `var` is `consum`, `income` or `wealth`.
- `survey` is `psid` or `scf`. The SCF has no consumption, and the non-detrended PSID
  series start when the PSID measures do: income in 1968, wealth in 1984, consumption in 1999.
- `pctl` is `p10`…`p90` or `top`. Each names a decile group: `p10` is the 0–10th percentile,
  `p90` the 80–90th, and `top` the top 10%.
- `stat` is one of:
  - `quantile`: the group's mean per household, in 2019 USD. (It is a group mean, not a
    point quantile.)
  - `level`: the group's total, in 2019 USD.
  - `share`: the group's share of the total. The ten shares sum to 1.

Items columns are `date, item, grouping, group, mean_per_hh, share`. `mean_per_hh` is in
2019Q4 USD. There are three groupings: `wealth_groups` (bottom50/next40/top10),
`income_quintiles` and `own_deciles`.

## Scripts (`julia/`, `python/`)

The two languages share the same API and agree to machine precision. The pipeline is
factors → coefficients → moments / micro data.

| Script | Does |
|---|---|
| `reconstruct` | Core library (`Reconstruction`, `quantile_at`, `copula_density_at`, `copula_pmf_grid`, `FactorMap`) |
| `factors_to_coefficients` | Factors → coefficients (`_average` convention); perturb a factor for counterfactuals |
| `coefficients_to_moments` | Decile cuts and averages, shares, Ginis, copula density |
| `coefficients_to_micro_data` | Weighted synthetic cross-sections, relative to the per-household mean. Writes `*_output.csv` next to the script; use `--out` / `OUT=` to write elsewhere |
| `posterior_bands`, `plot_factor_bands` | Factor draws → credible bands on a moment or on the factors |
| `cross_conditional_means` | Means of one variable within groups of another (2D maximum entropy) |
| `validate_reconstruction` (Python only), `tools_make_stata.py` | Consistency check against the functional data; rebuilds the `dd_*` files |

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'                      # once
python python/posterior_bands.py --measure wealth --date 2008-Q3          # numpy, pandas, matplotlib
python python/cross_conditional_means.py --dataset SCF --date all --groups 50,40,10
```

## Notes

- **Items:** these are model-based between SCF waves, and stocks and private business are
  not published. See [`docs/items_notes.md`](docs/items_notes.md) for the caveats.
- **Uncertainty:** `smoothed_factors.csv` is the posterior mode. For bands on any moment,
  push the draws through `FactorMap` (`posterior_bands` does this).
- **Method:** see [`docs/methods.md`](docs/methods.md) (basis, factor layout, trends) and
  [`docs/maxent_note.pdf`](docs/maxent_note.pdf) (maximum-entropy step).

## Cite

```bibtex
@techreport{BayerCalderonKuhn2026,
  author      = {Bayer, Christian and Calderon, Luis and Kuhn, Moritz},
  title       = {Distributional Dynamics},
  institution = {Centre for Economic Policy Research},
  type        = {CEPR Discussion Paper},
  number      = {19829},
  year        = {2026}
}
```

MIT license ([`LICENSE`](LICENSE)). The data are model-generated estimates. Estimation code,
sources and the replication package are in [DD_replication](https://github.com/LuisCald/DD_replication).
