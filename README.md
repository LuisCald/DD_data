# Distributional Dynamics — Data

Quarterly estimates of the joint distribution of consumption, income and wealth for U.S.
households, 1962-Q3 to 2024-Q1, combining five household surveys (PSID, SCF, CEX, SIPP, CPS)
with aggregate data. From Bayer, Calderon and Kuhn, *Distributional Dynamics*, CEPR DP 19829.

**Vintage:** October 2026 ([CHANGELOG](CHANGELOG.md)) ·
[distributionaldynamics.com](https://www.distributionaldynamics.com)

## Which file do I need?

| I want … | File |
|---|---|
| Group means and shares in a tidy table (bottom 50 / next 40 / top 10, quintiles, deciles; also consumption or wealth by income group) | `data/dd_tidy_long.csv` |
| Decile-group series in one row per quarter | `data/dd_series.{csv,dta}` |
| Household-level (cell) data to merge or compute my own statistics | `data/dd_microdata.{csv,dta}` |
| Uncertainty for statistics computed from the micro data | `data/implicates/` (5 files) |
| Balance sheets by DFA wealth group: real estate, mortgages, consumer credit, pensions, stocks | `data/items/dd_items_SCF.{csv,dta}` |
| Posterior bands on any moment, or the model objects | factor draws, coefficients ([Advanced](#advanced-files)) |

Column definitions for every file: [`docs/codebook.md`](docs/codebook.md).
Everything in plain English: the [data guide (PDF)](docs/DD_data_guide.pdf).

## Quickstart

```python
import pandas as pd
u = "https://www.distributionaldynamics.com/data/"
t = pd.read_csv(u + "dd_tidy_long.csv")
top10 = t.query("survey == 'scf' and variable == 'wealth' and ranked_by == 'wealth' and group == 'top10'")
top10.plot(x="time", y="share", title="Top-10% wealth share (SCF view)")
```

```r
u <- "https://www.distributionaldynamics.com/data/"
t <- read.csv(paste0(u, "dd_tidy_long.csv"))
c_by_inc <- subset(t, survey == "psid" & variable == "consum" & ranked_by == "income" & grouping == "quintiles")
mic <- read.csv(paste0(u, "dd_microdata.csv"))
with(subset(mic, time == "2019-Q4"), weighted.mean(income, cop_share))
```

```stata
use "https://www.distributionaldynamics.com/data/dd_microdata.dta", clear
collapse (mean) income wealth [aw=cop_share], by(year quarter)
use "https://www.distributionaldynamics.com/data/items/dd_items_SCF.dta", clear
```

## Units, timing and groups

- **Dollars:** 2019Q4 USD, deflated with the CPI-U (FRED `CPIAUCSL`). Consumption and income
  are annual flows per household (measured each quarter); wealth is the stock of net worth.
  Income is pre-tax and includes transfers.
- **Timing:** quarterly, `time` = `YYYY-Qq`. PSID series start when the PSID measures do:
  income 1968, wealth 1984 (1984-Q2), consumption 1999. The SCF view has income and wealth only.
- **Survey views:** `psid` and `scf` are the model's distribution in each survey's measurement
  (the SCF captures the top tail better, so its top shares are higher).
- **Groups:** households are ranked by `ranked_by` within the quarter. `bottom50/next40/top10`
  are the Fed DFA wealth groups; `0-20 … 80-100` quintiles; `d1 … d10` deciles. In
  `dd_series`, decile groups are named `p10 … p90, top` (`p10` = 0–10th percentile).
- **Group means across variables** (e.g. consumption by income quintile) come from a 2D
  maximum-entropy step. That is the recommended method; group means from the micro data
  cells are close but coarser ([`docs/methods.md`](docs/methods.md)).

## Synthetic micro data and implicates

`dd_microdata` has one row per cell of the 10×10×10 grid of consumption, income and wealth
deciles, 1,000 rows per quarter, 1984-Q2 to 2024-Q1. `consum/income/wealth` are means among
the cell's households; `cop_share` is the cell's population weight and sums to 1 each quarter.
Always use the weight. There is no household ID: to merge, rank your own households into
deciles and join on `time, consumgrid, incomegrid, wealthgrid`.

Before 1999-Q2 the PSID has no consumption: the 100 (income, wealth) cells sit in the
`consumgrid == 1` rows, and the other 900 rows have weight 0.

The five implicates repeat the file at five posterior draws. Compute your statistic on each
and combine them as with the SCF implicates. The cells are population quantities, so there is
no within-implicate sampling variance, and the variance is (1 + 1/5) × the between-implicate
variance:

```python
import numpy as np
imp = pd.concat(pd.read_csv(f"{u}implicates/PSID_synthetic_microdata_imp{m}.csv") for m in range(1, 6))
est = imp.query("time == '2019-Q4'").groupby("implicate").apply(lambda d: np.average(d.wealth, weights=d.cop_share))
print(est.mean(), np.sqrt(1.2 * est.var()))   # estimate, standard error
```

Five draws give rough bands; for precise bands use the 200 factor draws (`posterior_bands`).

## Balance-sheet items

`items/dd_items_SCF` gives real estate, mortgages, consumer credit, pensions (from 1983-Q3) and
stocks by wealth group (`bottom50/next40/top10`), income quintile and own decile. Stocks are
published for `next40` and `top10` only. Private business is not published. See
[`docs/items_notes.md`](docs/items_notes.md) for definitions and caveats; compare with the
[Fed DFA](https://www.federalreserve.gov/releases/z1/dataviz/dfa/).

## Advanced files

| File | Contents |
|---|---|
| `smoothed_factors.csv` | Posterior mode of the 43 states: `x1`–`x8` distributional factors, `x9`–`x32` their lags 1–3, `x33`–`x43` the 11 aggregate factors. The factors are latent and have no economic labels. |
| `smoothed_factor_draws.csv`, `smoothed_factors_bands.csv` | 200 posterior draws of the state path (`draw, time, x1…x43`) and 5/50/95% bands. |
| `{PSID,SCF,CEX,CPS,CPS2,SIPP1,SIPP2,SIPP3}_coefficients_{normal,average}.csv` | Legendre coefficients, 1,730 per quarter, for each survey view. `_normal` includes the HP trend (use in sample); `_average` uses the time-averaged trend (use for counterfactuals and with `FactorMap`). Blocks a survey does not measure are empty: the CPS views have income only (CPS2 = income years 1967–1992, CPS = 1993 on); the SIPP views income and wealth (SIPP1 = panels 1984–1993, SIPP2 = panels 1996–2008, SIPP3 = the 2014 panel and later releases); only the PSID has the 3D copula. |
| `{PSID,SCF}_functional_data{,_detrended}.csv` | The coefficients evaluated: copula cell masses `ciw_ijk` and decile profiles. |
| `aggregate_anchors.csv` | Per-household consumption, income and wealth and the number of households (`tot_hhs`), 1962–2024. Dates are quarter starts (`1962-07-01` = 1962-Q3). |
| `PSID_synthetic_microdata.csv` | `dd_microdata` without `year`/`quarter` (the layout the scripts write). |

**Scripts** (`python/`, `julia/`; same API in both): `reconstruct` (core library),
`factors_to_coefficients`, `coefficients_to_moments`, `coefficients_to_micro_data`,
`posterior_bands` / `plot_factor_bands`, `cross_conditional_means`. `tools_make_stata.py` and
`tools_make_tidy.py` rebuild the `dd_*` files.

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'                       # once
python python/posterior_bands.py --measure wealth --date 2008-Q3           # numpy, pandas, matplotlib
python python/cross_conditional_means.py --dataset SCF --date all --groups 20,20,20,20,20
```

Method details: [`docs/methods.md`](docs/methods.md).

## Glossary

- **Factor:** a latent time series that drives the whole distribution; 8 distributional and 11 aggregate factors.
- **Legendre coefficients:** the numbers that describe each quarter's marginal distributions and copula.
- **Copula:** the joint ranking of households in consumption, income and wealth, separate from the dollar amounts.
- **Survey view:** the distribution as each survey would measure it.
- **Implicate:** one complete version of the micro data drawn from the posterior.
- **Maximum entropy:** the method used for group means of one variable by another.

## Cite, license, contact

```bibtex
@techreport{BayerCalderonKuhn2026,
  author      = {Bayer, Christian and Calderon, Luis and Kuhn, Moritz},
  title       = {Distributional Dynamics},
  institution = {Centre for Economic Policy Research},
  type        = {CEPR Discussion Paper},
  number      = {19829},
  year        = {2026},
  month       = mar
}
```

Data and documentation: [CC BY 4.0](LICENSE) (use freely, with attribution: cite the paper).
Code: MIT ([`LICENSE-CODE`](LICENSE-CODE)). The data are model-generated estimates. Questions and
problems: [GitHub issues](https://github.com/LuisCald/DD_data/issues) or
luis.calderon@uni-bonn.de. Estimation code, sources and the processed survey data are in
[DD_replication](https://github.com/LuisCald/DD_replication).
