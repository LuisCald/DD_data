# Codebook

Dollar values are 2019Q4 USD (CPI-U) per household unless noted. `time` is `YYYY-Qq`.

## dd_tidy_long.csv

| Variable | Type | Description |
|---|---|---|
| `time`, `year`, `quarter` | str, int, int | Quarter |
| `survey` | str | `psid` or `scf` (survey view) |
| `variable` | str | Outcome: `consum`, `income`, `wealth` |
| `ranked_by` | str | Variable households are ranked by (equal to `variable` for own-ranked groups) |
| `grouping` | str | `deciles`, `quintiles`, `wealth_groups` |
| `group` | str | `d1`…`d10`; `0-20`…`80-100`; `bottom50`, `next40`, `top10` |
| `mean` | float | Group mean per household |
| `share` | float | Group share of the total; shares in a grouping sum to 1 |

## dd_series.{csv,dta}

| Variable | Type | Description |
|---|---|---|
| `year`, `quarter`, `time` | int, int, str | Quarter |
| `<var>_quantile_<g>_<survey>[_dt]` | float | Mean per household of decile group `g` (`p10` = 0–10th percentile … `p90`, `top` = top 10%); `_dt` = detrended |
| `<var>_level_<g>_<survey>[_dt]` | float | Group total (mean × households / 10), USD |
| `<var>_share_<g>_<survey>[_dt]` | float | Group share of the total; the ten shares sum to 1 |
| `factor_1`…`factor_43` | float | Smoothed states (= `x1`…`x43` in `smoothed_factors.csv`) |

`var` is `consum`, `income` or `wealth`; `survey` is `psid` or `scf` (no SCF consumption).

## dd_microdata.{csv,dta}, implicates/PSID_synthetic_microdata_imp{1..5}.csv

| Variable | Type | Description |
|---|---|---|
| `year`, `quarter`, `time` | int, int, str | Quarter (1984-Q2 to 2024-Q1) |
| `implicate`, `draw` | int | Implicate 1–5 and the posterior draw it uses (implicate files only) |
| `cop_share` | float | Cell population weight; sums to 1 per quarter |
| `grid_point` | int | Cell id: consumption, income, wealth decile concatenated (`372` = 3, 7, 2) |
| `consum`, `income`, `wealth` | float | Mean among the cell's households (`consum` empty before 1999-Q2) |
| `consumgrid`, `incomegrid`, `wealthgrid` | int | Decile of the cell, 1–10 |

## items/dd_items_SCF.{csv,dta} and items/<item>_SCF.csv

| Variable | Type | Description |
|---|---|---|
| `time` (+ `year`, `quarter` in .dta) | str | Quarter |
| `item` | str | `real_estate`, `mortgages`, `consumer_credit`, `pension`, `stocks` (combined file only) |
| `grouping` | str | `wealth_groups`, `income_quintiles`, `own_deciles` |
| `group` | str | `bottom50/next40/top10`; `0-20`…`80-100`; `d1`…`d10` |
| `mean_per_hh` | float | Group mean holding per household |
| `share` | float | Group share of the item total |

## Model files

| File | Variables |
|---|---|
| `smoothed_factors.csv` | `time`, `x1`…`x43` (1–8 distributional factors, 9–32 their lags 1–3, 33–43 aggregate factors) |
| `smoothed_factor_draws.csv` | `draw` (1–200), `time`, `x1`…`x43` |
| `smoothed_factors_bands.csv` | `time`, `xK_p05`, `xK_p50`, `xK_p95` |
| `{PSID,SCF,CEX,CPS,CPS2,SIPP1,SIPP2,SIPP3}_coefficients_{normal,average}.csv` | `time`, `x1`…`x1730`: 3 × 12 marginal coefficients (asinh scale), then 12³ − 34 copula coefficients |
| `{PSID,SCF}_functional_data{,_detrended}.csv` | `time`; `ciw_ijk` copula cell masses; `{quantiles,levels,shares}<var>_<u>` decile profiles relative to the per-household mean (`u` = upper edge, `0.999999` = top) |
| `aggregate_anchors.csv` | `time` (quarter start date), `consum_per_hh`, `income_per_hh`, `wealth_per_hh`, `tot_hhs` (number of households) |
