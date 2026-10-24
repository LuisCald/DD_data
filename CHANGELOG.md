# Changelog

## 2026-10-24 (balance-sheet items and products re-estimated)

- Fixes a bug in the measurement equation (DD_replication 7108924): in the annual surveys (SCF, PSID, CPS) only
  consumption, income and wealth had measurement rows, so any other measure in a run (the balance-sheet items and
  the liquidity and debt measures of the products) was never linked to the model's state in those surveys. Its
  distribution followed a trend through the survey waves, with quarterly movement only through the aggregate
  anchor and the copula. Balance-sheet measures now load on the current state as point-in-time stocks.
- The five items (`data/items/`: real estate, mortgages, consumer credit, pensions, stocks) and the four products
  (`data/products/`) are re-estimated with the fix and the settings of the main data (tight aggregate
  measurement-error prior, income rows at Q4 of the reference year, 1,600 sampler iterations).
- Item wealth-group shares move by at most 2–5 points, mostly between SCF waves (the 1970s gap); survey-year values
  change by less than half a point except pensions after 2023 (2023-Q4 top-10% share 49% vs 53%) and mortgages and
  real estate in 2023-24 (1–2 points). Quarter-to-quarter variation rises for most items. Product shares of
  households with no debt move by less than one point.

## 2026-10-24 (October 2026 vintage, corrected aggregate priors)

- Re-estimated after fixing the specification of the aggregate priors in the measurement equation (the tight
  measurement-error prior on the aggregate factors was not being applied) and with the SCF and PSID income
  rows placed at the fourth quarter of their reference year; 1,600 sampler iterations (400 before). All CIW
  files rebuilt from the new run: coefficients, functional data, factors, posterior draws and bands, micro data
  and implicates, `dd_series`, `dd_tidy_long`, CPS views. Published decile-group shares move by at most 1.2
  points (SCF-view wealth), 0.8 points (income) and 0.1 points (consumption); the 2022-Q3 SCF-view top-10%
  wealth share is 72.9% (73.0% before). The factors are re-identified (correlation with the previous vintage
  0.73–0.86 per factor), so `smoothed_factors`, the draws and `FactorMap` fits are not comparable across vintages.
- Balance-sheet items and products: re-estimated in the entry above.

## 2026-10-19 (corrected October 2026 vintage)

- Re-estimated with the corrected SCF 2022 wave (DD_replication issue #9): retirement accounts removed from
  2022 wealth (as in 1950–2019), CPI-U factors for 2022, and 2022 weights on the SCF+ scale (the Forbes 400
  now carry their weight in 2022). All CIW files rebuilt from the new run: coefficients, functional data,
  factors, posterior draws and bands, micro data and implicates, `dd_series`, `dd_tidy_long`, CPS views.
  Main change: SCF-view wealth from 2019 on (2022-Q3 top-10% share 73.0% vs 70.7% before); other series move
  by less than 0.4 points. Cross-tabulations (maximum entropy) move by up to about 3 points in some quarters,
  including early ones, from the re-estimation itself.
- New `data/products/`: four further joint distributions (consumption with gross liquid assets, total assets,
  consumer, total and mortgage debt), built by `tools_make_products.py`.
- Balance-sheet items re-estimated on the corrected SCF 2022 wave (2026-10-20). Through 2019 they match the
  earlier files; the changes come after 2020. Main change: pensions by wealth group (2023-Q4 top-10% share
  53% vs 70% before, which had jumped from 51% in 2019). The stocks next-40 share stays near 9% (12% before
  in 2023–24). Own-decile series are unchanged.

## 2026-10 (October 2026 vintage)

- CPS survey views (income only) in `dd_series` and `dd_tidy_long`: `cps` (1993-Q4 on) and `cps1967`
  (1967-Q4 to 1993-Q3), built from the published CPS and CPS2 coefficients by `tools_make_cps_views.py`.
  Existing PSID and SCF rows are unchanged.

- Estimates from the paper's replication specification: 8 distributional and 11 aggregate factors.
- Synthetic micro data carry the means among each cell's households; five implicates from joint posterior draws.
- New `dd_tidy_long.csv`: group means and shares, including groups across variables (maximum entropy).
- Balance-sheet items: real estate, mortgages, consumer credit, pensions, and stocks (next 40 / top 10).
- Cleanup: integer cell keys, weight 0 (not missing) on the empty consumption cells before 1999-Q2,
  `year`/`quarter` in the implicates, empty SCF consumption columns dropped from `dd_series`,
  anchors start in 1962, items use `time` like the other files; small negative item means and
  two near-zero consumption cells set to 0.
- `dd_master` renamed `dd_microdata` (same content).
- License: data CC BY 4.0, code MIT.
- Legendre coefficients for the CPS (CPS, CPS2) and SIPP (SIPP1–3) views; the scripts accept these datasets.
