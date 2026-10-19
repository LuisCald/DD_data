# Changelog

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
- The balance-sheet items are unchanged (estimated on the earlier SCF 2022 wave).

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
