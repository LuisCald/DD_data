# Changelog

## 2026-10 (October 2026 vintage)

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
