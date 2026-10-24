# Methods and technical notes

Full methodology: Bayer, Calderon, Kuhn, "Distributional Dynamics", CEPR DP 19829
(2026). Estimation code: [DD_replication](https://github.com/LuisCald/DD_replication).

## Specification (October 2026 vintage)

These are the estimates from the paper's replication specification: 8 distributional
factors and 11 aggregate factors, with the full innovation covariance in the likelihood.
That covariance includes the contemporaneous correlation between distributional and
aggregate innovations. The 24 October release re-estimates the model with the tight
measurement-error prior on the aggregate factors (which an earlier release had not applied),
the SCF and PSID income rows placed at the fourth quarter of their reference year, and 1,600
sampler iterations; see the CHANGELOG for what moved.

`smoothed_factors.csv` has 43 state columns:

- `x1`–`x8`: the distributional factors.
- `x9`–`x32`: lags 1–3 of `x1`–`x8`. For example, `x9` is `x1` lagged once and `x17`
  is `x1` lagged twice. The lags are there because the model conditions annual surveys
  on 4-quarter averages.
- `x33`–`x43`: the 11 aggregate factors.

## Representation

- Basis: orthonormal Legendre polynomials on [0,1], $Q_o(u) = \sqrt{2o+1}\,P_o(2u-1)$.
  There are 12 orders per marginal and a 12³ copula tensor. The 34 immutable copula
  entries are not stored; the helpers reconstruct them.
- Each coefficient file has 1730 columns: 3 × 12 marginal coefficients plus
  12³ − 34 copula coefficients.
- Marginal coefficients are stored on an asinh scale, and the helpers apply the sinh
  back-transform. Returned values are therefore already in natural units, relative to
  the per-household mean. The copula block is stored as is.
- `_normal` files add the HP trend back (use these in sample). `_average` files use a
  time-averaged trend (use these for extrapolation and with `FactorMap`).
- `FactorMap` re-estimates the factors → coefficients map from the published files.
  The fit is exact (R² = 1) under the `_average` convention. The `_normal` HP trend is
  not a function of the factors.

## Posterior uncertainty

`smoothed_factors.csv` is the posterior mode. `smoothed_factor_draws.csv` holds 200
draws. Each draw combines parameter uncertainty (a posterior θ draw) with state
uncertainty: given that θ, the whole state path is drawn from p(x_1, …, x_T | data) by
forward-filtering backward-sampling. Draws are therefore coherent over time, so changes
between dates within one draw carry the right uncertainty. To get a band on any moment,
push each draw through `FactorMap` and compute the moment; `posterior_bands` does this
for you.

## Synthetic micro data

`dd_microdata` / `PSID_synthetic_microdata.csv` follow the replication pipeline's
`construct_micro_dataset` layout: one row per cell of the 10×10×10 decile copula grid.
The weight `cop_share` is the copula mass of the cell. Each of `consum`, `income` and
`wealth` is the mean of that variable among the households in the cell,
E[m | u_c ∈ U_i, u_y ∈ U_j, u_w ∈ U_k], integrated from the copula density and the
quantile functions — the cell version of the group means in the paper's recession figure.
Two safeguards: cells with copula mass below 5e-5 (about 12% of cells, 0.1% of the
weight) carry the decile average, because their conditional mean is a ratio of
near-zero integrals; and values the truncated copula pushes outside their own decile are
clipped to that decile's range. In 1984–1999 the PSID has no consumption; the file then
conditions on income and wealth only. `coefficients_to_micro_data` reproduces the file
(`--values decile` gives the earlier layout, where each cell carried the unconditional
decile average). In those quarters the (income, wealth) cells sit in the rows with
`consumgrid == 1`; the nine other copies of each cell carry weight 0 and `consum` is empty.

The five `implicates/` files repeat this at five posterior draws (weights: the mode
weights plus that draw's change in cell mass). The 12-term (degree-11) basis truncates the extreme
wealth tail, so the weighted mean of wealth is typically 0.94–0.98 of the per-household
aggregate (median 0.96). The Forbes-type top is not representable in this basis.

## Group means across variables

`cross_conditional_means` computes group means of one variable within groups of another
from a 2D maximum-entropy density fitted to the copula coefficients. This avoids the bias
of the truncated Legendre series for heavy-tailed outcomes. `METHOD=series` gives the
series version for comparison. This is the recommended method for group means across
variables, and it is what `dd_tidy_long.csv` and the balance-sheet items use. Group means
computed from the synthetic micro data cells are close but not identical, because the
cells are coarser (deciles) and the cell values inherit the truncated series.
