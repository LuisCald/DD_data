# Balance-sheet items: notes and caveats

Files: `data/items/dd_items_SCF.{csv,dta}` (all items) and
`data/items/{real_estate,mortgages,consumer_credit,pension,stocks}_SCF.csv` (one item each).

## What they are

The distribution of five balance-sheet items from the model's SCF reconstruction.
Each item is estimated jointly with income and wealth, using the same specification
as the main data. Coverage is 1962-Q3 to 2024-Q1; pensions start in 1983-Q3.
Stocks are published for the `next40` and `top10` wealth groups only (see Caveats).

Columns: `time, item, grouping, group, mean_per_hh, share`. The per-item CSVs have no
`item` column. The `.dta` file also has `year` and `quarter`.

Groupings:

- `wealth_groups`: households ranked by net worth: `bottom50`, `next40`, `top10`
  (the grouping used in the Fed's Distributional Financial Accounts, DFA).
- `income_quintiles`: households ranked by income: `0-20` … `80-100`.
- `own_deciles`: households ranked by the item itself: `d1` … `d10`.

`mean_per_hh` is the group's average holding per household in 2019Q4 dollars (CPI-U). Item
totals are anchored to the Financial Accounts. `share` is the group's share of the
item total.

## Method

Zero holdings are modelled explicitly, as a point mass at zero plus the holders'
distribution. Group means of one variable by another (for example, real estate by
wealth group) come from the 2D maximum-entropy step (see [`methods.md`](methods.md)).

## Caveats

- Real estate includes non-owner-occupied property. The DFA books part of this under
  business equity.
- The pension total is anchored to DB+DC entitlements, but the survey item is DC/IRA
  balances.
- Between SCF waves the series are model-based interpolations. This matters most before
  1989, when there are few waves. The series step when a new wave enters.
- Small negative group means (deciles inside the mass of zero holders, at most a few hundred
  dollars) are set to 0, and those shares renormalized.
- Stocks: only the `next40` and `top10` wealth groups are published. Bottom-50% stock
  holdings are below the resolution of the representation, and the shares are of the
  model's item total. Private business is not published.
