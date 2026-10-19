"""Build data/CPS_functional_data.csv and data/CPS1967_functional_data.csv: the CPS survey views.

The model reports the estimated income distribution as each survey would measure it. The CPS changed
how it measures income in 1993 (new interview method, higher top codes), so the model treats the two
eras as two surveys, with two published coefficient files:

  cps       CPS_coefficients_normal.csv,  income years 1993 on   -> 1993-Q4 .. 2024-Q1
  cps1967   CPS2_coefficients_normal.csv, income years 1967-1992 -> 1967-Q4 .. 1993-Q3

Each view is kept on its own survey's coverage (the CPS2 file continues after 1993 only as a model
extrapolation of the old measurement, which is not published). Income only: the CPS has no
consumption or wealth, so there are no cross-tabulations.

Columns follow {PSID,SCF}_functional_data.csv: `quantilesincome_<u>` = decile-group mean per household
(2019Q4 USD), `levelsincome_<u>` = group total (mean x households / 10), `sharesincome_<u>` = share of
the total, for decile groups ending at u = 0.1 .. 0.9 and 0.999999 (top 10%). Decile means come from
the coefficients exactly as in python/coefficients_to_moments.py; on the PSID view this reproduces the
published shares to within 3e-4.

Usage: python3 tools_make_cps_views.py   (then tools_make_stata.py and tools_make_tidy.py)
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "python"))
from reconstruct import Reconstruction  # noqa: E402
from coefficients_to_moments import decile_averages, has_measure  # noqa: E402

VIEWS = [("CPS", "CPS", "1993-Q4", "2024-Q1"),
         ("CPS1967", "CPS2", "1967-Q4", "1993-Q3")]
EDGES = [f"{u / 10:.1f}" for u in range(1, 10)] + ["0.999999"]


def quarter_key(t):
    y, q = t.split("-Q")
    return int(y) * 4 + int(q)


def anchors():
    a = pd.read_csv(os.path.join(HERE, "data", "aggregate_anchors.csv"))
    ts = pd.to_datetime(a.time)
    a["time"] = ts.dt.year.astype(str) + "-Q" + ts.dt.quarter.astype(str)
    return a.set_index("time")[["income_per_hh", "tot_hhs"]]


def build(out_name, coef_name, first, last, anc):
    r = Reconstruction(os.path.join(HERE, "data", f"{coef_name}_coefficients_normal.csv"))
    lo, hi = quarter_key(first), quarter_key(last)
    rows = []
    for t in r.available_dates():
        if not (lo <= quarter_key(t) <= hi) or not has_measure(r, t, "income"):
            continue
        rel = decile_averages(r, t, "income")              # decile means relative to the per-HH mean
        if not np.all(np.isfinite(rel)):
            continue
        mean = rel * anc.loc[t, "income_per_hh"]
        row = {"time": t}
        row.update({f"quantilesincome_{e}": v for e, v in zip(EDGES, mean)})
        row.update({f"levelsincome_{e}": v * anc.loc[t, "tot_hhs"] / 10 for e, v in zip(EDGES, mean)})
        row.update({f"sharesincome_{e}": v for e, v in zip(EDGES, rel / rel.sum())})
        rows.append(row)
    df = pd.DataFrame(rows)
    path = os.path.join(HERE, "data", f"{out_name}_functional_data.csv")
    df.to_csv(path, index=False)
    print(f"wrote {path}: {len(df)} quarters, {df.time.iloc[0]}..{df.time.iloc[-1]}")


def main():
    anc = anchors()
    for out_name, coef_name, first, last in VIEWS:
        build(out_name, coef_name, first, last, anc)


if __name__ == "__main__":
    main()
