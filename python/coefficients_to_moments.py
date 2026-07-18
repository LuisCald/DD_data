#!/usr/bin/env python3
"""
coefficients_to_moments.py — turn the published coefficient files into
moments: decile cut points and averages, levels, shares, Ginis, and the
copula density at the joint median. Python twin of
coefficients_to_moments.jl (see that header for the definitions).

Per measure and date:
  * decile cut points  Xi^-1(u), u = 0.1 ... 0.9   (relative to per-HH mean)
  * decile averages    32-node quadrature within each decile
  * levels             group average x <measure>_per_hh (aggregate_anchors.csv)
  * shares             bottom 50% / next 40% / top 10% of the total
  * Gini               robust Gini (Raffinetti, Siletti & Vernizzi 2015) on
                       the 10 decile averages, equal weights — handles
                       negative wealth

    python python/coefficients_to_moments.py
    python python/coefficients_to_moments.py --dataset SCF --date 2020-Q3
    python python/coefficients_to_moments.py --date all       # full time series

Output: python/coefficients_to_moments_output.csv, long format
  date, measure, stat, index, value

Dependencies: numpy, pandas.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)
from reconstruct import Reconstruction  # noqa: E402

MEASURES = ["consum", "income", "wealth"]
GRID = 10
NODES_PER_DECILE = 32
GROUPS = [("bot50", slice(0, 5)), ("mid40", slice(5, 9)), ("top10", slice(9, 10))]


def has_measure(r, date, m):
    try:
        r.quantile_at(date, m, [0.5])
        return True
    except Exception:
        return False


def decile_averages(r, date, m):
    avgs = np.zeros(GRID)
    for d in range(GRID):
        us = d / GRID + (np.arange(NODES_PER_DECILE) + 0.5) / (GRID * NODES_PER_DECILE)
        avgs[d] = float(np.mean(r.quantile_at(date, m, us)))
    return avgs


def gini_robust(x, w=None):
    """Generalized (robust) Gini of Raffinetti, Siletti & Vernizzi (2015);
    handles negatives. Ported from Distributional_Counterfactuals."""
    x = np.asarray(x, dtype=float)
    w = np.ones_like(x) if w is None else np.asarray(w, dtype=float)
    t_plus = float(np.sum(w * np.maximum(x, 0)))
    t_minus = float(np.sum(w * -np.minimum(x, 0)))
    num = float(np.sum(np.outer(w, w) * np.abs(np.subtract.outer(x, x))))
    denom = 2 * (w.sum() - w.min()) * (t_plus + t_minus)
    return num / denom


def load_anchors():
    df = pd.read_csv(os.path.join(REPO, "data", "aggregate_anchors.csv"))
    ts = pd.to_datetime(df["time"])
    keys = ts.dt.year.astype(str) + "-Q" + ts.dt.quarter.astype(str)
    return {k: {m: float(row[f"{m}_per_hh"]) for m in MEASURES}
            for k, (_, row) in zip(keys, df.iterrows())}


def moments_for_date(r, anchors, date, verbose=True):
    rows = []
    anchor = anchors.get(date)
    us = np.round(np.arange(0.1, 1.0, 0.1), 2)
    measures = [m for m in MEASURES if has_measure(r, date, m)]

    for m in measures:
        cuts = r.quantile_at(date, m, us)
        avgs = decile_averages(r, date, m)
        tot = avgs.sum()
        g = gini_robust(avgs, np.full(GRID, 0.1))

        rows += [(date, m, "decile_cut", str(u), v) for u, v in zip(us, cuts)]
        rows += [(date, m, "decile_avg", str(d + 1), avgs[d]) for d in range(GRID)]
        for name, ids in GROUPS:
            rows.append((date, m, "share", name, avgs[ids].sum() / tot))
            if anchor is not None:
                rows.append((date, m, "level_per_hh", name, avgs[ids].mean() * anchor[m]))
        rows.append((date, m, "gini", "", g))

        if verbose:
            print(f"{m:>7} deciles:  " + "  ".join(f"{v:6.3f}" for v in avgs))
            shares = " / ".join(f"{avgs[ids].sum() / tot:.3f}" for _, ids in GROUPS)
            print(f"{'':>7} shares (bot50/mid40/top10): {shares}   gini: {g:.3f}")

    if len(measures) == 3:
        cmed = float(r.copula_density_at(date, 0.5, 0.5, 0.5))
        rows.append((date, "joint", "copula_density_median", "", cmed))
        if verbose:
            print(f"\ncopula density at joint median: {cmed:.4f}")
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", default="PSID", choices=["PSID", "SCF", "CEX"])
    ap.add_argument("--trend", default="normal", choices=["normal", "average"])
    ap.add_argument("--date", default="2008-Q3", help="a quarter, or 'all'")
    args = ap.parse_args()

    csv = os.path.join(REPO, "data", f"{args.dataset}_coefficients_{args.trend}.csv")
    if not os.path.exists(csv):
        sys.exit(f"Coefficient file not found: {csv}")
    r = Reconstruction(csv)
    anchors = load_anchors()

    rows = []
    if args.date.lower() == "all":
        for n, d in enumerate(r.available_dates(), 1):
            rows += moments_for_date(r, anchors, d, verbose=False)
            if n % 50 == 0:
                print(f"  {n} quarters done")
    else:
        print(f"Dataset {args.dataset} ({args.trend} trend), date {args.date}\n")
        rows = moments_for_date(r, anchors, args.date)

    out = os.path.join(HERE, "coefficients_to_moments_output.csv")
    pd.DataFrame(rows, columns=["date", "measure", "stat", "index", "value"]).to_csv(out, index=False)
    print(f"\nWrote {len(rows)} rows -> {out}")


if __name__ == "__main__":
    main()
