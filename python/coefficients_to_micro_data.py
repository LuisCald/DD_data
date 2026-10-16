#!/usr/bin/env python3
"""
coefficients_to_micro_data.py — synthetic micro data from the published
coefficient files. Python twin of coefficients_to_micro_data.jl; see that
file's header for the full documentation (row layout, how it follows
CreateTimeSeries.jl / construct_micro_dataset, and units).

Short version: one weighted synthetic cross-section per quarter — one row per
cell of the 10x10x10 decile copula grid. `cop_share` is the cell's household
weight share; consum/income/wealth are the MEANS OF EACH MEASURE AMONG THE
HOUSEHOLDS IN THE CELL (copula-weighted, `--values cell`, default since Oct 2026),
relative to the per-household mean. `--values decile` gives the earlier layout,
where each cell carries the unconditional decile average of each measure. Scale to dollar levels with
data/aggregate_anchors.csv (level = value * <measure>_per_hh).
(The published data/PSID_synthetic_microdata.csv and dd_master are already
scaled this way, i.e. in 2019 USD per household.)

    python coefficients_to_micro_data.py                        # PSID, all quarters
    python coefficients_to_micro_data.py --dataset PSID --date 2008-Q3
    python coefficients_to_micro_data.py --out my_micro.csv     # write elsewhere

Default output: python/coefficients_to_micro_data_output.csv (git-ignored;
never overwrites the published data/ files).

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
from reconstruct import Reconstruction, Q, I, GRID_COP  # noqa: E402

MEASURES = ["consum", "income", "wealth"]   # alphabetical — model convention
GRID = 10
NODES_PER_DECILE = 32


def decile_averages(r: Reconstruction, date: str, measure: str) -> np.ndarray:
    """Mean of the quantile function within each decile (mirrors the quadgk
    decile integral in CreateTimeSeries.jl; exact for the polynomial basis)."""
    avgs = np.zeros(GRID)
    for d in range(GRID):
        us = d / GRID + (np.arange(NODES_PER_DECILE) + 0.5) / (GRID * NODES_PER_DECILE)
        avgs[d] = float(np.mean(r.quantile_at(date, measure, us)))
    return avgs


# Gauss–Legendre nodes per decile for the cell means (degree-11 quantile
# polynomial under sinh and degree-11 copula basis: 64 nodes is ample)
_GL_X, _GL_W = np.polynomial.legendre.leggauss(64)
_EDGES = np.arange(GRID + 1) / GRID
# below this copula mass a cell's conditional mean is numerically unstable (ratio of
# near-zero integrals); it carries the decile average instead (≈12% of cells, ≈0.1% of weight)
MIN_CELL_MASS = 5e-5


def cell_means(r: Reconstruction, date: str) -> dict:
    """Mean of each measure among the households in each 10x10x10 cell.

    E[m | u_c in U_i, u_y in U_j, u_w in U_k] = (1/w_ijk) * integral over the cell of
    sinh(Xi^-1_m(u_m)) * c(u) du, with w_ijk the cell's copula mass. The copula
    density is a Legendre polynomial, so this factors into one-dimensional
    integrals per decile (closed form on the conditioning side, Gauss–Legendre on
    the outcome side) — the cell version of RecessionPanels.box_conditional_means.
    A quarter without the consumption marginal conditions on (income, wealth)
    only. Cells with (near-)zero mass, and values a truncated copula pushes outside
    their own decile, fall back to the decile average / are clipped at the decile bounds.
    Returns {measure: (10,10,10) array or None}.
    """
    kap = r._kappa(date)
    U = (_EDGES[:-1, None] + _EDGES[1:, None]) / 2 + _GL_X[None, :] / (2 * GRID)   # (10, n)
    Wq = np.broadcast_to(_GL_W / (2 * GRID), U.shape)
    Qn = np.stack([Q(o, U) for o in range(GRID_COP)], 0)                            # (12, 10, n)
    Iint = np.stack([I(o, _EDGES[1:]) - I(o, _EDGES[:-1]) for o in range(GRID_COP)], 0)  # (12, 10)
    has = {m: not np.all(np.isnan(r._xi(date, m))) for m in MEASURES}
    if not has["consum"]:
        kap = kap.copy(); kap[1:] = 0.0     # only the (income, wealth) slice is defined
    # consumption unobserved: integrate it out (∫_0^1 Q_o = 1{o=0}) in every cell
    Ic = Iint if has["consum"] else np.tile(np.eye(GRID_COP)[:, :1], (1, GRID))
    side = {"consum": Ic, "income": Iint, "wealth": Iint}
    mass = np.einsum("abc,ai,bj,ck->ijk", kap, side["consum"], side["income"], side["wealth"])
    out = {}
    for mi, m in enumerate(MEASURES):
        if not has[m]:
            out[m] = None
            continue
        f = r.quantile_at(date, m, U.ravel()).reshape(U.shape)                       # relative to mean
        dec = (f * Wq).sum(1) * GRID                                                  # decile averages
        F = np.einsum("aGn,Gn,Gn->aG", Qn, f, Wq)                                     # ∫_{U_d} f Q_a
        mats = [side["consum"], side["income"], side["wealth"]]
        mats[mi] = F
        num = np.einsum("abc,ai,bj,ck->ijk", kap, *mats)
        ok = mass > MIN_CELL_MASS
        v = np.where(ok, num / np.where(ok, mass, 1.0), np.nan)
        # range the quantile function spans within each decile (the fitted function
        # need not be monotone near an atom, e.g. zero wealth)
        lo, hi = f.min(1), f.max(1)
        sh = [1, 1, 1]; sh[mi] = GRID
        v = np.clip(v, lo.reshape(sh), hi.reshape(sh))
        out[m] = np.where(np.isnan(v), np.broadcast_to(dec.reshape(sh), v.shape), v)
    return out


def micro_rows(r: Reconstruction, date: str, digits: int = 6, values: str = "cell") -> pd.DataFrame:
    pmf = r.copula_pmf_grid(date)                      # 10x10x10, matches ciw_*
    if values == "cell":
        cm = cell_means(r, date)
        val = lambda m, i, j, k: np.nan if cm[m] is None else round(float(cm[m][i, j, k]), digits)
    else:
        q = {m: np.round(decile_averages(r, date, m), digits) for m in MEASURES}
        idx = {"consum": 0, "income": 1, "wealth": 2}
        val = lambda m, i, j, k: q[m][(i, j, k)[idx[m]]]
    rows = []
    for i in range(1, GRID + 1):                        # same order as construct_micro_dataset
        for j in range(1, GRID + 1):
            for k in range(1, GRID + 1):
                rows.append((round(float(pmf[i - 1, j - 1, k - 1]), digits),
                             int(f"{i}{j}{k}"),
                             *[val(m, i - 1, j - 1, k - 1) for m in MEASURES],
                             i, j, k, date))
    return pd.DataFrame(rows, columns=["cop_share", "grid_point", *MEASURES,
                                       "consumgrid", "incomegrid", "wealthgrid", "time"])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", default="PSID", choices=["PSID", "SCF", "CEX"])
    ap.add_argument("--trend", default="normal", choices=["normal", "average"])
    ap.add_argument("--date", default=None, help="one quarter, e.g. 2008-Q3; omit for all")
    ap.add_argument("--values", default="cell", choices=["cell", "decile"],
                    help="cell = mean among the cell's households (default); decile = unconditional decile averages")
    ap.add_argument("--out", default=os.path.join(HERE, "coefficients_to_micro_data_output.csv"),
                    help="output CSV path (default: python/coefficients_to_micro_data_output.csv)")
    args = ap.parse_args()

    csv = os.path.join(REPO, "data",
                       f"{args.dataset}_coefficients_{args.trend}.csv")
    if not os.path.exists(csv):
        sys.exit(f"Coefficient file not found: {csv}")
    r = Reconstruction(csv)

    dates = [args.date] if args.date else r.available_dates()
    out = pd.concat([micro_rows(r, d, values=args.values) for d in dates], ignore_index=True)

    dest = args.out
    out.to_csv(dest, index=False)
    print(f"Wrote {len(out)} rows ({len(dates)} quarters x {GRID**3} cells) -> {dest}")


if __name__ == "__main__":
    main()
