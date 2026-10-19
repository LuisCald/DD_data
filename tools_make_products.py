"""Build data/products/<product>/: additional joint distributions (beyond consumption x income x wealth).

Each product is a separate estimation of the model with three measures (October 2026 vintage):

  consum_grossliq_debt      consumption x gross liquid assets x consumer debt
  consum_totassets_totdebt  consumption x total assets x total debt
  consum_grossliq_mortdebt  consumption x gross liquid assets x mortgage debt
  consum_income_mortdebt    consumption x income x mortgage debt

Debt measures have a point mass at zero (households with no debt); the model carries the share of
zeros in its state, so it is estimated every quarter.

Usage: python3 tools_make_products.py <product> <results_dir>
  <results_dir> = the run's from_mcmc/data folder (or a copy of the files listed in INPUTS).

Writes, for the SCF and PSID survey views:
  <product>_tidy_long.csv   group means and shares, same columns as dd_tidy_long.csv:
                            own-ranked deciles, quintiles and bottom 50 / next 40 / top 10, and each measure
                            by quintiles and bottom 50 / next 40 / top 10 of the other two (cross groups)
  <product>_share_zero.csv  share of households with zero holdings of the debt measure, per quarter
  {PSID,SCF}_coefficients_{normal,average}.csv and coefficients_meta.csv (column layout; the debt block's
                            last column is the logit of the share of zeros, not a Legendre coefficient)
"""
import os
import shutil
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SURVEYS = ["PSID", "SCF"]
EDGES = [f"{u / 10:.1f}" for u in range(1, 10)] + ["0.999999"]
DECILES = [f"d{i}" for i in range(1, 11)]
GROUPINGS = {
    "deciles": [(f"d{i + 1}", [i]) for i in range(10)],
    "quintiles": [(f"{20 * i}-{20 * (i + 1)}", [2 * i, 2 * i + 1]) for i in range(5)],
    "wealth_groups": [("bottom50", list(range(5))), ("next40", list(range(5, 9))), ("top10", [9])],
}
CROSS = {"bot50": ("wealth_groups", "bottom50", 0.5), "mid40": ("wealth_groups", "next40", 0.4),
         "top10": ("wealth_groups", "top10", 0.1),
         **{f"q{i + 1}": ("quintiles", f"{20 * i}-{20 * (i + 1)}", 0.2) for i in range(5)}}
NONNEG = {"consum", "grossliq", "debt", "totassets", "totdebt", "mortdebt"}
SUFFIX = "_A non-diag_.csv"


def measures_of(product):
    return product.split("_")


def floor_small_negatives(df, measure, col="mean"):
    """Nonnegative measures, own-ranked groups: set small negative group means to 0. These are bottom groups that lie
    inside the zero mass (or at the bottom of the support), where the reconstruction is a few hundred dollars below
    zero; same rule as the published items. Larger negatives stop the build."""
    if measure not in NONNEG:
        return df
    neg = df[col] < 0
    if (df.loc[neg, col] < -1000).any():
        raise ValueError(f"{measure}: own-ranked group mean below -1000 USD")
    df.loc[neg, col] = 0.0
    return df


def own_ranked(f, survey, measures):
    rows = []
    for m in measures:
        q = [f"quantiles{m}_{e}" for e in EDGES]
        s = [f"shares{m}_{e}" for e in EDGES]
        if q[0] not in f or f[q].isna().all().all():
            continue
        for grouping, groups in GROUPINGS.items():
            for g, idx in groups:
                rows.append(pd.DataFrame({
                    "time": f.time, "survey": survey.lower(), "variable": m, "ranked_by": m,
                    "grouping": grouping, "group": g,
                    "mean": f[[q[i] for i in idx]].mean(axis=1, skipna=False),
                    "share": f[[s[i] for i in idx]].sum(axis=1, min_count=len(idx))}))
    return rows


def cross_groups(c, survey, measures):
    rows = []
    for outcome in measures:
        for cond in measures:
            if cond == outcome:
                continue
            for code, (grouping, g, width) in CROSS.items():
                col = f"{outcome}_by_{cond}_{code}"
                if col not in c:
                    continue
                rows.append(pd.DataFrame({
                    "time": c.time, "survey": survey.lower(), "variable": outcome, "ranked_by": cond,
                    "grouping": grouping, "group": g, "mean": c[col], "_w": width}))
    if not rows:
        return []
    df = pd.concat(rows, ignore_index=True)
    df["share"] = df["mean"] * df["_w"]
    key = ["time", "survey", "variable", "ranked_by", "grouping"]
    df["share"] /= df.groupby(key)["share"].transform("sum")
    return [df.drop(columns="_w")]


def main(product, src):
    measures = measures_of(product)
    out = os.path.join(HERE, "data", "products", product)
    os.makedirs(out, exist_ok=True)
    rows, zeros = [], []
    for s in SURVEYS:
        f = pd.read_csv(os.path.join(src, f"{s}_functional_data{SUFFIX}"))
        rows += own_ranked(f, s, measures)
        cpath = os.path.join(src, f"{s}_cross_conditional_dfa{SUFFIX}")
        if os.path.exists(cpath):
            rows += cross_groups(pd.read_csv(cpath), s, measures)
        for m in measures:
            zpath = os.path.join(src, f"{s}_participation_{m}{SUFFIX}")
            if os.path.exists(zpath):
                z = pd.read_csv(zpath)
                zeros.append(pd.DataFrame({"time": z.time, "survey": s.lower(), "variable": m,
                                           "share_zero": z.share_zero}))
        for kind in ["normal", "average"]:
            shutil.copyfile(os.path.join(src, f"{s}_coefficients_{kind}.csv"),
                            os.path.join(out, f"{s}_coefficients_{kind}.csv"))
    meta = os.path.join(src, "coefficients_meta.csv")
    if os.path.exists(meta):
        shutil.copyfile(meta, os.path.join(out, "coefficients_meta.csv"))

    df = pd.concat(rows, ignore_index=True).dropna(subset=["mean", "share"])
    own = df.variable == df.ranked_by
    for m in measures:
        sel = own & (df.variable == m)
        df.loc[sel] = floor_small_negatives(df.loc[sel].copy(), m)
    # ... and their shares: clip at 0 and rescale so each own-ranked grouping still sums to 1
    sel = own & df.variable.isin(NONNEG)
    df.loc[sel, "share"] = df.loc[sel, "share"].clip(lower=0.0)
    key = ["time", "survey", "variable", "grouping"]
    df.loc[sel, "share"] /= df[sel].groupby(key)["share"].transform("sum")
    # Cross groups: a negative mean of a nonnegative measure is a failed cross-tab (the 2D step's truncated-series
    # fallback), not small noise. Leave those cells out instead of flooring them.
    bad = ~own & df.variable.isin(NONNEG) & (df["mean"] < 0)
    if bad.any():
        print(f"  left out {bad.sum()} cross-group cells with a negative mean: "
              + ", ".join(f"{k}: {v}" for k, v in df[bad].groupby(["survey", "variable", "ranked_by", "group"]).size().items()))
        df = df[~bad]
    yq = df.time.str.extract(r"(\d{4})-Q(\d)").astype(int)
    df.insert(1, "year", yq[0])
    df.insert(2, "quarter", yq[1])
    order = {k: i for i, k in enumerate(["deciles", "quintiles", "wealth_groups"])}
    df = df.sort_values(["survey", "variable", "ranked_by", "grouping", "time"],
                        key=lambda x: x.map(order) if x.name == "grouping" else x, kind="stable")
    df.to_csv(os.path.join(out, f"{product}_tidy_long.csv"), index=False, float_format="%.6g")
    msg = f"{product}: {len(df)} tidy rows, {df.time.min()}..{df.time.max()}"
    if zeros:
        z = pd.concat(zeros, ignore_index=True).dropna(subset=["share_zero"])
        z.to_csv(os.path.join(out, f"{product}_share_zero.csv"), index=False, float_format="%.6g")
        msg += f"; share_zero rows {len(z)}"
    print(msg)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
