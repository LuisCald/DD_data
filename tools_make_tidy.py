"""Build data/dd_tidy_long.csv: group means and shares in one long table.

One row per (quarter, survey view, variable, ranking variable, grouping, group):

  time, year, quarter   quarter (YYYY-Qq) and its parts
  survey                psid | scf | cps | cps1967 (the model's PSID, SCF or CPS view; the CPS
                        views are income only: cps = 1993 on, cps1967 = 1967-1993)
  variable              consum | income | wealth (the outcome)
  ranked_by             the variable households are ranked by. Equal to `variable` for
                        the usual own-ranked groups (e.g. the wealth share of the top 10% of
                        wealth); different for cross groups (e.g. consumption of the
                        bottom income quintile)
  grouping              deciles | quintiles | wealth_groups (bottom 50 / next 40 / top 10)
  group                 d1..d10 | 0-20..80-100 | bottom50, next40, top10
  mean                  group mean per household, 2019Q4 USD (CPI-U)
  share                 group's share of the variable's total (cross groups: of the
                        groups' combined total)

Own-ranked groups come from dd_series (decile groups, aggregated for quintiles and
50/40/10). Cross groups come from the 2D maximum-entropy step
(python/cross_conditional_means.py), which is the recommended method for group means of
one variable by another. Rows with no estimate (e.g. PSID wealth before 1984,
consumption before 1999, SCF consumption) are left out.

Usage: python3 tools_make_tidy.py      (runs cross_conditional_means.py; a few minutes)
"""
import os
import subprocess
import sys
import tempfile

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DECILES = ["p10", "p20", "p30", "p40", "p50", "p60", "p70", "p80", "p90", "top"]
GROUPINGS = {
    "deciles": [(f"d{i + 1}", [d]) for i, d in enumerate(DECILES)],
    "quintiles": [(f"{20 * i}-{20 * (i + 1)}", DECILES[2 * i:2 * i + 2]) for i in range(5)],
    "wealth_groups": [("bottom50", DECILES[:5]), ("next40", DECILES[5:9]), ("top10", DECILES[9:])],
}
CROSS_LABELS = {"bot50": "bottom50", "mid40": "next40", "top10": "top10",
                **{f"p{20 * i}-{20 * (i + 1)}": f"{20 * i}-{20 * (i + 1)}" for i in range(5)}}
CROSS_RUNS = {"wealth_groups": ("50,40,10", [0.5, 0.4, 0.1]), "quintiles": ("20,20,20,20,20", [0.2] * 5)}


def own_ranked(series):
    rows = []
    for survey in ["psid", "scf", "cps", "cps1967"]:
        for var in ["consum", "income", "wealth"]:
            mcol = {d: f"{var}_quantile_{d}_{survey}" for d in DECILES}
            scol = {d: f"{var}_share_{d}_{survey}" for d in DECILES}
            if mcol["p10"] not in series:
                continue
            for grouping, groups in GROUPINGS.items():
                for g, ds in groups:
                    rows.append(pd.DataFrame({
                        "time": series.time, "survey": survey, "variable": var, "ranked_by": var,
                        "grouping": grouping, "group": g,
                        # deciles are equal-sized, so the group mean is the mean of decile means;
                        # min_count keeps a group missing when any decile is missing
                        "mean": series[[mcol[d] for d in ds]].mean(axis=1, skipna=False),
                        "share": series[[scol[d] for d in ds]].sum(axis=1, min_count=len(ds))}))
    return pd.concat(rows, ignore_index=True)


def cross_groups(tmp):
    rows = []
    for dataset in ["PSID", "SCF"]:
        for grouping, (sizes, widths) in CROSS_RUNS.items():
            out = os.path.join(tmp, f"{dataset}_{grouping}.csv")
            subprocess.run([sys.executable, os.path.join(HERE, "python", "cross_conditional_means.py"),
                            "--dataset", dataset, "--date", "all", "--groups", sizes, "--out", out],
                           check=True, stdout=subprocess.DEVNULL)
            c = pd.read_csv(out)
            c = c[c.outcome != c.conditioning]
            w = dict(zip(sorted(c.group.unique(), key=list(CROSS_LABELS).index), widths))
            rows.append(pd.DataFrame({
                "time": c.date, "survey": dataset.lower(), "variable": c.outcome, "ranked_by": c.conditioning,
                "grouping": grouping, "group": c.group.map(CROSS_LABELS),
                "mean": c.level, "share": c.rel_mean * c.group.map(w)}))
    df = pd.concat(rows, ignore_index=True)
    # shares of the groups' combined total (the max-ent groups add up to the mean the
    # basis can represent, which misses part of the extreme wealth tail)
    key = ["time", "survey", "variable", "ranked_by", "grouping"]
    df["share"] /= df.groupby(key).share.transform("sum")
    return df


def main():
    os.chdir(HERE)
    series = pd.read_csv("data/dd_series.csv")
    with tempfile.TemporaryDirectory() as tmp:
        df = pd.concat([own_ranked(series), cross_groups(tmp)], ignore_index=True)
    df = df.dropna(subset=["mean", "share"])
    yq = df.time.str.extract(r"(\d{4})-Q(\d)").astype(int)
    df.insert(1, "year", yq[0])
    df.insert(2, "quarter", yq[1])
    order = {k: i for i, k in enumerate(["deciles", "quintiles", "wealth_groups"])}
    df = df.sort_values(["survey", "variable", "ranked_by", "grouping", "time"],
                        key=lambda s: s.map(order) if s.name == "grouping" else s, kind="stable")
    df.to_csv("data/dd_tidy_long.csv", index=False, float_format="%.6g")
    print(f"wrote data/dd_tidy_long.csv: {len(df)} rows, {df.time.min()}..{df.time.max()}")


if __name__ == "__main__":
    main()
