"""Generate Stata (.dta, v118) twins of the published data products.

  data/dd_master.dta  — synthetic microdata (the master file)
  data/dd_series.dta  — quantiles/levels/shares (PSID & SCF views, normal +
                        detrended) + smoothed factors, one row per quarter

Coefficient files stay CSV-only (pipeline users). Factor DRAWS deliberately
excluded — dataset notes point to them for proper uncertainty treatment.
"""
import pandas as pd
import re

CITE = "Bayer, Calderon & Kuhn: Distributional Dynamics. https://www.distributionaldynamics.com"

def add_tq(df):
    yq = df["time"].str.extract(r"(\d{4})-Q(\d)")
    df.insert(0, "quarter", yq[1].astype(int))
    df.insert(0, "year", yq[0].astype(int))
    return df

def grid_name(sfx):
    return "top" if abs(sfx - 0.999999) < 1e-9 else f"p{int(round(sfx*100))}"

def series_label(meas, stat, sfx, vlab):
    """Variable label. All three stats are DECILE-GROUP values (2019 USD):
    `quantiles` = mean per household of the group, `levels` = group total,
    `shares` = group share of the total (the ten shares sum to 1).
    p10 = group 0-10th pct, ..., p90 = 80-90th pct, top = top 10%."""
    name = {"consum": "Consumption", "income": "Income", "wealth": "Wealth"}[meas]
    grp = ("top 10%" if abs(sfx - 0.999999) < 1e-9
           else f"decile group ending at {int(round(sfx*100))}th pct")
    what = {"quantiles": "mean/HH", "levels": "total", "shares": "share"}[stat]
    return f"{name} {what} of {grp} ({vlab})"


# ── dd_series ────────────────────────────────────────────────────────────────
frames = []
labels = {}
for view, vtag, vlab in [("PSID_functional_data", "psid", "PSID view"),
                         ("PSID_functional_data_detrended", "psid_dt", "PSID view, detrended"),
                         ("SCF_functional_data", "scf", "SCF view"),
                         ("SCF_functional_data_detrended", "scf_dt", "SCF view, detrended")]:
    df = pd.read_csv(f"data/{view}.csv")
    keep = ["time"] + [c for c in df.columns if not c.startswith("ciw_") and c != "time"]
    df = df[keep]
    ren = {}
    for c in df.columns:
        if c == "time":
            continue
        m = re.match(r"(quantiles|levels|shares)(consum|income|wealth)_([0-9.]+)", c)
        stat, meas, sfx = m.group(1), m.group(2), float(m.group(3))
        new = f"{meas}_{stat[:-1] if stat != 'levels' else 'level'}_{grid_name(sfx)}_{vtag}"
        ren[c] = new
        labels[new] = series_label(meas, stat, sfx, vlab)
    df = df.rename(columns=ren)
    frames.append(df.set_index("time"))

fac = pd.read_csv("data/smoothed_factors.csv").set_index("time")
fac.columns = [f"factor_{c[1:]}" for c in fac.columns]
for c in fac.columns:
    k = int(c.split('_')[1])
    labels[c] = (f"Distributional factor {k}" if k <= 8 else
                 f"Lag {(k-9)//8+1} of distributional factor {(k-9)%8+1}" if k <= 32 else
                 f"Aggregate factor {k-32}")
frames.append(fac)

series = pd.concat(frames, axis=1).reset_index()
series = add_tq(series)
labels.update({"year": "Year", "quarter": "Quarter (1-4)", "time": "Quarter label (YYYY-Qq)"})
series.to_stata("data/dd_series.dta", write_index=False, version=118,
                variable_labels={k: v[:80] for k, v in labels.items() if k in series.columns},
                data_label="Distributional Dynamics: reconstructed series"[:80])

# ── dd_master (microdata) ────────────────────────────────────────────────────
micro = pd.read_csv("data/PSID_synthetic_microdata.csv")
micro = add_tq(micro)
mlab = {
    "year": "Year", "quarter": "Quarter (1-4)", "time": "Quarter label (YYYY-Qq)",
    "cop_share": "Copula cell probability mass (weight)",
    "grid_point": "Copula grid cell index",
    "consum": "Consumption, decile-group mean (2019 USD per household)",
    "income": "Income, decile-group mean (2019 USD per household)",
    "wealth": "Wealth, decile-group mean (2019 USD per household)",
    "consumgrid": "Consumption grid coordinate (rank cell)",
    "incomegrid": "Income grid coordinate (rank cell)",
    "wealthgrid": "Wealth grid coordinate (rank cell)",
}
micro.to_stata("data/dd_master.dta", write_index=False, version=118,
               variable_labels={k: v[:80] for k, v in mlab.items() if k in micro.columns},
               data_label="Distributional Dynamics: synthetic microdata"[:80])
print("written:", "dd_series.dta", series.shape, "| dd_master.dta", micro.shape)

# ── CSV twins (same harmonized names; URL-readable from Julia/Python/R) ──────
series.to_csv("data/dd_series.csv", index=False)
micro.to_csv("data/dd_master.csv", index=False)
print("csv twins written")
