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

def grid_label(sfx):
    return "top (99.9999th pct)" if abs(sfx - 0.999999) < 1e-9 else f"{int(round(sfx*100))}th pct"

def grid_name(sfx):
    return "top" if abs(sfx - 0.999999) < 1e-9 else f"p{int(round(sfx*100))}"

STAT_LABEL = {"quantiles": "quantile", "levels": "decile-group mean level", "shares": "cumulative share"}

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
        labels[new] = f"{meas.capitalize()}, {STAT_LABEL[stat]} at {grid_label(sfx)} ({vlab})"
    df = df.rename(columns=ren)
    frames.append(df.set_index("time"))

fac = pd.read_csv("data/smoothed_factors.csv").set_index("time")
fac.columns = [f"factor_{c[1:]}" for c in fac.columns]
for c in fac.columns:
    labels[c] = f"Smoothed latent factor {c.split('_')[1]} (see repository docs for ordering)"
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
    "consum": "Consumption (relative to per-household average)",
    "income": "Income (relative to per-household average)",
    "wealth": "Wealth (relative to per-household average)",
    "consumgrid": "Consumption grid coordinate (rank cell)",
    "incomegrid": "Income grid coordinate (rank cell)",
    "wealthgrid": "Wealth grid coordinate (rank cell)",
}
micro.to_stata("data/dd_master.dta", write_index=False, version=118,
               variable_labels={k: v[:80] for k, v in mlab.items() if k in micro.columns},
               data_label="Distributional Dynamics: synthetic microdata"[:80])
print("written:", "dd_series.dta", series.shape, "| dd_master.dta", micro.shape)
