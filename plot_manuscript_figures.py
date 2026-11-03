"""Manuscript Figs. 2-4 (Test 13) and 6-8 (Test 20) from the final model runs.

Inputs : outputs/extra/t13_base/summary.csv and outputs/extra/t20_base/summary.csv (written by extra_runs2.jl,
         1 s dense output of the full cycle) and the digitised experimental data in inputs/.
Outputs: outputs/manuscript_figs/fig{2,3,4,6,7,8}_*.{pdf,png}
Run    : python plot_manuscript_figures.py
Temperature points carry the +-2.2 K type-K thermocouple uncertainty; the pressure uncertainty (+-2.5 kPa) is
smaller than the symbols; the reconstructed phase masses carry no independent uncertainty.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.abspath(__file__))
INP = os.path.join(ROOT, "inputs")
OUT = os.path.join(ROOT, "outputs", "manuscript_figs")
os.makedirs(OUT, exist_ok=True)
V, eps, rho_s, M = 2.4946e-3, 0.49, 517.6, 2.0159e-3
Ms = rho_s * (1 - eps) * V
T_ERR = 2.2  # K

plt.rcParams.update({"font.family": "sans-serif", "font.size": 11, "axes.labelsize": 12, "legend.fontsize": 9,
                     "axes.grid": True, "grid.linestyle": ":", "grid.alpha": 0.6, "figure.dpi": 150})
C_CENTRE, C_MID, C_WALL = "#D62728", "#9467BD", "#7F7F7F"
C_GAS, C_ADS, C_TOT, C_P = "#1F77B4", "#FF7F0E", "#2CA02C", "#1F77B4"

def read_exp(name):
    return pd.read_csv(os.path.join(INP, name), skipinitialspace=True)

def load(tag):
    d = pd.read_csv(os.path.join(ROOT, "outputs", "extra", tag, "summary.csv"))
    return d.drop_duplicates(subset="t", keep="last")

def stage_lines(ax, bounds):
    for b in bounds[1:-1]:
        ax.axvline(b, color="k", lw=0.6, ls="--", alpha=0.5)

def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("saved", name)

# ---------------------------------------------------------------- Test 13
d = load("t13_base"); b13 = [0, 1042, 3189, 4131]
# Fig 2: pressure
fig, ax = plt.subplots(figsize=(6.5, 4.2))
ax.plot(d.t, d.P / 1e6, color=C_P, lw=2, label="1-D model")
for suf in ["", "_dormant", "_discharge"]:
    e = read_exp(f"xiao_pressure_finite_element{suf}.csv")
    ax.plot(e.Time, e.Pressure, "o", ms=5, mfc="white", mec=C_P, mew=1.2, label="Experiment (Richard et al.)" if suf == "" else None)
stage_lines(ax, b13); ax.set_xlabel("Time / s"); ax.set_ylabel("Pressure / MPa"); ax.set_xlim(0, 4200); ax.set_ylim(0, 10)
ax.legend(loc="lower center"); save(fig, "fig2_pressure_no13")

# Fig 3: masses
fig, ax = plt.subplots(figsize=(6.5, 4.2))
m_gas = d.rho_avg * V * eps * 1e3; m_ads = Ms * M * d.na_avg_cyl * 1e3
for y, c, lab, key in [(m_gas + m_ads, C_TOT, "Total", "total"), (m_ads, C_ADS, "Adsorbed", "adsorbed"), (m_gas, C_GAS, "Gas", "gas")]:
    ax.plot(d.t, y, color=c, lw=2, label=f"{lab}, model")
    for suf in ["", "_dormant", "_discharge"]:
        e = read_exp(f"xiao_{key}_mass_finite_element{suf}.csv")
        ax.plot(e.Time, e.Mass * 1e3, "o", ms=5, mfc="white", mec=c, mew=1.2, label=f"{lab}, reconstructed (Xiao et al.)" if suf == "" else None)
stage_lines(ax, b13); ax.set_xlabel("Time / s"); ax.set_ylabel("Hydrogen mass / g"); ax.set_xlim(0, 4200); ax.set_ylim(0, 24)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3, frameon=False); save(fig, "fig3_mass_no13")

# Fig 4: temperatures with error bars
fig, ax = plt.subplots(figsize=(6.5, 4.6))
for col, c, lab, fn in [("T_center", C_CENTRE, "Centre ($r=0$)", "xiao_center_temp_finite_element"),
                        ("T_mid", C_MID, "Half radius ($r=r_i/2$)", "xiao_middle_temp_finite_element"),
                        ("T_wall", C_WALL, "Wall", "xiao_wall_temp_finite_element")]:
    ax.plot(d.t, d[col], color=c, lw=2, label=f"{lab}, model")
    for suf in ["", "_dormant", "_discharge"]:
        e = read_exp(f"{fn}{suf}.csv")
        ax.errorbar(e.Time, e.Temperature, yerr=T_ERR, fmt="o", ms=5, mfc="white", mec=c, ecolor=c, mew=1.2, capsize=2.5, lw=1,
                    label=f"{lab}, experiment ($\\pm$2.2 K)" if suf == "" else None)
stage_lines(ax, b13); ax.set_xlabel("Time / s"); ax.set_ylabel("Temperature / K"); ax.set_xlim(0, 4200); ax.set_ylim(250, 312)
ax.legend(loc="lower left", ncol=1); save(fig, "fig4_temperature_no13")

# ---------------------------------------------------------------- Test 20
d = load("t20_base"); b20 = [0, 953, 3822, 4694]
fig, ax = plt.subplots(figsize=(6.5, 4.2))
ax.plot(d.t, d.P / 1e6, color=C_P, lw=2, label="1-D model")
e = read_exp("xiao19_pressure.csv")
ax.plot(e.Time, e.Pressure, "o", ms=5, mfc="white", mec=C_P, mew=1.2, label="Experiment (Richard et al.)")
stage_lines(ax, b20); ax.set_xlabel("Time / s"); ax.set_ylabel("Pressure / MPa"); ax.set_xlim(0, 4800); ax.set_ylim(0, 10)
ax.legend(loc="lower center"); save(fig, "fig6_pressure_no20")

fig, ax = plt.subplots(figsize=(6.5, 4.2))
m_gas = d.rho_avg * V * eps * 1e3; m_ads = Ms * M * d.na_avg_cyl * 1e3
for y, c, lab, fn in [(m_gas + m_ads, C_TOT, "Total", "xiao19_total_mass.csv"), (m_ads, C_ADS, "Adsorbed", "xiao19_ads_mass.csv"), (m_gas, C_GAS, "Gas", "xiao19_gas_mass.csv")]:
    ax.plot(d.t, y, color=c, lw=2, label=f"{lab}, model")
    e = read_exp(fn)
    ax.plot(e.Time, e.Mass * 1e3, "o", ms=5, mfc="white", mec=c, mew=1.2, label=f"{lab}, reconstructed (Xiao et al.)")
stage_lines(ax, b20); ax.set_xlabel("Time / s"); ax.set_ylabel("Hydrogen mass / g"); ax.set_xlim(0, 4800); ax.set_ylim(0, 22)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3, frameon=False); save(fig, "fig7_mass_no20")

fig, ax = plt.subplots(figsize=(6.5, 4.6))
ax.plot(d.t, d.T_avg_cyl, color=C_ADS, lw=2, label="Volume-averaged bed temperature, model")
for fn, lab in [("xiao19_avg_temperature.csv", "Average bed temperature, experiment ($\\pm$2.2 K)"), ("xiao19_avg_temperature_dormant.csv", None), ("xiao19_avg_temperature_discharge.csv", None)]:
    e = read_exp(fn)
    ax.errorbar(e.Time, e.Temperature, yerr=T_ERR, fmt="o", ms=5, mfc="white", mec=C_ADS, ecolor=C_ADS, mew=1.2, capsize=2.5, lw=1, label=lab)
stage_lines(ax, b20); ax.set_xlabel("Time / s"); ax.set_ylabel("Temperature / K"); ax.set_xlim(0, 4800); ax.set_ylim(282, 318)
ax.legend(loc="lower left"); save(fig, "fig8_temperature_no20")
