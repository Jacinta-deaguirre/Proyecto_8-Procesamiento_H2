"""Figure 5 of the manuscript: radial temperature and adsorption profiles during charging
(Test No. 13). Reads outputs/charge_temperature_profiles.csv and
outputs/charge_adsorption_profiles.csv written by finite_model_xiao.jl.
Run:  python plot_radial_profiles.py
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def create_plot(folder="outputs/"):
    df_temp = pd.read_csv(f"{folder}charge_temperature_profiles.csv")
    df_ads = pd.read_csv(f"{folder}charge_adsorption_profiles.csv")
    x_radius = df_temp["Radius_m"]
    t_cols = [c for c in df_temp.columns if c.startswith("t_")]
    num_profiles_to_plot, num_stationary_points = 16, 3
    num_linspace_points = num_profiles_to_plot - num_stationary_points
    if len(t_cols) < num_profiles_to_plot:
        indices = np.arange(len(t_cols))
    else:
        stop_index = len(t_cols) - 1 - num_stationary_points
        linspace_indices = np.linspace(0, stop_index, num=num_linspace_points, dtype=int)
        last_three = np.arange(len(t_cols) - num_stationary_points, len(t_cols))
        indices = np.concatenate((linspace_indices, last_three))
    cols_to_plot = [t_cols[i] for i in indices]
    colors = plt.cm.viridis_r(np.linspace(0, 1, num_profiles_to_plot))

    plt.style.use("seaborn-v0_8-talk")
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Inter", "DejaVu Sans", "Arial"]
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(18, 9), sharex=True)
    handles = []
    for i, col in enumerate(cols_to_plot):
        ax1.plot(x_radius, df_temp[col], color=colors[i], linestyle="-")
        ax2.plot(x_radius, df_ads[col], color=colors[i], linestyle="-")
        label = col.replace("t_", "").replace("s", " s")
        handles.append(plt.Line2D([0], [0], color=colors[i], lw=2.5, label=label))
    ax1.set_title("Temperature Profiles", fontsize=22); ax1.set_ylabel("Temperature (K)", fontsize=20)
    ax1.set_xlabel("Radius (m)", fontsize=20); ax1.grid(True, linestyle=":", alpha=0.6); ax1.tick_params(axis="both", which="major", labelsize=16)
    ax2.set_title("Adsorption Profiles", fontsize=22); ax2.set_ylabel("Adsorption (mol/kg)", fontsize=20)
    ax2.set_xlabel("Radius (m)", fontsize=20); ax2.grid(True, linestyle=":", alpha=0.6); ax2.tick_params(axis="both", which="major", labelsize=16)
    ax1.set_xlim(left=0)
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.01), title="Time (s)", fontsize=14,
               ncols=int(np.ceil(len(cols_to_plot) / 2)))
    fig.tight_layout(); plt.subplots_adjust(top=0.88, bottom=0.2, wspace=0.25)
    for ext in ("svg", "png"):
        plt.savefig(folder + f"charge_radial_profiles_subplots.{ext}", dpi=200, bbox_inches="tight", format=ext)
    print("saved", folder + "charge_radial_profiles_subplots.svg/.png")


if __name__ == "__main__":
    create_plot(sys.argv[1] if len(sys.argv) > 1 else "outputs/")
