"""Generate the M5 figures:

- figures/m5_regional_max_gap.png : regional max-gap heatmap
- figures/m5_access_fraction_map.png : regional access-fraction heatmap
- figures/m5_sensitivity_trade.png : 4-panel sensitivity trade study
- figures/m5_representative_site_timeline.png : optional 14-day site
  elevation/access timeline

Reads results/m5_grid_metrics.csv and results/m5_sensitivity.csv (written
by scripts/m5_verification_report.py — run that first) plus recomputes
the representative-site timeline directly.

Run:
    python scripts/m5_verification_report.py   # first, to produce CSVs
    python scripts/m5_figures.py
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from molniya_design.access import find_access_intervals, range_az_el_series
from molniya_design.constants import (
    BASELINE_ELEMENTS_DEG,
    MIN_ELEVATION_DEG,
    SITE_LAT_DEG,
    SITE_LON_DEG,
)
from molniya_design.coverage import build_satellite_ecef_trajectory
from molniya_design.j2 import SecularElements
from molniya_design.propagation import true_to_mean_anomaly

FIG_DIR = Path(__file__).resolve().parent.parent / "figures"
RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
BASE = BASELINE_ELEMENTS_DEG
M0_DEG = float(np.degrees(true_to_mean_anomaly(np.radians(BASE["nu0_deg"]), BASE["e"])))


def _load_grid_csv():
    rows = []
    with open(RESULTS_DIR / "m5_grid_metrics.csv") as f:
        for row in csv.DictReader(f):
            rows.append({k: float(v) for k, v in row.items()})
    return rows


def _load_sensitivity_csv():
    rows = []
    with open(RESULTS_DIR / "m5_sensitivity.csv") as f:
        for row in csv.DictReader(f):
            rows.append(row)
    return rows


def figure_regional_maps():
    rows = _load_grid_csv()
    lats = sorted(set(r["lat_deg"] for r in rows))
    lons = sorted(set(r["lon_deg"] for r in rows))
    lat_idx = {v: i for i, v in enumerate(lats)}
    lon_idx = {v: i for i, v in enumerate(lons)}

    gap_grid = np.full((len(lats), len(lons)), np.nan)
    access_grid = np.full((len(lats), len(lons)), np.nan)
    for r in rows:
        i, j = lat_idx[r["lat_deg"]], lon_idx[r["lon_deg"]]
        gap_grid[i, j] = r["max_gap_s"] / 3600.0
        access_grid[i, j] = r["access_fraction"] * 100.0

    worst_idx = np.unravel_index(np.nanargmax(gap_grid), gap_grid.shape)
    worst_lat, worst_lon = lats[worst_idx[0]], lons[worst_idx[1]]
    best_idx = np.unravel_index(np.nanargmax(access_grid), access_grid.shape)
    best_lat, best_lon = lats[best_idx[0]], lons[best_idx[1]]

    lon_edges = np.array(lons + [360.0]) - (lons[1] - lons[0]) / 2.0
    lat_step = lats[1] - lats[0]
    lat_edges = np.array(lats + [lats[-1] + lat_step]) - lat_step / 2.0

    # --- Figure 1: max gap ---
    fig, ax = plt.subplots(figsize=(12, 5.5))
    pcm = ax.pcolormesh(lon_edges, lat_edges, gap_grid, cmap="inferno_r", shading="flat")
    cbar = fig.colorbar(pcm, ax=ax)
    cbar.set_label("Maximum no-access gap [hours]")
    ax.scatter([SITE_LON_DEG], [SITE_LAT_DEG], color="cyan", marker="*", s=220,
               edgecolor="black", linewidth=0.6, zorder=5, clip_on=False,
               label=f"Representative site ({SITE_LAT_DEG:.0f}N,{SITE_LON_DEG:.0f}E)")
    ax.scatter([worst_lon], [worst_lat], color="red", marker="X", s=140,
               edgecolor="white", linewidth=0.8, zorder=5, clip_on=False,
               label=f"Worst point ({worst_lat:.1f}N,{worst_lon:.1f}E)")
    ax.set_xlabel("Longitude [deg]")
    ax.set_ylabel("Latitude [deg N]")
    ax.set_xlim(0, 360)
    ax.set_ylim(lat_edges[0], lat_edges[-1])
    ax.set_title(
        "M5: Regional maximum no-access gap, 60-75N, 14-day horizon\n"
        "Geometric access, first-order secular J2, spherical Earth — NOT RF coverage"
    )
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, fontsize=9)
    fig.tight_layout()
    out = FIG_DIR / "m5_regional_max_gap.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"wrote {out}")

    # --- Figure 2: access fraction ---
    fig, ax = plt.subplots(figsize=(12, 5.5))
    pcm = ax.pcolormesh(lon_edges, lat_edges, access_grid, cmap="viridis", shading="flat")
    cbar = fig.colorbar(pcm, ax=ax)
    cbar.set_label("Access fraction [%]")
    ax.scatter([SITE_LON_DEG], [SITE_LAT_DEG], color="red", marker="*", s=220,
               edgecolor="black", linewidth=0.6, zorder=5, label=f"Representative site ({SITE_LAT_DEG:.0f}N,{SITE_LON_DEG:.0f}E)")
    ax.scatter([best_lon], [best_lat], color="lime", marker="P", s=140,
               edgecolor="black", linewidth=0.8, zorder=5, label=f"Best point ({best_lat:.1f}N,{best_lon:.1f}E)")
    ax.set_xlabel("Longitude [deg]")
    ax.set_ylabel("Latitude [deg N]")
    ax.set_xlim(0, 360)
    ax.set_ylim(lat_edges[0], lat_edges[-1])
    ax.set_title(
        "M5: Regional geometric access fraction, 60-75N, 14-day horizon\n"
        "Elevation >= 10 deg line-of-sight access (NOT the M1 anomaly dwell proxy, NOT RF coverage)"
    )
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, fontsize=9)
    fig.tight_layout()
    out = FIG_DIR / "m5_access_fraction_map.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"wrote {out}")


def figure_sensitivity_trade():
    rows = _load_sensitivity_csv()
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    # panel 1: worst gap vs inclination
    incl = [r for r in rows if r["study"] == "inclination"]
    x = [float(r["value"]) for r in incl]
    y = [float(r["worst_max_gap_s"]) / 3600.0 for r in incl]
    ax = axes[0, 0]
    ax.plot(x, y, "o-", color="steelblue")
    ax.axvline(BASE["i_deg"], color="crimson", ls="--", lw=1.2, label=f"Baseline (critical) i={BASE['i_deg']:.4f} deg")
    ax.set_xlabel("Inclination [deg]")
    ax.set_ylabel("Worst regional max gap [h]")
    ax.set_title("Worst max gap vs inclination")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # panel 2: worst gap vs min elevation
    elev = [r for r in rows if r["study"] == "elevation_min"]
    x = [float(r["value"]) for r in elev]
    y = [float(r["worst_max_gap_s"]) / 3600.0 for r in elev]
    ax = axes[0, 1]
    ax.plot(x, y, "o-", color="darkorange")
    ax.axvline(MIN_ELEVATION_DEG, color="crimson", ls="--", lw=1.2, label=f"Baseline el_min={MIN_ELEVATION_DEG:.0f} deg")
    ax.set_xlabel("Minimum elevation [deg]")
    ax.set_ylabel("Worst regional max gap [h]")
    ax.set_title("Worst max gap vs minimum elevation")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # panel 3: access fraction vs perigee altitude
    peri = [r for r in rows if r["study"] == "perigee_altitude"]
    x = [float(r["value"]) for r in peri]
    y = [float(r["regional_mean_access_fraction"]) * 100.0 for r in peri]
    ax = axes[1, 0]
    ax.plot(x, y, "o-", color="darkgreen")
    ax.axvline(600.0, color="crimson", ls="--", lw=1.2, label="Baseline perigee alt = 600 km")
    ax.set_xlabel("Perigee altitude [km]")
    ax.set_ylabel("Regional mean access fraction [%]")
    ax.set_title("Access fraction vs perigee altitude (a fixed)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # panel 4: access metric vs argument of perigee (exclude omega=90 outlier from main view, annotate separately)
    argp = [r for r in rows if r["study"] == "argp"]
    x = sorted(set(float(r["value"]) for r in argp))
    x_north = [v for v in x if 200 <= v <= 340]
    y_north = [
        float(next(r["regional_mean_access_fraction"] for r in argp if float(r["value"]) == v)) * 100.0
        for v in x_north
    ]
    ax = axes[1, 1]
    ax.plot(x_north, y_north, "o-", color="purple")
    ax.axvline(270.0, color="crimson", ls="--", lw=1.2, label="Baseline argp = 270 deg (northern)")
    omega90_val = next(float(r["regional_mean_access_fraction"]) * 100.0 for r in argp if float(r["value"]) == 90.0)
    ax.set_xlabel("Argument of perigee [deg]")
    ax.set_ylabel("Regional mean access fraction [%]")
    ax.set_title(f"Access fraction vs arg. of perigee\n(omega=90 deg control case: {omega90_val:.2f}% -- off chart, catastrophic)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.suptitle("M5: Design sensitivity trade study (regional 60-75N band, 7-day horizon, J2 secular)", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = FIG_DIR / "m5_sensitivity_trade.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"wrote {out}")


def figure_site_timeline():
    e0 = SecularElements(a_km=BASE["a_km"], e=BASE["e"], i_deg=BASE["i_deg"],
                          raan_deg=BASE["raan_deg"], argp_deg=BASE["argp_deg"], m_deg=M0_DEG)
    t_span = (0.0, 14 * 86400.0)
    dt = 120.0
    n = int(round((t_span[1] - t_span[0]) / dt)) + 1
    t_s = np.linspace(t_span[0], t_span[1], n)
    traj = build_satellite_ecef_trajectory(e0, t_s)
    r_ecef_3n = traj.T
    _, _, el = range_az_el_series(r_ecef_3n, SITE_LAT_DEG, SITE_LON_DEG)
    intervals = find_access_intervals(t_s, el, threshold_deg=MIN_ELEVATION_DEG)

    fig, ax = plt.subplots(figsize=(14, 5.5))
    ax.plot(t_s / 86400.0, el, color="steelblue", lw=0.8)
    ax.axhline(MIN_ELEVATION_DEG, color="crimson", ls="--", lw=1.0, label=f"Min elevation = {MIN_ELEVATION_DEG:.0f} deg")
    for iv in intervals:
        ax.axvspan(iv.start_s / 86400.0, iv.end_s / 86400.0, color="steelblue", alpha=0.15)
    from matplotlib.patches import Patch
    access_proxy = Patch(facecolor="steelblue", alpha=0.15, label="Access interval")
    handles, labels = ax.get_legend_handles_labels()
    handles.append(access_proxy)
    labels.append(access_proxy.get_label())
    ax.legend(handles, labels, loc="upper right", fontsize=9)
    ax.set_xlim(0, 14)
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Elevation [deg]")
    ax.set_title(
        "M5: Representative-site (65N,40E) elevation/access timeline, 14 days\n"
        "First-order secular J2, geometric access only"
    )
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out = FIG_DIR / "m5_representative_site_timeline.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"wrote {out}")


def main():
    FIG_DIR.mkdir(exist_ok=True)
    figure_regional_maps()
    figure_sensitivity_trade()
    figure_site_timeline()


if __name__ == "__main__":
    main()
