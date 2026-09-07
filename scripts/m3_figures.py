"""Generate the M3 figures:

- figures/m3_ground_track.png : two-body + spherical-Earth-rotation
  ground track over one sidereal day, with apogee/perigee markers and the
  representative site. No J2.
- figures/m3_access_vs_time.png : elevation vs time over one sidereal day
  with the 10 deg threshold, shaded access intervals, apogee markers.
- figures/m3_range_vs_elevation.png : supporting figure, range vs
  elevation colored by time, showing the geometric relationship during
  each pass.

Run:
    python scripts/m3_figures.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from molniya_design.access import access_metrics, find_access_intervals, range_az_el_series
from molniya_design.constants import (
    BASELINE_ELEMENTS_DEG,
    M1_PERIOD_S,
    MIN_ELEVATION_DEG,
    SIDEREAL_DAY_S,
    SITE_LAT_DEG,
    SITE_LON_DEG,
)
from molniya_design.elements import coe_to_rv
from molniya_design.groundtrack import apsis_ground_points, compute_ground_track

FIG_DIR = Path(__file__).resolve().parent.parent / "figures"
BASE = BASELINE_ELEMENTS_DEG


def _split_on_wrap(lon, lat):
    """Split a (lon, lat) polyline into segments wherever longitude wraps
    across the +/-180 deg boundary, so plotting doesn't draw a fake
    line straight across the map."""
    segments = []
    seg_lon, seg_lat = [lon[0]], [lat[0]]
    for k in range(1, len(lon)):
        if abs(lon[k] - lon[k - 1]) > 180.0:
            segments.append((seg_lon, seg_lat))
            seg_lon, seg_lat = [], []
        seg_lon.append(lon[k])
        seg_lat.append(lat[k])
    segments.append((seg_lon, seg_lat))
    return segments


def figure_ground_track(r0, v0):
    t_span = (0.0, SIDEREAL_DAY_S)
    gt = compute_ground_track(r0, v0, t_span, dt_s=30.0)
    apsides = apsis_ground_points(r0, v0, (0.0, 2.05 * M1_PERIOD_S))

    fig, ax = plt.subplots(figsize=(12, 6.5))

    for seg_lon, seg_lat in _split_on_wrap(gt.lon_deg, gt.lat_deg):
        ax.plot(seg_lon, seg_lat, color="firebrick", lw=1.2, zorder=3)
    ax.plot([], [], color="firebrick", lw=1.2, label="Ground track (two-body + spherical-Earth rotation, 1 sidereal day)")

    perigees = [p for p in apsides if p.kind == "perigee"]
    apogees = [p for p in apsides if p.kind == "apogee"]
    ax.scatter(
        [p.lon_deg for p in apogees], [p.lat_deg for p in apogees],
        color="darkgreen", marker="s", s=70, zorder=5, label="Apogee",
    )
    ax.scatter(
        [p.lon_deg for p in perigees], [p.lat_deg for p in perigees],
        color="darkorange", marker="^", s=70, zorder=5, label="Perigee",
    )
    ax.scatter(
        [SITE_LON_DEG], [SITE_LAT_DEG], color="blue", marker="*", s=220,
        zorder=6, edgecolor="black", linewidth=0.5,
        label=f"Representative site ({SITE_LAT_DEG:.0f}N, {SITE_LON_DEG:.0f}E)",
    )

    ax.axhline(0, color="gray", lw=0.6, alpha=0.6)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xticks(np.arange(-180, 181, 30))
    ax.set_yticks(np.arange(-90, 91, 30))
    ax.grid(alpha=0.3)
    ax.set_xlabel("Longitude [deg]")
    ax.set_ylabel("Geocentric latitude [deg]")
    ax.set_title(
        "M3: Baseline Molniya ground track — two-body propagation + spherical-Earth rotation (no J2)\n"
        "Earth-fixed (ECEF) geocentric lat/lon, engineering-reference epoch (theta_G0=0 at t=0)"
    )
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2, fontsize=9)
    fig.tight_layout()
    out = FIG_DIR / "m3_ground_track.png"
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")
    return gt


def figure_access_vs_time(r0, v0, gt=None):
    t_span = (0.0, SIDEREAL_DAY_S)
    if gt is None:
        gt = compute_ground_track(r0, v0, t_span, dt_s=15.0)
    rng, az, el = range_az_el_series(gt.r_ecef, SITE_LAT_DEG, SITE_LON_DEG)
    intervals = find_access_intervals(gt.t_s, el, threshold_deg=MIN_ELEVATION_DEG)

    apsides = apsis_ground_points(r0, v0, (0.0, 2.05 * M1_PERIOD_S))
    apogee_times = [p.t_s for p in apsides if p.kind == "apogee" and p.t_s <= SIDEREAL_DAY_S + 1.0]

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(gt.t_s / 3600.0, el, color="steelblue", lw=1.2, label="Elevation (geometric, spherical Earth)")
    ax.axhline(MIN_ELEVATION_DEG, color="crimson", ls="--", lw=1.2, label=f"Min elevation = {MIN_ELEVATION_DEG:.0f} deg")

    for iv in intervals:
        ax.axvspan(iv.start_s / 3600.0, iv.end_s / 3600.0, color="steelblue", alpha=0.15)

    for t in apogee_times:
        ax.axvline(t / 3600.0, color="darkgreen", ls=":", lw=1.0)

    # legend proxies (avoid drawing the shaded/line artists twice, which
    # would double-stack alpha and make the labeled interval visibly
    # darker than the others)
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    access_proxy = Patch(facecolor="steelblue", alpha=0.15, label="Access interval (elevation >= threshold)")
    apogee_proxy = Line2D([0], [0], color="darkgreen", ls=":", lw=1.0, label="Apogee time")

    ax.set_xlim(0, SIDEREAL_DAY_S / 3600.0)
    ax.set_xlabel("Time [hours]")
    ax.set_ylabel("Elevation [deg]")
    ax.set_title(
        "M3: Geometric elevation vs time at representative site over one sidereal day\n"
        "(spherical-Earth line-of-sight geometry only — NOT a communications-link or availability result)"
    )
    handles, labels = ax.get_legend_handles_labels()
    handles += [access_proxy, apogee_proxy]
    labels += [access_proxy.get_label(), apogee_proxy.get_label()]
    ax.legend(handles, labels, loc="upper right", fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out = FIG_DIR / "m3_access_vs_time.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"wrote {out}")
    return gt, rng, az, el, intervals


def figure_range_vs_elevation(gt, rng, el):
    fig, ax = plt.subplots(figsize=(8, 6.5))
    sc = ax.scatter(el, rng, c=gt.t_s / 3600.0, cmap="viridis", s=6)
    ax.axvline(MIN_ELEVATION_DEG, color="crimson", ls="--", lw=1.0, label=f"Min elevation = {MIN_ELEVATION_DEG:.0f} deg")
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("Time [hours]")
    ax.set_xlabel("Elevation [deg]")
    ax.set_ylabel("Range [km]")
    ax.set_title("M3: Range vs elevation, colored by time\n(geometric line-of-sight only)")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out = FIG_DIR / "m3_range_vs_elevation.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"wrote {out}")


def main():
    r0, v0 = coe_to_rv(
        BASE["a_km"], BASE["e"], BASE["i_deg"], BASE["raan_deg"],
        BASE["argp_deg"], BASE["nu0_deg"],
    )
    FIG_DIR.mkdir(exist_ok=True)
    figure_ground_track(r0, v0)
    gt, rng, az, el, intervals = figure_access_vs_time(r0, v0)
    figure_range_vs_elevation(gt, rng, el)


if __name__ == "__main__":
    main()
