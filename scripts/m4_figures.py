"""Generate the M4 figures:

- figures/m4_j2_rate_sensitivity.png : Omega_dot and omega_dot vs
  inclination, marking the critical inclination and baseline.
- figures/m4_ground_track_two_body_vs_j2.png : 7-day two-body vs
  first-order secular J2 ground-track comparison.
- figures/m4_element_drift.png : Omega(t) and omega(t) over 14 days,
  baseline (critical) vs an off-critical case.

Run:
    python scripts/m4_figures.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from molniya_design.constants import BASELINE_ELEMENTS_DEG, SIDEREAL_DAY_S
from molniya_design.elements import coe_to_rv
from molniya_design.groundtrack import (
    apogee_ground_points_j2,
    apsis_ground_points,
    compute_ground_track,
    compute_ground_track_j2,
)
from molniya_design.j2 import SecularElements, compute_secular_rates, propagate_elements_j2_secular
from molniya_design.propagation import true_to_mean_anomaly

FIG_DIR = Path(__file__).resolve().parent.parent / "figures"
BASE = BASELINE_ELEMENTS_DEG
M0_DEG = float(np.degrees(true_to_mean_anomaly(np.radians(BASE["nu0_deg"]), BASE["e"])))


def figure_rate_sensitivity():
    i_deg = np.linspace(50.0, 75.0, 500)
    raan_dot = np.array([compute_secular_rates(BASE["a_km"], BASE["e"], i).raan_dot_deg_day for i in i_deg])
    argp_dot = np.array([compute_secular_rates(BASE["a_km"], BASE["e"], i).argp_dot_deg_day for i in i_deg])

    i_crit = BASE["i_deg"]

    fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True)

    ax0 = axes[0]
    ax0.plot(i_deg, raan_dot, color="steelblue", lw=1.6)
    ax0.axvline(i_crit, color="crimson", ls="--", lw=1.2, label=f"Critical/baseline i = {i_crit:.4f} deg")
    ax0.scatter([i_crit], [compute_secular_rates(BASE["a_km"], BASE["e"], i_crit).raan_dot_deg_day],
                color="crimson", zorder=5, s=50)
    ax0.set_ylabel("RAAN_dot [deg/day]")
    ax0.set_title("M4: J2 secular rate sensitivity vs inclination (a, e fixed at baseline)")
    ax0.legend(loc="best", fontsize=9)
    ax0.grid(alpha=0.3)

    ax1 = axes[1]
    ax1.plot(i_deg, argp_dot, color="darkgreen", lw=1.6)
    ax1.axhline(0.0, color="gray", lw=0.8)
    ax1.axvline(i_crit, color="crimson", ls="--", lw=1.2, label=f"Critical/baseline i = {i_crit:.4f} deg (argp_dot = 0)")
    ax1.scatter([i_crit], [0.0], color="crimson", zorder=5, s=50)
    ax1.set_xlabel("Inclination [deg]")
    ax1.set_ylabel("argp_dot (omega_dot) [deg/day]")
    ax1.legend(loc="best", fontsize=9)
    ax1.grid(alpha=0.3)

    fig.tight_layout()
    out = FIG_DIR / "m4_j2_rate_sensitivity.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"wrote {out}")


def _split_on_wrap(lon, lat):
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


def figure_ground_track_comparison(elements0):
    t_span = (0.0, 7 * 86400.0)
    r0, v0 = coe_to_rv(BASE["a_km"], BASE["e"], BASE["i_deg"], BASE["raan_deg"], BASE["argp_deg"], BASE["nu0_deg"])

    gt_2body = compute_ground_track(r0, v0, t_span, dt_s=120.0)
    gt_j2 = compute_ground_track_j2(elements0, t_span, dt_s=120.0)

    apogees_2body = [p for p in apsis_ground_points(r0, v0, (0.0, t_span[1] * 1.01)) if p.kind == "apogee"]
    apogees_j2 = apogee_ground_points_j2(elements0, (0.0, t_span[1] * 1.01))

    fig, ax = plt.subplots(figsize=(12, 6.5))

    for seg_lon, seg_lat in _split_on_wrap(gt_2body.lon_deg, gt_2body.lat_deg):
        ax.plot(seg_lon, seg_lat, color="steelblue", lw=1.0, alpha=0.8, zorder=2)
    ax.plot([], [], color="steelblue", lw=1.2, label="Two-body (no J2), 7 days")

    for seg_lon, seg_lat in _split_on_wrap(gt_j2.lon_deg, gt_j2.lat_deg):
        ax.plot(seg_lon, seg_lat, color="darkorange", lw=1.0, ls="--", alpha=0.9, zorder=3)
    ax.plot([], [], color="darkorange", lw=1.2, ls="--", label="First-order secular J2 (mean elements), 7 days")

    ax.scatter([p.lon_deg for p in apogees_2body], [p.lat_deg for p in apogees_2body],
               color="steelblue", marker="s", s=45, zorder=5, edgecolor="black", linewidth=0.4,
               label="Apogee (two-body)")
    ax.scatter([p.lon_deg for p in apogees_j2], [p.lat_deg for p in apogees_j2],
               color="darkorange", marker="D", s=45, zorder=5, edgecolor="black", linewidth=0.4,
               label="Apogee (J2 secular)")

    ax.axhline(0, color="gray", lw=0.6, alpha=0.6)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xticks(np.arange(-180, 181, 30))
    ax.set_yticks(np.arange(-90, 91, 30))
    ax.grid(alpha=0.3)
    ax.set_xlabel("Longitude [deg]")
    ax.set_ylabel("Geocentric latitude [deg]")
    ax.set_title(
        "M4: Two-body vs first-order secular-J2 ground track, 7 days\n"
        "Same spherical-Earth rotation model as M3; mean-element secular J2, not full ephemeris"
    )
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2, fontsize=9)

    # zoomed inset on the last (day-7) apogee, where the two tracks have
    # visibly separated -- at full-map scale the two curves are
    # indistinguishable (J2 drift over 7 days is genuinely small, ~1.17
    # deg -- see DESIGN.md M4), so the inset is needed to actually show it
    last_2body = apogees_2body[-1]
    last_j2 = apogees_j2[-1]
    axins = ax.inset_axes([0.66, 0.46, 0.32, 0.32])
    for seg_lon, seg_lat in _split_on_wrap(gt_2body.lon_deg, gt_2body.lat_deg):
        axins.plot(seg_lon, seg_lat, color="steelblue", lw=1.4, alpha=0.9)
    for seg_lon, seg_lat in _split_on_wrap(gt_j2.lon_deg, gt_j2.lat_deg):
        axins.plot(seg_lon, seg_lat, color="darkorange", lw=1.4, ls="--", alpha=0.9)
    axins.scatter([last_2body.lon_deg], [last_2body.lat_deg], color="steelblue", marker="s", s=55,
                  zorder=5, edgecolor="black", linewidth=0.5)
    axins.scatter([last_j2.lon_deg], [last_j2.lat_deg], color="darkorange", marker="D", s=55,
                  zorder=5, edgecolor="black", linewidth=0.5)
    lon_c = 0.5 * (last_2body.lon_deg + last_j2.lon_deg)
    axins.set_xlim(lon_c - 1.5, lon_c + 1.5)
    axins.set_ylim(last_2body.lat_deg - 0.3, last_2body.lat_deg + 0.3)
    axins.set_title(f"Day-7 apogee zoom (diff = {abs(last_j2.lon_deg-last_2body.lon_deg):.3f} deg lon)", fontsize=7.5, pad=3)
    axins.tick_params(labelsize=7)
    axins.grid(alpha=0.3)
    ax.indicate_inset_zoom(axins, edgecolor="black", alpha=0.5)

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = FIG_DIR / "m4_ground_track_two_body_vs_j2.png"
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def figure_element_drift(elements0):
    t_days = np.linspace(0.0, 14.0, 400)
    t_s = t_days * 86400.0

    e_baseline = propagate_elements_j2_secular(elements0, t_s, wrap=False)

    elements0_65 = SecularElements(
        a_km=elements0.a_km, e=elements0.e, i_deg=65.0,
        raan_deg=elements0.raan_deg, argp_deg=elements0.argp_deg, m_deg=elements0.m_deg,
    )
    e_65 = propagate_elements_j2_secular(elements0_65, t_s, wrap=False)

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    ax0 = axes[0]
    ax0.plot(t_days, e_baseline.raan_deg, color="steelblue", lw=1.6, label=f"RAAN(t), i={elements0.i_deg:.4f} deg (baseline)")
    ax0.set_ylabel("RAAN [deg] (unwrapped)")
    ax0.set_title("M4: Element drift over 14 days under first-order secular J2")
    ax0.legend(loc="best", fontsize=9)
    ax0.grid(alpha=0.3)

    ax1 = axes[1]
    ax1.plot(t_days, e_baseline.argp_deg, color="steelblue", lw=1.8,
              label=f"argp(t), i={elements0.i_deg:.4f} deg (critical/baseline) — effectively frozen")
    ax1.plot(t_days, e_65.argp_deg, color="darkorange", lw=1.6, ls="--",
              label="argp(t), i=65.0 deg (off-critical) — measurable drift")
    ax1.axhline(270.0, color="gray", lw=0.8, ls=":")
    ax1.set_xlabel("Time [days]")
    ax1.set_ylabel("Argument of perigee [deg] (unwrapped)")
    ax1.legend(loc="best", fontsize=9)
    ax1.grid(alpha=0.3)

    fig.tight_layout()
    out = FIG_DIR / "m4_element_drift.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"wrote {out}")


def main():
    elements0 = SecularElements(
        a_km=BASE["a_km"], e=BASE["e"], i_deg=BASE["i_deg"],
        raan_deg=BASE["raan_deg"], argp_deg=BASE["argp_deg"], m_deg=M0_DEG,
    )
    FIG_DIR.mkdir(exist_ok=True)
    figure_rate_sensitivity()
    figure_ground_track_comparison(elements0)
    figure_element_drift(elements0)


if __name__ == "__main__":
    main()
