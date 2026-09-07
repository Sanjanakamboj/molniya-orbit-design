"""Generate the two M2 diagnostic figures:

- figures/m2_orbit_geometry.png : 3D ECI two-body trajectory, Earth,
  perigee/apogee markers. NOT a ground track (no Earth rotation applied).
- figures/m2_conservation_and_apsides.png : radius/altitude vs time with
  detected apsides, plus energy/angular-momentum relative error vs time.

Run:
    python scripts/m2_figures.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (registers 3D projection)

from molniya_design.constants import (
    BASELINE_ELEMENTS_DEG,
    M1_PERIOD_S,
    R_EARTH,
)
from molniya_design.elements import coe_to_rv
from molniya_design.propagation import find_apsides, propagate
from molniya_design.twobody import specific_angular_momentum, specific_energy

FIG_DIR = Path(__file__).resolve().parent.parent / "figures"
BASE = BASELINE_ELEMENTS_DEG


def figure_orbit_geometry(r0, v0):
    t_eval = np.linspace(0.0, M1_PERIOD_S, 2000)
    sol = propagate(r0, v0, (0.0, M1_PERIOD_S), t_eval=t_eval, rtol=1e-13, atol=1e-13)
    xs, ys, zs = sol.y[0], sol.y[1], sol.y[2]

    apsides = find_apsides(r0, v0, (0.0, M1_PERIOD_S * 1.01), rtol=1e-13, atol=1e-13)
    perigee = next(a for a in apsides if a.kind == "perigee")
    apogee = next(a for a in apsides if a.kind == "apogee" and a.t_s > 1.0)

    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111, projection="3d")

    # Earth sphere
    u, v = np.mgrid[0 : 2 * np.pi : 40j, 0 : np.pi : 20j]
    ex = R_EARTH * np.cos(u) * np.sin(v)
    ey = R_EARTH * np.sin(u) * np.sin(v)
    ez = R_EARTH * np.cos(v)
    ax.plot_surface(ex, ey, ez, color="steelblue", alpha=0.35, linewidth=0, zorder=1)

    ax.plot(xs, ys, zs, color="firebrick", lw=1.6, label="Two-body ECI trajectory (1 period)", zorder=3)
    ax.scatter(*r0, color="black", s=40, marker="o", label="Epoch state (apogee, nu0=180 deg)", zorder=5)
    ax.scatter(
        perigee.state[0], perigee.state[1], perigee.state[2],
        color="darkorange", s=60, marker="^", label=f"Perigee (r={perigee.r_km:.0f} km)", zorder=5,
    )
    ax.scatter(
        apogee.state[0], apogee.state[1], apogee.state[2],
        color="darkgreen", s=60, marker="s", label=f"Apogee (r={apogee.r_km:.0f} km)", zorder=5,
    )

    max_r = 1.05 * np.linalg.norm(apogee.state[0:3])
    ax.set_xlim(-max_r, max_r)
    ax.set_ylim(-max_r, max_r)
    ax.set_zlim(-max_r, max_r)
    ax.set_box_aspect([1, 1, 1])

    ax.set_xlabel("ECI X [km]")
    ax.set_ylabel("ECI Y [km]")
    ax.set_zlabel("ECI Z [km]")
    ax.set_title(
        "M2: Baseline Molniya orbit — inertial (ECI) two-body trajectory\n"
        "(Earth-centered inertial frame; this is NOT an Earth-fixed ground track — see M3)"
    )
    ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout()
    out = FIG_DIR / "m2_orbit_geometry.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"wrote {out}")


def figure_conservation_and_apsides(r0, v0):
    t_eval = np.linspace(0.0, 2.0 * M1_PERIOD_S, 3000)
    sol = propagate(r0, v0, (0.0, 2.0 * M1_PERIOD_S), t_eval=t_eval, rtol=1e-13, atol=1e-13)
    r_mag = np.linalg.norm(sol.y[0:3], axis=0)
    alt = r_mag - R_EARTH

    energies = np.array(
        [specific_energy(sol.y[0:3, k], sol.y[3:6, k]) for k in range(sol.y.shape[1])]
    )
    h_mags = np.array(
        [
            np.linalg.norm(specific_angular_momentum(sol.y[0:3, k], sol.y[3:6, k]))
            for k in range(sol.y.shape[1])
        ]
    )
    e_rel_err = np.abs(energies - energies[0]) / np.abs(energies[0])
    h_rel_err = np.abs(h_mags - h_mags[0]) / h_mags[0]
    e_rel_err = np.clip(e_rel_err, 1e-16, None)  # avoid log(0)
    h_rel_err = np.clip(h_rel_err, 1e-16, None)

    apsides = find_apsides(r0, v0, (0.0, 2.05 * M1_PERIOD_S), rtol=1e-13, atol=1e-13)

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=False)

    ax0 = axes[0]
    ax0.plot(sol.t / 3600.0, alt, color="steelblue", lw=1.3, label="Altitude (two-body propagation)")
    for a in apsides:
        color = "darkorange" if a.kind == "perigee" else "darkgreen"
        marker = "^" if a.kind == "perigee" else "s"
        ax0.scatter(a.t_s / 3600.0, a.r_km - R_EARTH, color=color, marker=marker, s=50, zorder=5)
    ax0.scatter([], [], color="darkorange", marker="^", label="Detected perigee (r.v=0 crossing)")
    ax0.scatter([], [], color="darkgreen", marker="s", label="Detected apogee (r.v=0 crossing)")
    ax0.set_xlabel("Time [hours]")
    ax0.set_ylabel("Altitude [km]")
    ax0.set_title("M2: Altitude vs time over 2 orbital periods, with detected apsides")
    ax0.legend(loc="upper right", fontsize=9)
    ax0.grid(alpha=0.3)

    ax1 = axes[1]
    ax1.semilogy(sol.t / 3600.0, e_rel_err, color="crimson", lw=1.2, label="|specific energy relative error|")
    ax1.semilogy(sol.t / 3600.0, h_rel_err, color="purple", lw=1.2, label="|angular momentum relative error|")
    ax1.set_xlabel("Time [hours]")
    ax1.set_ylabel("Relative error (log scale)")
    ax1.set_title("M2: Conservation-quantity relative error vs time (rtol=atol=1e-13)")
    ax1.legend(loc="upper right", fontsize=9)
    ax1.grid(alpha=0.3, which="both")

    fig.tight_layout()
    out = FIG_DIR / "m2_conservation_and_apsides.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"wrote {out}")


def main():
    r0, v0 = coe_to_rv(
        BASE["a_km"], BASE["e"], BASE["i_deg"], BASE["raan_deg"],
        BASE["argp_deg"], BASE["nu0_deg"],
    )
    FIG_DIR.mkdir(exist_ok=True)
    figure_orbit_geometry(r0, v0)
    figure_conservation_and_apsides(r0, v0)


if __name__ == "__main__":
    main()
