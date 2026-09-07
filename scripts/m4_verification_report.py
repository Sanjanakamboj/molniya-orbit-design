"""M4 verification report generator.

Runs the J2 secular rate / ground-track / Cartesian cross-check analysis
and prints a numeric summary used to populate the M4 section of
DESIGN.md.

Run:
    python scripts/m4_verification_report.py
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from molniya_design.constants import BASELINE_ELEMENTS_DEG, SIDEREAL_DAY_S
from molniya_design.elements import coe_to_rv, rv_to_coe
from molniya_design.groundtrack import (
    apogee_ground_points_j2,
    compute_ground_track,
    compute_ground_track_j2,
    longitude_separation_deg,
)
from molniya_design.j2 import SecularElements, compute_secular_rates, propagate_elements_j2_secular
from molniya_design.j2_cartesian import j2_eom
from molniya_design.propagation import true_to_mean_anomaly

BASE = BASELINE_ELEMENTS_DEG
M0_DEG = float(np.degrees(true_to_mean_anomaly(np.radians(BASE["nu0_deg"]), BASE["e"])))


def main():
    print("=" * 78)
    print("M4 VERIFICATION REPORT")
    print("=" * 78)

    elements0 = SecularElements(
        a_km=BASE["a_km"], e=BASE["e"], i_deg=BASE["i_deg"],
        raan_deg=BASE["raan_deg"], argp_deg=BASE["argp_deg"], m_deg=M0_DEG,
    )

    rates = compute_secular_rates(elements0.a_km, elements0.e, elements0.i_deg)
    print(f"\nBaseline secular rates (i={elements0.i_deg:.6f} deg):")
    print(f"  n        = {rates.n_rad_s:.10e} rad/s")
    print(f"  RAAN_dot = {rates.raan_dot_rad_s:.10e} rad/s = {rates.raan_dot_deg_day:.6f} deg/day")
    print(f"  argp_dot = {rates.argp_dot_rad_s:.10e} rad/s = {rates.argp_dot_deg_day:.6e} deg/day")
    print(f"  M_dot    = {rates.m_dot_rad_s:.10e} rad/s = {rates.m_dot_deg_day:.6f} deg/day")
    print(f"  M1 target RAAN_dot = -0.145135 deg/day")
    print(f"  5cos^2(i)-1 = {5*np.cos(np.radians(elements0.i_deg))**2-1:.3e}")

    print(f"\nOff-critical inclination comparison:")
    for i_deg in [60.0, 63.0, elements0.i_deg, 64.0, 65.0]:
        r = compute_secular_rates(elements0.a_km, elements0.e, i_deg)
        print(f"  i={i_deg:12.6f} deg  5cos2i-1={5*np.cos(np.radians(i_deg))**2-1: .6f}  "
              f"argp_dot={r.argp_dot_deg_day: .6f} deg/day  RAAN_dot={r.raan_dot_deg_day: .6f} deg/day")

    # --- ground track drift at multiple horizons ---
    print(f"\nApogee ground-track drift under secular J2 (same-side apogee comparison):")
    for days in [1, 3, 7, 14]:
        t_span = (0.0, days * 86400.0)
        pts = apogee_ground_points_j2(elements0, t_span)
        first, last = pts[0], pts[-1]
        sep = longitude_separation_deg(first.lon_deg, last.lon_deg)
        lat_drift = last.lat_deg - first.lat_deg
        print(f"  {days:2d}d: n_apogees={len(pts):3d}  first=(lat={first.lat_deg:.5f}, lon={first.lon_deg:.5f})  "
              f"last=(lat={last.lat_deg:.5f}, lon={last.lon_deg:.5f})  lon_drift={sep:.5f} deg  lat_drift={lat_drift:.3e} deg")

    # decomposition of the 1-day drift
    pts1 = apogee_ground_points_j2(elements0, (0.0, 86400.0))
    dt = pts1[-1].t_s
    raan_component = rates.raan_dot_deg_day * dt / 86400.0
    extra_rot_deg = (dt - SIDEREAL_DAY_S) / SIDEREAL_DAY_S * 360.0
    print(f"\n1-day apogee-drift decomposition (dt={dt:.3f} s to same-side apogee):")
    print(f"  RAAN-regression component = {raan_component:.6f} deg")
    print(f"  extra-Earth-rotation component (dt != exactly 1 sidereal day) = {-extra_rot_deg:.6f} deg")
    print(f"  sum = {raan_component - extra_rot_deg:.6f} deg  (observed = {longitude_separation_deg(pts1[0].lon_deg, pts1[-1].lon_deg):.6f} deg)")

    # --- one-day repeat error under J2 (should NOT repeat exactly, unlike M3 two-body) ---
    gt_j2_1day = compute_ground_track_j2(elements0, (0.0, SIDEREAL_DAY_S), dt_s=60.0)
    lon_start, lon_end = gt_j2_1day.lon_deg[0], gt_j2_1day.lon_deg[-1]
    print(f"\nOne-sidereal-day repeat error under J2 (start vs end of window, NOT an apogee-to-apogee comparison):")
    print(f"  lon(t=0) = {lon_start:.6f} deg, lon(t=sidereal_day) = {lon_end:.6f} deg, "
          f"diff = {longitude_separation_deg(lon_start, lon_end):.6f} deg")

    # --- omega stability over 14 days: baseline vs off-critical ---
    print(f"\nArgument-of-perigee stability over 14 days:")
    for i_deg, label in [(elements0.i_deg, "baseline (critical)"), (65.0, "off-critical i=65")]:
        e0 = SecularElements(a_km=elements0.a_km, e=elements0.e, i_deg=i_deg,
                              raan_deg=elements0.raan_deg, argp_deg=elements0.argp_deg, m_deg=elements0.m_deg)
        e14 = propagate_elements_j2_secular(e0, 14 * 86400.0, wrap=False)
        r_apo, _ = coe_to_rv(e14.a_km, e14.e, e14.i_deg, e14.raan_deg % 360.0, e14.argp_deg % 360.0, 180.0)
        from molniya_design.propagation import geocentric_latitude_deg
        lat_apo = geocentric_latitude_deg(r_apo)
        r_apo0, _ = coe_to_rv(e0.a_km, e0.e, e0.i_deg, e0.raan_deg, e0.argp_deg, 180.0)
        lat_apo0 = geocentric_latitude_deg(r_apo0)
        print(f"  {label}: omega(0)={e0.argp_deg:.6f} deg, omega(14d)={e14.argp_deg:.6f} deg, "
              f"delta={e14.argp_deg-e0.argp_deg:.6f} deg | apogee_lat(0)={lat_apo0:.6f} deg, "
              f"apogee_lat(14d)={lat_apo:.6f} deg, delta={lat_apo-lat_apo0:.6e} deg")

    # --- two-body vs J2 ground track divergence ---
    r0, v0 = coe_to_rv(BASE["a_km"], BASE["e"], BASE["i_deg"], BASE["raan_deg"], BASE["argp_deg"], BASE["nu0_deg"])
    print(f"\nTwo-body vs J2 ground-track divergence:")
    for days in [1, 3, 7, 14]:
        t_span = (0.0, days * 86400.0)
        gt_2body = compute_ground_track(r0, v0, t_span, dt_s=300.0)
        gt_j2 = compute_ground_track_j2(elements0, t_span, dt_s=300.0)
        lon_diff_end = longitude_separation_deg(gt_2body.lon_deg[-1], gt_j2.lon_deg[-1])
        print(f"  {days:2d}d: longitude diff at window end = {lon_diff_end:.6f} deg")

    # --- Cartesian J2 cross-check ---
    print(f"\nCartesian J2 cross-check (osculating propagation, 7 days, perigee-sampled fit):")
    y0 = np.concatenate([r0, v0])

    def rdotv(t, y):
        return np.dot(y[0:3], y[3:6])

    rdotv.direction = 0
    sol = solve_ivp(j2_eom, (0.0, 7 * 86400.0), y0, method="DOP853", rtol=1e-13, atol=1e-13,
                     dense_output=True, events=rdotv)
    event_times = sol.t_events[0]
    states = [sol.sol(t) for t in event_times]
    radii = [np.linalg.norm(s[0:3]) for s in states]
    r_mid = 0.5 * (min(radii) + max(radii))
    perigee_t, perigee_raan, perigee_argp = [], [], []
    for t, s, r in zip(event_times, states, radii):
        if r < r_mid:
            coe = rv_to_coe(s[0:3], s[3:6])
            perigee_t.append(t)
            perigee_raan.append(coe["raan_deg"])
            perigee_argp.append(coe["argp_deg"])
    perigee_t = np.array(perigee_t)
    perigee_raan_u = np.degrees(np.unwrap(np.radians(perigee_raan)))
    perigee_argp_u = np.degrees(np.unwrap(np.radians(perigee_argp)))
    print(f"  {len(perigee_t)} perigee samples over 7 days")
    A = np.vstack([perigee_t, np.ones_like(perigee_t)]).T
    raan_slope = np.linalg.lstsq(A, perigee_raan_u, rcond=None)[0][0] * 86400.0
    argp_slope = np.linalg.lstsq(A, perigee_argp_u, rcond=None)[0][0] * 86400.0
    print(f"  fitted RAAN_dot = {raan_slope:.6f} deg/day  (theory {rates.raan_dot_deg_day:.6f} deg/day, "
          f"rel err {abs(raan_slope-rates.raan_dot_deg_day)/abs(rates.raan_dot_deg_day):.3e})")
    print(f"  fitted argp_dot = {argp_slope:.6e} deg/day  (theory {rates.argp_dot_deg_day:.3e} deg/day)")

    # J2=0 reduction check
    from molniya_design.twobody import two_body_eom
    d1 = j2_eom(0.0, y0, J2=0.0)
    d2 = two_body_eom(0.0, y0)
    print(f"\n  J2=0 Cartesian EOM matches two-body EOM exactly: {np.allclose(d1, d2, atol=1e-15)}")

    print("\n" + "=" * 78)
    print("END M4 VERIFICATION REPORT")
    print("=" * 78)


if __name__ == "__main__":
    main()
