"""M3 verification report generator.

Runs the frame/ground-track/access checks and prints a numeric summary
used to populate the M3 section of DESIGN.md.

Run:
    python scripts/m3_verification_report.py
"""

from __future__ import annotations

import numpy as np

from molniya_design.access import (
    access_metrics,
    elevation_independent_check,
    find_access_intervals,
    range_az_el_series,
)
from molniya_design.constants import (
    BASELINE_ELEMENTS_DEG,
    M1_PERIOD_S,
    MIN_ELEVATION_DEG,
    SIDEREAL_DAY_S,
    SITE_LAT_DEG,
    SITE_LON_DEG,
)
from molniya_design.elements import coe_to_rv
from molniya_design.frames import EARTH_ROTATION_RATE_RAD_S
from molniya_design.groundtrack import apsis_ground_points, compute_ground_track, longitude_separation_deg

BASE = BASELINE_ELEMENTS_DEG


def main():
    print("=" * 78)
    print("M3 VERIFICATION REPORT")
    print("=" * 78)

    r0, v0 = coe_to_rv(
        BASE["a_km"], BASE["e"], BASE["i_deg"], BASE["raan_deg"],
        BASE["argp_deg"], BASE["nu0_deg"],
    )

    print(f"\nEarth rotation rate omega_E = {EARTH_ROTATION_RATE_RAD_S:.10e} rad/s")
    print(f"omega_E * T (baseline period) = {np.degrees(EARTH_ROTATION_RATE_RAD_S*M1_PERIOD_S):.10f} deg (expect 180)")
    print(f"omega_E * sidereal_day = {np.degrees(EARTH_ROTATION_RATE_RAD_S*SIDEREAL_DAY_S):.10f} deg (expect 360)")

    # --- apsis ground points ---
    t_span = (0.0, 2.05 * M1_PERIOD_S)
    pts = apsis_ground_points(r0, v0, t_span)
    print(f"\nApsis ground points over 2.05 periods:")
    for p in pts:
        print(f"  t={p.t_s:12.3f} s  {p.kind:8s} lat={p.lat_deg:10.5f} deg  lon={p.lon_deg:10.5f} deg  r={p.r_km:.3f} km")

    apogees = [p for p in pts if p.kind == "apogee"]
    print(f"\nSuccessive apogee longitude separations:")
    for k in range(1, len(apogees)):
        sep = longitude_separation_deg(apogees[k - 1].lon_deg, apogees[k].lon_deg)
        print(f"  apogee[{k-1}] -> apogee[{k}]: {sep:.8f} deg (expect ~180 or ~-180)")

    if len(apogees) > 2:
        repeat_err = longitude_separation_deg(apogees[0].lon_deg, apogees[2].lon_deg)
        print(f"\nOne-sidereal-day repeat error (apogee[0] vs apogee[2]): {repeat_err:.3e} deg")
        print(f"  time diff = {apogees[2].t_s - apogees[0].t_s:.6f} s vs sidereal day {SIDEREAL_DAY_S:.6f} s")

    # --- site definition ---
    print(f"\nRepresentative site: lat={SITE_LAT_DEG} deg N, lon={SITE_LON_DEG} deg E, "
          f"min elevation={MIN_ELEVATION_DEG} deg")

    # --- ground track + access over one sidereal day, dt=15s (tightest) ---
    dt = 15.0
    t_span_day = (0.0, SIDEREAL_DAY_S)
    gt = compute_ground_track(r0, v0, t_span_day, dt_s=dt)
    rng, az, el = range_az_el_series(gt.r_ecef, SITE_LAT_DEG, SITE_LON_DEG)
    intervals = find_access_intervals(gt.t_s, el, threshold_deg=MIN_ELEVATION_DEG)
    m = access_metrics(intervals, t_span_day[0], t_span_day[1])

    print(f"\nAccess geometry over one sidereal day (dt={dt}s, {len(gt.t_s)} samples):")
    print(f"  number of passes       = {m.num_passes}")
    print(f"  total access time      = {m.total_access_time_s:.3f} s ({m.total_access_time_s/3600:.4f} h)")
    print(f"  access fraction        = {m.access_fraction:.6f}")
    print(f"  longest pass           = {m.longest_pass_s:.3f} s ({m.longest_pass_s/3600:.4f} h)")
    print(f"  max no-access gap      = {m.max_gap_s:.3f} s ({m.max_gap_s/3600:.4f} h)")
    peak_el = max((iv.peak_elevation_deg for iv in intervals), default=float("nan"))
    print(f"  peak elevation (best pass) = {peak_el:.4f} deg")
    for k, iv in enumerate(intervals):
        print(f"    pass {k}: start={iv.start_s:10.2f}s end={iv.end_s:10.2f}s dur={iv.duration_s:8.2f}s "
              f"peak_el={iv.peak_elevation_deg:7.3f} deg @ t={iv.peak_time_s:10.2f}s")

    # --- independent elevation cross-check over the same trajectory ---
    max_el_err, max_rng_err = 0.0, 0.0
    for k, t in enumerate(gt.t_s):
        rng2, el2 = elevation_independent_check(gt.r_ecef[:, k], SITE_LAT_DEG, SITE_LON_DEG)
        max_el_err = max(max_el_err, abs(el[k] - el2))
        max_rng_err = max(max_rng_err, abs(rng[k] - rng2))
    print(f"\nIndependent elevation cross-check (ENU vs triangle-geometry), over {len(gt.t_s)} samples:")
    print(f"  max |elevation error| = {max_el_err:.3e} deg")
    print(f"  max |range error|     = {max_rng_err:.3e} km")

    # --- timestep convergence study ---
    print(f"\nTimestep convergence study (one sidereal day):")
    dts = [120.0, 60.0, 15.0]
    results = []
    for dtc in dts:
        gtc = compute_ground_track(r0, v0, t_span_day, dt_s=dtc)
        rngc, azc, elc = range_az_el_series(gtc.r_ecef, SITE_LAT_DEG, SITE_LON_DEG)
        ivsc = find_access_intervals(gtc.t_s, elc, threshold_deg=MIN_ELEVATION_DEG)
        mc = access_metrics(ivsc, t_span_day[0], t_span_day[1])
        results.append((dtc, mc, ivsc))
        print(f"  dt={dtc:6.1f}s  n_samples={len(gtc.t_s):6d}  total_access={mc.total_access_time_s:10.3f}s  "
              f"max_gap={mc.max_gap_s:10.3f}s  longest_pass={mc.longest_pass_s:10.3f}s  n_passes={mc.num_passes}")

    print(f"\n  Convergence deltas (|value(dt) - value(15s)|):")
    ref_total = results[-1][1].total_access_time_s
    ref_gap = results[-1][1].max_gap_s
    ref_longest = results[-1][1].longest_pass_s
    for dtc, mc, _ in results:
        print(f"    dt={dtc:6.1f}s  d_total={abs(mc.total_access_time_s-ref_total):.4f}s  "
              f"d_gap={abs(mc.max_gap_s-ref_gap):.4f}s  d_longest={abs(mc.longest_pass_s-ref_longest):.4f}s")

    print("\n" + "=" * 78)
    print("END M3 VERIFICATION REPORT")
    print("=" * 78)


if __name__ == "__main__":
    main()
