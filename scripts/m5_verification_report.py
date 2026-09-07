"""M5 verification report generator: regional coverage, sensitivity
studies, convergence, and independent cross-checks.

Saves machine-readable results under results/ and prints a full summary
used to populate the M5 section of DESIGN.md.

Run:
    python scripts/m5_verification_report.py
"""

from __future__ import annotations

import csv
import json
import time
from pathlib import Path

import numpy as np

from molniya_design.access import elevation_independent_check, range_az_el
from molniya_design.constants import (
    BASELINE_ELEMENTS_DEG,
    MIN_ELEVATION_DEG,
    R_EARTH,
    SITE_LAT_DEG,
    SITE_LON_DEG,
)
from molniya_design.coverage import (
    build_satellite_ecef_trajectory,
    point_coverage,
    regional_coverage,
    summarize_regional,
)
from molniya_design.j2 import SecularElements, compute_secular_rates
from molniya_design.propagation import true_to_mean_anomaly

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
BASE = BASELINE_ELEMENTS_DEG
M0_DEG = float(np.degrees(true_to_mean_anomaly(np.radians(BASE["nu0_deg"]), BASE["e"])))

LAT_MIN, LAT_MAX = 60.0, 75.0


def baseline_elements(i_deg=None, argp_deg=None, raan_deg=None, a_km=None, e=None):
    return SecularElements(
        a_km=a_km if a_km is not None else BASE["a_km"],
        e=e if e is not None else BASE["e"],
        i_deg=i_deg if i_deg is not None else BASE["i_deg"],
        raan_deg=raan_deg if raan_deg is not None else BASE["raan_deg"],
        argp_deg=argp_deg if argp_deg is not None else BASE["argp_deg"],
        m_deg=M0_DEG,
    )


def regional_summary_dict(summ):
    return {
        "point_weighted_mean_access_fraction": summ.point_weighted_mean_access_fraction,
        "area_weighted_mean_access_fraction": summ.area_weighted_mean_access_fraction,
        "min_access_fraction": summ.min_access_fraction,
        "max_access_fraction": summ.max_access_fraction,
        "worst_max_gap_s": summ.worst_max_gap_s,
        "worst_max_gap_h": summ.worst_max_gap_s / 3600.0,
        "median_max_gap_s": summ.median_max_gap_s,
        "worst_gap_point_lat_deg": summ.worst_gap_point[0],
        "worst_gap_point_lon_deg": summ.worst_gap_point[1],
        "best_access_point_lat_deg": summ.best_access_point[0],
        "best_access_point_lon_deg": summ.best_access_point[1],
    }


def main():
    report = {}
    print("=" * 78)
    print("M5 VERIFICATION REPORT")
    print("=" * 78)

    # -----------------------------------------------------------------
    # 1. Baseline regional coverage (authoritative: 14 days, 2.5 deg grid)
    # -----------------------------------------------------------------
    print("\n--- Baseline regional coverage (14 days, dt=60s, 2.5 deg grid) ---")
    t0 = time.time()
    e0 = baseline_elements()
    t_span_14 = (0.0, 14 * 86400.0)
    results_14 = regional_coverage(e0, t_span_14, 60.0, LAT_MIN, LAT_MAX, 2.5, 2.5, MIN_ELEVATION_DEG)
    summ_14 = summarize_regional(results_14)
    print(f"  (computed in {time.time()-t0:.2f} s, {len(results_14)} grid points)")
    print(f"  point-weighted mean access fraction = {summ_14.point_weighted_mean_access_fraction:.6f}")
    print(f"  area-weighted mean access fraction  = {summ_14.area_weighted_mean_access_fraction:.6f}")
    print(f"  min/max access fraction = {summ_14.min_access_fraction:.6f} / {summ_14.max_access_fraction:.6f}")
    print(f"  worst max gap = {summ_14.worst_max_gap_s:.3f} s ({summ_14.worst_max_gap_s/3600:.4f} h) at {summ_14.worst_gap_point}")
    print(f"  median max gap = {summ_14.median_max_gap_s:.3f} s ({summ_14.median_max_gap_s/3600:.4f} h)")
    print(f"  best access point = {summ_14.best_access_point}")
    report["baseline_14day"] = regional_summary_dict(summ_14)

    # save full grid metrics
    with open(RESULTS_DIR / "m5_grid_metrics.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["lat_deg", "lon_deg", "access_fraction", "max_gap_s", "num_complete_passes",
                    "longest_pass_s", "peak_elevation_deg"])
        for r in results_14:
            w.writerow([r.lat_deg, r.lon_deg, r.access_fraction, r.max_gap_s, r.num_complete_passes,
                        r.longest_pass_s, r.peak_elevation_deg])
    print(f"  wrote results/m5_grid_metrics.csv ({len(results_14)} rows)")

    # representative site (14-day, non-periodic boundary handling)
    n14 = int(round((t_span_14[1] - t_span_14[0]) / 60.0)) + 1
    t_s_14 = np.linspace(t_span_14[0], t_span_14[1], n14)
    traj_14 = build_satellite_ecef_trajectory(e0, t_s_14)
    site_14 = point_coverage(traj_14, t_s_14, SITE_LAT_DEG, SITE_LON_DEG, MIN_ELEVATION_DEG)
    print(f"\n  Representative site (65N,40E), 14 days: access_fraction={site_14.access_fraction:.6f}, "
          f"max_gap={site_14.max_gap_s:.3f}s ({site_14.max_gap_s/3600:.4f}h), "
          f"longest_pass={site_14.longest_pass_s:.3f}s ({site_14.longest_pass_s/3600:.4f}h), "
          f"complete_passes={site_14.num_complete_passes}, leading_partial={site_14.has_leading_partial}, "
          f"trailing_partial={site_14.has_trailing_partial}, peak_el={site_14.peak_elevation_deg:.4f} deg")
    report["representative_site_14day"] = {
        "lat_deg": SITE_LAT_DEG, "lon_deg": SITE_LON_DEG,
        "access_fraction": site_14.access_fraction, "max_gap_s": site_14.max_gap_s,
        "longest_pass_s": site_14.longest_pass_s, "num_complete_passes": site_14.num_complete_passes,
        "has_leading_partial": site_14.has_leading_partial, "has_trailing_partial": site_14.has_trailing_partial,
        "peak_elevation_deg": site_14.peak_elevation_deg,
    }

    # -----------------------------------------------------------------
    # 2. Horizon stability (1/3/7/14 days)
    # -----------------------------------------------------------------
    print("\n--- Horizon stability (5 deg grid for speed, dt=120s) ---")
    horizon_rows = []
    for days in [1, 3, 7, 14]:
        t_span = (0.0, days * 86400.0)
        results = regional_coverage(e0, t_span, 120.0, LAT_MIN, LAT_MAX, 5.0, 5.0, MIN_ELEVATION_DEG)
        summ = summarize_regional(results)
        print(f"  {days:2d}d: worst_gap={summ.worst_max_gap_s:9.3f}s ({summ.worst_max_gap_s/3600:.4f}h) "
              f"at {summ.worst_gap_point}  ptmean={summ.point_weighted_mean_access_fraction:.4f} "
              f"areamean={summ.area_weighted_mean_access_fraction:.4f}")
        horizon_rows.append({"days": days, **regional_summary_dict(summ)})
    report["horizon_stability"] = horizon_rows

    # -----------------------------------------------------------------
    # 3. RAAN / longitude dependence
    # -----------------------------------------------------------------
    print("\n--- RAAN dependence: regional aggregate (should be ~invariant) vs site-specific (should NOT be) ---")
    raan_rows = []
    t_span_raan = (0.0, 7 * 86400.0)
    n_raan = int(round((t_span_raan[1] - t_span_raan[0]) / 120.0)) + 1
    t_s_raan = np.linspace(t_span_raan[0], t_span_raan[1], n_raan)
    for raan in [0.0, 60.0, 120.0, 180.0]:
        eR = baseline_elements(raan_deg=raan)
        results = regional_coverage(eR, t_span_raan, 120.0, LAT_MIN, LAT_MAX, 5.0, 10.0, MIN_ELEVATION_DEG)
        summ = summarize_regional(results)
        traj = build_satellite_ecef_trajectory(eR, t_s_raan)
        site = point_coverage(traj, t_s_raan, SITE_LAT_DEG, SITE_LON_DEG, MIN_ELEVATION_DEG)
        print(f"  RAAN={raan:6.1f}: regional ptmean={summ.point_weighted_mean_access_fraction:.6f}  "
              f"site_access_fraction={site.access_fraction:.6f}  site_max_gap={site.max_gap_s:.3f}s")
        raan_rows.append({
            "raan_deg": raan,
            "regional_point_weighted_mean_access_fraction": summ.point_weighted_mean_access_fraction,
            "site_access_fraction": site.access_fraction,
            "site_max_gap_s": site.max_gap_s,
        })
    report["raan_dependence"] = raan_rows

    # -----------------------------------------------------------------
    # 4. Minimum-elevation sensitivity
    # -----------------------------------------------------------------
    print("\n--- Minimum-elevation sensitivity (7 days, 5 deg grid) ---")
    elev_rows = []
    t_span_7 = (0.0, 7 * 86400.0)
    n7 = int(round((t_span_7[1] - t_span_7[0]) / 120.0)) + 1
    t_s_7 = np.linspace(t_span_7[0], t_span_7[1], n7)
    traj_7 = build_satellite_ecef_trajectory(e0, t_s_7)
    for elev_min in [5.0, 10.0, 15.0, 20.0]:
        results = regional_coverage(e0, t_span_7, 120.0, LAT_MIN, LAT_MAX, 5.0, 5.0, elev_min)
        summ = summarize_regional(results)
        site = point_coverage(traj_7, t_s_7, SITE_LAT_DEG, SITE_LON_DEG, elev_min)
        print(f"  el_min={elev_min:5.1f} deg: regional_mean_access={summ.point_weighted_mean_access_fraction:.6f}  "
              f"worst_gap={summ.worst_max_gap_s:.3f}s  site_access={site.access_fraction:.6f}  "
              f"site_longest_pass={site.longest_pass_s:.3f}s")
        elev_rows.append({
            "elevation_min_deg": elev_min,
            "regional_mean_access_fraction": summ.point_weighted_mean_access_fraction,
            "worst_max_gap_s": summ.worst_max_gap_s,
            "site_access_fraction": site.access_fraction,
            "site_longest_pass_s": site.longest_pass_s,
        })
    report["elevation_sensitivity"] = elev_rows

    # -----------------------------------------------------------------
    # 5. Inclination sensitivity
    # -----------------------------------------------------------------
    print("\n--- Inclination sensitivity (7 days, 5 deg grid, M4 J2 secular model) ---")
    incl_rows = []
    for i_deg in [60.0, 62.0, BASE["i_deg"], 65.0, 70.0]:
        eI = baseline_elements(i_deg=i_deg)
        results = regional_coverage(eI, t_span_7, 120.0, LAT_MIN, LAT_MAX, 5.0, 5.0, MIN_ELEVATION_DEG)
        summ = summarize_regional(results)
        traj = build_satellite_ecef_trajectory(eI, t_s_7)
        site = point_coverage(traj, t_s_7, SITE_LAT_DEG, SITE_LON_DEG, MIN_ELEVATION_DEG)
        rates = compute_secular_rates(eI.a_km, eI.e, eI.i_deg)
        print(f"  i={i_deg:9.4f} deg: regional_mean_access={summ.point_weighted_mean_access_fraction:.6f}  "
              f"worst_gap={summ.worst_max_gap_s:9.3f}s  site_access={site.access_fraction:.6f}  "
              f"argp_dot={rates.argp_dot_deg_day: .6f} deg/day")
        incl_rows.append({
            "i_deg": i_deg,
            "regional_mean_access_fraction": summ.point_weighted_mean_access_fraction,
            "worst_max_gap_s": summ.worst_max_gap_s,
            "site_access_fraction": site.access_fraction,
            "argp_dot_deg_day": rates.argp_dot_deg_day,
        })
    report["inclination_sensitivity"] = incl_rows

    # -----------------------------------------------------------------
    # 6. Argument-of-perigee sensitivity
    # -----------------------------------------------------------------
    print("\n--- Argument-of-perigee sensitivity (7 days, 5 deg grid) ---")
    argp_rows = []
    for argp_deg in [240.0, 255.0, 270.0, 285.0, 300.0, 90.0]:
        eW = baseline_elements(argp_deg=argp_deg)
        results = regional_coverage(eW, t_span_7, 120.0, LAT_MIN, LAT_MAX, 5.0, 5.0, MIN_ELEVATION_DEG)
        summ = summarize_regional(results)
        traj = build_satellite_ecef_trajectory(eW, t_s_7)
        site = point_coverage(traj, t_s_7, SITE_LAT_DEG, SITE_LON_DEG, MIN_ELEVATION_DEG)
        print(f"  argp={argp_deg:6.1f} deg: regional_mean_access={summ.point_weighted_mean_access_fraction:.6f}  "
              f"worst_gap={summ.worst_max_gap_s:9.3f}s  site_access={site.access_fraction:.6f}  "
              f"site_longest_pass={site.longest_pass_s:.3f}s")
        argp_rows.append({
            "argp_deg": argp_deg,
            "regional_mean_access_fraction": summ.point_weighted_mean_access_fraction,
            "worst_max_gap_s": summ.worst_max_gap_s,
            "site_access_fraction": site.access_fraction,
            "site_longest_pass_s": site.longest_pass_s,
        })
    report["argp_sensitivity"] = argp_rows

    # -----------------------------------------------------------------
    # 7. Eccentricity / perigee-altitude sensitivity (a fixed)
    # -----------------------------------------------------------------
    print("\n--- Perigee-altitude sensitivity (a fixed at half-sidereal-day value, 7 days, 5 deg grid) ---")
    perigee_rows = []
    for perigee_alt_km in [300.0, 600.0, 1000.0, 2000.0]:
        rp = R_EARTH + perigee_alt_km
        e_ecc = 1.0 - rp / BASE["a_km"]
        ra = BASE["a_km"] * (1.0 + e_ecc)
        eP = baseline_elements(e=e_ecc)
        results = regional_coverage(eP, t_span_7, 120.0, LAT_MIN, LAT_MAX, 5.0, 5.0, MIN_ELEVATION_DEG)
        summ = summarize_regional(results)
        traj = build_satellite_ecef_trajectory(eP, t_s_7)
        site = point_coverage(traj, t_s_7, SITE_LAT_DEG, SITE_LON_DEG, MIN_ELEVATION_DEG)
        rates = compute_secular_rates(eP.a_km, eP.e, eP.i_deg)
        print(f"  hp={perigee_alt_km:6.1f} km (e={e_ecc:.6f}, ha={ra-R_EARTH:9.2f} km): "
              f"regional_mean_access={summ.point_weighted_mean_access_fraction:.6f}  "
              f"worst_gap={summ.worst_max_gap_s:9.3f}s  site_access={site.access_fraction:.6f}  "
              f"RAAN_dot={rates.raan_dot_deg_day: .6f} deg/day")
        perigee_rows.append({
            "perigee_alt_km": perigee_alt_km, "eccentricity": e_ecc, "apogee_alt_km": ra - R_EARTH,
            "regional_mean_access_fraction": summ.point_weighted_mean_access_fraction,
            "worst_max_gap_s": summ.worst_max_gap_s,
            "site_access_fraction": site.access_fraction,
            "raan_dot_deg_day": rates.raan_dot_deg_day,
        })
    report["perigee_sensitivity"] = perigee_rows

    # -----------------------------------------------------------------
    # 8. Dwell-time cross-check: M1 proxy vs M3 site vs M5 regional
    # -----------------------------------------------------------------
    print("\n--- Dwell-time cross-check ---")
    print(f"  M1 anomaly-based dwell proxy (+/-60 deg true anomaly of apogee): 84.01% of orbital period")
    print(f"  M3 single-site (65N,40E) access fraction, 1 sidereal day, two-body: 80.30%")
    print(f"  M5 single-site (65N,40E) access fraction, 14 days, J2 secular: {site_14.access_fraction*100:.2f}%")
    print(f"  M5 regional point-weighted mean access fraction, 14 days: {summ_14.point_weighted_mean_access_fraction*100:.2f}%")
    report["dwell_comparison"] = {
        "m1_anomaly_dwell_proxy_fraction": 0.8401,
        "m3_site_access_fraction_1day_twobody": 0.803027,
        "m5_site_access_fraction_14day_j2": site_14.access_fraction,
        "m5_regional_point_weighted_mean_14day": summ_14.point_weighted_mean_access_fraction,
    }

    # -----------------------------------------------------------------
    # 9. Time convergence
    # -----------------------------------------------------------------
    print("\n--- Time-step convergence (representative site, 3 days) ---")
    t_span_3 = (0.0, 3 * 86400.0)
    conv_time_rows = []
    for dt in [120.0, 60.0, 30.0, 15.0]:
        n = int(round((t_span_3[1] - t_span_3[0]) / dt)) + 1
        t_s = np.linspace(t_span_3[0], t_span_3[1], n)
        traj = build_satellite_ecef_trajectory(e0, t_s)
        site = point_coverage(traj, t_s, SITE_LAT_DEG, SITE_LON_DEG, MIN_ELEVATION_DEG)
        results = regional_coverage(e0, t_span_3, dt, LAT_MIN, LAT_MAX, 5.0, 5.0, MIN_ELEVATION_DEG)
        summ = summarize_regional(results)
        print(f"  dt={dt:5.1f}s: site_access={site.access_fraction:.6f}  site_longest_pass={site.longest_pass_s:.3f}s  "
              f"regional_worst_gap={summ.worst_max_gap_s:.3f}s")
        conv_time_rows.append({
            "dt_s": dt, "site_access_fraction": site.access_fraction,
            "site_longest_pass_s": site.longest_pass_s, "regional_worst_gap_s": summ.worst_max_gap_s,
        })
    report["time_convergence"] = conv_time_rows

    # -----------------------------------------------------------------
    # 10. Spatial-grid convergence
    # -----------------------------------------------------------------
    print("\n--- Spatial-grid convergence (3 days, dt=120s) ---")
    conv_grid_rows = []
    for grid_step in [5.0, 2.5, 1.0]:
        t0g = time.time()
        results = regional_coverage(e0, t_span_3, 120.0, LAT_MIN, LAT_MAX, grid_step, grid_step, MIN_ELEVATION_DEG)
        summ = summarize_regional(results)
        print(f"  grid={grid_step:4.1f} deg ({len(results):5d} pts, {time.time()-t0g:.2f}s): "
              f"worst_gap={summ.worst_max_gap_s:9.3f}s at {summ.worst_gap_point}  "
              f"ptmean={summ.point_weighted_mean_access_fraction:.6f}")
        conv_grid_rows.append({
            "grid_step_deg": grid_step, "n_points": len(results),
            "worst_max_gap_s": summ.worst_max_gap_s,
            "worst_gap_point_lat_deg": summ.worst_gap_point[0], "worst_gap_point_lon_deg": summ.worst_gap_point[1],
            "point_weighted_mean_access_fraction": summ.point_weighted_mean_access_fraction,
        })
    report["spatial_convergence"] = conv_grid_rows

    # -----------------------------------------------------------------
    # 11. Independent elevation cross-check at 3 points
    # -----------------------------------------------------------------
    print("\n--- Independent elevation cross-check (ENU vs triangle-geometry) ---")
    cross_check_rows = []
    check_points = [
        ("representative_site", SITE_LAT_DEG, SITE_LON_DEG),
        ("worst_regional_point", summ_14.worst_gap_point[0], summ_14.worst_gap_point[1]),
        ("band_edge_75N", 75.0, 180.0),
    ]
    n2 = int(round((t_span_3[1] - t_span_3[0]) / 300.0)) + 1
    t_s_check = np.linspace(t_span_3[0], t_span_3[1], n2)
    traj_check = build_satellite_ecef_trajectory(e0, t_s_check)
    for label, lat, lon in check_points:
        max_err = 0.0
        for k in range(traj_check.shape[0]):
            _, _, el1 = range_az_el(traj_check[k], lat, lon)
            _, el2 = elevation_independent_check(traj_check[k], lat, lon)
            max_err = max(max_err, abs(el1 - el2))
        print(f"  {label} ({lat}N,{lon}E): max |elevation error| = {max_err:.3e} deg")
        cross_check_rows.append({"label": label, "lat_deg": lat, "lon_deg": lon, "max_elevation_error_deg": max_err})
    report["independent_crosscheck"] = cross_check_rows

    # -----------------------------------------------------------------
    # save results
    # -----------------------------------------------------------------
    report["baseline_orbit_elements"] = {
        "a_km": BASE["a_km"], "e": BASE["e"], "i_deg": BASE["i_deg"],
        "raan_deg": BASE["raan_deg"], "argp_deg": BASE["argp_deg"], "m0_deg": M0_DEG,
    }
    report["grid_definition"] = {"lat_min_deg": LAT_MIN, "lat_max_deg": LAT_MAX, "lat_step_deg": 2.5, "lon_step_deg": 2.5}
    report["elevation_min_deg"] = MIN_ELEVATION_DEG
    report["propagation_horizon_days"] = 14

    with open(RESULTS_DIR / "m5_baseline_regional_metrics.json", "w") as f:
        json.dump(report, f, indent=2, default=float)
    print(f"\nwrote {RESULTS_DIR / 'm5_baseline_regional_metrics.json'}")

    # sensitivity CSV (combined)
    with open(RESULTS_DIR / "m5_sensitivity.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["study", "parameter", "value", "regional_mean_access_fraction", "worst_max_gap_s", "site_access_fraction"])
        for row in elev_rows:
            w.writerow(["elevation_min", "elevation_min_deg", row["elevation_min_deg"],
                        row["regional_mean_access_fraction"], row["worst_max_gap_s"], row["site_access_fraction"]])
        for row in incl_rows:
            w.writerow(["inclination", "i_deg", row["i_deg"],
                        row["regional_mean_access_fraction"], row["worst_max_gap_s"], row["site_access_fraction"]])
        for row in argp_rows:
            w.writerow(["argp", "argp_deg", row["argp_deg"],
                        row["regional_mean_access_fraction"], row["worst_max_gap_s"], row["site_access_fraction"]])
        for row in perigee_rows:
            w.writerow(["perigee_altitude", "perigee_alt_km", row["perigee_alt_km"],
                        row["regional_mean_access_fraction"], row["worst_max_gap_s"], row["site_access_fraction"]])
    print(f"wrote {RESULTS_DIR / 'm5_sensitivity.csv'}")

    with open(RESULTS_DIR / "m5_convergence.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["study", "parameter", "value", "metric", "result"])
        for row in conv_time_rows:
            w.writerow(["time", "dt_s", row["dt_s"], "site_access_fraction", row["site_access_fraction"]])
            w.writerow(["time", "dt_s", row["dt_s"], "regional_worst_gap_s", row["regional_worst_gap_s"]])
        for row in conv_grid_rows:
            w.writerow(["spatial", "grid_step_deg", row["grid_step_deg"], "worst_max_gap_s", row["worst_max_gap_s"]])
            w.writerow(["spatial", "grid_step_deg", row["grid_step_deg"], "point_weighted_mean_access_fraction", row["point_weighted_mean_access_fraction"]])
    print(f"wrote {RESULTS_DIR / 'm5_convergence.csv'}")

    print("\n" + "=" * 78)
    print("END M5 VERIFICATION REPORT")
    print("=" * 78)


if __name__ == "__main__":
    main()
