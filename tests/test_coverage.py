"""Milestone 5 test suite: regional coverage/revisit engine, window-
boundary handling, sensitivity behavior, and cross-checks.

Covers checklist items A-R from the M5 task description.
"""

from __future__ import annotations

import numpy as np
import pytest

from molniya_design.access import (
    AccessInterval,
    access_metrics_boundary_aware,
    elevation_independent_check,
    find_access_intervals,
    merge_cyclic_boundary_intervals,
    range_az_el,
    summarize_boundary_passes,
)
from molniya_design.constants import (
    BASELINE_ELEMENTS_DEG,
    M1_PERIOD_S,
    MIN_ELEVATION_DEG,
    R_EARTH,
    SIDEREAL_DAY_S,
    SITE_LAT_DEG,
    SITE_LON_DEG,
)
from molniya_design.coverage import (
    PointCoverageResult,
    build_grid,
    build_satellite_ecef_trajectory,
    point_coverage,
    regional_coverage,
    summarize_regional,
)
from molniya_design.elements import coe_to_r_array, coe_to_rv
from molniya_design.frames import eci_array_to_ecef_array, eci_to_ecef
from molniya_design.j2 import SecularElements, compute_secular_rates, propagate_elements_j2_secular
from molniya_design.propagation import true_to_mean_anomaly

BASE = BASELINE_ELEMENTS_DEG
M0_DEG = float(np.degrees(true_to_mean_anomaly(np.radians(BASE["nu0_deg"]), BASE["e"])))


@pytest.fixture(scope="module")
def elements0():
    return SecularElements(
        a_km=BASE["a_km"], e=BASE["e"], i_deg=BASE["i_deg"],
        raan_deg=BASE["raan_deg"], argp_deg=BASE["argp_deg"], m_deg=M0_DEG,
    )


# ---------------------------------------------------------------------------
# vectorized-primitive cross-checks (needed before trusting anything else)
# ---------------------------------------------------------------------------


def test_vectorized_coe_to_r_matches_scalar():
    raan_arr = np.array([0.0, 45.0, 270.0])
    argp_arr = np.array([270.0, 10.0, 300.0])
    nu_arr = np.array([0.0, 90.0, 180.0])
    r_vec = coe_to_r_array(BASE["a_km"], BASE["e"], BASE["i_deg"], raan_arr, argp_arr, nu_arr)
    for k in range(3):
        r_scalar, _ = coe_to_rv(BASE["a_km"], BASE["e"], BASE["i_deg"], raan_arr[k], argp_arr[k], nu_arr[k])
        assert np.allclose(r_vec[k], r_scalar, atol=1e-9)


def test_vectorized_ecef_matches_scalar():
    r_eci = np.array([[1000.0, 2000.0, 3000.0], [-4000.0, 5000.0, -6000.0]])
    t_s = np.array([0.0, 12345.0])
    r_ecef_vec = eci_array_to_ecef_array(r_eci, t_s)
    for k in range(2):
        r_ecef_scalar = eci_to_ecef(r_eci[k], t_s[k])
        assert np.allclose(r_ecef_vec[k], r_ecef_scalar, atol=1e-9)


# ---------------------------------------------------------------------------
# A. M3 access result reproduces in the J2=0 / 1-day limit
# ---------------------------------------------------------------------------


def test_A_j2_zero_reproduces_m3_access():
    """With J2=0, the coverage engine's regional-trajectory path should
    match M3's per-point access fraction at the representative site to
    high precision."""
    elements0_j2zero = SecularElements(
        a_km=BASE["a_km"], e=BASE["e"], i_deg=BASE["i_deg"],
        raan_deg=BASE["raan_deg"], argp_deg=BASE["argp_deg"], m_deg=M0_DEG,
    )
    t_span = (0.0, SIDEREAL_DAY_S)
    dt = 15.0
    n = int(round((t_span[1] - t_span[0]) / dt)) + 1
    t_s = np.linspace(t_span[0], t_span[1], n)
    traj = build_satellite_ecef_trajectory(elements0_j2zero, t_s, J2=0.0)
    res = point_coverage(traj, t_s, SITE_LAT_DEG, SITE_LON_DEG, MIN_ELEVATION_DEG, periodic=True, period_s=SIDEREAL_DAY_S)

    # M3 committed value: access_fraction = 0.803027 (see results/m3_verification_report.txt)
    assert res.access_fraction == pytest.approx(0.803027, abs=2e-4)


# ---------------------------------------------------------------------------
# B. threshold-crossing interpolation test
# ---------------------------------------------------------------------------


def test_B_threshold_crossing_interpolation():
    t = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    el = np.array([5.0, 9.0, 11.0, 8.0, 4.0])  # crosses 10 between idx1-2, and idx2-3
    intervals = find_access_intervals(t, el, threshold_deg=10.0)
    assert len(intervals) == 1
    # linear interp: between (1,9) and (2,11): frac=(10-9)/(11-9)=0.5 -> t=1.5
    assert intervals[0].start_s == pytest.approx(1.5, abs=1e-9)
    # between (2,11) and (3,8): frac=(10-11)/(8-11)=1/3 -> t=2+1/3
    assert intervals[0].end_s == pytest.approx(2.0 + 1.0 / 3.0, abs=1e-9)


# ---------------------------------------------------------------------------
# C. boundary-window interval merge test (periodic case)
# ---------------------------------------------------------------------------


def test_C_periodic_boundary_merge():
    # synthetic: access at both start and end of window, known to be cyclic
    t = np.linspace(0, 10, 11)
    el = np.array([15, 12, 5, 3, 2, 3, 5, 8, 12, 14, 16], dtype=float)  # above 10 at t=0,1 and t=8,9,10
    intervals = find_access_intervals(t, el, threshold_deg=10.0)
    assert len(intervals) == 2
    assert intervals[0].is_boundary_start
    assert intervals[-1].is_boundary_end

    merged = merge_cyclic_boundary_intervals(intervals, 0.0, 10.0, period_s=10.0)
    assert len(merged) == 1
    assert merged[0].duration_s == pytest.approx(intervals[0].duration_s + intervals[-1].duration_s)

    m = access_metrics_boundary_aware(intervals, 0.0, 10.0, periodic=True, period_s=10.0)
    assert m["num_passes"] == 1
    assert m["boundary_treatment"] == "merged_cyclic"


# ---------------------------------------------------------------------------
# D. non-periodic lead/trail gap handling test
# ---------------------------------------------------------------------------


def test_D_non_periodic_boundary_not_merged():
    t = np.linspace(0, 10, 11)
    el = np.array([15, 12, 5, 3, 2, 3, 5, 8, 12, 14, 16], dtype=float)
    intervals = find_access_intervals(t, el, threshold_deg=10.0)

    m = access_metrics_boundary_aware(intervals, 0.0, 10.0, periodic=False)
    # leading and trailing are NOT assumed to be the same pass
    assert m["num_passes"] == 0  # zero *complete* interior passes in this synthetic signal
    assert m["has_leading_partial"] is True
    assert m["has_trailing_partial"] is True
    assert m["boundary_treatment"] == "labeled_partial_lead_trail"
    # total access time / access fraction are unaffected by the merge question
    assert m["total_access_time_s"] == pytest.approx(intervals[0].duration_s + intervals[-1].duration_s)


def test_D_non_periodic_with_interior_complete_pass():
    t = np.linspace(0, 20, 21)
    el = np.array(
        [15, 12, 5, 3, 2, 3, 5, 3, 2, 3, 12, 15, 12, 3, 2, 3, 5, 8, 12, 14, 16], dtype=float
    )
    intervals = find_access_intervals(t, el, threshold_deg=10.0)
    assert len(intervals) == 3  # leading partial, one interior complete, trailing partial
    m = access_metrics_boundary_aware(intervals, 0.0, 20.0, periodic=False)
    assert m["num_passes"] == 1
    assert m["num_complete_passes"] == 1
    assert m["has_leading_partial"] and m["has_trailing_partial"]


# ---------------------------------------------------------------------------
# E. access fraction bounded [0,1]
# ---------------------------------------------------------------------------


def test_E_access_fraction_bounded(elements0):
    t_span = (0.0, 3 * 86400.0)
    dt = 300.0
    n = int(round((t_span[1] - t_span[0]) / dt)) + 1
    t_s = np.linspace(t_span[0], t_span[1], n)
    traj = build_satellite_ecef_trajectory(elements0, t_s)
    for lat, lon in [(60.0, 0.0), (75.0, 180.0), (65.0, 40.0), (60.0, 270.0)]:
        res = point_coverage(traj, t_s, lat, lon, MIN_ELEVATION_DEG)
        assert 0.0 <= res.access_fraction <= 1.0


# ---------------------------------------------------------------------------
# F. higher elevation threshold never increases access fraction
# ---------------------------------------------------------------------------


def test_F_higher_threshold_never_increases_access(elements0):
    t_span = (0.0, 3 * 86400.0)
    dt = 120.0
    n = int(round((t_span[1] - t_span[0]) / dt)) + 1
    t_s = np.linspace(t_span[0], t_span[1], n)
    traj = build_satellite_ecef_trajectory(elements0, t_s)

    thresholds = [5.0, 10.0, 15.0, 20.0]
    fractions = [point_coverage(traj, t_s, SITE_LAT_DEG, SITE_LON_DEG, th).access_fraction for th in thresholds]
    for k in range(1, len(fractions)):
        assert fractions[k] <= fractions[k - 1] + 1e-12


# ---------------------------------------------------------------------------
# G. higher elevation threshold does not reduce max gap (controlled case)
# ---------------------------------------------------------------------------


def test_G_higher_threshold_never_decreases_max_gap(elements0):
    t_span = (0.0, 3 * 86400.0)
    dt = 120.0
    n = int(round((t_span[1] - t_span[0]) / dt)) + 1
    t_s = np.linspace(t_span[0], t_span[1], n)
    traj = build_satellite_ecef_trajectory(elements0, t_s)

    thresholds = [5.0, 10.0, 15.0, 20.0]
    gaps = [point_coverage(traj, t_s, SITE_LAT_DEG, SITE_LON_DEG, th).max_gap_s for th in thresholds]
    for k in range(1, len(gaps)):
        assert gaps[k] >= gaps[k - 1] - 1e-6


# ---------------------------------------------------------------------------
# H. omega=270 outperforms omega=90 for northern target
# ---------------------------------------------------------------------------


def test_H_omega_270_beats_omega_90_northern():
    t_span = (0.0, 7 * 86400.0)
    dt = 120.0
    n = int(round((t_span[1] - t_span[0]) / dt)) + 1
    t_s = np.linspace(t_span[0], t_span[1], n)

    e_270 = SecularElements(a_km=BASE["a_km"], e=BASE["e"], i_deg=BASE["i_deg"],
                             raan_deg=BASE["raan_deg"], argp_deg=270.0, m_deg=M0_DEG)
    e_90 = SecularElements(a_km=BASE["a_km"], e=BASE["e"], i_deg=BASE["i_deg"],
                            raan_deg=BASE["raan_deg"], argp_deg=90.0, m_deg=M0_DEG)

    traj_270 = build_satellite_ecef_trajectory(e_270, t_s)
    traj_90 = build_satellite_ecef_trajectory(e_90, t_s)

    res_270 = point_coverage(traj_270, t_s, SITE_LAT_DEG, SITE_LON_DEG, MIN_ELEVATION_DEG)
    res_90 = point_coverage(traj_90, t_s, SITE_LAT_DEG, SITE_LON_DEG, MIN_ELEVATION_DEG)

    assert res_270.access_fraction > res_90.access_fraction
    assert res_270.max_gap_s < res_90.max_gap_s
    # dramatic difference expected
    assert res_270.access_fraction > 2.0 * max(res_90.access_fraction, 1e-6)


# ---------------------------------------------------------------------------
# I. critical-inclination omega remains fixed under M4 model
# ---------------------------------------------------------------------------


def test_I_critical_inclination_argp_fixed(elements0):
    e_14 = propagate_elements_j2_secular(elements0, 14 * 86400.0, wrap=False)
    assert e_14.argp_deg == pytest.approx(270.0, abs=1e-9)


# ---------------------------------------------------------------------------
# J. off-critical orientation changes with time
# ---------------------------------------------------------------------------


def test_J_off_critical_argp_changes():
    e0 = SecularElements(a_km=BASE["a_km"], e=BASE["e"], i_deg=65.0,
                          raan_deg=BASE["raan_deg"], argp_deg=270.0, m_deg=M0_DEG)
    e_14 = propagate_elements_j2_secular(e0, 14 * 86400.0, wrap=False)
    assert e_14.argp_deg != pytest.approx(270.0, abs=1e-6)


# ---------------------------------------------------------------------------
# K. RAAN rotation invariance of full-longitude regional aggregate
# ---------------------------------------------------------------------------


def test_K_raan_invariance_regional_aggregate():
    t_span = (0.0, 3 * 86400.0)
    dt = 300.0
    means = []
    for raan in [0.0, 90.0, 200.0]:
        e0 = SecularElements(a_km=BASE["a_km"], e=BASE["e"], i_deg=BASE["i_deg"],
                              raan_deg=raan, argp_deg=BASE["argp_deg"], m_deg=M0_DEG)
        results = regional_coverage(e0, t_span, dt, 60.0, 75.0, 5.0, 10.0, MIN_ELEVATION_DEG)
        summ = summarize_regional(results)
        means.append(summ.point_weighted_mean_access_fraction)
    assert means[0] == pytest.approx(means[1], abs=1e-6)
    assert means[0] == pytest.approx(means[2], abs=1e-6)


# ---------------------------------------------------------------------------
# L. site-specific RAAN timing shift behaves consistently (NOT invariant)
# ---------------------------------------------------------------------------


def test_L_site_specific_raan_dependence():
    t_span = (0.0, 3 * 86400.0)
    dt = 120.0
    n = int(round((t_span[1] - t_span[0]) / dt)) + 1
    t_s = np.linspace(t_span[0], t_span[1], n)

    fractions = []
    for raan in [0.0, 90.0]:
        e0 = SecularElements(a_km=BASE["a_km"], e=BASE["e"], i_deg=BASE["i_deg"],
                              raan_deg=raan, argp_deg=BASE["argp_deg"], m_deg=M0_DEG)
        traj = build_satellite_ecef_trajectory(e0, t_s)
        res = point_coverage(traj, t_s, SITE_LAT_DEG, SITE_LON_DEG, MIN_ELEVATION_DEG)
        fractions.append(res.access_fraction)
    # a fixed site's metrics should differ (not be invariant) under RAAN shift
    assert fractions[0] != pytest.approx(fractions[1], abs=1e-4)


# ---------------------------------------------------------------------------
# M. regional worst-point calculation matches direct point evaluation
# ---------------------------------------------------------------------------


def test_M_worst_point_matches_direct_evaluation(elements0):
    t_span = (0.0, 3 * 86400.0)
    dt = 300.0
    results = regional_coverage(elements0, t_span, dt, 60.0, 75.0, 5.0, 15.0, MIN_ELEVATION_DEG)
    summ = summarize_regional(results)

    n = int(round((t_span[1] - t_span[0]) / dt)) + 1
    t_s = np.linspace(t_span[0], t_span[1], n)
    traj = build_satellite_ecef_trajectory(elements0, t_s)
    direct = point_coverage(traj, t_s, summ.worst_gap_point[0], summ.worst_gap_point[1], MIN_ELEVATION_DEG)
    assert direct.max_gap_s == pytest.approx(summ.worst_max_gap_s, rel=1e-9)


# ---------------------------------------------------------------------------
# N. time-step convergence
# ---------------------------------------------------------------------------


def test_N_timestep_convergence(elements0):
    t_span = (0.0, 3 * 86400.0)
    n_dt = int(round((t_span[1] - t_span[0]) / 15.0)) + 1
    t_s_fine = np.linspace(t_span[0], t_span[1], n_dt)
    traj_fine_src = build_satellite_ecef_trajectory(elements0, t_s_fine)

    fractions = []
    for dt in [120.0, 60.0, 15.0]:
        n = int(round((t_span[1] - t_span[0]) / dt)) + 1
        t_s = np.linspace(t_span[0], t_span[1], n)
        traj = build_satellite_ecef_trajectory(elements0, t_s)
        res = point_coverage(traj, t_s, SITE_LAT_DEG, SITE_LON_DEG, MIN_ELEVATION_DEG)
        fractions.append(res.access_fraction)

    d1 = abs(fractions[0] - fractions[2])
    d2 = abs(fractions[1] - fractions[2])
    assert d2 <= d1 + 1e-9  # finer step -> closer to the finest reference


# ---------------------------------------------------------------------------
# O. grid refinement behavior
# ---------------------------------------------------------------------------


def test_O_grid_refinement_stability(elements0):
    t_span = (0.0, 3 * 86400.0)
    dt = 300.0
    coarse = summarize_regional(regional_coverage(elements0, t_span, dt, 60.0, 75.0, 5.0, 10.0, MIN_ELEVATION_DEG))
    fine = summarize_regional(regional_coverage(elements0, t_span, dt, 60.0, 75.0, 2.5, 5.0, MIN_ELEVATION_DEG))
    # aggregate means should be reasonably close between grid resolutions
    assert abs(coarse.point_weighted_mean_access_fraction - fine.point_weighted_mean_access_fraction) < 0.05


def test_O_build_grid_shape():
    lats, lons = build_grid(60.0, 75.0, 2.5, 2.5)
    assert lats[0] == 60.0 and lats[-1] == 75.0
    assert lons[0] == 0.0 and lons[-1] < 360.0
    assert len(lats) == 7
    assert len(lons) == 144


# ---------------------------------------------------------------------------
# P. independent elevation geometry cross-check
# ---------------------------------------------------------------------------


def test_P_independent_elevation_crosscheck(elements0):
    t_span = (0.0, 2 * 86400.0)
    dt = 600.0
    n = int(round((t_span[1] - t_span[0]) / dt)) + 1
    t_s = np.linspace(t_span[0], t_span[1], n)
    traj = build_satellite_ecef_trajectory(elements0, t_s)

    for lat, lon in [(65.0, 40.0), (60.0, 270.0), (75.0, 180.0)]:  # site, worst-ish, band-edge
        max_err = 0.0
        for k in range(0, len(t_s), 20):
            _, _, el1 = range_az_el(traj[k], lat, lon)
            _, el2 = elevation_independent_check(traj[k], lat, lon)
            max_err = max(max_err, abs(el1 - el2))
        assert max_err < 1e-8


# ---------------------------------------------------------------------------
# Q. M1-M4 regression tests remain green (spot check via production code)
# ---------------------------------------------------------------------------


def test_Q_regression_m1_m4_still_holds():
    rates = compute_secular_rates(BASE["a_km"], BASE["e"], BASE["i_deg"])
    assert rates.raan_dot_deg_day == pytest.approx(-0.145135, abs=1e-5)
    assert abs(rates.argp_dot_deg_day) < 1e-12


# ---------------------------------------------------------------------------
# R. dimensional / hour / day conversions
# ---------------------------------------------------------------------------


def test_R_unit_conversions():
    hours_per_day = 24.0
    seconds_per_hour = 3600.0
    assert SIDEREAL_DAY_S / seconds_per_hour == pytest.approx(23.9344695, abs=1e-4)
    # gap in seconds -> hours round trip
    gap_s = 8551.051
    gap_h = gap_s / seconds_per_hour
    assert gap_h * seconds_per_hour == pytest.approx(gap_s, rel=1e-12)
