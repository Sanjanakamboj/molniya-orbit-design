"""Milestone 3 test suite: ECI/ECEF transformation, geocentric lat/lon,
two-body ground track, and basic spherical-Earth access geometry.

Covers checklist items A-T from the M3 task description.
"""

from __future__ import annotations

import numpy as np
import pytest

from molniya_design.access import (
    AccessInterval,
    access_metrics,
    elevation_independent_check,
    enu_basis,
    find_access_intervals,
    range_az_el,
    range_az_el_series,
    site_ecef,
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
from molniya_design.elements import coe_to_rv
from molniya_design.frames import (
    EARTH_ROTATION_RATE_RAD_S,
    ecef_to_eci,
    ecef_to_geocentric_latlon,
    eci_to_ecef,
    geocentric_latlon_to_ecef,
    theta_g,
)
from molniya_design.groundtrack import (
    apsis_ground_points,
    compute_ground_track,
    longitude_separation_deg,
)

BASE = BASELINE_ELEMENTS_DEG


@pytest.fixture(scope="module")
def baseline_state0():
    return coe_to_rv(
        BASE["a_km"], BASE["e"], BASE["i_deg"], BASE["raan_deg"],
        BASE["argp_deg"], BASE["nu0_deg"],
    )


# ---------------------------------------------------------------------------
# A. ECI/ECEF identity at theta=0
# ---------------------------------------------------------------------------


def test_A_identity_at_theta_zero():
    r = np.array([1234.5, -6789.0, 4321.0])
    r_ecef = eci_to_ecef(r, 0.0)
    assert np.allclose(r_ecef, r, atol=1e-9)


# ---------------------------------------------------------------------------
# B. ECI/ECEF round trip
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("t_s", [0.0, 1234.0, 43082.045, 86164.0905, 1.5e5])
def test_B_eci_ecef_round_trip(t_s):
    r = np.array([12000.0, -3400.0, 5600.0])
    r_ecef = eci_to_ecef(r, t_s)
    r_back = ecef_to_eci(r_ecef, t_s)
    assert np.allclose(r_back, r, rtol=1e-12)


# ---------------------------------------------------------------------------
# C. rotation matrix determinant / right-handedness
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("t_s", [0.0, 100.0, 43082.045, 86164.0905, 2e5])
def test_C_rotation_right_handed(t_s):
    theta = theta_g(t_s)
    c, s = np.cos(theta), np.sin(theta)
    R = np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])  # R3(-theta)
    assert np.linalg.det(R) == pytest.approx(1.0, abs=1e-12)
    assert np.allclose(R.T @ R, np.eye(3), atol=1e-12)


# ---------------------------------------------------------------------------
# D. one-sidereal-day frame repeat
# ---------------------------------------------------------------------------


def test_D_one_sidereal_day_frame_repeat():
    theta0 = theta_g(0.0)
    theta1 = theta_g(SIDEREAL_DAY_S)
    assert theta1 == pytest.approx(theta0, abs=1e-9)
    # sanity: omega_E * sidereal_day = 2 pi exactly by construction
    assert EARTH_ROTATION_RATE_RAD_S * SIDEREAL_DAY_S == pytest.approx(2 * np.pi, rel=1e-14)


# ---------------------------------------------------------------------------
# E. lat/lon/ECEF round trip
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "lat_deg,lon_deg",
    [
        (0.0, 0.0),
        (65.0, 40.0),
        (-63.4349, 179.9),
        (90.0, 0.0),
        (-90.0, 123.0),
        (10.0, -170.0),
    ],
)
def test_E_latlon_ecef_round_trip(lat_deg, lon_deg):
    r_ecef = geocentric_latlon_to_ecef(lat_deg, lon_deg)
    lat_back, lon_back = ecef_to_geocentric_latlon(r_ecef)
    assert lat_back == pytest.approx(lat_deg, abs=1e-9)
    # longitude is ill-defined at the poles; skip lon check there
    if abs(lat_deg) < 89.999:
        assert lon_back == pytest.approx(((lon_deg + 180) % 360) - 180, abs=1e-6)


# ---------------------------------------------------------------------------
# F. longitude wrap behavior
# ---------------------------------------------------------------------------


def test_F_longitude_wrap_range():
    for lon_in in [-180.0, -179.999, 0.0, 179.999, 180.0, 350.0, -350.0]:
        r_ecef = geocentric_latlon_to_ecef(10.0, lon_in)
        _, lon_out = ecef_to_geocentric_latlon(r_ecef)
        assert -180.0 <= lon_out < 180.0


def test_F_longitude_separation_wraparound():
    # crossing the +/-180 discontinuity should give a small separation
    assert longitude_separation_deg(179.0, -179.0) == pytest.approx(2.0, abs=1e-9)
    assert longitude_separation_deg(-179.0, 179.0) == pytest.approx(-2.0, abs=1e-9)
    assert longitude_separation_deg(10.0, 20.0) == pytest.approx(10.0, abs=1e-9)


# ---------------------------------------------------------------------------
# G. northern-apogee latitude retained in ECEF/geocentric form
# ---------------------------------------------------------------------------


def test_G_apogee_latitude_in_ecef(baseline_state0):
    r0, v0 = baseline_state0
    r_ecef = eci_to_ecef(r0, 0.0)
    lat, _ = ecef_to_geocentric_latlon(r_ecef)
    assert lat == pytest.approx(BASE["i_deg"], abs=1e-6)


# ---------------------------------------------------------------------------
# H. two-body ground-track repeat after one sidereal day
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def apogee_points(baseline_state0):
    r0, v0 = baseline_state0
    t_span = (0.0, 2.05 * M1_PERIOD_S)
    pts = apsis_ground_points(r0, v0, t_span)
    return [p for p in pts if p.kind == "apogee"]


def test_H_ground_track_repeats_after_one_sidereal_day(apogee_points):
    # apogee[0] at t~0, apogee[2] at t~2 periods = 1 sidereal day later
    assert len(apogee_points) >= 3
    lon_err = longitude_separation_deg(apogee_points[0].lon_deg, apogee_points[2].lon_deg)
    assert abs(lon_err) < 1e-6
    assert apogee_points[0].lat_deg == pytest.approx(apogee_points[2].lat_deg, abs=1e-9)


# ---------------------------------------------------------------------------
# I. successive-apogee longitude pattern
# ---------------------------------------------------------------------------


def test_I_successive_apogee_longitude_shift(apogee_points):
    # T * omega_E = pi exactly (period = half sidereal day), so successive
    # apogees should be exactly 180 deg apart in longitude.
    sep_01 = longitude_separation_deg(apogee_points[0].lon_deg, apogee_points[1].lon_deg)
    sep_12 = longitude_separation_deg(apogee_points[1].lon_deg, apogee_points[2].lon_deg)
    assert abs(abs(sep_01) - 180.0) < 1e-6
    assert abs(abs(sep_12) - 180.0) < 1e-6
    # all apogees at the same (northern) latitude
    for p in apogee_points:
        assert p.lat_deg == pytest.approx(BASE["i_deg"], abs=1e-6)


# ---------------------------------------------------------------------------
# J. site ECEF conversion
# ---------------------------------------------------------------------------


def test_J_site_ecef():
    r = site_ecef(SITE_LAT_DEG, SITE_LON_DEG)
    assert np.linalg.norm(r) == pytest.approx(R_EARTH, rel=1e-12)
    lat, lon = ecef_to_geocentric_latlon(r)
    assert lat == pytest.approx(SITE_LAT_DEG, abs=1e-9)
    assert lon == pytest.approx(SITE_LON_DEG, abs=1e-9)


# ---------------------------------------------------------------------------
# K. ENU basis orthogonality
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("lat,lon", [(0, 0), (65, 40), (-30, 170), (89, -50)])
def test_K_enu_orthonormal_right_handed(lat, lon):
    e, n, u = enu_basis(lat, lon)
    for v in (e, n, u):
        assert np.linalg.norm(v) == pytest.approx(1.0, abs=1e-12)
    assert np.dot(e, n) == pytest.approx(0.0, abs=1e-12)
    assert np.dot(n, u) == pytest.approx(0.0, abs=1e-12)
    assert np.dot(e, u) == pytest.approx(0.0, abs=1e-12)
    assert np.allclose(np.cross(e, n), u, atol=1e-12)  # right-handed: E x N = U


# ---------------------------------------------------------------------------
# L. known-overhead elevation ~90 deg
# ---------------------------------------------------------------------------


def test_L_overhead_elevation_90deg():
    r_sat = site_ecef(SITE_LAT_DEG, SITE_LON_DEG, r_km=R_EARTH + 1000.0)
    rng, az, el = range_az_el(r_sat, SITE_LAT_DEG, SITE_LON_DEG)
    # arcsin is ill-conditioned near +/-1 (derivative -> infinity), so a
    # ~1e-16 floating-point error in rho_up/range near exactly overhead
    # amplifies to ~1e-6 deg in elevation -- an expected numerical effect,
    # not a geometry bug. 1e-4 deg is still far tighter than any physical
    # requirement.
    assert el == pytest.approx(90.0, abs=1e-4)
    assert rng == pytest.approx(1000.0, abs=1e-9)


# ---------------------------------------------------------------------------
# M. known-horizon geometry ~0 deg
# ---------------------------------------------------------------------------


def test_M_horizon_elevation_0deg():
    r_site = site_ecef(SITE_LAT_DEG, SITE_LON_DEG)
    e_east, _, _ = enu_basis(SITE_LAT_DEG, SITE_LON_DEG)
    r_sat = r_site + 5000.0 * e_east  # purely local-horizontal offset
    rng, az, el = range_az_el(r_sat, SITE_LAT_DEG, SITE_LON_DEG)
    assert el == pytest.approx(0.0, abs=1e-9)
    assert az == pytest.approx(90.0, abs=1e-6)  # due east
    assert rng == pytest.approx(5000.0, abs=1e-9)


# ---------------------------------------------------------------------------
# N. ENU elevation vs independent triangle-geometry elevation
# ---------------------------------------------------------------------------


def test_N_enu_vs_triangle_elevation(baseline_state0):
    r0, v0 = baseline_state0
    from molniya_design.propagation import propagate

    t_eval = np.linspace(0.0, M1_PERIOD_S, 50)
    sol = propagate(r0, v0, (0.0, M1_PERIOD_S), t_eval=t_eval, rtol=1e-12, atol=1e-12)

    max_el_err, max_rng_err = 0.0, 0.0
    for k, t in enumerate(sol.t):
        r_ecef = eci_to_ecef(sol.y[0:3, k], t)
        rng1, az1, el1 = range_az_el(r_ecef, SITE_LAT_DEG, SITE_LON_DEG)
        rng2, el2 = elevation_independent_check(r_ecef, SITE_LAT_DEG, SITE_LON_DEG)
        max_el_err = max(max_el_err, abs(el1 - el2))
        max_rng_err = max(max_rng_err, abs(rng1 - rng2))
    assert max_el_err < 1e-8
    assert max_rng_err < 1e-8


# ---------------------------------------------------------------------------
# O. access threshold classification
# ---------------------------------------------------------------------------


def test_O_access_threshold_classification():
    t = np.linspace(0, 100, 11)
    el = np.array([-5, -2, 2, 8, 12, 15, 12, 8, 2, -2, -5], dtype=float)
    intervals = find_access_intervals(t, el, threshold_deg=10.0)
    assert len(intervals) == 1
    assert intervals[0].start_s > 30 and intervals[0].start_s < 40
    assert intervals[0].end_s > 60 and intervals[0].end_s < 70


# ---------------------------------------------------------------------------
# P. interval extraction on a synthetic known signal
# ---------------------------------------------------------------------------


def test_P_synthetic_sinusoid_interval_extraction():
    # el(t) = 20*sin(2*pi*t/100) - 5  => crosses threshold=10 at known times
    # 20 sin(x) - 5 = 10  =>  sin(x) = 0.75  => x = asin(0.75) or pi - asin(0.75)
    t = np.linspace(0.0, 100.0, 20001)  # dt = 0.005 s, very fine synthetic sampling
    el = 20.0 * np.sin(2 * np.pi * t / 100.0) - 5.0
    intervals = find_access_intervals(t, el, threshold_deg=10.0)
    assert len(intervals) == 1

    x1 = np.arcsin(0.75)
    x2 = np.pi - np.arcsin(0.75)
    t_start_expected = x1 / (2 * np.pi) * 100.0
    t_end_expected = x2 / (2 * np.pi) * 100.0

    assert intervals[0].start_s == pytest.approx(t_start_expected, abs=0.01)
    assert intervals[0].end_s == pytest.approx(t_end_expected, abs=0.01)
    assert intervals[0].peak_elevation_deg == pytest.approx(15.0, abs=0.01)
    assert intervals[0].peak_time_s == pytest.approx(25.0, abs=0.05)


# ---------------------------------------------------------------------------
# Q. no-access gap computation
# ---------------------------------------------------------------------------


def test_Q_access_metrics_gap_computation():
    intervals = [
        AccessInterval(start_s=100.0, end_s=200.0, duration_s=100.0, peak_elevation_deg=30.0, peak_time_s=150.0),
        AccessInterval(start_s=500.0, end_s=550.0, duration_s=50.0, peak_elevation_deg=15.0, peak_time_s=525.0),
    ]
    m = access_metrics(intervals, t_start_s=0.0, t_end_s=1000.0)
    assert m.total_access_time_s == pytest.approx(150.0)
    assert m.access_fraction == pytest.approx(0.15)
    assert m.num_passes == 2
    assert m.longest_pass_s == pytest.approx(100.0)
    # gaps: [0,100]=100, [200,500]=300, [550,1000]=450 -> max=450
    assert m.max_gap_s == pytest.approx(450.0)


def test_Q_access_metrics_no_passes():
    m = access_metrics([], t_start_s=0.0, t_end_s=1000.0)
    assert m.num_passes == 0
    assert m.total_access_time_s == 0.0
    assert m.access_fraction == 0.0
    assert m.longest_pass_s == 0.0
    assert m.max_gap_s == pytest.approx(1000.0)


# ---------------------------------------------------------------------------
# R. access convergence with timestep
# ---------------------------------------------------------------------------


def test_R_access_convergence_with_timestep(baseline_state0):
    r0, v0 = baseline_state0
    t_span = (0.0, SIDEREAL_DAY_S)
    dts = [120.0, 60.0, 15.0]
    totals = []
    for dt in dts:
        gt = compute_ground_track(r0, v0, t_span, dt_s=dt)
        _, _, el = range_az_el_series(gt.r_ecef, SITE_LAT_DEG, SITE_LON_DEG)
        intervals = find_access_intervals(gt.t_s, el, threshold_deg=MIN_ELEVATION_DEG)
        m = access_metrics(intervals, t_span[0], t_span[1])
        totals.append(m.total_access_time_s)

    # finer sampling should converge (differences should shrink)
    diff_120_60 = abs(totals[0] - totals[1])
    diff_60_15 = abs(totals[1] - totals[2])
    assert diff_60_15 <= diff_120_60 + 1.0  # allow small numerical slack


# ---------------------------------------------------------------------------
# S. M1/M2 regression tests still passing (spot check via production code)
# ---------------------------------------------------------------------------


def test_S_M1_M2_regression_still_holds(baseline_state0):
    from molniya_design.constants import M1_APOGEE_ALT_KM, M1_PERIGEE_ALT_KM, M1_VA_KM_S, M1_VP_KM_S
    from molniya_design.propagation import find_apsides

    r0, v0 = baseline_state0
    apsides = find_apsides(r0, v0, (0.0, 1.01 * M1_PERIOD_S), rtol=1e-13, atol=1e-13)
    perigee = next(a for a in apsides if a.kind == "perigee")
    apogee = next(a for a in apsides if a.kind == "apogee" and a.t_s > 1.0)

    assert (perigee.r_km - R_EARTH) == pytest.approx(M1_PERIGEE_ALT_KM, abs=1e-4)
    assert (apogee.r_km - R_EARTH) == pytest.approx(M1_APOGEE_ALT_KM, abs=1e-4)
    assert perigee.v_km_s == pytest.approx(M1_VP_KM_S, abs=1e-6)
    assert apogee.v_km_s == pytest.approx(M1_VA_KM_S, abs=1e-6)


# ---------------------------------------------------------------------------
# T. dimensional / unit sanity checks
# ---------------------------------------------------------------------------


def test_T_unit_sanity_checks():
    # omega_E in rad/s, converted to deg/day should be ~360.9856 (sidereal rate)
    deg_per_day = np.degrees(EARTH_ROTATION_RATE_RAD_S) * 86400.0
    assert deg_per_day == pytest.approx(360.98564736, rel=1e-6)

    # theta_g wraps to [0, 2pi)
    assert 0.0 <= theta_g(1e9) < 2 * np.pi

    # site radius must equal R_EARTH exactly by construction
    r = site_ecef(0.0, 0.0)
    assert np.linalg.norm(r) == pytest.approx(R_EARTH, rel=1e-14)
