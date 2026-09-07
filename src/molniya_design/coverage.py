"""Regional high-latitude coverage / revisit analysis (M5).

Adds **no new orbital dynamics**: reuses M2's Kepler-solving machinery,
M3's frame/access geometry, and M4's secular-J2 propagation unchanged.
The only new production logic here is (a) vectorized primitives needed
for regional-grid performance (the satellite trajectory is independent
of ground point, so it is computed **once** and reused for every grid
point — see :func:`build_satellite_ecef_trajectory`, which in turn uses
`elements.coe_to_r_array` and `frames.eci_array_to_ecef_array`), and (b)
grid construction / regional aggregation.

Scope: geometric access only (spherical Earth, first-order secular J2,
no refraction, no link budget) — see DESIGN.md limitations. This module
answers "how much line-of-sight geometric access does one spacecraft
give," not "how much communications capacity."
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .access import (
    access_metrics_boundary_aware,
    find_access_intervals,
    range_az_el_series,
)
from .constants import J2_EARTH, MU_EARTH, R_EARTH, SIDEREAL_DAY_S
from .elements import coe_to_r_array
from .frames import eci_array_to_ecef_array
from .j2 import SecularElements, compute_secular_rates
from .propagation import eccentric_to_true_anomaly, mean_to_eccentric_anomaly_array


def build_satellite_ecef_trajectory(
    elements0: SecularElements,
    t_s: np.ndarray,
    J2: float = J2_EARTH,
    Re: float = R_EARTH,
    mu: float = MU_EARTH,
    theta_g0_rad: float = 0.0,
) -> np.ndarray:
    """Vectorized satellite ECEF trajectory under the M4 secular-J2
    model, shape (N, 3) km. Computed once per (elements0, t_s) and reused
    across every ground point in a regional grid, since satellite
    dynamics do not depend on any ground point."""
    t_s = np.asarray(t_s, dtype=float)
    rates = compute_secular_rates(elements0.a_km, elements0.e, elements0.i_deg, J2, Re, mu)

    raan_deg = elements0.raan_deg + np.degrees(rates.raan_dot_rad_s) * t_s
    argp_deg = elements0.argp_deg + np.degrees(rates.argp_dot_rad_s) * t_s
    m_rad = np.radians(elements0.m_deg) + rates.m_dot_rad_s * t_s

    E_rad = mean_to_eccentric_anomaly_array(m_rad, elements0.e)
    nu_rad = eccentric_to_true_anomaly(E_rad, elements0.e)
    nu_deg = np.degrees(nu_rad)

    r_eci = coe_to_r_array(elements0.a_km, elements0.e, elements0.i_deg, raan_deg, argp_deg, nu_deg, mu=mu)
    r_ecef = eci_array_to_ecef_array(r_eci, t_s, theta_g0_rad)
    return r_ecef


# ---------------------------------------------------------------------------
# Regional grid
# ---------------------------------------------------------------------------


def build_grid(lat_min_deg: float, lat_max_deg: float, lat_step_deg: float, lon_step_deg: float):
    """Regular lat/lon grid: latitudes in [lat_min, lat_max] inclusive at
    lat_step spacing, longitudes covering the full [0, 360) range at
    lon_step spacing (no duplicate 0/360 point)."""
    n_lat = int(round((lat_max_deg - lat_min_deg) / lat_step_deg)) + 1
    lats = np.linspace(lat_min_deg, lat_max_deg, n_lat)
    n_lon = int(round(360.0 / lon_step_deg))
    lons = np.linspace(0.0, 360.0, n_lon, endpoint=False)
    return lats, lons


# ---------------------------------------------------------------------------
# Per-point coverage
# ---------------------------------------------------------------------------


@dataclass
class PointCoverageResult:
    lat_deg: float
    lon_deg: float
    total_access_time_s: float
    access_fraction: float
    num_passes: int
    num_complete_passes: int
    has_leading_partial: bool
    has_trailing_partial: bool
    longest_pass_s: float
    max_gap_s: float
    peak_elevation_deg: float


def point_coverage(
    r_ecef_traj: np.ndarray,
    t_s: np.ndarray,
    lat_deg: float,
    lon_deg: float,
    elevation_min_deg: float,
    periodic: bool = False,
    period_s: float | None = None,
) -> PointCoverageResult:
    """Access metrics for one ground point, given a precomputed (N,3)
    satellite ECEF trajectory (reuses M3's `range_az_el_series` and
    `find_access_intervals` unchanged)."""
    r_ecef_traj_3n = r_ecef_traj.T  # (3, N), matches access.range_az_el_series's convention
    _, _, el = range_az_el_series(r_ecef_traj_3n, lat_deg, lon_deg)
    intervals = find_access_intervals(t_s, el, threshold_deg=elevation_min_deg)
    m = access_metrics_boundary_aware(intervals, t_s[0], t_s[-1], periodic=periodic, period_s=period_s)
    peak_el = max((iv.peak_elevation_deg for iv in intervals), default=float(np.max(el)))

    return PointCoverageResult(
        lat_deg=lat_deg,
        lon_deg=lon_deg,
        total_access_time_s=m["total_access_time_s"],
        access_fraction=m["access_fraction"],
        num_passes=m["num_passes"],
        num_complete_passes=m["num_complete_passes"],
        has_leading_partial=m["has_leading_partial"],
        has_trailing_partial=m["has_trailing_partial"],
        longest_pass_s=m["longest_pass_s"],
        max_gap_s=m["max_gap_s"],
        peak_elevation_deg=peak_el,
    )


def regional_coverage(
    elements0: SecularElements,
    t_span: tuple[float, float],
    dt_s: float,
    lat_min_deg: float,
    lat_max_deg: float,
    lat_step_deg: float,
    lon_step_deg: float,
    elevation_min_deg: float,
    J2: float = J2_EARTH,
    Re: float = R_EARTH,
    mu: float = MU_EARTH,
    theta_g0_rad: float = 0.0,
    periodic: bool = False,
    period_s: float | None = None,
) -> list[PointCoverageResult]:
    """Regional coverage over a lat/lon grid. The satellite trajectory is
    built once (§ build_satellite_ecef_trajectory) and reused for every
    grid point — the dominant cost is O(n_grid * n_time) vectorized ENU
    projections (`range_az_el_series`), not satellite dynamics."""
    n = int(round((t_span[1] - t_span[0]) / dt_s)) + 1
    t_s = np.linspace(t_span[0], t_span[1], n)
    r_ecef_traj = build_satellite_ecef_trajectory(elements0, t_s, J2, Re, mu, theta_g0_rad)

    lats, lons = build_grid(lat_min_deg, lat_max_deg, lat_step_deg, lon_step_deg)
    results = []
    for lat in lats:
        for lon in lons:
            results.append(point_coverage(r_ecef_traj, t_s, lat, lon, elevation_min_deg, periodic, period_s))
    return results


# ---------------------------------------------------------------------------
# Regional aggregation
# ---------------------------------------------------------------------------


@dataclass
class RegionalSummary:
    point_weighted_mean_access_fraction: float
    area_weighted_mean_access_fraction: float
    min_access_fraction: float
    max_access_fraction: float
    worst_max_gap_s: float
    median_max_gap_s: float
    worst_gap_point: tuple[float, float]
    best_access_point: tuple[float, float]


def summarize_regional(results: list[PointCoverageResult]) -> RegionalSummary:
    """Point-weighted mean = simple average over grid points (each point
    counted equally regardless of the shrinking physical area it
    represents at high latitude). Area-weighted mean uses weight ∝
    cos(latitude), the standard spherical-Earth area element weighting,
    and is the more physically meaningful "how much of the actual band
    area is covered" statistic. Both are reported and explicitly
    labeled, per M5 §2."""
    af = np.array([r.access_fraction for r in results])
    gaps = np.array([r.max_gap_s for r in results])
    lats = np.array([r.lat_deg for r in results])

    weights_area = np.cos(np.radians(lats))
    point_weighted_mean = float(af.mean())
    area_weighted_mean = float(np.sum(af * weights_area) / np.sum(weights_area))

    worst_gap_idx = int(np.argmax(gaps))
    best_access_idx = int(np.argmax(af))

    return RegionalSummary(
        point_weighted_mean_access_fraction=point_weighted_mean,
        area_weighted_mean_access_fraction=area_weighted_mean,
        min_access_fraction=float(af.min()),
        max_access_fraction=float(af.max()),
        worst_max_gap_s=float(gaps[worst_gap_idx]),
        median_max_gap_s=float(np.median(gaps)),
        worst_gap_point=(results[worst_gap_idx].lat_deg, results[worst_gap_idx].lon_deg),
        best_access_point=(results[best_access_idx].lat_deg, results[best_access_idx].lon_deg),
    )
