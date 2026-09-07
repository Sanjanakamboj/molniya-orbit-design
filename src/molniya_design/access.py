"""Basic spherical-Earth topocentric access geometry.

Scope (M3): geometry only — range, azimuth, elevation, and threshold-based
access intervals from a spherical-Earth line-of-sight test. This is
**not** a communications-availability or RF link-budget result (no
antenna gain, no refraction, no link margin) — see DESIGN.md limitations.

Site and ENU/SEZ convention
----------------------------
The ground site is a fixed point on a sphere of radius ``R_EARTH`` at
geocentric (lat, lon). Its ECEF position and the local East-North-Up
(ENU) basis (all unit vectors, expressed in ECEF components) are:

    r_site = R_E * [cos(lat)cos(lon), cos(lat)sin(lon), sin(lat)]

    e_east  = [-sin(lon),            cos(lon),           0        ]
    e_north = [-sin(lat)cos(lon), -sin(lat)sin(lon), cos(lat)]
    e_up    = [ cos(lat)cos(lon),  cos(lat)sin(lon), sin(lat)]

which is a right-handed, orthonormal triad (tested).

For a satellite ECEF position ``r_sat``, the topocentric line-of-sight
vector is ``rho = r_sat - r_site``. Projecting onto the ENU basis:

    rho_e = rho . e_east,   rho_n = rho . e_north,   rho_u = rho . e_up

    range     = |rho|
    elevation = asin(rho_u / range)                      [-90, +90] deg
    azimuth   = atan2(rho_e, rho_n) mod 360                [0, 360) deg
                (0 = north, 90 = east — clockwise from north)

No atmospheric refraction is modeled (elevation is the pure geometric
line-of-sight angle above the local horizontal plane).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import R_EARTH


def site_ecef(lat_deg: float, lon_deg: float, r_km: float = R_EARTH) -> np.ndarray:
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    return r_km * np.array(
        [np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)]
    )


def enu_basis(lat_deg: float, lon_deg: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (e_east, e_north, e_up) unit vectors in ECEF components."""
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    e_east = np.array([-np.sin(lon), np.cos(lon), 0.0])
    e_north = np.array([-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)])
    e_up = np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
    return e_east, e_north, e_up


def range_az_el(
    r_sat_ecef: np.ndarray, lat_deg: float, lon_deg: float, r_site_km: float = R_EARTH
) -> tuple[float, float, float]:
    """ENU-based (range_km, azimuth_deg, elevation_deg) of a satellite ECEF
    position as seen from a spherical-Earth ground site."""
    r_site = site_ecef(lat_deg, lon_deg, r_site_km)
    e_east, e_north, e_up = enu_basis(lat_deg, lon_deg)

    rho = np.asarray(r_sat_ecef, dtype=float) - r_site
    rng = np.linalg.norm(rho)

    rho_e = np.dot(rho, e_east)
    rho_n = np.dot(rho, e_north)
    rho_u = np.dot(rho, e_up)

    el = np.degrees(np.arcsin(rho_u / rng))
    az = np.degrees(np.arctan2(rho_e, rho_n)) % 360.0
    return rng, az, el


def range_az_el_series(
    r_sat_ecef_traj: np.ndarray, lat_deg: float, lon_deg: float, r_site_km: float = R_EARTH
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Vectorized :func:`range_az_el` over a (3, N) ECEF trajectory.
    Returns (range_km, azimuth_deg, elevation_deg), each shape (N,)."""
    r_sat_ecef_traj = np.asarray(r_sat_ecef_traj, dtype=float)
    r_site = site_ecef(lat_deg, lon_deg, r_site_km)
    e_east, e_north, e_up = enu_basis(lat_deg, lon_deg)

    rho = r_sat_ecef_traj - r_site[:, None]
    rng = np.linalg.norm(rho, axis=0)

    rho_e = e_east @ rho
    rho_n = e_north @ rho
    rho_u = e_up @ rho

    el = np.degrees(np.arcsin(rho_u / rng))
    az = np.degrees(np.arctan2(rho_e, rho_n)) % 360.0
    return rng, az, el


def elevation_independent_check(
    r_sat_ecef: np.ndarray, lat_deg: float, lon_deg: float, r_site_km: float = R_EARTH
) -> tuple[float, float]:
    """Independent (range_km, elevation_deg) via the Earth-center/site/
    satellite triangle, using the geocentric angle psi between the site
    and satellite direction vectors:

        cos(psi) = r_site_hat . r_sat_hat

    Standard closed-form triangle-geometry relations (law of cosines for
    range, and the direct psi -> elevation identity, e.g. Vallado eq.
    5-x): with r_sat = |r_sat_ecef| and R_E = r_site_km,

        rho = sqrt(r_sat^2 + R_E^2 - 2 R_E r_sat cos(psi))
        tan(el) = (cos(psi) - R_E/r_sat) / sin(psi)

    This does **not** use the ENU basis or dot-products against e_east/
    e_north/e_up at all — a genuinely separate derivation path from
    :func:`range_az_el`, used purely as a cross-check.
    """
    r_site = site_ecef(lat_deg, lon_deg, r_site_km)
    r_sat = np.asarray(r_sat_ecef, dtype=float)

    r_sat_mag = np.linalg.norm(r_sat)
    r_site_mag = np.linalg.norm(r_site)

    cos_psi = np.dot(r_site, r_sat) / (r_site_mag * r_sat_mag)
    cos_psi = np.clip(cos_psi, -1.0, 1.0)
    psi = np.arccos(cos_psi)
    sin_psi = np.sin(psi)

    rho = np.sqrt(r_sat_mag**2 + r_site_mag**2 - 2.0 * r_site_mag * r_sat_mag * cos_psi)

    if abs(sin_psi) < 1e-12:
        # satellite directly overhead (psi ~ 0): elevation -> 90 deg
        el = 90.0 if cos_psi > 0 else -90.0
    else:
        el = np.degrees(np.arctan2(cos_psi - r_site_mag / r_sat_mag, sin_psi))
    return rho, el


# ---------------------------------------------------------------------------
# Access interval extraction
# ---------------------------------------------------------------------------


@dataclass
class AccessInterval:
    start_s: float
    end_s: float
    duration_s: float
    peak_elevation_deg: float
    peak_time_s: float
    is_boundary_start: bool = False
    """True if this interval's access run was already in progress at the
    very first sample of the series (elevation >= threshold at t[0]) —
    i.e. we do not know from this window alone when it actually began.
    See M5 §3 / `summarize_boundary_passes` for how this is used."""
    is_boundary_end: bool = False
    """True if this interval's access run is still in progress at the
    very last sample of the series — we do not know when it actually
    ends. See M5 §3 / `summarize_boundary_passes`."""


def find_access_intervals(
    t_s: np.ndarray, elevation_deg: np.ndarray, threshold_deg: float = 10.0
) -> list[AccessInterval]:
    """Extract access intervals (elevation >= threshold) from a sampled
    elevation time series.

    Threshold crossings are refined by **linear interpolation** between
    the bracketing samples (not simply the nearest sample), and the peak
    elevation within each interval is refined with a local quadratic
    (parabolic) fit around the sampled maximum when interior points are
    available. Requires small enough sample spacing that elevation is
    well-approximated as monotonic between crossing samples and locally
    parabolic near the peak — validated by the M3 convergence study.
    """
    t_s = np.asarray(t_s, dtype=float)
    el = np.asarray(elevation_deg, dtype=float)
    above = el >= threshold_deg

    intervals: list[AccessInterval] = []
    n = len(t_s)
    i = 0
    while i < n:
        if not above[i]:
            i += 1
            continue
        # start of an access run
        start_idx = i
        # refine start crossing time (if not the very first sample)
        if start_idx > 0 and not above[start_idx - 1]:
            t0, t1 = t_s[start_idx - 1], t_s[start_idx]
            e0, e1 = el[start_idx - 1], el[start_idx]
            frac = (threshold_deg - e0) / (e1 - e0)
            start_t = t0 + frac * (t1 - t0)
        else:
            start_t = t_s[start_idx]

        j = start_idx
        while j < n and above[j]:
            j += 1
        end_idx = j - 1
        if end_idx + 1 < n:
            t0, t1 = t_s[end_idx], t_s[end_idx + 1]
            e0, e1 = el[end_idx], el[end_idx + 1]
            frac = (threshold_deg - e0) / (e1 - e0)
            end_t = t0 + frac * (t1 - t0)
        else:
            end_t = t_s[end_idx]

        # peak within [start_idx, end_idx], refined with a parabolic fit
        local_idx = np.arange(start_idx, end_idx + 1)
        peak_local = local_idx[np.argmax(el[local_idx])]
        if 0 < peak_local < n - 1:
            y0, y1, y2 = el[peak_local - 1], el[peak_local], el[peak_local + 1]
            denom = y0 - 2 * y1 + y2
            if abs(denom) > 1e-14:
                delta = 0.5 * (y0 - y2) / denom
                delta = np.clip(delta, -1.0, 1.0)
            else:
                delta = 0.0
            dt = t_s[peak_local + 1] - t_s[peak_local]
            peak_t = t_s[peak_local] + delta * dt
            peak_el = y1 - 0.25 * (y0 - y2) * delta
        else:
            peak_t = t_s[peak_local]
            peak_el = el[peak_local]

        intervals.append(
            AccessInterval(
                start_s=start_t,
                end_s=end_t,
                duration_s=end_t - start_t,
                peak_elevation_deg=peak_el,
                peak_time_s=peak_t,
                is_boundary_start=(start_idx == 0),
                is_boundary_end=(end_idx == n - 1),
            )
        )
        i = j
    return intervals


@dataclass
class AccessMetrics:
    total_access_time_s: float
    access_fraction: float
    num_passes: int
    longest_pass_s: float
    max_gap_s: float


def access_metrics(
    intervals: list[AccessInterval], t_start_s: float, t_end_s: float
) -> AccessMetrics:
    total_duration = t_end_s - t_start_s
    total_access = sum(iv.duration_s for iv in intervals)
    num_passes = len(intervals)
    longest = max((iv.duration_s for iv in intervals), default=0.0)

    if num_passes == 0:
        max_gap = total_duration
    else:
        gaps = [intervals[0].start_s - t_start_s]
        for k in range(1, num_passes):
            gaps.append(intervals[k].start_s - intervals[k - 1].end_s)
        gaps.append(t_end_s - intervals[-1].end_s)
        max_gap = max(gaps)

    return AccessMetrics(
        total_access_time_s=total_access,
        access_fraction=total_access / total_duration,
        num_passes=num_passes,
        longest_pass_s=longest,
        max_gap_s=max_gap,
    )


# ---------------------------------------------------------------------------
# M5: window-boundary pass handling (fixes the M3 §5 "3 passes when there
# are physically 2" issue, generalized to both the periodic two-body case
# and the non-periodic multi-day J2 case)
# ---------------------------------------------------------------------------


def merge_cyclic_boundary_intervals(
    intervals: list[AccessInterval], t_start_s: float, t_end_s: float, period_s: float
) -> list[AccessInterval]:
    """For a **genuinely periodic/cyclic** analysis window (e.g. exactly
    one sidereal day of two-body dynamics, where the ground track is
    verified to repeat — M3 §M3.3), merge a trailing boundary interval
    with a leading boundary interval into one physical pass, since they
    are known to be the same continuous access run wrapping around the
    window edge.

    Do **not** call this for the M4/M5 J2 case — the ground track does
    not exactly repeat under J2 (M4 §M4.5), so a window-edge interval
    pair is not reliably the same physical pass. Use
    :func:`summarize_boundary_passes` for that case instead.
    """
    if len(intervals) < 2:
        return intervals
    first, last = intervals[0], intervals[-1]
    if first is last or not (first.is_boundary_start and last.is_boundary_end):
        return intervals

    merged_duration = first.duration_s + last.duration_s
    if first.peak_elevation_deg >= last.peak_elevation_deg:
        peak_el, peak_t = first.peak_elevation_deg, first.peak_time_s
    else:
        peak_el, peak_t = last.peak_elevation_deg, last.peak_time_s

    merged = AccessInterval(
        start_s=last.start_s,
        end_s=first.end_s + period_s,  # "unwrapped" bookkeeping time; duration_s is the meaningful field
        duration_s=merged_duration,
        peak_elevation_deg=peak_el,
        peak_time_s=peak_t,
        is_boundary_start=False,
        is_boundary_end=False,
    )
    return [merged] + intervals[1:-1]


def summarize_boundary_passes(intervals: list[AccessInterval]) -> dict:
    """For a **non-periodic** (e.g. multi-day J2) analysis window: report
    the physical pass count honestly, distinguishing complete interior
    passes from leading/trailing partial passes truncated by the window
    edge. Unlike the periodic case, a leading and trailing boundary
    interval are *not* assumed to be the same physical pass (the ground
    track has drifted, M4 §M4.5) — they are reported separately as
    "partial" rather than merged or silently counted as full passes."""
    n = len(intervals)
    if n == 0:
        return {"complete": 0, "leading_partial": False, "trailing_partial": False, "physical_count_reported": 0}

    if n == 1 and intervals[0].is_boundary_start and intervals[0].is_boundary_end:
        return {
            "complete": 0, "leading_partial": True, "trailing_partial": True,
            "physical_count_reported": 0,
            "note": "single interval spans the full window; both ends truncated, not a complete observed pass",
        }

    leading_partial = intervals[0].is_boundary_start
    trailing_partial = intervals[-1].is_boundary_end
    complete = n - int(leading_partial) - int(trailing_partial)
    return {
        "complete": complete,
        "leading_partial": leading_partial,
        "trailing_partial": trailing_partial,
        "physical_count_reported": complete,
    }


def access_metrics_boundary_aware(
    intervals: list[AccessInterval],
    t_start_s: float,
    t_end_s: float,
    periodic: bool,
    period_s: float | None = None,
) -> dict:
    """Boundary-aware wrapper around :func:`access_metrics`. total_access_
    time/access_fraction/longest_pass/max_gap are unaffected by the
    boundary-merge question (those quantities are correct either way from
    the raw interval list); only ``num_passes`` bookkeeping differs:

    - ``periodic=True`` (period_s required): boundary intervals are
      merged first (:func:`merge_cyclic_boundary_intervals`), then
      ``num_passes`` is the physically correct merged count.
    - ``periodic=False``: ``num_passes`` reports only *complete* interior
      passes (:func:`summarize_boundary_passes`); leading/trailing
      partial passes are flagged separately rather than counted as whole
      passes or silently merged.
    """
    if periodic:
        if period_s is None:
            raise ValueError("period_s is required when periodic=True")
        working = merge_cyclic_boundary_intervals(intervals, t_start_s, t_end_s, period_s)
        base = access_metrics(working, t_start_s, t_end_s)
        return {
            "total_access_time_s": base.total_access_time_s,
            "access_fraction": base.access_fraction,
            "num_passes": base.num_passes,
            "num_complete_passes": base.num_passes,
            "has_leading_partial": False,
            "has_trailing_partial": False,
            "longest_pass_s": base.longest_pass_s,
            "max_gap_s": base.max_gap_s,
            "boundary_treatment": "merged_cyclic",
        }
    else:
        base = access_metrics(intervals, t_start_s, t_end_s)
        summary = summarize_boundary_passes(intervals)
        return {
            "total_access_time_s": base.total_access_time_s,
            "access_fraction": base.access_fraction,
            "num_passes": summary["physical_count_reported"],
            "num_complete_passes": summary["complete"],
            "has_leading_partial": summary["leading_partial"],
            "has_trailing_partial": summary["trailing_partial"],
            "longest_pass_s": base.longest_pass_s,
            "max_gap_s": base.max_gap_s,
            "boundary_treatment": "labeled_partial_lead_trail",
        }
