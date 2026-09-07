"""Two-body ECI propagation -> ECEF -> geocentric ground track, plus an
M4 first-order-secular-J2 ground track built from the same frame code.

The two-body path (unchanged since M3) chains M2's two-body propagator
(:mod:`molniya_design.propagation`) with M3's frame transformation
(:mod:`molniya_design.frames`) and the M2 apsis detector — it adds no new
dynamics. The J2 path (new in M4) chains :mod:`molniya_design.j2`'s
secular mean-element propagator + state reconstruction with the *same*
M3 frame transformation code, so any difference between the two ground
tracks is attributable only to the J2 secular dynamics, not to a
different Earth-rotation/frame model.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import J2_EARTH, MU_EARTH, R_EARTH
from .frames import ecef_to_geocentric_latlon, eci_to_ecef
from .j2 import SecularElements, mean_anomaly_crossing_times, propagate_elements_j2_secular, reconstruct_eci_state
from .propagation import Apsis, find_apsides, propagate


@dataclass
class GroundTrack:
    t_s: np.ndarray
    r_eci: np.ndarray  # (3, N) km
    r_ecef: np.ndarray  # (3, N) km
    lat_deg: np.ndarray  # (N,)
    lon_deg: np.ndarray  # (N,)


def compute_ground_track(
    r0: np.ndarray,
    v0: np.ndarray,
    t_span: tuple[float, float],
    dt_s: float,
    mu: float = MU_EARTH,
    theta_g0_rad: float = 0.0,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> GroundTrack:
    """Propagate (two-body, no J2) and compute the Earth-fixed ground
    track at uniform time steps ``dt_s`` over ``t_span``."""
    n = int(round((t_span[1] - t_span[0]) / dt_s)) + 1
    t_eval = np.linspace(t_span[0], t_span[1], n)
    sol = propagate(r0, v0, t_span, t_eval=t_eval, mu=mu, rtol=rtol, atol=atol)

    r_eci = sol.y[0:3, :]
    r_ecef = np.zeros_like(r_eci)
    lat = np.zeros(sol.t.shape)
    lon = np.zeros(sol.t.shape)
    for k, t in enumerate(sol.t):
        r_ecef[:, k] = eci_to_ecef(r_eci[:, k], t, theta_g0_rad)
        lat[k], lon[k] = ecef_to_geocentric_latlon(r_ecef[:, k])

    return GroundTrack(t_s=sol.t, r_eci=r_eci, r_ecef=r_ecef, lat_deg=lat, lon_deg=lon)


@dataclass
class ApsisGroundPoint:
    t_s: float
    kind: str
    lat_deg: float
    lon_deg: float
    r_km: float


def apsis_ground_points(
    r0: np.ndarray,
    v0: np.ndarray,
    t_span: tuple[float, float],
    mu: float = MU_EARTH,
    theta_g0_rad: float = 0.0,
    rtol: float = 1e-13,
    atol: float = 1e-13,
) -> list[ApsisGroundPoint]:
    """Detect apsides (M2 method) and report their Earth-fixed lat/lon."""
    apsides: list[Apsis] = find_apsides(r0, v0, t_span, mu=mu, rtol=rtol, atol=atol)
    points = []
    for a in apsides:
        r_ecef = eci_to_ecef(a.state[0:3], a.t_s, theta_g0_rad)
        lat, lon = ecef_to_geocentric_latlon(r_ecef)
        points.append(ApsisGroundPoint(t_s=a.t_s, kind=a.kind, lat_deg=lat, lon_deg=lon, r_km=a.r_km))
    return points


def longitude_separation_deg(lon1_deg: float, lon2_deg: float) -> float:
    """Signed shortest angular separation lon2 - lon1, wrapped to
    (-180, 180] deg, avoiding the +180/-180 discontinuity when comparing
    successive ground-track passes."""
    d = (lon2_deg - lon1_deg + 180.0) % 360.0 - 180.0
    if d == -180.0:
        d = 180.0
    return d


# ---------------------------------------------------------------------------
# M4: first-order secular J2 ground track (same frame/rotation model as M3)
# ---------------------------------------------------------------------------


def compute_ground_track_j2(
    elements0: SecularElements,
    t_span: tuple[float, float],
    dt_s: float,
    mu: float = MU_EARTH,
    J2: float = J2_EARTH,
    Re: float = R_EARTH,
    theta_g0_rad: float = 0.0,
) -> GroundTrack:
    """Secular-J2 mean-element ground track: propagate mean elements
    (:func:`molniya_design.j2.propagate_elements_j2_secular`), reconstruct
    ECI state at each sample (:func:`molniya_design.j2.reconstruct_eci_state`),
    then reuse the *same* M3 ECI->ECEF->lat/lon chain used by
    :func:`compute_ground_track`."""
    n = int(round((t_span[1] - t_span[0]) / dt_s)) + 1
    t_eval = np.linspace(t_span[0], t_span[1], n)

    r_eci = np.zeros((3, n))
    r_ecef = np.zeros((3, n))
    lat = np.zeros(n)
    lon = np.zeros(n)

    for k, t in enumerate(t_eval):
        elements_t = propagate_elements_j2_secular(elements0, t, J2=J2, Re=Re, mu=mu)
        r, v = reconstruct_eci_state(elements_t, mu=mu)
        r_eci[:, k] = r
        r_ecef[:, k] = eci_to_ecef(r, t, theta_g0_rad)
        lat[k], lon[k] = ecef_to_geocentric_latlon(r_ecef[:, k])

    return GroundTrack(t_s=t_eval, r_eci=r_eci, r_ecef=r_ecef, lat_deg=lat, lon_deg=lon)


def apogee_ground_points_j2(
    elements0: SecularElements,
    t_span: tuple[float, float],
    mu: float = MU_EARTH,
    J2: float = J2_EARTH,
    Re: float = R_EARTH,
    theta_g0_rad: float = 0.0,
) -> list[ApsisGroundPoint]:
    """Apogee (M=180 deg) ground points under secular J2, using the exact
    analytical crossing times from
    :func:`molniya_design.j2.mean_anomaly_crossing_times` (M4 §9's
    preferred robust method — not nearest-sample search)."""
    times = mean_anomaly_crossing_times(elements0, 180.0, t_span, J2=J2, Re=Re, mu=mu)
    points = []
    for t in times:
        elements_t = propagate_elements_j2_secular(elements0, t, J2=J2, Re=Re, mu=mu)
        r, v = reconstruct_eci_state(elements_t, mu=mu)
        r_ecef = eci_to_ecef(r, t, theta_g0_rad)
        lat, lon = ecef_to_geocentric_latlon(r_ecef)
        points.append(ApsisGroundPoint(t_s=t, kind="apogee", lat_deg=lat, lon_deg=lon, r_km=np.linalg.norm(r)))
    return points


def perigee_ground_points_j2(
    elements0: SecularElements,
    t_span: tuple[float, float],
    mu: float = MU_EARTH,
    J2: float = J2_EARTH,
    Re: float = R_EARTH,
    theta_g0_rad: float = 0.0,
) -> list[ApsisGroundPoint]:
    """Perigee (M=0 deg) ground points under secular J2 (see
    :func:`apogee_ground_points_j2`)."""
    times = mean_anomaly_crossing_times(elements0, 0.0, t_span, J2=J2, Re=Re, mu=mu)
    points = []
    for t in times:
        elements_t = propagate_elements_j2_secular(elements0, t, J2=J2, Re=Re, mu=mu)
        r, v = reconstruct_eci_state(elements_t, mu=mu)
        r_ecef = eci_to_ecef(r, t, theta_g0_rad)
        lat, lon = ecef_to_geocentric_latlon(r_ecef)
        points.append(ApsisGroundPoint(t_s=t, kind="perigee", lat_deg=lat, lon_deg=lon, r_km=np.linalg.norm(r)))
    return points
