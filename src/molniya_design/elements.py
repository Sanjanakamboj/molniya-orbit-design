"""Classical orbital elements <-> Cartesian ECI state conversion.

Conventions
-----------
- **Frame:** ECI (Earth-Centered Inertial), right-handed, equatorial:
  X toward the reference equinox/reference direction, Z along Earth's
  rotation axis (north), Y completes the right-handed set
  (Y = Z cross X). This is the standard geocentric equatorial inertial
  frame used throughout orbital mechanics (e.g. Vallado, *Fundamentals of
  Astrodynamics and Applications*). No specific equinox epoch (e.g. J2000)
  is attached to this frame at M2 — see DESIGN.md §0.4; it is simply "the"
  inertial frame RAAN is measured in, consistent with the M1 convention.

- **Perifocal frame (PQW):** P toward perigee, W along the orbit normal
  (angular momentum direction), Q = W cross P completing the right-handed
  set, in the orbital plane, 90 deg ahead of perigee in the direction of
  motion.

- **Angles:** all classical elements (i, RAAN, argument of perigee, true
  anomaly) are documented and stored in **degrees** at the public API
  boundary and converted to **radians** immediately for internal
  trigonometry. All outputs of ``rv_to_coe`` are wrapped to [0, 360) deg
  for angles.

- **Rotation:** PQW -> ECI is the standard 3-1-3 Euler sequence, expressed
  as the *active* rotation matrix

      R = R3(RAAN) @ R1(i) @ R3(argp)

  applied to a PQW column vector, where R3(theta) is a right-handed
  rotation of a vector by +theta about the Z axis and R1(theta) is a
  right-handed rotation of a vector by +theta about the X axis:

      R3(theta) = [[cos theta, -sin theta, 0],
                   [sin theta,  cos theta, 0],
                   [0,          0,         1]]

      R1(theta) = [[1, 0,          0        ],
                   [0, cos theta, -sin theta],
                   [0, sin theta,  cos theta]]

  This is the standard Vallado/Curtis transformation and rotates a vector
  expressed in the PQW basis into the same vector expressed in the ECI
  basis (an *active* rotation of the vector's coordinates, equivalently a
  passive rotation of the frame from ECI to PQW applied in reverse). It is
  right-handed and orientation-preserving (det R = +1) for all i, RAAN,
  argp, which is confirmed by a unit test.

- **Units:** distances in km, velocities in km/s, mu in km^3/s^2 (see
  ``constants.py``), angles in degrees at the API boundary.
"""

from __future__ import annotations

import numpy as np

from .constants import MU_EARTH


def _rot3(theta_rad: float) -> np.ndarray:
    """Right-handed active rotation of a vector by +theta about Z."""
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _rot1(theta_rad: float) -> np.ndarray:
    """Right-handed active rotation of a vector by +theta about X."""
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


def coe_to_rv(
    a_km: float,
    e: float,
    i_deg: float,
    raan_deg: float,
    argp_deg: float,
    nu_deg: float,
    mu: float = MU_EARTH,
) -> tuple[np.ndarray, np.ndarray]:
    """Classical orbital elements -> ECI Cartesian state.

    Parameters
    ----------
    a_km : semi-major axis, km
    e : eccentricity, dimensionless (0 <= e < 1 for the ellipse used here)
    i_deg, raan_deg, argp_deg, nu_deg : inclination, RAAN, argument of
        perigee, true anomaly, all in degrees
    mu : gravitational parameter, km^3/s^2

    Returns
    -------
    (r_eci, v_eci) : each a (3,) numpy array, km and km/s respectively.
    """
    i = np.radians(i_deg)
    raan = np.radians(raan_deg)
    argp = np.radians(argp_deg)
    nu = np.radians(nu_deg)

    p = a_km * (1.0 - e**2)
    cos_nu, sin_nu = np.cos(nu), np.sin(nu)
    r_mag = p / (1.0 + e * cos_nu)

    r_pqw = r_mag * np.array([cos_nu, sin_nu, 0.0])
    v_pqw = np.sqrt(mu / p) * np.array([-sin_nu, e + cos_nu, 0.0])

    R = _rot3(raan) @ _rot1(i) @ _rot3(argp)

    r_eci = R @ r_pqw
    v_eci = R @ v_pqw
    return r_eci, v_eci


def rv_to_coe(
    r_eci: np.ndarray, v_eci: np.ndarray, mu: float = MU_EARTH
) -> dict:
    """ECI Cartesian state -> classical orbital elements.

    Independent of ``coe_to_rv`` (does not call it internally): derives
    elements directly from angular-momentum, eccentricity-vector, and
    node-vector geometry, per the standard algorithm (e.g. Vallado
    Algorithm 9 / Curtis Algorithm 4.1), restricted to non-circular,
    non-equatorial orbits (the Molniya baseline is neither), matching the
    stated M2 scope.

    Returns
    -------
    dict with keys: ``a_km, e, i_deg, raan_deg, argp_deg, nu_deg, p_km``.
    All angles wrapped to [0, 360) deg.
    """
    r_eci = np.asarray(r_eci, dtype=float)
    v_eci = np.asarray(v_eci, dtype=float)

    r = np.linalg.norm(r_eci)
    v = np.linalg.norm(v_eci)

    h_vec = np.cross(r_eci, v_eci)
    h = np.linalg.norm(h_vec)

    n_vec = np.cross(np.array([0.0, 0.0, 1.0]), h_vec)
    n = np.linalg.norm(n_vec)

    e_vec = (np.cross(v_eci, h_vec) / mu) - (r_eci / r)
    e = np.linalg.norm(e_vec)

    energy = v**2 / 2.0 - mu / r
    a_km = -mu / (2.0 * energy)
    p_km = h**2 / mu

    i_rad = np.arccos(np.clip(h_vec[2] / h, -1.0, 1.0))

    raan_rad = np.arccos(np.clip(n_vec[0] / n, -1.0, 1.0))
    if n_vec[1] < 0:
        raan_rad = 2 * np.pi - raan_rad

    argp_rad = np.arccos(np.clip(np.dot(n_vec, e_vec) / (n * e), -1.0, 1.0))
    if e_vec[2] < 0:
        argp_rad = 2 * np.pi - argp_rad

    nu_rad = np.arccos(np.clip(np.dot(e_vec, r_eci) / (e * r), -1.0, 1.0))
    if np.dot(r_eci, v_eci) < 0:
        nu_rad = 2 * np.pi - nu_rad

    return {
        "a_km": a_km,
        "e": e,
        "i_deg": np.degrees(i_rad) % 360.0,
        "raan_deg": np.degrees(raan_rad) % 360.0,
        "argp_deg": np.degrees(argp_rad) % 360.0,
        "nu_deg": np.degrees(nu_rad) % 360.0,
        "p_km": p_km,
    }
