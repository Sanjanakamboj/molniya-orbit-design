"""Cartesian two-body equations of motion (no perturbations).

State vector convention
------------------------
    y = [rx, ry, rz, vx, vy, vz]

Units: km, km/s (ECI frame, see elements.py for the frame convention).
Time in seconds.

Dynamics
--------
    rddot = -mu * r / |r|^3

J2 and any other perturbation is explicitly **not** included here (M2
scope guard); it is added starting at M4.
"""

from __future__ import annotations

import numpy as np

from .constants import MU_EARTH


def two_body_eom(t: float, y: np.ndarray, mu: float = MU_EARTH) -> np.ndarray:
    """Right-hand side of the two-body Cartesian ODE, y' = f(t, y).

    Parameters
    ----------
    t : time, s (unused — the two-body field is time-invariant, but the
        signature matches what scipy.integrate.solve_ivp requires).
    y : state [rx, ry, rz, vx, vy, vz], km and km/s.
    mu : gravitational parameter, km^3/s^2.

    Returns
    -------
    dydt : [vx, vy, vz, ax, ay, az], km/s and km/s^2.
    """
    r = y[0:3]
    v = y[3:6]
    r_mag = np.linalg.norm(r)
    a = -mu * r / r_mag**3
    return np.concatenate([v, a])


def specific_energy(r: np.ndarray, v: np.ndarray, mu: float = MU_EARTH) -> float:
    """epsilon = v^2/2 - mu/r, km^2/s^2."""
    r_mag = np.linalg.norm(r)
    v_mag = np.linalg.norm(v)
    return v_mag**2 / 2.0 - mu / r_mag


def specific_angular_momentum(r: np.ndarray, v: np.ndarray) -> np.ndarray:
    """h = r x v, km^2/s (vector)."""
    return np.cross(r, v)


def eccentricity_vector(r: np.ndarray, v: np.ndarray, mu: float = MU_EARTH) -> np.ndarray:
    """e_vec = (v x h)/mu - r/|r|, dimensionless (vector, points to perigee)."""
    h = specific_angular_momentum(r, v)
    r_mag = np.linalg.norm(r)
    return np.cross(v, h) / mu - r / r_mag
