"""Milestone 1 placeholder tests.

These confirm the package imports and that later-milestone APIs (ground
track, J2 propagation, coverage) do not exist yet. Real numerical two-body
propagation/element-conversion tests (Milestone 2) live in
``tests/test_twobody.py``.

NOTE (M2 update, transparent): the original M1 file also guarded against
"accidentally starting M2 early" by asserting
``molniya_design.propagate``/``two_body_propagate`` were absent. M2 has now
been explicitly approved and implemented (see ``src/molniya_design/
propagation.py``), so that specific guard is obsolete and has been removed
here; the M3/M4/M5 guards below are unchanged and still enforced.
"""

import molniya_design


def test_package_imports():
    assert molniya_design is not None


def test_version_string_present():
    assert isinstance(molniya_design.__version__, str)
    assert molniya_design.__version__ == "0.1.0"


def test_ground_track_api_not_yet_implemented():
    """Guards against accidentally starting M3 early."""
    assert not hasattr(molniya_design, "ground_track")
    assert not hasattr(molniya_design, "eci_to_ecef")


def test_coverage_api_not_yet_implemented():
    """Guards against accidentally starting M5 early."""
    assert not hasattr(molniya_design, "coverage")
    assert not hasattr(molniya_design, "access_windows")


def test_j2_api_not_yet_implemented():
    """Guards against accidentally starting M4 early."""
    assert not hasattr(molniya_design, "j2_acceleration")
    assert not hasattr(molniya_design, "j2_secular_propagate")
