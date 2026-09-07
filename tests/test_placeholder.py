"""Milestone 1 placeholder tests.

These only confirm the package imports and that later-milestone APIs
(propagation, ground track, coverage) do not exist yet. Real numerical
tests against the DESIGN.md verification plan begin in Milestone 2.
"""

import molniya_design


def test_package_imports():
    assert molniya_design is not None


def test_version_string_present():
    assert isinstance(molniya_design.__version__, str)
    assert molniya_design.__version__ == "0.1.0"


def test_propagation_api_not_yet_implemented():
    """Guards against accidentally starting M2 early."""
    assert not hasattr(molniya_design, "propagate")
    assert not hasattr(molniya_design, "two_body_propagate")


def test_ground_track_api_not_yet_implemented():
    """Guards against accidentally starting M3 early."""
    assert not hasattr(molniya_design, "ground_track")
    assert not hasattr(molniya_design, "eci_to_ecef")


def test_coverage_api_not_yet_implemented():
    """Guards against accidentally starting M5 early."""
    assert not hasattr(molniya_design, "coverage")
    assert not hasattr(molniya_design, "access_windows")
