"""Milestone 1 placeholder tests.

These confirm the package imports. Real numerical tests for each
milestone live in their own files: ``tests/test_twobody.py`` (M2),
``tests/test_frames_access.py`` (M3), ``tests/test_j2.py`` (M4),
``tests/test_coverage.py`` (M5).

NOTE (transparent updates as milestones were approved and implemented):
the original M1 file guarded against "accidentally starting M2/M3/M4/M5
early" by asserting specific later-milestone attributes were absent from
the top-level package. Each guard was removed, one at a time, in the
commit where that milestone was explicitly approved and implemented:
M2 (``propagate``/``two_body_propagate``), M3 (``ground_track``/
``eci_to_ecef``), M4 (``j2_acceleration``/``j2_secular_propagate``), and
now M5 (``coverage``/``access_windows`` — ``coverage`` in particular is
the real name of the new ``src/molniya_design/coverage.py`` module, so
the guard would otherwise fail simply because the module exists and was
imported elsewhere in the test session, not because of any API misuse).
No M6-specific guard existed to remove; M6 is final polish, not new
production API surface.
"""

import molniya_design


def test_package_imports():
    assert molniya_design is not None


def test_version_string_present():
    assert isinstance(molniya_design.__version__, str)
    # bumped to 1.0.0 at M6 final packaging (was 0.1.0 through M1-M5)
    assert molniya_design.__version__ == "1.0.0"
