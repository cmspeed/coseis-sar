"""Mosaic of per-path optical offsets: path alignment, inverse-variance weights, uncertainty."""

from typing import Any, Tuple

import pytest

import coseis_api as api

np = pytest.importorskip("numpy")
pytest.importorskip("scipy")
mosaic = api.optical_mosaic

ROWS, COLS = 60, 90


def two_paths(
    offset: Tuple[float, float] = (0.4, -0.2), noise: Tuple[float, float] = (0.1, 0.1)
) -> Tuple[Any, Any, Any, Any]:
    """
    Two paths over a fault (a 1 m dx step at column 45): path 0 covers columns 0-59 and
    path 1 columns 30-89, so they overlap in columns 30-59 (1800 nodes). Path 1 has a
    constant (dx, dy) offset and the given noise; path 0 has noise 0.1 m.
    """
    rng = np.random.default_rng(0)
    signal = np.where(np.arange(COLS) >= 45, 1.0, 0.0)[None, :].repeat(ROWS, axis=0)
    dx = np.full((2, ROWS, COLS), np.nan, np.float32)
    dy = np.full((2, ROWS, COLS), np.nan, np.float32)
    for k, (cols, shift, sigma) in enumerate(
        ((slice(0, 60), (0.0, 0.0), 0.1), (slice(30, 90), offset, noise[1]))
    ):
        shape = (ROWS, cols.stop - cols.start)
        dx[k][:, cols] = signal[:, cols] + shift[0] + rng.normal(0, sigma, shape)
        dy[k][:, cols] = shift[1] + rng.normal(0, sigma, shape)
    covered = np.isfinite(dx)
    interp = np.where(covered, 0.0, np.nan).astype(np.float32)
    chip = np.where(covered, 4.0, np.nan).astype(np.float32)
    return dx, dy, interp, chip


def test_local_sigma_ignores_gaps() -> None:
    a = np.ones((20, 20), np.float32)
    a[:, 10:] = np.nan
    s = mosaic.local_sigma(a)
    assert np.isnan(s[:, 10:]).all()
    np.testing.assert_allclose(s[:, :10], mosaic.SIGMA_FLOOR_M, atol=1e-6)


def test_aligned_weighted_removes_the_offset_between_paths() -> None:
    dx, dy, interp, chip = two_paths()
    dx[1][:5, 30:35] = np.nan  # path 0 has more valid nodes: it is the reference
    out = mosaic.merge(dx, dy, interp, chip, "aligned_weighted")
    for comp, true in (("dx", 0.4), ("dy", -0.2)):
        record = out["alignment"][comp]
        assert record["reference_path_index"] == 0
        assert record["offsets_m"][0] == 0.0
        assert record["offsets_m"][1] == pytest.approx(true, abs=0.02)
        assert record["overlaps"][0]["paths"] == [0, 1]
    # no step where path 0 ends (column 60), unlike the plain mean
    mean = mosaic.merge(dx, dy, interp, chip, "mean")
    for result, step in ((mean, 0.2), (out, 0.0)):
        jump = np.nanmedian(result["dx"][:, 60:63]) - np.nanmedian(result["dx"][:, 57:60])
        assert jump == pytest.approx(step, abs=0.05)
    # the fault step is kept
    assert np.nanmedian(out["dx"][:, 46:50] - out["dx"][:, 40:44]) == pytest.approx(1, abs=0.05)


def test_uncertainty_falls_where_paths_overlap() -> None:
    dx, dy, interp, chip = two_paths()
    out = mosaic.merge(dx, dy, interp, chip, "aligned_weighted")
    sigma = out["sigma_h"]
    # columns away from the fault, whose step raises the local scatter of dx
    single, both = np.nanmedian(sigma[:, 5:25]), np.nanmedian(sigma[:, 35:41])
    floor = mosaic.SIGMA_FLOOR_M
    assert single == pytest.approx(np.hypot(0.1 + floor, 0.1 + floor), abs=0.02)
    assert both == pytest.approx(single / np.sqrt(2), rel=0.1)
    assert np.isfinite(sigma).sum() == ROWS * COLS


def test_noisier_path_gets_less_weight() -> None:
    dx, dy, interp, chip = two_paths(offset=(0.0, 0.0), noise=(0.1, 0.4))
    out = mosaic.merge(dx, dy, interp, chip, "aligned_weighted")
    overlap = slice(30, 44)  # west of the fault, where the truth is 0
    weighted = np.nanstd(out["dy"][:, overlap])
    plain = np.nanstd(mosaic.merge(dx, dy, interp, chip, "mean")["dy"][:, overlap])
    assert weighted < 0.8 * plain


def test_paths_without_a_usable_overlap_keep_their_values() -> None:
    dx, dy, interp, chip = two_paths()
    interp[1][:] = np.where(np.isfinite(dx[1]), 1.0, np.nan)  # all interpolated: no clean nodes
    out = mosaic.merge(dx, dy, interp, chip, "aligned_weighted")
    assert out["alignment"]["dx"]["offsets_m"] == [0.0, 0.0]
    assert out["alignment"]["dx"]["overlaps"] == []


def test_mean_and_unknown_methods() -> None:
    dx, dy, interp, chip = two_paths(offset=(0.0, 0.0))
    out = mosaic.merge(dx, dy, interp, chip, "mean")
    assert set(out) == {"dx", "dy"}
    np.testing.assert_allclose(out["dx"][:, 35], (dx[0][:, 35] + dx[1][:, 35]) / 2, atol=1e-6)
    with pytest.raises(ValueError, match="unknown mosaic method"):
        mosaic.merge(dx, dy, interp, chip, "median")
