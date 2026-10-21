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


def test_aligned_mean_aligns_like_weighted_but_weighs_equally() -> None:
    dx, dy, interp, chip = two_paths(offset=(0.4, -0.2), noise=(0.1, 0.4))
    weighted = mosaic.merge(dx, dy, interp, chip, "aligned_weighted")
    plain = mosaic.merge(dx, dy, interp, chip, "aligned_mean")
    assert plain["alignment"]["dx"]["offsets_m"][1] == pytest.approx(0.4, abs=0.05)
    overlap = slice(30, 44)  # truth 0 west of the fault, after alignment
    assert np.nanstd(plain["dy"][:, overlap]) > np.nanstd(weighted["dy"][:, overlap])
    assert np.isfinite(plain["sigma_h"]).all()


def patchy_pair(seed: int = 3) -> Tuple[Any, Any, Any, Any, Any, Any]:
    """
    Two members over the same 60 x 90 km at 1 km nodes, true motion 0 everywhere. Member 0 is
    rough (0.3 m node noise) but unbiased; member 1 is smooth (0.1 m noise) but carries patches:
    one constant per 30 km square, between -1 and 1 m. Returns the stacks, the patch values and
    the node centres.
    """
    rng = np.random.default_rng(seed)
    patches = rng.uniform(-1.0, 1.0, (2, 3))
    patch_field = np.kron(patches, np.ones((30, 30)))
    dx = np.stack([rng.normal(0, 0.3, (60, 90)), patch_field + rng.normal(0, 0.1, (60, 90))])
    dy = np.stack([rng.normal(0, 0.3, (60, 90)), patch_field + rng.normal(0, 0.1, (60, 90))])
    interp = np.zeros_like(dx)
    chip = np.full_like(dx, 4.0)
    x = 400_000.0 + (np.arange(90) + 0.5) * 1000
    y = 4_000_000.0 - (np.arange(60) + 0.5) * 1000
    return dx.astype(np.float32), dy.astype(np.float32), interp, chip, patches, (x, y)


def test_structured_error_measures_patches_not_noise() -> None:
    dx, _, _, _, patches, (x, y) = patchy_pair()
    far = np.ones((60, 90), bool)
    patchy = mosaic.structured_error(dx[1:2], far, x, y)
    rough = mosaic.structured_error(dx[0:1], far, x, y)
    assert patchy["blocks"] == 54  # 6 x 9 blocks of 10 km
    assert patchy["s_m"] == pytest.approx(np.std(patches), abs=0.05)
    assert rough["s_m"] < 0.06  # 0.3 m node noise averages out in 100-node blocks
    small = np.zeros((60, 90), bool)
    small[:5, :5] = True  # 25 far-field nodes
    assert mosaic.structured_error(dx[1:2], small, x, y)["s_m"] is None


def test_aligned_total_down_weights_the_patchy_member() -> None:
    dx, dy, interp, chip, _, (x, y) = patchy_pair()
    structure = {"far": np.ones((60, 90), bool), "x": x, "y": y, "groups": ["rough", "patchy"]}
    local = mosaic.merge(dx, dy, interp, chip, "aligned_weighted")
    total = mosaic.merge(dx, dy, interp, chip, "aligned_total", structure)
    for comp in ("dx", "dy"):  # the truth is 0: the merged values are the error
        assert np.sqrt(np.nanmean(total[comp] ** 2)) < 0.7 * np.sqrt(np.nanmean(local[comp] ** 2))
        record = total["alignment"][comp]["structured_error"]
        assert record["patchy"]["s_m"] > 5 * record["rough"]["s_m"]
    # the uncertainty now includes the structured error
    assert np.nanmedian(total["sigma_h"]) > np.nanmedian(local["sigma_h"])
    assert "structured_error" not in local["alignment"]["dx"]


def test_aligned_total_without_a_far_field_is_aligned_weighted() -> None:
    dx, dy, interp, chip, _, (x, y) = patchy_pair()
    structure = {"far": np.zeros((60, 90), bool), "x": x, "y": y, "groups": ["a", "b"]}
    local = mosaic.merge(dx, dy, interp, chip, "aligned_weighted")
    total = mosaic.merge(dx, dy, interp, chip, "aligned_total", structure)
    for name in ("dx", "dy", "sigma_h"):
        np.testing.assert_allclose(total[name], local[name], atol=1e-6)
    assert total["alignment"]["dx"]["structured_error"]["a"]["s_m"] is None
    with pytest.raises(ValueError, match="far-field nodes"):
        mosaic.merge(dx, dy, interp, chip, "aligned_total")
