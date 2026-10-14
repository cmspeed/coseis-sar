"""Optical autoRIFT processing profiles: pixel geometry, autoRIFT's constraints, driver options."""

import pytest

from coseis_api import (
    AUTORIFT_PROFILES,
    aligned_bounds,
    clipped_chip_sizes,
    driver_options,
    pixel_geometry,
    resolve_profile,
)


def test_legacy90_is_the_v0_product_geometry() -> None:
    profile = AUTORIFT_PROFILES["legacy90"]
    geom = pixel_geometry(profile, 10.0)
    assert geom.grid_px == 9
    assert geom.chip_sizes_px == (36, 72)
    assert geom.ratio == 4
    assert geom.search_px == 4 and geom.effective_search_px == 6  # autoRIFT's minSearch
    assert geom.oversample_px is None  # driver default (128)
    assert geom.coarse_filter_span_m == {36: 24480.0, 72: 48960.0}
    assert any("minSearch" in w for w in geom.warnings)
    # v0 used ITS_LIVE tiles and reference velocities, the tile-crop grid and no sign fix
    assert profile.geogrid_inputs == "its_live" and not profile.null_reference_velocity
    assert not profile.align_grid
    assert driver_options(profile, geom) == {
        "ChipSizeMinX": 36,
        "ChipSizeMinY": 36,
        "ChipSizeMaxX": 72,
        "ChipSizeMaxY": 72,
        "SearchLimitX": 4,
        "SearchLimitY": 4,
        "CorrectOpticalDy0Sign": False,
    }


def test_g90_is_v0_without_reference_velocities_or_parameter_tiles() -> None:
    profile = AUTORIFT_PROFILES["g90"]
    geom = pixel_geometry(profile, 10.0)
    legacy = pixel_geometry(AUTORIFT_PROFILES["legacy90"], 10.0)
    assert (geom.grid_px, geom.chip_sizes_px) == (legacy.grid_px, legacy.chip_sizes_px)
    assert geom.effective_search_px == legacy.effective_search_px == 6
    assert profile.geogrid_inputs == "flat" and profile.align_grid
    assert profile.prefilter == AUTORIFT_PROFILES["legacy90"].prefilter
    options = driver_options(profile, geom)
    assert options["NullReferenceVelocity"] is True
    assert options["SearchLimitX"] == 6 and "OverSampleRatio" not in options


def test_s2_30m_geometry() -> None:
    profile = AUTORIFT_PROFILES["s2_30m"]
    geom = pixel_geometry(profile, 10.0)
    assert geom.grid_px == 3
    assert geom.chip_sizes_px == (12, 24, 48)
    assert geom.ratio == 4
    assert geom.oversample_px == {12: 64, 24: 64, 48: 128}
    assert geom.warnings == ()
    assert driver_options(profile, geom)["OverSampleRatio"] == {12: 64, 24: 64, 48: 128}
    assert profile.parameters()["oversample"] == {"120": 64, "240": 64, "480": 128}


def test_legacy90_on_landsat_pixels() -> None:
    geom = pixel_geometry(AUTORIFT_PROFILES["legacy90"], 30.0)
    assert (geom.grid_px, geom.chip_sizes_px, geom.ratio) == (3, (12, 24), 4)


@pytest.mark.parametrize(
    "overrides, message",
    [
        ({"grid_m": 90, "chip_min_m": 240, "chip_max_m": 480}, "multiple of the grid"),
        ({"chip_min_m": 130, "chip_max_m": 260}, "multiples of 4 px"),
        ({"chip_max_m": 1000}, "power of two"),
        ({"grid_m": 25}, "whole number of pixels"),
        ({"prefilter": "median"}, "prefilter"),
        ({"oversample": {360: 64}}, "every chip size"),
    ],
)
def test_illegal_combinations_raise(overrides: dict, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        pixel_geometry(resolve_profile("legacy90", **overrides), 10.0)


def test_high_chip_to_grid_ratio_warns() -> None:
    geom = pixel_geometry(resolve_profile("legacy90", grid_m=30), 10.0)
    assert geom.ratio == 12
    assert any("ratio 12" in w for w in geom.warnings)


def test_resolve_profile_overrides() -> None:
    # Values equal to the profile's own keep the named profile
    same = resolve_profile("g90", prefilter="wallis_hps", grid_m=90)
    assert same is AUTORIFT_PROFILES["g90"]
    changed = resolve_profile("g90", prefilter="hps")
    assert changed.name == "g90+custom" and changed.prefilter == "hps"
    assert changed.parameters() != AUTORIFT_PROFILES["g90"].parameters()
    with pytest.raises(ValueError, match="unknown profile"):
        resolve_profile("s2_10m")


def test_clipped_chip_sizes() -> None:
    legacy = pixel_geometry(AUTORIFT_PROFILES["legacy90"], 10.0)
    assert clipped_chip_sizes(legacy, 42_000, 113_000) == [72]
    assert clipped_chip_sizes(legacy, 100_000, 100_000) == []
    assert (
        clipped_chip_sizes(pixel_geometry(AUTORIFT_PROFILES["s2_30m"], 10.0), 42_000, 42_000) == []
    )


def test_aligned_bounds() -> None:
    assert aligned_bounds(375183.0, 4000001.0, 399999.0, 4024001.0, 30) == (
        375180,
        3999990,
        400020,
        4024020,
    )
    assert aligned_bounds(100.0, 100.0, 200.0, 200.0, 90, pad_m=5000) == (-4950, -4950, 5220, 5220)
