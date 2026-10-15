"""Processing profiles for optical autoRIFT: grid, chip and search parameters and their checks.

Free of GDAL and hyp3_autorift imports so the parameter logic can be tested anywhere. The pixel
arithmetic mirrors what hyp3-autorift's vendored drivers and autoRIFT 2.x do with these values:

- ChipSize0X = ceil(chip_min_m / pixel_m / 4) * 4    (Geogrid / testautoRIFT)
- GridSpacingX = int(ChipSize0X * grid_m / chip_min_m) (testautoRIFT)
- chip sizes are ChipSize0X * 2**k; ChipSize0X must be a multiple of GridSpacingX (autoRIFT)
- non-zero search limits below autoRIFT's minSearch (6 px) are raised to 6 px
- the coarse ("sparse") pass is sampled every 4 * ChipSize0X pixels (sparseSearchSampleRate 4) and
  filtered with a window of 4 * r + 1 coarse nodes, r = ChipSize0X / GridSpacingX. autoRIFT clips
  that window to the array size, and a clipped even-sized window compares each node with the wrong
  neighbour, so chip sizes whose window spans more than the image are mostly rejected.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass, replace

# autoRIFT 2.x defaults this module has to agree with
AUTORIFT_MIN_SEARCH_PX = 6
AUTORIFT_SPARSE_SEARCH_SAMPLE_RATE = 4
AUTORIFT_FILT_WIDTH = 5

PREFILTERS = ("wallis_hps", "wallis", "hps")
GEOGRID_INPUTS = ("its_live", "flat")


@dataclass(frozen=True)
class AutoriftProfile:
    """
    Parameters that define an optical autoRIFT product.

    grid_m: output posting (the Geogrid DEM grid), m.
    chip_min_m, chip_max_m: smallest and largest chip, m. autoRIFT tries chip_min first and doubles
        the chip where it fails, up to chip_max; the smallest successful chip sets the local
        resolution (offset band 4).
    search_px: search half-width around the a-priori offset, image pixels (raised to 6 by autoRIFT).
    oversample: correlation-surface oversampling: an int for every chip, a {chip size in m: ratio}
        dict, or None for the hyp3-autorift driver default (128 for every chip).
    null_reference_velocity: zero the a-priori offsets. ITS_LIVE reference velocities describe
        glacier flow and are meaningless (and mis-signed in N-S for optical) for coseismic offsets.
    geogrid_inputs: "its_live" crops and resamples the ITS_LIVE parameter tile; "flat" writes a zero
        DEM, zero slopes and zero reference velocity on the target grid (no parameter file needed;
        optical Geogrid only uses the DEM grid).
    prefilter: "wallis_hps" = Wallis fill here, then autoRIFT's 5x5 high-pass (the 90 m products);
        "wallis" = Wallis fill only, nodata from its zero mask; "hps" = autoRIFT's high-pass only,
        as ITS_LIVE does for Sentinel-2.
    align_grid: put the output grid on multiples of grid_m. False reproduces the v0 products,
        whose grid is the ITS_LIVE DEM crop (scene box plus 5 km), off the grid_m lattice.
    correct_optical_dy0_sign: the hyp3-autorift driver's fix for the north-south sign of optical
        a-priori offsets (fork branch `coseis`). False reproduces v0. Irrelevant without
        reference velocities.
    """

    name: str
    grid_m: float
    chip_min_m: float
    chip_max_m: float
    search_px: int
    oversample: int | dict[float, int] | None = None
    null_reference_velocity: bool = True
    geogrid_inputs: str = "flat"
    prefilter: str = "wallis_hps"
    align_grid: bool = True
    correct_optical_dy0_sign: bool = True

    def parameters(self) -> dict[str, object]:
        """All parameters, for VERSION.json (oversample dict keys become strings)."""
        params = asdict(self)
        if isinstance(self.oversample, dict):
            params["oversample"] = {f"{k:g}": v for k, v in self.oversample.items()}
        return params


PROFILES: dict[str, AutoriftProfile] = {
    # The v0 products (October 2026): 90 m posting, 360-720 m chips, ITS_LIVE tile inputs and
    # reference velocities, Wallis fill + high-pass, oversampling 128.
    "legacy90": AutoriftProfile(
        name="legacy90",
        grid_m=90,
        chip_min_m=360,
        chip_max_m=720,
        search_px=4,
        oversample=None,
        null_reference_velocity=False,
        geogrid_inputs="its_live",
        prefilter="wallis_hps",
        align_grid=False,
        correct_optical_dy0_sign=False,
    ),
    # v0 with one change, the start of the 90 m testbed: no reference velocities and no ITS_LIVE
    # parameter tiles (flat Geogrid inputs on a grid aligned to multiples of 90 m). The search is
    # set to the 6 px that v0 used in effect (autoRIFT raised its 4 px to its minSearch).
    "g90": AutoriftProfile(
        name="g90",
        grid_m=90,
        chip_min_m=360,
        chip_max_m=720,
        search_px=6,
        oversample=None,
        null_reference_velocity=True,
        geogrid_inputs="flat",
        prefilter="wallis_hps",
    ),
    # 30 m posting with 120 m smallest chips (same 4x chip/grid ratio as legacy90, so autoRIFT's
    # outlier filters behave the same), falling back to 240 and 480 m chips where 120 m fails.
    "s2_30m": AutoriftProfile(
        name="s2_30m",
        grid_m=30,
        chip_min_m=120,
        chip_max_m=480,
        search_px=6,
        oversample={120: 64, 240: 64, 480: 128},
        null_reference_velocity=True,
        geogrid_inputs="flat",
        prefilter="wallis_hps",
    ),
}


@dataclass(frozen=True)
class PixelGeometry:
    """A profile expressed in image pixels, as hyp3-autorift and autoRIFT will use it."""

    pixel_m: float
    grid_px: int
    chip_sizes_px: tuple[int, ...]
    ratio: int
    search_px: int
    effective_search_px: int
    oversample_px: int | dict[int, int] | None
    coarse_filter_span_m: dict[int, float]
    warnings: tuple[str, ...]

    @property
    def chip_min_px(self) -> int:
        return self.chip_sizes_px[0]

    @property
    def chip_max_px(self) -> int:
        return self.chip_sizes_px[-1]


def _whole(value: float, what: str) -> int:
    rounded = round(value)
    if abs(value - rounded) > 1e-9:
        raise ValueError(f"{what} must be a whole number of pixels, got {value:g}")
    return int(rounded)


def pixel_geometry(profile: AutoriftProfile, pixel_m: float) -> PixelGeometry:
    """
    Convert a profile to image pixels and check it against autoRIFT's constraints.

    Raises ValueError for combinations autoRIFT would reject or silently change; returns warnings
    for combinations that run but behave differently from the tested regime.
    """
    if profile.prefilter not in PREFILTERS:
        raise ValueError(f"prefilter must be one of {PREFILTERS}, got {profile.prefilter!r}")
    if profile.geogrid_inputs not in GEOGRID_INPUTS:
        raise ValueError(
            f"geogrid_inputs must be one of {GEOGRID_INPUTS}, got {profile.geogrid_inputs!r}"
        )

    grid_px = _whole(
        profile.grid_m / pixel_m, f"grid_m={profile.grid_m:g} at {pixel_m:g} m pixels"
    )
    chip0 = _whole(
        profile.chip_min_m / pixel_m, f"chip_min_m={profile.chip_min_m:g} at {pixel_m:g} m pixels"
    )
    if chip0 % 4:
        raise ValueError(
            f"chip_min_m={profile.chip_min_m:g} is {chip0} px; autoRIFT chips must be multiples of 4 px "
            f"(Geogrid would silently round it up to {math.ceil(chip0 / 4) * 4} px)"
        )
    chip_max = _whole(
        profile.chip_max_m / pixel_m, f"chip_max_m={profile.chip_max_m:g} at {pixel_m:g} m pixels"
    )
    levels = math.log2(chip_max / chip0) if chip_max >= chip0 else -1
    if levels < 0 or abs(levels - round(levels)) > 1e-9:
        raise ValueError(
            f"chip_max_m must be chip_min_m times a power of two, got {profile.chip_max_m:g}"
        )
    chips = tuple(chip0 * 2**k for k in range(int(round(levels)) + 1))

    if chip0 > grid_px and chip0 % grid_px:
        raise ValueError(
            f"the smallest chip ({chip0} px = {profile.chip_min_m:g} m) must be a multiple of the grid "
            f"({grid_px} px = {profile.grid_m:g} m); allowed smallest chips at this grid: "
            f"{_allowed_chips(grid_px, pixel_m)}"
        )
    ratio = max(chip0 // grid_px, 1)

    if profile.search_px < 1:
        raise ValueError("search_px must be at least 1")
    effective_search = max(profile.search_px, AUTORIFT_MIN_SEARCH_PX)

    oversample_px: int | dict[int, int] | None
    if isinstance(profile.oversample, dict):
        oversample_px = {}
        for size_m, value in profile.oversample.items():
            size_px = _whole(size_m / pixel_m, f"oversample chip size {size_m:g} m")
            if size_px not in chips:
                raise ValueError(
                    f"oversample entry for {size_m:g} m is not one of the chip sizes {chips} px"
                )
            oversample_px[size_px] = int(value)
        missing = [c for c in chips if c not in oversample_px]
        if missing:
            raise ValueError(
                f"oversample needs an entry for every chip size; missing {missing} px"
            )
    else:
        oversample_px = profile.oversample

    # Coarse pass: every sparseSearchSampleRate * ratio nodes of the (nested) grid, i.e. every
    # 4 * chip pixels; its filter window is (FiltWidth - 1) * ratio + 1 coarse nodes.
    window = (AUTORIFT_FILT_WIDTH - 1) * ratio + 1
    spans = {c: window * AUTORIFT_SPARSE_SEARCH_SAMPLE_RATE * c * pixel_m for c in chips}

    warnings = []
    if ratio > 4:
        widest_km = max(spans.values()) / 1000
        warnings.append(
            f"chip/grid ratio {ratio} > 4: autoRIFT's overlap-adjusted outlier filter needs "
            f"{_frac_valid(ratio):.0%} of a {window}-node window to agree, and the coarse filter "
            f"spans up to {widest_km:.0f} km, wider than most composites (clipped-window bug); "
            "expect much lower coverage. A smaller chip_min_m keeps the ratio at 2-4."
        )
    if profile.search_px < AUTORIFT_MIN_SEARCH_PX:
        warnings.append(
            f"search_px={profile.search_px} is raised to autoRIFT's minSearch of {AUTORIFT_MIN_SEARCH_PX} px"
        )
    return PixelGeometry(
        pixel_m=pixel_m,
        grid_px=grid_px,
        chip_sizes_px=chips,
        ratio=ratio,
        search_px=profile.search_px,
        effective_search_px=effective_search,
        oversample_px=oversample_px,
        coarse_filter_span_m=spans,
        warnings=tuple(warnings),
    )


def _frac_valid(ratio: int, frac_valid: float = 8 / 25) -> float:
    overlap = 1 - 1 / ratio
    return frac_valid * (1 - overlap) + overlap**2


def _allowed_chips(grid_px: int, pixel_m: float, up_to_px: int = 64) -> list[str]:
    allowed = [c for c in range(4, up_to_px + 1, 4) if c <= grid_px or c % grid_px == 0]
    return [f"{c * pixel_m:g} m" for c in allowed]


def clipped_chip_sizes(geometry: PixelGeometry, width_m: float, height_m: float) -> list[int]:
    """Chip sizes (px) whose coarse outlier-filter window is wider than the image (see module doc)."""
    smaller = min(width_m, height_m)
    return [c for c, span in geometry.coarse_filter_span_m.items() if span > smaller]


def aligned_bounds(
    minx: float, miny: float, maxx: float, maxy: float, grid_m: float, pad_m: float = 0.0
) -> tuple[float, float, float, float]:
    """Expand a box (plus pad_m) outward to multiples of grid_m."""
    return (
        math.floor((minx - pad_m) / grid_m) * grid_m,
        math.floor((miny - pad_m) / grid_m) * grid_m,
        math.ceil((maxx + pad_m) / grid_m) * grid_m,
        math.ceil((maxy + pad_m) / grid_m) * grid_m,
    )


def resolve_profile(name: str, **overrides: object) -> AutoriftProfile:
    """A named profile with any non-None overrides applied; the name becomes 'custom' if changed."""
    if name not in PROFILES:
        raise ValueError(f"unknown profile {name!r}; choose from {sorted(PROFILES)}")
    base = PROFILES[name]
    changes = {
        key: value
        for key, value in overrides.items()
        if value is not None and getattr(base, key) != value
    }
    if not changes:
        return base
    return replace(base, name=f"{name}+custom", **changes)


def driver_options(profile: AutoriftProfile, geometry: PixelGeometry) -> dict[str, object]:
    """
    The hyp3-autorift generateAutoriftProduct keyword arguments a profile sets (besides the
    Geogrid/autoRIFT file inputs), as used for every path and recorded with each product version.
    """
    options: dict[str, object] = {
        "ChipSizeMinX": geometry.chip_min_px,
        "ChipSizeMinY": geometry.chip_min_px,
        "ChipSizeMaxX": geometry.chip_max_px,
        "ChipSizeMaxY": geometry.chip_max_px,
        "SearchLimitX": geometry.search_px,
        "SearchLimitY": geometry.search_px,
        "CorrectOpticalDy0Sign": profile.correct_optical_dy0_sign,
    }
    if profile.null_reference_velocity:
        options["NullReferenceVelocity"] = True
    if geometry.oversample_px is not None:
        options["OverSampleRatio"] = geometry.oversample_px
    if profile.prefilter == "wallis":
        # Already Wallis-filtered by aria_coseis: no further filtering, nodata from the zero masks
        options["PreprocessingMethods"] = ("wallis_fill", "wallis_fill")
    return options
