"""
Combine per-path offsets into one event mosaic (dx and dy are treated separately).

Methods:
- "mean": the mean of the paths that have a value at each node.
- "aligned_mean": as aligned_weighted, but the aligned paths are averaged with equal weights.
- "aligned_total": as aligned_weighted, but each path (member) is weighted by its total error,
  1 / (local variance + s^2). s is the structured error of the member's group (its composite
  recipe; all members form one group without groups): the spread of the 10 km block medians of the
  group's aligned values on the far-field nodes. Each block median pools thousands of nodes, so
  node noise cancels and, far from the ruptures, real motion is small: what remains is the smooth
  composite error (shadows, snow, surface change), which local scatter does not see. A group with
  too little far field (structured_error) gets s = 0, so its weights are those of aligned_weighted.
  The uncertainty is 1 / sqrt(sum of weights) per component, so it includes s.
- "aligned_weighted":
  1. Align the paths. One constant offset per path is fitted by least squares to the median
     differences in all overlaps. The differences are measured on clean nodes: both paths
     valid, not interpolated, each path's smallest chip, and local scatter at most its 75th
     percentile, trimmed at 3 NMAD. The path with the most valid nodes is the reference
     (offset 0).
  2. Each node is the mean of the aligned paths weighted by 1 / local variance, the scatter in
     a WINDOW x WINDOW-node window.
  3. The uncertainty is 1 / sqrt(sum of weights) per component. The deliverable reports
     sigma_h = sqrt(sigma_dx^2 + sigma_dy^2).

Ground motion is the same for every path, so the overlap differences hold only errors, and the
alignment needs no deformation mask. It aligns the paths to each other, not to stable ground: a
bias shared by all paths remains. For aligned_weighted and aligned_mean, sigma reflects noise only,
not smooth systematic errors such as the snow and sun biases of the composites; aligned_total adds
each group's structured error, assuming the members' errors are independent (members that share a
composite are not, so its sigma is then somewhat too small).
"""

from __future__ import annotations

from typing import Any

import numpy as np

METHODS = ("mean", "aligned_weighted", "aligned_mean", "aligned_total")
WINDOW = 7
SIGMA_FLOOR_M = 0.05
MIN_OVERLAP_NODES = 500
BLOCK_M = 10_000.0  # structured error: block size
MIN_BLOCK_NODES = 21  # nodes for a block median
MIN_BLOCKS = 4  # blocks for a spread
MIN_FAR_NODES = 1000  # far-field nodes of a group


def local_sigma(
    values: np.ndarray, window: int = WINDOW, floor: float = SIGMA_FLOOR_M
) -> np.ndarray:
    """Standard deviation in a window x window neighbourhood, ignoring NaN; NaN where no value."""
    from scipy import ndimage

    valid = np.isfinite(values)
    data = np.where(valid, values, 0.0).astype(np.float64)
    count = ndimage.uniform_filter(valid.astype(np.float64), window)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = ndimage.uniform_filter(data, window) / count
        var = ndimage.uniform_filter(data * data, window) / count - mean**2
    sigma = np.sqrt(np.clip(var, 0.0, None)) + floor
    sigma[~valid] = np.nan
    return sigma.astype(np.float32)


def path_offsets(
    values: np.ndarray,
    interp: np.ndarray,
    chip: np.ndarray,
    sigmas: np.ndarray,
    min_nodes: int = MIN_OVERLAP_NODES,
) -> tuple[np.ndarray, int, list[dict[str, Any]]]:
    """
    Constant offset of each path (stacks of shape (paths, rows, cols)) for one component.

    Returns (offsets with the reference at 0, reference index, the pairwise measurements used).
    Paths without a usable overlap keep offset 0.
    """
    n = values.shape[0]
    valid = np.isfinite(values)
    reference = int(np.argmax(valid.reshape(n, -1).sum(axis=1)))
    rows, rhs, weights, pairs = [], [], [], []
    for i in range(n):
        for j in range(i + 1, n):
            both = valid[i] & valid[j]
            if both.sum() < min_nodes:
                continue
            clean = both & (interp[i] == 0) & (interp[j] == 0)
            clean &= (chip[i] == np.nanmin(chip[i][both])) & (chip[j] == np.nanmin(chip[j][both]))
            if clean.sum() < min_nodes:
                continue
            clean &= sigmas[i] <= np.nanpercentile(sigmas[i][clean], 75)
            clean &= sigmas[j] <= np.nanpercentile(sigmas[j][clean], 75)
            diff = values[i][clean] - values[j][clean]
            if diff.size < min_nodes:
                continue
            median = np.median(diff)
            nmad = 1.4826 * np.median(np.abs(diff - median))
            diff = diff[np.abs(diff - median) <= 3 * max(nmad, SIGMA_FLOOR_M)]
            if diff.size < min_nodes:
                continue
            row = np.zeros(n)
            row[i], row[j] = 1.0, -1.0
            rows.append(row)
            rhs.append(float(np.median(diff)))
            weights.append(np.sqrt(diff.size))
            pairs.append(
                {"paths": [i, j], "median_difference_m": round(rhs[-1], 4), "nodes": diff.size}
            )
    offsets = np.zeros(n)
    if rows:
        design, target, w = np.array(rows), np.array(rhs), np.array(weights)
        free = [k for k in range(n) if k != reference]
        solution, *_ = np.linalg.lstsq(design[:, free] * w[:, None], target * w, rcond=None)
        offsets[free] = solution
    return offsets, reference, pairs


def structured_error(
    values: np.ndarray, far: np.ndarray, x: np.ndarray, y: np.ndarray
) -> dict[str, Any]:
    """
    Spread (standard deviation) of the 10 km block medians of values (stack (members, rows, cols)
    of one group, aligned) on the far-field nodes; x, y are the node centres (m). The blocks lie
    on multiples of BLOCK_M. s is None with fewer than MIN_FAR_NODES far-field values or
    MIN_BLOCKS blocks of at least MIN_BLOCK_NODES.
    """
    X, Y = np.meshgrid(x, y)
    keep = far[None] & np.isfinite(values)
    v = values[keep].astype(np.float64)
    if v.size < MIN_FAR_NODES:
        return {"s_m": None, "far_values": int(v.size), "blocks": 0}
    xs = np.broadcast_to(X, values.shape)[keep]
    ys = np.broadcast_to(Y, values.shape)[keep]
    key = np.floor(xs / BLOCK_M).astype(np.int64) * 10_000_000 + np.floor(ys / BLOCK_M).astype(
        np.int64
    )
    order = np.argsort(key, kind="stable")
    _, starts, counts = np.unique(key[order], return_index=True, return_counts=True)
    vs = v[order]
    medians = [np.median(vs[a : a + n]) for a, n in zip(starts, counts) if n >= MIN_BLOCK_NODES]
    if len(medians) < MIN_BLOCKS:
        return {"s_m": None, "far_values": int(v.size), "blocks": len(medians)}
    return {"s_m": float(np.std(medians)), "far_values": int(v.size), "blocks": len(medians)}


def merge(
    dx: np.ndarray,
    dy: np.ndarray,
    interp: np.ndarray,
    chip: np.ndarray,
    method: str,
    structure: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Mosaic stacks of per-path dx, dy (m), interp_mask and chip size (shape (paths, rows, cols)).

    aligned_total needs structure: "far" (bool (rows, cols), the far-field nodes), "x" and "y"
    (node centres, m) and optionally "groups" (one label per path; default: one group).

    Returns dx, dy and, for the aligned methods, sigma_h (m) and the alignment record (for
    aligned_total with each group's structured error).
    """
    if method not in METHODS:
        raise ValueError(f"unknown mosaic method {method!r}; choose from {METHODS}")
    if method == "aligned_total" and (structure is None or "far" not in structure):
        raise ValueError("aligned_total needs the far-field nodes (structure['far'])")
    out: dict[str, Any] = {}
    if method == "mean":
        for name, values in (("dx", dx), ("dy", dy)):
            valid = np.isfinite(values)
            count = valid.sum(axis=0)
            with np.errstate(invalid="ignore", divide="ignore"):
                mean = np.where(valid, values, 0.0).astype(np.float64).sum(axis=0) / count
            mean[count == 0] = np.nan
            out[name] = mean.astype(np.float32)
        return out

    variances = []
    alignment: dict[str, Any] = {}
    for name, values in (("dx", dx), ("dy", dy)):
        sigmas = np.stack([local_sigma(v) for v in values])
        offsets, reference, pairs = path_offsets(values, interp, chip, sigmas)
        aligned = values - offsets[:, None, None]
        valid = np.isfinite(aligned)
        sig = np.where(valid, sigmas, 1.0).astype(np.float64)
        structured: dict[str, Any] = {}
        if method == "aligned_total":  # add each group's structured error to its members' variance
            groups = list(structure.get("groups") or ["all"] * len(values))
            for g in sorted(set(groups)):
                members = [k for k, name in enumerate(groups) if name == g]
                structured[g] = structured_error(
                    aligned[members], structure["far"], structure["x"], structure["y"]
                )
            s2 = np.array([(structured[g]["s_m"] or 0.0) ** 2 for g in groups])
            sig = np.sqrt(sig**2 + s2[:, None, None])
        if method == "aligned_mean":  # equal weights; variance of the mean of the members
            w = valid.astype(np.float64)
            n = w.sum(axis=0)
            with np.errstate(invalid="ignore", divide="ignore"):
                variance = (w * sig**2).sum(axis=0) / n**2
        else:
            w = np.where(valid, 1.0 / np.square(sig), 0.0)
            with np.errstate(invalid="ignore", divide="ignore"):
                variance = 1.0 / w.sum(axis=0)
        total = w.sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            merged = (np.where(valid, aligned, 0.0) * w).sum(axis=0) / total
        merged[total == 0] = np.nan
        variance[total == 0] = np.nan
        out[name] = merged.astype(np.float32)
        variances.append(variance)
        alignment[name] = {
            "reference_path_index": reference,
            "offsets_m": [round(float(o), 4) for o in offsets],
            "overlaps": pairs,
        }
        if structured:
            alignment[name]["structured_error"] = {
                g: {**r, "s_m": None if r["s_m"] is None else round(r["s_m"], 4)}
                for g, r in structured.items()
            }
    out["sigma_h"] = np.sqrt(variances[0] + variances[1]).astype(np.float32)
    out["alignment"] = alignment
    return out
