"""Explicit spatial estimands and two conditional random-labelling nulls."""

import numpy as np
from scipy.spatial import cKDTree


def base_metrics(xy, tumor, cd8, eligible, radius=80):
    xy = np.asarray(xy, float)
    tumor = np.asarray(tumor, bool) & eligible
    cd8 = np.asarray(cd8, bool) & eligible & ~tumor
    pool = np.flatnonzero(eligible & ~tumor)
    out = {
        "eligible_cells": int(eligible.sum()),
        "tumor_cells": int(tumor.sum()),
        "cd8_cells": int(cd8.sum()),
        "radius_px": radius,
    }
    if tumor.sum() < 20 or cd8.sum() < 10 or len(pool) < 20:
        return {
            **out,
            "eligible": False,
            "reason": "Requires at least 20 tumor, 10 CD8 and 20 eligible non-tumor cells",
            "access_fraction": None,
            "available_fraction": None,
            "access_ratio": None,
            "median_distance_px": None,
        }
    d = cKDTree(xy[tumor]).query(xy[pool], workers=1)[0]
    is_cd8 = cd8[pool]
    near = d <= radius
    access = float(near[is_cd8].mean())
    available = float(near.mean())
    return {
        **out,
        "eligible": True,
        "access_fraction": access,
        "available_fraction": available,
        "access_ratio": access / available if available else None,
        "median_distance_px": float(np.median(d[is_cd8])),
    }


def analyze(r, radius=80, permutations=199, seed=17):
    if radius not in (40, 80, 120):
        raise ValueError("Radius must be 40, 80 or 120 native pixels")
    if not isinstance(permutations, int) or not 0 <= permutations <= 999:
        raise ValueError("Invalid null count")
    result = base_metrics(r["xy"], r["tumor"], r["cd8"], r["eligible"], radius)
    result.update(region=r["id"], patient=r["patient"])
    if not result["eligible"] or not permutations:
        return result
    rng = np.random.default_rng(seed)
    pool = np.flatnonzero(r["eligible"] & ~r["tumor"])
    xy = r["xy"][pool]
    near = (
        cKDTree(r["xy"][r["tumor"] & r["eligible"]]).query(xy, workers=1)[0] <= radius
    )
    cd8 = r["cd8"][pool]
    n = int(cd8.sum())
    observed = int(near[cd8].sum())
    block = np.floor((xy - r["xy"].min(axis=0)) / 400).astype(int)
    strata = []
    for b in np.unique(block, axis=0):
        ix = np.flatnonzero((block == b).all(axis=1))
        strata.append((ix, int(cd8[ix].sum())))
    global_counts = []
    local_counts = []
    for _ in range(permutations):
        global_counts.append(int(near[rng.choice(len(pool), n, replace=False)].sum()))
        local_counts.append(
            sum(
                int(near[rng.choice(ix, k, replace=False)].sum())
                for ix, k in strata
                if k
            )
        )
    expected = sum(k * float(near[ix].mean()) for ix, k in strata) / n
    result.update(
        global_null_p=(1 + int(np.sum(np.array(global_counts) <= observed)))
        / (permutations + 1),
        block_null_p=(1 + int(np.sum(np.array(local_counts) <= observed)))
        / (permutations + 1),
        block_expected_fraction=expected,
        block_access_ratio=result["access_fraction"] / expected if expected else None,
        global_null_counts=global_counts,
        block_null_counts=local_counts,
        observed_near_count=observed,
        null_replicates=permutations,
    )
    return result


def bh_adjust(values):
    p = np.asarray(values, float)
    n = len(p)
    if not n:
        return []
    order = np.argsort(p)
    sorted_q = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.clip(sorted_q, 0, 1)
    return q.tolist()


def terrain(r, points, k=20):
    use = r["eligible"]
    n = int(use.sum())
    if n < 1:
        raise ValueError("No eligible cells")
    d, ix = cKDTree(r["xy"][use]).query(points, k=min(k, n), workers=1)
    if d.ndim == 1:
        d = d[:, None]
        ix = ix[:, None]
    return r["tumor"][use][ix].mean(axis=1), d[:, 0] <= 80


def grid(r, width=64, height=48):
    lo = r["xy"].min(axis=0)
    hi = r["xy"].max(axis=0)
    xx, yy = np.meshgrid(
        np.linspace(lo[0], hi[0], width), np.linspace(lo[1], hi[1], height)
    )
    fraction, support = terrain(r, np.c_[xx.ravel(), yy.ravel()])
    return {
        "width": width,
        "height": height,
        "bounds": [lo.tolist(), hi.tolist()],
        "values": np.where(support, fraction, -1)
        .reshape(height, width)
        .round(4)
        .tolist(),
    }


def sensitivity(r, seed=17):
    rng = np.random.default_rng(seed)
    rows = []
    for radius in (40, 80, 120):
        rows.append(
            dict(
                setting=f"radius {radius} px",
                **base_metrics(r["xy"], r["tumor"], r["cd8"], r["eligible"], radius),
            )
        )
    for sd in (3, 8):
        p = r["xy"] + rng.normal(0, sd, r["xy"].shape)
        rows.append(
            dict(
                setting=f"coordinate jitter SD {sd} px",
                **base_metrics(p, r["tumor"], r["cd8"], r["eligible"]),
            )
        )
    keep = r["eligible"] & (rng.random(len(r["xy"])) >= 0.2)
    rows.append(
        dict(
            setting="retain random 80% cells",
            **base_metrics(r["xy"], r["tumor"], r["cd8"], keep),
        )
    )
    tumor = r["tumor"] | (r["labels"] == "tumor cells / immune cells")
    eligible = r["eligible"] | (r["labels"] == "tumor cells / immune cells")
    rows.append(
        dict(
            setting="mixed tumor annotations included",
            **base_metrics(r["xy"], tumor, r["cd8"], eligible),
        )
    )
    return rows
