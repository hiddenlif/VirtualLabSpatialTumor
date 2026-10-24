"""Read the frozen, reduced author-deposited cell table."""

import functools
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FEATURES = [
    "CD3",
    "CD4",
    "CD45",
    "CD2",
    "CD5",
    "CD7",
    "FOXP3",
    "CD20",
    "CD68",
    "Cytokeratin",
    "aSMA",
    "CD45RO",
]
MARKERS = ["CD8", *FEATURES, "Collagen IV", "CD31", "CDX2"]
EXCLUDED = {
    "dirt",
    "undefined",
    "tumor cells / immune cells",
    "immune cells / vasculature",
}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def digest(x):
    return hashlib.sha256(
        json.dumps(x, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


@functools.lru_cache(maxsize=1)
def load():
    meta = json.loads((ROOT / "data/derived/manifest.json").read_text())
    path = ROOT / "data/derived/cells.npz"
    if sha(path) != meta["derived_sha256"]:
        raise ValueError("Derived data digest changed")
    with np.load(path, allow_pickle=False) as z:
        arrays = {k: z[k] for k in z.files}
    for a in arrays.values():
        a.flags.writeable = False
    return arrays, meta


def region(region_id):
    a, m = load()
    if region_id not in m["regions"]:
        raise ValueError("Unknown region")
    r = m["region_names"].index(region_id)
    ix = np.flatnonzero(a["region"] == r)
    labels = np.array(m["type_names"])[a["type"][ix]]
    return {
        "id": region_id,
        "patient": int(a["patient"][ix[0]]),
        "cell_id": a["cell_id"][ix],
        "xy": a["xy"][ix].copy(),
        "markers": a["markers"][ix].copy(),
        "labels": labels,
        "eligible": ~np.isin(labels, list(EXCLUDED)),
        "tumor": labels == "tumor cells",
        "cd8": labels == "CD8+ T cells",
        "neighborhood": a["neighborhood"][ix],
        "metadata": m["regions"][region_id],
    }
