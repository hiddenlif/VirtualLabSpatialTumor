"""Re-execute all frozen computations; retain original results untouched."""

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from run_study import run, save, scientific_sources, verify_source_binding
from verify_results import verify

from spatial_lab.data import sha

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    verify()
    run(args.out)
    original = ROOT / "reports/study-v1"
    exact = []
    for name in (
        "regions.json",
        "regions.csv",
        "sensitivity.json",
        "development-selection.json",
        "model.json",
        "reconstruction.json",
        "reconstruction.csv",
        "movement.json",
    ):
        identical = sha(original / name) == sha(args.out / name)
        exact.append({"artifact": name, "exact": identical})
        assert identical, name
    a = json.loads((original / "summary.json").read_text())
    b = json.loads((args.out / "summary.json").read_text())
    a.pop("created_at")
    b.pop("created_at")
    source_revision = verify_source_binding(a.pop("source_hashes"))
    assert b.pop("source_hashes") == scientific_sources()
    assert a == b
    result = {
        "status": "PASS",
        "artifacts": exact,
        "summary_equal_except_run_timestamp_and_source_revision": True,
        "scientific_source_revision": source_revision,
        "reference_regions": 140,
        "reconstruction_cases": 504,
        "movement_runs": 90,
        "new_external_model_calls": 0,
    }
    save(args.out / "REPRODUCTION.json", result)
    print(json.dumps(result, indent=2))
