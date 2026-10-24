# Reproduce the investigation

The complete reduced input is included under `data/derived/`; initial reproduction needs no network download. Install the three pinned dependencies in `requirements.txt` with Python 3.12+. The retained run used Python 3.13.15, NumPy 2.5.2, SciPy 1.18.0 and pandas 2.3.3 on macOS ARM64. Other platforms can differ in floating-point arithmetic; the strict replay command reports a mismatch rather than silently accepting changed results.

```sh
python3 -m unittest discover -s tests -v
python3 scripts/verify_results.py
python3 scripts/reproduce.py --out reports/my-reproduction
```

Use a new output path for each run. `run_study.py` and `reproduce.py` refuse existing output directories. `verify_results.py` checks source/data/artifact bindings, the reporting roster, and recomputes patient-weighted summaries without refitting a model.

`reproduce.py` reruns all 140 regional analyses, 980 sensitivity conditions, development selection, 504 reporting reconstructions and 90 movement runs. It checks eight substantive artifacts byte for byte and compares the summary excluding its execution timestamp. It makes zero external model calls. Ridge is genuinely refitted and its predictions are re-executed; this is not a cached model demonstration.

`PLAN.json` predates reporting outputs and specifies the source hash, analysis and population splits. `reports/study-v1/RUN.json` binds scientific source files and all initial study artifacts. Later frontend, API and documentation changes are separately covered by delivery checksums. `reports/reproduction-v1/REPRODUCTION.json` preserves the completed fresh replay.

To reconstruct the reduced input from the original CSV, use `python3 scripts/prepare_data.py --out NEW_DIRECTORY`. If absent, the original is downloaded from the publisher's retained file URL. The full 223 MB CSV must match its known SHA-256. Existing reduced data are not overwritten. New preparation timestamps differ; numeric arrays, IDs and categorical mappings should be compared independently of that audit timestamp.

Interactive computations are saved separately under `.runtime/`. Their dossier contains source/model bindings, input settings, result digests and complete outputs. They do not change the frozen study. Source-cell CSV exports contain processed observations and author annotations, not the generated replacement channel. JSON dossiers preserve those reconstruction and simulation layers separately.
