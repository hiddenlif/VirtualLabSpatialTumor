# Development verification

Verified 2026-10-24 14:56 PDT · local first release 0.1.0.

| Check | Result | Retained evidence |
|---|---|---|
| Scientific and local HTTP tests | 16 passed; no failures | [Final test output](../reports/final-tests.txt) |
| Scientific-source, data and study artifact hashes | Passed; reporting roster and patient-weighted means agree | [Integrity record](../reports/final-integrity.json) |
| Fresh reconstruction of the complete study | Eight substantive artifacts identical; summary identical except execution timestamp | [Reproduction record](../reports/reproduction-v1/REPRODUCTION.json) |
| JavaScript syntax | Passed | [Output](../reports/javascript-check.txt) |
| Real browser exports | JSON dossier matched API; all fields in 2,434 CSV rows matched the reduced source | [Export check](../reports/browser-export-check.json) |
| Browser workflow | Five views, both reconstruction patterns/method examples, truth toggle, sensitivity, movement, source-cell lookup, notes, history reload and exports checked | [QA record](../reports/development-verification.json) |
| Documented scientific failure | P35/reg070_A, band seed 101, ridge: 533 hidden cells; RMSE 2.044; proximity error 20.88 pp | [Reconstruction screenshot](../reports/reconstruction-screen.png) |
| Mobile breakpoint | Atlas, reconstruction and study fit a 390 × 844 px viewport without page overflow; override reset | [Mobile study screenshot](../reports/mobile-study-screen.png) |

The numerical distance implementation was compared with a direct pairwise-distance formulation. Tests also check source joins and missing annotations, patient split separation, masked-target noninterference, a small independent interpolation example, random-labelling bounds, the known BH example, unsupported terrain, reproducible movement, source preservation, export bindings and rejected requests without history mutations. This is bounded engineering and numerical evidence, not independent biological validation.

Browser logs included repeated browser-extension asynchronous-listener messages. Functional checks and actual exported files supplied the verification evidence; no claim of a clean browser-console certification is made. Browser QA was scoped to the available macOS Chrome/in-app browsers and the recorded mobile override, not an exhaustive accessibility or cross-browser audit.

All 140 original regions, 980 sensitivity rows, 504 reporting reconstructions and 90 hypothetical movement runs remain retained. The source and inference limits are documented in [the investigation](INVESTIGATION.md) and [source audit](SOURCE_AUDIT.md).

`SHA256SUMS` covers the delivered code, documentation, reduced input and retained public evidence. It excludes Git/runtime metadata, interpreter caches, itself and the large original CSV; the original CSV's publisher-verified checksum remains in the data manifest. Runtime records are separate from frozen study results. The owned preview is stopped at handoff; `run.command` restarts it with saved records intact. No repository creation or publishing was performed.
