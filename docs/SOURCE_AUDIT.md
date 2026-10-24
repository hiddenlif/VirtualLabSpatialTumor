# Source audit

Audited 2026-10-24 14:52 PDT.

The actual author-deposited CSV was downloaded, parsed in full, and matched to the publisher-provided SHA-256. The audit inspected 258,385 rows and 99 columns. CellID is unique. All 140 `File Name` regions join to exactly one of 35 patients; each patient contributes four regions. Required identity, coordinate and cell-type fields have no missing values.

The retained representation includes all rows, original cell IDs, patient/region joins, X/Y coordinates, 16 selected protein channels, author cell types and final neighborhood annotations. Numeric values are stored as the parsed float64 values. The full original remains in `data/source/CRC_clusters_neighborhoods_markers.csv`; its 222,999,273 bytes are ignored by future Git to keep a large source file out of the repository.

## Exclusions and missingness

| Author annotation | Rows excluded from primary metrics |
|---|---:|
| dirt | 7,357 |
| undefined | 6,524 |
| immune cells / vasculature | 2,153 |
| tumor cells / immune cells | 1,797 |
| Total | 17,831 |

All these rows remain inspectable. A separate sensitivity condition includes the mixed tumor annotation as tumor.

There are **4,566 missing final neighborhood names**, including 726 tumor cells and 39 CD8 T cells. They remain in the data with an explicit unassigned-neighborhood sentinel. This missingness is distinct from the deliberately masked CD8 channel in the reconstruction experiment. Neighborhood labels are not model features or primary proximity inputs.

## Verified and unverified evidence

The [author deposit](https://data.mendeley.com/datasets/mpjzbtfgfr/1) supplies the actual cell table, annotations, data description and CC BY 4.0 licence. The [author repository](https://github.com/nolanlab/NeighborhoodCoordination) confirms the intended patient/core/phenotype table and explains the accompanying analysis functions. The [primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7479520/) supplies study context; no original scientific code was copied into this project.

The [TCIA collection](https://www.cancerimagingarchive.net/collection/crc_ffpe-codex_cellneighs/) describes the original multiplexed images and acquisition resolution. Those raw images, masks and cell-to-image registration were **not downloaded or independently reviewed** here. Native CSV pixel units are retained; no unverified micron conversion is applied. A point map is not displayed as a microscopy image, and agreement with source annotations is not called independent biological accuracy.

`data/derived/manifest.json` contains source and derived hashes, channel-column mappings, complete joins, exclusions, units and attribution. `scripts/prepare_data.py` downloads from the retained publisher file URL when needed, checks its hash, and refuses to overwrite derived artifacts.
