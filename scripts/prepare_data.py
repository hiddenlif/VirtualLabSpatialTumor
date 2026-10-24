"""Verify the original CSV and create a lossless selected-column representation."""
import json
import subprocess
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from spatial_lab.data import MARKERS, EXCLUDED, sha

URL = 'https://data.mendeley.com/public-files/datasets/mpjzbtfgfr/files/c24351b3-76d7-444f-9edf-0246356b0c78/file_downloaded'
EXPECTED = '416cc3926a7a900ce3b22a33be699535ca35ae85fca4585a8ba5dad6d7f3c677'

def prepare(out):
    raw = ROOT / 'data/source/CRC_clusters_neighborhoods_markers.csv'
    if not raw.exists():
        subprocess.run(['curl', '--fail', '--location', '--max-time', '300', URL, '--output', str(raw)], check=True)
    if sha(raw) != EXPECTED: raise ValueError('Source checksum mismatch')
    if (out / 'manifest.json').exists() or (out / 'cells.npz').exists(): raise FileExistsError(out)
    out.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(raw)
    assert len(df) == 258385 and len(df.columns) == 99
    keys = ['CellID', 'patients', 'File Name', 'X:X', 'Y:Y', 'ClusterName']
    assert not df[keys].isna().any().any() and df.CellID.is_unique
    assert (df.groupby('File Name').patients.nunique() == 1).all()
    assert df.patients.nunique() == 35 and df['File Name'].nunique() == 140
    region_names = sorted(df['File Name'].unique())
    type_names = sorted(df.ClusterName.unique())
    missing_neighborhoods = int(df['neighborhood name'].isna().sum())
    neighborhood_labels = df['neighborhood name'].fillna('Not assigned by authors')
    neighborhood_names = sorted(neighborhood_labels.unique())
    cols = []
    for marker in MARKERS:
        hits = [c for c in df if c.startswith(marker + ' -')]
        assert len(hits) == 1, (marker, hits)
        cols.append(hits[0])
    markers = df[cols].to_numpy(dtype=np.float64)
    assert np.isfinite(markers).all() and (markers >= 0).all()
    arrays = dict(cell_id=df.CellID.to_numpy(np.int32), patient=df.patients.to_numpy(np.int16),
                  region=pd.Categorical(df['File Name'], categories=region_names).codes.astype(np.int16),
                  type=pd.Categorical(df.ClusterName, categories=type_names).codes.astype(np.int16),
                  neighborhood=pd.Categorical(neighborhood_labels, categories=neighborhood_names).codes.astype(np.int16),
                  xy=df[['X:X', 'Y:Y']].to_numpy(np.float64), markers=markers)
    np.savez_compressed(out / 'cells.npz', **arrays)
    regions = {}
    for name, x in df.groupby('File Name', sort=True):
        eligible = ~x.ClusterName.isin(EXCLUDED)
        regions[name] = dict(patient=int(x.patients.iloc[0]), group=int(x.groups.iloc[0]), spot=str(x.spots.iloc[0]),
                             cell_count=len(x), eligible_count=int(eligible.sum()),
                             tumor_count=int((x.ClusterName == 'tumor cells').sum()), cd8_count=int((x.ClusterName == 'CD8+ T cells').sum()),
                             bounds=[[float(x['X:X'].min()), float(x['Y:Y'].min())], [float(x['X:X'].max()), float(x['Y:Y'].max())]])
    meta = dict(created_at=subprocess.check_output(['date', '+%Y-%m-%d %H:%M %Z'], text=True).strip(),
                data_doi='10.17632/mpjzbtfgfr.1', source_url=URL, source_sha256=EXPECTED, source_bytes=raw.stat().st_size,
                derived_sha256=sha(out / 'cells.npz'), rows=len(df), original_columns=99,
                marker_names=MARKERS, marker_source_columns=dict(zip(MARKERS, cols)), type_names=type_names,
                region_names=region_names, neighborhood_names=neighborhood_names, missing_final_neighborhoods=missing_neighborhoods, regions=regions,
                excluded_annotations=sorted(EXCLUDED), exclusions=df.loc[df.ClusterName.isin(EXCLUDED), 'ClusterName'].value_counts().to_dict(),
                units={'xy':'native source pixels; no registered micron calibration applied', 'markers':'author-processed fluorescence intensity, arbitrary units'},
                representation='Original row order, IDs, coordinates, selected float64 intensities and annotations; no rounding, imputation or relabelling',
                licence='CC BY 4.0', licence_url='https://creativecommons.org/licenses/by/4.0/',
                attribution='Schürch, Christian (2020), Coordinated cellular neighborhoods orchestrate antitumoral immunity at the colorectal cancer invasive front, Mendeley Data, V1, doi:10.17632/mpjzbtfgfr.1')
    (out / 'manifest.json').write_text(json.dumps(meta, indent=2))
    return meta

if __name__ == '__main__':
    import argparse
    p=argparse.ArgumentParser(); p.add_argument('--out', type=Path, default=ROOT/'data/derived'); args=p.parse_args()
    m=prepare(args.out)
    print(json.dumps({k:m[k] for k in ('rows','source_sha256','derived_sha256','exclusions')}, indent=2))
