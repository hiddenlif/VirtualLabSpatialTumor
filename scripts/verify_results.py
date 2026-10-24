"""Check source binding and independently rebuild retained descriptive aggregates."""
import csv
import json
import math
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from spatial_lab.data import sha,digest,load
from run_study import scientific_sources

def verify():
    p=ROOT/'reports/study-v1';run=json.loads((p/'RUN.json').read_text());s=json.loads((p/'summary.json').read_text())
    assert scientific_sources()==run['scientific_sources']==s['source_hashes']
    for name,h in run['artifacts'].items():assert sha(p/name)==h,name
    a,m=load();regions=json.loads((p/'regions.json').read_text());rec=json.loads((p/'reconstruction.json').read_text());motion=json.loads((p/'movement.json').read_text())
    assert len(regions)==140 and {r['region'] for r in regions}==set(m['regions'])
    assert sum(r['eligible'] for r in regions)==123
    assert sum(m['exclusions'].values())==17831 and m['missing_final_neighborhoods']==4566
    test={5,10,15,20,25,30,35}
    assert {r['patient'] for r in rec}==test and len(rec)==504
    assert len({(r['region'],r['pattern'],r['method'],r['seed']) for r in rec})==504
    for row in s['reconstruction']:
        subset=[r for r in rec if r['method']==row['method'] and r['pattern']==row['pattern']]
        for key in ('rmse_log1p','mae_log1p','access_absolute_error','interval90_coverage'):
            patient_means=[]
            for patient in sorted(test):
                values=[r[key] for r in subset if r['patient']==patient and r[key] is not None]
                if values:patient_means.append(sum(values)/len(values))
            if not patient_means:assert row[key] is None
            else:assert math.isclose(sum(patient_means)/len(patient_means),row[key],abs_tol=1e-12)
    assert len(motion)==90 and all(r['eligible'] for r in motion)
    return dict(status='PASS',reference_regions=140,reconstruction_cases=504,movement_runs=90,
                reporting_patients=sorted(test),source_and_artifact_hashes_match=True,patient_weighted_summaries_match=True,
                scope='Integrity and aggregate recalculation; not an independent biological validation')

if __name__=='__main__':print(json.dumps(verify(),indent=2))
