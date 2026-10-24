"""Re-execute all frozen computations; retain original results untouched."""
import argparse
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from spatial_lab.data import sha
from run_study import run,save
from verify_results import verify

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    verify();run(args.out)
    original=ROOT/'reports/study-v1';exact=[]
    for name in ('regions.json','regions.csv','sensitivity.json','development-selection.json','model.json','reconstruction.json','reconstruction.csv','movement.json'):
        identical=sha(original/name)==sha(args.out/name);exact.append(dict(artifact=name,exact=identical));assert identical,name
    a=json.loads((original/'summary.json').read_text());b=json.loads((args.out/'summary.json').read_text())
    a.pop('created_at');b.pop('created_at');assert a==b
    result=dict(status='PASS',artifacts=exact,summary_equal_except_run_timestamp=True,reference_regions=140,
                reconstruction_cases=504,movement_runs=90,new_external_model_calls=0)
    save(args.out/'REPRODUCTION.json',result);print(json.dumps(result,indent=2))
