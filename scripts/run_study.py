"""Run the frozen descriptive analysis, development selection, and reporting cases."""
import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import spatial_lab
import numpy as np
import scipy
import pandas as pd
from spatial_lab.data import load,region,sha,digest
from spatial_lab.analysis import analyze,bh_adjust,sensitivity
from spatial_lab.reconstruction import fit,reconstruct
from spatial_lab.movement import simulate

def scientific_sources():
    paths=['PLAN.json','data/derived/manifest.json',*['spatial_lab/'+x for x in ('__init__.py','data.py','analysis.py','reconstruction.py','movement.py')],
           'scripts/prepare_data.py','scripts/run_study.py']
    return {p:sha(ROOT/p) for p in paths}

def save(path,value):
    with path.open('x') as f:json.dump(value,f,indent=2,allow_nan=False)

def csv_save(path,rows):
    if not rows:return
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

def patient_means(rows,key):
    return {p:float(np.mean([r[key] for r in rows if r['patient']==p and r.get(key) is not None]))
            for p in sorted(set(r['patient'] for r in rows)) if any(r.get(key) is not None for r in rows if r['patient']==p)}

def run(out):
    if out.exists():raise FileExistsError('Choose a new output directory: '+str(out))
    out.mkdir(parents=True)
    a,meta=load();plan=json.loads((ROOT/'PLAN.json').read_text())
    source=scientific_sources();baseline=[];sensitive=[]
    for index,name in enumerate(meta['region_names']):
        r=region(name);baseline.append(analyze(r,seed=17+index))
        sensitive.extend(dict(region=name,patient=r['patient'],**x) for x in sensitivity(r,seed=17+index))
    valid=[r for r in baseline if r['eligible']]
    for field in ('global','block'):
        q=bh_adjust([r[field+'_null_p'] for r in valid])
        for row,value in zip(valid,q):row[field+'_null_q']=value
    save(out/'regions.json',baseline);csv_save(out/'regions.csv',[{k:v for k,v in r.items() if not k.endswith('_counts')} for r in baseline])
    save(out/'sensitivity.json',sensitive)
    print('Reference analysis complete:',len(baseline),'regions',flush=True)
    dev=plan['reconstruction']['development_patients'];test=plan['reconstruction']['test_patients']
    dev_regions=[k for k,v in meta['regions'].items() if v['patient'] in dev]
    selection=[];candidates={alpha:fit(alpha) for alpha in (1,10,100)}
    for alpha,model in candidates.items():
        scores=[reconstruct(region(name),'ridge',pattern,17,model=model) for name in dev_regions for pattern in ('random','band')]
        selection.append(dict(method='ridge',parameter=alpha,patient_mean_rmse=float(np.mean(list(patient_means(scores,'rmse_log1p').values())))))
    for k in (5,15,30):
        scores=[reconstruct(region(name),'spatial',pattern,17,model=candidates[10],k=k) for name in dev_regions for pattern in ('random','band')]
        selection.append(dict(method='spatial',parameter=k,patient_mean_rmse=float(np.mean(list(patient_means(scores,'rmse_log1p').values())))))
    best={m:min([x for x in selection if x['method']==m],key=lambda x:x['patient_mean_rmse'])['parameter'] for m in ('ridge','spatial')}
    model=candidates[best['ridge']];model['spatial_k']=best['spatial'];model['development_patients']=dev;model['test_patients']=test
    save(out/'development-selection.json',selection);save(out/'model.json',model)
    cases=[]
    for name,v in meta['regions'].items():
        if v['patient'] not in test:continue
        r=region(name)
        for pattern in ('random','band'):
            for seed in (101,107,113):
                for method in ('median','spatial','ridge'):
                    cases.append(reconstruct(r,method,pattern,seed,model=model))
    save(out/'reconstruction.json',cases);csv_save(out/'reconstruction.csv',cases)
    print('Reconstruction reporting complete:',len(cases),'cases',flush=True)
    motion=[];motion_regions=[]
    for patient in (5,20,35):
        names=sorted(k for k,v in meta['regions'].items() if v['patient']==patient)[:2];motion_regions.extend(names)
        for name in names:
            r=region(name)
            for permeability in (.1,.5,1.):
                for seed in (101,107,113,127,131):
                    motion.append(simulate(r,permeability,seed,paths=False))
    save(out/'movement.json',motion)
    rng=np.random.default_rng(17);pm=patient_means(baseline,'access_ratio');values=np.array(list(pm.values()))
    boot=rng.choice(values,(2000,len(values)),replace=True).mean(axis=1)
    reconstruction_summary=[]
    for pattern in ('random','band'):
        for method in ('median','spatial','ridge'):
            subset=[r for r in cases if r['pattern']==pattern and r['method']==method]
            summary=dict(pattern=pattern,method=method,cases=len(subset))
            for key in ('rmse_log1p','mae_log1p','access_absolute_error','interval90_coverage'):
                means=patient_means(subset,key);summary[key]=float(np.mean(list(means.values()))) if means else None
            reconstruction_summary.append(summary)
    summary=dict(status='EXECUTED',created_at=subprocess.check_output(['date','+%Y-%m-%d %H:%M %Z'],text=True).strip(),
                 source_hashes=source,plan_sha256=digest(plan),data_sha256=meta['derived_sha256'],
                 rows=meta['rows'],eligible_cells=sum(r['eligible_count'] for r in meta['regions'].values()),
                 patients=35,regions=140,eligible_regions=len(valid),patients_with_eligible_regions=len(pm),
                 region_global_depletion_q05=sum(r['global_null_q']<.05 for r in valid),
                 region_block_depletion_q05=sum(r['block_null_q']<.05 for r in valid),
                 patient_mean_access_ratio=float(values.mean()),patient_bootstrap_mean_ci95=np.quantile(boot,[.025,.975]).tolist(),
                 patient_access_ratios=pm,reconstruction_cases=len(cases),model_selection=best,reconstruction=reconstruction_summary,
                 movement_runs=len(motion),movement_eligible=sum(x['eligible'] for x in motion),movement_regions=motion_regions,
                 environment=dict(python=sys.version.split()[0],numpy=np.__version__,scipy=scipy.__version__,pandas=pd.__version__),
                 interpretation='Descriptive sampled-tissue analysis; exploratory conditional nulls; patient-level reconstruction reporting; hypothetical movement')
    save(out/'summary.json',summary)
    save(out/'RUN.json',dict(command='python3 scripts/run_study.py --out NEW_DIRECTORY',scientific_sources=source,
                           plan_sha256=digest(plan),artifacts={p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()},
                           model_calls=0))
    print(json.dumps(summary,indent=2))
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=ROOT/'reports/study-v1');args=p.parse_args();run(args.out)
