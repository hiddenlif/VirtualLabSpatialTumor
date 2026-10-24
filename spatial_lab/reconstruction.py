"""Recover a deliberately masked processed protein channel; truth stays separate."""
import functools
import json
import numpy as np
from scipy.spatial import cKDTree
from .data import ROOT, load, MARKERS, digest

def split_masks(a):
    p=a['patient']; return np.isin(p%5,[1,2,3]),p%5==4,p%5==0

def fit(alpha=10):
    a,m=load(); train,_,_=split_masks(a)
    valid=~np.isin(np.array(m['type_names'])[a['type']],m['excluded_annotations']); train &= valid
    x=np.log1p(a['markers'][train,1:13]); y=np.log1p(a['markers'][train,0])
    mu=x.mean(axis=0); sd=x.std(axis=0); sd[sd<1e-9]=1
    z=(x-mu)/sd; ym=float(y.mean())
    coef=np.linalg.solve(z.T@z+float(alpha)*np.eye(z.shape[1]),z.T@(y-ym))
    residual=y-(ym+z@coef)
    return dict(alpha=alpha,mean=mu.tolist(),scale=sd.tolist(),coefficients=coef.tolist(),intercept=ym,
                median=float(np.median(y)),high_threshold=float(np.quantile(y,.9)),
                residual_q05=float(np.quantile(residual,.05)),residual_q95=float(np.quantile(residual,.95)),
                training_patients=sorted(np.unique(a['patient'][train]).tolist()),features=MARKERS[1:13],
                target='log1p(CD8 intensity)',data_sha256=m['derived_sha256'])

@functools.lru_cache(maxsize=1)
def frozen_model():
    return json.loads((ROOT/'reports/study-v1/model.json').read_text())

def mask_channel(n,xy,pattern,seed):
    rng=np.random.default_rng(seed); count=max(1,int(n*.3)); mask=np.zeros(n,bool)
    if pattern=='random': mask[rng.choice(n,count,replace=False)]=True
    elif pattern=='band':
        start=int(rng.uniform(.2,.5)*n); mask[np.argsort(xy[:,0],kind='stable')[start:start+count]]=True
    else:raise ValueError('Unknown missingness pattern')
    return mask

def spatial_predict(xy,y,mask,k):
    observed=np.flatnonzero(~mask)
    if not len(observed):raise ValueError('No observed marker values')
    d,ix=cKDTree(xy[observed]).query(xy[mask],k=min(k,len(observed)),workers=1)
    if d.ndim==1:d=d[:,None];ix=ix[:,None]
    weights=1/np.maximum(d,1)**2
    return (y[observed[ix]]*weights).sum(axis=1)/weights.sum(axis=1)

def reconstruct(r,method='ridge',pattern='random',seed=101,model=None,k=None,include_cells=False):
    if method not in ('median','spatial','ridge'):raise ValueError('Unknown reconstruction method')
    if not isinstance(seed,int) or not 0<=seed<=1000000:raise ValueError('Seed must be an integer in 0..1000000')
    model=frozen_model() if model is None else model
    use=r['eligible']; xy=r['xy'][use]; marker=r['markers'][use]; y=np.log1p(marker[:,0]); mask=mask_channel(len(y),xy,pattern,seed)
    if method=='median': yp=np.full(mask.sum(),model['median'])
    elif method=='spatial':yp=spatial_predict(xy,y,mask,k or model.get('spatial_k',15))
    else:
        z=(np.log1p(marker[mask,1:13])-np.array(model['mean']))/np.array(model['scale'])
        yp=np.maximum(0,model['intercept']+z@np.array(model['coefficients']))
    error=yp-y[mask]; filled=y.copy();filled[mask]=yp
    tumor=r['tumor'][use]; nt=~tumor
    if tumor.any(): near=cKDTree(xy[tumor]).query(xy,workers=1)[0]<=80
    else:near=np.zeros(len(xy),bool)
    truth_high=(y>=model['high_threshold'])&nt; pred_high=(filled>=model['high_threshold'])&nt
    true_access=float(near[truth_high].mean()) if truth_high.any() else None
    pred_access=float(near[pred_high].mean()) if pred_high.any() else None
    cover=None
    if method=='ridge':
        low=yp+model['residual_q05'];high=yp+model['residual_q95']
        cover=float(((y[mask]>=low)&(y[mask]<=high)).mean())
    result=dict(region=r['id'],patient=r['patient'],method=method,pattern=pattern,seed=seed,hidden_cells=int(mask.sum()),
                rmse_log1p=float(np.sqrt(np.mean(error**2))),mae_log1p=float(np.mean(np.abs(error))),
                marker_high_threshold_log1p=model['high_threshold'],true_marker_high_access=true_access,
                reconstructed_marker_high_access=pred_access,
                access_absolute_error=abs(true_access-pred_access) if true_access is not None and pred_access is not None else None,
                interval90_coverage=cover,model_sha256=digest(model),
                target='processed CD8 intensity; marker-high is not a CD8 T-cell identity call')
    if include_cells:
        result.update(cell_ids=r['cell_id'][use].tolist(),mask=mask.tolist(),observed_log1p=np.where(mask,np.nan,y).tolist(),
                      reconstructed_log1p=filled.tolist(),validation_truth_log1p=y.tolist(),
                      uncertainty90_low=(filled+model['residual_q05']).tolist() if method=='ridge' else None,
                      uncertainty90_high=(filled+model['residual_q95']).tolist() if method=='ridge' else None)
        result['observed_log1p']=[None if x is None or not np.isfinite(x) else x for x in result['observed_log1p']]
    return result
