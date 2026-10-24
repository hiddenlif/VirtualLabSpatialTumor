"""Conditional random walks, in pixels/step, on an inferred two-dimensional field."""
import numpy as np
from scipy.spatial import cKDTree
from .analysis import terrain

def simulate(r,permeability=.5,seed=101,agents=64,steps=120,step_px=8,bias=.5,paths=True):
    if permeability not in (.1,.5,1.):raise ValueError('Permeability must be 0.1, 0.5 or 1.0')
    if not isinstance(seed,int) or not 0<=seed<=1000000:raise ValueError('Invalid seed')
    if not 1<=agents<=128 or not 1<=steps<=200:raise ValueError('Simulation exceeds local limits')
    if not r['tumor'].any():return dict(region=r['id'],eligible=False,reason='No annotated tumor cells')
    rng=np.random.default_rng(seed); eligible_xy=r['xy'][r['eligible']]
    fraction,support=terrain(r,eligible_xy); candidates=eligible_xy[(fraction<.2)&support]
    if len(candidates)<5:return dict(region=r['id'],eligible=False,reason='Insufficient non-tumor entry support')
    candidates=candidates[candidates[:,0]<=np.quantile(candidates[:,0],.2)]
    pos=candidates[rng.choice(len(candidates),agents,replace=True)].copy()
    low=r['xy'].min(axis=0);high=r['xy'].max(axis=0);extent=np.maximum(high-low,1)
    tumor_tree=cKDTree(r['xy'][r['tumor']]); entered=np.zeros(agents,bool); contact=np.zeros(agents,bool)
    first=np.full(agents,-1,int); travel=np.zeros(agents); history=[pos.copy()]
    for tick in range(steps):
        _,nearest=tumor_tree.query(pos,workers=1); vector=r['xy'][r['tumor']][nearest]-pos
        vector/=np.maximum(np.linalg.norm(vector,axis=1,keepdims=True),1e-12)
        theta=rng.uniform(0,2*np.pi,agents); direction=np.c_[np.cos(theta),np.sin(theta)]+bias*vector
        direction/=np.maximum(np.linalg.norm(direction,axis=1,keepdims=True),1e-12)
        proposed=pos+step_px*direction
        u=(proposed-low)%(2*extent);proposed=low+np.where(u<=extent,u,2*extent-u)
        terrain_fraction,on_tissue=terrain(r,proposed)
        admit=on_tissue&((terrain_fraction<.5)|(rng.random(agents)<permeability))
        travel+=np.linalg.norm(proposed-pos,axis=1)*admit;pos=np.where(admit[:,None],proposed,pos)
        actual,_=terrain(r,pos);entered|=actual>=.5
        now=tumor_tree.query(pos,workers=1)[0]<=20
        first[(first<0)&now]=tick+1;contact|=now
        if paths:history.append(pos.copy())
    result=dict(region=r['id'],patient=r['patient'],eligible=True,seed=seed,permeability=permeability,
                agents=agents,steps=steps,step_px=step_px,tumor_bias=bias,entry_fraction=float(entered.mean()),
                centroid_contact_fraction=float(contact.mean()),median_first_contact_step=float(np.median(first[first>=0])) if contact.any() else None,
                mean_travel_px=float(travel.mean()),interpretation='Hypothetical trajectories; no measured speed, permeability or therapeutic effect')
    if paths:result['paths']=np.stack(history,axis=1).round(3).tolist()
    return result
