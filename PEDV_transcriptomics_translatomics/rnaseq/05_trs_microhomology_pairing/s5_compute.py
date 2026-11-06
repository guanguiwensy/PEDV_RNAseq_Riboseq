import os as _os, sys as _sys
_d = _os.path.dirname(_os.path.abspath(__file__))
while _d != _os.path.dirname(_d) and not _os.path.isdir(_os.path.join(_d, "config")):
    _d = _os.path.dirname(_d)
_sys.path.insert(0, _os.path.join(_d, "config"))
from pedv_paths import PEDV_P, PEDV_K, PEDV_R, PEDV_COV, PEDV_SERVER_ROOT, REPO_ROOT, SCRIPT_DIRS, ORIGIN
_sys.path.extend(_x for _x in SCRIPT_DIRS if _x not in _sys.path)
__REPO_FILE__ = globals().get("__REPO_FILE__") or _os.path.abspath(__file__)
_o = ORIGIN.get(_os.path.relpath(_os.path.abspath(__file__), REPO_ROOT).replace(chr(92), "/"))
if _o:
    __file__ = _o
import sys; sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from s3_lib import *
RES=W+'res/'
m=load_all()
null=np.load(W+'sim_null_dinuc.npy')*12; ids=pd.read_csv(W+'sim_null_dinuc_ids.tsv',sep='\t').junction_id.tolist(); idx={j:i for i,j in enumerate(ids)}
d=m[m.pass_filter&m.cls.isin(CLS3)]
obs=d.maxid.values*12; ii=[idx[j] for j in d.junction_id]; nl=null[ii]; n=len(d)
rng=np.random.default_rng(SEED)
draws=nl[np.arange(n)[:,None],rng.integers(0,nl.shape[1],(n,NPERM))].mean(0)
bo=boot_mean(obs,2000); dl=bo-draws[rng.integers(0,len(draws),len(bo))]
r=dict(group='noncanonical (dep + distant + local)',n=n,obs_mean_matches_of_12=obs.mean(),bg_mean_matches=draws.mean(),bg_p2_5=np.percentile(draws,2.5),bg_p97_5=np.percentile(draws,97.5),
       delta_matches=obs.mean()-draws.mean(),delta_ci_lo=np.percentile(dl,2.5),delta_ci_hi=np.percentile(dl,97.5),P_two_sided=p2(obs.mean(),draws))
pd.DataFrame([r]).to_csv(RES+'sim_merged_effect.tsv',sep='\t',index=False); print(pd.DataFrame([r]).round(4).to_string())
# one background draw per junction (pooled) for the box plot
rr=np.random.default_rng(SEED+11)
bg=nl[np.arange(n),rr.integers(0,nl.shape[1],n)]
pd.DataFrame(dict(junction_id=d.junction_id.values,cls=d.cls.values,best_shift_matches=obs,background_one_shuffle=bg)).to_csv(RES+'sim_merged_points.tsv',sep='\t',index=False)
