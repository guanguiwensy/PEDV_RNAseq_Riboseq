"""Synonymous-substitution null ingredients (same definitions as Sup L p5a: single-nt synonymous changes, no stop/Met/Trp changes, ts/tv kappa from Sup L null_params.json)."""
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
from r6lib import *
kappa=json.load(open(SL+'work/null_params.json'))['kappa']
def _ts(a,b): return {a,b} in ({'A','G'},{'C','T'})
alts={}
for s in code:
    if s in STOPS or tbl.get(s) in ('M','W'): continue
    Lx=[]
    for k in range(3):
        for nb in 'ACGT':
            if nb==s[k]: continue
            t=s[:k]+nb+s[k+1:]
            if t in STOPS or tbl.get(t)!=tbl.get(s): continue
            Lx.append((code[t],kappa if _ts(s[k],nb) else 1.0))
    if Lx: tot=sum(w for _,w in Lx); alts[code[s]]=[(c,w/tot) for c,w in Lx]
_z=np.load(OUT+'work/proj.npz',allow_pickle=True); P=_z['P']; contig=_z['contig']; ids=list(_z['ids'])
_cl=pd.read_csv(SL+'work/dedup_clusters.tsv',sep='\t')
cluster=_cl.set_index('accession').loc[ids,'cluster'].values
rep=_cl.set_index('accession').loc[ids,'is_cluster_representative'].eq('yes').values
K={}; NCAND={}
for n,(a,b) in CDS.items():
    nco=(b-a+1)//3-1; p0=np.arange(a-1,a-1+3*nco,3)
    rc=np.array([code[ref[p:p+3]] for p in p0])
    ok=contig[:,p0]&contig[:,p0+1]&(P[:,p0]<4)&(P[:,p0+1]<4)&(P[:,p0+2]<4)
    gc=P[:,p0].astype(int)*16+P[:,p0+1]*4+P[:,p0+2]
    gcs=np.where(ok,gc,0)
    syn=ok&(gc!=rc[None,:])&(aa_of[gcs]==aa_of[rc][None,:])&(~isstop[gcs])
    K[n]=syn.sum(1); NCAND[n]=sum(1 for c in rc if c in alts)
def kbin(k): return 5*np.round(np.asarray(k)/5).astype(int)
