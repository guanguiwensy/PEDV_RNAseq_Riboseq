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
from common import *
S0=np.array(['ACGT'.find(c) for c in seq],dtype=np.int8)
CANON_SITES={'S':20597,'ORF3':24718,'E':25342,'M':25597,'N':26285}
GENES=['S','ORF3','E','M','N']
CLS3=['TRSL_dep_noncanonical','TRSL_indep_distant','TRSL_indep_local']
CONS={'TRSL_dep_noncanonical':lambda d,a:np.ones(len(d),bool),
      'TRSL_indep_distant':lambda d,a:(a-d-1)>=1000,
      'TRSL_indep_local':lambda d,a:((a-d-1)>0)&((a-d-1)<1000)}
def mh_vec(d,a):
    d=np.asarray(d); a=np.asarray(a); n=len(d); tot=np.zeros(n,int)
    for direction in (1,-1):
        alive=np.ones(n,bool)
        for s in range(80):
            if direction==1: i1=d+s; i2=a-1+s
            else: i1=d-s-1; i2=a-2-s
            ok=alive&(i1>=0)&(i2>=0)&(i1<GL)&(i2<GL)
            eq=np.zeros(n,bool); eq[ok]=S0[i1[ok]]==S0[i2[ok]]
            alive=eq; tot+=eq
            if not alive.any(): break
    return tot
def chain(don,acc,cons,stat_fn,nsamp,sweeps=2,burn=20,seed=SEED):
    rng=np.random.default_rng(seed); don=np.asarray(don); cur=np.asarray(acc).copy(); n=len(don)
    out=[]; h=n//2
    for it in range(burn+nsamp*sweeps):
        p=rng.permutation(n); i=p[:h]; j=p[h:2*h]
        ai=cur[i]; aj=cur[j]
        ok=cons(don[i],aj)&cons(don[j],ai)&(ai!=aj)
        cur[i[ok]]=aj[ok]; cur[j[ok]]=ai[ok]
        if it>=burn and (it-burn)%sweeps==sweeps-1: out.append(stat_fn(don,cur))
    return np.array(out),cur
def boot_mean(x,nb=2000,seed=SEED):
    rng=np.random.default_rng(seed); x=np.asarray(x,float); n=len(x); out=[]
    ch=max(1,int(2e6//max(n,1)))
    for s in range(0,nb,ch):
        k=min(ch,nb-s); out.append(x[rng.integers(0,n,(k,n))].mean(1))
    return np.concatenate(out)
def p2(obs,nl):
    mu=nl.mean(); return float((1+np.sum(np.abs(nl-mu)>=abs(obs-mu)))/(len(nl)+1))
def holm(ps):
    ps=np.asarray(ps,float); o=np.argsort(ps); m=len(ps); adj=np.empty(m); run=0
    for r,i in enumerate(o):
        run=max(run,(m-r)*ps[i]); adj[i]=min(1,run)
    return adj
def load_all():
    m=pd.read_csv(W+'rri_mfe_per_junction_ALL.tsv',sep='\t')
    for c in ['pass_filter','canonical_primary','hotspot_27620_27682']: m[c]=m[c].astype(bool)
    m['near_leader_donor_le100']=m.donor<=100
    site=np.array(list(CANON_SITES.values()))
    m['canonical_site_acceptor']=[bool((abs(a-site)<=5).any()) for a in m.acceptor]
    m['flag_leaderdonor_at_canon_site_in_distant']=(m.cls=='TRSL_indep_distant')&m.canonical_site_acceptor&m.near_leader_donor_le100
    return m
