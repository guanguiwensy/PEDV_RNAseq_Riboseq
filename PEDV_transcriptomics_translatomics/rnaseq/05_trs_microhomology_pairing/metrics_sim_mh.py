# >>> release-candidate path shim (added during repository curation; no analytical change) >>>
import os as _os, sys as _sys
_d = _os.path.dirname(_os.path.abspath(__file__))
while _d != _os.path.dirname(_d) and not _os.path.isdir(_os.path.join(_d, "config")):
    _d = _os.path.dirname(_d)
_sys.path.insert(0, _os.path.join(_d, "config"))
from pedv_paths import PEDV_P, PEDV_K, PEDV_R, PEDV_COV, PEDV_SERVER_ROOT, REPO_ROOT, SCRIPT_DIRS, ORIGIN
_sys.path.extend(_x for _x in SCRIPT_DIRS if _x not in _sys.path)  # sibling modules of other repo modules (own directory keeps priority)
__REPO_FILE__ = globals().get("__REPO_FILE__") or _os.path.abspath(__file__)  # real location of this script inside the repository (sibling-script loading must use this, not the virtual __file__)
_o = ORIGIN.get(_os.path.relpath(_os.path.abspath(__file__), REPO_ROOT).replace(chr(92), "/"))
if _o:
    __file__ = _o  # virtual ORIGINAL location inside the working tree: keeps the scripts' own dirname(__file__)-relative path logic unchanged
# <<< release-candidate path shim <<<
import sys; sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
import random
rng=np.random.default_rng(SEED); random.seed(SEED)
REF='CAACTAAACGAA'  # GI-b 64-75 = CTAAAC(67-72) +-3
assert seq[63:75]==REF, seq[63:75]
enc={c:i for i,c in enumerate('ACGT')}
S=np.array([enc.get(c,4) for c in seq],dtype=np.int8); RE=np.array([enc[c] for c in REF],dtype=np.int8)
KS=np.arange(-20,21)
def mm_vec(arr):   # matches of every 12-mer start vs REF
    n=len(arr)-11
    return sum((arr[i:i+n]==RE[i]).astype(np.int16) for i in range(12))
MM=mm_vec(S)                      # index = 0-based start
def maxid_at(a):                  # a 1-based acceptor, window body[a+k-3 .. a+k+8]
    st=(a-1)+KS-3
    ok=(st>=0)&(st<len(MM))
    v=np.full(len(KS),np.nan); v[ok]=MM[st[ok]]/12
    return v
# genome-wide precompute
ALLV=np.array([maxid_at(a) for a in range(1,GL+1)])
MAXID_G=np.nanmax(ALLV,axis=1)
print('genome-wide maxid mean',np.nanmean(MAXID_G))
core_mm=lambda s6:sum(x!=y for x,y in zip(s6,'CTAAAC'))
def best_core(a):
    best=9
    for k in KS:
        s=a-1+k
        if 0<=s<=GL-6: best=min(best,core_mm(seq[s:s+6]))
    return best
def longest_run(a,k):
    s=(a-1)+k-3; run=best=0
    for i in range(12):
        if 0<=s+i<GL and seq[s+i]==REF[i]: run+=1; best=max(best,run)
        else: run=0
    return best
m=load_pool()
rows=[]
for r in m.itertuples():
    v=ALLV[r.acceptor-1]; kb=int(KS[np.nanargmax(v)])
    rows.append((np.nanmax(v),kb,v[20],best_core(r.acceptor),longest_run(r.acceptor,kb)))
m[['maxid','best_k','id_k0','core_min_mm','longest_contig_at_best']]=pd.DataFrame(rows,index=m.index)
# ---------------- microhomology
def LR(d,a):
    R=0
    while a+R<=GL and d+1+R<=GL and seq[d+R]==seq[a-1+R]: R+=1
    L=0
    while d-1-L>=0 and a-2-L>=0 and seq[d-1-L]==seq[a-2-L]: L+=1
    return L,R
lr=[LR(d,a) for d,a in zip(m.donor,m.acceptor)]
m['L']=[x[0] for x in lr]; m['R']=[x[1] for x in lr]
m['mh_len']=m.L+m.R; m['ambig_window']=m.L+m.R+1
# ---------------- composition-preserving shuffled acceptor context (dinucleotide) null
def dinuc_shuffle(s):
    last=s[-1]; nxt={}
    for a,b in zip(s,s[1:]): nxt.setdefault(a,[]).append(b)
    verts=[v for v in nxt if v!=last]
    while True:
        le={v:random.choice(nxt[v]) for v in verts}
        ok=True
        for v in verts:
            x=v; seen=0
            while x!=last and seen<10: x=le[x] if x in le else last; seen+=1
            if x!=last: ok=False; break
        if ok: break
    lists={}
    for v,l in nxt.items():
        l=list(l)
        if v in le: l.remove(le[v])
        random.shuffle(l)
        if v in le: l.append(le[v])
        lists[v]=l
    out=[s[0]]; x=s[0]
    for _ in range(len(s)-1):
        x=lists[x].pop(0); out.append(x)
    return ''.join(out)
def ctx_maxid(ctx):   # ctx 52 nt, acceptor at idx 23
    arr=np.array([enc.get(c,4) for c in ctx],dtype=np.int8)
    mm=mm_vec(arr)
    return mm[20+KS].max()/12
NSH=200
use=m[m.pass_filter|m.canonical_primary].copy()
null=np.zeros((len(use),NSH))
for i,a in enumerate(use.acceptor.values):
    ctx=seq[a-24:a+28]
    assert len(ctx)==52
    for j in range(NSH): null[i,j]=ctx_maxid(dinuc_shuffle(ctx))
np.save(W+'sim_null_dinuc.npy',null); use[['junction_id']].to_csv(W+'sim_null_dinuc_ids.tsv',sep='\t',index=False)
m.to_csv(W+'per_junction_sim_mh_ALL.tsv',sep='\t',index=False)
print('done')
