"""Seed-chain-extend homology mapper: reference position -> genome position (per genome, alignment based; no subtraction).
refmap[i] (0-based ref) = 0-based genome pos, -1 = internal deletion in genome, -2 = outside genome extent (end-truncated)."""
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
import numpy as np
from bisect import bisect_left
from Bio import Align
K=15
_al=Align.PairwiseAligner(mode='global',match_score=2,mismatch_score=-3,open_gap_score=-8,extend_gap_score=-1,
                          end_insertion_score=0,end_deletion_score=0)
def build_index(ref,k=K):
    d={}
    for i in range(len(ref)-k+1):
        s=ref[i:i+k]; d[s]=-1 if s in d else i
    return {s:i for s,i in d.items() if i>=0}
def _lis(pairs):
    # pairs sorted by ref pos; longest chain strictly increasing in genome pos
    tails=[];tail_idx=[];prev=[-1]*len(pairs)
    for n,(r,g) in enumerate(pairs):
        j=bisect_left(tails,g)
        if j==len(tails): tails.append(g); tail_idx.append(n)
        else: tails[j]=g; tail_idx[j]=n
        prev[n]=tail_idx[j-1] if j>0 else -1
    out=[];n=tail_idx[-1] if tail_idx else -1
    while n>=0: out.append(pairs[n]); n=prev[n]
    return out[::-1]
def map_genome(ref,idx,g,k=K,max_diag_jump=2500):
    L=len(ref);G=len(g)
    anc=[]
    for j in range(G-k+1):
        i=idx.get(g[j:j+k])
        if i is not None: anc.append((i,j))
    anc.sort()
    ch=_lis(anc)
    # drop diagonal outliers (isolated spurious anchors)
    if len(ch)>=3:
        keep=[]
        for n,(r,j) in enumerate(ch):
            d=j-r
            nb=[ch[m][1]-ch[m][0] for m in (n-1,n+1) if 0<=m<len(ch)]
            if all(abs(d-x)>max_diag_jump for x in nb): continue
            keep.append((r,j))
        ch=keep
    refmap=np.full(L,-2,dtype=np.int32)
    if not ch: return refmap,0.0
    # merge into blocks (rs,re,diag)
    blocks=[]
    for r,j in ch:
        d=j-r
        if blocks and blocks[-1][2]==d and r<=blocks[-1][1]: blocks[-1][1]=max(blocks[-1][1],r+k)
        else: blocks.append([r,r+k,d])
    # trim overlapping blocks of different diag
    for b0,b1 in zip(blocks,blocks[1:]):
        if b1[0]<b0[1]: b1[0]=b0[1] if b0[1]<b1[1] else b1[1]
    def put(r0,r1,d):
        for r in range(r0,r1):
            p=r+d
            refmap[r]= p if 0<=p<G else -2
    # 5' extension
    put(0,blocks[0][0],blocks[0][2])
    for n,(rs,re,d) in enumerate(blocks):
        put(rs,re,d)
        if n+1<len(blocks):
            nrs,_,nd=blocks[n+1]
            if nd==d: put(re,nrs,d)
            else:
                rg=ref[re:nrs]; gg=g[re+d:nrs+nd]
                if len(rg)==0:
                    continue
                if len(gg)==0: refmap[re:nrs]=-1; continue
                if len(rg)>6000 or len(gg)>6000:
                    refmap[re:nrs]=-1; continue
                a=_al.align(rg,gg)[0]
                for (t0,t1),(q0,q1) in zip(*a.aligned):
                    refmap[re+t0:re+t1]=(re+d)+np.arange(q0,q1)
                # unaligned ref positions inside gap -> deletion
                seg=refmap[re:nrs]
                # positions not covered by aligned blocks remain -2 initially; set to -1
                cov=np.zeros(nrs-re,bool)
                for (t0,t1),_ in zip(*a.aligned): cov[t0:t1]=True
                seg[~cov]=-1
    last=blocks[-1]; put(last[1],L,last[2])
    # identity on mapped positions
    m=refmap>=0
    idxs=np.nonzero(m)[0]
    ident=sum(1 for i in idxs[::7] if ref[i]==g[refmap[i]])/max(1,len(idxs[::7]))
    return refmap,ident
