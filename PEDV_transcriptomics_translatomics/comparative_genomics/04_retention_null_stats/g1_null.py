"""Composition null (synonymous codons resampled per ORF from that ORF's own codon usage; amino-acid sequence fixed).
One simulation feeds G1 (stop/AUG/near-cognate density per frame), G2 (stop-anchored ORF length distribution) and G4 (long-ORF capacity).
Observed values are computed by the SAME code applied to the reference (R=1).
Resampled: codons wholly inside the CDS body (excluding the canonical stop codon and the ORF1a/ORF1b overlap 12592-12622 with its flanking codons).
Fixed: everything else (UTRs, stop codons, Met/Trp codons, overlap, codons already claimed by an upstream ORF)."""
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
import sys, time
from r6lib import *
NREP=int(sys.argv[1]) if len(sys.argv)>1 else 10000
CH=200; SEED=20261019
L=LREF; NT=L-2
OVL=(12592,12622)
body={}
for n,(a,b) in CDS.items():
    lo,hi=a,b-3
    if n=='ORF1a': hi=OVL[0]-1
    if n=='ORF1b': lo=OVL[1]+1
    body[n]=(lo,hi)
CSTOP={b-2 for a,b in CDS.values()}
# ---- slots
claimed=np.zeros(L+3,bool); slots={}
for n,(a,b) in CDS.items():
    lo,hi=body[n]; sl=[]
    p=a
    while p+2<=b-3:
        if p>=lo and p+2<=hi and not claimed[p:p+3].any():
            c=ref[p-1:p+2]
            if tbl.get(c) is not None and len(c)==3 and set(c)<=set("ACGT"): sl.append(p)
        p+=3
    for p in sl: claimed[p:p+3]=True
    slots[n]=np.array(sl)
usage={}   # per ORF, per aa: (codon ids, cumulative probs, slot array)
plan=[]
for n,sl in slots.items():
    cods=np.array([code[ref[p-1:p+2]] for p in sl]); aas=aa_of[cods]
    for aa in np.unique(aas):
        m=aas==aa; cc=np.unique(cods[m])
        if len(cc)<2: continue
        cnt=np.array([(cods[m]==c).sum() for c in cc],float); cum=np.cumsum(cnt/cnt.sum())
        plan.append((n,aa,sl[m],cc,cum))
print('resampled codon slots',{n:len(s) for n,s in slots.items()},'aa groups',len(plan))
# ---- masks
def trip_mask(lo,hi,f,a):  # triplet 1-based first base p in [lo,hi-2], frame (p-a)%3==f
    m=np.zeros(NT,bool); p=np.arange(lo,hi-1); p=p[(p-a)%3==f]; m[p-1]=True; return m
groups=[(n,f) for n in ORDER for f in range(3)]
GM=np.stack([trip_mask(*body[n],f,CDS[n][0]) for n,f in groups]).astype(np.float32)
GMALL=np.stack([GM[[i for i,(n,ff) in enumerate(groups) if ff==f]].sum(0) for f in range(3)])
GMAT=np.vstack([GM,GMALL])            # 21 + 3 groups
gnames=[f'{n}|{f}' for n,f in groups]+[f'ALL|{f}' for f in range(3)]
# start-region masks for family counting (start codon p1 in CDS interval, overlap excluded)
def startreg(n,rf):
    a,b=CDS[n]; p=np.arange(a,b-1)  # start triplet fully inside CDS interval
    p=p[~((p+2>=OVL[0])&(p<=OVL[1]))]
    p=p[(p-a)%3==rf]; m=np.zeros(NT,bool); m[p-1]=True; return m
SR={(n,rf):startreg(n,rf) for n in ORDER for rf in (1,2)}
TH=[8,20,41]   # tiers: aa>=8, >=20, >40 (aa>=41)
MAXAA=600
def process(A):
    """A: (R,L) int8 nucleotide codes. returns dict of arrays."""
    R=A.shape[0]
    codes=(A[:,:-2].astype(np.int16)*16+A[:,1:-1]*4+A[:,2:])
    S=isstop[codes]; Au=isaug[codes]; Nc=isnc[codes]; St=isstart[codes]
    out={'stop':S.astype(np.float32)@GMAT.T,'aug':Au.astype(np.float32)@GMAT.T,'nc':Nc.astype(np.float32)@GMAT.T}
    nfam=np.zeros((R,len(ORDER),2,len(TH)),np.int32); niso=np.zeros((R,len(ORDER),2,len(TH)),np.int32)
    hist=np.zeros((len(ORDER),2,MAXAA+1),np.int64)
    for gf in range(3):
        idx=np.arange(gf,NT,3); n=len(idx)
        sx=S[:,idx]; st=St[:,idx]; ar=np.arange(n)
        BIG=n+5
        ns=np.where(sx,ar[None,:],BIG)
        nxt=np.minimum.accumulate(ns[:,::-1],axis=1)[:,::-1]
        valid=nxt<BIG
        aa=np.where(valid,nxt-ar[None,:],0)
        nxtc=np.minimum(nxt,n-1)
        p1stop=idx[nxtc]+1
        notcanon=~np.isin(p1stop,list(CSTOP))
        cumS=np.cumsum(st,axis=1)
        prevcum=np.maximum.accumulate(np.where(sx,cumS,0),axis=1)
        prevcum=np.concatenate([np.zeros((R,1),prevcum.dtype),prevcum[:,:-1]],axis=1)
        firststart=st&((cumS-prevcum)==1)&valid
        for ri,rn in enumerate(ORDER):
            for rfi,rf in enumerate((1,2)):
                sm=SR[(rn,rf)][idx]
                if not sm.any(): continue
                ok=sm[None,:]&valid&notcanon
                for ti,t in enumerate(TH):
                    Q=st&ok&(aa>=t)
                    niso[:,ri,rfi,ti]+=Q.sum(1)
                    # families: stop-runs containing >=1 qualifying start
                    cq=np.cumsum(Q,axis=1); pq=np.maximum.accumulate(np.where(sx,cq,0),axis=1)
                    pq=np.concatenate([np.zeros((R,1),pq.dtype),pq[:,:-1]],axis=1)
                    nfam[:,ri,rfi,ti]+=(sx&((cq-pq)>0)).sum(1)
                # length histogram: family max aa = from first start of run; assigned to region of that first start
                fs=firststart&ok&(aa>=8)
                v=np.minimum(aa[fs],MAXAA)
                hist[ri,rfi]+=np.bincount(v,minlength=MAXAA+1)
    out.update(nfam=nfam,niso=niso,hist=hist); return out
def build(rng,R):
    A=np.tile(refarr.astype(np.int8),(R,1))
    for n,aa,sl,cc,cum in plan:
        u=rng.random((R,len(sl)))
        k=(u[...,None]>cum[:-1][None,None,:]).sum(-1)
        cod=cc[k]
        for j in range(3):
            A[:,sl-1+j]=np.array([enc[ord(c[j])] for c in codon_list],np.int8)[cod]
    return A
# observed
obs=process(refarr.astype(np.int8)[None,:])
np.savez(OUT+'work/null_observed.npz',gnames=gnames,**obs)
if NREP>0:
    ss=np.random.SeedSequence(SEED); kids=ss.spawn(NREP//CH)
    acc={k:[] for k in ('stop','aug','nc','nfam','niso')}; hist=0; t0=time.time()
    for i,kd in enumerate(kids):
        A=build(np.random.default_rng(kd),CH); o=process(A)
        for k in acc: acc[k].append(o[k])
        hist=hist+o['hist']
        if i%5==0: print(i,len(kids),round(time.time()-t0),'s',flush=True)
    np.savez_compressed(OUT+'work/null_replicates.npz',gnames=gnames,hist=hist,seed=SEED,nrep=NREP,chunk=CH,**{k:np.concatenate(v) for k,v in acc.items()})
    # sanity: protein unchanged
    print('done',time.time()-t0)
