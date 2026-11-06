"""G5: per-family retention probability under k random synonymous substitutions (k = each genome's own synonymous codon changes in the parent CDS), R replicates per (family, k-bin).
Disruption = 1 - retention. Window simulation: only codons overlapping [first start, stop+2] can matter, each candidate codon is substituted independently with probability q = k/N_candidate (Bernoulli approximation of choosing k of N_candidate; documented)."""
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
import sys,time
from synlib import *
R=int(sys.argv[1]) if len(sys.argv)>1 else 10000
SEED=20261020
f=families(); iso=isoforms()
ia=f[f.primary_category=='internal alternative-frame ORF'].copy()
ia['role']=np.where(ia.max_peptide_aa>40,'target_gt40','control_8to40')
ia=ia[ia.fixed_parent_assignment.isin(CDS)]
print(ia.role.value_counts().to_dict(),flush=True)
isog=iso.groupby('family_id')
cs=np.array([[ord(c) for c in cc] for cc in codon_list]); ENC=np.array([[enc[ord(c)] for c in cc] for cc in codon_list],np.int8)
def sim_family(fid,role,parent,qs,rng):
    d=isog.get_group(fid)
    thr=41 if role=='target_gt40' else 8
    st=sorted(d[d.peptide_length_aa>=thr].start_position.tolist())
    stop=int(d.stop_position.iloc[0]); smin=st[0]
    wlo=smin; whi=stop+2
    W0=refarr[wlo-1:whi].astype(np.int8)
    a,b=CDS[parent]; 
    codons=[]
    p=a
    while p+2<=b-3:
        if p+2>=wlo and p<=whi:
            c=code.get(ref[p-1:p+2])
            if c is not None and c in alts:
                al=alts[c]; codons.append((p-wlo,np.array([x for x,_ in al]),np.cumsum([w for _,w in al])))
        p+=3
    s0=stop-smin; nc=s0//3+1
    stpos=np.array([(x-smin)//3 for x in st])
    res={}
    for q in qs:
        if q<=0: res[q]=1.0; continue
        W=np.tile(W0,(R,1))
        for off,ac,cw in codons:
            u=rng.random(R); sel=u<q
            if not sel.any(): continue
            j=np.searchsorted(cw,u[sel]/q,side='right').clip(0,len(ac)-1)
            newc=ac[j]; rows=np.nonzero(sel)[0]
            for t in range(3):
                pos=off+t
                if 0<=pos<W.shape[1]: W[rows,pos]=ENC[newc,t]
        cod=W[:,0:3*nc].reshape(R,nc,3); cv=cod[:,:,0].astype(np.int16)*16+cod[:,:,1]*4+cod[:,:,2]
        stop_ok=isstop[cv[:,nc-1]]
        stp=isstop[cv[:,:nc-1]]
        idx=np.where(stp,np.arange(nc-1)[None,:],-1).max(1) if nc>1 else np.full(R,-1)
        ok=np.zeros(R,bool)
        for k in stpos:
            if k>=nc-1: continue
            ok|=(k>idx)&isstart[cv[:,k]]
        res[q]=float((ok&stop_ok).mean())
    return res
t0=time.time(); out=[]
ss=np.random.SeedSequence(SEED); kids=ss.spawn(len(ia))
for i,(_,F) in enumerate(ia.reset_index(drop=True).iterrows()):
    par=F.fixed_parent_assignment
    kcs=sorted(set(kbin(K[par]).tolist())); qs=[min(kc,NCAND[par])/NCAND[par] for kc in kcs]
    res=sim_family(F.family_id,F.role,par,qs,np.random.default_rng(kids[i]))
    for kc,q in zip(kcs,qs): out.append(dict(family_id=F.family_id,role=F.role,parent=par,kbin=kc,q=q,p_retain=res[q]))
    if i%20==0: print(i,len(ia),round(time.time()-t0),flush=True)
pd.DataFrame(out).to_csv(OUT+'work/g5_null_pretain_by_family_kbin.csv',index=False)
json.dump(dict(seed=SEED,R=R,kappa=kappa,n_families=len(ia),method='Bernoulli(q) per candidate codon in window; q=min(kbin,Ncand)/Ncand; k-bin width 5 (as Sup L)'),open(OUT+'work/g5_null_params.json','w'))
print('done',time.time()-t0)
