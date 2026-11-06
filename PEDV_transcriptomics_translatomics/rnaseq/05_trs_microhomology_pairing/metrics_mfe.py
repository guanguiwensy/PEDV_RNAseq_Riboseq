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
import RNA, random, time
from multiprocessing import Pool
from scipy.stats import spearmanr
SU=seq.replace('T','U'); COMP=str.maketrans('ACGU','UGCA')
rc=lambda s:s.translate(COMP)[::-1]
HW=15   # window 30 nt: donor [d-14,d+15], acceptor [a-15,a+14]
def dwin(d): return SU[d-15:d+15]
def awin(a): return SU[a-16:a+14]
def energy(mode,D,A):
    return RNA.duplexfold(rc(D) if mode=='a' else D, A).energy
def row(t):
    d,accs,mode=t; D=dwin(d)
    return [energy(mode,D,awin(a)) for a in accs]
def matrix(donors,accs,mode,pool):
    donors=list(donors); accs=list(accs)
    out=pool.map(row,[(d,accs,mode) for d in donors],chunksize=4)
    return np.array(out)
def gc(s): return (s.count('G')+s.count('C'))/len(s)
if __name__=='__main__':
    NM=int(sys.argv[1]) if len(sys.argv)>1 else NPERM
    random.seed(SEED); rng=np.random.default_rng(SEED)
    m=pd.read_csv(W+'per_junction_sim_mh_ALL.tsv',sep='\t')
    m['pass_filter']=m.pass_filter.astype(bool); m['canonical_primary']=m.canonical_primary.astype(bool)
    t0=time.time()
    with Pool(26) as pool:
        # observed, all events
        for mode in ['b','a']:
            u=m[['donor','acceptor']].drop_duplicates()
            ds=sorted(u.donor.unique())
            # pair list rows by donor
            grp=u.groupby('donor').acceptor.apply(list)
            res=pool.map(row,[(d,grp[d],mode) for d in ds],chunksize=16)
            mp={(d,a):e for d,accs,r in zip(ds,[grp[d] for d in ds],res) for a,e in zip(accs,r)}
            m['mfe_'+mode]=[mp[(d,a)] for d,a in zip(m.donor,m.acceptor)]
        print('observed done',time.time()-t0,flush=True)
        m['gc_donor']=[gc(dwin(d)) for d in m.donor]; m['gc_acceptor']=[gc(awin(a)) for a in m.acceptor]
        # cofold cross-check (filtered + canonical primary)
        sel=m.pass_filter|m.canonical_primary
        def cf(dd,aa):
            fc=RNA.fold_compound(dwin(dd)+'&'+awin(aa)); return fc.mfe_dimer()[1]
        m.loc[sel,'mfe_cofold_b']=[cf(d,a) for d,a in zip(m.donor[sel],m.acceptor[sel])]
        m.to_csv(W+'rri_mfe_per_junction_ALL.tsv',sep='\t',index=False)
        print('cofold done',time.time()-t0,flush=True)
        # ---------- pair-shuffle nulls (swap chain on precomputed matrices)
        cons={'TRSL_dep_noncanonical':lambda d,a:True,'TRSL_indep_distant':lambda d,a:(a-d-1)>=1000,'TRSL_indep_local':lambda d,a:0<(a-d-1)<1000}
        res=[]; 
        def two_sided(obs,nl):
            mu=nl.mean(); return float((1+np.sum(np.abs(nl-mu)>=abs(obs-mu)))/(len(nl)+1))
        nulls={}
        for g in cons:
            d=m[(m.cls==g)&m.pass_filter]
            dn=d.donor.values; ac=d.acceptor.values; n=len(d)
            ud=sorted(set(dn)); ua=sorted(set(ac)); di={x:i for i,x in enumerate(ud)}; ai={x:i for i,x in enumerate(ua)}
            for mode in ['b','a']:
                M=matrix(ud,ua,mode,pool)
                obs=d['mfe_'+mode].mean()
                cur=list(ac); e=[M[di[x],ai[y]] for x,y in zip(dn,cur)]
                assert abs(np.mean(e)-obs)<1e-6
                stats=[]; burn=50
                rr=random.Random(SEED)
                for it in range(burn+NM):
                    for _ in range(n):
                        i=rr.randrange(n); j=rr.randrange(n)
                        if i==j or cur[i]==cur[j]: continue
                        if cons[g](dn[i],cur[j]) and cons[g](dn[j],cur[i]):
                            cur[i],cur[j]=cur[j],cur[i]
                            e[i]=M[di[dn[i]],ai[cur[i]]]; e[j]=M[di[dn[j]],ai[cur[j]]]
                    if it>=burn: stats.append(np.mean(e))
                stats=np.array(stats); nulls[(g,mode)]=stats
                res.append(dict(group=g,mode=mode,n=n,observed_mean_mfe=obs,null_mean=stats.mean(),null_p2_5=np.percentile(stats,2.5),null_p97_5=np.percentile(stats,97.5),P=two_sided(obs,stats),n_perm=NM,
                                gc_donor=d.gc_donor.mean(),gc_acceptor=d.gc_acceptor.mean()))
                print(g,mode,res[-1],time.time()-t0,flush=True)
        # ---------- canonical: position nulls
        pool_dep=np.unique(m[(m.cls=='TRSL_dep_noncanonical')&m.pass_filter].acceptor.values)
        uni=np.arange(1000,27900)
        for gname,d in [('canonical_primary5',m[m.canonical_primary]),('canonical_all118',m[m.cls=='canonical'])]:
            ud=sorted(set(d.donor)); di={x:i for i,x in enumerate(ud)}
            for mode in ['b','a']:
                Mu=matrix(ud,uni,mode,pool)       # donors x all positions
                obs=d['mfe_'+mode].mean()
                for nm,src in [('acceptor pool of TRS-L-dep noncanonical',pool_dep),('uniform genome positions 1000-27900',uni)]:
                    idx=src-1000 if nm.startswith('uniform') else np.searchsorted(uni,src)
                    dj=np.array([di[x] for x in d.donor])
                    stats=np.array([Mu[dj, rng.choice(idx,size=len(d))].mean() for _ in range(NM)])
                    res.append(dict(group=gname,mode=mode,null_type=nm,n=len(d),observed_mean_mfe=obs,null_mean=stats.mean(),null_p2_5=np.percentile(stats,2.5),null_p97_5=np.percentile(stats,97.5),P=two_sided(obs,stats),n_perm=NM,
                                gc_donor=d.gc_donor.mean(),gc_acceptor=d.gc_acceptor.mean()))
                    print(res[-1],flush=True)
        R=pd.DataFrame(res); R.to_csv(W+'separability_mfe.tsv',sep='\t',index=False)
        # Spearman MFE vs reads (filtered), and canonical-5
        sp=[]
        for g in ['canonical_primary5','TRSL_dep_noncanonical','TRSL_indep_distant','TRSL_indep_local']:
            d=m[m.canonical_primary] if g=='canonical_primary5' else m[(m.cls==g)&m.pass_filter]
            for mode in ['b','a']:
                r,p=spearmanr(d['mfe_'+mode],d.reads); sp.append(dict(group=g,mode=mode,n=len(d),rho=r,P=p))
        pd.DataFrame(sp).to_csv(W+'spearman_mfe_reads.tsv',sep='\t',index=False)
        print(pd.DataFrame(sp))
