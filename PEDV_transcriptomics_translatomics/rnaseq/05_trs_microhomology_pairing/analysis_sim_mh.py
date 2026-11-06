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
import random, json
rng=np.random.default_rng(SEED); random.seed(SEED)
m=pd.read_csv(W+'per_junction_sim_mh_ALL.tsv',sep='\t')
m['pass_filter']=m.pass_filter.astype(bool); m['canonical_primary']=m.canonical_primary.astype(bool)
print(m[m.pass_filter].groupby('cls').size().to_dict(), m.canonical_primary.sum())
CDS={'ORF1a':(293,12622),'ORF1b':(12592,20613),'S':(20610,24758),'ORF3':(24758,25033),'E':(25361,25600),'M':(25608,26288),'N':(26300,27625)}
def region(pos):
    h=[(s,g) for g,(s,e) in CDS.items() if s<=pos<=e]
    if h: return max(h)[1]
    return "5'UTR" if pos<293 else ("3'UTR" if pos>27625 else 'intergenic')
regpos={}
for p in range(1,GL+1): regpos.setdefault(region(p),[]).append(p)
regpos={k:np.array(v) for k,v in regpos.items()}
# recompute per-position maxid (exact same as pipeline)
enc={c:i for i,c in enumerate('ACGT')}
S=np.array([enc.get(c,4) for c in seq],dtype=np.int8); RE=np.array([enc[c] for c in 'CAACTAAACGAA'],dtype=np.int8)
n=len(S)-11; MM=sum((S[i:i+n]==RE[i]).astype(np.int16) for i in range(12))
KS=np.arange(-20,21)
def mid(a):
    st=(a-1)+KS-3; ok=(st>=0)&(st<len(MM)); v=np.full(41,np.nan); v[ok]=MM[st[ok]]/12; return v
ALLV=np.array([mid(a) for a in range(1,GL+1)]); MAXG=np.nanmax(ALLV,axis=1); K0G=ALLV[:,20]
null=np.load(W+'sim_null_dinuc.npy'); ids=pd.read_csv(W+'sim_null_dinuc_ids.tsv',sep='\t').junction_id.tolist()
idx={j:i for i,j in enumerate(ids)}
NR=NPERM
def two_sided(obs,nulld):
    mu=nulld.mean(); return float((1+np.sum(np.abs(nulld-mu)>=abs(obs-mu)))/(len(nulld)+1))
res=[]
groups={'canonical_primary5':m[m.canonical_primary],
        'canonical_all118':m[m.cls=='canonical'],
        'TRSL_dep_noncanonical':m[(m.cls=='TRSL_dep_noncanonical')&m.pass_filter],
        'TRSL_indep_distant':m[(m.cls=='TRSL_indep_distant')&m.pass_filter],
        'TRSL_indep_local':m[(m.cls=='TRSL_indep_local')&m.pass_filter]}
# ---- similarity: null A dinuc-shuffled context, null B region-matched genome positions
for g,d in groups.items():
    obs=d.maxid.mean(); obs0=d.id_k0.mean()
    ii=[idx[j] for j in d.junction_id] if all(j in idx for j in d.junction_id) else None
    if ii is not None:
        nl=null[ii]  # n x NSH
        draws=nl[np.arange(len(ii))[:,None], rng.integers(0,nl.shape[1],(len(ii),NR))].mean(0)
    rp=np.array([rng.choice(regpos[r],size=NR) for r in d.acceptor_region]) # n x NR
    if g.startswith('canonical'):   # exclude canonical +-10
        pass
    dB=MAXG[rp-1].mean(0)
    dB0=K0G[rp-1]; dB0=np.nanmean(dB0,axis=0)
    res.append(dict(metric='TRS-L 12nt best-shift identity (max over -20..+20)',group=g,n=len(d),observed_mean=obs,
        null_dinucshuffle_mean=draws.mean(),null_dinuc_p2_5=np.percentile(draws,2.5),null_dinuc_p97_5=np.percentile(draws,97.5),P_dinuc=two_sided(obs,draws),
        null_regionpos_mean=dB.mean(),null_reg_p2_5=np.percentile(dB,2.5),null_reg_p97_5=np.percentile(dB,97.5),P_regionpos=two_sided(obs,dB),
        observed_id_at_k0=obs0,null_reg_k0_mean=np.nanmean(dB0)))
# per-canonical site
pc=groups['canonical_primary5'][['canonical_gene','donor','acceptor','reads','maxid','best_k','id_k0','core_min_mm','longest_contig_at_best']]
print(pc.to_string())
print('genome-wide maxid: mean %.3f  sd %.3f  q95 %.3f'%(np.nanmean(MAXG),np.nanstd(MAXG),np.nanpercentile(MAXG,95)))
for _,r in pc.iterrows():
    print(r.canonical_gene,'percentile of maxid in genome-wide:',round((MAXG<r.maxid).mean()*100+(MAXG==r.maxid).mean()*50,1),' id_k0 pct:',round((K0G<r.id_k0).mean()*100,1))
# ---- microhomology: swap-chain pair permutation within class
def Rlen(d,a):
    R=0
    while a+R<=GL and d+1+R<=GL and seq[d+R]==seq[a-1+R]: R+=1
    L=0
    while d-1-L>=0 and a-2-L>=0 and seq[d-1-L]==seq[a-2-L]: L+=1
    return L+R
def swap_chain(don,acc,constraint,nperm,thin=1,burn=50):
    n=len(don); acc=list(acc); don=list(don)
    stats=[]; cur=acc[:]
    mh=[Rlen(d,a) for d,a in zip(don,cur)]
    total=burn+nperm*thin
    for it in range(total):
        for _ in range(n):
            i=random.randrange(n); j=random.randrange(n)
            if i==j: continue
            ai,aj=cur[i],cur[j]
            if constraint(don[i],aj) and constraint(don[j],ai) and (ai!=aj):
                cur[i],cur[j]=aj,ai
                mh[i]=Rlen(don[i],aj); mh[j]=Rlen(don[j],ai)
        if it>=burn and (it-burn)%thin==0:
            stats.append((np.mean(mh),np.mean(np.array(mh)>=3),np.mean(np.array(mh)>=4)))
    return np.array(stats)
cons={'TRSL_dep_noncanonical':lambda d,a:True,
      'TRSL_indep_distant':lambda d,a:(a-d-1)>=1000,
      'TRSL_indep_local':lambda d,a:0<(a-d-1)<1000}
NMH=int(sys.argv[1]) if len(sys.argv)>1 else 2000
for g in ['TRSL_dep_noncanonical','TRSL_indep_distant','TRSL_indep_local']:
    d=groups[g]
    st=swap_chain(d.donor.values,d.acceptor.values,cons[g],NMH)
    obs=(d.mh_len.mean(),(d.mh_len>=3).mean(),(d.mh_len>=4).mean())
    for k,nm in enumerate(['mean_mh','frac_mh>=3','frac_mh>=4']):
        res.append(dict(metric='microhomology '+nm,group=g,n=len(d),observed_mean=obs[k],null_pairshuffle_mean=st[:,k].mean(),
           null_ps_p2_5=np.percentile(st[:,k],2.5),null_ps_p97_5=np.percentile(st[:,k],97.5),P_pairshuffle=two_sided(obs[k],st[:,k]),n_perm=len(st)))
# canonical MH: position null. donor kept, acceptor from (i) dep-noncanonical acceptor pool (ii) uniform body positions
pool=np.unique(groups['TRSL_dep_noncanonical'].acceptor.values)
for g in ['canonical_primary5','canonical_all118']:
    d=groups[g]
    for nm,src in [('acceptor pool of TRS-L-dep noncanonical',pool),('uniform genome positions 1000-27900',np.arange(1000,27900))]:
        stat=[]
        for _ in range(NR):
            a=rng.choice(src,size=len(d)); stat.append(np.mean([Rlen(dd,aa) for dd,aa in zip(d.donor.values,a)]))
        stat=np.array(stat); obs=d.mh_len.mean()
        res.append(dict(metric='microhomology mean_mh (position null: %s)'%nm,group=g,n=len(d),observed_mean=obs,null_pos_mean=stat.mean(),
             null_pos_p2_5=np.percentile(stat,2.5),null_pos_p97_5=np.percentile(stat,97.5),P_position=two_sided(obs,stat)))
R=pd.DataFrame(res); R.to_csv(W+'separability_sim_mh.tsv',sep='\t',index=False)
pd.set_option('display.width',250); pd.set_option('display.max_columns',30)
print(R.round(4).to_string())
