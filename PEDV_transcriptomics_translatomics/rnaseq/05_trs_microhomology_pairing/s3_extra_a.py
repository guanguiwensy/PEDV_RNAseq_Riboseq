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
from s3_lib import *
import random, itertools
from scipy.stats import spearmanr, rankdata
import RNA
rng=np.random.default_rng(SEED); random.seed(SEED)
m=load_all()
NR=NPERM
RES=W+'res/'
REF='CAACTAAACGAA'
enc={c:i for i,c in enumerate('ACGT')}; RE=np.array([enc[c] for c in REF],dtype=np.int8)
KS=np.arange(-20,21)
def mm_vec(arr):
    n=len(arr)-11
    return sum((arr[i:i+n]==RE[i]).astype(np.int16) for i in range(12))
MM=mm_vec(S0)
def mid(a):
    st=(a-1)+KS-3; ok=(st>=0)&(st<len(MM)); v=np.full(41,np.nan); v[ok]=MM[st[ok]]; return v
ALLV=np.array([mid(a) for a in range(1,GL+1)]); MAXG=np.nanmax(ALLV,axis=1)
G={'canonical':m[m.canonical_primary],'TRSL_dep_noncanonical':m[(m.cls=='TRSL_dep_noncanonical')&m.pass_filter],
   'TRSL_indep_distant':m[(m.cls=='TRSL_indep_distant')&m.pass_filter],'TRSL_indep_local':m[(m.cls=='TRSL_indep_local')&m.pass_filter]}
null=np.load(W+'sim_null_dinuc.npy')*12; ids=pd.read_csv(W+'sim_null_dinuc_ids.tsv',sep='\t').junction_id.tolist(); idx={j:i for i,j in enumerate(ids)}
CDS={'ORF1a':(293,12622),'ORF1b':(12592,20613),'S':(20610,24758),'ORF3':(24758,25033),'E':(25361,25600),'M':(25608,26288),'N':(26300,27625)}
def region(pos):
    h=[(s,g) for g,(s,e) in CDS.items() if s<=pos<=e]
    return max(h)[1] if h else ("5'UTR" if pos<293 else ("3'UTR" if pos>27625 else 'intergenic'))
regpos={}
for p in range(1,GL+1): regpos.setdefault(region(p),[]).append(p)
regpos={k:np.array(v) for k,v in regpos.items()}
def dinuc_shuffle(s):
    last=s[-1]; nxt={}
    for a,b in zip(s,s[1:]): nxt.setdefault(a,[]).append(b)
    verts=[v for v in nxt if v!=last]
    while True:
        le={v:random.choice(nxt[v]) for v in verts}; ok=True
        for v in verts:
            x=v; c=0
            while x!=last and c<10: x=le.get(x,last); c+=1
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
effect=[]; curve_obs={}; curve_null={}
NSH_CURVE=50
for g,d in G.items():
    n=len(d); obs=d.maxid.values*12
    ii=[idx[j] for j in d.junction_id]; nl=null[ii]
    draws=nl[np.arange(n)[:,None],rng.integers(0,nl.shape[1],(n,NR))].mean(0)
    rp=np.array([rng.choice(regpos[r],size=NR) for r in d.acceptor_region])
    dB=MAXG[rp-1].mean(0)
    bo=boot_mean(obs,2000)
    for nm,dr in [('composition-preserving (dinucleotide-shuffled context)',draws),('region-matched random genomic positions',dB)]:
        dl=bo-dr[rng.integers(0,len(dr),len(bo))]
        effect.append(dict(group=g,background=nm,n=n,obs_mean_matches_of_12=obs.mean(),bg_mean_matches=dr.mean(),bg_p2_5=np.percentile(dr,2.5),bg_p97_5=np.percentile(dr,97.5),
            delta_matches=obs.mean()-dr.mean(),delta_ci_lo=np.percentile(dl,2.5),delta_ci_hi=np.percentile(dl,97.5),P_two_sided=p2(obs.mean(),dr)))
    ctxs=[seq[a-24:a+28] for a in d.acceptor.values]
    cur=np.zeros((n,41))
    for i,c in enumerate(ctxs):
        acc=np.zeros(41)
        for _ in range(NSH_CURVE):
            arr=np.array([enc.get(x,4) for x in dinuc_shuffle(c)],dtype=np.int8); acc+=mm_vec(arr)[20+KS]
        cur[i]=acc/NSH_CURVE
    curve_null[g]=cur
    curve_obs[g]=np.array([mid(a) for a in d.acceptor.values])
    print('sim done',g,flush=True)
E=pd.DataFrame(effect)
for bg,sub in E.groupby('background'):
    E.loc[sub.index,'P_holm_across_4_groups']=holm(sub.P_two_sided.values)
E['P_floor']=1/(NR+1)
E.to_csv(RES+'sim_effect_matches.tsv',sep='\t',index=False)
cr=[]
for g in G:
    for kind,arr in [('observed',curve_obs[g]),('composition-shuffled background',curve_null[g])]:
        n=len(arr); r=np.random.default_rng(SEED)
        bs=np.array([np.nanmean(arr[r.integers(0,n,n)],axis=0) for _ in range(1000)])
        for j,k in enumerate(KS):
            cr.append(dict(group=g,kind=kind,offset=k,mean_matches=np.nanmean(arr[:,j]),ci_lo=np.nanpercentile(bs[:,j],2.5),ci_hi=np.nanpercentile(bs[:,j],97.5),n=n))
pd.DataFrame(cr).to_csv(RES+'sim_curves_D1.tsv',sep='\t',index=False)
D2=[]; r=np.random.default_rng(SEED+1)
for g,d in G.items():
    ii=[idx[j] for j in d.junction_id]
    bgv=null[ii,r.integers(0,null.shape[1],len(ii))]
    for (jid,co,mx,kb,bgx) in zip(d.junction_id,d.canonical_gene,d.maxid*12,d.best_k,bgv):
        D2.append(dict(group=g,junction_id=jid,canonical_gene=co,best_shift_matches=mx,best_k=kb,background_one_shuffle=bgx))
pd.DataFrame(D2).to_csv(RES+'sim_D2_points.tsv',sep='\t',index=False)
# ===== Panel A
core_start={'S':20602,'E':25351,'M':25597,'N':26285,'ORF3':None}
core_mm_tab={'S':1,'E':1,'M':1,'N':0,'ORF3':np.nan}
canon=m[m.cls=='canonical']; reads_g=canon.groupby('canonical_gene').reads.sum()
prim=m[m.canonical_primary].set_index('canonical_gene')
A=[]
for g in GENES:
    a0=core_start[g] if core_start[g] else CANON_SITES[g]
    body=seq[a0-1-3:a0-1+9]; assert len(body)==12
    match=[x==y for x,y in zip(body,REF)]
    run=best=0
    for t in match:
        run=run+1 if t else 0; best=max(best,run)
    A.append(dict(gene=g,alignment_anchor=a0,anchor_type='TRS-B core start' if core_start[g] else 'observed breakpoint (TRS-B core undetermined)',
        body_12mer=body,ref_12mer=REF,matched_of_12=sum(match),longest_contiguous=best,identity=sum(match)/12,core_hexamer=seq[a0-1:a0+5],
        core_mismatch_vs_CTAAAC=core_mm_tab[g],primary_donor=int(prim.loc[g,'donor']),primary_acceptor=int(prim.loc[g,'acceptor']),
        mfe_a_kcal=prim.loc[g,'mfe_a'],mfe_b_kcal=prim.loc[g,'mfe_b'],canonical_reads=int(reads_g[g]),canonical_share=reads_g[g]/reads_g.sum(),
        maxid_matches=prim.loc[g,'maxid']*12,best_k=int(prim.loc[g,'best_k']),longest_contig_at_best_shift=int(prim.loc[g,'longest_contig_at_best'])))
A=pd.DataFrame(A); A.to_csv(RES+'panelA_table.tsv',sep='\t',index=False)
y=A.canonical_reads.values.astype(float)
def exact_spearman(x,y):
    rx=rankdata(x); ry=rankdata(y); obs=np.corrcoef(rx,ry)[0,1]; c=0; t=0
    for p in itertools.permutations(range(5)):
        v=np.corrcoef(rx,ry[list(p)])[0,1]; t+=1; c+= abs(v)>=abs(obs)-1e-12
    return obs,c/t
Bt=[]
for nm,col in [('best-shift identity (matches of 12)','maxid_matches'),('longest contiguous match (at best shift)','longest_contig_at_best_shift'),
               ('MFE (a) nascent-strand re-pairing','mfe_a_kcal'),('MFE (b) direct positive-positive','mfe_b_kcal'),
               ('matched bases at aligned TRS-B core (Panel A)','matched_of_12')]:
    r_,p_=exact_spearman(A[col].values,y); Bt.append(dict(metric=nm,column=col,rho=r_,P_exact_two_sided=p_,n=5,P_min_possible=2/120))
pd.DataFrame(Bt).to_csv(RES+'panelB_spearman.tsv',sep='\t',index=False)
SU=seq.replace('T','U'); COMP=str.maketrans('ACGU','UGCA'); rc=lambda s:s.translate(COMP)[::-1]
dw=lambda d:SU[d-15:d+15]; aw=lambda a:SU[a-16:a+14]
sens=[]
for g in GENES:
    a=CANON_SITES[g]; sens.append(dict(gene=g,mfe_a_donor66=RNA.duplexfold(rc(dw(66)),aw(a)).energy,mfe_b_donor66=RNA.duplexfold(dw(66),aw(a)).energy))
sens=pd.DataFrame(sens); sens['reads']=y
rr=exact_spearman(sens.mfe_a_donor66.values,y); rb=exact_spearman(sens.mfe_b_donor66.values,y)
sens['rho_a_fixed_donor66']=rr[0]; sens['P_a']=rr[1]; sens['rho_b_fixed_donor66']=rb[0]; sens['P_b']=rb[1]
sens.to_csv(RES+'panelB_sens_fixed_donor66.tsv',sep='\t',index=False)
def freq(accs):
    M=np.zeros((41,4)); n=0
    for a in accs:
        s=seq[a-26:a+15]
        if len(s)!=41: continue
        n+=1
        for i,c in enumerate(s): M[i,'ACGT'.find(c)]+=1
    return M/n,n
CF=[]
regs=m[m.pass_filter|m.canonical_primary].acceptor_region.values
bg_acc=np.concatenate([rng.choice(regpos[r_],size=20) for r_ in regs])
for g,acc in [('canonical',G['canonical'].acceptor.values),('TRSL_dep_noncanonical',G['TRSL_dep_noncanonical'].acceptor.values),
              ('TRSL_indep_distant',G['TRSL_indep_distant'].acceptor.values),('TRSL_indep_local',G['TRSL_indep_local'].acceptor.values),
              ('background_region_matched_positions',bg_acc)]:
    F,n=freq(acc)
    for i,pos in enumerate(range(-25,16)):
        for j,nt in enumerate('ACGT'): CF.append(dict(group=g,n=n,position=pos,nucleotide=nt.replace('T','U'),frequency=F[i,j]))
pd.DataFrame(CF).to_csv(RES+'panelC_frequencies.tsv',sep='\t',index=False)
pd.set_option('display.width',250); pd.set_option('display.max_columns',30)
print(E.round(4).to_string()); print(A.round(3).to_string()); print(pd.DataFrame(Bt).round(4)); print(sens.round(3).to_string())
