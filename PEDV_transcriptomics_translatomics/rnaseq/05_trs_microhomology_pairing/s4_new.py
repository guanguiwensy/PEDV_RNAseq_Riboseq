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
from metrics_mfe import matrix
from multiprocessing import Pool
import itertools
from scipy.stats import rankdata
RES=W+'res/'
def chain_span(don,acc,edges,stat_fn,nsamp,sweeps=2,burn=30,seed=SEED):
    rng=np.random.default_rng(seed); don=np.asarray(don); cur=np.asarray(acc).copy(); n=len(don); h=n//2
    span=cur-don-1; b=np.digitize(span,edges); out=[]; acc_rate=[]
    for it in range(burn+nsamp*sweeps):
        p=rng.permutation(n); i=p[:h]; j=p[h:2*h]
        ai=cur[i]; aj=cur[j]
        ok=(np.digitize(aj-don[i]-1,edges)==b[i])&(np.digitize(ai-don[j]-1,edges)==b[j])&(ai!=aj)
        cur[i[ok]]=aj[ok]; cur[j[ok]]=ai[ok]; acc_rate.append(ok.mean())
        if it>=burn and (it-burn)%sweeps==sweeps-1: out.append(stat_fn(don,cur))
    return np.array(out),cur,float(np.mean(acc_rate))
def edges_of(d):
    sp=(d.acceptor-d.donor-1).values; return np.percentile(sp,[25,50,75])
def stat3(d,a):
    x=mh_vec(d,a); return np.array([x.mean(),(x>=3).mean(),(x>=4).mean()])
if __name__=='__main__':
    m=load_all()
    old=pd.read_csv(RES+'mh_main.tsv',sep='\t')
    rows=[]
    sets=[]
    for g in CLS3:
        d=m[(m.cls==g)&m.pass_filter]; sets.append((g,'filtered (main)',d,d))
    # sensitivity sets use the quartile edges of the parent set
    dd=m[(m.cls=='TRSL_indep_distant')&m.pass_filter]; sets.append(('TRSL_indep_distant','excl. donor<=100',dd[~dd.near_leader_donor_le100],dd))
    ll=m[(m.cls=='TRSL_indep_local')&m.pass_filter]; sets.append(('TRSL_indep_local','excl. 27620->27682',ll[~ll.hotspot_27620_27682],ll))
    for g,lab,d,parent in sets:
        ed=edges_of(parent)
        st,cur,ar=chain_span(d.donor.values,d.acceptor.values,ed,stat3,NPERM)
        obs=d.mh_len.values; o=np.array([obs.mean(),(obs>=3).mean(),(obs>=4).mean()])
        bo=boot_mean(obs,1000); dl=bo-st[np.random.default_rng(SEED).integers(0,len(st),len(bo)),0]
        r=dict(group=g,subset=lab,n=len(d),span_quartile_edges='/'.join('%d'%e for e in ed),swap_acceptance=ar,obs_mean_mh=o[0],null_mean_mh=st[:,0].mean(),null_p2_5=np.percentile(st[:,0],2.5),null_p97_5=np.percentile(st[:,0],97.5),
               delta_nt=o[0]-st[:,0].mean(),delta_ci_lo=np.percentile(dl,2.5),delta_ci_hi=np.percentile(dl,97.5),fold_mean=o[0]/st[:,0].mean(),P_mean=p2(o[0],st[:,0]),
               obs_frac_ge3=o[1],null_frac_ge3=st[:,1].mean(),fold_ge3=o[1]/st[:,1].mean(),P_ge3=p2(o[1],st[:,1]))
        ol=old[(old.group==g)&(old.subset==('filtered (main)' if lab=='filtered (main)' else lab))]
        if len(ol): r['unstratified_null_mean_mh']=ol.null_mean_mh.iloc[0]; r['unstratified_delta_nt']=ol.delta_nt.iloc[0]
        rows.append(r); print(g,lab,'obs %.3f null %.3f delta %.3f acc %.3f'%(o[0],st[:,0].mean(),o[0]-st[:,0].mean(),ar),flush=True)
        if lab=='filtered (main)':
            xs=mh_vec(d.donor.values,cur)
            pd.DataFrame(dict(group=g,kind='span-stratified pair-shuffled (one draw)',mh_len=xs)).to_csv(RES+'mh_E1_spanshuf_%s.tsv'%g,sep='\t',index=False)
    M=pd.DataFrame(rows); base=M[M.subset=='filtered (main)']
    for c in ['P_mean','P_ge3']: M.loc[base.index,c+'_holm_across_3_classes']=holm(base[c].values)
    M.to_csv(RES+'mh_main_spanstrat.tsv',sep='\t',index=False)
    # ---------------- MFE span-stratified
    out=[]; pts=[]
    with Pool(26) as pool:
        for g in CLS3:
            d=m[(m.cls==g)&m.pass_filter]; dn=d.donor.values; ac=d.acceptor.values; ed=edges_of(d)
            ud=sorted(set(dn)); ua=sorted(set(ac)); di={x:i for i,x in enumerate(ud)}
            dI=np.array([di[x] for x in dn]); lut=np.full(GL+5,-1); lut[np.array(ua)]=np.arange(len(ua))
            for mode in ['b','a']:
                Mx=matrix(ud,ua,mode,pool)
                subs=[('filtered (main)',np.ones(len(d),bool))]
                if g=='TRSL_indep_distant': subs.append(('excl. donor<=100',~d.near_leader_donor_le100.values))
                for lab,sel in subs:
                    fn=lambda don,cur,sel=sel:Mx[dI[sel],lut[cur]].mean()
                    st,cur,ar=chain_span(dn[sel],ac[sel],ed,fn,NPERM)
                    obs=d['mfe_'+mode].values[sel]; bo=boot_mean(obs,1000); dl=bo-st[np.random.default_rng(SEED).integers(0,len(st),len(bo))]
                    out.append(dict(group=g,mode=mode,subset=lab,n=int(sel.sum()),swap_acceptance=ar,obs_mean_mfe=obs.mean(),null_mean=st.mean(),null_p2_5=np.percentile(st,2.5),null_p97_5=np.percentile(st,97.5),
                                    delta_kcal=obs.mean()-st.mean(),delta_ci_lo=np.percentile(dl,2.5),delta_ci_hi=np.percentile(dl,97.5),P=p2(obs.mean(),st)))
                    if lab=='filtered (main)':
                        e_sh=Mx[dI,lut[cur]]
                        for v in e_sh: pts.append(dict(group=g,mode=mode,kind='span-stratified pair-shuffled (one draw)',mfe=v))
                    print(g,mode,lab,out[-1]['obs_mean_mfe'],out[-1]['null_mean'],out[-1]['swap_acceptance'],flush=True)
    F=pd.DataFrame(out)
    for mode in ['a','b']:
        ix=F[(F['mode']==mode)&(F.subset=='filtered (main)')].index; F.loc[ix,'P_holm_across_3_classes']=holm(F.loc[ix,'P'].values)
    F.to_csv(RES+'mfe_main_spanstrat.tsv',sep='\t',index=False); pd.DataFrame(pts).to_csv(RES+'mfe_F1_points_spanstrat.tsv',sep='\t',index=False)
    # ---------------- Panel B: same-source metrics, ties, per-sample reads
    A=pd.read_csv(RES+'panelA_table.tsv',sep='\t')
    ps=pd.read_csv(B+'figures/Fig2_work/per_sample_events.tsv.gz',sep='\t')
    man=pd.read_csv(B+'sample_manifest_GIb.tsv',sep='\t'); s20=man[man.batch=='batch3'].sample_id.tolist()
    ps=ps[ps['sample'].isin(s20)&ps.donor.between(55,71)]
    rows=[]; rng=np.random.default_rng(SEED)
    for g in GENES:
        site=CANON_SITES[g]; x=ps[(ps.acc-site).abs()<=5].groupby('sample').reads.sum().reindex(s20).fillna(0).values
        bs=np.array([x[rng.integers(0,len(x),len(x))].mean() for _ in range(5000)])
        rows.append(dict(gene=g,mean_reads_per_sample=x.mean(),ci_lo=np.percentile(bs,2.5),ci_hi=np.percentile(bs,97.5),total_reads=x.sum(),n_samples=len(x),min_sample=x.min(),max_sample=x.max()))
    PS=pd.DataFrame(rows); PS.to_csv(RES+'panelB_reads_per_sample.tsv',sep='\t',index=False)
    print(PS.round(1)); print(A[['gene','canonical_reads']])
    y=PS.mean_reads_per_sample.values
    def exact_spearman(x,y):
        rx=rankdata(x); ry=rankdata(y); obs=np.corrcoef(rx,ry)[0,1]; c=t=0
        for p in itertools.permutations(range(5)):
            v=np.corrcoef(rx,ry[list(p)])[0,1]; t+=1; c+=abs(v)>=abs(obs)-1e-12
        return obs,c/t
    rows=[]
    for nm,col in [('Matched bases (of 12), aligned at TRS-B core (= Panel A)','matched_of_12'),('Longest contiguous match (nt) (= Panel A)','longest_contiguous'),
                   ('MFE (a) nascent-strand re-pairing','mfe_a_kcal'),('MFE (b) direct plus-strand pairing','mfe_b_kcal'),('[supplement] best-shift matched bases','maxid_matches')]:
        v=A[col].values; ties=int(len(v)-len(set(v)))
        r_,p_=exact_spearman(v,y); rows.append(dict(metric=nm,column=col,values=';'.join(str(x) for x in v),rho=r_,P_exact_two_sided=p_,n=5,n_tied_values=ties,P_min_possible=2/120))
    pd.DataFrame(rows).to_csv(RES+'panelB_spearman_v3.tsv',sep='\t',index=False); print(pd.DataFrame(rows).round(3).to_string())
    # rho possible values / exact P table for rho=0.9 (n=5, no ties)
    perm=[np.corrcoef(np.arange(5),np.array(p))[0,1] for p in itertools.permutations(range(5))]
    pe=pd.Series(np.round(perm,3)).value_counts().sort_index(); print(pe.to_string())
    tot=len(perm); print('P(|rho|>=0.9)=',np.mean(np.abs(perm)>=0.9-1e-9),'one-sided',np.mean(np.array(perm)>=0.9-1e-9))
