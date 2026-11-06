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
from metrics_mfe import matrix
from multiprocessing import Pool
from scipy.stats import spearmanr
RES=W+'res/'
if __name__=='__main__':
    m=load_all()
    out=[]; pts=[]; f2=[]
    with Pool(26) as pool:
        for g in CLS3:
            d=m[(m.cls==g)&m.pass_filter]
            dn=d.donor.values; ac=d.acceptor.values
            ud=sorted(set(dn)); ua=sorted(set(ac)); di={x:i for i,x in enumerate(ud)}; ai={x:i for i,x in enumerate(ua)}
            dI=np.array([di[x] for x in dn]); aI_all=np.array([ai[x] for x in ac])
            for mode in ['b','a']:
                Mx=matrix(ud,ua,mode,pool)
                # chain operates on acceptor VALUES; stat maps values->index
                lut=np.full(GL+5,-1); lut[np.array(ua)]=np.arange(len(ua)); amap=lambda cur,lut=lut:lut[cur]
                def stat_fn(don,cur,dI=dI,Mx=Mx,amap=amap,subset=None):
                    e=Mx[dI,amap(cur)]
                    return e if subset is None else e[subset]
                sub_sets=[('filtered (main)',None)]
                if g=='TRSL_indep_distant': sub_sets.append(('excl. donor<=100',~d.near_leader_donor_le100.values))
                for lab,subset in sub_sets:
                    sel=np.ones(len(d),bool) if subset is None else subset
                    dd=dn[sel]; aa=ac[sel]
                    fn=lambda don,cur,sel=sel:Mx[dI[sel],amap(cur)].mean()
                    st,cur=chain(dd,aa,CONS[g],fn,NPERM)
                    obs=d['mfe_'+mode].values[sel]
                    bo=boot_mean(obs,1000); dl=bo-st[np.random.default_rng(SEED).integers(0,len(st),len(bo))]
                    out.append(dict(group=g,mode=mode,subset=lab,n=int(sel.sum()),obs_mean_mfe=obs.mean(),null_mean=st.mean(),null_p2_5=np.percentile(st,2.5),null_p97_5=np.percentile(st,97.5),
                                    delta_kcal=obs.mean()-st.mean(),delta_ci_lo=np.percentile(dl,2.5),delta_ci_hi=np.percentile(dl,97.5),P=p2(obs.mean(),st),
                                    gc_donor=d.gc_donor.values[sel].mean(),gc_acceptor=d.gc_acceptor.values[sel].mean(),n_perm=len(st)))
                    if subset is None:
                        e_sh=Mx[dI,amap(cur)]
                        for v in obs: pts.append(dict(group=g,mode=mode,kind='observed',mfe=v))
                        for v in e_sh: pts.append(dict(group=g,mode=mode,kind='pair-shuffled (one draw)',mfe=v))
                    print(g,mode,lab,out[-1]['obs_mean_mfe'],out[-1]['null_mean'],flush=True)
    M=pd.DataFrame(out)
    for mode in ['a','b']:
        idx=M[(M['mode']==mode)&(M.subset=='filtered (main)')].index
        M.loc[idx,'P_holm_across_3_classes']=holm(M.loc[idx,'P'].values)
    M.to_csv(RES+'mfe_main.tsv',sep='\t',index=False)
    pd.DataFrame(pts).to_csv(RES+'mfe_F1_points.tsv',sep='\t',index=False)
    # F2: Spearman MFE vs reads (with/without flagged events; canonical primary 5 shown as points)
    for g in CLS3:
        d=m[(m.cls==g)&m.pass_filter]
        for mode in ['b','a']:
            r,p=spearmanr(d['mfe_'+mode],d.reads); f2.append(dict(group=g,mode=mode,subset='filtered (main)',n=len(d),rho=r,P=p))
            if g=='TRSL_indep_distant':
                dd=d[~d.flag_leaderdonor_at_canon_site_in_distant]; r,p=spearmanr(dd['mfe_'+mode],dd.reads)
                f2.append(dict(group=g,mode=mode,subset='excl. %d leader-donor events at canonical-site acceptors'%d.flag_leaderdonor_at_canon_site_in_distant.sum(),n=len(dd),rho=r,P=p))
    c=m[m.canonical_primary]
    for mode in ['b','a']:
        r,p=spearmanr(c['mfe_'+mode],c.reads); f2.append(dict(group='canonical_primary5',mode=mode,subset='primary events',n=5,rho=r,P=p))
    F=pd.DataFrame(f2); F.to_csv(RES+'mfe_F2_spearman.tsv',sep='\t',index=False)
    pd.set_option('display.width',250); pd.set_option('display.max_columns',40)
    print(M.round(4).to_string()); print(F.round(4).to_string())
