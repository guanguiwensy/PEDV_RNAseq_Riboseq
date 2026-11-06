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
RES=W+'res/'
m=load_all()
def stat3(d,a):
    x=mh_vec(d,a); return np.array([x.mean(),(x>=3).mean(),(x>=4).mean()])
def analyse(d,cons,nsamp,label,extra=None,seed=SEED):
    don=d.donor.values; acc=d.acceptor.values; n=len(d)
    obs_x=d.mh_len.values
    st,cur=chain(don,acc,cons,stat3,nsamp,seed=seed)
    o=np.array([obs_x.mean(),(obs_x>=3).mean(),(obs_x>=4).mean()])
    bo=boot_mean(obs_x,1000); bo3=boot_mean((obs_x>=3).astype(float),1000)
    dl=bo-st[np.random.default_rng(seed).integers(0,len(st),len(bo)),0]
    r=dict(label=label,n=n,obs_mean_mh=o[0],null_mean_mh=st[:,0].mean(),null_p2_5=np.percentile(st[:,0],2.5),null_p97_5=np.percentile(st[:,0],97.5),
           delta_nt=o[0]-st[:,0].mean(),delta_ci_lo=np.percentile(dl,2.5),delta_ci_hi=np.percentile(dl,97.5),fold_mean=o[0]/st[:,0].mean(),P_mean=p2(o[0],st[:,0]),
           obs_frac_ge3=o[1],null_frac_ge3=st[:,1].mean(),fold_ge3=o[1]/st[:,1].mean(),P_ge3=p2(o[1],st[:,1]),
           obs_frac_ge4=o[2],null_frac_ge4=st[:,2].mean(),fold_ge4=o[2]/max(st[:,2].mean(),1e-9),P_ge4=p2(o[2],st[:,2]),n_perm=len(st))
    if extra: r.update(extra)
    return r,cur
if __name__=='__main__':
    main=[]; e1=[]
    for g in CLS3:
        d=m[(m.cls==g)&m.pass_filter]
        r,cur=analyse(d,CONS[g],NPERM,g); r['group']=g; r['subset']='filtered (main)'; main.append(r)
        xs=mh_vec(d.donor.values,cur)
        for v in d.mh_len.values: e1.append(dict(group=g,kind='observed',mh_len=v))
        for v in xs: e1.append(dict(group=g,kind='pair-shuffled (one draw)',mh_len=v))
        print(g,r['obs_mean_mh'],r['null_mean_mh'],flush=True)
        if g=='TRSL_indep_distant':
            d2=d[~d.near_leader_donor_le100]; r2,_=analyse(d2,CONS[g],NPERM,g+' excl donor<=100'); r2['group']=g; r2['subset']='excl. donor<=100'; main.append(r2)
            print('distant excl',r2['obs_mean_mh'],r2['null_mean_mh'])
        if g=='TRSL_indep_local':
            d2=d[~d.hotspot_27620_27682]; r2,_=analyse(d2,CONS[g],NPERM,g+' excl hotspot'); r2['group']=g; r2['subset']='excl. 27620->27682'; main.append(r2)
    M=pd.DataFrame(main)
    base=M[M.subset=='filtered (main)']
    for c in ['P_mean','P_ge3','P_ge4']:
        M.loc[base.index,c+'_holm_across_3_classes']=holm(base[c].values)
    M.to_csv(RES+'mh_main.tsv',sep='\t',index=False)
    pd.DataFrame(e1).to_csv(RES+'mh_E1_points.tsv',sep='\t',index=False)
    # canonical: position nulls
    rng=np.random.default_rng(SEED)
    dep_pool=np.unique(m[(m.cls=='TRSL_dep_noncanonical')&m.pass_filter].acceptor.values); uni=np.arange(1000,27900)
    cn=[]
    for gname,d in [('canonical_primary5',m[m.canonical_primary]),('canonical_all118',m[m.cls=='canonical'])]:
        obs=d.mh_len.values
        for nm,src in [('acceptor pool of TRS-L-dep noncanonical',dep_pool),('uniform genome positions 1000-27900',uni)]:
            dd=np.tile(d.donor.values,NPERM); aa=rng.choice(src,size=len(dd))
            x=mh_vec(dd,aa).reshape(NPERM,len(d)).mean(1)
            cn.append(dict(group=gname,n=len(d),background=nm,obs_mean_mh=obs.mean(),null_mean=x.mean(),null_p2_5=np.percentile(x,2.5),null_p97_5=np.percentile(x,97.5),P=p2(obs.mean(),x)))
    pd.DataFrame(cn).to_csv(RES+'mh_canonical_position_null.tsv',sep='\t',index=False)
    print(pd.DataFrame(cn))
    # ---------- strata by positive samples (all events of class, stratum-specific background)
    bins=[('1',lambda s:s==1),('2',lambda s:s==2),('3-4',lambda s:(s>=3)&(s<=4)),('>=5',lambda s:s>=5)]
    rows=[]; cells=[]
    for g in CLS3:
        dall=m[m.cls==g].copy()
        dall['sbin']=np.select([b[1](dall.n_samples) for b in bins],[b[0] for b in bins],default='?')
        # reads quartiles within class (duplicates collapsed)
        dall['rbin']=np.select([dall.reads==1,dall.reads==2,dall.reads<10],['1','2','3-9'],default='>=10')
        print(g,dall.rbin.value_counts().to_dict(),flush=True)
        for bname,_ in bins:
            d=dall[dall.sbin==bname]
            if len(d)<10: rows.append(dict(group=g,samples_bin=bname,n=len(d))); continue
            r,_=analyse(d,CONS[g],3000,g+bname); r.update(group=g,samples_bin=bname); rows.append(r)
            for rb,dc in d.groupby('rbin'):
                cell=dict(group=g,samples_bin=bname,reads_bin=rb,n=len(dc),obs_mean_mh=dc.mh_len.mean())
                if len(dc)>=10:
                    r2,_=analyse(dc,CONS[g],1500,g,seed=SEED+7)
                    cell.update(null_mean_mh=r2['null_mean_mh'],delta_nt=r2['delta_nt'],delta_ci_lo=r2['delta_ci_lo'],delta_ci_hi=r2['delta_ci_hi'],P_mean=r2['P_mean'])
                cells.append(cell)
        print('strata done',g,flush=True)
    pd.DataFrame(rows).to_csv(RES+'mh_strata_samples.tsv',sep='\t',index=False)
    pd.DataFrame(cells).to_csv(RES+'mh_2d_cells.tsv',sep='\t',index=False)
    pd.set_option('display.width',250); pd.set_option('display.max_columns',40)
    print(M[['label','n','obs_mean_mh','null_mean_mh','delta_nt','delta_ci_lo','delta_ci_hi','fold_mean','fold_ge3','P_mean','P_mean_holm_across_3_classes']].round(4).to_string())
    print(pd.DataFrame(rows)[['group','samples_bin','n','obs_mean_mh','null_mean_mh','delta_nt','delta_ci_lo','delta_ci_hi','P_mean']].round(3).to_string())
    print(pd.DataFrame(cells).round(3).to_string())
