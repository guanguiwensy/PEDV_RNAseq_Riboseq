"""G5 tests: observed per-family disruption (homology mapping, Sup L) vs synonymous-substitution null (g5_null.py), exact Poisson-binomial tails, BH,
full genome set (1,083) and cluster representatives (579). Targets = >40 aa internal alt-frame +1/+2 families (n=71); controls = 8-40 aa internal alt-frame families."""
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
import sys
SUF=sys.argv[1] if len(sys.argv)>1 else ''
from synlib import *
from scipy import stats
f=families()
nul=pd.read_csv(OUT+'work/g5_null'+SUF+'_pretain_by_family_kbin.csv')
pr={(r.family_id,r.kbin):r.p_retain for r in nul.itertuples()}
m=pd.read_csv(SL+'PEDV_family_homology_mapping.csv',encoding='utf-8-sig',usecols=['accession','family_id','hom_status_gt40','status_ge8'])
gi={a:i for i,a in enumerate(ids)}; m['gi']=m.accession.map(gi)
fam=f[f.family_id.isin(nul.family_id.unique())].set_index('family_id')
role=nul.drop_duplicates('family_id').set_index('family_id').role
def pb_tails(p,d):
    pmf=np.zeros(len(p)+1); pmf[0]=1
    for x in p:
        new=pmf*(1-x); new[1:]+=pmf[:-1]*x; pmf=new
    c=np.cumsum(pmf); lo=min(1,c[d]); hi=min(1,1-(c[d-1] if d>0 else 0)); return lo,hi,pmf
def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); q=np.empty(n); r=p[o]*n/(np.arange(n)+1); r=np.minimum.accumulate(r[::-1])[::-1]; q[o]=np.minimum(r,1); return q
rows=[]
for scope in ('full_1083','cluster_reps_579'):
    keep=np.ones(len(ids),bool) if scope=='full_1083' else rep
    for crit in ('primary','secondary_ge8'):
        for fid,F in fam.iterrows():
            r=role[fid]
            if crit=='secondary_ge8' and r!='target_gt40': continue
            col='hom_status_gt40' if (r=='target_gt40' and crit=='primary') else 'status_ge8'
            d=m[m.family_id==fid]; st=pd.Series(d[col].values,index=d.gi.values)
            st=st[st.isin(['retained','lost'])]; st=st[keep[st.index.values]]
            par=F.fixed_parent_assignment
            kb=kbin(K[par][st.index.values])
            pd_=np.array([1-pr[(fid,int(k))] for k in kb])
            dobs=int((st=='lost').sum()); n=len(st)
            lo,hi,_=pb_tails(pd_,dobs)
            rows.append(dict(scope=scope,criterion=crit,family_id=fid,role=r,region=par,rel_frame=f'+{int(F.relative_frame_to_parent)}',max_aa=F.max_peptide_aa,n_starts=F.n_start_variants,
              n_assessable=n,n_disrupted=dobs,obs_rate=dobs/max(n,1),null_expected_disrupted=pd_.sum(),null_rate=pd_.mean() if n else np.nan,excess_rate=dobs/max(n,1)-(pd_.mean() if n else np.nan),
              p_lower=lo,p_upper=hi,p_two=min(1,2*min(lo,hi))))
T=pd.DataFrame(rows)
T['null_degenerate']=T.null_expected_disrupted<1   # null predicts <1 disrupted genome: upper-tail p is uninformative
for (sc,cr,ro),idx in T.groupby(['scope','criterion','role']).groups.items():
    T.loc[idx,'q_lower']=bh(T.loc[idx,'p_lower']); T.loc[idx,'q_upper']=bh(T.loc[idx,'p_upper'])
T.to_csv(OUT+'source_data/G5_per_family_disruption_vs_null'+SUF+'.csv',index=False,encoding='utf-8-sig')
# summaries
def summ(d):
    n=len(d)
    return dict(n_families=n,n_obs_rate_below_null=int((d.obs_rate<d.null_rate).sum()),n_obs_rate_equal=int((d.obs_rate==d.null_rate).sum()),n_obs_rate_above_null=int((d.obs_rate>d.null_rate).sum()),
      n_p_lower_lt_005=int((d.p_lower<.05).sum()),n_BH_lower_q_lt_005=int((d.q_lower<.05).sum()),n_p_upper_lt_005=int((d.p_upper<.05).sum()),n_BH_upper_q_lt_005=int((d.q_upper<.05).sum()),
      n_null_degenerate=int(d.null_degenerate.sum()),n_BH_upper_nondegenerate=int(((d.q_upper<.05)&~d.null_degenerate).sum()),median_obs_rate=d.obs_rate.median(),median_null_rate=d.null_rate.median(),median_excess=d.excess_rate.median(),n_zero_disruption=int((d.n_disrupted==0).sum()))
S=[]
for (sc,cr,ro),d in T.groupby(['scope','criterion','role']):
    S.append(dict(scope=sc,criterion=cr,role=ro,region='ALL',rel_frame='ALL',**summ(d)))
    for rg,dd in d.groupby('region'): S.append(dict(scope=sc,criterion=cr,role=ro,region=rg,rel_frame='ALL',**summ(dd)))
    for rf,dd in d.groupby('rel_frame'): S.append(dict(scope=sc,criterion=cr,role=ro,region='ALL',rel_frame=rf,**summ(dd)))
S=pd.DataFrame(S); S.to_csv(OUT+'source_data/G5_summary_by_scope_role_region'+SUF+'.csv',index=False,encoding='utf-8-sig')
# targets vs controls (excess rate), pooled and region-matched; length-stratified controls
C=[]
for sc in ('full_1083','cluster_reps_579'):
    d=T[(T.scope==sc)&(T.criterion=='primary')]
    t=d[d.role=='target_gt40']; c=d[d.role=='control_8to40']
    C.append(dict(scope=sc,stratum='all regions',n_target=len(t),n_control=len(c),median_excess_target=t.excess_rate.median(),median_excess_control=c.excess_rate.median(),mwu_p_two_sided=stats.mannwhitneyu(t.excess_rate,c.excess_rate).pvalue))
    for rg in ORDER:
        tt=t[t.region==rg]; cc=c[c.region==rg]
        if len(tt)>=1 and len(cc)>=1: C.append(dict(scope=sc,stratum=rg,n_target=len(tt),n_control=len(cc),median_excess_target=tt.excess_rate.median(),median_excess_control=cc.excess_rate.median(),mwu_p_two_sided=stats.mannwhitneyu(tt.excess_rate,cc.excess_rate).pvalue if len(tt)>1 and len(cc)>1 else np.nan))
    for lo,hi in ((8,15),(16,25),(26,40)):
        cc=c[(c.max_aa>=lo)&(c.max_aa<=hi)]
        C.append(dict(scope=sc,stratum=f'control {lo}-{hi} aa',n_target=len(t),n_control=len(cc),median_excess_target=t.excess_rate.median(),median_excess_control=cc.excess_rate.median(),mwu_p_two_sided=stats.mannwhitneyu(t.excess_rate,cc.excess_rate).pvalue))
Cdf=pd.DataFrame(C); Cdf.to_csv(OUT+'source_data/G5_targets_vs_controls'+SUF+'.csv',index=False,encoding='utf-8-sig')
pd.set_option('display.width',260); pd.set_option('display.max_columns',40)
print(S[(S.region=='ALL')&(S.rel_frame=='ALL')].round(4).to_string()); print(Cdf.round(4).to_string())
