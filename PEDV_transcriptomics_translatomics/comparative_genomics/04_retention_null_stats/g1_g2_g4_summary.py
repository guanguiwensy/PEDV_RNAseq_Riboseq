"""Summaries of the composition null (10,000 replicates, seed 20261019): G1 densities, G2/G3 family counts, G4 region capacity."""
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
from r6lib import *
from scipy import stats
D=OUT+'work/'; o=np.load(D+'null_observed.npz'); z=np.load(D+'null_replicates.npz')
N=int(z['nrep']); g=list(z['gnames'])
OVL=(12592,12622)
body={}
for n,(a,b) in CDS.items():
    lo,hi=a,b-3
    if n=='ORF1a': hi=OVL[0]-1
    if n=='ORF1b': lo=OVL[1]+1
    body[n]=(lo,hi)
def garwood(k):
    return (0 if k==0 else stats.chi2.ppf(.025,2*k)/2, stats.chi2.ppf(.975,2*k+2)/2)
def pemp(x,v):  # two-sided empirical p with +1 correction
    lo=(1+(x<=v).sum())/(N+1); hi=(1+(x>=v).sum())/(N+1); return min(1,2*min(lo,hi))
rows=[]
for k,label in (('stop','stop (UAA/UAG/UGA)'),('aug','AUG'),('nc','near-cognate x7'),('start8','AUG+7 near-cognate')):
    for gi,gn in enumerate(g):
        reg,f=gn.split('|')
        if k=='start8': ov=o['aug'][0,gi]+o['nc'][0,gi]; x=z['aug'][:,gi]+z['nc'][:,gi]
        else: ov=o[k][0,gi]; x=z[k][:,gi]
        if reg=='ALL': kb=sum(body[n][1]-body[n][0]+1 for n in ORDER)/1000
        else: kb=(body[reg][1]-body[reg][0]+1)/1000
        lo,hi=garwood(int(ov))
        rows.append(dict(codon_class=label,region=reg,rel_frame={'0':'0','1':'+1','2':'+2'}[f],region_body_kb=round(kb,3),observed_count=int(ov),observed_per_kb=round(ov/kb,2),
            obs_poisson95_lo_per_kb=round(lo/kb,2),obs_poisson95_hi_per_kb=round(hi/kb,2),null_mean_per_kb=round(x.mean()/kb,2),null_p2_5_per_kb=round(np.percentile(x,2.5)/kb,2),null_p97_5_per_kb=round(np.percentile(x,97.5)/kb,2),
            obs_over_null=round(ov/x.mean(),3) if x.mean()>0 else np.nan,empirical_p_two_sided=round(pemp(x,ov),4),n_replicates=N,seed=int(z['seed'])))
G1=pd.DataFrame(rows); G1.to_csv(OUT+'source_data/G1_codon_density_by_frame_region.csv',index=False,encoding='utf-8-sig')
# ratio +2/+1 stops
rr=[]
for reg in ORDER+['ALL']:
    i1,i2=g.index(f'{reg}|1'),g.index(f'{reg}|2')
    ob=(o['stop'][0,i2],o['stop'][0,i1]); x1,x2=z['stop'][:,i1],z['stop'][:,i2]
    r=(x2+0.5)/(x1+0.5)
    rr.append(dict(region=reg,stop_plus1=int(ob[1]),stop_plus2=int(ob[0]),obs_ratio_p2_over_p1=round(ob[0]/max(ob[1],1),2),null_ratio_mean=round(r.mean(),2),null_ratio_p2_5=round(np.percentile(r,2.5),2),null_ratio_p97_5=round(np.percentile(r,97.5),2),
        obs_diff=int(ob[0]-ob[1]),null_diff_mean=round((x2-x1).mean(),1),emp_p_diff_ge_obs=round((1+((x2-x1)>=ob[0]-ob[1]).sum())/(N+1),4),emp_p_diff_le_obs=round((1+((x2-x1)<=ob[0]-ob[1]).sum())/(N+1),4)))
pd.DataFrame(rr).to_csv(OUT+'source_data/G1_stop_plus2_vs_plus1_ratio.csv',index=False,encoding='utf-8-sig')
# families by region, rf, tier
TN=['ge8','ge20','gt40']; rows=[]
for ri,rn in enumerate(ORDER+['ALL']):
    for rfi,rf in enumerate(('+1','+2')):
        for ti,t in enumerate(TN):
            for kind,key in (('families','nfam'),('start_isoforms','niso')):
                if rn=='ALL': ov=o[key][0,:,rfi,ti].sum(); x=z[key][:,:,rfi,ti].sum(1)
                else: ov=o[key][0,ri,rfi,ti]; x=z[key][:,ri,rfi,ti]
                rows.append(dict(region=rn,rel_frame=rf,tier=t,unit=kind,observed=int(ov),null_mean=round(x.mean(),2),null_p2_5=int(np.percentile(x,2.5)),null_p97_5=int(np.percentile(x,97.5)),obs_over_null=round(ov/x.mean(),3) if x.mean()>0 else np.nan,empirical_p_two_sided=round(pemp(x,ov),4),n_replicates=N))
G3n=pd.DataFrame(rows); G3n.to_csv(OUT+'source_data/G2G3_family_counts_vs_composition_null.csv',index=False,encoding='utf-8-sig')
# G4: capacity rates vs null; M+N vs rest
LEN={n:b-a+1 for n,(a,b) in CDS.items()}
rows=[]
ti=2; rfi=0
tot_o=o['nfam'][0,:,rfi,ti]; tot_n=z['nfam'][:,:,rfi,ti]
for ri,rn in enumerate(ORDER):
    x=tot_n[:,ri]; ov=tot_o[ri]
    rows.append(dict(region=rn,length_nt=LEN[rn],obs_families=int(ov),obs_per_kb=round(1000*ov/LEN[rn],2),null_mean_families=round(x.mean(),2),null_mean_per_kb=round(1000*x.mean()/LEN[rn],2),null_p2_5=int(np.percentile(x,2.5)),null_p97_5=int(np.percentile(x,97.5)),obs_over_null=round(ov/x.mean(),3),emp_p_two_sided=round(pemp(x,ov),4),
      obs_iso=int(o['niso'][0,ri,rfi,ti]),null_mean_iso=round(z['niso'][:,ri,rfi,ti].mean(),1),emp_p_iso_two_sided=round(pemp(z['niso'][:,ri,rfi,ti],o['niso'][0,ri,rfi,ti]),4)))
pd.DataFrame(rows).to_csv(OUT+'source_data/G4_capacity_vs_composition_null.csv',index=False,encoding='utf-8-sig')
# M+N vs rest (observed rates, uniform and composition null)
iM,iN=ORDER.index('M'),ORDER.index('N'); mn=[iM,iN]; rest=[i for i in range(7) if i not in mn]
kMN=tot_o[mn].sum(); LMN=LEN['M']+LEN['N']; kR=tot_o[rest].sum(); LR=sum(LEN[ORDER[i]] for i in rest)
res={'M+N families':int(kMN),'M+N length':LMN,'M+N per kb':1000*kMN/LMN,'rest families':int(kR),'rest length':LR,'rest per kb':1000*kR/LR,'rate ratio':(kMN/LMN)/(kR/LR),
     'binom p one-sided (>=)':stats.binomtest(int(kMN),int(kMN+kR),LMN/(LMN+LR),alternative='greater').pvalue,'binom p two-sided':stats.binomtest(int(kMN),int(kMN+kR),LMN/(LMN+LR)).pvalue}
# composition null for the ratio of per-kb rates M+N vs rest
nm=tot_n[:,mn].sum(1)/LMN; nr=tot_n[:,rest].sum(1)/LR; ratio=(nm+1e-9)/(nr+1e-9)
res['null ratio mean']=ratio.mean(); res['null ratio p2.5']=np.percentile(ratio,2.5); res['null ratio p97.5']=np.percentile(ratio,97.5)
res['emp p (null ratio >= obs)']=(1+(ratio>=res['rate ratio']).sum())/(N+1)
# ranks of N and M among the 7 by per-kb rate (observed), and under null rate-normalized
per=1000*tot_o/np.array([LEN[n] for n in ORDER]); res['rank of N (1=highest) among 7 incl E,ORF3']=int((per>per[iN]).sum()+1); res['rank of M']=int((per>per[iM]).sum()+1)
json.dump({k:(float(v) if not isinstance(v,str) else v) for k,v in res.items()},open(OUT+'source_data/G4_M_N_vs_rest.json','w'),indent=1,ensure_ascii=False)
pd.set_option('display.width',250); pd.set_option('display.max_columns',30)
print(pd.DataFrame(rr).to_string()); print(pd.DataFrame(rows).to_string()); print(json.dumps(res,indent=1))
print(G3n[(G3n.region=='ALL')].to_string())
print(G1[(G1.region=='ALL')&(G1.codon_class!='near-cognate x7')].to_string())
