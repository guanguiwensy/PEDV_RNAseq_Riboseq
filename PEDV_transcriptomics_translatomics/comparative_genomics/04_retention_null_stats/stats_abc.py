"""Chi-square tests requested for panels A, B, C (+ translation-efficiency correlation, panel E).
A: codon class x frame (hit vs non-hit among triplets in the CDS body, ORF1a/1b overlap excluded; same triplet sets as g1_null.py).
B: family length class (8-20, 21-40, >40 aa) x frame (+1 vs +2); also >40 vs <=40.
C: >40-aa +1/+2 families per kb across ORF1a, ORF1b, S, M, N (ORF3 and E excluded as requested): chi-square goodness of fit against rate proportional to length (and M+N vs the rest).
E: Spearman/Pearson between translation efficiency (Ribo P-site density / RNA relative depth) and families per kb, n=5 ORFs, exact permutation p.
Note: triplets of different frames overlap in sequence (not independent); chi-square p-values are therefore descriptive; the 10,000-draw composition null (G1_codon_density_by_frame_region.csv) is the stronger test."""
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
import itertools
from r6lib import *
from scipy import stats
SD=OUT+'source_data/'
OVL=(12592,12622)
body={}
for n,(a,b) in CDS.items():
    lo,hi=a,b-3
    if n=='ORF1a': hi=OVL[0]-1
    if n=='ORF1b': lo=OVL[1]+1
    body[n]=(lo,hi)
# ---------- A
cod=np.array([code[ref[i:i+3]] if set(ref[i:i+3])<=set('ACGT') else 0 for i in range(LREF-2)])
cls={'stop':isstop[cod],'AUG':isaug[cod],'AUG+7 near-cognate':(isaug|isnc)[cod]}
ntrip={}; hits={k:{} for k in cls}
for fr in (0,1,2):
    n=0; h={k:0 for k in cls}
    for g,(lo,hi) in body.items():
        a=CDS[g][0]
        for p in range(lo,hi-1):
            if (p-a)%3==fr:
                n+=1
                for k in cls: h[k]+=int(cls[k][p-1])
    ntrip[fr]=n
    for k in cls: hits[k][fr]=h[k]
rows=[]
def chi(a,b,na,nb,k,fa,fb):
    tab=np.array([[a,na-a],[b,nb-b]]); c,p,_,e=stats.chi2_contingency(tab,correction=False)
    cy,py,_,_=stats.chi2_contingency(tab,correction=True)
    odds=(a/(na-a))/(b/(nb-b)) if b>0 and a>0 else np.nan
    return dict(codon_class=k,comparison=f'frame {fa} vs frame {fb}',hits_a=a,triplets_a=na,rate_a_per_kb=round(1000*a/(3*na),2),hits_b=b,triplets_b=nb,rate_b_per_kb=round(1000*b/(3*nb),2),chi2=round(c,2),p_chi2=p,p_chi2_yates=py,odds_ratio=round(odds,3) if odds==odds else np.nan,min_expected=round(e.min(),1))
for k in cls:
    for fa,fb in ((1,2),(0,1),(0,2)):
        if k=='stop' and 0 in (fa,fb): continue   # frame 0 contains no stop by construction
        rows.append(chi(hits[k][fa],hits[k][fb],ntrip[fa],ntrip[fb],k,('+%d'%fa if fa else '0'),('+%d'%fb if fb else '0')))
    tab=np.array([[hits[k][f],ntrip[f]-hits[k][f]] for f in (0,1,2)] if k!='stop' else [[hits[k][f],ntrip[f]-hits[k][f]] for f in (1,2)])
    c,p,dof,_=stats.chi2_contingency(tab,correction=False)
    rows.append(dict(codon_class=k,comparison='all frames' if k!='stop' else '+1 and +2 (same as above)',chi2=round(c,2),p_chi2=p,dof=dof))
A=pd.DataFrame(rows); A.to_csv(SD+'R6-A_chi_square.csv',index=False,encoding='utf-8-sig')
# ---------- B
fam=pd.read_csv(SD+'G2_family_lengths_plus1_plus2.csv',dtype={'rel_frame':str}); iso=isoforms(); iso=iso[iso.family_id.isin(fam.family_id)].merge(fam[['family_id','rel_frame']],on='family_id')
def cls3(v): return np.array([((v>=8)&(v<=20)).sum(),((v>20)&(v<=40)).sum(),(v>40).sum()])
rowsB=[]
for unit,df,col in (('family',fam,'max_aa'),('start isoform',iso,'peptide_length_aa')):
    t=np.array([cls3(df[df.rel_frame=='+1'][col].values),cls3(df[df.rel_frame=='+2'][col].values)])
    c,p,dof,e=stats.chi2_contingency(t,correction=False)
    rowsB.append(dict(unit=unit,test='chi-square 2x3 (8-20, 21-40, >40 aa) x (+1, +2)',table=str(t.tolist()),chi2=round(c,2),dof=dof,p=p,min_expected=round(e.min(),1)))
    t2=np.array([[t[0,:2].sum(),t[0,2]],[t[1,:2].sum(),t[1,2]]]); c2,p2,_,e2=stats.chi2_contingency(t2,correction=False)
    rowsB.append(dict(unit=unit,test='chi-square 2x2 (<=40 vs >40 aa) x (+1, +2)',table=str(t2.tolist()),chi2=round(c2,2),dof=1,p=p2,min_expected=round(e2.min(),1),fisher_p=stats.fisher_exact(t2)[1]))
B=pd.DataFrame(rowsB); B.to_csv(SD+'R6-B_chi_square.csv',index=False,encoding='utf-8-sig')
# ---------- C
f=families(); ia=f[(f.primary_category=='internal alternative-frame ORF')&(f.max_peptide_aa>40)&f.relative_frame_to_parent.isin([1,2])]
LEN={n:b-a+1 for n,(a,b) in CDS.items()}; ORFS5=['ORF1a','ORF1b','S','M','N']
rowsC=[]; cnt={}
for label,sub in (('+1 and +2 (as drawn; 71 families)',ia),('+1 only (68 families)',ia[ia.relative_frame_to_parent==1])):
    k=np.array([(sub.fixed_parent_assignment==r).sum() for r in ORFS5]); L=np.array([LEN[r] for r in ORFS5]); cnt[label]=k
    exp=k.sum()*L/L.sum(); c,p=stats.chisquare(k,exp)
    mc=(1+sum(1 for _ in [0]))  # placeholder
    rng=np.random.default_rng(20261023); sim=rng.multinomial(k.sum(),L/L.sum(),size=100000); cs=(((sim-exp)**2)/exp).sum(1); pmc=(1+(cs>=c).sum())/100001
    mn=k[3:].sum(); rest=k[:3].sum(); t=np.array([[mn,L[3:].sum()/1000],[rest,L[:3].sum()/1000]])
    pb=stats.binomtest(int(mn),int(k.sum()),L[3:].sum()/L.sum()).pvalue
    tab=np.array([[mn,k.sum()-mn],[L[3:].sum()/L.sum()*k.sum(),(1-L[3:].sum()/L.sum())*k.sum()]])
    rowsC.append(dict(definition=label,counts=dict(zip(ORFS5,k.tolist())),per_kb=dict(zip(ORFS5,np.round(1000*k/L,2).tolist())),chi2_GOF_vs_length_proportional=round(c,2),dof=4,p_chi2=p,p_monte_carlo_100000=pmc,min_expected=round(exp.min(),2),M_plus_N_rate_per_kb=round(1000*mn/L[3:].sum(),2),rest_rate_per_kb=round(1000*rest/L[:3].sum(),2) if False else round(1000*k[:3].sum()/L[:3].sum(),2),p_exact_binomial_M_plus_N_vs_length_share=pb))
C=pd.DataFrame(rowsC); C.to_csv(SD+'R6-C_chi_square.csv',index=False,encoding='utf-8-sig')
# ---------- E: translation efficiency vs capacity (5 ORFs)
act=pd.read_csv(SD+'G4_activity_vs_capacity_table.csv').set_index('region').loc[ORFS5]
k71=cnt['+1 and +2 (as drawn; 71 families)']; act['fam_per_kb_71']=1000*k71/act.length_nt.values
act['TE_CHX']=act.CHX_density_ppm_per_nt/act.RNA_rel_depth_batch3_5_8_16h_median; act['TE_HAR']=act.HAR_P_per_nt/act.RNA_rel_depth_batch3_5_8_16h_median
rowsE=[]
for y,yl in (('fam_per_kb_71','>40-aa +1/+2 families per kb'),('fam_per_kb_ruleA','>40-aa +1 families per kb')):
    for x in ('TE_CHX','TE_HAR','CHX_density_ppm_per_nt','RNA_rel_depth_batch3_5_8_16h_median'):
        xv=act[x].values; yv=act[y].values
        rho=stats.spearmanr(xv,yv).statistic; r=stats.pearsonr(np.log10(xv),yv).statistic
        perms=[stats.spearmanr(xv,p).statistic for p in itertools.permutations(yv)]
        pex=float(np.mean(np.abs(perms)>=abs(rho)-1e-12));
        rowsE.append(dict(capacity=yl,predictor=x,n=5,spearman_rho=round(rho,3),p_spearman_exact_permutation_two_sided=round(pex,4),pearson_r_log10x=round(r,3),p_pearson=round(stats.pearsonr(np.log10(xv),yv).pvalue,4),min_attainable_p=round(2/120,4)))
E=pd.DataFrame(rowsE); E.to_csv(SD+'R6-E_TE_correlation.csv',index=False,encoding='utf-8-sig'); act.reset_index().to_csv(SD+'R6-E_panel_data.csv',index=False,encoding='utf-8-sig')
pd.set_option('display.width',250); pd.set_option('display.max_columns',30); pd.set_option('display.max_colwidth',80)
print(ntrip); print(A.drop(columns=['min_expected']).to_string()); print(B.to_string()); print(C.T.to_string()); print(E.to_string()); print(act[['length_nt','fam_per_kb_71','TE_CHX','TE_HAR']].round(2))
