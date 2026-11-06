"""G3 (counts, reconciliation of 68/69, two attribution rules) and G4 (capacity per canonical ORF; Poisson/binomial tests).
Rule A = Sup L (fixed_parent_assignment of the family's first isoform; stop shared with a canonical CDS stop -> that CDS, frame 0).
Rule B = Figure 6 (parent = CDS containing the start, +1 frame only)."""
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
from r6lib import *
from scipy import stats
f=families(); iso=isoforms(); f6=fig6fam()
f['stop']=f.shared_stop_position
f['relf']=f.relative_frame_to_parent
iso=iso.merge(f[['family_id','fixed_parent_assignment','relf']].rename(columns={'fixed_parent_assignment':'famparent'}),on='family_id')
tiers={'ge8':8,'ge20':20,'gt40':40}
rows=[]
def cnt(frame,reg,thr):
    ff=f[(f.relf==frame)&(f.fixed_parent_assignment==reg)&(f.max_peptide_aa>thr if thr>=40 else f.max_peptide_aa>=thr)]
    ii=iso[(iso.relf==frame)&(iso.famparent==reg)&(iso.family_id.isin(ff.family_id))&(iso.peptide_length_aa>thr if thr>=40 else iso.peptide_length_aa>=thr)]
    return len(ff),len(ii)
for frame in (1,2):
    for reg in ORDER:
        for t,thr in tiers.items():
            nf,ni=cnt(frame,reg,thr)
            rows.append(dict(rule='A_SupL',rel_frame=f'+{frame}',region=reg,tier=t,n_families=nf,n_start_isoforms=ni))
G3=pd.DataFrame(rows)
G3.to_csv(OUT+'source_data/G3_counts_by_frame_region_tier_ruleA.csv',index=False,encoding='utf-8-sig')
# Rule B / Figure 6 and crosswalk by stop coordinate
f6['stop']=f6.shared_stop
A=f[(f.relf==1)&(f.max_peptide_aa>40)][['family_id','stop','fixed_parent_assignment','n_start_variants','max_peptide_aa']].rename(columns={'family_id':'SupL_family_id','fixed_parent_assignment':'ruleA_region','n_start_variants':'SupL_n_starts_ge8','max_peptide_aa':'SupL_max_aa'})
Bt=f6[['family_id','parent_ORF','stop','n_start_variants','max_ORF_length_aa']].rename(columns={'family_id':'Fig6_family_id','parent_ORF':'ruleB_region','n_start_variants':'Fig6_n_starts_gt40','max_ORF_length_aa':'Fig6_max_aa'})
assert A.stop.is_unique and Bt.stop.is_unique
X=A.merge(Bt,on='stop',how='outer')
# SupL families sharing the stop of Fig6 families but not +1 (E-x1)
extra=f[f.stop.isin(X[X.SupL_family_id.isna()].stop)][['stop','family_id','fixed_parent_assignment','relf','n_start_variants','max_peptide_aa']]
X.loc[X.SupL_family_id.isna(),'note']='Fig6-only; in Sup L this stop-anchored family is '+extra.set_index('stop').reindex(X[X.SupL_family_id.isna()].stop).apply(lambda r:f"{r.family_id}: parent {r.fixed_parent_assignment}, rel frame {int(r.relf)}, {r.n_start_variants} starts (>=8aa)",axis=1).values
X['region_agree']=(X.ruleA_region==X.ruleB_region)
X.sort_values('stop').to_csv(OUT+'source_data/G3_crosswalk_ruleA_vs_ruleB_plus1_gt40.csv',index=False,encoding='utf-8-sig')
print('A +1>40:',len(A),'B:',len(Bt),'disagree:',X[~X.region_agree].to_dict('records'))
# isoform reconciliation 567
E1=iso[iso.family_id=='SDF_f2_stop26286']
print('Fig6 start variants total',f6.n_start_variants.sum(),'| rule A +1 >40 isoforms',iso[(iso.relf==1)&(iso.family_id.isin(A.SupL_family_id))&(iso.peptide_length_aa>40)].shape[0])
print('E-x1 isoforms >40 in SupL family',(E1.peptide_length_aa>40).sum(),'of',len(E1),'; those with start in E:',E1[(E1.start_position>=25361)&(E1.start_position<=25600)&(E1.peptide_length_aa>40)].shape[0])
# G4 capacity
LEN={n:b-a+1 for n,(a,b) in CDS.items()}  # Sup L / Fig6 lengths incl. ORF1a/1b overlap
def garwood(k,n=None):
    lo=0 if k==0 else stats.chi2.ppf(0.025,2*k)/2; hi=stats.chi2.ppf(0.975,2*k+2)/2; return lo,hi
def cap(famdf,isodf,label,lengths,rate_ref=None):
    out=[]
    tot=len(famdf); L=sum(lengths.values())
    for r in ORDER:
        k=int((famdf.region==r).sum()); Lr=lengths[r]
        lo,hi=garwood(k)
        # nt coverage: union of [first start, stop+2] intervals intersected with region
        sub=famdf[famdf.region==r]; iv=[]
        a,b=CDS[r]
        for _,x in sub.iterrows():
            s0=int(x.first_start); iv.append((max(s0,a),min(int(x.stop)+2,b)))
        iv=sorted(v for v in iv if v[0]<=v[1]); u=0; cur=None
        for s,e in iv:
            if cur is None or s>cur[1]+1:
                if cur: u+=cur[1]-cur[0]+1
                cur=[s,e]
            else: cur[1]=max(cur[1],e)
        if cur: u+=cur[1]-cur[0]+1
        rest_k=tot-k; rest_L=L-Lr
        pv=stats.binomtest(k,tot,Lr/L).pvalue if tot>0 else np.nan   # region vs uniform allocation of all families
        ni=int((isodf.region==r).sum())
        out.append(dict(rule=label,region=r,length_nt=Lr,n_families=k,fam_per_kb=round(1000*k/Lr,3),fam_per_kb_CI95_lo=round(1000*lo/Lr,3),fam_per_kb_CI95_hi=round(1000*hi/Lr,3),
          n_start_isoforms=ni,iso_per_kb=round(1000*ni/Lr,2),covered_nt=u,covered_fraction=round(u/Lr,3),
          rest_fam_per_kb=round(1000*rest_k/rest_L,3),rate_ratio_vs_rest=round((k/Lr)/(rest_k/rest_L),3) if rest_k else np.nan,
          expected_if_uniform=round(tot*Lr/L,2),exact_binomial_p_two_sided=round(pv,4)))
    return pd.DataFrame(out)
# frame +1 >40 families under rule A and rule B, with first-start coordinate
first_start=iso.groupby('family_id').start_position.min()
fa=f[(f.relf==1)&(f.max_peptide_aa>40)].copy(); fa['region']=fa.fixed_parent_assignment; fa['first_start']=fa.family_id.map(first_start)
ia=iso[(iso.family_id.isin(fa.family_id))&(iso.peptide_length_aa>40)].copy(); ia['region']=ia.family_id.map(fa.set_index('family_id').region)
# rule B
fb=f6.copy(); fb['region']=fb.parent_ORF; fb['first_start']=fb.all_candidate_starts.str.split('|').map(lambda x:min(int(v) for v in x))
ib=fb.copy(); ib['n']=ib.n_start_variants
ibi=pd.DataFrame({'region':np.repeat(fb.region.values,fb.n_start_variants.values)})
CA=cap(fa,ia,'A_SupL_68',LEN); CB=cap(fb,ibi,'B_Fig6_69',LEN)
LEN2=dict(LEN); LEN2['ORF1a']-=31; LEN2['ORF1b']-=31
CA2=cap(fa,ia,'A_SupL_68_overlap31nt_removed_from_lengths',LEN2)
cap_all=pd.concat([CA,CB,CA2]); cap_all.to_csv(OUT+'source_data/G4_capacity_by_ORF.csv',index=False,encoding='utf-8-sig')
pd.set_option('display.width',250)
print(CA.to_string()); print(CB.to_string())
# same for +2 and for 0 frame counts (small)
print(G3[(G3.tier=='gt40')&((G3.n_families>0))].to_string())
