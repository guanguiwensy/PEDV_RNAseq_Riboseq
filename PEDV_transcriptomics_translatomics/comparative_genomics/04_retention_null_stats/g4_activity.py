"""G4 'translation activity' descriptors vs long +1 ORF capacity. n=7 canonical ORFs: Spearman is DESCRIPTIVE ONLY (no significance statement).
Ribo: Fig6_canonical_background_audit.tsv (CHX_P_total, HAR_P_total; strict_core 29/30 nt, offset +12, canonical-frame normalisation of Fig6). RNA: GIb_rnaseq_reanalysis/coverage/*.depth.tsv (23 libraries, non-strand-specific)."""
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
import sample_registry as SR
from r6lib import *
from scipy import stats
BT=PEDV_P+""
au=pd.read_csv(BT+"SupM_alt_ORF_evidence/Fig6_reaudit/tables/Fig6_canonical_background_audit.tsv",sep='\t').set_index('parent_ORF')
man=pd.read_csv(BT+"GIb_rnaseq_reanalysis/sample_manifest_GIb.tsv",sep='\t')
man=man[man.sample_id.isin(SR.WHITELIST_LIBRARIES)].reset_index(drop=True); SR.assert_exactly(man.sample_id, SR.RNA_VERO_ALIASES, "Sup15 RNA depth")
rel={}
for s in man.sample_id:
    d=np.loadtxt(BT+f"GIb_rnaseq_reanalysis/coverage/{s}.depth.tsv",usecols=2); g=d.mean()
    rel[s]={n:d[a-1:b].mean()/g for n,(a,b) in CDS.items()}
rel=pd.DataFrame(rel).T
rel_all=rel.median(); b3=man[(man.batch=='batch3')&man.time_h.isin([5,8,16])].sample_id; rel_b3=rel.loc[b3].median()
cap=pd.read_csv(OUT+'source_data/G4_capacity_by_ORF.csv')
rows=[]
for n in ORDER:
    L=CDS[n][1]-CDS[n][0]+1
    ca=cap[(cap.rule=='A_SupL_68')&(cap.region==n)].iloc[0]; cb=cap[(cap.rule=='B_Fig6_69')&(cap.region==n)].iloc[0]
    rows.append(dict(region=n,length_nt=L,CHX_P_total=au.loc[n,'CHX_P_total'],CHX_density_ppm_per_nt=au.loc[n,'CHX_canonical_density_ppm_per_nt'],HAR_P_total=au.loc[n,'HAR_P_total'],HAR_P_per_nt=au.loc[n,'HAR_P_total']/L,
        RNA_rel_depth_all20_median=rel_all[n],RNA_rel_depth_batch3_5_8_16h_median=rel_b3[n],CHX_ppm_over_RNA=au.loc[n,'CHX_canonical_density_ppm_per_nt']/rel_b3[n],
        fam_per_kb_ruleA=ca.fam_per_kb,fam_per_kb_ruleB=cb.fam_per_kb,iso_per_kb_ruleA=ca.iso_per_kb,covered_fraction_ruleA=ca.covered_fraction))
A=pd.DataFrame(rows); A.to_csv(OUT+'source_data/G4_activity_vs_capacity_table.csv',index=False,encoding='utf-8-sig')
C=[]
for act in ('CHX_density_ppm_per_nt','HAR_P_per_nt','RNA_rel_depth_all20_median','RNA_rel_depth_batch3_5_8_16h_median','CHX_ppm_over_RNA'):
    for capc in ('fam_per_kb_ruleA','fam_per_kb_ruleB','iso_per_kb_ruleA','covered_fraction_ruleA'):
        r=stats.spearmanr(A[act],A[capc]); C.append(dict(activity=act,capacity=capc,n=7,spearman_rho=round(r.statistic,3),descriptive_only_p_not_to_be_reported=round(r.pvalue,3)))
pd.DataFrame(C).to_csv(OUT+'source_data/G4_activity_spearman_descriptive.csv',index=False,encoding='utf-8-sig')
pd.set_option('display.width',250); print(A.round(3).to_string()); print(pd.DataFrame(C).to_string())
