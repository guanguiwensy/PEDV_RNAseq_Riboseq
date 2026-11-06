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
import csv,json,collections,os
from openpyxl import Workbook
from openpyxl.styles import Font,Alignment
B=PEDV_K+r""; O=B+"Supplementary_Figure_L/"; W=O+"work/"
T=PEDV_P+r""
TS=O+"SupTable_TSV/"; os.makedirs(TS,exist_ok=True)
def rd(fn,delim=None,enc='utf-8-sig'):
    d=delim or ('\t' if fn.endswith('.tsv') else ',')
    return list(csv.DictReader(open(fn,encoding=enc),delimiter=d))
U=rd(O+"SupL_candidate_universe_reconciliation.tsv")
alt=[u for u in U if u['included_in_SupL_A']=='yes']          # 73 alternative-frame families (>40 aa)
ALT={u['family_id'] for u in alt}
feat=json.load(open(W+"featured.json")); C12=feat['C12_list']; REAS=feat['C12_reason']; ORIG=set(feat['orig_AC']+feat['orig_D'])
sel={r['family_id']:r for r in rd(O+"SupL_C1C2_family_selection.csv")}
name6={r['original_family_id']:r['displayed_name'] for r in rd(T+"global_perturbation_screen/plus1_orf_capacity/panelC_family_name_mapping.tsv",delim='\t',enc='utf-8')}
fig6=feat['fig6']; sdf2p1={v:k for k,v in fig6.items()}
v2={r['family']:r for r in rd(T+"SupM_alt_ORF_evidence/Fig6_checkpoint3/tables/Fig6_joint_evidence_classification_v2.tsv")}
LR=[r for r in rd(O+"SupL_family_retention_by_lineage.tsv") if r['family_id'] in ALT]
L={(r['family_id'],r['lineage_group']):r for r in LR}
lin_groups=['GI-a','GI-b','GII-a','GII-b','GII-c','unresolved','ALL','N1_historical_G1b_label']
# ---------- 1 Family_summary
S1=[]
for u in sorted(alt,key=lambda u:int(u['reference_start_min'])):
    f=u['family_id']; p1=sdf2p1.get(f,''); nm=name6.get(p1,'')
    arch=v2.get(nm,{})
    row=dict(family_id=f,Figure6_family_id=p1,Figure6_name=nm,parent_ORF=u['parent_ORF_fixed_assignment'],absolute_genomic_frame=u['genomic_frame_static'],parent_genomic_frame=u['parent_genomic_frame'],
        relative_alt_frame=u['relative_frame_to_parent'],category=u['category'],stop_anchor_first_nt=u['stop_codon_first_nt'],start_min=u['reference_start_min'],start_max=u['reference_start_max'],
        max_candidate_length_aa=u['max_peptide_aa'],n_start_variants_gt40=u['n_start_variants_gt40'])
    for g in lin_groups:
        r=L[(f,g)]; row[f'retention_{g}']=r['retention_rate']; row[f'n_assessable_{g}']=r['assessed_genomes']
    row['retention_ALL_dedup_clusters']=L[(f,'ALL')]['dedup_retention_rate']
    row['retention_ALL_anyStart_sensitivity']=L[(f,'ALL')]['sensitivity_anyStart_retention_rate']
    row['in_Figure6_Sup_M']=u['included_in_Figure6']; row['in_SupN_main']=u['included_in_SupN'][:3]
    row['Figure6_original_featured']='yes' if f in ORIG else 'no'
    row['Figure6_provisional_featured_C1C2']='yes' if f in C12 else 'no'
    row['C1C2_selection_basis']=REAS.get(f,'')
    row['Figure6_evidence_architecture_v2_label_only']=arch.get('evidence_architecture','')
    row['note']='Ribo-seq labels are descriptive and were not used in any sequence-retention calculation'
    S1.append(row)
# ---------- 2 Start_variant_conservation
sc=[r for r in rd(O+"SupL_start_codon_conservation.tsv") if r['family_id'] in ALT]
cls=['exact_codon_retention','AUG_retention','near_cognate_retention','start_class_switching','loss_of_eligible_start','not_assessable']
all_=collections.defaultdict(lambda:collections.Counter()); meta={}
for r in sc:
    meta[(r['family_id'],r['start_position'])]=r
    for c in cls: all_[(r['family_id'],r['start_position'])][c]+=int(r[c])
S2=[]
for r in sc:
    d=dict(family_id=r['family_id'],start_position=r['start_position'],reference_start_codon=r['reference_start_codon'],start_class='AUG' if r['reference_start_codon']=='AUG' else 'near-cognate',
           candidate_length_aa=r['reference_peptide_aa'],lineage=r['lineage'],n_genomes_in_lineage=r['n_genomes_in_lineage'])
    tot=sum(int(r[c]) for c in cls)
    for c in cls: d[c]=r[c]; d[c+'_fraction']=round(int(r[c])/tot,4) if tot else ''
    d['n_assessable']=tot-int(r['not_assessable']); S2.append(d)
for (f,s),c in all_.items():
    r=meta[(f,s)]; tot=sum(c.values())
    d=dict(family_id=f,start_position=s,reference_start_codon=r['reference_start_codon'],start_class='AUG' if r['reference_start_codon']=='AUG' else 'near-cognate',
           candidate_length_aa=r['reference_peptide_aa'],lineage='ALL',n_genomes_in_lineage=1083)
    for k in cls: d[k]=c[k]; d[k+'_fraction']=round(c[k]/tot,4) if tot else ''
    d['n_assessable']=tot-c['not_assessable']; S2.append(d)
# ---------- 3 Family_retention_by_lineage
S3=[]
for r in LR:
    S3.append(dict(family_id=r['family_id'],lineage_group=r['lineage_group'],n_genomes=r['n_genomes'],assessable_genomes=r['assessed_genomes'],retained_genomes=r['retained_genomes'],
        lost_eligible_architecture=r['lost_eligible_architecture'],not_assessable_genomes=r['not_assessable'],retention_fraction=r['retention_rate'],wilson95_lo=r['wilson95_lo'],wilson95_hi=r['wilson95_hi'],
        evidence_flag=r['evidence_flag'],dedup_clusters_assessed=r['dedup_clusters_assessed'],dedup_retention_fraction=r['dedup_retention_rate'],
        anyStart_sensitivity_fraction=r['sensitivity_anyStart_retention_rate'],coding_potential_ge8aa_fraction=r['coding_potential_ge8aa_retention_rate']))
# ---------- 4 Synonymous_null_comparison
NL=rd(O+"SupL_synonymous_null_comparison.tsv")
S4=[r for r in NL if r['family_id'] in ALT]
inn={r['family_id'] for r in S4 if r.get('observed_retention')}
for u in alt:
    if u['family_id'] not in inn and not any(r['family_id']==u['family_id'] for r in S4):
        S4.append(dict(family_id=u['family_id'],parent_CDS=u['parent_ORF_fixed_assignment'],model_assumptions='EXCLUDED: '+u['exclusion_reason']))
for r in S4: r['Figure6_name']=name6.get(sdf2p1.get(r['family_id'],''),'')
# ---------- 5 Genome_level_retention
hm=rd(O+"PEDV_family_homology_mapping.csv")
typ={r['accession']:r for r in rd(W+"typing_S_based.tsv",delim='\t',enc='utf-8')}
S5=[]
for r in hm:
    if r['family_id'] not in ALT: continue
    st=r['hom_status_gt40']; st={'NA':'not_assessable'}.get(st,st)
    S5.append(dict(accession=r['accession'],S_gene_lineage=r['S_type'],historical_label=typ[r['accession']]['N1_genotype_historical'],family_id=r['family_id'],
        retention_status=st,anyStart_status=r['status_gt40'],reason=r['reason'],max_peptide_aa_in_genome=r['max_peptide_aa_in_genome'],n_reference_starts_gt40=r['n_ref_starts_gt40'],
        n_reference_starts_retained=r['n_hom_starts_retained_gt40'],**{k:r[k+'_gt40'] for k in cls}))
# ---------- 6 Definitions_and_reconciliation
cnt=rd(O+"SupL_candidate_universe_counts_by_step.tsv")
S6=[]
def sec(k,v,n=''): S6.append(dict(section=k,item=v,value=n,detail=''))
for c in cnt: S6.append(dict(section='Counts by step',item=c['step'],value=f"{c['n_families']} families / {c['n_start_variants']} start variants",detail=c['note']))
excl=collections.defaultdict(list)
for u in U:
    if u['length_tier']!='gt40': continue
    if u['canonical_stop_family']=='yes': excl['canonical-stop (not alternative-frame; shares stop with annotated CDS)'].append(u['family_id']+'('+u['sharing_stop_with']+')')
    elif u['included_in_SupL_C3_null']!='yes': excl['alternative-frame but not in synonymous null (C3): '+u['category']].append(u['family_id'])
    if u['included_in_SupL_A']=='yes' and u['included_in_Figure6']!='yes': excl['in L (73) but not in Figure 6/M (69)'].append(u['family_id']+' ['+u['relative_frame_to_parent']+', '+u['category']+']')
    if u['included_in_Figure6']=='yes' and u['canonical_stop_family']=='yes': excl['in Figure 6/M (69) but not in Sup N main heatmap (68)'].append(u['family_id']+' ('+u['Figure6_name']+')')
for k,v in excl.items(): S6.append(dict(section='Excluded / not shared families',item=k,value=len(v),detail='; '.join(v)))
defs=[("Sequence-defined family","All AUG/CUG/GUG/UUG/AUA/AUU/AUC/ACG starts that share one in-frame stop codon on refdata_GIb_v1.1; a family is not a translated ORF and many starts are not many proteins."),
("Alternative-frame family (73)","Family with >40 aa isoforms whose stop codon is not the stop of an annotated CDS."),
("Absolute vs relative frame","Absolute frame=(position-1) mod 3; relative frame = (absolute - parent CDS frame) mod 3. ORF1b +1 = absolute frame 1 = ORF1a coding frame."),
("Start-codon classes","exact non-AUG codon; AUG; different near-cognate; start-class switch (AUG<->near-cognate); loss of eligible start; not assessable. Denominator for fractions = all reference start variant x genome pairs in the lineage."),
("Family retention","Stop codon intact in frame at the alignment-mapped position AND >=1 reference start variant (>40 aa) remains an eligible codon in the same frame with no in-frame stop between start and stop. Not-assessable (end-incomplete/ambiguous) excluded from denominator. AnyStart sensitivity: any eligible start in the open stretch."),
("Not assessable","Mapped stop or start codon missing because the genome is end-incomplete/unaligned, or contains ambiguous bases; never counted as loss."),
("Lineage (curated S-gene)","Nearest S-gene panel member (SGP v1.0, Guo 2019 anchors); p-distance<=0.03; >=3,400 S positions; ties within 0.003 across labels -> unresolved. GI-a n=23, GI-b n=30, GII-a 752, GII-b 96, GII-c 158, unresolved 24."),
("Historical label","Database label 'G1b' (n=168; older nomenclature, mostly S-INDEL/GII-c strains). Not the same classification as curated GI-b; shown only as a separate sensitivity group."),
("Synonymous null","Parent-CDS protein fixed; per-genome number of synonymous changes; single-nt synonymous substitutions weighted by empirical ts/tv; 40 replicates per bin; no control of overlapping RNA constraints, codon re-fitting or phylogeny."),
("Figure 6 featured","Original Figure 6 featured = ORF1a-x22, N-x1, N-x2 (A/C), ORF1a-x3, ORF1b-x13, ORF1b-x8 (D). C1/C2 list = provisional checkpoint-3 evidence-architecture list (selection basis column); Figure 6 not frozen."),
("Independence","Sequence retention was computed without any Ribo-seq quantity (HAR D, CHX); Figure 6 labels were attached afterwards.")]
for k,v in defs: S6.append(dict(section='Definitions',item=k,value='',detail=v))
# ---------- write
sheets=[("Family_summary",S1),("Start_variant_conservation",S2),("Family_retention_by_lineage",S3),("Synonymous_null_comparison",S4),("Genome_level_retention",S5),("Definitions_and_reconciliation",S6)]
wb=Workbook(write_only=True)
for name,rows in sheets:
    cols=[]
    for r in rows:
        for k in r:
            if k not in cols: cols.append(k)
    ws=wb.create_sheet(name)
    ws.append(cols)
    for r in rows: ws.append([r.get(c,'') for c in cols])
    with open(TS+f"{name}.tsv","w",encoding='utf-8',newline='') as fh:
        w=csv.writer(fh,delimiter='\t'); w.writerow(cols); [w.writerow([r.get(c,'') for c in cols]) for r in rows]
    print(name,len(rows),len(cols))
wb.save(O+"Supplementary_Table_PEDV_Alternative_ORF_Conservation.xlsx")
