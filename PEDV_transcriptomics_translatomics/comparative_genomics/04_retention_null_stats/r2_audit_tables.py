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
import csv,json,collections,numpy as np
B=PEDV_K+r""; O=B+"Supplementary_Figure_L/"; W=O+"work/"
T=PEDV_P+r""
pool=json.load(open(W+"refpool.json"))['fams']; fm={F['fid']:F for F in pool}
feat=json.load(open(W+"featured.json")); from pedv_paths import check_sup16_featured; check_sup16_featured(feat); NAME=feat['names']; FEATC=set(feat['featured']); ORIG=set(feat['orig_AC']+feat['orig_D']); REV=set(feat['rev_rule'])
fig6=feat['fig6']; inv6={v:k for k,v in fig6.items()}
CST={12620:'ORF1a',20611:'ORF1b',24756:'S',25031:'ORF3',25598:'E',26286:'M',27623:'N'}
PF={"ORF1a":(293,12622),"ORF1b":(12592,20613),"S":(20610,24758),"ORF3":(24758,25033),"E":(25361,25600),"M":(25608,26288),"N":(26300,27625)}
supm={r['family_id']:r for r in csv.DictReader(open(T+"SupM_alt_ORF_evidence/tables/SupM_candidate_family_universe.csv",encoding='utf-8-sig'))}
tier={r['family_id']:r for r in csv.DictReader(open(T+"SupM_alt_ORF_evidence/tables/SupM_candidate_tier_assignments.csv",encoding='utf-8-sig'))}
f6=list(csv.DictReader((l for l in open(T+"global_perturbation_screen/plus1_orf_capacity/tables/plus1_ORF_families_gt40aa.tsv",encoding='utf-8') if not l.startswith('#')),delimiter='\t'))
det={r['family_id']:r['detectability_status'] for r in f6}
def wr(fn,rows,cols=None):
    cols=cols or list(rows[0])
    with open(O+fn,"w",encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t'); w.writeheader(); w.writerows(rows)
gcds={n:(a-1)%3 for n,(a,b) in PF.items()}
nullr={r['family_id']:r for r in csv.DictReader(open(W+"synonymous_null_by_family.csv",encoding='utf-8-sig'))}
VER="refdata_GIb_v1.1 (27,959 nt; PEDV_GIb.annotation.gtf)"
U=[];AR=[]
for F in pool:
    f0=F['isos'][0]; st=[x for x in F['isos'] if x['aa']>40]
    gt40=F['maxaa']>40; canon=F['stop'] in CST
    n6=inv6.get(F['fid']); in6=n6 is not None
    nonc=gt40 and not canon
    why=[]
    if not gt40: why.append("max peptide %d aa (<=40): outside >40 aa main set; table only"%F['maxaa'])
    if canon: why.append("shares stop codon with annotated %s (canonical-stop family): not an alternative-frame family; kept in tables, not drawn in A/B/D"%CST[F['stop']])
    if nonc and f0['category'] in('upstream ORF','intergenic-UTR ORF'): why.append("start outside every CDS (UTR/gap): drawn in A, excluded from per-CDS statistics and null")
    if nonc and not in6:
        w="not in Figure 6: "
        if f0['relative_frame']!='+1': w+="relative frame to parent "+(f0['relative_frame'] or 'n/a')+" (Figure 6 enumerated +1 only); "
        if f0['parent'] in('-',''): w+="start outside parent CDS"
        why.append(w)
    if canon and in6: why.append("in Figure 6 (E-x1): Figure 6 enumerated only starts inside parent E")
    nullin=F['fid'] in nullr
    if nonc and not nullin: why.append("excluded from synonymous null: not an internal alternative-frame ORF of a parent CDS")
    U.append(dict(family_id=F['fid'],reference_start_min=min(x['start'] for x in F['isos']),reference_start_max=max(x['start'] for x in F['isos']),stop_codon_first_nt=F['stop'],
        genomic_frame_static=F['frame'],parent_ORF_fixed_assignment=f0['parent'],parent_genomic_frame=gcds.get(f0['parent'],''),relative_frame_to_parent=f0['relative_frame'],
        max_peptide_aa=F['maxaa'],n_start_variants_all_ge8=F['n'],n_start_variants_gt40=len(st),length_tier=F['tier'],canonical_stop_family='yes' if canon else 'no',
        sharing_stop_with=CST.get(F['stop'],''),category=f0['category'],
        included_in_SupL_A='yes' if nonc else 'no',included_in_SupL_B='yes' if nonc else 'no',
        included_in_SupL_C_featured='yes' if F['fid'] in FEATC else 'no',included_in_SupL_C3_null='yes' if nullin else 'no',included_in_SupL_D='yes' if nonc else 'no',
        included_in_Figure6='yes' if in6 else 'no',Figure6_name=NAME.get(F['fid'],''),Figure6_original_featured_provisional='yes' if F['fid'] in ORIG else 'no',
        Figure6_revised_rule_pending='yes' if F['fid'] in REV else 'no',
        included_in_SupM='yes' if (in6 and n6 in supm) else 'no',
        included_in_SupN='yes (derived: Figure 6 set minus canonical-stop families; SupN source table not examined)' if (in6 and not canon) else 'no',
        Figure6_CHX_HAR_evaluability=det.get(n6,'') if in6 else '',SupM_evidence_tier=tier.get(n6,{}).get('evidence_tier','') if in6 else '',
        exclusion_reason='; '.join(why),coordinate_annotation_version=VER))
wr("SupL_candidate_universe_reconciliation.tsv",U)
for u in U:
    F=fm[u['family_id']]
    AR.append(dict(family_id=u['family_id'],parent_CDS=u['parent_ORF_fixed_assignment'],absolute_genomic_frame=u['genomic_frame_static'],parent_genomic_frame=u['parent_genomic_frame'],
      relative_alt_frame=u['relative_frame_to_parent'],stop_anchor_first_nt=u['stop_codon_first_nt'],stop_codon=F['isos'][0]['stop_codon'],
      candidate_starts_gt40='|'.join(str(x['start']) for x in F['isos'] if x['aa']>40),candidate_starts_all_ge8='|'.join(str(x['start']) for x in F['isos']),
      start_codons_gt40='|'.join(x['start_codon'] for x in F['isos'] if x['aa']>40),n_AUG_gt40=sum(1 for x in F['isos'] if x['aa']>40 and x['aug']=='yes'),
      n_near_cognate_gt40=sum(1 for x in F['isos'] if x['aa']>40 and x['aug']=='no'),predicted_lengths_aa_gt40='|'.join(str(x['aa']) for x in F['isos'] if x['aa']>40),
      overlapping_annotations=F['isos'][0]['overlapping'],category=u['category'],region_of_first_start=F['isos'][0]['region'],
      prf_or_ORF1ab_overlap='yes' if any(x['prf_overlap']=='yes' for x in F['isos']) else 'no',length_tier=u['length_tier'],
      canonical_stop_family=u['canonical_stop_family'],included_in_main_set=u['included_in_SupL_A'],exclusion_reason=u['exclusion_reason']))
wr("SupL_family_start_architecture.tsv",AR)
lin=list(csv.DictReader(open(O+"PEDV_ORF_conservation_by_lineage.csv",encoding='utf-8-sig')))
R=[]
for r in lin:
    if r['length_tier']!='gt40': continue
    F=fm[r['family_id']]
    R.append(dict(family_id=r['family_id'],canonical_stop_family='yes' if F['stop'] in CST else 'no',lineage_group=r['lineage_group'],n_genomes=r['n_genomes'],
      assessed_genomes=r['n_assessable'],retained_genomes=r['n_family_retained'],lost_eligible_architecture=r['n_family_lost'],
      not_assessable=r['n_not_assessable'],retention_rate=r['fraction_retained'],wilson95_lo=r['wilson95_lo'],wilson95_hi=r['wilson95_hi'],
      evidence_flag=r['evidence_flag'],dedup_clusters_assessed=r['n_clusters_assessable_dedup'],dedup_retention_rate=r['fraction_retained_dedup'],
      sensitivity_anyStart_retention_rate=r['fraction_retained_anyStart_sensitivity'],coding_potential_ge8aa_retention_rate=r['fraction_retained_ge8'],
      definition="stop codon at mapped position intact AND in frame AND >=1 reference start variant (>40 aa) still an eligible codon with no in-frame stop between start and stop; not-assessable (end-incomplete / ambiguous) excluded from denominator"))
wr("SupL_family_retention_by_lineage.tsv",R)
N=[]
ASSUME="parent-CDS amino-acid sequence fixed; single-nt synonymous substitutions only (no stop/Met/Trp changes); transitions weighted by empirical ts/tv (ORF1a); substitutions per replicate = each genome's observed synonymous codon changes in that CDS (binned by 5); 40 replicates per bin; same retention rule on the reference sequence; changes outside parent CDS and indels not simulated; codon usage not re-fitted; dinucleotide/RNA-structure and phylogenetic relatedness not controlled"
for fid,r in nullr.items():
    F=fm[fid]; f0=F['isos'][0]; pa,pb=PF[f0['parent']]
    N.append(dict(family_id=fid,parent_CDS=r['parent_cds'] if 'parent_cds' in r else r['parent'],n_assessable_genomes=r['n_assessable'],observed_retained=r['n_retained_observed'],observed_retention=r['observed_fraction'],
        expected_retention_null=r['expected_fraction_synonymous_null'],expected_count=r['expected_count'],null_sd_count=r.get('null_sd_count',''),
        empirical_percentile_of_observed=r.get('empirical_percentile_of_observed',''),z=r['z'],null_replicates_per_bin=40,
        crosses_parent_CDS_boundary='yes' if (F['stop']+2>pb or min(x['start'] for x in F['isos'] if x['aa']>40)<pa) else 'no',model_assumptions=ASSUME))
excl=[dict(family_id=u['family_id'],parent_CDS=u['parent_ORF_fixed_assignment'],model_assumptions="EXCLUDED from null: "+(u['exclusion_reason'] or 'not an internal alternative-frame ORF of a parent CDS')) for u in U if u['length_tier']=='gt40' and u['family_id'] not in nullr]
wr("SupL_synonymous_null_comparison.tsv",N+excl,list(N[0]))
typ=list(csv.DictReader(open(W+"typing_S_based.tsv",encoding='utf-8'),delimiter='\t'))
qc={r['accession_as_in_1083_file']:r for r in csv.DictReader(open(O+"PEDV_genome_database_QC.csv",encoding='utf-8-sig'))}
meta={}
for r in csv.DictReader(open(B+"N1_conservation_analysis/data_external/PEDV_metadata.tsv",encoding='utf-8'),delimiter='\t'): meta[r['accession']]=r
pg={r['accession']:r for r in csv.DictReader(open(O+"SourceData_PanelD_per_genome_completeness.csv",encoding='utf-8-sig'))}
LA=[]
for t in typ:
    a=t['accession']; q=qc[a]; m=meta.get(a.split('.')[0],{})
    LA.append(dict(accession=a,S_type_this_study=t['S_type'],assignment_basis=t['S_type_basis'],nearest_panel_member=t['nearest_panel_member'],nearest_p_distance=t['nearest_p_distance'],
      S_gene_coverage=t['S_coverage'],unresolved='yes' if t['S_type']=='unresolved' else 'no',historical_label_original_database=t['N1_genotype_historical'],historical_status=t['N1_status'],
      recombination_flag_source=m.get('recombination_flag',''),genome_length=q['genome_length'],end_completeness=q['end_completeness'],non_ACGT_fraction=q['non_ACGT_fraction'],
      near_duplicate_cluster=q['near_duplicate_cluster_le0p1pct'],n_alt_families_assessable=pg.get(a,{}).get('n_assessable',''),per_genome_retention_fraction=pg.get(a,{}).get('fraction',''),
      country=t['country'],collection_date=t['collection_date'],strain=t['strain']))
wr("SupL_lineage_curation_audit.tsv",LA)
cnt=[]
def add(step,fam,iso,note): cnt.append(dict(step=step,n_families=fam,n_start_variants=iso,note=note))
allp=[(F,x) for F in pool for x in F['isos']]
add("1 all-frame complete ORFs >=8 aa, 8 start codons",len(pool),sum(F['n'] for F in pool),"sequence-defined universe (stop-anchored) on refdata_GIb_v1.1")
g40=[F for F in pool if F['maxaa']>40]; add("2 families with max peptide >40 aa",len(g40),sum(1 for F,x in allp if x['aa']>40),"start variants counted at >40 aa")
cs=[F for F in g40 if F['stop'] in CST]; add("3 of which share a stop with an annotated CDS",len(cs),sum(1 for F in cs for x in F['isos'] if x['aa']>40),"canonical-stop families; includes E-x1 (= M in-frame extension), an S-frame extension through ORF1b, ORF1b-frame extension through ORF1a")
nc=[F for F in g40 if F['stop'] not in CST]; add("4 non-canonical-stop families >40 aa (Panels A, B, D)",len(nc),sum(1 for F in nc for x in F['isos'] if x['aa']>40),"supersedes earlier 76 (3 canonical-stop families had been mislabelled) and 77/798 caption (upstream ORF counted twice; 798 summed variants of mislabelled families)")
inn=[F for F in nc if F['isos'][0]['category']=='internal alternative-frame ORF']; add("5 internal alternative-frame families",len(inn),sum(1 for F in inn for x in F['isos'] if x['aa']>40),"first start inside a parent CDS")
p1=[F for F in inn if F['isos'][0]['relative_frame']=='+1']; add("6 of which relative +1 to fixed parent",len(p1),sum(1 for F in p1 for x in F['isos'] if x['aa']>40),"Figure 6 additionally required every start inside the parent (567 variants)")
add("7 families in synonymous null (C3)",len(nullr),'',"internal alternative-frame families of a parent CDS (earlier 74 included 3 canonical-stop families)")
add("8 Figure 6 / Sup M families",len(fig6),567,"parent +1 frame, start inside parent CDS, >40 aa, AUG/near-cognate (Figure 6 table)")
add("9 Figure 6 families that are non-canonical-stop in L",sum(1 for k in fig6.values() if fm[k]['stop'] not in CST),'',"= Sup N main analysis (68); E-x1 excluded as the M-frame in-frame extension")
add("10 Figure 6 CHX/HAR evaluable families",sum(1 for r in f6 if r['detectability_status']=='evaluable'),'',"Figure 6 table detectability_status (label only; not used in L)")
add("11 L non-canonical families not in Figure 6",len([F for F in nc if F['fid'] not in inv6]),'',"relative +2 frame, upstream/UTR or start outside parent: not enumerated in Figure 6")
wr("SupL_candidate_universe_counts_by_step.tsv",cnt)
for c in cnt: print(c['step'],'|',c['n_families'],c['n_start_variants'])
