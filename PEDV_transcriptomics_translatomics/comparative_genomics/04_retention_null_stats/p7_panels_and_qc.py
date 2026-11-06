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
import sys,csv,json,collections,math,numpy as np
from scipy import stats
B=PEDV_K+r""; O=B+"Supplementary_Figure_L/"; W=O+"work/"
pool=json.load(open(W+"refpool.json"))['fams']; fm={F['fid']:F for F in pool}
feat=json.load(open(W+"featured.json")); from pedv_paths import check_sup16_featured; check_sup16_featured(feat); FEAT=feat['featured']; NAME=feat['names']; ORIG=set(feat['orig_AC']+feat['orig_D']); REV=set(feat['rev_rule'])
typing={r['accession']:r for r in csv.DictReader(open(W+"typing_S_based.tsv",encoding='utf-8'),delimiter='\t')}
lin=list(csv.DictReader(open(O+"PEDV_ORF_conservation_by_lineage.csv",encoding='utf-8-sig')))
L={(r['family_id'],r['lineage_group']):r for r in lin}
nullr={r['family_id']:r for r in csv.DictReader(open(W+"synonymous_null_by_family.csv",encoding='utf-8-sig'))}
CDS={"ORF1a":(293,12622),"ORF1b":(12592,20613),"S":(20610,24758),"ORF3":(24758,25033),"E":(25361,25600),"M":(25608,26288),"N":(26300,27625)}
def wr(fn,rs,cols=None):
    cols=cols or list(rs[0])
    with open(O+fn,"w",encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols); w.writeheader(); w.writerows(rs)
CSTOPS={12620,20611,24756,25031,25598,26286,27623}
INFR={'in-frame extension','in-frame internal initiation','annotated canonical CDS start','in-frame within CDS (stop differs)'}
G40=[F for F in pool if F['maxaa']>40]
for F in G40: F['inframe']=F['stop'] in CSTOPS
# ---------- Panel A
A=[]
for F in sorted(G40,key=lambda F:F['isos'][0]['start']):
    st=[x for x in F['isos'] if x['aa']>40]; f0=F['isos'][0]
    A.append(dict(family_id=F['fid'],figure6_name=NAME.get(F['fid'],''),genomic_frame_static=F['frame'],start_min=min(x['start'] for x in st),start_max=max(x['start'] for x in st),
        stop_position=F['stop'],n_start_variants_gt40=len(st),n_AUG=sum(1 for x in st if x['aug']=='yes'),n_near_cognate=sum(1 for x in st if x['aug']=='no'),
        contains_AUG='yes' if any(x['aug']=='yes' for x in st) else 'no',max_aa=max(x['aa'] for x in st),
        parent=f0['parent'],relative_frame=f0['relative_frame'],category=f0['category'],canonical_frame_family='yes' if F['inframe'] else 'no',
        featured='yes' if F['fid'] in FEAT else 'no',fig6_original_featured_provisional='yes' if F['fid'] in ORIG else 'no',fig6_revised_rule_pending='yes' if F['fid'] in REV else 'no',overlaps_PRF_or_ORF1ab_overlap='yes' if any(x['prf_overlap']=='yes' for x in F['isos']) else 'no'))
wr("SourceData_PanelA_families_gt40.csv",A)
# density tracks (non-canonical-frame families only), window 1000 / step 100
isoall=[(x['start'],F['fid']) for F in G40 if not F['inframe'] for x in F['isos'] if x['aa']>40]
famstop=[F['stop'] for F in G40 if not F['inframe']]
famfirst=[min(x['start'] for x in F['isos'] if x['aa']>40) for F in G40 if not F['inframe']]
tr=[]
for c in range(500,27959-500+1,100):
    a,b=c-500,c+500
    tr.append(dict(window_center=c,families_per_kb=sum(1 for s in famfirst if a<=s<b)/1.0,start_sites_per_kb=sum(1 for s,_ in isoall if a<=s<b)/1.0))
wr("SourceData_PanelA_density_tracks.csv",tr)
# ---------- Panel B
NC=[F for F in G40 if not F['inframe']]
regions=list(CDS)+['UTR/gap']
reg_len={k:b-a+1 for k,(a,b) in CDS.items()}
tot_cds=sum(reg_len.values())
Brow=[];spf={}
for r in regions:
    fs=[F for F in NC if (F['isos'][0]['parent']==r and F['isos'][0]['category'] not in('upstream ORF','intergenic-UTR ORF') if r!='UTR/gap' else F['isos'][0]['category'] in('upstream ORF','intergenic-UTR ORF') or F['isos'][0]['parent'] in('-',''))]
    # UTR/gap: families whose first start is outside every CDS
    ns=[sum(1 for x in F['isos'] if x['aa']>40) for F in fs]
    spf[r]=ns
    ln=reg_len.get(r,None)
    loc=sum(1 for F in NC for x in F['isos'] if x['aa']>40 and ((x['parent']==r and x['category'] not in('upstream ORF','intergenic-UTR ORF')) if r!='UTR/gap' else (x['category'] in('upstream ORF','intergenic-UTR ORF') or x['parent'] in('-',''))))
    Brow.append(dict(region=r,length_nt=ln if ln else '',n_families=len(fs),n_start_isoforms=sum(ns),n_start_sites_by_start_location=loc,mean_starts_per_family=round(float(np.mean(ns)),2) if ns else '',
        median_starts_per_family=float(np.median(ns)) if ns else '',max_starts_per_family=max(ns) if ns else '',
        families_per_kb=round(len(fs)/ln*1000,3) if ln else '',start_sites_per_kb=round(loc/ln*1000,2) if ln else '',
        max_family_share_of_region_starts=round(max(ns)/sum(ns),3) if ns else ''))
# genome-wide rates over CDS regions
tf=sum(b['n_families'] for b in Brow if b['region'] in CDS); ti=sum(b['n_start_isoforms'] for b in Brow if b['region'] in CDS)
for b in Brow:
    if b['length_nt']:
        ef=tf*b['length_nt']/tot_cds; ei=ti*b['length_nt']/tot_cds
        b['expected_families_if_uniform']=round(ef,2); b['family_count_vs_uniform_poisson_p_two_sided']=round(float(stats.poisson.sf(b['n_families']-1,ef)*2 if b['n_families']>ef else stats.poisson.cdf(b['n_families'],ef)*2),4)
        b['expected_starts_if_uniform']=round(ei,1)
wr("SourceData_PanelB_region_density.csv",Brow)
bs=[]
for r,ns in spf.items():
    for n in ns: bs.append(dict(region=r,n_start_variants_gt40=n))
wr("SourceData_PanelB_starts_per_family.csv",bs)
# N vs rest
rest=[n for r,ns in spf.items() if r not in('N','UTR/gap') for n in ns]
Nn=spf['N']
mw=stats.mannwhitneyu(Nn,rest,alternative='two-sided')
rng=np.random.default_rng(1)
allf=np.array(Nn+rest); perm=[]
obs=np.mean(Nn)-np.mean(rest)
for _ in range(20000):
    rng.shuffle(allf); perm.append(allf[:len(Nn)].mean()-allf[len(Nn):].mean())
perm_p=float((np.abs(perm)>=abs(obs)).mean())
famrate_N=len(Nn)/reg_len['N']*1000; famrate_rest=sum(len(spf[r]) for r in CDS if r!='N')/(tot_cds-reg_len['N'])*1000
# exact conditional test for rate ratio (Poisson): binomial
nN=len(Nn); nR=sum(len(spf[r]) for r in CDS if r!='N')
pbin=stats.binomtest(nN,nN+nR,reg_len['N']/tot_cds).pvalue
locd={b['region']:b['n_start_sites_by_start_location'] for b in Brow}
isorate_N=locd['N']/reg_len['N']*1000; isorate_rest=sum(locd[r] for r in CDS if r!='N')/(tot_cds-reg_len['N'])*1000
Bst=dict(N_families=nN,N_start_isoforms=int(locd['N']),N_starts_per_family=Nn,N_fam_per_kb=round(famrate_N,3),rest_fam_per_kb=round(famrate_rest,3),
    fam_rate_ratio_N_over_rest=round(famrate_N/famrate_rest,3),fam_rate_binom_p=round(float(pbin),4),
    N_iso_per_kb=round(isorate_N,2),rest_iso_per_kb=round(isorate_rest,2),iso_rate_ratio=round(isorate_N/isorate_rest,3),
    starts_per_family_mean_N=round(float(np.mean(Nn)),2),starts_per_family_mean_rest=round(float(np.mean(rest)),2),
    starts_per_family_median_N=float(np.median(Nn)),starts_per_family_median_rest=float(np.median(rest)),
    mannwhitney_p=round(float(mw.pvalue),4),permutation_p_mean_diff=round(perm_p,4),n_rest_families=len(rest))
json.dump(Bst,open(W+"panelB_stats.json","w"),indent=1)
print(json.dumps(Bst,indent=1))
# ---------- featured summary + Panel C
cls=['exact_codon_retention','AUG_retention','near_cognate_retention','start_class_switching','loss_of_eligible_start','not_assessable']
groups=['GI-b','GI-a','GII-a','GII-b','GII-c','ALL','N1_historical_G1b_label']
fs=[];pc=[]
for fid in FEAT:
    F=fm[fid]; st=[x for x in F['isos'] if x['aa']>40]; f0=F['isos'][0]
    d=dict(family_id=fid,figure6_name_label=NAME[fid],start_positions_gt40='|'.join(str(x['start']) for x in st),start_codons='|'.join(x['start_codon'] for x in st),
           stop_position=F['stop'],max_aa=F['maxaa'],parent=f0['parent'],category=f0['category'],n_starts_gt40=len(st))
    labs=[]
    for g in groups:
        r=L[(fid,g)]; d[f'{g}_n_assessable']=r['n_assessable']; d[f'{g}_n_retained']=r['n_family_retained']; d[f'{g}_fraction']=r['fraction_retained']
        d[f'{g}_wilson95']=f"{r['wilson95_lo']}-{r['wilson95_hi']}"; d[f'{g}_fraction_dedup']=r['fraction_retained_dedup']; d[f'{g}_n_clusters']=r['n_clusters_assessable_dedup']
    gb=L[(fid,'GI-b')]
    fr=lambda r: float(r['fraction_retained']) if r['fraction_retained']!='' and int(r['n_assessable'])>=10 else None
    others={g:fr(L[(fid,g)]) for g in ['GI-a','GII-a','GII-b','GII-c']}
    if fr(gb) is None: labs.append('insufficient sample size in GI-b (n_assessable<10)')
    else:
        if fr(gb)>=0.9: labs.append('highly conserved in GI-b (>=0.90)')
        vals=[v for v in others.values() if v is not None]
        if fr(gb)>=0.9 and vals and all(v>=0.9 for v in vals): labs.append('highly conserved across all PEDV lineages with n>=10')
        if vals and any(fr(gb)-v>=0.10 for v in vals): labs.append('GI-b completeness higher than >=1 other lineage (diff>=0.10)')
        if fr(gb)<0.9: labs.append('not highly conserved in GI-b')
    d['conservation_statement_rule_based']='; '.join(labs)
    n=nullr.get(fid)
    d['synonymous_null_observed']=n['observed_fraction'] if n else ''; d['synonymous_null_expected']=n['expected_fraction_synonymous_null'] if n else ''; d['synonymous_null_z']=(n['z'] if float(n['expected_fraction_synonymous_null'])<0.999 else 'NA(degenerate: expected~1)') if n else ''
    d['functional_interpretation']='none (sequence retention only; hypothesis to be tested)'
    fs.append(d)
    for g in ['GI-b','ALL']:
        r=L[(fid,g)]; tot=sum(int(r['isoform_'+c]) for c in cls)
        row=dict(family_id=fid,name=NAME[fid],group=g,n_isoform_genome_pairs=tot)
        for c in cls: row[c]=r['isoform_'+c]; row[c+'_frac']=round(int(r['isoform_'+c])/tot,4) if tot else ''
        pc.append(row)
wr("PEDV_featured_family_GIb_summary.csv",fs)
wr("SourceData_PanelC_start_class_counts.csv",pc)
# ---------- background controls
ctl=[]
def frac(fid,g): r=L[(fid,g)]; return float(r['fraction_retained']) if r['n_assessable'] not in('0',0) and r['fraction_retained']!='' else None
for fid in FEAT:
    F=fm[fid]; f0=F['isos'][0]
    same_parent=[G for G in NC if G['fid']!=fid and G['fid'] not in FEAT and G['isos'][0]['parent']==f0['parent']]
    matched=[G for G in NC if G['fid']!=fid and G['fid'] not in FEAT and G['isos'][0]['parent']==f0['parent'] and G['isos'][0]['relative_frame']==f0['relative_frame']
             and abs(G['maxaa']-F['maxaa'])<=max(10,0.25*F['maxaa'])]
    for lab,S in (('same_parent_non_featured',same_parent),('region_length_frame_matched',matched)):
        for g in ['GI-b','ALL']:
            v=[frac(G['fid'],g) for G in S]; v=[x for x in v if x is not None]
            ctl.append(dict(featured_family=fid,figure6_name=NAME[fid],control_set=lab,lineage_group=g,n_control_families=len(v),
                control_median_fraction=round(float(np.median(v)),4) if v else '',control_q25=round(float(np.percentile(v,25)),4) if v else '',control_q75=round(float(np.percentile(v,75)),4) if v else '',
                featured_fraction=frac(fid,g),featured_percentile_rank_in_control=round(float(np.mean([x<=frac(fid,g) for x in v])),3) if v and frac(fid,g) is not None else ''))
    n=nullr.get(fid)
    ctl.append(dict(featured_family=fid,figure6_name=NAME[fid],control_set='synonymous_substitution_null(parent-CDS protein held fixed)',lineage_group='ALL',
        n_control_families='',control_median_fraction=n['expected_fraction_synonymous_null'] if n else '',featured_fraction=n['observed_fraction'] if n else '',
        featured_percentile_rank_in_control=(('z='+n['z']) if float(n['expected_fraction_synonymous_null'])<0.999 else 'z=NA(degenerate null)') if n else 'not applicable'))
wr("PEDV_ORF_background_controls.csv",ctl,['featured_family','figure6_name','control_set','lineage_group','n_control_families','control_median_fraction','control_q25','control_q75','featured_fraction','featured_percentile_rank_in_control'])
# ---------- Panel D matrix
D=[]
order=[F['fid'] for F in sorted(G40,key=lambda F:F['isos'][0]['start'])]
for g in ['GI-a','GI-b','GII-a','GII-b','GII-c','unresolved','ALL','N1_historical_G1b_label']:
    for fid in order:
        r=L[(fid,g)]; D.append(dict(lineage_group=g,family_id=fid,n_assessable=r['n_assessable'],fraction_retained=r['fraction_retained'],evidence_flag=r['evidence_flag']))
wr("SourceData_PanelD_lineage_family_matrix.csv",D)
# ---------- genome database QC
qc=list(csv.DictReader(open(W+"genome_map_qc.tsv",encoding='utf-8'),delimiter='\t'))
cl={r['accession']:r for r in csv.DictReader(open(W+"dedup_clusters.tsv",encoding='utf-8'),delimiter='\t')}
f1=set()
for l in open(B+"N1_conservation_analysis/data_external/PEDV.fasta",encoding='utf-8'):
    if l[0]=='>': f1.add(l[1:].split()[0])
comp={r['accession']:r for r in csv.DictReader(open(O+"SourceData_PanelD_per_genome_completeness.csv",encoding='utf-8-sig'))}
Q=[]
for r in qc:
    a=r['accession']; t=typing[a]
    inc=int(r['ref_pos_outside_genome'])
    Q.append(dict(accession_as_in_1083_file=a,accession_version_normalized=a if '.' in a else a+'.1',in_PEDV_fasta_1041=('yes' if (a in f1 or a+'.1' in f1) else 'no'),
        accession_unique='yes',pedv_by_homology=('yes' if float(r['identity_mapped_sampled'])>=0.9 else 'no'),mapped_identity_to_refdata_GIb_v1_1=r['identity_mapped_sampled'],
        genome_length=r['length'],ref_positions_outside_genome_5p3p=inc,end_completeness=('end_incomplete(>=30 ref nt not covered)' if inc>=30 else 'complete_ends'),
        ref_positions_deleted_internal=r['ref_pos_deleted'],non_ACGT_bases=r['non_ACGT'],N_bases=r['N'],
        non_ACGT_fraction=round(int(r['non_ACGT'])/int(r['length']),6),near_duplicate_cluster_le0p1pct=cl[a]['cluster'],cluster_representative=cl[a]['is_cluster_representative'],
        strain=t['strain'],country=t['country'],collection_date=t['collection_date'],
        S_type_this_study=t['S_type'],S_type_basis=t['S_type_basis'],N1_historical_genotype_label=t['N1_genotype_historical'],N1_status=t['N1_status'],
        n_gt40_families_assessable=comp[a]['n_assessable'] if a in comp else '',
        decision='included_main_analysis',
        reason='PEDV by homology; ends/ambiguity handled as not-assessable per family, never as loss' if inc<30 else 'included; 3/5-prime incomplete -> families overlapping uncovered ends are not-assessable (excluded from denominators)'))
wr("PEDV_genome_database_QC.csv",Q)
print('QC rows',len(Q),collections.Counter(q['end_completeness'] for q in Q))
