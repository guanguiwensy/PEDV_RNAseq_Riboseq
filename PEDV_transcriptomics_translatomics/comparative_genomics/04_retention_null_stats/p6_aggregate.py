"""Aggregation of sequence-only results into the 8 source-data tables + panel data. Figure 6 columns are attached LAST and only as labels."""
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
T=PEDV_P+r"global_perturbation_screen/plus1_orf_capacity/"
pool=json.load(open(W+"refpool.json"))['fams']; fm={F['fid']:F for F in pool}
LIN=['GI-a','GI-b','GII-a','GII-b','GII-c','unresolved']
typing={r['accession']:r for r in csv.DictReader(open(W+"typing_S_based.tsv",encoding='utf-8'),delimiter='\t')}
cl={}
for r in csv.DictReader(open(W+"dedup_clusters.tsv",encoding='utf-8'),delimiter='\t'): cl[r['accession']]=int(r['cluster'])
rows=list(csv.DictReader(open(O+"PEDV_family_homology_mapping.csv",encoding='utf-8-sig')))
def wilson(k,n,z=1.96):
    if n==0: return ('','')
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return (round(max(0,c-h),4),round(min(1,c+h),4))
def wr(fn,rs,cols=None):
    cols=cols or list(rs[0])
    with open(O+fn,"w",encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols); w.writeheader(); w.writerows(rs)
CSTOPS={12620,20611,24756,25031,25598,26286,27623}
def stat(r):
    t=fm[r['family_id']]['tier']
    if t=='gt40':
        s=r['hom_status_gt40']; return {'NA':'not_assessable'}.get(s,s)   # homologous-start definition (primary)
    return r['status_ge8']
def stat_any(r):
    return r['status_gt40'] if fm[r['family_id']]['tier']=='gt40' else r['status_ge8']
groups={'ALL':lambda a:True,**{l:(lambda a,l=l:typing[a]['S_type']==l) for l in LIN},
        'N1_historical_G1b_label':lambda a:typing[a]['N1_genotype_historical']=='GIb'}
# index
by=collections.defaultdict(list)
for r in rows: by[r['family_id']].append(r)
cls=['exact_codon_retention','AUG_retention','near_cognate_retention','start_class_switching','loss_of_eligible_start','not_assessable']
out=[]; ret={}   # ret[(fid,group)] = (k,n,k_dedup,n_dedup)
for fid,rs in by.items():
    tier=fm[fid]['tier']
    for g,f in groups.items():
        sub=[r for r in rs if f(r['accession'])]
        c=collections.Counter(stat(r) for r in sub)
        k=c['retained']; lost=c['lost']; n=k+lost
        # cluster-level (majority among assessable members; clusters with no assessable member dropped)
        cc=collections.defaultdict(lambda:[0,0])
        for r in sub:
            s=stat(r)
            if s in('retained','lost'): cc[cl[r['accession']]][0 if s=='retained' else 1]+=1
        kd=sum(1 for a,b in cc.values() if a>b); nd=sum(1 for a,b in cc.values() if a!=b)
        nd_all=len(cc); kd=sum(1 for a,b in cc.values() if a>=b)   # ties counted retained
        lo,hi=wilson(k,n)
        c8=collections.Counter(r['status_ge8'] for r in sub); k8=c8['retained']; n8=k8+c8['lost']
        ca=collections.Counter(stat_any(r) for r in sub); ka=ca['retained']; na=ka+ca['lost']
        ret[(fid,g)]=(k,n,kd,nd_all)
        d=dict(family_id=fid,length_tier=tier,lineage_group=g,n_genomes=len(sub),n_not_assessable=c['not_assessable'],n_unresolved=c['unresolved'],
               n_assessable=n,n_family_retained=k,n_family_lost=lost,fraction_retained=(round(k/n,4) if n else ''),wilson95_lo=lo,wilson95_hi=hi,
               n_assessable_anyStart=na,n_retained_anyStart=ka,fraction_retained_anyStart_sensitivity=(round(ka/na,4) if na else ''),
               n_assessable_ge8=n8,n_retained_coding_potential_ge8=k8,fraction_retained_ge8=(round(k8/n8,4) if n8 else ''),
               evidence_flag=('insufficient_n(<10 assessable)' if n<10 else 'ok'),
               n_clusters_assessable_dedup=nd_all,n_clusters_retained_dedup=kd,fraction_retained_dedup=(round(kd/nd_all,4) if nd_all else ''))
        for q in cls: d['isoform_'+q]=sum(int(r[q+'_gt40']) for r in sub)
        out.append(d)
wr("PEDV_ORF_conservation_by_lineage.csv",out)
print('lineage table rows',len(out))
# ---- per-genome completeness
comp=collections.defaultdict(lambda:[0,0])
for r in rows:
    if fm[r['family_id']]['tier']!='gt40' or fm[r['family_id']]['stop'] in CSTOPS: continue
    s=stat(r)
    if s in('retained','lost'): comp[r['accession']][0 if s=='retained' else 1]+=1
pg=[dict(accession=a,S_type=typing[a]['S_type'],cluster=cl[a],n_gt40_families_retained=v[0],n_lost=v[1],n_assessable=v[0]+v[1],fraction=round(v[0]/(v[0]+v[1]),4)) for a,v in comp.items()]
wr("SourceData_PanelD_per_genome_completeness.csv",pg)
# ---- Figure 6 labels and crosswalk (labels attached only after all sequence results)
fam6=list(csv.DictReader((l for l in open(T+"tables/plus1_ORF_families_gt40aa.tsv",encoding='utf-8') if not l.startswith('#')),delimiter='\t'))
name6={r['original_family_id']:r['displayed_name'] for r in csv.DictReader(open(T+"panelC_family_name_mapping.tsv",encoding='utf-8'),delimiter='\t')}
top5=[r['displayed_name'] for r in csv.DictReader(open(T+"PanelC_top5_family_selection.tsv",encoding='utf-8'),delimiter='\t')]
FEATURED_NAMES=['ORF1a-x3','ORF1a-x4','ORF1a-x22','S-x1','S-x2','N-x1','N-x2']   # examined in Panel C (user-specified)
FIG6_ORIG_AC=['ORF1a-x22','N-x1','N-x2']; FIG6_ORIG_D=['ORF1a-x3','ORF1b-x13','ORF1b-x8']
FIG6_REV_RULE=['N-x2','ORF1a-x9','ORF1a-x23','S-x1','ORF1a-x3','S-x2','N-x5']
CST={12620:'ORF1a',20611:'ORF1b',24756:'S',25031:'ORF3',25598:'E',26286:'M',27623:'N'}
pc={'ORF1a':293,'ORF1b':12592,'S':20610,'ORF3':24758,'E':25361,'M':25608,'N':26300}
cw=[]; used=set(); match={}
bystop=collections.defaultdict(list)
for F in pool: bystop[F['stop']].append(F)
for r in fam6:
    fid6=r['family_id']; par=r['parent_ORF']; stop=int(r['shared_stop']); st6=[int(x) for x in r['all_candidate_starts'].split('|')]
    frame=((pc[par]-1)%3+1)%3
    cand=[F for F in bystop.get(stop,[]) if F['frame']==frame]
    if len(cand)==1:
        F=cand[0]; s_new=[x['start'] for x in F['isos'] if x['aa']>40]
        ok=sorted(s_new)==sorted(st6); used.add(F['fid']); match[fid6]=F['fid']
        cw.append(dict(figure6_family_id=fid6,figure6_displayed_name=name6.get(fid6,''),figure6_parent=par,figure6_shared_stop=stop,figure6_n_start_variants=len(st6),
            sequence_defined_family_id=F['fid'],one_to_one='yes' if ok else 'no',n_start_variants_sequence_scan_gt40=len(s_new),
            n_start_variants_sequence_scan_ge8=F['n'],reason_if_not_one_to_one='' if ok else ('Figure 6 enumerated only starts inside parent %s; the sequence-defined family shares the stop codon of annotated %s (in-frame N-terminal extension reaching into the parent), %d vs %d start variants'%(par,CST.get(F['stop'],'-'),len(s_new),len(st6)) if F['stop'] in CSTOPS else 'start sets differ'),
            fig6_max_D=r['max_D'],fig6_n_positive_starts=r['n_positive_starts'],fig6_n_cross_time_supported_starts=r['n_cross_time_supported_starts'],
            fig6_best_local_percentile=r['best_local_percentile'],fig6_HAR_response_architecture=r['HAR_response_architecture'],fig6_detectability_status=r['detectability_status']))
    else:
        cw.append(dict(figure6_family_id=fid6,figure6_displayed_name=name6.get(fid6,''),figure6_parent=par,figure6_shared_stop=stop,figure6_n_start_variants=len(st6),
            sequence_defined_family_id='UNMATCHED',one_to_one='no',reason_if_not_one_to_one='no family with this stop/frame in sequence scan'))
fcols=list(cw[0])
for F in pool:
    if F['maxaa']>40 and F['fid'] not in used:
        f0=F['isos'][0]
        why=[]
        if f0['category']!='internal alternative-frame ORF': why.append('category='+f0['category']+' (not an internal alternative-frame ORF)')
        if f0['relative_frame'] not in('+1',): why.append('relative frame to parent '+(f0['relative_frame'] or 'n/a')+' (Figure 6 analysed +1 only)')
        if f0['parent'] in('-',''): why.append('start outside every parent CDS (Figure 6 enumerated starts inside parents only)')
        cw.append(dict(figure6_family_id='NONE',figure6_displayed_name='',figure6_parent='',figure6_shared_stop='',figure6_n_start_variants='',
            sequence_defined_family_id=F['fid'],one_to_one='sequence-only family',n_start_variants_sequence_scan_gt40=sum(1 for x in F['isos'] if x['aa']>40),
            n_start_variants_sequence_scan_ge8=F['n'],reason_if_not_one_to_one='; '.join(why) or 'unclassified'))
wr("PEDV_figure6_family_crosswalk.csv",cw,fcols)
print('crosswalk: matched',len(used),'one_to_one yes',sum(1 for c in cw if c['one_to_one']=='yes'),'seq-only',sum(1 for c in cw if c['one_to_one']=='sequence-only family'))
featured=[match[k] for k,v in name6.items() if v in FEATURED_NAMES and k in match]
inv={v:k for k,v in name6.items()}
json.dump(dict(featured=featured,names={match[k]:name6[k] for k in match},fig6=match,orig_AC=[match[inv[n]] for n in FIG6_ORIG_AC if inv[n] in match],orig_D=[match[inv[n]] for n in FIG6_ORIG_D if inv.get(n) in match],rev_rule=[match[inv[n]] for n in FIG6_REV_RULE if inv[n] in match]),open(W+"featured.json","w"))
print('featured',featured)
# ---- Figure 6 label vs architecture retention (correlation only)
lab={r['family_id']:r for r in fam6}
xs=[];Dm=[];Pc=[];Cs=[];ids_=[]
for k,fid in match.items():
    r=lab[k]; k_,n_,kd,nd=ret[(fid,'ALL')]
    if n_<10: continue
    ids_.append(fid); xs.append(k_/n_)
    Dm.append(float(r['max_D']) if r['max_D'] not in('',) else np.nan); Pc.append(float(r['best_local_percentile']) if r['best_local_percentile'] not in('',) else np.nan)
    Cs.append(1 if int(r['n_cross_time_supported_starts'])>=1 else 0)
xs=np.array(xs);Dm=np.array(Dm);Pc=np.array(Pc);Cs=np.array(Cs)
res={}
for nm,v in (('max_D',Dm),('best_local_percentile',Pc)):
    ok=~np.isnan(v); rho,p=stats.spearmanr(v[ok],xs[ok]); res[nm]=dict(rho=round(float(rho),3),p=round(float(p),4),n=int(ok.sum()))
u=stats.mannwhitneyu(xs[Cs==1],xs[Cs==0],alternative='two-sided')
res['cross_time_supported_vs_not']=dict(n1=int((Cs==1).sum()),n0=int((Cs==0).sum()),median1=round(float(np.median(xs[Cs==1])),4),median0=round(float(np.median(xs[Cs==0])),4),U=float(u.statistic),p=round(float(u.pvalue),4))
# same on GI-b only
xb=np.array([ret[(f,'GI-b')][0]/ret[(f,'GI-b')][1] if ret[(f,'GI-b')][1]>=10 else np.nan for f in ids_])
ok=~np.isnan(xb)
res['GIb_only_cross_time_supported_vs_not']=dict(n1=int(((Cs==1)&ok).sum()),n0=int(((Cs==0)&ok).sum()),median1=round(float(np.nanmedian(xb[(Cs==1)&ok])),4),median0=round(float(np.nanmedian(xb[(Cs==0)&ok])),4),
    p=round(float(stats.mannwhitneyu(xb[(Cs==1)&ok],xb[(Cs==0)&ok]).pvalue),4))
res['frac_families_retained_ge_0.99_ALL']=dict(n=len(xs),n_ge_0_99=int((xs>=0.99).sum()))
json.dump(res,open(W+"fig6_label_vs_retention.json","w"),indent=1)
print(json.dumps(res,indent=1))
