"""New C1/C2 family list (provisional, from Figure 6 checkpoint-3 evidence architecture labels - NOT from conservation) + panel data."""
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
import csv,json,collections
B=PEDV_K+r""; O=B+"Supplementary_Figure_L/"; W=O+"work/"
T=PEDV_P+r""
pool=json.load(open(W+"refpool.json"))['fams']; fm={F['fid']:F for F in pool}
feat=json.load(open(W+"featured.json"))
if 'featured_old7' not in feat: feat['featured_old7']=feat['featured']
fig6=feat['fig6']                         # P1 id -> SDF id
name6={r['original_family_id']:r['displayed_name'] for r in csv.DictReader(open(T+"global_perturbation_screen/plus1_orf_capacity/panelC_family_name_mapping.tsv",encoding='utf-8'),delimiter='\t')}
v2=list(csv.DictReader(open(T+"SupM_alt_ORF_evidence/Fig6_checkpoint3/tables/Fig6_joint_evidence_classification_v2.tsv",encoding='utf-8-sig'),delimiter='\t'))
V={r['family']:r for r in v2}
SEL=[("evidence architecture 'CHX candidate-frame signal elevated above matched background' (Fig6 checkpoint 3 v2)",["N-x1"]),
     ("evidence architecture 'Complementary CHX/HAR evidence' (Fig6 checkpoint 3 v2); Figure 6D case",["ORF1a-x3"]),
     ("evidence architecture 'HAR differential-response informative' (Fig6 checkpoint 3 v2)",["N-x2","ORF1a-x9","ORF1a-x23","S-x1","S-x2","N-x5"]),
     ("alternative Figure 6D case proposed in checkpoint 3 notes (partial CHX evidence); provisional",["S-x4"])]
reason={}
for why,names in SEL:
    for n in names: reason[n]=why
# consistency check against v2 table
for n in reason:
    assert n in V, n
inv={v:k for k,v in name6.items()}
ids=[fig6[inv[n]] for n in reason]
order=sorted(ids,key=lambda f:min(x['start'] for x in fm[f]['isos']))   # genomic order
feat['C12_list']=order; feat['C12_reason']={f:reason[name6[[k for k,v in fig6.items() if v==f][0]]] for f in order}
import json as _js  # curation patch
_fixed=_js.load(open(REPO_ROOT+'/config/sup16_featured_families.json',encoding='utf-8'))['featured_family_ids']
assert order==_fixed, ('recomputed featured list differs from the fixed Sup16 config', order, _fixed)
feat['featured']=list(_fixed)  # identical to 'order'; fixed config is authoritative
json.dump(feat,open(W+"featured.json","w"))
names={f:name6[[k for k,v in fig6.items() if v==f][0]] for f in order}
print([names[f] for f in order])
# ---- data
sc=list(csv.DictReader(open(O+"SupL_start_codon_conservation.tsv",encoding='utf-8'),delimiter='\t'))
cls=['exact_codon_retention','AUG_retention','near_cognate_retention','start_class_switching','loss_of_eligible_start','not_assessable']
agg=collections.defaultdict(lambda:collections.Counter())
for r in sc:
    for c in cls:
        agg[(r['family_id'],r['lineage'])][c]+=int(r[c])
        agg[(r['family_id'],'ALL')][c]+=int(r[c])
lin=list(csv.DictReader(open(O+"PEDV_ORF_conservation_by_lineage.csv",encoding='utf-8-sig')))
L={(r['family_id'],r['lineage_group']):r for r in lin}
rows1=[];rows2=[]
for f in order:
    for g in ['GI-a','GI-b','GII-a','GII-b','GII-c','ALL']:
        c=agg[(f,g)]; tot=sum(c.values())
        d=dict(family=names[f],family_id=f,lineage=g,n_start_variants_gt40=sum(1 for x in fm[f]['isos'] if x['aa']>40),n_start_variant_genome_pairs=tot)
        for k in cls: d[k]=c[k]; d[k+'_fraction']=round(c[k]/tot,4) if tot else ''
        d['n_assessable_pairs']=tot-c['not_assessable']
        rows1.append(d)
        r=L[(f,g)]
        rows2.append(dict(family=names[f],family_id=f,lineage=g,n_genomes=r['n_genomes'],n_assessable=r['n_assessable'],n_not_assessable=r['n_not_assessable'],
            n_retained=r['n_family_retained'],n_lost=r['n_family_lost'],retention_fraction=r['fraction_retained'],wilson95_lo=r['wilson95_lo'],wilson95_hi=r['wilson95_hi'],evidence_flag=r['evidence_flag']))
def wr(fn,rows):
    with open(O+fn,"w",encoding='utf-8-sig',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
wr("SourceData_PanelC1_start_codon_by_lineage.csv",rows1)
wr("SourceData_PanelC2_family_retention_by_lineage.csv",rows2)
wr("SupL_C1C2_family_selection.csv",[dict(order=i+1,family=names[f],family_id=f,parent=fm[f]['isos'][0]['parent'],start_min=min(x['start'] for x in fm[f]['isos']),stop=fm[f]['stop'],
    selection_basis=feat['C12_reason'][f],Fig6_evidence_architecture_v2=V[names[f]]['evidence_architecture'],Fig6_CHX_status=V[names[f]]['CHX_status'],Fig6_HAR_status=V[names[f]]['HAR_status'],
    original_Figure6_featured=V[names[f]]['featured_fig6'],selection_status="provisional (Figure 6 not frozen)") for i,f in enumerate(order)])
