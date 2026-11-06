"""Build Supplementary Figure 16 C1/C2 data for the configured families."""
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
import csv,json,collections
B=PEDV_K+r""; O=B+"Supplementary_Figure_L/"; W=O+"work/"
T=PEDV_P+r""
pool=json.load(open(W+"refpool.json"))['fams']; fm={F['fid']:F for F in pool}
feat=json.load(open(W+"featured.json"))
fig6=feat['fig6']                         # P1 id -> SDF id
name6={r['original_family_id']:r['displayed_name'] for r in csv.DictReader(open(T+"global_perturbation_screen/plus1_orf_capacity/panelC_family_name_mapping.tsv",encoding='utf-8'),delimiter='\t')}
order=json.load(open(REPO_ROOT+'/config/sup16_featured_families.json',encoding='utf-8'))['featured_family_ids']
assert set(order) <= set(fig6.values()), "configured featured families are absent from the Figure 6 crosswalk"
feat['C12_list']=order
feat['featured']=list(order)
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
    selection_status="configured") for i,f in enumerate(order)])
