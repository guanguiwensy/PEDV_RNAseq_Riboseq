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
import sys,csv,json,collections,numpy as np,warnings
warnings.filterwarnings('ignore')
sys.path[0:0] = [REPO_ROOT + '/comparative_genomics/02_alignment_homology_lineage', REPO_ROOT + '/comparative_genomics/03_orf_scan', REPO_ROOT + '/comparative_genomics/04_retention_null_stats']
from scan_core import *
B=PEDV_K+r""; O=B+"Supplementary_Figure_L/"; W=O+"work/"
ref="".join(l.strip() for l in open(PEDV_P+r"refbuild_GIb/v1.1/virus/PEDV_GIb.fa") if l[0]!='>').upper()
z=np.load(W+"refmaps.npz"); ids=list(z['ids']); M=z['M']
gen={};cur=None
for l in open(B+"N1_conservation_analysis/results/PEDV_genomes_all.fasta",encoding='utf-8'):
    if l[0]=='>': cur=l[1:].split()[0]; gen[cur]=[]
    else: gen[cur].append(l.strip().upper())
gen={k:''.join(v) for k,v in gen.items()}
typing={r['accession']:r['S_type'] for r in csv.DictReader(open(W+"typing_S_based.tsv",encoding='utf-8'),delimiter='\t')}
pool=json.load(open(W+"refpool.json"))['fams']
refstops={F['stop'] for F in pool}; refstops_gt40={F['stop'] for F in pool if F['maxaa']>40}
agg=collections.defaultdict(lambda: collections.Counter())
per=[]
for i,a in enumerate(ids):
    g=gen[a]; inv=np.full(len(g)+3,-1,np.int32)
    m=M[i]; ok=m>=0; inv[m[ok]]=np.nonzero(ok)[0]
    comp,trunc=scan_seq(g,min_aa=8)
    fam=collections.defaultdict(list)
    for r in comp: fam[(r['frame'],r['stop'])].append(r)
    n40=0;match=0;gain=0;ctx=0;unm=0
    for (f,s),v in fam.items():
        mx=max(x['aa'] for x in v)
        if mx<=40: continue
        n40+=1
        s0=s-1
        r0=inv[s0]; r1=inv[s0+1]; r2=inv[s0+2]
        if r0<0 or r1!=r0+1 or r2!=r0+2: cat='unmappable_or_inserted';unm+=1; key='NA'
        else:
            key=int(r0)+1
            if key in refstops: cat='matches_reference_family';match+=1
            elif ref[r0:r0+3] in STOPS: cat='start_gained_at_reference_stop';gain+=1
            else: cat='new_stop_context_not_reference_stop';ctx+=1
        agg[(key,cat)][typing[a]]+=1; agg[(key,cat)]['ALL']+=1
    per.append(dict(accession=a,S_type=typing[a],n_families_gt40_independent=n40,n_match_reference_family=match,n_start_gained_at_reference_stop=gain,
                    n_new_stop_context=ctx,n_unmappable_or_inserted=unm,n_boundary_truncated_starts=len(trunc)))
with open(O+"PEDV_independent_scan_per_genome.csv","w",encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(per[0])); w.writeheader(); w.writerows(per)
lin=['GI-a','GI-b','GII-a','GII-b','GII-c','unresolved','ALL']
rows=[]
for (key,cat),c in sorted(agg.items(),key=lambda kv:(str(kv[0][1]),kv[0][0] if isinstance(kv[0][0],int) else 10**9)):
    rows.append(dict(reference_stop_coordinate=key,category=cat,**{f"n_genomes_{l}":c.get(l,0) for l in lin}))
with open(O+"PEDV_independent_scan_families_by_reference_coordinate.csv","w",encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print(collections.Counter(r['category'] for r in rows))
import statistics
print('median families>40 per genome',statistics.median(p['n_families_gt40_independent'] for p in per),'min/max',min(p['n_families_gt40_independent'] for p in per),max(p['n_families_gt40_independent'] for p in per))
nov=[r for r in rows if r['category']!='matches_reference_family' and r['n_genomes_ALL']>=20]
print(len(nov)); 
for r in nov[:20]: print(r)
