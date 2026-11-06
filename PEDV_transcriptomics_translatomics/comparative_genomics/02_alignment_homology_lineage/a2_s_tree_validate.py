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
import csv,re,collections,sys
sys.stdout.reconfigure(encoding='utf-8')
from Bio import Phylo
P=PEDV_P+r"pedv_s_gene_phylogeny/"
N=PEDV_K+r"N1_conservation_analysis/results/"
O=PEDV_K+r"Supplementary_Figure_L/audit/"
lab={r['accession']:r for r in csv.DictReader(open(P+"intermediate/panel_curated_labels.tsv",encoding='utf-8'),delimiter='\t')}
anch={a:r['label_Guo2019'] for a,r in lab.items() if r['label_Guo2019'].startswith('G')}
print(collections.Counter(anch.values()))
# cdhit clusters
clu={};rep={}
cur=None;members=collections.defaultdict(list)
for l in open(P+"intermediate/r1/refpanel/cdhit_S_0.985.fasta.clstr"):
    if l[0]=='>': cur=l.strip()[1:]; continue
    m=re.search(r">(\S+?)\.\.\.",l); a=m.group(1); members[cur].append(a)
    if l.rstrip().endswith('*'): rep[cur]=a
a2rep={a:rep[c] for c,ms in members.items() for a in ms}
allS=[r['accession'] for r in csv.DictReader(open(P+"intermediate/r1/refpanel/all_S_qc.tsv",encoding='utf-8'),delimiter='\t')]
print('allS',len(allS),'clustered',len(a2rep))
t=Phylo.read(P+"tree/S_main.treefile","newick")
leaves={c.name for c in t.get_terminals()}
print('tree leaves',len(leaves))
# nearest labelled anchor by patristic distance
anch_in=[a for a in anch if a in leaves]
print('anchors in tree',len(anch_in))
def nearest(leaf):
    best=None
    for a in anch_in:
        if a==leaf: return a,0.0
        d=t.distance(leaf,a)
        if best is None or d<best[1]: best=(a,d)
    return best
cache={}
def call(acc):
    r=a2rep.get(acc) or a2rep.get(acc+'.1')
    if r is None: return ('NA','not_in_S_panel_input','NA')
    if r not in leaves: return ('NA','rep_not_in_tree',r)
    if r not in cache: cache[r]=nearest(r)
    a,d=cache[r]
    return (anch[a],('direct' if r==acc else 'via_cdhit_rep'),f"{r}|{a}|{d:.4f}")
gn=list(csv.DictReader(open(N+"genotype_normalized.tsv",encoding='utf-8'),delimiter='\t'))
out=[]
for r in gn:
    acc=r['accession']; g,how,info=call(acc); out.append((acc,r['genotype_normalized'],r['genotype_status'],g,how,info))
with open(O+"s_tree_genotype_crosscheck.tsv","w",encoding='utf-8',newline='') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['accession','N1_genotype','N1_status','S_tree_nearest_anchor_label','method','rep|anchor|patristic']); w.writerows(out)
c=collections.Counter((o[1],o[3]) for o in out)
for k,v in sorted(c.items()): print(k,v)
print(collections.Counter(o[4] for o in out))
