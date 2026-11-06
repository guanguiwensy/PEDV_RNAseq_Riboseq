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
import sys,csv,collections,numpy as np,warnings
warnings.filterwarnings('ignore')
sys.path[0:0] = [REPO_ROOT + '/comparative_genomics/02_alignment_homology_lineage', REPO_ROOT + '/comparative_genomics/03_orf_scan', REPO_ROOT + '/comparative_genomics/04_retention_null_stats']
from mapper import *
from Bio import Phylo
B=PEDV_K+r""; W=B+"Supplementary_Figure_L/work/"; P=PEDV_P+r"pedv_s_gene_phylogeny/"
ref="".join(l.strip() for l in open(PEDV_P+r"refbuild_GIb/v1.1/virus/PEDV_GIb.fa") if l[0]!='>').upper()
z=np.load(W+"refmaps.npz"); ids=list(z['ids']); M=z['M']
gen={};cur=None
for l in open(B+"N1_conservation_analysis/results/PEDV_genomes_all.fasta",encoding='utf-8'):
    if l[0]=='>': cur=l[1:].split()[0]; gen[cur]=[]
    else: gen[cur].append(l.strip().upper())
gen={k:''.join(v) for k,v in gen.items()}
enc=np.full(256,4,np.uint8)
for i,c in enumerate('ACGT'): enc[ord(c)]=i
def proj(g,m):
    a=enc[np.frombuffer(g.encode(),np.uint8)]
    out=np.full(len(m),5,np.uint8)      # 5 = gap(deleted) ; 6 = outside
    ok=m>=0; out[ok]=a[m[ok]]; out[m==-2]=6; return out   # 4 = ambiguous
Gp=np.stack([proj(gen[a],M[i]) for i,a in enumerate(ids)])
# ---- S projection
S0,S1=20609,24758
Sref=ref[S0:S1]; sidx=build_index(Sref)
# panel members
al={};cur=None
for l in open(P+"alignment/S_gene_alignment.fasta"):
    if l[0]=='>': cur=l[1:].split()[0]; al[cur]=[]
    else: al[cur].append(l.strip().upper())
al={k:''.join(v).replace('-','') for k,v in al.items()}
# labels from SGP anchors via tree
lab={r['accession']:r['label_Guo2019'] for r in csv.DictReader(open(P+"intermediate/panel_curated_labels.tsv",encoding='utf-8'),delimiter='\t') if r['label_Guo2019'].startswith('G')}
t=Phylo.read(P+"tree/S_main.treefile","newick")
leaves=[c.name for c in t.get_terminals()]
anch=[a for a in lab if a in leaves]
plabel={}
for lf in leaves:
    if lf not in al or lf.startswith(('CV777_','LZ202401_')): continue
    if lf in lab: plabel[lf]=lab[lf]; continue
    best=min(anch,key=lambda a:t.distance(lf,a)); plabel[lf]=lab[best]
print(collections.Counter(plabel.values()),len(plabel))
pnames=list(plabel)
def sproj(s):
    m,_=map_genome(Sref,sidx,s); return proj(s,m)
Pp=np.stack([sproj(al[n]) for n in pnames])
Sg=Gp[:,S0:S1]
def pdist(a,B):
    ok=(a<4)[None,:]&(B<4); mm=(a[None,:]!=B)&ok
    return mm.sum(1)/np.maximum(ok.sum(1),1),ok.sum(1)
rows=[]
for i,a in enumerate(ids):
    d,n=pdist(Sg[i],Pp)
    # exclude self (panel member identical accession)
    for j,nm in enumerate(pnames):
        if nm==a or nm+'.1'==a or nm.split('.')[0]==a.split('.')[0]: d[j]=9
    cov=(Sg[i]<4).sum()/ (S1-S0)
    ok=n>=3400
    d=np.where(ok,d,9)
    j=int(np.argmin(d)); d1=float(d[j])
    if d1>=9: call='unresolved';why='low_S_overlap'
    elif d1>0.03: call='unresolved';why='nearest_panel_distance>0.03'
    else:
        near={plabel[pnames[k]] for k in range(len(pnames)) if d[k]<=d1+0.003}
        if len(near)==1: call=near.pop();why='nearest_panel_member'
        else: call='unresolved';why='ambiguous_within_0.003:'+'/'.join(sorted(near))
    rows.append((a,call,why,pnames[j],round(d1,4),round(cov,3)))
gn={r['accession']:r for r in csv.DictReader(open(B+"N1_conservation_analysis/results/genotype_normalized.tsv",encoding='utf-8'),delimiter='\t')}
meta={}
for r in csv.DictReader(open(B+"N1_conservation_analysis/data_external/PEDV_metadata.tsv",encoding='utf-8'),delimiter='\t'):
    meta[r['accession']+'.'+r['version']]=r; meta[r['accession']]=r
rows2=[]
for r in rows:
    g=gn[r[0]]; m=meta.get(r[0]) or meta.get(r[0].split('.')[0]) or {}
    rows2.append(list(r)+[g['genotype_normalized'],g['genotype_status'],m.get('country',''),m.get('collection_date',''),m.get('strain','')])
with open(W+"typing_S_based.tsv","w",encoding='utf-8',newline='') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['accession','S_type','S_type_basis','nearest_panel_member','nearest_p_distance','S_coverage','N1_genotype_historical','N1_status','country','collection_date','strain']); w.writerows(rows2)
c=collections.Counter((r[6],r[1]) for r in rows2)
for k,v in sorted(c.items()): print(k,v)
print(collections.Counter(r[1] for r in rows2))
# ---- dedup clustering on whole-genome projection
valid=(Gp<4)
cl=-np.ones(len(ids),int); reps=[]
order=np.argsort([-int(valid[i].sum()) for i in range(len(ids))])   # most complete first
for i in order:
    if cl[i]>=0: continue
    k=len(reps); reps.append(i); cl[i]=k
    rest=np.nonzero(cl<0)[0]
    if len(rest)==0: break
    ok=valid[rest]&valid[i][None,:]
    mm=((Gp[rest]!=Gp[i][None,:])&ok).sum(1); n=ok.sum(1)
    close=(n>=20000)&(mm/np.maximum(n,1)<=0.001)
    cl[rest[close]]=k
print('clusters at <=0.1% divergence:',len(reps))
np.save(W+"cluster_id.npy",cl)
with open(W+"dedup_clusters.tsv","w",encoding='utf-8',newline='') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['accession','cluster','is_cluster_representative'])
    for i,a in enumerate(ids): w.writerow([a,int(cl[i]),'yes' if reps[cl[i]]==i else 'no'])
