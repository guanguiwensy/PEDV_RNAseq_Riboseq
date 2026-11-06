"""Project the 1,083 genomes onto reference coordinates using the Sup L homology map (refmaps.npz; read-only).
state codes: 0-3 = A,C,G,T ; 4 = ambiguous/IUPAC ; 5 = internal deletion (-1) ; 6 = outside genome (-2)."""
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
from r6lib import *
z=np.load(SL+'work/refmaps.npz'); ids=list(z['ids']); M=z['M']
gen={};cur=None
for l in open(B+'N1_conservation_analysis/results/PEDV_genomes_all.fasta',encoding='utf-8'):
    if l[0]=='>': cur=l[1:].split()[0]; gen[cur]=[]
    else: gen[cur].append(l.strip().upper())
P=np.zeros((len(ids),LREF),np.uint8); 
for gi,a in enumerate(ids):
    e=enc[np.frombuffer(''.join(gen[a]).encode(),np.uint8)]; m=M[gi]
    ok=m>=0; P[gi,ok]=e[m[ok]]; P[gi,m==-1]=5; P[gi,m==-2]=6
# contiguity of consecutive reference positions (False at indel junctions)
contig=np.zeros((len(ids),LREF),bool); contig[:,:-1]=(M[:,1:]==M[:,:-1]+1)&(M[:,:-1]>=0)
np.savez_compressed(OUT+'work/proj.npz',ids=np.array(ids),P=P,contig=contig)
cl=pd.read_csv(SL+'work/dedup_clusters.tsv',sep='\t'); print(cl.is_cluster_representative.value_counts().to_dict(), cl.cluster.nunique())
print((P==5).mean(),(P==6).mean(),(P==4).mean())
