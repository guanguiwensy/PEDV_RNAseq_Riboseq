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
import sys,time,numpy as np,csv,hashlib
sys.path[0:0] = [REPO_ROOT + '/comparative_genomics/02_alignment_homology_lineage', REPO_ROOT + '/comparative_genomics/03_orf_scan', REPO_ROOT + '/comparative_genomics/04_retention_null_stats']
from mapper import *
B=PEDV_K+r""; W=B+"Supplementary_Figure_L/work/"
ref="".join(l.strip() for l in open(PEDV_P+r"refbuild_GIb/v1.1/virus/PEDV_GIb.fa") if l[0]!='>').upper()
gen={};cur=None
for l in open(B+"N1_conservation_analysis/results/PEDV_genomes_all.fasta",encoding='utf-8'):
    if l[0]=='>': cur=l[1:].split()[0]; gen[cur]=[]
    else: gen[cur].append(l.strip().upper())
gen={k:''.join(v) for k,v in gen.items()}
idx=build_index(ref)
ids=list(gen); M=np.zeros((len(ids),len(ref)),dtype=np.int32); rows=[]
for n,a in enumerate(ids):
    m,i=map_genome(ref,idx,gen[a]); M[n]=m
    g=gen[a]; nonacgt=sum(1 for c in g if c not in 'ACGT')
    rows.append([a,len(g),round(i,4),int((m==-1).sum()),int((m==-2).sum()),nonacgt,g.count('N')])
np.savez_compressed(W+"refmaps.npz",ids=np.array(ids),M=M)
with open(W+"genome_map_qc.tsv","w",encoding='utf-8',newline='') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['accession','length','identity_mapped_sampled','ref_pos_deleted','ref_pos_outside_genome','non_ACGT','N']); w.writerows(rows)
import statistics
print('n',len(ids),'min ident',min(r[2] for r in rows),'median',statistics.median(r[2] for r in rows))
print('low identity',[r for r in rows if r[2]<0.93][:10])
print('out>200',sum(1 for r in rows if r[4]>200),'del>100',sum(1 for r in rows if r[3]>100))
