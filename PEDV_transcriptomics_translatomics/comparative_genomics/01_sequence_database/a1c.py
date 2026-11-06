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
import sys,re,collections
sys.stdout.reconfigure(encoding='utf-8')
exec(open(__import__('os').path.join(__import__('os').path.dirname(__REPO_FILE__), 'a1_source_audit.py'),encoding='utf-8').read().split("print(\"fasta1\"")[0])
diff=[k for k in set(f1)&set(f2) if f1[k]!=f2[k]]
ok=sum(1 for k in diff if re.sub('[^ACGT]','',f2[k])==f1[k]); print('f1 == f2 with non-ACGT deleted:',ok,'/',len(diff))
tot=collections.Counter(); 
for k,s in f2.items():
    n=sum(1 for c in s if c not in 'ACGT'); tot['genomes_with_nonACGT']+= n>0; tot['nonACGT_total']+=n
    tot['N_total']+=s.count('N')
print(tot)
# versionless ids
vl=[k for k in f2 if '.' not in k]; print(len(vl))
hdr={}
for l in open(N+"results/PEDV_genomes_all.fasta",encoding='utf-8'):
    if l[0]=='>': p=l[1:].split(); hdr[p[0]]=l.strip()
print([hdr[k] for k in vl[:2]])
# f1 vs f2 for accessions in f1 whose f2 id is versionless
miss=[k for k in set(f1)-set(f2)]
print(sum(1 for k in miss if k.split('.')[0] in f2), len(miss))
ok2=sum(1 for k in miss if k.split('.')[0] in f2 and (f1[k]==f2[k.split('.')[0]] or re.sub('[^ACGT]','',f2[k.split('.')[0]])==f1[k])); print('versionless matches',ok2)
# genomes in f2 with ambiguity at positions in N / candidate regions later
