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
import sys,collections
sys.stdout.reconfigure(encoding='utf-8')
exec(open(__import__('os').path.join(__import__('os').path.dirname(__REPO_FILE__), 'a1_source_audit.py'),encoding='utf-8').read().split("print(\"fasta1\"")[0])
diff=[k for k in set(f1)&set(f2) if f1[k]!=f2[k]]
import itertools
cat=collections.Counter()
ex=[]
for k in diff:
    a,b=f1[k],f2[k]
    if a.replace('N','')==b.replace('N',''): cat['same_after_removing_N']+=1
    elif len(a)==len(b):
        mm=sum(1 for x,y in zip(a,b) if x!=y); cat['same_len_mismatch']+=1; ex.append((k,len(a),len(b),mm))
    else:
        cat['len_diff']+=1; ex.append((k,len(a),len(b),None))
print(cat); print(ex[:15])
# is one a substring of other / ambiguity chars
nonacgt=collections.Counter()
for k in diff[:40]:
    for s,n in ((f1[k],'f1'),(f2[k],'f2')):
        nonacgt[(n,''.join(sorted(set(s)-set('ACGT'))))]+=1
print(nonacgt)
# are the diffs only IUPAC vs N?
c2=0
for k in diff:
    a,b=f1[k],f2[k]
    if len(a)==len(b) and all(x==y or x not in 'ACGT' or y not in 'ACGT' for x,y in zip(a,b)): c2+=1
print('diff only at ambiguity positions',c2)
# header of f2 for one
for l in open(N+"results/PEDV_genomes_all.fasta",encoding='utf-8'):
    if l.startswith('>KT199103') or l.startswith('>KR265789'): print(l.strip()[:200])
