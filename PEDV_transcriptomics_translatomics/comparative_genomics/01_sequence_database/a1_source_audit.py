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
import csv,re,hashlib,collections,sys
sys.stdout.reconfigure(encoding='utf-8')
B=PEDV_K+r""; N=B+"N1_conservation_analysis/"
def fa(p):
    d=collections.OrderedDict(); cur=None; dup=[]
    for l in open(p,encoding='utf-8'):
        if l[0]=='>':
            cur=l[1:].split()[0]
            if cur in d: dup.append(cur)
            d[cur]=[]
        else: d[cur].append(l.strip().upper())
    return {k:''.join(v) for k,v in d.items()},dup,l
def tsv(p): return list(csv.DictReader(open(p,encoding='utf-8'),delimiter='\t'))
f1,d1,_=fa(N+"data_external/PEDV.fasta")
f2,d2,_=fa(N+"results/PEDV_genomes_all.fasta")
meta=tsv(N+"data_external/PEDV_metadata.tsv")
gn=tsv(N+"results/genotype_normalized.tsv"); gi=tsv(N+"results/genotype_inference.tsv"); mf=tsv(N+"results/PEDV_target_manifest.tsv")
print("fasta1",len(f1),"dup",len(d1),"fasta2",len(f2),"dup",len(d2))
print("meta rows",len(meta),"gn",len(gn),"gi",len(gi),"manifest",len(mf))
S1=set(f1);S2=set(f2);Sm={r['accession']+'.'+r['version'] for r in meta};Sg={r['accession'] for r in gn};Sgi={r['accession'] for r in gi};Smf={r['accession'] for r in mf}
print("meta ids==fasta1?",Sm==S1)
print("fasta1 not in fasta2",sorted(S1-S2)); print("fasta2 not in fasta1",len(S2-S1),sorted(S2-S1))
print("gn==fasta2",Sg==S2,"gi==gn",Sg==Sgi,"mf==gn",Smf==Sg,"gn-f2",sorted(Sg-S2),"f2-gn",sorted(S2-Sg))
# same seq identical between files
diff=[k for k in S1&S2 if f1[k]!=f2[k]]; print("shared acc w/ different seq",len(diff),diff[:5])
# PEDV.txt CDS accessions
acc=set(); n=0
for l in open(B+"PEDV.txt",encoding='utf-8'):
    if l[0]=='>':
        n+=1; m=re.match(r">lcl\|([A-Z0-9_]+\.\d+)_cds",l)
        if m: acc.add(m.group(1))
print("PEDV.txt cds",n,"distinct genomes",len(acc),"in f2",len(acc&S2),"not in f2",len(acc-S2),sorted(acc-S2)[:10],"f2 not in txt",len(S2-acc))
# manifest columns
c=collections.Counter(r['sequence_source'] for r in mf); print(c)
print(collections.Counter(r['note'] for r in mf).most_common(8))
print(collections.Counter(r['qc_status'] for r in mf))
for k in ['KP162057.1','MT783684.1','PP999774.1','AF353511.1']: print(k,k in S1,k in S2, [ (r['genotype_normalized'],r['genotype_status']) for r in gn if r['accession']==k])
# ids of GIb
gb=[r['accession'] for r in gn if r['genotype_normalized']=='GIb']; print('GIb',len(gb))
import statistics
L=[len(v) for v in f2.values()]; print('len min/med/max',min(L),statistics.median(L),max(L))
ids=[r for r in gi if r['accession'] in S2]
print('mf-f1 differ', len(Smf-S1))
open(PEDV_K+r"Supplementary_Figure_L/audit/f2_minus_f1.txt","w").write("\n".join(sorted(S2-S1)))
