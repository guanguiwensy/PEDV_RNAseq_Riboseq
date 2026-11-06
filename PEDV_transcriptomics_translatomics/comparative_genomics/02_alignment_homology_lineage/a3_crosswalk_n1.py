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
import csv,sys
sys.stdout.reconfigure(encoding='utf-8')
from Bio import Align
from Bio.Seq import Seq
G=PEDV_P+r"refbuild_GIb/v1.1/"
gib="".join(l.strip() for l in open(G+"virus/PEDV_GIb.fa") if l[0]!='>').upper()
print('GIb len',len(gib))
n1="".join(l.strip() for l in open(PEDV_K+r"N1_conservation_analysis/results/N1_reference.nt.fa") if l[0]!='>').upper()
print('N1 len',len(n1), 'aa',len(str(Seq(n1).translate()).rstrip('*')))
# independent local alignment against GIb 26000..27959 (no coordinate subtraction)
al=Align.PairwiseAligner(mode='local',match_score=2,mismatch_score=-3,open_gap_score=-6,extend_gap_score=-2)
reg0=25500
reg=gib[reg0:]
a=al.align(reg,n1)[0]
(ts,te),(qs,qe)=a.aligned[0][0],a.aligned[0][-1]
tb=[(int(x[0])+reg0+1,int(x[1])+reg0) for x in a.aligned[0]]
print('local aln score',a.score,'GIb span',tb[0][0],tb[-1][1],'query span',int(a.aligned[1][0][0])+1,int(a.aligned[1][-1][1]))
# count matches
m=mm=0
for (t0,t1),(q0,q1) in zip(a.aligned[0],a.aligned[1]):
    for i in range(int(t1-t0)):
        m+= reg[int(t0)+i]==n1[int(q0)+i]; mm+= reg[int(t0)+i]!=n1[int(q0)+i]
print('matches',m,'mismatch',mm,'gap blocks',len(a.aligned[0])-1)
s,e=tb[0][0],tb[-1][1]
print('GIb codons at start/stop:',gib[s-1:s+2], gib[e-3:e], 'len',e-s+1,'aa',(e-s+1)//3-1)
aa=str(Seq(gib[s-1:e]).translate()); print(aa)
print('N1 aa',str(Seq(n1).translate()))
# project's coordinate map
mp={}
for r in csv.DictReader(open(G+"coordinate_map/LZ202401ref_to_GIb.tsv",encoding='utf-8'),delimiter='\t'):
    if r['LZ202401_GIIa_pos']: mp[int(r['LZ202401_GIIa_pos'])]=r
for p in (26434,26772,26379,27704,26300+79):
    r=mp.get(p); print('map LZ',p,'->',r and (r['GIb_pos'],r['status']))
print('subtraction 79:',26434-79,26772-79)
# offsets along genome for LZ->GIb
import collections
off=collections.Counter()
last=None; segs=[]
for p in sorted(mp):
    r=mp[p]
    if r['GIb_pos'] in ('','-','NA'): continue
    try: o=p-int(r['GIb_pos'])
    except: continue
    if o!=last: segs.append((p,o)); last=o
print('offset change points (LZ pos, LZ-GIb):',segs)
# ORF1a-x22 check
st,sp=8445,8568
print('x22 start codon',gib[st-1:st+2],'stop',gib[sp-1:sp+2],'ncodons',(sp-st)//3,'frame',(st-1)%3,(sp-1)%3)
cods=[gib[i:i+3] for i in range(st-1,sp+2,3)]
print(len(cods),cods[0],cods[-1],[i for i,c in enumerate(cods[:-1]) if c in('TAA','TAG','TGA')])
print(str(Seq(gib[st-1:sp-1]).translate()))
