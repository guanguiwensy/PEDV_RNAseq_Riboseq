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
import csv,sys,collections
sys.stdout.reconfigure(encoding='utf-8')
sys.path[0:0] = [REPO_ROOT + '/comparative_genomics/02_alignment_homology_lineage', REPO_ROOT + '/comparative_genomics/03_orf_scan', REPO_ROOT + '/comparative_genomics/04_retention_null_stats']
from scan_core import *
G=PEDV_P+r"refbuild_GIb/v1.1/virus/PEDV_GIb.fa"
T=PEDV_P+r"global_perturbation_screen/plus1_orf_capacity/tables/"
O=PEDV_K+r"Supplementary_Figure_L/audit/"
seq="".join(l.strip() for l in open(G) if l[0]!='>').upper()
comp,trunc=scan_seq(seq)
print('all-frame complete isoforms (>=1aa)',len(comp),'truncated',len(trunc))
ORFS=[("ORF1a",293,12622),("ORF1b",12592,20613),("S",20610,24758),("ORF3",24758,25033),("E",25361,25600),("M",25608,26288),("N",26300,27625)]
# --- reproduction of Fig6-side table
rep=[]
for name,s,e in ORFS:
    c=((s-1)%3+1)%3
    for r in comp:
        if r['frame']==c and s<=r['start']<=e and r['aa']>40:
            rep.append((f"P1_{name}_stop{r['stop']}",r['start'],r['start_codon'].replace('T','U'),r['stop'],r['stop_codon'].replace('T','U'),r['aa']))
old=list(csv.DictReader((l for l in open(T+"plus1_ORF_start_variants_gt40aa.tsv",encoding='utf-8') if not l.startswith('#')),delimiter='\t'))
oldfam=list(csv.DictReader((l for l in open(T+"plus1_ORF_families_gt40aa.tsv",encoding='utf-8') if not l.startswith('#')),delimiter='\t'))
oldset={(r['family_id'],int(r['start_position']),r['start_codon'],int(r['stop_position']),r['stop_codon'],int(r['ORF_length_aa'])) for r in old}
newset=set(rep)
print('old isoforms',len(old),len(oldset),'new',len(rep),len(newset))
print('only old',sorted(oldset-newset)[:10],'only new',sorted(newset-oldset)[:10])
print('families old',len(oldfam),'new',len({r[0] for r in rep}),'same ids',{r['family_id'] for r in oldfam}=={r[0] for r in rep})
print('sum n_start_variants',sum(int(r['n_start_variants']) for r in oldfam))
fo={r['family_id']:(int(r['n_start_variants']),r['all_candidate_starts']) for r in oldfam}
fn=collections.defaultdict(list)
for r in rep: fn[r[0]].append(r[1])
bad=[k for k in fo if fo[k][0]!=len(fn[k]) or fo[k][1]!="|".join(map(str,sorted(fn[k])))]
print('family start-list mismatches',bad)
print('by parent families',collections.Counter(k.split('_')[1] for k in fn),'isoforms',collections.Counter(r[0].split('_')[1] for r in rep))
# does any start of an old family sit in two parents' ranges?  (overlap rule)
# --- scale: all frames, thresholds
for lab,fn_ in (('>=8',lambda a:a>=8),('>=20',lambda a:a>=20),('>40',lambda a:a>40),('>=40',lambda a:a>=40)):
    for aug_only in (False,True):
        sel=[r for r in comp if fn_(r['aa']) and (r['start_codon']=='ATG' or not aug_only)]
        fam=collections.defaultdict(int)
        for r in sel: fam[(r['frame'],r['stop'])]+=1
        pf=collections.Counter(f for f,_ in fam); pi=collections.Counter(r['frame'] for r in sel)
        print(lab,'AUG-only' if aug_only else '8 starts','families',len(fam),dict(pf),'isoforms',len(sel),dict(pi))
print('truncated by frame',collections.Counter(r['frame'] for r in trunc))
print('truncated starts >=8 aa open',sum(1 for r in trunc if r['aa_open']>=8))
# stops per frame
for f in (0,1,2): print('frame',f,'stops',sum(1 for p in range(f,len(seq)-2,3) if seq[p:p+3] in STOPS))
# density per region family counts for >40 all frames by stop position region (rough)
reg=[("5UTR",1,292),("ORF1a",293,12591),("ORF1a/1b ovl",12592,12622),("ORF1b",12623,20609),("S",20610,24757),("ORF3",24758,25033),("gap",25034,25360),("E",25361,25600),("gap2",25601,25607),("M",25608,26288),("gap3",26289,26299),("N",26300,27625),("3UTR",27626,27959)]
for lab in ('>40',):
    fam={}
    for r in comp:
        if r['aa']>40: fam.setdefault((r['frame'],r['stop']),[]).append(r['start'])
    print('region  len  fam_stop_in_region  starts_in_region  per_kb(fam)  per_kb(starts)')
    for n,a,b in reg:
        nf=sum(1 for (f,s) in fam if a<=s<=b); ns=sum(1 for v in fam.values() for x in v if a<=x<=b)
        print(n,b-a+1,nf,ns,round(nf/((b-a+1)/1000),1),round(ns/((b-a+1)/1000),1))
