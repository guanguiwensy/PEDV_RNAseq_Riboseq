"""Shared inputs for Result 6 genomics (read-only inputs; outputs only into Result6_genomics_figures)."""
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
import numpy as np, pandas as pd, json, collections
from Bio.Data.CodonTable import unambiguous_dna_by_id
B=PEDV_K+""; SL=B+"Supplementary_Figure_L/"; OUT=B+"Result6_genomics_figures/"
REFFA=PEDV_P+"refbuild_GIb/v1.1/virus/PEDV_GIb.fa"
ref="".join(l.strip() for l in open(REFFA) if l[0]!='>').upper(); LREF=len(ref)
CDS={"ORF1a":(293,12622),"ORF1b":(12592,20613),"S":(20610,24758),"ORF3":(24758,25033),"E":(25361,25600),"M":(25608,26288),"N":(26300,27625)}
ORDER=list(CDS)
COL={"ORF1a":"#1F3A93","ORF1b":"#7FB3E6","S":"#1B9E8A","ORF3":"#B03A7A","E":"#E8833A","M":"#D9A441","N":"#3C9A3C"}
tbl=unambiguous_dna_by_id[1].forward_table
STOPS={"TAA","TAG","TGA"}; AUG={"ATG"}; NC={"CTG","GTG","TTG","ATA","ATT","ATC","ACG"}; STARTS=AUG|NC
codon_list=[a+b+c for a in "ACGT" for b in "ACGT" for c in "ACGT"]
code={s:i for i,s in enumerate(codon_list)}
aa_of=np.array([tbl.get(c,'*') for c in codon_list])
isstop=np.array([c in STOPS for c in codon_list]); isaug=np.array([c in AUG for c in codon_list]); isnc=np.array([c in NC for c in codon_list])
isstart=isaug|isnc
enc=np.full(256,4,np.uint8)
for i,c in enumerate("ACGT"): enc[ord(c)]=i
refarr=enc[np.frombuffer(ref.encode(),np.uint8)].astype(np.int64)
def families(): return pd.read_csv(SL+"PEDV_reference_ORF_families.csv",encoding='utf-8-sig')
def isoforms(): return pd.read_csv(SL+"PEDV_reference_ORF_isoforms.csv",encoding='utf-8-sig')
def fig6fam(): return pd.read_csv(PEDV_P+"global_perturbation_screen/plus1_orf_capacity/tables/plus1_ORF_families_gt40aa.tsv",sep='\t',comment='#')
def region_of(pos):
    return [n for n,(a,b) in CDS.items() if a<=pos<=b]
