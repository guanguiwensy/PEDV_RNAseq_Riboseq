"""Stop-anchored scanning core (shared by audit and analysis). Coordinates 1-based. aa length excludes the stop codon.
stop_position = first base of the stop codon (matches Figure-6-side 'shared_stop')."""
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
STARTS={"ATG","CTG","GTG","TTG","ATA","ATT","ATC","ACG"}   # AUG + 7 near-cognate (frozen set)
STOPS={"TAA","TAG","TGA"}
def scan_seq(seq, starts=STARTS, min_aa=1):
    """Return (complete, truncated). complete: list of dict(start,stop,frame,aa,codon,stop_codon,)
    truncated: starts with no in-frame stop before the 3' end."""
    L=len(seq); comp=[]; trunc=[]
    for f in (0,1,2):
        pend=[]                      # starts since last stop (5'->3')
        for p in range(f, L-2, 3):   # 0-based
            c=seq[p:p+3]
            if c in STOPS:
                for s in pend:
                    aa=(p-s)//3
                    if aa>=min_aa: comp.append(dict(start=s+1,stop=p+1,frame=f,aa=aa,start_codon=seq[s:s+3],stop_codon=c))
                pend=[]
            elif c in starts:
                pend.append(p)
        for s in pend:
            aa=(L-s)//3
            trunc.append(dict(start=s+1,stop=None,frame=f,aa_open=aa,start_codon=seq[s:s+3]))
    return comp,trunc
