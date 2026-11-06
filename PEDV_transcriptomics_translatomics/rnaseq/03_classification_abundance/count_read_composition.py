#!/usr/bin/env python3
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
import sys, pysam
V="PEDV_GIb"
def cls(name): return "V" if name==V else "H"
def main(sample, base):
    bam=f"{base}/align/{sample}/Aligned.sorted.bam"
    log=f"{base}/align/{sample}/Log.final.out"
    total=None
    for l in open(log):
        if "Number of input reads" in l: total=int(l.split("|")[1])
    classes={}
    n_pairs_seen=0
    with pysam.AlignmentFile(bam,"rb") as b:
        for r in b.fetch(until_eof=True):
            if r.is_supplementary: continue
            s=classes.setdefault(r.query_name,set())
            s.add(cls(r.reference_name))
    c={"V":0,"H":0,"both":0}
    for s in classes.values():
        if len(s)==2: c["both"]+=1
        elif "V" in s: c["V"]+=1
        else: c["H"]+=1
    mapped=len(classes)
    print("\t".join(map(str,[sample,total,mapped,c["V"],c["H"],c["both"],total-mapped])),flush=True)
if __name__=="__main__": main(sys.argv[1],sys.argv[2])
