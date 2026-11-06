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
import os, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); G = os.path.dirname(HERE)
ESV = os.path.dirname(G); TASKS = os.path.dirname(os.path.dirname(ESV))
ANN = os.path.join(ESV, "strain_assembly_phylo_v1", "tables", "annot"); ASM = os.path.join(ESV, "strain_assembly_phylo_v1", "assembly")
GIB_FA = os.path.join(ESV, "phase5", "_in5", "_work", "ref", "PEDV_GIb.fa")
def read_fa(p): return "".join(l.strip() for l in open(p) if not l.startswith(">")).upper()
GIB = read_fa(GIB_FA)
DATASETS = {"LJX": dict(fa=os.path.join(ASM, "LJX_merged_denovo.fa"), lift=os.path.join(ANN, "LJXmerged_lift_to_GIb_v1.1.tsv"), ref="LJXm", label="LJX_merged"),
            "LLC-PK": dict(fa=os.path.join(ASM, "LLCPK_denovo.fa"), lift=os.path.join(ANN, "LLCPK_lift_to_GIb_v1.1.tsv"), ref="LLCPK", label="LLC-PK")}
SEQ = {k: read_fa(v["fa"]) for k, v in DATASETS.items()}
def load_maps():
    q2g = {}; g2q = {}
    for k, v in DATASETS.items():
        lf = pd.read_csv(v["lift"], sep="\t", dtype=str); L = len(lf)
        a = np.zeros(L + 2, dtype=np.int64); b = np.zeros(len(GIB) + 2, dtype=np.int64)
        for q, g in zip(lf.new_pos.astype(int), lf.GIb_pos):
            if g not in ("-", ""):
                g = int(g); a[q] = g; b[g] = q
        q2g[k] = a; g2q[k] = b
    return q2g, g2q
Q2G, G2Q = load_maps()
CODON = {}
_b = "TCAG"; _aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(_b):
    for j, b in enumerate(_b):
        for k, c in enumerate(_b): CODON[a + b + c] = _aa[16 * i + 4 * j + k]
def translate(s):
    return "".join(CODON.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))
def orf_from(seq, start1):
    """seq 1-based start -> (aa string without stop, stop_first_base_1based or None)"""
    i = start1 - 1; aa = []
    while i + 3 <= len(seq):
        c = CODON.get(seq[i:i + 3], "X")
        if c == "*": return "".join(aa), i + 1
        aa.append(c); i += 3
    return "".join(aa), None
def write_tsv(df, name):
    p = os.path.join(G, "tables", name); df.to_csv(p, sep="\t", index=False); print("wrote", name, df.shape)
