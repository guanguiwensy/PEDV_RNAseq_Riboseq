# -*- coding: utf-8 -*-
"""PK1_psite_v1 建表公共函数与常量（改自 riboseq_GIb_mapping_psite/scripts/common.py RGP_v1.3；改动见 psite/CHANGES_vs_RGP.md）。"""
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
import sample_registry as SR
import csv, os, gzip, datetime, collections
import numpy as np

VERSION = "PK1_psite_v1"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)                         # R1: analysis_PK1/psite
WORK = f"{OUT}/_work/server"                        # R1 上的统计输出（stats/、host_psite/、samples_PK1.tsv）
R1_ROOT = PEDV_SERVER_ROOT+"/analysis_PK1"
BAM_TEMPLATE = R1_ROOT + "/align_m25/ribo/{sid}/star.Aligned.sortedByCoord.out.bam"
VCHR = "PEDV_GIb"; VLEN = 27959
ORFS = [("ORF1a", 293, 12622), ("ORF1b", 12592, 20613), ("S", 20610, 24758), ("ORF3", 24758, 25033),
        ("E", 25361, 25600), ("M", 25608, 26288), ("N", 26300, 27625)]   # PEDV_GIb.orf_check.tsv（1-based，含终止密码子）
J_LO, J_HI = 7, 19; TRIM = 45; LONG_MIN = 35
CORE, EXT, SHORT = (29, 30), (31, 32, 33, 34), (23, 24, 25, 26, 27, 28)
# Vero（ChlSab1.1）冻结值，仅作对照列：strict_core 29/30 nt 的 offset = 12（Module 1 M1_v1.0）；其余读长无冻结值
VERO_FROZEN = {29: 12, 30: 12}
# 深度判据（与 Module 1 同）：宿主起始密码子覆盖 reads >= 500（C1）；PEDV ORF 唯一区 5′端 reads >= 1000（C2）
MIN_HOST_START = 500; MIN_PEDV_C2 = 1000

rd = lambda p: list(csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"))

INCLUDED = {SR.lib("RIBO_PK1_CHX_1MOI_16h"), SR.lib("RIBO_PK1_HAR_1MOI_16h")}   # the two PK1 Ribo-seq libraries of the final sample table
SAMPLES = [r for r in rd(f"{WORK}/samples_PK1.tsv") if r["sample_id"] in INCLUDED]
for s in SAMPLES: s["time_h"] = int(s["time_h"])
SAMPLES.sort(key=lambda s: (s["time_h"], {"CHX": 0, "HAR": 1}[s["drug"]], s["sample_id"]))
IDS = [s["sample_id"] for s in SAMPLES]
META = {s["sample_id"]: s for s in SAMPLES}
def label(sid): m = META[sid]; return f"{m['drug']} {m['time_h']} h"
def mcols(sid): m = META[sid]; return {"sample_id": sid, "time_h": m["time_h"], "drug": m["drug"]}

def zone(L):
    if L < 23: return "lt23"
    if L <= 28: return "short_23_28"
    if L <= 30: return "core_29_30"
    if L <= 34: return "ext_31_34"
    return "LPF_ge35"

def star_final(p):
    d = {}
    for l in open(p, encoding="utf-8"):
        if "|" not in l: continue
        k, v = l.split("|", 1); d[k.strip()] = v.strip()
    return d

def provenance(title, extra):
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    lines = [f"# {title}", f"# version: {VERSION}", f"# generated: {now} (R1, Python {__import__('sys').version.split()[0]})",
             f"# server_root: {OUT}", "# coordinates: 1-based, PEDV_GIb (refdata_GIb_v1.1, 27,959 nt); host Sscrofa11.1 (Ensembl 110 GTF)"]
    lines += [f"# {x}" for x in extra]
    return "\n".join(lines) + "\n"

def write(path, rows, header_lines, cols=None):
    keys = cols or list(dict.fromkeys(k for r in rows for k in r))
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "wt", encoding="utf-8", newline="") as f:
        f.write(header_lines)
        w = csv.DictWriter(f, fieldnames=keys, delimiter="\t", restval="", extrasaction="ignore", lineterminator="\n")
        w.writeheader(); w.writerows(rows)

def bam_prov():
    rows = []
    for sid in IDS:
        q = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(f"{WORK}/stats/{sid}/qc.tsv", encoding="utf-8") if "\t" in l}
        rows.append(f"input_bam {sid}: {BAM_TEMPLATE.format(sid=sid)} size={q.get('bam_size','?')} bytes")
    return rows
