# -*- coding: utf-8 -*-
"""riboseq_GIb_mapping_psite 本地建表的公共函数与常量（版本见 VERSION）。"""
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
import csv, os, gzip, datetime, collections
import numpy as np

VERSION = "RGP_v1.3"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)                         # 本地交付根目录
WORK = f"{OUT}/_work/server"                        # 从 R1 同步回来的原始统计
M1 = os.path.join(os.path.dirname(OUT), "module1_rpf_qc")
RB = os.path.join(os.path.dirname(OUT), "refbuild_GIb")
R1_ROOT = PEDV_SERVER_ROOT+"/analysis_GIb/Riboseq_GIb_psite"
VCHR = "PEDV_GIb"; VLEN = 27959
ORFS = [("ORF1a", 293, 12622), ("ORF1b", 12592, 20613), ("S", 20610, 24758), ("ORF3", 24758, 25033),
        ("E", 25361, 25600), ("M", 25608, 26288), ("N", 26300, 27625)]   # PEDV_GIb.orf_check.tsv（1-based，含终止密码子）
J_LO, J_HI = 7, 19; TRIM = 45; LONG_MIN = 35
CORE, EXT, SHORT = (29, 30), (31, 32, 33, 34), (23, 24, 25, 26, 27, 28)

rd = lambda p: list(csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"))

# Analysis samples = the six Vero Ribo-seq libraries of the FINAL sample table (metadata/sample.txt, via config/sample_registry.py).
# Libraries that were aligned earlier but are not in the table (e.g. replicates/other time points) never enter the results.
import sample_registry as SR
INCLUDED = set(SR.LIB_RIBO_VERO)
_ALL = rd(f"{WORK}/samples_GIb.tsv")
SR.assert_whitelisted([r["sample_id"] for r in _ALL if r["sample_id"] in INCLUDED], "samples_GIb.tsv")
SAMPLES = [r for r in _ALL if r["sample_id"] in INCLUDED]
for s in SAMPLES:
    s["time_h"] = int(s["time_h"])
    _a = SR.SAMPLES[SR.alias_of(s["sample_id"])]          # drug/time must agree with the registry
    assert (s["drug"], s["time_h"]) == (_a["treatment"], int(_a["hpi"])), (s["sample_id"], s["drug"], s["time_h"])
SAMPLES.sort(key=lambda s: (s["time_h"], {"CHX": 0, "HAR": 1}[s["drug"]], s["sample_id"]))
assert {s["sample_id"] for s in SAMPLES} == INCLUDED, "samples_GIb.tsv lacks some registry libraries"
IDS = [s["sample_id"] for s in SAMPLES]
META = {s["sample_id"]: s for s in SAMPLES}
def label(sid): m = META[sid]; return f"{m['drug']} {m['time_h']}h ({sid})"
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
    lines = [f"# {title}", f"# version: {VERSION}", f"# generated: {now} (local, Python {__import__('sys').version.split()[0]})",
             f"# server_root: {R1_ROOT}", "# coordinates: 1-based, PEDV_GIb (refdata_GIb_v1.1, 27,959 nt)"]
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
        rows.append(f"input_bam {sid}: {R1_ROOT}/align/{sid}/{sid}.Aligned.sortedByCoord.out.bam size={q.get('bam_size','?')} bytes")
    return rows
