#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""宿主 P-site 汇总：逐基因 CDS 计数（psite_counts_host_summary.tsv）与起始/终止密码子 P-site metagene（host_metagene_GIb.tsv）。
输入：R1 host_psite_pass.py 输出（host_psite/<id>/）。宿主侧统计，不涉及 PEDV reads。"""
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
import os, sys, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

SETROWS = [r for r in csv.DictReader((l for l in open(f"{OUT}/readlength_sets_PK1.tsv", encoding="utf-8") if not l.startswith("#")), delimiter="\t")]
INSET = collections.defaultdict(dict)
for r in SETROWS:
    if r["read_length"] != "≥150": INSET[r["sample_id"]][int(r["read_length"])] = r
def member(sid, L, k):
    r = INSET[sid].get(L)
    if k == "all_lengths": return True     # host_psite_pass 只对有 offset 的读长给 P-site
    return bool(r) and r[f"in_{k}"] == "True"

gene_rows, meta_rows = [], []
for sid in IDS:
    G = collections.defaultdict(lambda: collections.Counter()); info = {}
    for r in rd(f"{WORK}/host_psite/{sid}/host_gene_cds_psite.tsv"):
        L = int(r["L"]); t = r["transcript_id"]; info[t] = (r["gene_id"], r["gene_name"], r["cds_len_nt"])
        f = [int(r["f0"]), int(r["f1"]), int(r["f2"])]
        for k in ("strict_core", "fixed_29_30", "expanded", "all_lengths"):
            if member(sid, L, k):
                for i in range(3): G[t][f"{k}_f{i}"] += f[i]
    for t, c in G.items():
        g, gn, cl = info[t]
        row = {**mcols(sid), "transcript_id": t, "gene_id": g, "gene_name": gn, "cds_len_nt": cl}
        for k in ("strict_core", "fixed_29_30", "expanded", "all_lengths"):
            tot = sum(c[f"{k}_f{i}"] for i in range(3)); row[f"{k}_total"] = tot
            if k == "strict_core":
                for i in range(3): row[f"strict_core_f{i}"] = c[f"{k}_f{i}"]
        gene_rows.append(row)
    M = collections.defaultdict(lambda: collections.Counter())
    for r in rd(f"{WORK}/host_psite/{sid}/host_psite_metagene.tsv"):
        L = int(r["L"]); key = (r["anchor"], int(r["rel"])); n = int(r["n"])
        M[("L%d" % L)][key] += n
        if member(sid, L, "strict_core"): M["strict_core"][key] += n
    tot_sc = sum(c["strict_core_total"] for c in (x for x in gene_rows if x["sample_id"] == sid))
    for grp, c in M.items():
        if grp != "strict_core" and not (20 <= int(grp[1:]) <= 40): continue
        for (an, rel), n in sorted(c.items()):
            meta_rows.append({**mcols(sid), "length_group": grp, "anchor": an, "rel_nt": rel, "psite_count": n,
                              "per_million_strict_core_CDS_psites": f"{1e6*n/tot_sc:.3f}" if tot_sc else "NA"})
prov = lambda t: provenance(t, [
    "host: Ensembl canonical protein_coding 转录本（猪 Sscrofa11.1.110 GTF，带 CDS 的 14,944 条）；P-site = 5′端起第 offset+1 个比对碱基（offset 来自 psite_offset_PK1.tsv）",
    "CDS 计数：P-site 位于 GTF CDS 特征内（ATG 第1位起，至终止密码子前），frame = (P-site 转录本坐标 − CDS 起点) mod 3",
    "metagene：start rel 0 = ATG 的 A；stop rel 0 = 终止密码子第 1 位（CDS 末位 +1）；窗口 −60..+60；包含 UTR（按转录本外显子坐标）",
    "scripts: psite/scripts/host_psite_pass.py；psite/scripts/build_D_host_tables.py"] + bam_prov())
write(f"{OUT}/psite_counts_host_summary_PK1.tsv", gene_rows, prov("psite_counts_host_summary：宿主逐基因（canonical 转录本）CDS P-site 计数"))
write(f"{OUT}/host_metagene_PK1.tsv", meta_rows, prov("host_metagene_PK1：宿主起始/终止密码子 P-site metagene（逐读长 20–40 nt + strict_core 合并）"))
print(len(gene_rows), len(meta_rows))
