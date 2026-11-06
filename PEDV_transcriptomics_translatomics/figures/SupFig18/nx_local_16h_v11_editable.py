#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""PK1_VERO_ORF_VALIDATION_v1 第二步：N-x1、N-x2 的 Vero/PK1 局部真实数据对照（图 3）与描述性证据表。
底层数据处理与既有口径一致：29 nt P-site（5′端 + 12 nt）、+ 链、frame(p) = (p - ORF_start) mod 3（ORF_start = 父 ORF N 起点 26300）、PEDV_GIb（refdata_GIb_v1.1）坐标。
Vero 数据：SupM_alt_ORF_evidence/cache/psite_by_length.pkl（Fig5/6 审计用的同一缓存，含 RIBO_VERO_CHX_1MOI_16h）；PK1 数据：PK1_psite 交付的长表（29 nt 子集）。
不计算任何 P 值；表中全是计数、份额和比值（描述性）。"""
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
import os, sys, gzip, csv, pickle, collections, json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib as _mpl
_mpl.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "font.family": "sans-serif",
                      "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"], "axes.unicode_minus": False})

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpecFromSubplotSpec

ROOT = PEDV_P+""
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = HERE + "/data"
OUT = sys.argv[1]; os.makedirs(OUT + "/figures_editable", exist_ok=True)
VLEN = 27959; NSTART = 26300
FC = {0: "#7B3294", 1: "#56B4E9", 2: "#999999"}
COND = [("CHX", 16), ("HAR", 16)]   # 用户指示：只用 16 h
VERO_SID = {("CHX", 16): SR.lib("RIBO_VERO_CHX_1MOI_16h"), ("HAR", 16): SR.lib("RIBO_VERO_HAR_1MOI_16h")}
PK1_SID = {("CHX", 16): SR.lib("RIBO_PK1_CHX_1MOI_16h"), ("HAR", 16): SR.lib("RIBO_PK1_HAR_1MOI_16h")}
VERO_RNA = {("CHX", 16): SR.lib("RNA_VERO_CHX_1MOI_16h"), ("HAR", 16): SR.lib("RNA_VERO_HAR_1MOI_16h")}   # v1.1：按用户样本表配对（）
PK1_RNA = {("CHX", 16): SR.lib("RNA_PK1_CHX_1MOI_16h"), ("HAR", 16): SR.lib("RNA_PK1_HAR_1MOI_16h")}

# ---------- 候选（取自上一任务冻结的 candidate_list_frozen.tsv，即 Vero 既有结果）----------
cl = [l for l in open(f"{ROOT}/PK1_orf_evidence/candidate_list_frozen.tsv", encoding="utf-8") if not l.startswith("#")]
CAND = {}
for r in csv.DictReader(cl, delimiter="\t"):
    if r["alias_name"] in ("N-x1", "N-x2"):
        sv = [(int(x.split(":")[0]), x.split(":")[1]) for x in r["start_variants"].split(";") if x]
        CAND[r["alias_name"]] = dict(starts=sv, stop=int(r["stop"]), stop_last=int(r["stop_last_nt"]), best={k: r[k] for k in ("start_HARTIS_HAR_best", "start_HARTIS_CHX_best", "start_PRICE_best", "start_Fig6_plot")},
                                     length_aa=r["length_aa"], evidence=r["vero_strongest_evidence_type"], detail=r["vero_evidence_detail"])
print({k: (v["starts"], v["stop"]) for k, v in CAND.items()})

# ---------- 数据 ----------
cache = pickle.load(open(f"{ROOT}/SupM_alt_ORF_evidence/cache/psite_by_length.pkl", "rb"))
VERO = {c: np.asarray(cache["A"][VERO_SID[c]][29], float) for c in COND}
PK1 = {c: np.zeros(VLEN + 2) for c in COND}
rev = {v: k for k, v in PK1_SID.items()}
with gzip.open(f"{ROOT}/PK1_psite/psite_counts_PEDV_GIb_PK1_29nt.tsv.gz", "rt", encoding="utf-8") as f:
    for line in f:
        if line.startswith("#") or line.startswith("sample_id"): continue
        p = line.rstrip("\n").split("\t")
        if p[3] != "+" or p[5] != "29": continue
        if rev[p[0]] not in COND: continue
        PK1[rev[p[0]]][int(p[4])] += int(p[8])
def load_depth(path):
    a = np.zeros(VLEN + 2)
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        for l in f:
            c = l.rstrip("\n").split("\t")
            if len(c) == 3: a[int(c[1])] = float(c[2])   # contig pos depth
            elif len(c) == 2: a[int(c[0])] = float(c[1])  # pos depth
    return a
RNA_V = {c: load_depth(f"{DATA}/vero_rna_depth_samecond_{VERO_RNA[c]}.tsv.gz") for c in COND}   # v1.1：与 PK1 同条件重算（唯一比对、未去重）
RNA_P = {c: load_depth(f"{DATA}/rna_depth_{PK1_RNA[c]}.tsv.gz") for c in COND}
HOST_V = {SR.lib("RIBO_VERO_CHX_1MOI_8h"): (0.841, 0.0775, 0.0815), SR.lib("RIBO_VERO_CHX_1MOI_16h"): (0.8454, 0.099, 0.0555), SR.lib("RIBO_VERO_HAR_1MOI_8h"): (0.6701, 0.1714, 0.1585), SR.lib("RIBO_VERO_HAR_1MOI_16h"): (0.8637, 0.0815, 0.0548)}  # QC_summary.tsv body_F0/F1/F2
HOST_P = {}
for r in csv.DictReader(open(f"{DATA}/pk1_host_body_frame_shares_29nt.tsv"), delimiter="\t"): HOST_P[r["sample_id"]] = (float(r["F0"]), float(r["F1"]), float(r["F2"]), int(r["body_n"]))
frame = lambda p: (p - NSTART) % 3

# ---------- 描述性证据表 ----------
rows = []
def region_frames(a, lo, hi):
    n = [0.0, 0.0, 0.0]; pos1 = []
    for p in range(lo, hi + 1):
        v = a[p]
        if v: n[frame(p)] += v
        if v and frame(p) == 1: pos1.append(v)
    return n, sorted(pos1, reverse=True)
for name, c in CAND.items():
    s0 = min(s for s, _ in c["starts"]); lo = s0 + 30; hi = c["stop"] - 30
    for cell, D in (("Vero", VERO), ("PK1", PK1)):
        for cond in COND:
            a = D[cond]; T = a.sum()
            n, f1 = region_frames(a, lo, hi); tot = sum(n)
            # N 主 ORF 内、不含两个候选区的 F1 份额（病毒内部对照）
            n_ref, _ = region_frames(a, 26991, 27580)
            tr = sum(n_ref)
            sid = VERO_SID[cond] if cell == "Vero" else PK1_SID[cond]
            host = HOST_V[sid] if cell == "Vero" else HOST_P[sid][:3]
            nzf1 = len(f1); npos1 = sum(1 for p in range(lo, hi + 1) if frame(p) == 1)
            top1 = f1[0] if f1 else 0; top3 = sum(f1[:3])
            f1_share = n[1] / tot if tot else float("nan")
            r = dict(candidate=name, cell=cell, condition="%s %d h" % cond, library=sid, total_29nt_Psites_plus=int(T), body_window="%d-%d" % (lo, hi),
                     body_n0=int(n[0]), body_n1=int(n[1]), body_n2=int(n[2]), body_total=int(tot),
                     body_F1_share=round(f1_share, 4) if tot else "NA", body_F0_share=round(n[0] / tot, 4) if tot else "NA", body_F2_share=round(n[2] / tot, 4) if tot else "NA",
                     host_body_F1_share_29nt=round(host[1], 4), F1_over_host_F1=round(f1_share / host[1], 2) if tot else "NA",
                     Nbody_noncand_F1_share=round(n_ref[1] / tr, 4) if tr else "NA", F1_over_Nbody_F1=round(f1_share / (n_ref[1] / tr), 2) if tot and tr else "NA",
                     F1_positions_nonzero=nzf1, F1_positions_total=npos1, F1_nonzero_frac=round(nzf1 / npos1, 3), F1_top1_frac=round(top1 / n[1], 3) if n[1] else "NA",
                     F1_share_without_top3_positions=round((n[1] - top3) / (tot - top3), 4) if tot - top3 > 0 else "NA",
                     CPM29_body=round(tot / T * 1e6, 1) if T else "NA")
            rows.append(r)
# HAR/CHX 起始邻域（描述性；每个起点变体一行）
srows = []
for name, c in CAND.items():
    for s, codon in c["starts"]:
        for cell, D in (("Vero", VERO), ("PK1", PK1)):
            for cond in COND:
                a = D[cond]; T = a.sum()
                w5 = [0.0, 0.0, 0.0]
                for p in range(s, s + 13): w5[frame(p)] += a[p]
                loc = [a[p] for p in range(s - 30, s + 31) if p != s]
                f1pos = [a[p] for p in range(s - 30, s + 31) if frame(p) == 1 and p != s]
                srows.append(dict(candidate=name, start=s, codon=codon, start_frame_rel_N=frame(s), cell=cell, condition="%s %d h" % cond, total_29nt_Psites_plus=int(T),
                                  Psites_at_start=int(a[s]), W5_n0=int(w5[0]), W5_n1=int(w5[1]), W5_n2=int(w5[2]), W5_m=int(w5[1] + w5[2]),
                                  local_mean_F1_per_pos_pm30=round(float(np.mean(f1pos)), 2) if f1pos else "NA",
                                  start_over_local_F1_mean=round(a[s] / np.mean(f1pos), 2) if f1pos and np.mean(f1pos) > 0 else "NA",
                                  start_rank_among_F1_pm30=1 + sum(1 for x in f1pos if x > a[s]), CPM29_at_start=round(a[s] / T * 1e6, 1) if T else "NA"))
def wtsv(name, rr):
    with open(f"{OUT}/{name}", "w", encoding="utf-8", newline="") as f:
        f.write("# N-x1 / N-x2 局部描述性证据（Vero vs PK1，29 nt P-site，+ 链；frame 相对 N 起点 26300；无 P 值）；PK1_VERO_ORF_VALIDATION_v1 第二步\n")
        w = csv.DictWriter(f, fieldnames=list(rr[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(rr)
wtsv("Nx_body_frame_summary_16h.tsv", rows); wtsv("Nx_start_window_summary_16h.tsv", srows)

# ---------- 图 3 ----------
def draw_candidate(name):
    c = CAND[name]; s0 = min(s for s, _ in c["starts"]); stop = c["stop"]
    lo, hi = s0 - 40, stop + 40
    pos = np.arange(lo, hi + 1)
    fig = plt.figure(figsize=(17, 11), dpi=150)
    gs = fig.add_gridspec(5, 2, height_ratios=[0.7, 2.2, 2.2, 0.35, 2.4], hspace=0.55, wspace=0.12)
    # 基因模型
    ax = fig.add_subplot(gs[0, :]); ax.set_xlim(lo, hi); ax.set_ylim(0, 4); ax.axis("off")
    ax.add_patch(plt.Rectangle((max(lo, 26300), 2.6), min(hi, 27625) - max(lo, 26300), 0.7, fc="#dddddd", ec="k", lw=0.5)); ax.text(lo + 3, 3.65, "N canonical ORF (26300-27625), frame 0", fontsize=8)
    for k, (s, codon) in enumerate(c["starts"]):
        ax.plot([s, stop + 2], [1.9 - 0.45 * (k % 4), 1.9 - 0.45 * (k % 4)], color=FC[1], lw=2)
        ax.text(s, 2.05 - 0.45 * (k % 4), "%d %s" % (s, codon), fontsize=6, ha="center", va="bottom")
    ax.text(stop + 4, 1.0, "stop %d" % stop, fontsize=7)
    ax.set_title("%s (+1 frame vs N; Vero evidence: %s)\nPEDV_GIb refdata_GIb_v1.1; 29 nt RPF; P-site = 5' end + 12 nt, + strand; unit = P-sites per million 29-nt PEDV(+) P-sites; frame = (pos - 26300) mod 3" % (name, "PRICE nominal only (exploratory, uncorrected)" if "PRICE nominal only" in c["evidence"] else "see candidate_list_frozen.tsv"), fontsize=9, loc="left")
    def panel(spec, a, rna, label, cell, T, ymax, rmax, show_x):
        sub = GridSpecFromSubplotSpec(4, 1, subplot_spec=spec, hspace=0.05)
        ap = fig.add_subplot(sub[:3, 0]); ar = fig.add_subplot(sub[3, 0], sharex=ap)
        y = a[lo:hi + 1] / T * 1e6 if T else a[lo:hi + 1] * 0
        cols = [FC[frame(p)] for p in pos]
        ap.vlines(pos, 0, y, colors=cols, lw=1.0)
        for s, _ in c["starts"]: ap.axvline(s, color="k", lw=0.3, ls=":")
        ap.axvline(stop, color="r", lw=0.6)
        ap.set_xlim(lo, hi); ap.set_ylim(0, ymax); ap.tick_params(labelbottom=False, labelsize=7)
        ap.set_title("%s | %s | 29 nt P-sites in window: %d (genome total %d)" % (cell, label, int(a[lo:hi + 1].sum()), int(T)), fontsize=7.5, loc="left")
        ar.fill_between(pos, 0, rna[lo:hi + 1] / rna.sum() * 1e6, color="#bbbbbb", lw=0)
        ar.set_ylim(0, rmax); ar.tick_params(labelsize=6); ar.set_ylabel("RNA", fontsize=6)
        if not show_x: ar.tick_params(labelbottom=False)
        return ap
    for ri, cond in enumerate(COND):
        sidv, sidp = VERO_SID[cond], PK1_SID[cond]
        lab = "%s %d h" % cond
        Tv, Tp = VERO[cond].sum(), PK1[cond].sum()
        ymax = max((VERO[cond][lo:hi + 1] / Tv * 1e6).max(), (PK1[cond][lo:hi + 1] / Tp * 1e6).max()) * 1.08
        rmax = max((RNA_V[cond][lo:hi + 1] / RNA_V[cond].sum() * 1e6).max(), (RNA_P[cond][lo:hi + 1] / RNA_P[cond].sum() * 1e6).max()) * 1.08
        a1 = panel(gs[1 + ri, 0], VERO[cond], RNA_V[cond], lab + " [%s; RNA %s]" % (sidv, VERO_RNA[cond]), "Vero (ChlSab1.1 host)", Tv, ymax, rmax, ri == 1)
        a2 = panel(gs[1 + ri, 1], PK1[cond], RNA_P[cond], lab + " [%s; RNA %s]" % (sidp, PK1_RNA[cond]), "PK1 (Sscrofa11.1 host)", Tp, ymax, rmax, ri == 1)
        a2.set_yticklabels([])
        if ri == 0:
            from matplotlib.lines import Line2D
            a1.legend([Line2D([], [], color=FC[k], lw=2) for k in range(3)], ["frame 0 (N frame)", "frame +1 (candidate)", "frame +2"], fontsize=6, loc="upper left", ncol=3, frameon=False)
        a1.set_ylabel("P-sites / 10^6", fontsize=7)
    # HAR 起始区放大
    ax = fig.add_subplot(gs[3, :]); ax.axis("off"); ax.text(0, 0.2, "HAR start-region zoom: positions %d to %d (first start variant -15 to +45); dotted lines = start variants" % (s0 - 15, s0 + 45), fontsize=8)
    for k, cond in enumerate([("HAR", 16)]):
        lab = "%s %d h" % cond
        sub = gs[4, k]
        for cell_i, (cell, D) in enumerate((("Vero", VERO), ("PK1", PK1))):
            pass
    zlo, zhi = s0 - 15, s0 + 45; zpos = np.arange(zlo, zhi + 1)
    inner = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[4, :], wspace=0.28)
    k = 0
    for cond in [("HAR", 16)]:
        for cell, D in (("Vero", VERO), ("PK1", PK1)):
            T = D[cond].sum(); a = D[cond]
            axz = fig.add_subplot(inner[0, k]); k += 1
            y = a[zlo:zhi + 1] / T * 1e6
            axz.vlines(zpos, 0, y, colors=[FC[frame(p)] for p in zpos], lw=2)
            for s, codon in c["starts"]:
                if zlo <= s <= zhi: axz.axvline(s, color="k", lw=0.4, ls=":")
            axz.set_title("%s HAR %d h  (window n=%d)" % (cell, cond[1], int(a[zlo:zhi + 1].sum())), fontsize=7.5, loc="left"); axz.tick_params(labelsize=6)
            axz.set_xlim(zlo, zhi)
    ratios = ", ".join("%s %dh Vero/PK1 = %.1fx" % (cd[0], cd[1], VERO[cd].sum() / PK1[cd].sum()) for cd in COND)
    fig.text(0.01, 0.005, "Depth caveat (29-nt PEDV(+) P-site totals): " + ratios + ". Panels show descriptive signal only, not significance. MOI = 1 (PK1; confirmed by the user, not read from the data) / 1 (Vero). PK1 CHX 29-nt periodicity is weak (C2 fail).", fontsize=7)
    for ext in ("pdf", "svg"):
        fig.savefig(f"{OUT}/figures_editable/Fig3_{name}_local_Psite_Vero_vs_PK1_16h_v11.{ext}", dpi=200 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
for name in CAND: draw_candidate(name)
print("done")
