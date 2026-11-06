#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
import sample_registry as SR
import os, sys, csv
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib as _mpl
_mpl.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "font.family": "sans-serif",
                      "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"], "axes.unicode_minus": False})

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

OUT = sys.argv[1]; os.makedirs(OUT + "/figures_editable", exist_ok=True)
z = np.load(f"{OUT}/_plotdata.npz")
COND = [("CHX", 16), ("HAR", 16)]
CL = lambda c: "%s %d h" % c
VERO = {c: z["V29_%s_%d" % c] for c in COND}; PK1 = {c: z["P29_%s_%d" % c] for c in COND}
rows = list(csv.DictReader([l for l in open(f"{OUT}/PK1_VERO_CORRELATION.tsv", encoding="utf-8") if not l.startswith("#")], delimiter="\t"))
def val(level, target, cond, metric):
    for r in rows:
        if r["level"] == level and r["target"] == target and r["condition"] == cond and r["metric"] == metric:
            return (float(r["value"]) if r["value"] != "" else np.nan, r["ci95_lo"], r["ci95_hi"])
    return (np.nan, "", "")
f2 = lambda s: ("%.2f" % float(s)) if s not in ("", None) else "NA"
REG = [("5UTR", 1, 292, "#e0e0e0"), ("ORF1a", 293, 12622, "#cfe3f1"), ("ORF1b", 12623, 20609, "#f1e3cf"), ("S", 20610, 24757, "#d9efd0"), ("ORF3_E_M", 24758, 26299, "#f1d0d9"), ("N", 26300, 27625, "#e3d0f1"), ("3UTR", 27626, 27959, "#e0e0e0")]
SC = {"5UTR": "#888888", "ORF1a": "#4c8fc2", "ORF1b": "#c89b4c", "S": "#5aa84c", "ORF3_E_M": "#c8506f", "N": "#8a4cc8", "3UTR": "#555555"}
CV, CPK = "#0072B2", "#D55E00"
FOOT = (f"29-nt RPF; P-site = 5' end + 12 nt, + strand; PEDV_GIb (refdata_GIb_v1.1); unit = P-sites per million 29-nt PEDV(+) P-sites (CPM29); Vero host ChlSab1.1 ({SR.lib('RIBO_VERO_CHX_1MOI_16h')} CHX 16 h, {SR.lib('RIBO_VERO_HAR_1MOI_16h')} HAR 16 h), "
        f"PK1 host Sscrofa11.1 ({SR.lib('RIBO_PK1_CHX_1MOI_16h')}, {SR.lib('RIBO_PK1_HAR_1MOI_16h')}). MOI = 1 (PK1; confirmed by the user, not read from the data) / 1 (Vero). PK1 CHX 29-nt periodicity is weak (C2 fail). The two pairs share no PK1-internal comparison.")


fig, axs = plt.subplots(2, 1, figsize=(16, 7.5), sharex=True)
pos = np.arange(len(VERO[COND[0]]))
for ax, c in zip(axs, COND):
    for n, a, b, col in REG: ax.axvspan(a, b, color=col, alpha=0.55, lw=0)
    v = np.log10(1 + VERO[c]); p = np.log10(1 + PK1[c])
    ax.vlines(pos[1:], 0, v[1:], colors=CV, lw=0.5); ax.vlines(pos[1:], 0, -p[1:], colors=CPK, lw=0.5)
    ax.axhline(0, color="k", lw=0.5); ax.set_ylim(-4.6, 4.6)
    rs, lo, hi = val("region_L29", "whole_virus", CL(c), "Spearman_perpos_allpos")
    rb, lob, hib = val("region_L29", "whole_virus", CL(c), "Spearman_bin30nt")
    rp = val("region_L29", "whole_virus", CL(c), "Pearson_bin30nt")[0]
    ax.set_ylabel("%s\nlog10(1 + CPM29)" % CL(c), fontsize=9)
    ax.text(0.005, 0.93, "Vero (up) | PK1 (down, mirrored)   whole virus: Spearman per-position = %.2f [%s, %s]; Spearman 30-nt bins = %.2f [%s, %s]; Pearson 30-nt bins = %.2f" % (rs, f2(lo), f2(hi), rb, f2(lob), f2(hib), rp),
            transform=ax.transAxes, fontsize=8.5, va="top")
    ax.set_yticks([-4, -2, 0, 2, 4]); ax.set_yticklabels(["4", "2", "0", "2", "4"], fontsize=8)
for n, a, b, col in REG:
    axs[0].text((a + b) / 2, 4.45, n.replace("_", "/") if (b - a) > 600 else "", ha="center", va="top", fontsize=8, clip_on=True)
axs[-1].set_xlim(0, 27959); axs[-1].set_xlabel("PEDV_GIb position (nt, 1-based)")
axs[0].legend([Patch(color=CV), Patch(color=CPK)], ["Vero (ChlSab1.1)", "PK1 (Sscrofa11.1)"], loc="lower right", fontsize=8, ncol=2, frameon=False)
fig.suptitle("Fig 1 | PK1 vs Vero PEDV P-site profile, matched treatment and time (16 h; 29 nt)", fontsize=12, x=0.01, ha="left")
fig.text(0.01, 0.005, FOOT, fontsize=6.5, wrap=True)
fig.subplots_adjust(top=0.92, bottom=0.12, hspace=0.12, left=0.07, right=0.99)
for ext in ("pdf", "svg"): fig.savefig(f"{OUT}/figures_editable/Fig1_PK1_Vero_whole_virus_Psite_16h.{ext}", dpi=220 if ext == "png" else None)
plt.close(fig)


fig = plt.figure(figsize=(14, 12))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05], hspace=0.4, wspace=0.25)
for k, c in enumerate(COND):
    ax = fig.add_subplot(gs[0, k])
    for n, a, b, col in REG:
        v = VERO[c][a:b + 1]; p = PK1[c][a:b + 1]; m = len(v) // 30 * 30
        vb = v[:m].reshape(-1, 30).sum(1); pb = p[:m].reshape(-1, 30).sum(1)
        ax.scatter(np.log10(1 + vb), np.log10(1 + pb), s=7, color=SC[n], alpha=0.8, lw=0, label=n.replace("_", "/"))
    rs = val("region_L29", "whole_virus", CL(c), "Spearman_bin30nt"); rpv = val("region_L29", "whole_virus", CL(c), "Pearson_bin30nt")
    ax.set_title("%s (30-nt bins, n=%d)\nSpearman %.2f [%s, %s]; Pearson %.2f [%s, %s]" % (CL(c), 27959 // 30, rs[0], f2(rs[1]), f2(rs[2]), rpv[0], f2(rpv[1]), f2(rpv[2])), fontsize=9)
    ax.set_xlabel("Vero log10(1 + CPM29 per 30-nt bin)", fontsize=8); ax.set_ylabel("PK1 log10(1 + CPM29 per 30-nt bin)", fontsize=8)
    lim = max(ax.get_xlim()[1], ax.get_ylim()[1]); ax.plot([0, lim], [0, lim], color="k", lw=0.4, ls=":"); ax.tick_params(labelsize=7)
    if k == 0: ax.legend(fontsize=7, markerscale=2, frameon=False, loc="upper left", ncol=2)
REGN = ["whole_virus", "5UTR", "ORF1a", "ORF1b", "S", "ORF3_E_M", "N", "3UTR"]
def heat(ax, metric, title):
    M = np.array([[val("region_L29", r, CL(c), metric)[0] for c in COND] for r in REGN])
    im = ax.imshow(M, vmin=0, vmax=1, cmap="viridis", aspect="auto")
    for i, r in enumerate(REGN):
        for j, c in enumerate(COND):
            v, lo, hi = val("region_L29", r, CL(c), metric)
            ax.text(j, i, ("%.2f" % v) + (("\n[%s, %s]" % (f2(lo), f2(hi))) if lo != "" else ""), ha="center", va="center", fontsize=8, color="w" if v < 0.65 else "k")
    ax.set_xticks(range(len(COND))); ax.set_xticklabels([CL(c) for c in COND], fontsize=9); ax.set_yticks(range(len(REGN))); ax.set_yticklabels([r.replace("_", "/") for r in REGN], fontsize=8)
    ax.set_title(title, fontsize=9); return im
ax1 = fig.add_subplot(gs[1, 0]); im = heat(ax1, "Spearman_perpos_allpos", "Spearman rho, per-position (all positions incl. zeros)\n95% block-bootstrap CI")
ax2 = fig.add_subplot(gs[1, 1]); heat(ax2, "Spearman_bin30nt", "Spearman rho, 30-nt bins (sums of CPM29)\n95% block-bootstrap CI")
cax = fig.add_axes([0.93, 0.08, 0.012, 0.3]); cb = fig.colorbar(im, cax=cax); cb.ax.tick_params(labelsize=7); cb.set_label("Spearman rho", fontsize=7)
fig.suptitle("Fig 2 | PK1 vs Vero correlation (16 h, 29 nt): whole-virus scatter and region-wise heatmaps", fontsize=11.5, x=0.01, ha="left")
fig.text(0.01, 0.005, FOOT + " Regions are non-overlapping (first-come): 5'UTR 1-292, ORF1a 293-12622, ORF1b 12623-20609, S 20610-24757, ORF3/E/M 24758-26299, N 26300-27625, 3'UTR 27626-27959. Block bootstrap: 300-nt blocks, 500 resamples.", fontsize=6.5, wrap=True)
fig.subplots_adjust(top=0.9, bottom=0.08, left=0.07, right=0.92)
for ext in ("pdf", "svg"): fig.savefig(f"{OUT}/figures_editable/Fig2_PK1_Vero_correlation_scatter_heatmap_16h.{ext}", dpi=220 if ext == "png" else None)
plt.close(fig); print("figs done")
