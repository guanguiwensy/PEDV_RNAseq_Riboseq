#!/usr/bin/env python3
"""Supplementary Fig.1C: RNA-seq library and mapping statistics (20 libraries)."""
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
import os
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

D = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
F = f"{D}/figures"
d = pd.read_csv(f"{D}/figdata/suppfig1/Supplementary_Fig1C_RNAseq_mapping_statistics.tsv", sep="\t")
d["__o"] = d.condition.map({"CHX": 0, "HAR": 1})
d = d.sort_values(["__o", "time_hpi", "sample_id"]).reset_index(drop=True)

n = len(d)
y = np.arange(n)[::-1]
d["y"] = y
labels = [f"{r.condition} {int(r.time_hpi)}h ({r.sample_id[-4:]})" for _, r in d.iterrows()]
lab_by_y = {yy: l for yy, l in zip(y, labels)}

plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "font.size": 8})
M = 1e6
CNTC = {"total_raw_reads": "#D9D9D9", "clean_reads": "#9E9E9E", "uniquely_mapped_reads": "#404040"}
CNTL = {"total_raw_reads": "Total raw", "clean_reads": "Clean", "uniquely_mapped_reads": "Uniquely mapped"}
PEDVC, HOSTC = "#0072B2", "#009E73"


def panelC(ax1, ax2, ax3):
    h = 0.26
    for i, k in enumerate(["total_raw_reads", "clean_reads", "uniquely_mapped_reads"]):
        ax1.barh(d.y + (1 - i) * h, d[k] / M, height=h, color=CNTC[k], label=CNTL[k], zorder=3)
    ax1.set_xlabel("Read pairs (million)", fontsize=9)
    ax1.set_xlim(0, d.total_raw_reads.max() / M * 1.02)
    ax1.legend(handles=[Patch(color=CNTC[k], label=CNTL[k]) for k in CNTC], loc="lower center",
               bbox_to_anchor=(0.5, 1.0), ncol=1, frameon=False, fontsize=7.5, handlelength=1, handleheight=1)
    # mapping rate
    ax2.scatter(d.overall_mapping_rate_pct, d.y, s=20, marker="o", facecolor="#56B4E9", edgecolor="none", label="Overall", zorder=3)
    ax2.scatter(d.uniq_mapping_rate_pct, d.y, s=20, marker="D", facecolor="#000000", edgecolor="none", label="Unique", zorder=3)
    ax2.set_xlabel("Mapping rate\n(% of clean)", fontsize=9)
    ax2.set_xlim(88, 99)
    ax2.set_xticks([90, 95])
    ax2.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=1, frameon=False, fontsize=7.5, handletextpad=0.2)
    # PEDV vs host uniquely mapped reads (log)
    ax3.scatter(d.PEDV_mapped_reads / M, d.y, s=20, marker="o", facecolor=PEDVC, edgecolor="none", label="PEDV", zorder=3)
    ax3.scatter(d.host_mapped_reads / M, d.y, s=20, marker="s", facecolor=HOSTC, edgecolor="none", label="Host", zorder=3)
    ax3.set_xscale("log")
    ax3.set_xlabel("Uniquely mapped\nreads (million)", fontsize=9)
    ax3.set_xlim(0.08, 40)
    ax3.set_xticks([0.1, 1, 10])
    ax3.set_xticklabels(["0.1", "1", "10"])
    ax3.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=1, frameon=False, fontsize=7.5, handletextpad=0.2)
    for ax in (ax1, ax2, ax3):
        ax.set_ylim(-0.6, n - 0.4)
        ax.set_yticks(y)
        ax.grid(axis="x", color="#EEEEEE", lw=0.5, zorder=0)
        for s_ in ["top", "right"]:
            ax.spines[s_].set_visible(False)
        ax.tick_params(labelsize=7.5)
    ax1.set_yticklabels(labels, fontsize=7)
    for ax in (ax2, ax3):
        ax.set_yticklabels([])

    sep = (d[d.condition == "HAR"].y.max() + d[d.condition == "CHX"].y.min()) / 2
    for ax in (ax1, ax2, ax3):
        ax.axhline(sep, color="#BBBBBB", lw=0.8, ls="--", zorder=1)


mm = 1 / 25.4
fig = plt.figure(figsize=(180 * mm, 95 * mm))
gs = fig.add_gridspec(1, 3, width_ratios=[2.0, 0.85, 0.95], left=0.165, right=0.985, bottom=0.14, top=0.84, wspace=0.18)
ax1, ax2, ax3 = [fig.add_subplot(g) for g in gs]
panelC(ax1, ax2, ax3)
fig.savefig(f"{F}/Supplementary_Figure1C_mapping_QC.png", dpi=300, facecolor="white")
fig.savefig(f"{F}/Supplementary_Figure1C_mapping_QC.pdf", facecolor="white")
plt.close(fig)
print("ranges: raw %.1f-%.1f M, clean %.1f-%.1f M, uniq %.1f-%.1f M" % (
    d.total_raw_reads.min()/M, d.total_raw_reads.max()/M, d.clean_reads.min()/M, d.clean_reads.max()/M,
    d.uniquely_mapped_reads.min()/M, d.uniquely_mapped_reads.max()/M))
print("uniq map rate %.2f-%.2f; overall %.2f-%.2f" % (
    d.uniq_mapping_rate_pct.min(), d.uniq_mapping_rate_pct.max(), d.overall_mapping_rate_pct.min(), d.overall_mapping_rate_pct.max()))
print("PEDV reads %.3f-%.3f M; host %.2f-%.2f M" % (
    d.PEDV_mapped_reads.min()/M, d.PEDV_mapped_reads.max()/M, d.host_mapped_reads.min()/M, d.host_mapped_reads.max()/M))
