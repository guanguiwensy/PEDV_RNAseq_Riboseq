#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Supplementary Figure 13 — translation-associated features of PEDV footprints (panels A–D).
全部输入为 source_data/ 中的 CSV；数据在绘图前已落在坐标轴范围内，所有 artist 关闭 clip，热图为 pcolormesh 矢量网格，不栅格化。"""
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
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.dirname(HERE); SD = f"{OUT}/source_data"
READ_SET = sys.argv[1] if len(sys.argv) > 1 else "host_criteria"
NAME = "Supplementary_Figure_13" if READ_SET == "host_criteria" else f"Supplementary_Figure_13_{READ_SET.replace(' ', '_')}_sensitivity"
MM = 1 / 25.4
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.linewidth": 0.7, "svg.fonttype": "none", "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
                     "xtick.labelsize": 7, "ytick.labelsize": 7})
rc = lambda f: pd.read_csv(f"{SD}/{f}", comment="#")
A = rc("panel_A_strand_counts.csv"); Bd = rc("panel_B_viral_frame_periodicity.csv"); C = rc("panel_C_ORF_frame_enrichment.csv")
D = rc("panel_D_host_PEDV_length_profiles.csv"); Q = rc("panel_D_length_QC_summary.csv")
SAMPLES = [(SR.lib("RIBO_VERO_CHX_1MOI_5h"), "CHX", 5), (SR.lib("RIBO_VERO_CHX_1MOI_8h"), "CHX", 8), (SR.lib("RIBO_VERO_CHX_1MOI_16h"), "CHX", 16),
           (SR.lib("RIBO_VERO_HAR_1MOI_5h"), "HAR", 5), (SR.lib("RIBO_VERO_HAR_1MOI_8h"), "HAR", 8), (SR.lib("RIBO_VERO_HAR_1MOI_16h"), "HAR", 16)]
IDS = [s for s, _, _ in SAMPLES]; DRUG = {s: d for s, d, _ in SAMPLES}; TIME = {s: t for s, _, t in SAMPLES}
DRUGC = {"CHX": "#333333", "HAR": "#D55E00"}
PEDV, HOST = "#0072B2", "#009E73"; FC = ["#7B3294", "#56B4E9", "#999999"]; MINUS = "#E8A0A0"
SHADE_P = {5: "#9ECAE1", 8: "#4292C6", 16: "#08306B"}; SHADE_H = {5: "#A1D99B", 8: "#41AB5D", 16: "#00441B"}
ORF_ORDER = ["ORF1a", "ORF1b", "S", "ORF3", "E", "M", "N", "all_7ORF"]; ORF_LAB = {"all_7ORF": "All ORFs"}
W, H = 185.0, 252.0
fig = plt.figure(figsize=(W * MM, H * MM))
def axmm(l, t, w, h): return fig.add_axes([l / W, 1 - (t + h) / H, w / W, h / H])
def letter(x, y, s): fig.text(x / W, 1 - y / H, s, fontsize=14, fontweight="bold", ha="left", va="top")
xlab = lambda s: f"{DRUG[s]}\n{TIME[s]} h"

# ---------------- A ----------------
letter(1, 1, "A")
ax = axmm(15, 10, 56, 44)
a = A.set_index("sample").loc[IDS]
x = np.arange(6)
ax.bar(x, a.plus_fraction, color=PEDV, width=0.7, linewidth=0)
ax.bar(x, a.minus_fraction, bottom=a.plus_fraction, color=MINUS, width=0.7, linewidth=0)
for i, s in enumerate(IDS): ax.text(i, 1.02, f"{a.loc[s, 'strand_assignable_reads']/1e6:.2f}M", ha="center", va="bottom", fontsize=7)
ax.axvline(2.5, color="#BDBDBD", lw=0.7, ls="--", zorder=0)
ax.set_ylim(0, 1); ax.set_xlim(-0.6, 5.6); ax.set_xticks(x); ax.set_xticklabels([xlab(s) for s in IDS], fontsize=7)
for t_, s in zip(ax.get_xticklabels(), IDS): t_.set_color(DRUGC[DRUG[s]])
ax.set_yticks([0, 0.25, 0.5, 0.75, 1]); ax.set_yticklabels(["0", "0.25", "0.5", "0.75", "1"])
ax.set_ylabel("Fraction of strand-assignable\nPEDV reads", fontsize=7.5)
ax.legend([Rectangle((0, 0), 1, 1, fc=PEDV), Rectangle((0, 0), 1, 1, fc=MINUS)], ["(+) strand", "(−) strand"], frameon=False, fontsize=7, loc="upper center",
          bbox_to_anchor=(0.5, -0.20), ncol=2, handlelength=1.0, columnspacing=1.2)
ax.text(0.0, 1.14, "reads (M)", transform=ax.transAxes, fontsize=7, color="#555555")
ax2 = axmm(92, 10, 18, 44)
y = np.arange(6)
ax2.scatter(a.minus_fraction * 100, y, s=22, color=MINUS, edgecolor="#B05050", linewidth=0.5, zorder=3, label="PEDV (−) reads")
ax2.scatter(a.host_antisense_fraction_of_host_psites * 100, y, s=20, marker="D", facecolor="white", edgecolor="#555555", linewidth=0.8, zorder=3, label="host antisense P-sites")
for yy, m_, h_ in zip(y, a.minus_fraction * 100, a.host_antisense_fraction_of_host_psites * 100): ax2.plot([m_, h_], [yy, yy], color="#BDBDBD", lw=0.7, zorder=1)
ax2.set_ylim(5.6, -0.6); ax2.set_yticks(y); ax2.set_yticklabels([f"{DRUG[s]} {TIME[s]} h" for s in IDS], fontsize=7)
for t_, s in zip(ax2.get_yticklabels(), IDS): t_.set_color(DRUGC[DRUG[s]])
ax2.set_xlim(0, 11); ax2.set_xticks([0, 5, 10]); ax2.set_xlabel("% of reads / P-sites", fontsize=7)
ax2.legend(frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.45, -0.20), handletextpad=0.3, borderaxespad=0)

# ---------------- B ----------------
letter(116, 1, "B")
ax = axmm(130, 10, 52, 44)
b = Bd[Bd.read_set == READ_SET].set_index("sample").loc[IDS]
bot = np.zeros(6)
for k, (col, fc, nm) in enumerate(zip(["F0", "F1", "F2"], FC, ["Frame 0", "Frame 1", "Frame 2"])):
    ax.bar(x, b[col].values, bottom=bot, color=fc, width=0.7, linewidth=0, label=nm); bot += b[col].values
for i, s in enumerate(IDS): ax.text(i, 1.02, f"{b.loc[s, 'usable_psites']/1e6:.2f}M", ha="center", va="bottom", fontsize=7)
ax.axvline(2.5, color="#BDBDBD", lw=0.7, ls="--", zorder=0)
ax.set_ylim(0, 1); ax.set_xlim(-0.6, 5.6); ax.set_xticks(x); ax.set_xticklabels([xlab(s) for s in IDS], fontsize=7)
for t_, s in zip(ax.get_xticklabels(), IDS): t_.set_color(DRUGC[DRUG[s]])
ax.set_yticks([0, 0.25, 0.5, 0.75, 1]); ax.set_yticklabels(["0", "0.25", "0.5", "0.75", "1"])
ax.set_ylabel("Fraction of usable viral\nORF P-sites", fontsize=7.5)
ax.legend(frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=3, handlelength=1.0, columnspacing=1.0)
ax.text(0.0, 1.14, "usable P-sites (M)", transform=ax.transAxes, fontsize=7, color="#555555")

# ---------------- C ----------------
letter(1, 68, "C")
axC = axmm(24, 80, 128, 54)
cc = C[C.read_set == READ_SET]
cmap = plt.get_cmap("RdBu_r"); norm = TwoSlopeNorm(vmin=-0.5, vcenter=0.0, vmax=1.25)
nrow, ncol = len(ORF_ORDER), len(IDS)
for i, orf in enumerate(ORF_ORDER):
    for j, s in enumerate(IDS):
        r = cc[(cc["sample"] == s) & (cc.orf == orf)].iloc[0]
        x0, y0 = j, i
        if r.status == "valid":
            axC.add_patch(Rectangle((x0, y0), 1, 1, fc=cmap(norm(r.ribosome_fraction_f)), ec="white", lw=0.8))
            txt = f"{100*r.F0:.0f}"; col = "white" if (r.ribosome_fraction_f > 0.85 or r.ribosome_fraction_f < -0.3) else "black"
            axC.text(x0 + 0.5, y0 + 0.5, txt, ha="center", va="center", fontsize=7.5, color=col)
        elif r.status == "low_depth":
            axC.add_patch(Rectangle((x0, y0), 1, 1, fc="#E0E0E0", ec="white", lw=0.8, hatch="////", linewidth=0.8))
            axC.text(x0 + 0.5, y0 + 0.5, f"{100*r.F0:.0f}", ha="center", va="center", fontsize=7.5, color="#444444", style="italic")
        else:
            axC.add_patch(Rectangle((x0, y0), 1, 1, fc="#BDBDBD", ec="white", lw=0.8))
            axC.text(x0 + 0.5, y0 + 0.5, "NA", ha="center", va="center", fontsize=7.5, color="#444444")
axC.set_xlim(0, ncol); axC.set_ylim(nrow, 0); axC.axvline(3, color="white", lw=2.0); axC.axhline(nrow - 1, color="white", lw=2.0)
axC.set_xticks(np.arange(ncol) + 0.5); axC.set_xticklabels([f"{DRUG[s]} {TIME[s]} h" for s in IDS], fontsize=7.5)
for t_, s in zip(axC.get_xticklabels(), IDS): t_.set_color(DRUGC[DRUG[s]])
axC.set_yticks(np.arange(nrow) + 0.5); axC.set_yticklabels([ORF_LAB.get(o, o) for o in ORF_ORDER], fontsize=7.5)
for sp in axC.spines.values(): sp.set_visible(False)
axC.tick_params(length=0); axC.xaxis.tick_top()
axC.set_ylabel("PEDV ORF", fontsize=8)
for j, s in enumerate(IDS):
    n = cc[(cc["sample"] == s) & (cc.orf == "all_7ORF")].usable_psites.iloc[0]
    axC.text(j + 0.5, nrow + 0.25, f"{n/1e6:.2f}M", ha="center", va="top", fontsize=7)
axC.text(-0.1, nrow + 0.25, "usable P-sites", ha="right", va="top", fontsize=7, color="#555555")
cax = axmm(158, 86, 3.5, 40)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); cb = fig.colorbar(sm, cax=cax)
cb.set_label("f = (F0 − 1/3) / (0.80 − 1/3)", fontsize=7); cb.ax.tick_params(labelsize=7); cb.set_ticks([-0.5, 0, 0.5, 1.0])
cb.solids.set_rasterized(False); cb.solids.set_edgecolor("face")
leg = [Rectangle((0, 0), 1, 1, fc=cmap(norm(0.6)), ec="none"), Rectangle((0, 0), 1, 1, fc="#E0E0E0", hatch="////", ec="#888888", lw=0),
       Rectangle((0, 0), 1, 1, fc="#BDBDBD")]
fig.legend(leg, ["valid (≥1,000 reads); number = frame-0 %", "low depth (100–999 reads)", "missing (<100 reads)"], loc="upper left",
           bbox_to_anchor=(24 / W, 1 - 143 / H), ncol=3, frameon=False, fontsize=7, handlelength=1.2, columnspacing=1.2)

# ---------------- D ----------------
letter(1, 151, "D")
profiles = [("PEDV (+) all uniquely mapped", "Fraction of PEDV (+)\nuniquely mapped reads", SHADE_P),
            ("PEDV (+) 5' end in ORF-unique region", "Fraction of PEDV footprints\nin ORF-unique regions", SHADE_P),
            ("Host CDS-body (5' end)", "Fraction of host CDS-body\nfootprints", SHADE_H)]
titles = ["PEDV, all (+) reads", "PEDV, ORF-unique regions", "Host CDS reference"]
for k, (pn, ylab, sh) in enumerate(profiles):
    ax = axmm(20 + k * 58, 162, 44, 34)
    d = D[D.profile == pn]
    ax.axvspan(28.5, 31.5, color="#EDEDED", zorder=0, lw=0)
    for s in IDS:
        g = d[(d["sample"] == s) & d.aligned_length.between(20, 40)].sort_values("aligned_length")
        ax.plot(g.aligned_length, g.fraction, color=sh[TIME[s]], lw=1.1, ls="-" if DRUG[s] == "CHX" else (0, (3, 1.6)))
    ax.set_xlim(20, 40); ax.set_ylim(0, 0.4); ax.set_xticks([20, 25, 30, 35, 40]); ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4])
    ax.set_title(titles[k], fontsize=7.5, pad=3); ax.set_xlabel("Aligned footprint length (nt)", fontsize=7.5)
    ax.set_ylabel(ylab, fontsize=7)
hl = [Line2D([], [], color=PEDV, lw=1.2, ls="-"), Line2D([], [], color=PEDV, lw=1.2, ls=(0, (3, 1.6))),
      Line2D([], [], color=SHADE_P[5], lw=3), Line2D([], [], color=SHADE_P[8], lw=3), Line2D([], [], color=SHADE_P[16], lw=3)]
fig.legend(hl, ["CHX", "HAR", "5 h", "8 h", "16 h  (shade; PEDV blue, host green)"], loc="upper left", bbox_to_anchor=(20 / W, 1 - 208 / H), ncol=5, frameon=False,
           fontsize=7, handlelength=2.0, columnspacing=1.2)
# D4：16 hpi 专项（随时间）
fig.text(1 / W, 1 - 243 / H, "", fontsize=7)
metrics = [("median_length", "Median length (nt)", (27, 32)), ("fraction_29_31", "29–31 nt fraction", (0.3, 0.75)), ("frame0", "Frame-0 fraction", (0.3, 0.9))]
for k, (m, lab, yl) in enumerate(metrics):
    ax = axmm(22 + k * 58, 226, 42, 14)
    for drug, ls_, mk in (("CHX", "-", "o"), ("HAR", (0, (3, 1.6)), "s")):
        for pn, col in (("PEDV (+) 5' end in ORF-unique region", PEDV), ("Host CDS-body (5' end)", HOST)):
            t_, v_ = [], []
            for s in [s for s in IDS if DRUG[s] == drug]:
                q = Q[(Q["sample"] == s) & (Q.profile == pn)].iloc[0]
                if m == "frame0":
                    v = q["frame0_viral_ORF_unique_host_criteria"] if pn.startswith("PEDV") else q["frame0_host_CDS_body_host_criteria"]
                else: v = q[m]
                t_.append(TIME[s]); v_.append(v)
            ax.plot(t_, v_, color=col, ls=ls_, marker=mk, ms=3.5, lw=1.0)
    ax.set_xlim(3, 18); ax.set_xticks([5, 8, 16]); ax.set_ylim(*yl); ax.set_yticks({"median_length": [28, 30, 32], "fraction_29_31": [0.3, 0.5, 0.7], "frame0": [0.3, 0.6, 0.9]}[m]); ax.set_xlabel("hpi", fontsize=7, labelpad=1); ax.set_ylabel(lab, fontsize=7)
fig.legend([Line2D([], [], color=PEDV, lw=1.2), Line2D([], [], color=HOST, lw=1.2), Line2D([], [], color="#555555", lw=1.0, ls="-", marker="o", ms=3.5),
            Line2D([], [], color="#555555", lw=1.0, ls=(0, (3, 1.6)), marker="s", ms=3.5)], ["PEDV ORF footprints", "host CDS-body footprints", "CHX", "HAR"],
           loc="upper left", bbox_to_anchor=(20 / W, 1 - 217.0 / H), ncol=4, frameon=False, fontsize=7, handlelength=1.8, columnspacing=1.2)

# ---------------- 关闭 clip，导出 ----------------
for ax_ in fig.axes:
    for art in [ax_, *ax_.get_children()]:
        try: art.set_clip_on(False)
        except Exception: pass
for art in fig.legends + fig.texts + fig.patches: art.set_clip_on(False)
LOCKED = []
for ext, kw in (("pdf", {}), ("svg", {}), ("png", dict(dpi=600)), ("jpg", dict(dpi=600, pil_kwargs=dict(quality=95)))):
    try: fig.savefig(f"{OUT}/{NAME}.{ext}", facecolor="white", **kw)
    except PermissionError: LOCKED.append(ext)
plt.close(fig)
if LOCKED: print("LOCKED:", LOCKED)
print("saved", NAME)
