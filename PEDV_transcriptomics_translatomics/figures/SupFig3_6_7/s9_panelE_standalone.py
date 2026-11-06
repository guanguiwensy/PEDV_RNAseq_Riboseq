#!/usr/bin/env python3
"""Panel E 单独成图: 5′(donor) 与 3′(acceptor) 各一张, 画板 830 × 720 pt (Illustrator 单位 = 1/72 in)。
竖向 5 行 = 5 个时间点(3/5/8/12/16 hpi), 2 列 = MOI 0.1 | 1; x = 病毒基因组坐标, y = 每百万病毒 fragment 的 junction reads(对数)。
数据同总图 panel E(data/panelE_*): 仅 noncanonical, 100-nt bin, 同一 MOI×hpi 的 CHX+HAR reads 相加后除以两库病毒 fragment 之和(描述性)。"""
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
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import NullLocator

R = PEDV_P+"Fig2_Supp_junction_dynamics/"
D, OUT = R + "data/", R + "figures/PanelE_standalone/"
os.makedirs(OUT, exist_ok=True)
WPT, HPT = 830, 720
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "svg.fonttype": "none", "font.size": 9,
                     "axes.linewidth": 0.9, "xtick.major.width": 0.9, "ytick.major.width": 0.9,
                     "xtick.major.size": 3.5, "ytick.major.size": 3.5})
TIMES = [3, 5, 8, 12, 16]
TCOL = {3: "#3A2448", 5: "#3A4E89", 8: "#288D85", 12: "#62BB5D", 16: "#F0DE3E"}
TTXT = {3: "#3A2448", 5: "#3A4E89", 8: "#288D85", 12: "#4C9A47", 16: "#B7A300"}   # 文字用略深的黄/绿以保证可读
GENES = {"ORF1a": (293, 12622), "ORF1b": (12592, 20613), "S": (20610, 24758), "ORF3": (24758, 25033), "E": (25361, 25600),
         "M": (25608, 26288), "N": (26300, 27625)}
GENE_FILL = {'ORF1a': '#F3E3EC', 'ORF1b': '#E1EBF5', 'S': '#E6F1E3', 'ORF3': '#FBF0D9', 'E': '#EAE3F3', 'M': '#F8E3E3', 'N': '#E0F0F0'}
GL, RB = 27959, 100
NB = int(np.ceil(GL / RB))
MOIS = [0.1, 1.0]

pe = pd.read_csv(D + "panelE_breakpoint_histogram_pooled_DESCRIPTIVE.csv.gz")
pe = pe[pe.set == "noncanonical"].copy()
pe["b"] = ((pe.bin_start - 1) // 20 * 20 // RB).astype(int)
agg = pe.groupby(["MOI", "hpi", "breakpoint", "b"]).agg(reads=("reads", "sum"), reads_per_M_viral_frag=("reads_per_M_viral_frag", "sum")).reset_index()
agg["bin_start_nt"] = agg.b * RB + 1
agg["bin_end_nt"] = (agg.b + 1) * RB


def draw(bp, title, fname):
    fig = plt.figure(figsize=(WPT / 72, HPT / 72))
    L, Rm, top, bot = 66, 14, 66, 52
    gap_c, gap_r = 22, 7
    colw = (WPT - L - Rm - gap_c) / 2
    rowh = (HPT - top - bot - 4 * gap_r) / 5

    def axp(x, y_top, w, h):
        return fig.add_axes([x / WPT, 1 - (y_top + h) / HPT, w / WPT, h / HPT])

    fig.text(14 / WPT, 1 - 8 / HPT, "E", fontsize=14, fontweight="bold", va="top")
    fig.text(40 / WPT, 1 - 10 / HPT, title, fontsize=11, fontweight="bold", va="top")
    if bp == "5prime_donor":
        fig.text(WPT - 14 / WPT * WPT, 1 - 10 / HPT, "nt 1–100 bin: leader donors (nt 55–72) of TRS-L-dependent junctions", fontsize=8,
                 ha="right", va="top", color="#444444")
    peaks = []
    for ci, moi in enumerate(MOIS):
        x0 = L + ci * (colw + gap_c)
        fig.text((x0 + colw / 2) / WPT, 1 - 29 / HPT, f"MOI = {moi:g}", fontsize=10, fontweight="bold", ha="center", va="top")
        g_ = axp(x0, 46, colw, 11)
        g_.set_xlim(0, GL); g_.set_ylim(0, 1); g_.axis("off")
        for name, (s, e) in GENES.items():
            g_.add_patch(Rectangle((s, 0), e - s, 1, fc=GENE_FILL[name], ec="#9A9A9A", lw=0.4))
            if name in ("ORF1a", "ORF1b", "S", "M", "N"):
                g_.text((s + e) / 2, 0.5, name, ha="center", va="center", fontsize=8)
        for ri, t in enumerate(TIMES):
            ax = axp(x0, top + ri * (rowh + gap_r), colw, rowh)
            d = agg[(agg.MOI == moi) & (agg.hpi == t) & (agg.breakpoint == bp)]
            y = np.full(NB, np.nan)
            y[d.b.values] = d.reads_per_M_viral_frag.values
            ax.plot((np.arange(NB) + 0.5) * RB, y, color=TCOL[t], lw=1.0, marker="o", ms=1.8, mew=0)
            ax.set_yscale("log"); ax.set_ylim(0.03, 3000); ax.set_xlim(0, GL)
            ax.set_yticks([0.1, 10, 1000]); ax.set_yticklabels(["0.1", "10", "1000"]); ax.yaxis.set_minor_locator(NullLocator())
            ax.tick_params(labelsize=9, pad=2)
            for s_ in ("top", "right"):
                ax.spines[s_].set_visible(False)
            if ri < 4:
                ax.set_xticklabels([])
            else:
                ax.set_xticks([0, 5000, 10000, 15000, 20000, 25000])
                ax.set_xticklabels(["0", "5,000", "10,000", "15,000", "20,000", "25,000"], fontsize=9)
                ax.set_xlabel("PEDV genomic coordinate (nt)", fontsize=10, labelpad=4)
            if ci > 0:
                ax.set_yticklabels([])
            ax.text(0.008, 0.93, f"{t} hpi", transform=ax.transAxes, fontsize=10, fontweight="bold", color=TTXT[t],
                    ha="left", va="top")
            # 标注该子图最高的 2 个 100-nt bin(相距 >= 1 kb)
            dd = d.sort_values("reads_per_M_viral_frag", ascending=False)
            picked = []
            for r in dd.itertuples():
                pos = (r.b + 0.5) * RB
                if all(abs(pos - p) > 1000 for p in picked):
                    picked.append(pos)
                    peaks.append(dict(MOI=moi, hpi=t, breakpoint=bp, bin_start_nt=int(r.bin_start_nt), bin_end_nt=int(r.bin_end_nt),
                                      reads=int(r.reads), reads_per_M_viral_frag=float(r.reads_per_M_viral_frag)))
                if len(picked) == 2:
                    break
            for p, r in zip(picked, dd.drop_duplicates("b").itertuples()):
                pass
            for k_, p in enumerate(sorted(picked)):
                bi = int(p // RB)
                val = float(y[bi])
                if p < 3000:
                    continue                      # 第一个 bin(1-100 nt)=TRS-L 依赖事件的 leader donor, 不标坐标
                elif k_ == 0:
                    ha = "right"
                else:
                    ha = "left" if p < 27200 else "right"
                ax.annotate(f"{int(p - RB / 2 + 1)}", (p, val), xytext=(0, 3), textcoords="offset points", fontsize=8,
                            ha=ha, va="bottom", color="#333333", clip_on=False)
    fig.text(14 / WPT, 1 - (top + 2.5 * rowh + 2 * gap_r) / HPT, "Reads per million viral fragments (CHX + HAR pooled)", rotation=90,
             fontsize=10, ha="center", va="center")
    base = OUT + fname
    fig.savefig(base + ".pdf", facecolor="white")
    fig.savefig(base + ".svg", facecolor="white")
    fig.savefig(base + "_300dpi.png", dpi=300, facecolor="white")
    fig.savefig(base + "_830x720px.png", dpi=72, facecolor="white")
    plt.close(fig)
    cols = ["MOI", "hpi", "breakpoint", "bin_start_nt", "bin_end_nt", "reads", "reads_per_M_viral_frag"]
    agg[agg.breakpoint == bp][cols].sort_values(["MOI", "hpi", "bin_start_nt"]).to_csv(base + "_source_data.csv", index=False)
    pd.DataFrame(peaks).to_csv(base + "_labelled_peaks.csv", index=False)


draw("5prime_donor", "5′ breakpoint (donor) of noncanonical junction reads", "SuppFig_PanelE_5prime_donor_by_hpi_830x720")
draw("3prime_acceptor", "3′ breakpoint (acceptor) of noncanonical junction reads", "SuppFig_PanelE_3prime_acceptor_by_hpi_830x720")
print("ok")
