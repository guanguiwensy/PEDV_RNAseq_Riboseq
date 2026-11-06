#!/usr/bin/env python3
"""Panel C 单独成图: 画板 670 × 350 (Illustrator 单位 = pt, 72 per inch → 236.4 × 123.5 mm)。
数据与总图 panel C 完全相同(data/panelC_*.csv), 配色/线型/标记沿用总图。"""
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
import shutil
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, NullLocator
import os

R = PEDV_P+"Fig2_Supp_junction_dynamics/"
D, OUT = R + "data/", R + "figures/PanelC_standalone/"
os.makedirs(OUT, exist_ok=True)
WPT, HPT = 670, 350
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "svg.fonttype": "none", "font.size": 9,
                     "axes.linewidth": 0.9, "xtick.major.width": 0.9, "ytick.major.width": 0.9,
                     "xtick.major.size": 3.5, "ytick.major.size": 3.5})
TIMES = [3, 5, 8, 12, 16]
NCC = ["TRSL_dep_noncanonical", "TRSL_indep_distant", "TRSL_indep_local"]
NCOL = {"TRSL_dep_noncanonical": "#B8860B", "TRSL_indep_distant": "#4DAF4A", "TRSL_indep_local": "#984EA3"}
NLAB = {"TRSL_dep_noncanonical": "TRS-L-dependent noncanonical", "TRSL_indep_distant": "TRS-L-independent distant",
        "TRSL_indep_local": "TRS-L-independent local"}
Cd = pd.read_csv(D + "panelC_noncanonical_per_sample.csv")
pC = pd.read_csv(D + "panelC_noncanonical_pooled_DESCRIPTIVE.csv")

fig = plt.figure(figsize=(WPT / 72, HPT / 72))


def axp(x, y, w, h):          # 以 pt 为单位, 原点左下
    return fig.add_axes([x / WPT, y / HPT, w / WPT, h / HPT])


panels = [("frac_all_", "Relative burden", "Fraction of all classified\njunction reads", False),
          ("per_M_viral_frag_", "Normalized abundance", "Junction-supporting reads per\nmillion viral fragments", True)]
for i, (col, title, ylab, logy) in enumerate(panels):
    ax = axp(72 + i * 300, 88, 255, 224)
    for n in NCC:
        c = NCOL[n]
        for moi in (0.1, 1.0):
            d = Cd[Cd.MOI == moi]
            filled = moi == 1.0
            for tr, m in (("CHX", "o"), ("HAR", "s")):
                dd = d[d.treatment == tr].sort_values("hpi")
                ax.scatter(dd.hpi + (-0.2 if tr == "CHX" else 0.2), dd[col + n], s=34, marker=m,
                           facecolor=c if filled else "white", edgecolor=c, linewidths=1.2, zorder=3)
            pp = pC[pC.MOI == moi].sort_values("hpi")
            ax.plot(pp.hpi, pp[col + n], color=c, lw=1.6, ls="-" if filled else "--", zorder=2)
    ax.set_xticks(TIMES); ax.set_xlim(2, 17)
    ax.set_xlabel("Hours post infection (hpi)", fontsize=10, labelpad=4)
    ax.set_ylabel(ylab, fontsize=10, labelpad=5)
    ax.set_title(title, fontsize=10.5, fontweight="bold", pad=8)
    ax.tick_params(labelsize=9, pad=3)
    if logy:
        ax.set_yscale("log"); ax.set_ylim(150, 1800)
        ax.set_yticks([200, 500, 1000]); ax.set_yticklabels(["200", "500", "1000"]); ax.yaxis.set_minor_locator(NullLocator())
    else:
        ax.set_ylim(0, 0.045); ax.set_yticks([0, 0.01, 0.02, 0.03, 0.04])
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.2f}"))
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

h1 = [Line2D([0], [0], color=NCOL[n], lw=2.0, label=NLAB[n]) for n in NCC]
h2 = [Line2D([0], [0], color="#444", lw=1.6, ls="-", marker="o", markerfacecolor="#444", markersize=6, label="MOI = 1"),
      Line2D([0], [0], color="#444", lw=1.6, ls="--", marker="o", markerfacecolor="white", markersize=6, label="MOI = 0.1"),
      Line2D([0], [0], marker="o", color="none", markerfacecolor="#888", markeredgecolor="#888", markersize=6, label="CHX"),
      Line2D([0], [0], marker="s", color="none", markerfacecolor="#888", markeredgecolor="#888", markersize=6, label="HAR")]
l1 = fig.legend(handles=h1, loc="lower center", bbox_to_anchor=(0.5, 26 / HPT), ncol=3, frameon=False, fontsize=9,
                handlelength=2.2, columnspacing=2.2, borderaxespad=0)
fig.legend(handles=h2, loc="lower center", bbox_to_anchor=(0.5, 5 / HPT), ncol=4, frameon=False, fontsize=9,
           handlelength=2.2, columnspacing=2.2, borderaxespad=0)
fig.text(14 / WPT, 1 - 8 / HPT, "C", fontsize=14, fontweight="bold", va="top")

base = OUT + "SuppFig_PanelC_noncanonical_burden_abundance_670x350"
fig.savefig(base + ".pdf", facecolor="white")
fig.savefig(base + ".svg", facecolor="white")
fig.savefig(base + "_300dpi.png", dpi=300 * 1.0, facecolor="white")     # 236 x 123 mm @300 dpi
fig.savefig(base + "_670x350px.png", dpi=72, facecolor="white")          # 恰好 670 x 350 像素
plt.close(fig)
shutil.copy(D + "panelC_noncanonical_per_sample.csv", OUT + "PanelC_source_per_sample.csv")
shutil.copy(D + "panelC_noncanonical_pooled_DESCRIPTIVE.csv", OUT + "PanelC_trend_line_pooled_DESCRIPTIVE.csv")
print("ok")
