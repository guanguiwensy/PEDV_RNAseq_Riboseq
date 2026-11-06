# -*- coding: utf-8 -*-
"""5′ 端三个窗口（3–45、75–135、255–325）× 六轨，仿用户给的 gRNA 5′ translation landscape 面板 A；着色沿用 3′ 小图方案：
读框相对 ORF1a 读框（青绿 = ORF1a 框，橙 = +1，紫 = +2）。P-site = 正式 6 库（CHX 16 h = RIBO_VERO_CHX_1MOI_16h）。"""
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
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    exec(open(__import__('os').path.join(REPO_ROOT+'/riboseq/06_gRNA_terminal_translation', 'audit3utr.py'), encoding="utf-8").read().split("# ================= 表1")[0])
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.gridspec as gs
from matplotlib.patches import Patch
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.linewidth": 0.6, "pdf.fonttype": 42, "svg.fonttype": "none"})
FC = {0: "#2A9D8F", 1: "#E76F51", 2: "#8B5CF6"}
NAMES = [("CHX", 5), ("CHX", 8), ("CHX", 16), ("HAR", 5), ("HAR", 8), ("HAR", 16)]
MST = 293
HOT = [(21, 23), (102, 104), (111, 114), (270, 272), (285, 296)]
cat = pd.read_csv(ROOT + r"/fig5_gRNA5prime_figure_package/input_data/gRNA5p_catalogue/gRNA_5prime_candidate_catalogue.tsv", sep="\t", comment="#")
cat = cat[cat.length_aa >= 8]
WINS = [(3, 45, 22, 22, "AUG", 3), (75, 135, 113, 71, "UUG", 20), (255, 325, 271, 271, "UUG", 9)]   # (窗口起, 窗口止, 热点位置, 候选起始, 密码子, aa)
fig = plt.figure(figsize=(7.4, 4.6)); outer = gs.GridSpec(1, 3, figure=fig, left=0.03, right=0.985, top=0.82, bottom=0.12, wspace=0.12)
for ci, (a, z, hp, cs, ccod, caa) in enumerate(WINS):
    xs = np.arange(a, z + 1); fr = (xs - cs) % 3
    sub = gs.GridSpecFromSubplotSpec(6, 1, subplot_spec=outer[ci], hspace=0.05)
    lab = f"{hp}" + chr(10) + f"{ccod}@{cs} ({caa} aa)"
    for k, (dr, t) in enumerate(NAMES):
        ax = fig.add_subplot(sub[k]); v = S_FORMAL[LIBS[(dr, t)]][xs]
        for f in (0, 1, 2):
            m = fr == f
            if m.any(): ax.bar(xs[m], v[m], width=1.0, color=FC[f], lw=0)
        for (ra, rb) in HOT:
            if rb >= a and ra <= z: ax.axvspan(ra - 0.5, rb + 0.5, color="#F7EBC3", lw=0, zorder=0)
        ax.set_xlim(a - 0.5, z + 0.5); ax.set_yticks([])
        for sp in ("left", "top", "right"): ax.spines[sp].set_visible(False)
        ax.spines["bottom"].set_color("#888888"); ax.tick_params(labelbottom=(k == 5), length=2 if k == 5 else 0, labelsize=7)
        ax.text(0.985, 0.78, f"{dr} {t} h", transform=ax.transAxes, ha="right", va="center", fontsize=7.5)
        ax.text(0.015, 0.78, f"max {int(v.max())}", transform=ax.transAxes, ha="left", va="center", fontsize=5.8, color="#666666")
        if k == 0:
            tr = ax.get_xaxis_transform()
            ax.annotate(lab, xy=(hp, 0.97), xycoords=tr, xytext=(hp, 1.28), textcoords=tr, ha="center", va="bottom", fontsize=7.2, color="#9C5A00",
                        annotation_clip=False, arrowprops=dict(arrowstyle="-", color="#D62728", lw=0.9))
fig.text(0.015, 0.975, "A", fontsize=13, fontweight="bold", va="top")
fig.text(0.5, 0.975, "gRNA 5′ translation landscape", fontsize=8.5, ha="center", va="top")
fig.legend(handles=[Patch(color=FC[0], label="candidate ORF frame (frame 0)"), Patch(color=FC[1], label="frame +1"), Patch(color=FC[2], label="frame +2"), Patch(color="#F7EBC3", label="recurrent peak region")],
           loc="lower center", ncol=4, fontsize=6.5, frameon=False, bbox_to_anchor=(0.5, 0.0), handlelength=1.1)
for ext, kw in (("png", dict(dpi=300)), ("pdf", {}), ("svg", {})): fig.savefig(OUT + rf"/figures/Fig5_5prime_landscape_panelA_candframe.{ext}", facecolor="white", **kw)
print("ok")
