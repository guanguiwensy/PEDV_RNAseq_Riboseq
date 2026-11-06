# -*- coding: utf-8 -*-
"""仿 Figure 5 5′ landscape 的 panel A 版式：多个窗口并排、六轨（CHX/HAR × 5/8/16 h）、黄带 = HAR 峰簇、顶部标注起始位点与 aa 长度；
着色沿用 3′ 小图方案：读框相对被标注候选自身（青绿 = 候选读框，橙 = +1，紫 = +2）。"""
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
import io, contextlib, pickle
with contextlib.redirect_stdout(io.StringIO()):
    exec(open(__import__('os').path.join(REPO_ROOT+'/riboseq/06_gRNA_terminal_translation', 'audit3utr.py'), encoding="utf-8").read().split("# ================= 表1")[0])
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.gridspec as gs
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.linewidth": 0.6, "pdf.fonttype": 42, "svg.fonttype": "none"})
b = pickle.load(open(W + r"/audit_b.pkl", "rb")); cd2 = b["cd2"]; regions = b["regions"]
FC = {0: "#2A9D8F", 1: "#E76F51", 2: "#8B5CF6"}
NAMES = [("CHX", 5), ("CHX", 8), ("CHX", 16), ("HAR", 5), ("HAR", 8), ("HAR", 16)]
SEL = ["U3-23", "U3-32", "U3-51"]            # 3′UTR 内满足「起始可归因 HAR 响应」的候选
WIN = {"U3-23": (27676, 27736), "U3-32": (27748, 27808), "U3-51": (27905, 27959)}
fig = plt.figure(figsize=(7.4, 4.6)); outer = gs.GridSpec(1, 3, figure=fig, left=0.03, right=0.985, top=0.82, bottom=0.12, wspace=0.12)
for ci, cid in enumerate(SEL):
    c = cd2[cd2.candidate_id == cid].iloc[0]; s = int(c.start); a, z = WIN[cid]; xs = np.arange(a, z + 1); fr = (xs - s) % 3
    sub = gs.GridSpecFromSubplotSpec(6, 1, subplot_spec=outer[ci], hspace=0.05)
    for k, (dr, t) in enumerate(NAMES):
        ax = fig.add_subplot(sub[k]); v = S_FORMAL[LIBS[(dr, t)]][xs]
        for f in (0, 1, 2):
            m = fr == f
            if m.any(): ax.bar(xs[m], v[m], width=1.0, color=FC[f], lw=0)
        for (ra, rb) in regions:
            if rb >= a and ra <= z: ax.axvspan(ra - 0.5, rb + 0.5, color="#F7EBC3", lw=0, zorder=0)
        ax.axvline(s, color="#555555", lw=0.5, ls="--")
        ax.set_xlim(a - 0.5, z + 0.5); ax.set_yticks([]); ax.spines["left"].set_visible(False); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_color("#888888"); ax.tick_params(labelbottom=(k == 5), length=2 if k == 5 else 0, labelsize=7)
        ax.text(0.985, 0.78, f"{dr} {t} h", transform=ax.transAxes, ha="right", va="center", fontsize=7.5)
        ax.text(0.015, 0.78, f"max {int(v.max())}", transform=ax.transAxes, ha="left", va="center", fontsize=5.8, color="#666666")
        if k == 0:
            tr = ax.get_xaxis_transform()
            lab = f"{s}" + chr(10) + f"({int(c.len_aa)} aa)" + ((chr(10) + "unclosed") if not c.closed else "")
            ax.annotate(lab, xy=(s, 0.97), xycoords=tr, xytext=(s, 1.28), textcoords=tr,
                        ha="center", va="bottom", fontsize=7.5, color="#9C5A00", annotation_clip=False, arrowprops=dict(arrowstyle="-", color="#D62728", lw=0.9))
    ax.set_xlabel(f"{c.start_codon}  ·  {cid}", fontsize=6.5, labelpad=1)
fig.text(0.015, 0.975, "A", fontsize=13, fontweight="bold", va="top")
fig.text(0.5, 0.975, "PEDV 3′UTR candidate start windows (start-attributable HAR response)", fontsize=8.5, ha="center", va="top")
from matplotlib.patches import Patch
fig.legend(handles=[Patch(color=FC[0], label="candidate frame"), Patch(color=FC[1], label="+1"), Patch(color=FC[2], label="+2"), Patch(color="#F7EBC3", label="HAR peak cluster")],
           loc="lower center", ncol=4, fontsize=6.5, frameon=False, bbox_to_anchor=(0.5, 0.0), handlelength=1.1)
for ext, kw in (("png", dict(dpi=300)), ("pdf", {}), ("svg", {})): fig.savefig(OUT + rf"/figures/Fig5_3UTR_candidate_landscape_panelA_style.{ext}", facecolor="white", **kw)
print("ok")
