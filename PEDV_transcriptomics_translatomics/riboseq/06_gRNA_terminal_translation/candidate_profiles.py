# -*- coding: utf-8 -*-
"""3′ 区每个候选 ORF 一张六轨 P-site 小图（仿 FigS5K_framecolor：按读框着色；读框相对候选自身起始，0 = 候选读框）。"""
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
import io, contextlib, pickle, os
with contextlib.redirect_stdout(io.StringIO()):
    exec(open(__import__('os').path.join(__import__('os').path.dirname(__REPO_FILE__), 'audit3utr.py'), encoding="utf-8").read().split("# ================= 表1")[0])
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.gridspec as gs
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.linewidth": 0.5, "pdf.fonttype": 42})
b = pickle.load(open(W + r"/audit_b.pkl", "rb")); cd2 = b["cd2"]; regions = b["regions"]
FC = {0: "#2A9D8F", 1: "#E76F51", 2: "#8B5CF6"}
LABC = {"CHX": "#0A4F7D", "HAR": "#9C3A00"}
OUTD = OUT + r"/figures/candidate_profiles"; os.makedirs(OUTD, exist_ok=True)
NAMES = [("CHX", 5), ("CHX", 8), ("CHX", 16), ("HAR", 5), ("HAR", 8), ("HAR", 16)]
from matplotlib.backends.backend_pdf import PdfPages
PP = PdfPages(OUT + "/figures/3UTR_all_candidate_profiles.pdf")
for c in cd2.sort_values("start").itertuples():
    s, e = c.start, c.stop_end
    a, z = max(s - 30, 1), min(e + 20, L)
    xs = np.arange(a, z + 1); fr = (xs - s) % 3
    fig = plt.figure(figsize=(6.2, 4.4))
    sub = gs.GridSpec(7, 1, figure=fig, hspace=0.12, height_ratios=[1] * 6 + [0.7], left=0.04, right=0.985, top=0.86, bottom=0.09)
    for k, (dr, t) in enumerate(NAMES):
        ax = fig.add_subplot(sub[k]); v = S_FORMAL[LIBS[(dr, t)]][xs]
        for f in (0, 1, 2):
            m = fr == f
            if m.any(): ax.bar(xs[m], v[m], width=1.0, color=FC[f], lw=0)
        for (ra, rb) in regions:
            if rb >= a and ra <= z: ax.axvspan(ra - 0.5, rb + 0.5, color="#F3E3A0", alpha=0.5, lw=0, zorder=0)
        ax.axvline(s, color="#555555", lw=0.6, ls="--"); ax.axvline(e, color="#555555", lw=0.6, ls=":")
        ax.set_xlim(a - 0.5, z + 0.5); ax.set_yticks([]); ax.tick_params(labelbottom=False, length=0)
        for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
        ax.text(0.01, 0.78, f"{dr} {t} h  max {int(v.max())}", transform=ax.transAxes, fontsize=6, color=LABC[dr])
    ax = fig.add_subplot(sub[6]); ax.set_xlim(a - 0.5, z + 0.5); ax.set_ylim(0, 1)
    for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
    ax.set_yticks([]); ax.tick_params(labelsize=6.5)
    ax.add_patch(plt.Rectangle((s, 0.15), e - s + 1, 0.3, fc="#9CC9E6" if c.start_codon != "AUG" else "#0072B2", ec="#0072B2", lw=0.6, ls="-" if c.closed else (0, (2, 1.5))))
    if c.start_codon == "AUG": ax.plot([s], [0.3], marker="D", ms=4.5, mfc="#222222", mec="k", mew=0.4, zorder=5, clip_on=False)
    else: ax.plot([s], [0.3], marker="o", ms=4.5, mfc="white", mec="#222222", mew=1.0, zorder=5, clip_on=False)
    ax.set_xlabel("genome position (nt)", fontsize=7, labelpad=1)
    flags = []
    if c.n_start_attributable_libs > 0: flags.append(f"start-attributable HAR in {int(c.n_start_attributable_libs)} lib")
    if c.n_CHX_candidate_periodicity_libs > 0: flags.append(f"CHX periodic (top-1 removed) in {int(c.n_CHX_candidate_periodicity_libs)} lib")
    if c.start_region == "N_terminal60" and c.frame_rel_N == 0: flags.append("in frame with N CDS")
    if not c.closed: flags.append("unclosed")
    fig.suptitle(f"{c.candidate_id}  {c.start_codon}@{s}–{e}  {c.len_aa} aa  [{c.start_region}]", x=0.04, ha="left", y=0.985, fontsize=8.5, fontweight="bold")
    fig.text(0.04, 0.915, "; ".join(flags) if flags else "no start-attributable HAR response; no CHX candidate periodicity", fontsize=6.5, color="#444444", ha="left")
    fig.legend(handles=[Patch(color=FC[0], label="frame 0 (candidate)"), Patch(color=FC[1], label="+1"), Patch(color=FC[2], label="+2"),
                        Line2D([0], [0], color="#555555", lw=0.8, ls="--", label="start"), Line2D([0], [0], color="#555555", lw=0.8, ls=":", label="stop end"), Patch(color="#F3E3A0", label="HAR peak cluster")],
               loc="upper right", bbox_to_anchor=(0.985, 0.955), ncol=6, fontsize=5.8, frameon=False, handlelength=1.1, columnspacing=0.8)
    fig.savefig(OUTD + rf"/{c.candidate_id}_{c.start_codon}_{s}.png", dpi=150, facecolor="white")
    fig.savefig(OUTD + rf"/{c.candidate_id}_{c.start_codon}_{s}.pdf", facecolor="white"); PP.savefig(fig, facecolor="white"); plt.close(fig)
PP.close(); print("done", len(cd2))
