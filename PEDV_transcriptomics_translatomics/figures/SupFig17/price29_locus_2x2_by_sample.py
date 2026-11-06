# -*- coding: utf-8 -*-
"""Four PRICE 29-nt candidates (N-x1, N-x2, ORF1b-x13, S-x1) in one 2 x 2 figure.
Each panel: ORF schematic on top, then one P-site track per sample from top to bottom
(CHX 5 h, CHX 8 h, CHX 16 h, HAR 5 h, HAR 8 h, HAR 16 h), with bars coloured by reading frame.
Data and colours are identical to price29_locus_figs.py (single-candidate locus figures):
  P-sites = 29 nt raw + 30 nt with host-derived K correction (+2/K; frame 0 and +1 unchanged);
  frames = genome-wide reference frame map; teal = main translation frame (0, blank outside coding),
  orange-red = +1, purple = +2.
Window per candidate = longest in-frame start - 60 nt .. stop + 25 nt (same as the single-candidate figures).
Each track has its own linear y-axis set by the +1/+2 frames only; main-frame bars above it are clipped (^).
The single-candidate figures are not modified; output goes to ../figures_2x2_by_sample/.
"""
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
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Patch

B = PEDV_P+""
FB = f"{B}/Supplementary_Figure_N/FrameBiasCorrection/tables"
OUT = f"{B}/Supplementary_Figure_N/PRICE29_locus/figures_2x2_by_sample"; os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "axes.linewidth": 0.5, "pdf.fonttype": 42, "ps.fonttype": 42,
                     "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2, "ytick.major.size": 2})

CANDS = [("A", "N-x1", 26693), ("B", "N-x2", 26948), ("C", "ORF1b-x13", 18436), ("D", "S-x1", 21096)]   # stop = last nt of stop codon
ROWS = [(SR.lib("RIBO_VERO_CHX_1MOI_5h"), "CHX 5 h"), (SR.lib("RIBO_VERO_CHX_1MOI_8h"), "CHX 8 h"), (SR.lib("RIBO_VERO_CHX_1MOI_16h"), "CHX 16 h"),
        (SR.lib("RIBO_VERO_HAR_1MOI_5h"), "HAR 5 h"), (SR.lib("RIBO_VERO_HAR_1MOI_8h"), "HAR 8 h"), (SR.lib("RIBO_VERO_HAR_1MOI_16h"), "HAR 16 h")]
C0, C1, C2 = "#349789", "#DA6C54", "#695CA5"        # main frame (0), +1, +2 — same as single-candidate figures
STOPC, PRC, SHADE = "#B2182B", "#E69F00", "#F7DCC4"   # stop tick, PRICE start line, PRICE start shading
CANON = [("ORF1a", 293, 12622), ("ORF1b", 12592, 20613), ("S", 20610, 24758), ("ORF3", 24758, 25033), ("E", 25361, 25600), ("M", 25608, 26288), ("N", 26300, 27625)]

cov = pd.read_csv(f"{FB}/viral_coverage_before_after_per_position.tsv.gz", sep="\t", dtype={"read_length": str})
ref = pd.read_csv(f"{FB}/genome_reference_frame_map.tsv.gz", sep="\t").set_index("position")
M = pd.read_csv(f"{B}/orf_calling_GIb/_price_final/annotated_29nt_resultlevel_merge_of_3_samples_named.tsv", sep="\t").set_index("stop")
SEQ = "".join(l.strip() for l in open(f"{B}/refbuild_GIb/v1.1/virus/PEDV_GIb.fa") if not l.startswith(">")).upper().replace("U", "T")
STARTSET = {"ATG", "CTG", "GTG", "TTG", "ATA", "ATT", "ATC", "ACG"}; STOPS = {"TAA", "TAG", "TGA"}


def enum_starts(stop):
    """In-frame start codons (8-codon set) between the first upstream in-frame stop and this ORF's stop (stop = last nt)."""
    out, p = [], stop - 5
    while p >= 1:
        c = SEQ[p - 1:p + 2]
        if c in STOPS: break
        if c in STARTSET: out.append((p, c))
        p -= 3
    return out


G = 27959 + 3
SIG = {}
for s, _ in ROWS:
    a29 = np.zeros(G); a30 = np.zeros(G)
    for L, arr, col in (("29", a29, "raw_Psite_count"), ("30", a30, "normalised_coverage_corrected")):
        x = cov[(cov["sample"] == s) & (cov.read_length == L)]
        arr[x.position.values] = x[col].values
    SIG[s] = a29 + a30
FRAME = np.full(G, -1); CODING = np.zeros(G, bool)
FRAME[ref.index.values] = ref.reference_frame.values
CODING[ref.index.values] = ref["class"].str.startswith(("coding_unique", "overlap_assigned")).values


def nice_tick(v):
    """Largest 1/2/5 x 10^k value not above v (y tick label for each track)."""
    if v <= 1: return 1
    e = 10 ** np.floor(np.log10(v))
    for m in (5, 2, 1):
        if m * e <= v: return int(m * e)


fig = plt.figure(figsize=(7.2, 8.6))
outer = fig.add_gridspec(2, 2, left=0.105, right=0.99, top=0.965, bottom=0.125, hspace=0.20, wspace=0.27)
log = []
for k, (tag, name, stop) in enumerate(CANDS):
    r = M.loc[stop]; assert r.ORF_name == name
    ES = enum_starts(stop); big = min(p for p, _ in ES)
    lib_p = {lib: r[f"p_{lib}"] for lib in ("CHX5", "CHX8", "CHX16") if pd.notna(r[f"p_{lib}"])}
    best = min(lib_p, key=lib_p.get); pstart = int(r[f"start_{best}"])
    x0, x1 = max(1, big - 60), min(27959, stop + 25); xs = np.arange(x0, x1 + 1)
    parent = ref.loc[stop - 2, "reference_region"]; assert int(FRAME[stop - 2]) == 1, name   # ORF lies in the +1 frame
    sub = outer[k // 2, k % 2].subgridspec(8, 1, height_ratios=[1.55, 0.22] + [1] * 3 + [1] * 3, hspace=0.18)
    # spacer between CHX and HAR blocks
    axs = fig.add_subplot(sub[0]); axes = [fig.add_subplot(sub[i], sharex=axs) for i in range(2, 8)]
    # ---------- schematic ----------
    axs.set_xlim(x0, x1); axs.set_ylim(0, 1); axs.axis("off")
    axs.axvspan(pstart - 1, pstart + 3, color=SHADE, zorder=0, lw=0)
    for n, a_, b_ in CANON:
        if n == parent:
            axs.add_patch(Rectangle((max(a_, x0), 0.64), min(b_, x1) - max(a_, x0), 0.26, facecolor=C0, edgecolor="none"))
            axs.text((max(a_, x0) + min(b_, x1)) / 2, 0.77, f"canonical {n} frame", ha="center", va="center", color="white", fontsize=6.5)
    axs.plot([big, stop - 2], [0.42, 0.42], color=C1, lw=1.6, solid_capstyle="butt", zorder=2)
    for p_, c_ in ES:
        axs.plot(p_, 0.42, "^" if c_ == "ATG" else "o", ms=3.4 if c_ == "ATG" else 2.8, mfc="white", mec="k", mew=0.5, zorder=3)
    axs.plot([stop - 2, stop - 2], [0.30, 0.54], color=STOPC, lw=1.4, zorder=3)
    axs.plot([pstart, stop - 2], [0.14, 0.14], color="k", lw=1.0, solid_capstyle="butt")
    axs.plot([pstart, pstart], [0.06, 0.58], color=PRC, lw=1.2, zorder=4)
    axs.text(0.5, 1.12, name, transform=axs.transAxes, ha="center", va="bottom", fontsize=8.5, fontweight="bold")
    axs.text(-0.16, 1.12, tag, transform=axs.transAxes, ha="left", va="bottom", fontsize=11, fontweight="bold")
    # ---------- per-sample tracks ----------
    for i, (s, lab) in enumerate(ROWS):
        a = axes[i]; v = SIG[s][xs]; fr = FRAME[xs]; cd = CODING[xs]
        m0 = (fr == 0) & cd; m1 = fr == 1; m2 = fr == 2
        ymx = max(v[m1].max() if m1.any() else 0, v[m2].max() if m2.any() else 0, 1) * 1.08
        a.axvspan(pstart - 1, pstart + 3, color=SHADE, zorder=0, lw=0)
        a.vlines(xs[m0], 0, np.minimum(v[m0], ymx), color=C0, lw=0.7, zorder=2)
        a.vlines(xs[m2], 0, v[m2], color=C2, lw=0.7, zorder=2); a.vlines(xs[m1], 0, v[m1], color=C1, lw=0.7, zorder=3)
        clip = xs[m0][v[m0] > ymx]
        if len(clip): a.plot(clip, np.full(len(clip), ymx), "^", ms=1.8, color=C0, ls="", clip_on=False, zorder=4)
        a.set_ylim(0, ymx); a.set_xlim(x0, x1)
        a.set_yticks([0, nice_tick(ymx / 1.08)]); a.tick_params(axis="y", labelsize=5.5, pad=1)
        a.set_ylabel(lab, rotation=0, ha="right", va="center", fontsize=6, labelpad=3)   # sample name left of the track (no overlap with bars)
        a.spines[["top", "right"]].set_visible(False)
        if i < 5: a.tick_params(labelbottom=False)
        log.append(dict(panel=tag, candidate=name, stop_last_nt=stop, sample=s, track=lab, window=f"{x0}-{x1}",
                        PRICE_start_shown=pstart, PRICE_start_library=best, PRICE_p_that_library=lib_p[best],
                        plus1_max=float(v[m1].max() if m1.any() else 0), plus2_max=float(v[m2].max() if m2.any() else 0),
                        main_frame_max=float(v[m0].max() if m0.any() else 0), n_main_frame_clipped=int(len(clip))))
    axes[-1].ticklabel_format(axis="x", style="plain", useOffset=False); axes[-1].tick_params(axis="x", labelsize=5.5)
    axes[-1].set_xlabel("genomic coordinate (nt)", fontsize=6.5, labelpad=1)
    # move HAR block down a little (visual gap between CHX and HAR)
    for a in axes[3:]:
        pos = a.get_position(); a.set_position([pos.x0, pos.y0 - 0.006, pos.width, pos.height])

leg = [Line2D([], [], color=C0, lw=2, label="main translation frame (0)"), Line2D([], [], color=C1, lw=2, label="+1 frame"),
       Line2D([], [], color=C2, lw=2, label="+2 frame"),
       Line2D([], [], marker="^", ls="", mfc="white", mec="k", ms=4, label="AUG start"),
       Line2D([], [], marker="o", ls="", mfc="white", mec="k", ms=3.5, label="near-cognate start"),
       Line2D([], [], color=STOPC, lw=1.5, label="stop codon"),
       Line2D([], [], color="k", lw=1.2, label="PRICE ORF (CHX library with min. P)"),
       Patch(facecolor=SHADE, edgecolor=PRC, lw=1, label="PRICE start codon"),
       Line2D([], [], marker="^", ls="", color=C0, ms=3, label="main-frame bar clipped")]
fig.text(0.53, 0.068, "y axis: P-sites (29 nt raw + 30 nt K-corrected); each track has its own linear scale set by the +1/+2 frames.", ha="center", fontsize=6, color="0.3")
fig.legend(handles=leg, loc="lower center", bbox_to_anchor=(0.53, 0.0), ncol=5, frameon=False, fontsize=6, columnspacing=1.2, handletextpad=0.4)
for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}/PRICE29_candidates_2x2_by_sample.{ext}", dpi=300)
pd.DataFrame(log).to_csv(f"{OUT}/PRICE29_candidates_2x2_by_sample_tracks.tsv", sep="\t", index=False, float_format="%.5g")
print("written", OUT)
