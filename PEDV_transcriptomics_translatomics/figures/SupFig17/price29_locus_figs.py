# -*- coding: utf-8 -*-
"""Locus figures for every ORF with min PRICE p < 0.05 (3-sample result-level merge, 29 nt PRICE runs).
Data: 29 nt P-sites raw (frame uncorrected) + 30 nt P-sites with host-derived K correction (+2 / K, from FrameBiasCorrection; frame 0 and +1 unchanged).
Frames = piecewise genome-wide reference frame (FrameBiasCorrection/genome_reference_frame_map); Fig. 6D colour scheme: teal = main translation frame (0, blank outside coding), orange-red = +1, purple = +2.
Six whitelisted Vero Ribo-seq libraries: CHX 5/8/16 h, HAR 5/8/16 h. CHX rows: all three frames on one linear axis. HAR rows: y-axis set by +1/+2 frames only; canonical clipped (^)."""
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
import os, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
B = PEDV_P+""
FB = f"{B}/Supplementary_Figure_N/FrameBiasCorrection/tables"
O = f"{B}/Supplementary_Figure_N/PRICE29_locus"; FG = f"{O}/figures"; TB = f"{O}/tables"
os.makedirs(FG, exist_ok=True); os.makedirs(TB, exist_ok=True)
cov = pd.read_csv(f"{FB}/viral_coverage_before_after_per_position.tsv.gz", sep="\t", dtype={"read_length": str})
ref = pd.read_csv(f"{FB}/genome_reference_frame_map.tsv.gz", sep="\t").set_index("position")
M = pd.read_csv(f"{B}/orf_calling_GIb/_price_final/annotated_29nt_resultlevel_merge_of_3_samples_named.tsv", sep="\t")
M = M[M.min_p < 0.05].sort_values("min_p").reset_index(drop=True)
ROWS = [(SR.lib("RIBO_VERO_CHX_1MOI_5h"), "CHX", 5, "CHX5"), (SR.lib("RIBO_VERO_CHX_1MOI_8h"), "CHX", 8, "CHX8"), (SR.lib("RIBO_VERO_CHX_1MOI_16h"), "CHX", 16, "CHX16"),
        (SR.lib("RIBO_VERO_HAR_1MOI_5h"), "HAR", 5, None), (SR.lib("RIBO_VERO_HAR_1MOI_8h"), "HAR", 8, None), (SR.lib("RIBO_VERO_HAR_1MOI_16h"), "HAR", 16, None)]
def seq_fa(p=B + "/refbuild_GIb/v1.1/virus/PEDV_GIb.fa"):
    return "".join(l.strip() for l in open(p) if not l.startswith(">")).upper().replace("U", "T")
SEQ = seq_fa()
PRX = pd.read_csv(f"{B}/ribosome_pause_alt_start/v2_locus_atlas/tables/altstart_RHO_proximity_by_start.tsv", sep="	")
HARPOS = {(int(a), int(b)): bool(c) for a, b, c in zip(PRX.start, PRX.stop, PRX.HAR_positive_start)}
FAMNAME = {(int(a), int(b)): n for a, b, n in zip(PRX.start, PRX.stop, PRX.family)}
sg = pd.read_csv(f"{B}/GIb_rnaseq_reanalysis/junction_atlas/sgRNA_atlas_GIb.tsv", sep="	", comment="#")
SGP = [int(x) for x in sg[sg.total_reads >= 100].body_cluster_start]
motif = [i for i in range(1, len(SEQ) - 6) if SEQ[i-1:i+2] in ("AAA", "TTT", "GGG", "CCC") and SEQ[i+2:i+5] in ("AAA", "TTT") and SEQ[i+5] != "G"]
STARTSET = {"ATG", "CTG", "GTG", "TTG", "ATA", "ATT", "ATC", "ACG"}; STOPS = {"TAA", "TAG", "TGA"}


def enum_starts(stop):
    """all in-frame start codons (8-codon set) between the first upstream in-frame stop and this ORF's stop; stop = last nt of the stop codon (1-based)"""
    out = []; p = stop - 5
    while p >= 1:
        c = SEQ[p - 1:p + 2]
        if c in STOPS: break
        if c in STARTSET: out.append((p, c))
        p -= 3
    return out


G = 27959 + 3
SIG = {}
for s, dr, t, _ in ROWS:
    a29 = np.zeros(G); a30 = np.zeros(G)
    for L, arr, col in (("29", a29, "raw_Psite_count"), ("30", a30, "normalised_coverage_corrected")):
        x = cov[(cov["sample"] == s) & (cov.read_length == L)]
        arr[x.position.values] = x[col].values
    SIG[s] = a29 + a30
FRAME = np.full(G, -1); CODING = np.zeros(G, bool)
FRAME[ref.index.values] = ref.reference_frame.values
CODING[ref.index.values] = ref["class"].str.startswith(("coding_unique", "overlap_assigned")).values
C0, C1, C2 = "#349789", "#DA6C54", "#695CA5"   # Fig. 6D scheme: main translation frame (0) teal, +1 orange-red, +2 purple
LC = {"CHX5": "#9DC3E6", "CHX8": "#5B9BD5", "CHX16": "#1F4E8C"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.linewidth": 0.6})
safe = lambda n: re.sub(r"[^A-Za-z0-9_.-]", "_", n.replace("+", "p").replace("#", "n").replace("'", ""))
log = []
for _, r in M.iterrows():
    stop = int(r.stop); starts = {t: int(r[f"start_{t}"]) for t in ("CHX5", "CHX8", "CHX16") if pd.notna(r[f"start_{t}"])}
    ES = enum_starts(stop); allst = sorted(set([p for p, _ in ES]) | set(starts.values()))
    big = min(allst); bigc = SEQ[big - 1:big + 2]; big_aa = (stop - big + 1) // 3 - 1
    x0 = max(1, big - 60); x1 = min(27959, stop + 25)
    xs = np.arange(x0, x1 + 1)
    orf_frame = int(FRAME[stop - 2]); cls = ref.loc[stop - 2, "reference_region"] if (stop - 2) in ref.index else "?"
    fig = plt.figure(figsize=(11, 12.3))
    gs = fig.add_gridspec(7, 1, height_ratios=[1] * 6 + [0.95], hspace=0.28, left=0.1, right=0.97, top=0.872, bottom=0.095)
    axs = [fig.add_subplot(gs[i]) for i in range(7)]
    for ai, (s, dr, t, key) in enumerate(ROWS):
        a = axs[ai]; v = SIG[s][xs]; fr = FRAME[xs]; cd = CODING[xs]
        m0 = (fr == 0) & cd; m1 = fr == 1; m2 = fr == 2
        ymx = max(v[m1].max() if m1.any() else 0, v[m2].max() if m2.any() else 0, 1) * 1.08
        a.vlines(xs[m0], 0, np.minimum(v[m0], ymx), color=C0, lw=0.8)
        a.vlines(xs[m2], 0, v[m2], color=C2, lw=0.8); a.vlines(xs[m1], 0, v[m1], color=C1, lw=0.8)
        a.set_ylim(0, ymx); a.set_xlim(x0, x1)
        if m0.any() and (v[m0] > ymx).any():
            cl = xs[m0][v[m0] > ymx]; a.plot(cl, np.full(len(cl), ymx), "^", ms=3, color=C0, ls="", clip_on=False)
            a.text(0.995, 0.97, f"main-frame max {int(v[m0].max()):,} (clipped, ^)", transform=a.transAxes, ha="right", va="top", fontsize=6.5, color=C0)
        a.axvline(stop - 2, color="k", lw=0.6, alpha=0.5)
        for tk, sx in starts.items():
            a.axvline(sx, color=LC[tk], lw=0.7, ls=":", alpha=0.9)
        a.set_ylabel(f"{dr} {t} h\nP-site count", fontsize=8)
        if key and not np.isnan(r[f"p_{key}"]):
            a.text(0.005, 0.95, f"PRICE p = {r[f'p_{key}']:.3g} (start {starts[key]})", transform=a.transAxes, ha="left", va="top", fontsize=7, color="#444444")
        elif key:
            a.text(0.005, 0.95, "PRICE: not called in this library", transform=a.transAxes, ha="left", va="top", fontsize=7, color="#888888")
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
        a.tick_params(labelbottom=False)
        log.append(dict(orf=r.ORF_name, stop=stop, sample=s, window=f"{x0}-{x1}", reads_in_window=float(v.sum()), plus1_max=float(v[m1].max() if m1.any() else 0), plus2_max=float(v[m2].max() if m2.any() else 0)))
    a = axs[6]; a.set_xlim(x0, x1); a.set_ylim(-0.3, 5.6)
    Y_ST, Y_PR, Y_STOP, Y_SLIP, Y_SG = 4.9, 3.9, 2.9, 2.0, 1.1
    nhar = 0; fam_here = any((p, stop - 2) in HARPOS for p, _ in ES)
    for p_, c_ in ES:
        har = HARPOS.get((p_, stop - 2), False); nhar += har; col = "#e07a1f" if har else "#7f7f7f"
        if c_ == "ATG": a.plot(p_, Y_ST, "^", color=col, ms=4.5, mew=0, clip_on=False)
        else: a.plot(p_, Y_ST, "o", mfc="white", mec=col, ms=3.8, mew=0.7, clip_on=False)
    for tk, sx in starts.items():
        a.plot(sx, Y_PR, "D", color=LC[tk], ms=4.5, mew=0, clip_on=False)
    a.plot([big, big], [Y_ST - 0.35, Y_ST + 0.35], color="k", lw=0.8, clip_on=False)
    a.plot([stop - 2, stop - 2], [Y_STOP - 0.3, Y_STOP + 0.3], color="black", lw=1.2, clip_on=False); a.text(stop - 6, Y_STOP, "stop of this ORF", ha="right", va="center", fontsize=7)
    for m in motif:
        if x0 <= m <= x1 - 6: a.plot([m, m + 6], [Y_SLIP, Y_SLIP], color="#d55e00", lw=2.2, solid_capstyle="butt", clip_on=False)
    for c0 in SGP:
        if x0 <= c0 <= x1: a.plot([c0, c0], [Y_SG - 0.25, Y_SG + 0.25], color="#9a9a9a", lw=0.6, clip_on=False)
    a.set_yticks([Y_ST, Y_PR, Y_STOP, Y_SLIP, Y_SG]); a.set_yticklabels(["starts", "PRICE start", "stop", "-1 slippery", "sgRNA 5'"], fontsize=7)
    a.ticklabel_format(axis="x", style="plain", useOffset=False)
    for sp in ("top", "right", "left"):
        a.spines[sp].set_visible(False)
    a.set_xlabel("PEDV_GIb genome position (nt)")
    pl = ", ".join(f"{k}: {r[f'p_{k}']:.3g}" if pd.notna(r[f"p_{k}"]) else f"{k}: not called" for k in ("CHX5", "CHX8", "CHX16"))
    fl = {0: "canonical-frame (0)", 1: "+1 frame", 2: "+2 frame"}.get(orf_frame, "?")
    fig.text(0.1, 0.972, f"{r.ORF_name}  (stop {stop}; largest ORF: start {big} {bigc}, {big_aa} aa; {len(ES)} candidate starts, {nhar} HAR-positive)", fontsize=10, fontweight="bold", ha="left")
    fig.text(0.1, 0.955, f"ORF lies in the {fl} of the reference region {cls}", fontsize=8, ha="left", color="#444444")
    fig.text(0.1, 0.938, f"PRICE 29-nt-only p: {pl};  min p = {r.min_p:.3g};  libraries with p<0.05: {int(r['n_libs_p<0.05'])}", fontsize=8, ha="left", color="#444444")
    fig.text(0.1, 0.921, "P-sites: 29 nt raw + 30 nt with host-derived K correction (+2/K; frame 0 and +1 unchanged). Raw counts for 29 nt; corrected values are normalised coverage, not read counts.", fontsize=7, ha="left", color="#444444")
    fig.text(0.1, 0.905, "Each row has its own linear y-scale. All rows: y-axis set by the +1/+2 frames only; main-frame (teal) bars above it are clipped (^, true maximum printed).", fontsize=7, ha="left", color="#444444")
    fig.text(0.1, 0.889, ("HAR-positive status exists only for Fig. 6 +1-family starts; all other starts are drawn grey (not assessed). " if not fam_here else "HAR-positive = frozen Fig. 6 / v2 atlas call.") + " Black tick = longest start; dotted line = PRICE start per CHX library.", fontsize=7, ha="left", color="#444444")
    fig.legend(handles=[Line2D([], [], color=C0, lw=1.6, label="main translation frame (0)"), Line2D([], [], color=C1, lw=1.6, label="+1 frame"), Line2D([], [], color=C2, lw=1.6, label="+2 frame"),
                        Line2D([], [], marker="^", color="#e07a1f", ls="", ms=5, label="AUG, HAR-positive"), Line2D([], [], marker="^", color="#7f7f7f", ls="", ms=5, label="AUG, not HAR-positive / not assessed"),
                        Line2D([], [], marker="o", mfc="white", mec="#7f7f7f", ls="", ms=4, label="near-cognate start (orange = HAR-positive)"),
                        Line2D([], [], marker="D", color=LC["CHX5"], ls="", ms=4.5, label="PRICE start CHX 5 h"), Line2D([], [], marker="D", color=LC["CHX8"], ls="", ms=4.5, label="PRICE start CHX 8 h"), Line2D([], [], marker="D", color=LC["CHX16"], ls="", ms=4.5, label="PRICE start CHX 16 h")],
               loc="lower center", ncol=3, frameon=False, fontsize=7.5)
    fig.savefig(f"{FG}/locus_{safe(r.ORF_name)}_stop{stop}.png", dpi=170); plt.close(fig)
pd.DataFrame(log).to_csv(f"{TB}/locus_window_summary.tsv", sep="\t", index=False, float_format="%.5g")
M[["stop", "ORF_name", "min_p", "n_libs_p<0.05", "start_CHX5", "start_CHX8", "start_CHX16", "p_CHX5", "p_CHX8", "p_CHX16"]].to_csv(f"{TB}/orfs_plotted.tsv", sep="\t", index=False, float_format="%.5g")
print("made", len(M))
