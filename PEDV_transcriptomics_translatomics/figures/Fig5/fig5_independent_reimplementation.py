"""Fig 5 (gRNA 5' / 3' translation landscape) -- INDEPENDENT REIMPLEMENTATION written after the final sample-table correction.

Status: NOT the program that drew the submitted Fig. 5 (that script is withheld, see docs/EXCLUDED_FROM_PUBLIC_PACKAGE.md) and NOT a pixel copy of it.
It redraws the four panels from tables that were recomputed with the whitelisted libraries only (final sample table, config/sample_registry.py):

  A / C  tracks      RNA-seq depth = mean per-position depth of the 20 whitelisted Vero RNA-seq libraries (log10(depth+1));
                     Ribo-seq tracks = P-site CPM of the six whitelisted Vero Ribo-seq libraries (strict_core read lengths, + strand,
                     CPM = P-sites at the position / all strict_core + strand PEDV P-sites of the same library x 1e6; this is the definition of the
                     archived 'gRNA_5prime_HAR_CHX_positional_CPM' table and was verified against it), y-axis = sqrt(CPM).
  B                  5' candidate ORFs (gRNA_5prime_candidate_evidence.tsv from g1/g2): line = start..stop, marker = start codon type.
  D                  3' candidate ORFs (3UTR_candidate_ORFs.tsv from audit3utr_b/c).

DIFFERENCES FROM THE SUBMITTED FIGURE (see tests/FINAL_SAMPLE_WHITELIST_AUDIT.tsv and docs/SAMPLE_WHITELIST_DIFF_REPORT.md):
  * the submitted y-scale of panels A/C ('sqrt CPM' up to ~500) cannot be reconstructed from the archived data with the CPM definition above
    (maximum sqrt(CPM) in the 5' region is ~100); this script plots the archived definition;
  * the red ('ORF with signal') / blue ('without signal') rule of the submitted figure could not be recovered from the archived tables. The rules used
    here are explicit constants below (RULE_5P, RULE_3P) and every call is written to the source-data table, so the authors can change them;
  * CHX 16 h is the library RIBO_VERO_CHX_1MOI_16h of the final sample table (the earlier tables used a library that is not in it).
Usage: python fig5_independent_reimplementation.py [OUTDIR]
"""
import os, sys, csv, json
_d = os.path.dirname(os.path.abspath(__file__))
while _d != os.path.dirname(_d) and not os.path.isdir(os.path.join(_d, "config")):
    _d = os.path.dirname(_d)
sys.path.insert(0, os.path.join(_d, "config"))
from pedv_paths import PEDV_P
import sample_registry as SR
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

OUT = sys.argv[1] if len(sys.argv) > 1 else "out_fig5_independent"
os.makedirs(OUT, exist_ok=True)
GLEN = 27959
# explicit calling rules (editable)
RULE_5P = "red if a jackpot-robust CHX body signal (CHX_body_signal_deTop1) is present at >=1 of 5/8/16 h"
RULE_3P = "red if >=1 whitelisted library shows a start-attributable HAR response (n_start_attributable_libs>=1) or CHX candidate-frame periodicity in >=2 libraries"

# ------------------------------------------------------------------ P-site CPM tracks (whitelisted Ribo-seq libraries only)
RGP = PEDV_P + "riboseq_GIb_mapping_psite/"
ps = pd.read_csv(RGP + "psite_counts_PEDV_GIb.tsv.gz", sep="\t", comment="#", low_memory=False)
rls = pd.read_csv(RGP + "readlength_sets_GIb.tsv", sep="\t", comment="#")
SR.assert_exactly(set(ps.sample_id), SR.RIBO_VERO_ALIASES, "P-site table")
ps = ps[ps.strand == "+"]
sc = rls[rls.in_strict_core.astype(str) == "True"][["sample_id", "read_length"]].copy()
sc["read_length"] = sc.read_length.astype(int)
m = ps.merge(sc, on=["sample_id", "read_length"])
CPM = {}
for al in SR.RIBO_VERO_ALIASES:
    s = SR.lib(al)
    g = m[m.sample_id == s].groupby("position").psite_count.sum()
    a = np.zeros(GLEN + 2)
    a[g.index.values.astype(int)] = g.values
    CPM[al] = a / g.sum() * 1e6

# ------------------------------------------------------------------ RNA depth (20 whitelisted Vero RNA-seq libraries)
RNA = np.zeros(GLEN + 2)
for al in SR.RNA_VERO_ALIASES:
    d = pd.read_csv(PEDV_P + "GIb_rnaseq_reanalysis/coverage/%s.depth.tsv" % SR.lib(al), sep="\t", header=None)
    RNA[d[1].values] += d[2].values
RNA /= len(SR.RNA_VERO_ALIASES)
LOGRNA = np.log10(RNA + 1)

# ------------------------------------------------------------------ candidates
E5 = pd.read_csv(PEDV_P + "fig5_terminal_ORF/gRNA_5prime/PEDV_gRNA_5prime_candidate_evidence.tsv", sep="\t", comment="#", keep_default_na=False)
flag5 = [c for c in E5.columns if c.startswith("exploratory_prev_CHX_body_signal_deTop1_")]
E5["signal_call"] = E5[flag5].astype(str).apply(lambda r: any(v == "True" for v in r), axis=1)
E5["start"] = E5.start.astype(int); E5["stop"] = E5["stop"].astype(int)
C3 = pd.read_csv(PEDV_P + "Fig5_3UTR_translation_audit/tables/3UTR_candidate_ORFs.tsv", sep="\t", comment="#", keep_default_na=False)
C3["start"] = C3.start.astype(int); C3["stop_end"] = C3.stop_end.astype(int)
C3["signal_call"] = (C3.n_start_attributable_libs.astype(float) >= 1) | (C3.n_CHX_candidate_periodicity_libs.astype(float) >= 2)

# ------------------------------------------------------------------ drawing helpers
plt.rcParams.update({"font.family": ["Arial", "DejaVu Sans"], "pdf.fonttype": 42, "svg.fonttype": "none", "font.size": 8})
GREY = {5: "#c8c8c8", 8: "#8c8c8c", 16: "#000000"}
ORNG = {5: "#f4b183", 8: "#e07b39", 16: "#843c0c"}


def tracks(axs, x0, x1):
    xs = np.arange(x0, x1 + 1)
    axs[0].fill_between(xs, 0, LOGRNA[x0:x1 + 1], color="#9a9a9a", lw=0)
    axs[0].set_ylabel("RNA-seq\nlog$_{10}$ depth", rotation=0, ha="right", va="center")
    for t in (5, 8, 16):
        axs[1].plot(xs, np.sqrt(CPM["RIBO_VERO_CHX_1MOI_%dh" % t][x0:x1 + 1]), color=GREY[t], lw=0.8, label="%d h" % t)
        axs[2].plot(xs, np.sqrt(CPM["RIBO_VERO_HAR_1MOI_%dh" % t][x0:x1 + 1]), color=ORNG[t], lw=0.8, label="%d h" % t)
    axs[1].set_ylabel("CHX\nsqrt CPM", rotation=0, ha="right", va="center"); axs[2].set_ylabel("HAR\nsqrt CPM", rotation=0, ha="right", va="center")
    axs[1].legend(frameon=False, ncol=3, loc="upper right", fontsize=7); axs[2].legend(frameon=False, ncol=3, loc="upper right", fontsize=7)
    for a in axs:
        a.set_xlim(x0, x1)
        for sp in ("top", "right"): a.spines[sp].set_visible(False)
    for a in axs[:-1]: a.tick_params(labelbottom=False)


def lanes(df, s, e):
    """greedy lane packing of intervals"""
    ends = []; out = []
    for _, r in df.sort_values([s, e]).iterrows():
        for i, en in enumerate(ends):
            if r[s] > en + 2:
                ends[i] = r[e]; out.append(i); break
        else:
            ends.append(r[e]); out.append(len(ends) - 1)
    d = df.sort_values([s, e]).copy(); d["lane"] = out
    return d


def orfs(ax, df, s, e, codon, x0, x1, ymax_label=None):
    d = lanes(df[(df[s] <= x1) & (df[e] >= x0)], s, e)
    for _, r in d.iterrows():
        y = -r.lane
        col = "#d62728" if r.signal_call else "#1f77b4"
        ax.plot([max(r[s], x0), min(r[e], x1)], [y, y], color=col, lw=2.2, solid_capstyle="butt")
        if x0 <= r[s] <= x1:
            if r[codon] == "AUG": ax.plot(r[s], y, marker="D", ms=4, color="#444444", ls="")
            else: ax.plot(r[s], y, marker="o", ms=4, mfc="white", mec="#444444", ls="")
    ax.set_xlim(x0, x1); ax.set_ylim(-(d.lane.max() + 1) if len(d) else -1, 1)
    for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
    ax.set_yticks([])


fig = plt.figure(figsize=(8.6, 11.5))
gs = fig.add_gridspec(8, 1, height_ratios=[0.9, 1.1, 1.1, 2.2, 0.9, 1.1, 1.1, 2.4], hspace=0.35)
# --- A
axA = [fig.add_subplot(gs[i]) for i in range(3)]
tracks(axA, 1, 443)
axA[0].set_title("gRNA 5' translation landscape (independent reimplementation; whitelisted libraries)", fontsize=9)
axA[2].set_xlabel("gRNA genomic coordinate (nt)")
for ax_ in axA:
    ax_.axvspan(283, 300, color="#fdebc0", lw=0, zorder=0)
# --- B
axB = fig.add_subplot(gs[3])
orfs(axB, E5, "start", "stop", "codon", 1, 443)
axB.set_xlabel("gRNA genomic coordinate (nt)")
axB.text(0.01, 0.02, "marker: diamond = AUG start, open circle = near-cognate start; red = ORF with signal (rule: %s)" % RULE_5P, transform=axB.transAxes, fontsize=5.5, color="#444")
# --- C
axC = [fig.add_subplot(gs[4 + i]) for i in range(3)]
tracks(axC, 27329, 27959)
axC[0].set_title("gRNA 3' translation landscape", fontsize=9)
axC[2].set_xlabel("gRNA genomic coordinate (nt)")
for ax_ in axC:
    ax_.axvline(27625.5, color="#777", ls="--", lw=0.6)
# --- D
axD = fig.add_subplot(gs[7])
orfs(axD, C3, "start", "stop_end", "start_codon", 27329, 27959)
axD.set_xlabel("gRNA genomic coordinate (nt)")
axD.text(0.01, 0.02, "red = ORF with signal (rule: %s)" % RULE_3P, transform=axD.transAxes, fontsize=5.5, color="#444")
for lab, a in zip("ABCD", (axA[0], axB, axC[0], axD)):
    a.text(-0.12, 1.05, lab, transform=a.transAxes, fontsize=13, fontweight="bold")
fig.subplots_adjust(left=0.14, right=0.97, top=0.96, bottom=0.05)
for ext in ("pdf", "svg", "png"):
    fig.savefig(os.path.join(OUT, "Fig5_independent." + ext), dpi=300)

# ------------------------------------------------------------------ source data
E5[["candidate_id", "start", "codon", "stop", "length_aa", "class", "signal_call"] + flag5].to_csv(os.path.join(OUT, "Fig5B_candidates_signal_calls.csv"), index=False)
C3[["candidate_id", "start", "start_codon", "stop_end", "len_aa", "n_start_attributable_libs", "n_CHX_candidate_periodicity_libs", "signal_call"]].to_csv(os.path.join(OUT, "Fig5D_candidates_signal_calls.csv"), index=False)
rows = []
for x in list(range(1, 444)) + list(range(27329, 27960)):
    r = {"position": x, "RNA_log10_mean_depth": round(float(LOGRNA[x]), 4)}
    for al in SR.RIBO_VERO_ALIASES:
        r[al + "_CPM"] = round(float(CPM[al][x]), 4)
    rows.append(r)
pd.DataFrame(rows).to_csv(os.path.join(OUT, "Fig5AC_track_source_data.csv"), index=False)
json.dump({"rule_5prime": RULE_5P, "rule_3prime": RULE_3P, "rna_libraries": list(SR.RNA_VERO_ALIASES), "ribo_libraries": list(SR.RIBO_VERO_ALIASES),
           "max_sqrtCPM_5prime_region": {al: round(float(np.sqrt(CPM[al][1:444].max())), 2) for al in SR.RIBO_VERO_ALIASES}},
          open(os.path.join(OUT, "Fig5_independent_run_info.json"), "w"), indent=1)
print("red 5' ORFs:", E5[E5.signal_call].candidate_id.tolist())
print("red 3' ORFs:", C3[C3.signal_call].candidate_id.tolist())
