#!/usr/bin/env python3
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
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from fig1_moi_data import SRC, OUT, TIMES, GENOME_LEN, REG, per_sample_table, smooth

FIG = f"{OUT}/figures"
os.makedirs(FIG, exist_ok=True)
MM = 1 / 25.4
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "svg.fonttype": "none", "font.size": 8,
                     "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8})


PEDV_C, HOST_C, OTHER_C = "#0072B2", "#009E73", "#BDBDBD"
CHX_C, HAR_C = "#000000", "#ED8F98"
CHX_TIME = {3: "#AF9CE0", 5: "#E84271", 8: "#000000", 12: "#236EB0", 16: "#608713"}
HAR_TIME = {3: "#FEC29B", 5: "#4BFC84", 8: "#AFAFAF", 12: "#D0801B", 16: "#9D320E"}
BAND = {"ORF1a": "#E1EFFA", "ORF1b": "#E5F2E1", "S": "#F9DEE5", "ORF3": "#FBE3C7", "E": "#E6DAE8",
        "M": "#E4F1DF", "N": "#D5E7F3", "3′UTR": "#ECEBE9"}
GENES = {"ORF1a": (1, 12591), "ORF1b": (12592, 20609), "S": (20610, 24757), "ORF3": (24758, 25033),
         "E": (25361, 25600), "M": (25608, 26288), "N": (26300, 27625), "3′UTR": (27626, 27959)}
CM = LinearSegmentedColormap.from_list("d", ["#2166AC", "#F7F7F7", "#B2182B"])
MOI_LAB = {0.1: "0.1", 1.0: "1"}
YLAB_C = "Relative coverage\n(per 10⁶ summed viral depth)"

meta, S = per_sample_table()
S = S.sort_values(["MOI", "hpi", "treatment"]).reset_index(drop=True)
STYLE_ROWS = []


def save(fig, name, w_mm, h_mm, note):
    for ext, kw in (("png", {"dpi": 300}), ("pdf", {}), ("svg", {})):
        fig.savefig(f"{FIG}/{name}.{ext}", facecolor="white", bbox_inches="tight", **kw)
    plt.close(fig)
    STYLE_ROWS.append((name, w_mm, h_mm, note))


def moi_tag(moi):
    return f"MOI = {MOI_LAB[moi]}"


def clean(ax):
    for s_ in ["top", "right"]:
        ax.spines[s_].set_visible(False)
    ax.tick_params(labelsize=8, length=3)


# =============================== Panel B ===============================
def panelB_bars(ax, groups, xpos):
    b = np.zeros(len(groups))
    for k, c, lab in (("pedv", PEDV_C, "PEDV"), ("host", HOST_C, "Host (Vero)"), ("other", OTHER_C, "Other")):
        v = np.array([g[k] for g in groups])
        ax.bar(xpos, v, bottom=b, color=c, width=0.72 if len(xpos) <= 5 else 0.8, edgecolor="white", linewidth=0.4)
        b += v
    ax.set_ylim(0, 1)
    ax.set_yticks([0, .25, .5, .75, 1])
    ax.set_yticklabels(["0.00", "0.25", "0.50", "0.75", "1.00"])
    clean(ax)
    ax.set_ylabel("Relative fraction of reads", fontsize=10)


def legendB(ax):
    ax.legend(handles=[Patch(color=c, label=l) for c, l in ((PEDV_C, "PEDV"), (HOST_C, "Host (Vero)"), (OTHER_C, "Other"))],
              loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, frameon=False, fontsize=8, handlelength=1,
              columnspacing=1.2)


def fracs(d):
    return {"pedv": d.B_pedv_frac.mean(), "host": d.B_host_frac.mean(), "other": d.B_other_frac.mean()}


def figB():
    cols = ["sample_id", "MOI", "treatment", "hpi", "B_total_QC_read_pairs", "B_pedv_pairs", "B_host_pairs", "B_other_pairs",
            "B_pedv_frac", "B_host_frac", "B_other_frac"]
    # main: 1 MOI, mean of CHX+HAR
    d1 = S[S.MOI == 1.0]
    g = [fracs(d1[d1.hpi == t]) for t in TIMES]
    fig, ax = plt.subplots(figsize=(88 * MM, 72 * MM))
    panelB_bars(ax, g, np.arange(5))
    ax.set_xticks(range(5)); ax.set_xticklabels(TIMES)
    ax.set_xlabel("Hours post infection (hpi)", fontsize=10)
    ax.set_title("Relative fraction of sequencing reads", fontsize=10, fontweight="bold", pad=24)
    legendB(ax)
    ax.text(0.0, 1.0, "", transform=ax.transAxes)
    ax.text(0.5, -0.27, "MOI = 1; each bar = mean of CHX and HAR (n = 2); no biological replicates",
            ha="center", va="top", fontsize=7, transform=ax.transAxes)
    save(fig, "Fig1B_main_MOI1", 88, 72, "stacked; y 0-1; bar=mean of CHX+HAR at MOI 1")
    out = d1[cols].copy(); out["plotted_bar_PEDV"] = out.hpi.map({t: g[i]["pedv"] for i, t in enumerate(TIMES)})
    out["plotted_bar_Host"] = out.hpi.map({t: g[i]["host"] for i, t in enumerate(TIMES)})
    out["plotted_bar_Other"] = out.hpi.map({t: g[i]["other"] for i, t in enumerate(TIMES)})
    out.to_csv(f"{FIG}/Fig1B_main_MOI1_source_data.csv", index=False)

    # supplement 1: by MOI, two bars per hpi
    fig, ax = plt.subplots(figsize=(110 * MM, 75 * MM))
    xs, gs, tl = [], [], []
    for i, t in enumerate(TIMES):
        for j, moi in enumerate((0.1, 1.0)):
            xs.append(i * 3 + j); gs.append(fracs(S[(S.hpi == t) & (S.MOI == moi)])); tl.append(MOI_LAB[moi])
    xs = np.array(xs)
    panelB_bars(ax, gs, xs)
    ax.set_xticks(xs); ax.set_xticklabels(tl, fontsize=7)
    for i, t in enumerate(TIMES):
        ax.text(i * 3 + 0.5, -0.17, str(t), ha="center", va="top", fontsize=8, transform=ax.get_xaxis_transform())
    ax.text(-0.08, -0.07, "MOI", ha="right", va="top", fontsize=7, transform=ax.transAxes)
    ax.text(-0.08, -0.17, "hpi", ha="right", va="top", fontsize=8, transform=ax.transAxes)
    ax.set_xlabel("Hours post infection (hpi)", fontsize=10, labelpad=34)
    ax.set_title("Relative fraction of sequencing reads, by MOI", fontsize=10, fontweight="bold", pad=24)
    legendB(ax)
    ax.text(0.5, -0.50, "each bar = mean of CHX and HAR within one MOI (n = 2); no biological replicates",
            ha="center", va="top", fontsize=7, transform=ax.transAxes)
    save(fig, "Fig1B_suppl_byMOI", 110, 75, "two bars per hpi (0.1 / 1)")
    S[cols].to_csv(f"{FIG}/Fig1B_suppl_byMOI_source_data.csv", index=False)

    # supplement 2: per-sample bars, two MOI panels
    fig, axes = plt.subplots(1, 2, figsize=(180 * MM, 80 * MM), sharey=True)
    for ax, moi in zip(axes, (0.1, 1.0)):
        d = S[S.MOI == moi].sort_values(["hpi", "treatment"])
        xs, labs, pos = [], [], 0
        for t in TIMES:
            for tr in ("CHX", "HAR"):
                xs.append(pos); labs.append(f"{tr} {t}h"); pos += 1
            pos += 0.6
        gs = [{"pedv": r.B_pedv_frac, "host": r.B_host_frac, "other": r.B_other_frac} for r in d.itertuples()]
        panelB_bars(ax, gs, np.array(xs))
        ax.set_xticks(xs); ax.set_xticklabels(labs, rotation=90, fontsize=7)
        ax.set_title(moi_tag(moi), fontsize=9, fontweight="bold")
        if moi == 1.0:
            ax.set_ylabel("")
    fig.legend(handles=[Patch(color=c, label=l) for c, l in ((PEDV_C, "PEDV"), (HOST_C, "Host (Vero)"), (OTHER_C, "Other"))],
               loc="upper center", ncol=3, frameon=False, fontsize=8, handlelength=1, columnspacing=1.2, bbox_to_anchor=(0.5, 1.04))
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.text(0.5, -0.02, "one bar per sample (n = 1 per condition)", ha="center", va="top", fontsize=7)
    save(fig, "Fig1B_suppl_perSample", 180, 80, "20 bars, one per sample")
    S[cols].to_csv(f"{FIG}/Fig1B_suppl_perSample_source_data.csv", index=False)


# =============================== Panel C ===============================
def band(ax):
    tr = ax.get_xaxis_transform()
    for g, (a, b) in GENES.items():
        ax.fill_between([a, b], 0.935, 1.0, color=BAND[g], transform=tr, lw=0, zorder=1)
        if g == "3′UTR":
            ax.text(b + 150, 0.967, g, transform=tr, ha="left", va="center", fontsize=7)
        elif g in ("S", "N", "ORF1a", "ORF1b", "E", "M"):
            ax.text((a + b) / 2, 0.967, g, transform=tr, ha="center", va="center", fontsize=8 if g in ("ORF1a", "ORF1b", "S", "N") else 7)
        if g in ("S", "ORF3", "E", "M", "N"):
            ax.fill_between([a, b], 0, 0.935, color=BAND[g], alpha=0.45, transform=tr, lw=0, zorder=0)


def covC(ax, moi, ver):
    x = np.arange(1, GENOME_LEN + 1)
    d = S[S.MOI == moi]
    rows, ymin, ymax = [], np.inf, 0
    for r in d.itertuples():
        dep = pd.read_csv(f"{SRC}/figdata/fig1_depth/{r.sample_id}.{ver}.depth.tsv", sep="\t", header=None)[2].values.astype(float)
        rel = smooth(dep / dep.sum() * 1e6)
        tot = {"raw": r.D_total_mapped_fragments, "dedup": r.D_total_mapped_fragments}[ver]
        absd = dep / (r.D_total_mapped_fragments / 1e6)
        col = (CHX_TIME if r.treatment == "CHX" else HAR_TIME)[r.hpi]
        ax.plot(x, np.where(rel > 0, rel, np.nan), color=col, lw=1.0, ls="-" if r.treatment == "CHX" else "--", zorder=3)
        ymin, ymax = min(ymin, rel[rel > 0].min()), max(ymax, rel.max())
        for p in range(0, GENOME_LEN, 5):
            rows.append((r.sample_id, moi, r.treatment, r.hpi, ver, p + 1, dep[p], rel[p], absd[p]))
    ax.set_yscale("log")
    lo = 10 ** np.floor(np.log10(ymin))
    ax.set_ylim(lo, ymax * 10 ** 0.6)
    ax.set_xlim(-1500, GENOME_LEN + 3200)
    ax.set_xticks([0, 5000, 10000, 15000, 20000, 25000])
    ax.set_xlabel("genome position (nt, PEDV_GIb)", fontsize=10)
    ax.set_ylabel(YLAB_C, fontsize=9)
    ax.tick_params(labelsize=8)
    band(ax)
    return pd.DataFrame(rows, columns=["sample_id", "MOI", "treatment", "hpi", "bam", "position", "depth",
                                       "relative_coverage_per_1e6_summed_viral_depth_101nt_smoothed",
                                       "depth_per_M_unique_fragments_raw_MAPQ255_denominator"])


def legendC(fig, ax):
    hs, ls = [], []
    hs.append(Line2D([], [], alpha=0)); ls.append("CHX")
    for t in TIMES:
        hs.append(Line2D([], [], color=CHX_TIME[t], lw=1.2)); ls.append(f"{t}h")
    hs.append(Line2D([], [], alpha=0)); ls.append("|")
    hs.append(Line2D([], [], alpha=0)); ls.append("HAR")
    for t in TIMES:
        hs.append(Line2D([], [], color=HAR_TIME[t], lw=1.2, ls="--")); ls.append(f"{t}h")
    ax.legend(hs, ls, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=13, frameon=False, fontsize=8,
              handlelength=1.6, handletextpad=0.3, columnspacing=0.7)


def figC():
    jobs = [("Fig1C_main_MOI1_raw", 1.0, "raw", "Genome-wide RNA-seq reads coverage"),
            ("Fig1C_suppl_MOI0.1_raw", 0.1, "raw", "Genome-wide RNA-seq reads coverage")]
    for name, moi, ver, title in jobs:
        fig, ax = plt.subplots(figsize=(180 * MM, 78 * MM))
        df = covC(ax, moi, ver)
        ax.set_title(title, fontsize=11, fontweight="bold", pad=26)
        legendC(fig, ax)
        ax.text(0.012, 0.86, moi_tag(moi), transform=ax.transAxes, fontsize=9, fontweight="bold")
        fig.tight_layout()
        save(fig, name, 180, 78, f"{moi_tag(moi)}; {ver} BAM, MAPQ=255; log y; 10 lines = 10 single samples")
        df.to_csv(f"{FIG}/{name}_source_data.csv", index=False)

    fig, axes = plt.subplots(2, 1, figsize=(180 * MM, 140 * MM), sharex=True)
    dfs = []
    for ax, moi in zip(axes, (1.0, 0.1)):
        dfs.append(covC(ax, moi, "dedup"))
        ax.text(0.012, 0.86, f"{moi_tag(moi)} (duplicate-marked BAM, sensitivity)", transform=ax.transAxes, fontsize=8, fontweight="bold")
        if ax is not axes[-1]:
            ax.set_xlabel("")
    axes[0].set_title("Genome-wide RNA-seq reads coverage, duplicates removed", fontsize=11, fontweight="bold", pad=26)
    legendC(fig, axes[0])
    fig.tight_layout()
    save(fig, "Fig1C_suppl_dedup_MOI1_MOI0.1", 180, 140, "dedup BAM sensitivity; each panel one MOI")
    pd.concat(dfs).to_csv(f"{FIG}/Fig1C_suppl_dedup_MOI1_MOI0.1_source_data.csv", index=False)


# =============================== Panel D ===============================
def kin(ax, series):
    """series: list of (label, df, color, marker, filled, ls, xoff)"""
    for lab, d, c, m, filled, ls, off in series:
        d = d.sort_values("hpi")
        ax.plot(d.hpi + off, d.D_pedv_pct_of_mapped, color=c, lw=1.4, ls=ls, marker=m, ms=5,
                mfc=c if filled else "white", mec=c, mew=1.2, label=lab, zorder=3)
    ax.set_yscale("log")
    ax.set_ylim(0.3, 100)
    yt = [0.3, 1, 3, 10, 30, 100]
    ax.set_yticks(yt); ax.set_yticklabels([f"{v:g}" for v in yt]); ax.minorticks_off()
    ax.set_xticks(TIMES); ax.set_xlim(2, 17)
    ax.set_xlabel("Hours post infection (hpi)", fontsize=10)
    ax.set_ylabel("PEDV (% of mapped reads)", fontsize=10)
    clean(ax)


def figD():
    cols = ["sample_id", "MOI", "treatment", "hpi", "D_total_mapped_fragments", "D_pedv_mapped_fragments",
            "D_pedv_pct_of_mapped", "B_total_QC_read_pairs"]
    d1 = S[S.MOI == 1.0]
    fig, ax = plt.subplots(figsize=(88 * MM, 72 * MM))
    kin(ax, [("CHX", d1[d1.treatment == "CHX"], CHX_C, "o", True, "-", -0.2),
             ("HAR", d1[d1.treatment == "HAR"], HAR_C, "s", True, "-", 0.2)])
    ax.legend(frameon=False, fontsize=8, loc="upper left", handlelength=1.4)
    ax.text(0.98, 0.04, "MOI = 1", transform=ax.transAxes, ha="right", fontsize=8, fontweight="bold")
    ax.set_title("Total viral RNA kinetics", fontsize=10, fontweight="bold", pad=8)
    fig.tight_layout()
    save(fig, "Fig1D_main_MOI1", 88, 72, "log y 0.3-100; one sample per point; no error bars")
    d1[cols].to_csv(f"{FIG}/Fig1D_main_MOI1_source_data.csv", index=False)

    fig, ax = plt.subplots(figsize=(95 * MM, 75 * MM))
    ser = []
    for moi, ls, filled, moff in ((1.0, "-", True, 0.12), (0.1, "--", False, -0.12)):
        for tr, c, m, toff in (("CHX", CHX_C, "o", -0.2), ("HAR", HAR_C, "s", 0.2)):
            ser.append((f"{tr} · MOI {MOI_LAB[moi]}", S[(S.MOI == moi) & (S.treatment == tr)], c, m, filled, ls, toff + moff))
    kin(ax, ser)
    ax.legend(frameon=False, fontsize=7.5, loc="lower right", handlelength=2.0, ncol=1)
    ax.set_title("Total viral RNA kinetics, by MOI", fontsize=10, fontweight="bold", pad=8)
    fig.tight_layout()
    save(fig, "Fig1D_suppl_bothMOI", 95, 75, "4 series; filled solid = MOI 1, open dashed = MOI 0.1; no error bars")
    S[cols].to_csv(f"{FIG}/Fig1D_suppl_bothMOI_source_data.csv", index=False)


# =============================== Panel F ===============================
REGS = list(REG)


def Fmatrix(moi):
    d = S[S.MOI == moi]
    m = np.array([[d[d.hpi == t][f"F_{g}_cov_per_M_mapped"].mean() for t in TIMES] for g in REGS])
    L = np.log2(m)
    return m, L - L.mean(axis=1, keepdims=True)


def heat(ax, Lc, show_y=True):
    im = ax.imshow(Lc, cmap=CM, norm=TwoSlopeNorm(0, -4, 4), aspect="auto")
    ax.set_xticks(range(5)); ax.set_xticklabels(TIMES, fontsize=8)
    ax.set_yticks(range(len(REGS)))
    ax.set_yticklabels(REGS if show_y else [], fontsize=8)
    ax.set_xlabel("Hours post infection (hpi)", fontsize=10)
    ax.tick_params(length=2)
    return im


def figF():
    rows = []
    for moi in (0.1, 1.0):
        m, Lc = Fmatrix(moi)
        for i, g in enumerate(REGS):
            for j, t in enumerate(TIMES):
                a, b, ex = REG[g]
                rows.append({"MOI": moi, "region": g, "start": a, "end": b,
                             "excluded": ";".join(f"{x}-{y}" for x, y in ex), "hpi": t,
                             "mean_cov_per_M_mapped(CHX,HAR mean)": m[i, j], "log2_rowcentered_within_MOI": Lc[i, j]})
    FD = pd.DataFrame(rows)
    per = S[["sample_id", "MOI", "treatment", "hpi"] + [f"F_{g}_cov_per_M_mapped" for g in REGS]]

    def one(moi, name, w, h):
        m, Lc = Fmatrix(moi)
        fig, ax = plt.subplots(figsize=(w * MM, h * MM))
        im = heat(ax, Lc)
        ax.set_title("Regional RNA accumulation", fontsize=10, fontweight="bold", pad=8)
        cb = fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.32, fraction=0.06, ticks=[-4, -2, 0, 2, 4])
        cb.set_label("Relative coverage (log2, row-centered)", fontsize=8); cb.ax.tick_params(labelsize=8)
        fig.tight_layout()
        fig.text(0.5, -0.02, f"{moi_tag(moi)}; each cell = mean of CHX and HAR (n = 2); centered within MOI", ha="center", va="top", fontsize=7)
        save(fig, name, w, h, f"{moi_tag(moi)}; diverging, zero-centred, range -4..4, no clipping needed")
        FD[FD.MOI == moi].to_csv(f"{FIG}/{name}_source_data.csv", index=False)
        per[per.MOI == moi].to_csv(f"{FIG}/{name}_per_sample_values.csv", index=False)

    one(1.0, "Fig1F_main_MOI1", 88, 78)
    one(0.1, "Fig1F_suppl_MOI0.1", 88, 78)
    fig, axes = plt.subplots(1, 2, figsize=(165 * MM, 80 * MM), gridspec_kw={"wspace": 0.75})
    for k, (ax, moi) in enumerate(zip(axes, (0.1, 1.0))):
        im = heat(ax, Fmatrix(moi)[1], show_y=True)
        ax.set_title(moi_tag(moi), fontsize=9, fontweight="bold")
        if k == 1:
            ax.set_xlabel("Hours post infection (hpi)", fontsize=10)
    fig.suptitle("Regional RNA accumulation", fontsize=10, fontweight="bold", y=1.0)
    cb = fig.colorbar(im, ax=axes, orientation="horizontal", pad=0.28, fraction=0.05, ticks=[-4, -2, 0, 2, 4], shrink=0.5)
    cb.set_label("Relative coverage (log2, row-centered within MOI)", fontsize=8); cb.ax.tick_params(labelsize=8)
    save(fig, "Fig1F_suppl_bothMOI_sideBySide", 150, 80, "two heatmaps, each centred within its own MOI, shared colour scale")
    FD.to_csv(f"{FIG}/Fig1F_suppl_bothMOI_sideBySide_source_data.csv", index=False)
    per.to_csv(f"{FIG}/Fig1F_suppl_bothMOI_sideBySide_per_sample_values.csv", index=False)


if __name__ == "__main__":
    S.to_csv(f"{FIG}/Fig1_per_sample_all_panels_source_data.csv", index=False)
    figB(); figC(); figD(); figF()
    with open(f"{FIG}/Fig1_MOI_style_compliance_note.md", "w", encoding="utf-8") as f:
        f.write("# Fig.1 B/C/D/F (MOI-regrouped) style compliance note\n\n")
        f.write("- Font Arial (≥7 pt), 300 dpi PNG + PDF (fonttype 42) + SVG (editable text); white background; sizes ≤180 mm wide.\n")
        f.write("- Colours: PEDV #0072B2 / Host #009E73 / Other #BDBDBD; CHX #000000, HAR #ED8F98; C uses a fixed 10-colour palette; diverging heatmap #2166AC/#F7F7F7/#B2182B, centred at 0, range -4..4.\n")
        f.write("- No error bars anywhere (n = 1 per MOI × treatment × hpi, no biological replicates); no inferential statistics.\n")
        f.write("- Log axes: D y 0.3–100 (covers all values 0.85–38.9); C y spans all plotted data; no truncated axes; no insets.\n\n")
        f.write("| file | width mm | height mm | note |\n|---|---|---|---|\n")
        for n, w, h, note in STYLE_ROWS:
            f.write(f"| {n} | {w} | {h} | {note} |\n")
    print("done", len(STYLE_ROWS))
