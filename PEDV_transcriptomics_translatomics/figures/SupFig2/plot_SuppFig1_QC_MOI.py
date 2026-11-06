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
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from fig1_moi_data import SRC, OUT, load_metadata

SF = f"{SRC}/figdata/suppfig1"
OD = f"{OUT}/figures/SuppFig1_QC_MOI"
os.makedirs(OD, exist_ok=True)
MM = 1 / 25.4
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "svg.fonttype": "none", "font.size": 7.5,
                     "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8})

TIMES = [3, 5, 8, 12, 16]
TCOL = {3: "#3A2448", 5: "#3A4E89", 8: "#288D85", 12: "#62BB5D", 16: "#F0DE3E"}
TRCOL = {"CHX": "#A5C3E7", "HAR": "#F5BED6"}
MOICOL = {0.1: "#D9D9D9", 1.0: "#6E6E6E"}
CNT = {"total_raw_reads": ("#E5E5E7", "Total raw"), "clean_reads": ("#B9B9B9", "Clean"),
       "uniquely_mapped_reads": ("#666465", "Uniquely mapped")}
OVL_C, UNQ_C, PEDV_C, HOST_C = "#56B4E9", "#000000", "#0072B2", "#009E73"
CM = LinearSegmentedColormap.from_list("corr", ["#F7F7F7", "#FDDBC7", "#D6604D", "#B2182B"])
MOI_LAB = {0.1: "0.1", 1.0: "1"}

meta = load_metadata()[["sample_id", "MOI", "treatment", "hpi"]]
st = pd.read_csv(f"{SF}/Supplementary_Fig1C_RNAseq_mapping_statistics.tsv", sep="\t").merge(meta, on="sample_id")
assert len(st) == 20
st["trk"] = st.treatment.map({"CHX": 0, "HAR": 1})
ORDER = st.sort_values(["MOI", "trk", "hpi"]).sample_id.tolist()
st = st.set_index("sample_id").loc[ORDER].reset_index()
st["label"] = st.treatment + " " + st.hpi.astype(str) + "h"
n = len(st)
st["y"] = np.arange(n)[::-1]


def R(fig, x, y, w, h):
    W, H = fig.get_size_inches() / MM
    return [x / W, y / H, w / W, h / H]


def save(fig, name):
    for ext, kw in (("png", {"dpi": 300}), ("pdf", {}), ("svg", {})):
        fig.savefig(f"{OD}/{name}.{ext}", facecolor="white", **kw)
    plt.close(fig)


# ---------------------------------- A + B ----------------------------------
def drawAB(fig, x, y, w, h, letters=True):
    lab_w = 19.5
    avail = w - lab_w
    wr = np.array([2.0, 0.85, 0.95]); gap = 4.0
    ws = (avail - 2 * gap) * wr / wr.sum()
    xs = [x + lab_w, x + lab_w + ws[0] + gap, x + lab_w + ws[0] + ws[1] + 2 * gap]
    plot_h = h - 12.0
    axs = [fig.add_axes(R(fig, xs[i], y + 11.0, ws[i], plot_h - 11.0 + 2)) for i in range(3)]
    a1, a2, a3 = axs
    hh = 0.26
    for i, k in enumerate(CNT):
        a1.barh(st.y + (1 - i) * hh, st[k] / 1e6, height=hh, color=CNT[k][0], zorder=3)
    a1.set_xlim(0, st.total_raw_reads.max() / 1e6 * 1.02)
    a1.set_xlabel("Read pairs (million)", fontsize=8)
    a2.scatter(st.overall_mapping_rate_pct, st.y, s=18, marker="o", facecolor=OVL_C, edgecolor="none", zorder=3)
    a2.scatter(st.uniq_mapping_rate_pct, st.y, s=18, marker="D", facecolor=UNQ_C, edgecolor="none", zorder=3)
    a2.set_xlim(88, 99); a2.set_xticks([90, 95])
    a2.set_xlabel("Mapping rate\n(% of clean)", fontsize=8)
    a3.scatter(st.PEDV_mapped_reads / 1e6, st.y, s=18, marker="o", facecolor=PEDV_C, edgecolor="none", zorder=3)
    a3.scatter(st.host_mapped_reads / 1e6, st.y, s=18, marker="s", facecolor=HOST_C, edgecolor="none", zorder=3)
    a3.set_xscale("log"); a3.set_xlim(0.08, 40)
    a3.set_xticks([0.1, 1, 10]); a3.set_xticklabels(["0.1", "1", "10"])
    a3.set_xlabel("Uniquely mapped\nreads (million)", fontsize=8)
    for ax in axs:
        ax.set_ylim(-0.6, n - 0.4); ax.set_yticks(st.y)
        ax.grid(axis="x", color="#EEEEEE", lw=0.5, zorder=0)
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
        ax.tick_params(labelsize=7, length=3)

        ax.axhline(n / 2 - 0.5, color="#555555", lw=0.9, zorder=1)
        for blk in (0, 1):
            ax.axhline(n - 1 - (blk * 10 + 4.5), color="#BBBBBB", lw=0.7, ls="--", zorder=1)
    a1.set_yticklabels(st.label, fontsize=7)
    a2.set_yticklabels([]); a3.set_yticklabels([])

    for moi, (ya, yb) in {0.1: (n - 1, n / 2), 1.0: (n / 2 - 1, 0)}.items():
        yc = (ya + yb) / 2
        Wmm = fig.get_size_inches()[0] / MM
        tr = matplotlib.transforms.blended_transform_factory(fig.transFigure, a1.transData)
        fig.text((x + 2.2) / Wmm, yc, f"MOI = {MOI_LAB[moi]}", rotation=90, ha="center", va="center", fontsize=7.5,
                 fontweight="bold", transform=tr)

    a1.legend(handles=[Patch(color=CNT[k][0], label=CNT[k][1]) for k in CNT], loc="lower left",
              bbox_to_anchor=(-0.02, 1.01), ncol=3, frameon=False, fontsize=7.5, handlelength=1, handleheight=1,
              columnspacing=0.8, handletextpad=0.4)
    a2.legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor=OVL_C, markersize=5, label="Overall"),
                       Line2D([0], [0], marker="D", color="none", markerfacecolor=UNQ_C, markersize=5, label="Unique")],
              loc="lower left", bbox_to_anchor=(-0.05, 1.01), ncol=2, frameon=False, fontsize=7.5,
              columnspacing=0.6, handletextpad=0.1)
    a3.legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor=PEDV_C, markersize=5, label="PEDV"),
                       Line2D([0], [0], marker="s", color="none", markerfacecolor=HOST_C, markersize=5, label="Host")],
              loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=2, frameon=False, fontsize=7.5,
              columnspacing=0.6, handletextpad=0.1)
    return axs


# ---------------------------------- C ----------------------------------
P = pd.read_csv(f"{SF}/Supplementary_Fig1D_RNAseq_Pearson_correlation.tsv", sep="\t", index_col=0)
PM = P.loc[ORDER, ORDER].values
VMIN = np.floor(PM[~np.eye(n, dtype=bool)].min() * 100) / 100


def strip(fig, rect, colors, horizontal):
    ax = fig.add_axes(rect)
    arr = np.array([matplotlib.colors.to_rgb(c) for c in colors])[None, :, :] if horizontal else\
        np.array([matplotlib.colors.to_rgb(c) for c in colors])[:, None, :]
    ax.imshow(arr, aspect="auto", interpolation="nearest")
    ax.axis("off")
    return ax


def drawC(fig, x, y, w, h):
    sw, g = 2.0, 0.5
    side = 3 * sw + 3 * g
    top_leg, bot = 11.0, 6 * 1.0 + 13.0
    s = min(w - side - 0.5, h - top_leg - side - 17.0 - 1.0)
    hx, hy = x + side + (w - side - s) / 2 + 2.0, y + 17.0 + side
    ax = fig.add_axes(R(fig, hx, hy, s, s))
    ax.imshow(PM, cmap=CM, vmin=VMIN, vmax=1, aspect="auto", interpolation="nearest")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_linewidth(0.8)
    ax.axhline(n / 2 - 0.5, color="#333333", lw=1.1); ax.axvline(n / 2 - 0.5, color="#333333", lw=1.1)
    for blk in (0, 1):
        p = blk * 10 + 4.5
        ax.axhline(p, color="#777777", lw=0.5); ax.axvline(p, color="#777777", lw=0.5)
    lab = {"trt": st.treatment.map(TRCOL).tolist(), "hpi": st.hpi.map(TCOL).tolist(), "moi": st.MOI.map(MOICOL).tolist()}

    strip(fig, R(fig, hx - (sw + g), hy, sw, s), lab["trt"], False)
    strip(fig, R(fig, hx - 2 * (sw + g), hy, sw, s), lab["hpi"], False)
    strip(fig, R(fig, hx - 3 * (sw + g), hy, sw, s), lab["moi"], False)

    strip(fig, R(fig, hx, hy - (sw + g), s, sw), lab["trt"], True)
    strip(fig, R(fig, hx, hy - 2 * (sw + g), s, sw), lab["hpi"], True)
    strip(fig, R(fig, hx, hy - 3 * (sw + g), s, sw), lab["moi"], True)

    ax.legend(handles=[Patch(facecolor=TCOL[t], label=f"{t} h") for t in TIMES], loc="lower left",
              bbox_to_anchor=(-0.02, 1.07), ncol=5, frameon=False, fontsize=7.5, handlelength=1.1, handleheight=0.9,
              columnspacing=0.8, handletextpad=0.3)
    ax.add_artist(ax.get_legend())
    ax.legend(handles=[Patch(facecolor=TRCOL["CHX"], label="CHX"), Patch(facecolor=TRCOL["HAR"], label="HAR"),
                       Patch(facecolor=MOICOL[0.1], edgecolor="#999999", linewidth=0.4, label="MOI 0.1"),
                       Patch(facecolor=MOICOL[1.0], label="MOI 1")], loc="lower left", bbox_to_anchor=(-0.02, 1.0),
              ncol=4, frameon=False, fontsize=7.5, handlelength=1.1, handleheight=0.9, columnspacing=0.8, handletextpad=0.3)

    cax = fig.add_axes(R(fig, hx, y + 10.0, s, 2.5))
    grad = np.linspace(VMIN, 1, 256)[None, :]
    cax.imshow(grad, cmap=CM, vmin=VMIN, vmax=1, aspect="auto", extent=[VMIN, 1, 0, 1])
    cax.set_yticks([])
    cax.set_xticks(np.round(np.arange(VMIN, 1.0001, 0.01), 2))
    cax.set_xticklabels([f"{v:.2f}" for v in np.round(np.arange(VMIN, 1.0001, 0.01), 2)], fontsize=7)
    cax.set_xlabel("Pearson r (host log2(TPM+1))", fontsize=8, labelpad=2)
    cax.tick_params(length=2)
    return ax


# ---------------------------------- D ----------------------------------
PC = pd.read_csv(f"{SF}/Supplementary_Fig1E_RNAseq_PCA_coordinates.tsv", sep="\t").merge(meta, on="sample_id")
VAR = pd.read_csv(f"{SF}/Supplementary_Fig1E_RNAseq_PCA_variance.tsv", sep="\t").variance_pct.values
SHAPE = {"CHX": "o", "HAR": "s"}


def cov_ellipse(ax, xy, color, level=0.95):
    from scipy.stats import chi2
    from matplotlib.patches import Ellipse
    c = np.cov(xy.T)
    vals, vecs = np.linalg.eigh(c)
    k = np.sqrt(chi2.ppf(level, 2))
    ang = np.degrees(np.arctan2(vecs[1, 1], vecs[0, 1]))
    kw = dict(width=2 * k * np.sqrt(vals[1]), height=2 * k * np.sqrt(vals[0]), angle=ang)
    ax.add_patch(Ellipse(xy.mean(0), facecolor=color, alpha=0.13, edgecolor="none", zorder=1, **kw))
    ax.add_patch(Ellipse(xy.mean(0), facecolor="none", edgecolor=color, lw=1.0, alpha=0.9, zorder=2, **kw))


def drawD(fig, x, y, w, h, ellipse=False):
    ax = fig.add_axes(R(fig, x + 14, y + 10, w - 16, h - 11))
    if ellipse:
        for t in TIMES:
            g = PC[PC.hpi == t]
            assert len(g) == 4
            cov_ellipse(ax, g[["PC1", "PC2"]].values, TCOL[t])
    for r in PC.itertuples():
        c = TCOL[r.hpi]
        if r.MOI == 1.0:
            ax.scatter(r.PC1, r.PC2, s=46, marker=SHAPE[r.treatment], facecolor=c, edgecolor="#333333", linewidths=0.6, zorder=3)
        else:
            ax.scatter(r.PC1, r.PC2, s=46, marker=SHAPE[r.treatment], facecolor="white", edgecolor=c, linewidths=1.7, zorder=3)
    ax.set_xlabel(f"PC1 ({VAR[0]:.1f}% variance)", fontsize=8.5)
    ax.set_ylabel(f"PC2 ({VAR[1]:.1f}% variance)", fontsize=8.5)
    ax.axhline(0, color="#EEEEEE", lw=0.6, zorder=0); ax.axvline(0, color="#EEEEEE", lw=0.6, zorder=0)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.tick_params(labelsize=7.5)
    ax.set_ylim(-32, 56) if ellipse else ax.set_ylim(-25, 42)
    hp = [Line2D([0], [0], marker="o", color="none", markerfacecolor=TCOL[t], markeredgecolor="#333", markersize=6.5,
                 label=f"{t} h") for t in TIMES]
    ex = [Line2D([0], [0], marker="o", color="none", markerfacecolor="#BBBBBB", markeredgecolor="#333", markersize=6.5, label="CHX"),
          Line2D([0], [0], marker="s", color="none", markerfacecolor="#BBBBBB", markeredgecolor="#333", markersize=6.5, label="HAR"),
          Line2D([0], [0], marker="o", color="none", markerfacecolor="#777777", markeredgecolor="#333", markersize=6.5, label="MOI 1"),
          Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor="#777777", markeredgewidth=1.6,
                 markersize=6.5, label="MOI 0.1"),
          Line2D([0], [0], color="none", label="")]
    hs = [v for pair in zip(hp, ex) for v in pair]
    lg = ax.legend(handles=hs, ncol=5, loc="upper right", bbox_to_anchor=(1.0, 1.0), fontsize=7, handletextpad=0.1,
                   columnspacing=0.45, labelspacing=0.5, borderpad=0.4, frameon=True, fancybox=False, edgecolor="#333333")
    lg.get_frame().set_linestyle("--"); lg.get_frame().set_linewidth(0.7)
    return ax



def main():
    f = plt.figure(figsize=(180 * MM, 98 * MM)); drawAB(f, 0, 0, 180, 98); save(f, "SuppFig1_AB_library_and_mapping_QC_MOI")
    f = plt.figure(figsize=(105 * MM, 118 * MM)); drawC(f, 0, 0, 105, 118); save(f, "SuppFig1_C_sample_correlation_MOI")
    f = plt.figure(figsize=(120 * MM, 85 * MM)); drawD(f, 0, 0, 120, 85); save(f, "SuppFig1_D_PCA_MOI")

    f = plt.figure(figsize=(180 * MM, 205 * MM))
    drawAB(f, 0, 107, 180, 98)
    drawC(f, 0, 0, 92, 105)
    drawD(f, 88, 8, 92, 80)
    for L, (xx, yy) in {"A": (0.5, 203.5), "B": (92, 203.5), "C": (0.5, 103), "D": (92, 91)}.items():
        W, H = f.get_size_inches() / MM
        f.text(xx / W, yy / H, L, fontsize=13, fontweight="bold", va="top", ha="left")
    save(f, "SuppFig1_composite_ABCD_MOI_preview")
    f = plt.figure(figsize=(120 * MM, 85 * MM)); drawD(f, 0, 0, 120, 85, ellipse=True); save(f, "SuppFig1_D_PCA_MOI_ellipse")
    PC[["sample_id", "MOI", "treatment", "hpi", "PC1", "PC2"]].sort_values(["hpi", "MOI", "treatment"]).to_csv(
        f"{OD}\SuppFig1_D_PCA_ellipse_groups.csv", index=False)

    cols = ["sample_id", "MOI", "treatment", "hpi", "total_raw_reads", "clean_reads", "uniquely_mapped_reads",
            "multi_mapped_reads", "unmapped_reads", "uniq_mapping_rate_pct", "overall_mapping_rate_pct",
            "PEDV_mapped_reads", "host_mapped_reads"]
    st[cols].to_csv(f"{OD}/SuppFig1_AB_source_data.csv", index=False)
    P.loc[ORDER, ORDER].rename(index=lambda s: s, columns=lambda s: s).to_csv(f"{OD}/SuppFig1_C_pearson_matrix_ordered.csv")
    PC[["sample_id", "MOI", "treatment", "hpi", "PC1", "PC2", "PC3", "PC4"]].sort_values(["MOI", "treatment", "hpi"]).to_csv(
        f"{OD}/SuppFig1_D_PCA_source_data.csv", index=False)
    meta.set_index("sample_id").loc[ORDER].reset_index().to_csv(f"{OD}/SuppFig1_row_order.csv", index=False)
    print("vmin", VMIN, "PC var", VAR[:2].round(1))


if __name__ == "__main__":
    main()
