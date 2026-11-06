#!/usr/bin/env python3
"""Step 4: Supplementary Figure (PEDV RNA junction dynamics) A-E, A4 portrait, 180 mm wide.
全部由 data/*.csv(s3/s5 产出)读入, 无手工数值。配色沿用 Fig.2 (canonical 蓝系 / dep / distant / local), 时间梯度沿用 Supp Fig.1 hpi 色。
figure_skill.md: Arial, 最小字号 7 pt, PNG 300 dpi + PDF + SVG, 宽 <= 180 mm, 高尽量 <= 240 mm。"""
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
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import FuncFormatter, NullLocator

D = PEDV_P+"Fig2_Supp_junction_dynamics/data/"
F = PEDV_P+"Fig2_Supp_junction_dynamics/figures/"
MM = 1 / 25.4
FS, FSL = 7, 8                      # 刻度/图例/注释 7 pt; 轴标题 8 pt
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "svg.fonttype": "none", "font.size": FS,
                     "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
                     "xtick.major.size": 2.5, "ytick.major.size": 2.5})
TIMES = [3, 5, 8, 12, 16]
CAN = ["S", "ORF3", "E", "M", "N"]
CANON_COL = {"S": "#9ECAE1", "ORF3": "#6BAED6", "E": "#4292C6", "M": "#2171B5", "N": "#08306B"}   # = Fig.2
CANON_ALL = "#2171B5"
NCC = ["TRSL_dep_noncanonical", "TRSL_indep_distant", "TRSL_indep_local"]
NCOL = {"TRSL_dep_noncanonical": "#B8860B", "TRSL_indep_distant": "#4DAF4A", "TRSL_indep_local": "#984EA3"}  # = Fig.2
NLAB = {"TRSL_dep_noncanonical": "TRS-L-dependent noncanonical", "TRSL_indep_distant": "TRS-L-independent distant",
        "TRSL_indep_local": "TRS-L-independent local"}
TCOL = {3: "#3A2448", 5: "#3A4E89", 8: "#288D85", 12: "#62BB5D", 16: "#F0DE3E"}      # Supp Fig.1 hpi gradient
TR_COL = {"CHX": "#A5C3E7", "HAR": "#F5BED6"}
GENES = {"ORF1a": (293, 12622), "ORF1b": (12592, 20613), "S": (20610, 24758), "ORF3": (24758, 25033), "E": (25361, 25600),
         "M": (25608, 26288), "N": (26300, 27625)}
GENE_FILL = {'ORF1a': '#F3E3EC', 'ORF1b': '#E1EBF5', 'S': '#E6F1E3', 'ORF3': '#FBF0D9', 'E': '#EAE3F3', 'M': '#F8E3E3', 'N': '#E0F0F0'}
GL = 27959
W, H = 180, 250          # mm
MOIS = [0.1, 1.0]
fig = plt.figure(figsize=(W * MM, H * MM))


def ax_mm(x, y_top, w, h, **kw):
    return fig.add_axes([x / W, 1 - (y_top + h) / H, w / W, h / H], **kw)


def txt(x, y_top, s, **kw):
    kw.setdefault("fontsize", FS)
    return fig.text(x / W, 1 - y_top / H, s, **kw)


def label(letter, x, y_top):
    txt(x, y_top, letter, fontsize=12, fontweight="bold", va="top", ha="left")


def clean(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


A = pd.read_csv(D + "panelA_per_sample_class_fractions.csv")
Bd = pd.read_csv(D + "panelB_canonical_sgRNA_per_sample.csv")
pB = pd.read_csv(D + "panelB_canonical_sgRNA_pooled_DESCRIPTIVE.csv")
Cd = pd.read_csv(D + "panelC_noncanonical_per_sample.csv")
pC = pd.read_csv(D + "panelC_noncanonical_pooled_DESCRIPTIVE.csv")
GX = {0.1: 17, 1.0: 102}
GW = 76

# =============================== A ===============================
label("A", 0, 0)
hA_leg = [Patch(color=CANON_ALL, label="TRS-L-dependent canonical")] + [Patch(color=NCOL[n], label=NLAB[n]) for n in NCC]
fig.legend(handles=hA_leg, loc="upper left", bbox_to_anchor=(8 / W, 1 - 0.0 / H), ncol=2, frameon=False, fontsize=FS,
           handlelength=0.9, columnspacing=1.0, handletextpad=0.4, labelspacing=0.15, borderaxespad=0)
layers = [(9.5, "Fraction of all\nclassified reads"), (33.0, "Fraction of\nnoncanonical")]
for layer, (ytop, ylab) in enumerate(layers):
    for moi in MOIS:
        ax = ax_mm(GX[moi], ytop, GW, 13)
        d = A[A.MOI == moi].sort_values(["hpi", "treatment"])
        xpos, labs, pos = [], [], 0
        for t in TIMES:
            for tr in ("CHX", "HAR"):
                xpos.append(pos); labs.append(tr[0]); pos += 1
            pos += 0.7
        xpos = np.array(xpos)
        bottoms = np.zeros(len(d))
        comps = ([("f_canonical", CANON_ALL)] + [("f_" + n, NCOL[n]) for n in NCC]) if layer == 0 else \
            [("nc_frac_" + n, NCOL[n]) for n in NCC]
        for k, c in comps:
            v = d[k].values
            ax.bar(xpos, v, bottom=bottoms, width=0.9, color=c, edgecolor="white", linewidth=0.3)
            bottoms += v
        ax.set_ylim(0, 1); ax.set_xlim(-0.8, xpos.max() + 0.8)
        ax.set_yticks([0, 1]); ax.set_yticklabels(["0", "1"], fontsize=FS)
        ax.set_xticks(xpos); ax.set_xticklabels(labs, fontsize=FS)
        ax.tick_params(length=2, pad=1.2)
        for i, t in enumerate(TIMES):
            ax.text((xpos[2 * i] + xpos[2 * i + 1]) / 2, -0.33, str(t), ha="center", va="top", fontsize=FS,
                    transform=ax.get_xaxis_transform())
        if layer == 1:
            ax.text(-0.012, -0.33, "hpi", ha="right", va="top", fontsize=FS, transform=ax.transAxes)
        if layer == 0:
            ax.set_title(f"MOI = {moi:g}", fontsize=FS, fontweight="bold", pad=2)
        if moi == 0.1:
            ax.set_ylabel(ylab, fontsize=FS, labelpad=2)
        else:
            ax.set_yticklabels([])
        clean(ax)

# =============================== B ===============================
yB = 55
label("B", 0, yB)
hB = [Line2D([0], [0], color=CANON_COL[g], lw=1.6, label=g) for g in CAN]
hB += [Line2D([0], [0], marker="o", color="none", markerfacecolor="#999999", markeredgecolor="#222", markersize=4, label="CHX"),
       Line2D([0], [0], marker="s", color="none", markerfacecolor="#999999", markeredgecolor="#222", markersize=4, label="HAR"),
       Line2D([0], [0], color="#888888", lw=1.0, label="CHX+HAR pooled")]
fig.legend(handles=hB, loc="upper left", bbox_to_anchor=(12 / W, 1 - yB / H), ncol=8, frameon=False, fontsize=FS,
           handlelength=1.2, columnspacing=0.8, handletextpad=0.3, borderaxespad=0)
for moi in MOIS:
    ax = ax_mm(GX[moi], yB + 9, GW, 20)
    d = Bd[Bd.MOI == moi]
    for g in CAN:
        c = CANON_COL[g]
        for tr, m in (("CHX", "o"), ("HAR", "s")):
            dd = d[d.treatment == tr].sort_values("hpi")
            ax.scatter(dd.hpi + (-0.15 if tr == "CHX" else 0.15), dd["frac_" + g], s=12, marker=m, facecolor=c,
                       edgecolor="#222222", linewidths=0.3, zorder=3)
        pp = pB[pB.MOI == moi].sort_values("hpi")
        ax.plot(pp.hpi, pp["frac_" + g], color=c, lw=1.0, zorder=2)
    ax.set_yscale("log"); ax.set_ylim(0.002, 1.3)
    ax.set_yticks([0.01, 0.1, 1]); ax.set_yticklabels(["0.01", "0.1", "1"]); ax.yaxis.set_minor_locator(NullLocator())
    ax.set_xticks(TIMES); ax.set_xlim(2, 17)
    ax.set_xlabel("Hours post infection (hpi)", fontsize=FSL, labelpad=1.5)
    ax.set_title(f"MOI = {moi:g}", fontsize=FS, fontweight="bold", pad=2)
    ax.tick_params(labelsize=FS, length=2.5, pad=1.5)
    if moi == 0.1:
        ax.set_ylabel("Fraction of canonical\nreads", fontsize=FS, labelpad=2)
    else:
        ax.set_yticklabels([])
    clean(ax)

# =============================== C ===============================
yC = 96
label("C", 0, yC)
for ci, (col, ylab, logy) in enumerate([("frac_all_", "Fraction of all\nclassified reads", False),
                                        ("per_M_viral_frag_", "Reads per million\nviral fragments", True)]):
    ax = ax_mm(19 + ci * 42, yC + 6, 25, 42)
    for n in NCC:
        c = NCOL[n]
        for moi in MOIS:
            d = Cd[Cd.MOI == moi]
            filled = moi == 1.0
            for tr, m in (("CHX", "o"), ("HAR", "s")):
                dd = d[d.treatment == tr].sort_values("hpi")
                ax.scatter(dd.hpi + (-0.18 if tr == "CHX" else 0.18), dd[col + n], s=11, marker=m,
                           facecolor=c if filled else "white", edgecolor=c, linewidths=0.8, zorder=3)
            pp = pC[pC.MOI == moi].sort_values("hpi")
            ax.plot(pp.hpi, pp[col + n], color=c, lw=1.0, ls="-" if filled else "--", zorder=2)
    ax.set_xticks(TIMES); ax.set_xlim(2, 17)
    ax.set_xticklabels(["3", "", "8", "", "16"])
    ax.set_xlabel("hpi", fontsize=FSL, labelpad=1.5)
    ax.set_ylabel(ylab, fontsize=FS, labelpad=2)
    ax.tick_params(labelsize=FS, length=2.5, pad=1.5)
    if logy:
        ax.set_yscale("log"); ax.set_ylim(150, 1800)
        ax.set_yticks([200, 500, 1000]); ax.set_yticklabels(["200", "500", "1000"]); ax.yaxis.set_minor_locator(NullLocator())
    else:
        ax.set_ylim(0, 0.045); ax.set_yticks([0, 0.02, 0.04])
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.2f}"))
    clean(ax)
hC = [Line2D([0], [0], color=NCOL[n], lw=1.6, label=NLAB[n]) for n in NCC]
hC += [Line2D([0], [0], color="#444", lw=1.0, ls="-", marker="o", markerfacecolor="#444", markersize=4, label="MOI = 1"),
       Line2D([0], [0], color="#444", lw=1.0, ls="--", marker="o", markerfacecolor="white", markersize=4, label="MOI = 0.1"),
       Line2D([0], [0], marker="o", color="none", markerfacecolor="#888", markeredgecolor="#888", markersize=4, label="CHX"),
       Line2D([0], [0], marker="s", color="none", markerfacecolor="#888", markeredgecolor="#888", markersize=4, label="HAR")]
fig.legend(handles=hC, loc="upper left", bbox_to_anchor=(2 / W, 1 - (yC + 55) / H), ncol=2, frameon=False, fontsize=FS,
           handlelength=1.6, labelspacing=0.2, columnspacing=1.0, borderaxespad=0)
hD = [Patch(color=NCOL[n], label=NLAB[n]) for n in NCC] + [Patch(color=TR_COL["CHX"], label="CHX"),
                                                            Patch(color=TR_COL["HAR"], label="HAR")]

# =============================== D ===============================
label("D", 92, yC)
sel = pd.read_csv(D + "panelD_selected_clusters.csv")
rate = pd.read_csv(D + "panelD_selected_cluster_x_sample_per_M_viral_frag.csv", index_col=0)
inf = pd.read_csv(D + "panelD_selected_informative_mask.csv", index_col=0).astype(bool)
meta = A[["sample_id", "MOI", "treatment", "hpi"]].sort_values(["MOI", "hpi", "treatment"])
cols = meta.sample_id.tolist()
sel["order"] = sel.junction_class.map({n: i for i, n in enumerate(NCC)})
sel = sel.sort_values(["order", "reads_total"], ascending=[True, False]).reset_index(drop=True)
Z = np.log2(rate.loc[sel.cluster_id, cols] + 1.0)
ok = inf.loc[sel.cluster_id, cols]
Zs = Z.copy()
for cid in sel.cluster_id:
    v = Z.loc[cid]; m = ok.loc[cid]
    mu, sd = v[m].mean(), v[m].std(ddof=0)
    Zs.loc[cid] = ((v - mu) / sd) if sd > 0 else 0.0
Zs = Zs.clip(-2.5, 2.5)
CM = LinearSegmentedColormap.from_list("d", ["#2166AC", "#F7F7F7", "#B2182B"])
nr = len(sel)
hx, hy, hw, hh = 124, yC + 7.5, 54, 2.45 * nr
axD = ax_mm(hx, hy, hw, hh)
axD.imshow(np.ma.array(Zs.values, mask=~ok.values), cmap=CM, vmin=-2.5, vmax=2.5, aspect="auto", interpolation="nearest")
for i in range(nr):
    for j in range(len(cols)):
        if not ok.values[i, j]:
            axD.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor="#EDEDED", edgecolor="#A0A0A0", hatch="////", linewidth=0))
axD.set_xticks([]); axD.set_yticks(range(nr))
rl = []
for r in sel.itertuples():
    s_ = f"{r.donor_rep}→{r.acc_rep}"
    if r.donor_rep == 27620 and r.acc_rep == 27682:
        s_ += " †"
    rl.append(s_)
axD.set_yticklabels(rl, fontsize=FS)
axD.tick_params(length=0, pad=7)
axD.axvline(9.5, color="#222", lw=0.9)
for bnd in np.cumsum(sel.groupby("order").size().values)[:-1]:
    axD.axhline(bnd - 0.5, color="#222", lw=0.6)
for sp in axD.spines.values():
    sp.set_linewidth(0.6)
side = ax_mm(hx - 2.4, hy, 1.8, hh)
side.imshow(np.array([matplotlib.colors.to_rgb(NCOL[c]) for c in sel.junction_class])[:, None, :], aspect="auto"); side.axis("off")
axT = ax_mm(hx, hy - 2.8, hw, 2.0)
axT.imshow(np.array([matplotlib.colors.to_rgb(TR_COL[t]) for t in meta.treatment])[None, :, :], aspect="auto"); axT.axis("off")
for gi, moi in enumerate(MOIS):
    txt(hx + hw * (0.25 + 0.5 * gi), hy - 3.2, f"MOI = {moi:g}", ha="center", va="bottom", fontweight="bold")
hp = ax_mm(hx, hy + hh + 1.0, hw, 2.0)
hp.imshow(np.array([matplotlib.colors.to_rgb(TCOL[t]) for t in meta.hpi])[None, :, :], aspect="auto"); hp.axis("off")
for i, t in enumerate(TIMES):
    for g in (0, 1):
        j = g * 10 + 2 * i
        txt(hx + hw * (j + 1) / 20, hy + hh + 3.4, str(t), ha="center", va="top", fontsize=FS)
txt(hx - 1.5, hy + hh + 3.4, "hpi", ha="right", va="top", fontsize=FS)
cb = fig.add_axes([(hx + 8) / W, 1 - (hy + hh + 11.0) / H, 38 / W, 2.2 / H])
cb.imshow(np.linspace(-2.5, 2.5, 256)[None, :], cmap=CM, aspect="auto", extent=[-2.5, 2.5, 0, 1])
cb.set_yticks([]); cb.set_xticks([-2, 0, 2]); cb.tick_params(labelsize=FS, length=2, pad=1)
cb.set_xlabel("Row z-score", fontsize=FS, labelpad=1)
# D 图例放在 C 下方左半版(含 hatched 与 † 说明)
fig.legend(handles=[Patch(color=TR_COL["CHX"], label="CHX (column strip)"), Patch(color=TR_COL["HAR"], label="HAR (column strip)")],
           loc="upper left", bbox_to_anchor=(2 / W, 1 - (yC + 70) / H), ncol=2, frameon=False, fontsize=FS,
           handlelength=1.0, columnspacing=1.0, borderaxespad=0)
txt(2, yC + 75, "Panel D: row colours as in C; hatched cells,", va="top")
txt(2, yC + 79, "not informative (low depth); † authenticity under audit.", va="top")

# =============================== E ===============================
yE = hy + hh + 20.0
label("E", 0, yE)
pe = pd.read_csv(D + "panelE_breakpoint_histogram_pooled_DESCRIPTIVE.csv.gz")
pe = pe[pe.set == "noncanonical"].copy()
RB = 100
pe["b"] = ((pe.bin_start - 1) // 20 * 20 // RB)
agg = pe.groupby(["MOI", "hpi", "breakpoint", "b"]).reads_per_M_viral_frag.sum().reset_index()
nb = int(np.ceil(GL / RB))
wE, xE = 76, {0: 17, 1: 102}
tops = []
for ri, moi in enumerate(MOIS):
    for ci, bp in enumerate(("5prime_donor", "3prime_acceptor")):
        ax = ax_mm(xE[ci], yE + 9.5 + ri * 15.8, wE, 13)
        for t in TIMES:
            d = agg[(agg.MOI == moi) & (agg.hpi == t) & (agg.breakpoint == bp)]
            y = np.full(nb, np.nan)
            y[d.b.values.astype(int)] = d.reads_per_M_viral_frag.values
            ax.plot((np.arange(nb) + 0.5) * RB, y, color=TCOL[t], lw=0.6, alpha=0.9, marker="o", ms=1.0, mew=0)
        ax.set_yscale("log"); ax.set_ylim(0.03, 1500); ax.set_xlim(0, GL)
        ax.set_yticks([0.1, 10, 1000]); ax.set_yticklabels(["0.1", "10", "1000"]); ax.yaxis.set_minor_locator(NullLocator())
        ax.tick_params(labelsize=FS, length=2, pad=1)
        if ri == 0:
            ax.set_xticklabels([])
        else:
            ax.set_xticks([0, 10000, 20000])
            ax.set_xticklabels(["0", "10,000", "20,000"], fontsize=FS)
            ax.set_xlabel("PEDV genomic coordinate (nt)", fontsize=FSL, labelpad=1.5)
        if ci == 0:
            ax.set_ylabel(f"MOI = {moi:g}", fontsize=FS, labelpad=2)
        clean(ax)
        dd = agg[(agg.MOI == moi) & (agg.breakpoint == bp)].sort_values("reads_per_M_viral_frag", ascending=False)
        picked = []
        for r in dd.itertuples():
            pos = (r.b + 0.5) * RB
            if all(abs(pos - p) > 1000 for p in picked):
                picked.append(pos)
                tops.append(dict(MOI=moi, breakpoint=bp, bin_start=int(r.b * RB + 1), bin_end=int((r.b + 1) * RB),
                                 hpi=int(r.hpi), reads_per_M_viral_frag=float(r.reads_per_M_viral_frag)))
            if len(picked) == 3:
                break
pd.DataFrame(tops).to_csv(D + "panelE_top_peaks_by_panel.csv", index=False)
for ci in (0, 1):
    g_ = ax_mm(xE[ci], yE + 5.2, wE, 3.2)
    g_.set_xlim(0, GL); g_.set_ylim(0, 1); g_.axis("off")
    for name, (s, e) in GENES.items():
        g_.add_patch(Rectangle((s, 0), e - s, 1, fc=GENE_FILL[name], ec="#9A9A9A", lw=0.3))
        if name in ("ORF1a", "ORF1b", "S", "N"):
            g_.text((s + e) / 2, 0.5, name, ha="center", va="center", fontsize=FS)
    txt(xE[ci] + wE / 2, yE + 0.2, "5′ breakpoint (donor)" if ci == 0 else "3′ breakpoint (acceptor)", ha="center", va="top",
        fontsize=FS, fontweight="bold")
hE = [Line2D([0], [0], color=TCOL[t], lw=1.8, label=f"{t} hpi") for t in TIMES]
fig.legend(handles=hE, loc="lower center", bbox_to_anchor=(0.54, 0.0), ncol=5, frameon=False, fontsize=FS, handlelength=1.4,
           columnspacing=1.3, handletextpad=0.4, borderaxespad=0)
print("E bottom (mm) ~", yE + 9.5 + 16.5 + 13 + 14)

for ext, kw in (("png", {"dpi": 300}), ("pdf", {}), ("svg", {})):
    fig.savefig(F + f"Supplementary_Figure_junction_dynamics.{ext}", facecolor="white", **kw)
print("saved; size mm", W, H)
