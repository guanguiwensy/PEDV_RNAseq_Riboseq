#!/usr/bin/env python3
"""Fig.1D-F: RNA-seq病毒RNA积累动力学。
输入: figdata/fig1_counts.tsv, figdata/fig1_depth/{sample}.{raw,dedup}.depth.tsv, sample_manifest_GIb.tsv"""
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
import os
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap

D = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
F = f"{D}/figures"
man = pd.read_csv(f"{D}/sample_manifest_GIb.tsv", sep="\t")
man = man[man.sample_id.isin(SR.LIB_RNA_VERO)][["sample_id", "treatment", "time_h"]]
assert len(man) == 20
cnt = pd.read_csv(f"{D}/figdata/fig1_counts.tsv", sep="\t",
                  names=["sample_id", "ver", "total_mapped", "pedv_mapped", "r1_fwd", "r1_rev"])
qc = pd.read_csv(f"{D}/figures/RNA00S_read_composition_by_sample.csv")[["sample_id", "total_pairs"]]

# ---- 窗口 (1-based闭区间, refdata_GIb_v1.1坐标) ----
DEL1 = (3298, 3499)    # ORF1a 24 nt缺失(3398|3399)两侧各100 nt
DEL2 = (24901, 25100)  # ORF3 49 nt缺失(25000|25001)两侧各100 nt
W = {"5p_main": (1101, 20613, [DEL1]),   # ORF1a+ORF1b, 排除leader及其下游至~1 kb(1-1100)
     "5p_sens": (12592, 20613, []),      # ORF1b
     "3p_main": (20610, 27625, [DEL2]),  # S起点-N终点
     "3p_sens": (24759, 27625, [DEL2])}  # S终点之后(ORF3-N)
REG = {"Leader": (10, 72, []), "ORF1a": (293, 12591, [DEL1]), "ORF1b": (12592, 20609, []),
       "S": (20610, 24757, []), "ORF3/E/M": (24758, 26288, [DEL2]), "N/3′UTR": (26300, 27959, [])}


def wstat(dep, w):
    a, b, ex = w
    m = np.zeros(len(dep) + 1, bool)
    m[a:b + 1] = True
    for x, y in ex:
        m[x:y + 1] = False
    v = dep[m[1:]]
    return v.mean(), np.median(v), int(m[1:].sum())


rows = []
for _, s in man.iterrows():
    r = {"sample_id": s.sample_id, "treatment": s.treatment, "time_h": s.time_h}
    q = int(qc[qc.sample_id == s.sample_id].total_pairs.iloc[0])
    r["qc_read_pairs"] = q
    for ver in ["raw", "dedup"]:
        c = cnt[(cnt.sample_id == s.sample_id) & (cnt.ver == ver)].iloc[0]
        tag = "" if ver == "raw" else "_dedup"
        r[f"total_mapped_pairs{tag}"] = c.total_mapped
        r[f"pedv_mapped_pairs{tag}"] = c.pedv_mapped
        r[f"pedv_pct_of_mapped{tag}"] = 100 * c.pedv_mapped / c.total_mapped
        r[f"pedv_RPM{tag}"] = 1e6 * c.pedv_mapped / q
        if ver == "raw":
            r["r1_fwd_frac_proper_pairs"] = c.r1_fwd / (c.r1_fwd + c.r1_rev)
        dep = pd.read_csv(f"{D}/figdata/fig1_depth/{s.sample_id}.{ver}.depth.tsv", sep="\t", header=None)[2].values.astype(float)
        norm = c.pedv_mapped / 1e6
        for k, w in W.items():
            mean, med, n = wstat(dep, w)
            r[f"{k}_mean_cov{tag}"] = mean
            r[f"{k}_median_cov{tag}"] = med
            r[f"{k}_mean_cov_per_M_viral{tag}"] = mean / norm
            r[f"{k}_nt{tag}"] = n
        for fp in ["5p_main", "5p_sens"]:
            for tp in ["3p_main", "3p_sens"]:
                r[f"ratio_mean_{tp}_over_{fp}{tag}"] = r[f"{tp}_mean_cov{tag}"] / r[f"{fp}_mean_cov{tag}"]
                r[f"ratio_median_{tp}_over_{fp}{tag}"] = r[f"{tp}_median_cov{tag}"] / r[f"{fp}_median_cov{tag}"]
        for k, w in REG.items():
            mean, _, _ = wstat(dep, w)
            if ver == "raw":
                r[f"region_{k}_mean_cov"] = mean
                r[f"region_{k}_cov_per_M_mapped"] = mean / (c.total_mapped / 1e6)
            else:
                r[f"region_{k}_cov_per_M_mapped_dedup"] = mean / (c.total_mapped / 1e6)
    rows.append(r)
S = pd.DataFrame(rows).sort_values(["time_h", "treatment", "sample_id"]).reset_index(drop=True)
S["E_ratio_main"] = S["ratio_mean_3p_main_over_5p_main"]
S.to_csv(f"{F}/Fig1DEF_source_data_per_sample.csv", index=False)


def summ(col):
    g = S.groupby(["treatment", "time_h"])[col]
    o = pd.DataFrame({"mean": g.mean(), "sd": g.std(ddof=1), "n": g.count()})
    o["sem"] = o.sd / np.sqrt(o.n)
    o["min"] = g.min()
    o["max"] = g.max()
    return o.reset_index()


Ecols = ["E_ratio_main"] + [c for c in S.columns if c.startswith("ratio_")]
parts = []
for col in ["pedv_pct_of_mapped", "pedv_pct_of_mapped_dedup", "pedv_RPM"] + Ecols:
    o = summ(col)
    o.insert(0, "metric", col)
    parts.append(o)
SUM = pd.concat(parts)
SUM.to_csv(f"{F}/Fig1DEF_source_data_summary_by_condition.csv", index=False)

plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "font.size": 8})
STY = {"CHX": dict(c="#000000", m="o", off=-0.22), "HAR": dict(c="#D55E00", m="s", off=0.22)}
TIMES = [3, 5, 8, 12, 16]


def kin(ax, col, ylabel, ylim, yt):
    for t, st in STY.items():
        o = SUM[(SUM.metric == col) & (SUM.treatment == t)].sort_values("time_h")
        d = S[S.treatment == t]
        ax.scatter(d.time_h + st["off"], d[col], s=14, marker=st["m"], color=st["c"], alpha=0.4, linewidths=0, zorder=2)
        ax.errorbar(o.time_h.values + st["off"], o["mean"], yerr=o["sem"], color=st["c"], marker=st["m"], ms=4.5,
                    lw=1.4, elinewidth=1, capsize=2, label=t, zorder=3)
    ax.set_yscale("log")
    ax.set_ylim(*ylim)
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{v:g}" for v in yt])
    ax.minorticks_off()
    ax.set_xticks(TIMES)
    ax.set_xlim(2, 17)
    ax.set_xlabel("Hours post infection (hpi)", fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.tick_params(labelsize=8)
    for s_ in ["top", "right"]:
        ax.spines[s_].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="upper left", handlelength=1.2)


def panelD(ax):
    kin(ax, "pedv_pct_of_mapped", "PEDV (% of mapped reads)", (0.3, 100), [0.3, 1, 3, 10, 30, 100])


elo, ehi = S.E_ratio_main.min(), S.E_ratio_main.max()
ETICKS = [50, 100, 200, 400]


def panelE(ax):
    kin(ax, "E_ratio_main", "3′/5′ coverage ratio", (min(ETICKS), max(ETICKS)), ETICKS)
    if min(ETICKS) <= 1 <= max(ETICKS):
        ax.axhline(1, color="#8C8C8C", lw=0.9, ls="--", zorder=1)


Freg = list(REG)
Fcols = [f"region_{k}_cov_per_M_mapped" for k in Freg]
fm = S.groupby("time_h")[Fcols].mean().loc[TIMES]
L = np.log2(fm.T.values)
Lc = L - L.mean(axis=1, keepdims=True)
fout = pd.DataFrame(Lc, index=Freg, columns=[f"log2_rowcentered_{t}hpi" for t in TIMES])
fout.index.name = "region"
coords = pd.DataFrame({"region": Freg, "start": [REG[k][0] for k in Freg], "end": [REG[k][1] for k in Freg],
                       "excluded": [";".join(f"{a}-{b}" for a, b in REG[k][2]) for k in Freg]})
fraw = pd.DataFrame(fm.T.values, index=Freg, columns=[f"mean_cov_per_M_mapped_{t}hpi" for t in TIMES])
fraw.index.name = "region"
coords.merge(fraw.reset_index(), on="region").merge(fout.reset_index(), on="region").to_csv(
    f"{F}/Fig1F_regional_rna_accumulation.csv", index=False)
vmax = np.ceil(np.abs(Lc).max() * 2) / 2
CM = LinearSegmentedColormap.from_list("d", ["#2166AC", "#F7F7F7", "#B2182B"])


def panelF(ax, fig):
    im = ax.imshow(Lc, cmap=CM, norm=TwoSlopeNorm(0, -vmax, vmax), aspect="auto")
    ax.set_xticks(range(5))
    ax.set_xticklabels(TIMES, fontsize=8)
    ax.set_yticks(range(len(Freg)))
    ax.set_yticklabels(Freg, fontsize=8)
    ax.set_xlabel("Hours post infection (hpi)", fontsize=10)
    ax.tick_params(length=2)
    for i in range(Lc.shape[0]):
        for j in range(5):
            ax.text(j, i, f"{Lc[i, j]:.1f}", ha="center", va="center", fontsize=7,
                    color="white" if abs(Lc[i, j]) > vmax * 0.7 else "black")
    cb = fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.3, fraction=0.06)
    cb.set_label("Relative coverage\n(log2, row-centered)", fontsize=8)
    cb.ax.tick_params(labelsize=8)


def save(fig, name):
    fig.savefig(f"{F}/{name}.png", dpi=300, facecolor="white")
    fig.savefig(f"{F}/{name}.pdf", facecolor="white")
    plt.close(fig)


mm = 1 / 25.4
fig, ax = plt.subplots(figsize=(88 * mm, 70 * mm)); panelD(ax); fig.tight_layout(); save(fig, "Fig1D_viral_rna_kinetics")
fig, ax = plt.subplots(figsize=(88 * mm, 70 * mm)); panelE(ax); fig.tight_layout(); save(fig, "Fig1E_3to5_coverage_ratio")
fig, ax = plt.subplots(figsize=(88 * mm, 75 * mm)); panelF(ax, fig); fig.tight_layout(); save(fig, "Fig1F_regional_rna_accumulation")
fig = plt.figure(figsize=(180 * mm, 72 * mm))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.15], left=0.085, right=0.955, bottom=0.2, top=0.86, wspace=0.85)
a1, a2, a3 = [fig.add_subplot(g) for g in gs]
panelD(a1); panelE(a2); panelF(a3, fig)
for a, l, x in zip([a1, a2, a3], "DEF", [-0.42, -0.42, -0.72]):
    a.text(x, 1.08, l, transform=a.transAxes, fontsize=13, fontweight="bold")
save(fig, "Fig1DEF_combined")

pd.set_option("display.width", 250)
print(SUM[SUM.metric.isin(["pedv_pct_of_mapped", "E_ratio_main"])].round(3).to_string())
print(S[["sample_id", "treatment", "time_h", "pedv_pct_of_mapped", "pedv_RPM", "E_ratio_main", "r1_fwd_frac_proper_pairs"]].round(3).to_string())
print("E limits", elo, ehi, ETICKS)
print(fout.round(2))
