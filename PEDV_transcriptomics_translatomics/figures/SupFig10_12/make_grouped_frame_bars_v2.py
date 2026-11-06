#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""逐读长三读框并排条形图（v2）：CHX 16 h = RIBO_VERO_CHX_1MOI_16h；不再标 strict_core，仅标"通过宿主判据"。
数据：read_length_frame_periodicity 的 CSV（五个库）+ extra_RIBO_VERO_CHX_1MOI_16h 的 CSV（813），均为宿主 CDS-body（两端各去 45 nt）P-site。"""
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
import os, datetime
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.dirname(HERE); TASK = os.path.dirname(OUT)
FRM = f"{TASK}/read_length_frame_periodicity"; FIG = f"{OUT}/figures"; SD = f"{OUT}/source_data"
MM = 1 / 25.4; NOW = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
SAMPLES = [(SR.lib("RIBO_VERO_CHX_1MOI_5h"), "CHX", 5), (SR.lib("RIBO_VERO_CHX_1MOI_8h"), "CHX", 8), (SR.lib("RIBO_VERO_CHX_1MOI_16h"), "CHX", 16),
           (SR.lib("RIBO_VERO_HAR_1MOI_5h"), "HAR", 5), (SR.lib("RIBO_VERO_HAR_1MOI_8h"), "HAR", 8), (SR.lib("RIBO_VERO_HAR_1MOI_16h"), "HAR", 16)]
DRUGC = {"CHX": "#000000", "HAR": "#D55E00"}; FC = ["#7B3294", "#56B4E9", "#999999"]; LS = list(range(25, 36))
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.linewidth": 0.7, "svg.fonttype": "none", "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
                     "xtick.labelsize": 7, "ytick.labelsize": 7})
a = pd.read_csv(f"{FRM}/read_length_frame_counts_and_fractions.csv", comment="#")
SR.assert_table_whitelisted(a, "sample", "read_length_frame_counts"); SR.assert_exactly(set(a["sample"]), SR.RIBO_VERO_ALIASES, "read_length_frame_counts")
a = a[a.region == "CDS_body_45_45"]
keep = ["sample", "treatment", "time_hpi", "read_length", "total_host_CDS_psites", "N0", "N1", "N2", "F0", "F1", "F2", "top_to_second_frame_ratio",
        "psite_offset_applied", "offset_confidence", "host_criteria_pass"]
D = a[keep].copy()
D = D[D.read_length.between(25, 35)].copy()
D["host_criteria_pass"] = D.host_criteria_pass.astype(bool)
head = [f"# Reading-frame distribution of host CDS-body P-sites by read length (Figure 4 v2 libraries; CHX 16 h = {SR.lib('RIBO_VERO_CHX_1MOI_16h')})", "# version: FIG4v2_1.0",
        f"# generated: {NOW} (local)", "# P-site = 5' end + calibrated offset of the same library and length; CDS-body = canonical CDS minus first/last 45 nt (host ChlSab1.1)",
        "# host_criteria_pass: offset confidence high & >=5,000 body P-sites & F0 >= 0.50 & top/second frame ratio >= 1.5 (no reads removed by the flag)",
        "# script: scripts/make_grouped_frame_bars_v2.py"]
with open(f"{SD}/frame_fractions_by_length_v2.csv", "w", encoding="utf-8", newline="") as f:
    f.write("\n".join(head) + "\n"); D.to_csv(f, index=False, float_format="%.6g")

def nl(n): return f"{n/1e6:.1f}M" if n >= 1e6 else (f"{n/1e3:.0f}k" if n >= 1e3 else str(int(n)))
W, H = 185.0, 200.0
fig = plt.figure(figsize=(W * MM, H * MM))
PW, PH, LM, GX, TOP = 78.0, 36.0, 15.0, 11.0, 7.0
for sid, drug, t in SAMPLES:
    r = [5, 8, 16].index(t); c = 0 if drug == "CHX" else 1
    ax = fig.add_axes([(LM + c * (PW + GX)) / W, 1 - (TOP + r * (PH + 24) + PH) / H, PW / W, PH / H])
    d = D[D["sample"] == sid].sort_values("read_length"); x = np.arange(len(LS))
    for k, (col, fc) in enumerate(zip(["F0", "F1", "F2"], FC)):
        ax.bar(x + (k - 1) * 0.27, d[col].values, width=0.25, color=fc, linewidth=0)
    ax.axhline(1 / 3, color="#555555", lw=0.8, ls="--", zorder=0)
    ax.set_xlim(-0.7, len(LS) - 0.3); ax.set_ylim(0, 1.0); ax.set_xticks(x); ax.set_xticklabels([str(L) for L in LS])
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0]); ax.set_yticklabels(["0", "0.25", "0.5", "0.75", "1"])
    ax.set_title(f"{drug} {t} h", color=DRUGC[drug], fontweight="bold", fontsize=8, pad=2)
    if c == 0: ax.set_ylabel("Fraction of host CDS-body P-sites", fontsize=7.5)
    else: ax.tick_params(labelleft=False)
    ax.set_xlabel("Aligned footprint length (nt)", fontsize=7.5, labelpad=1)
    for i, rr in enumerate(d.itertuples()):
        ax.text(i, -0.30, nl(rr.total_host_CDS_psites), fontsize=7, ha="center", va="top")
        if rr.host_criteria_pass: ax.plot([i], [0.95], ls="", marker="o", ms=3.4, mfc="black", mec="black", zorder=5)
    ax.text(-0.75, -0.30, "n", fontsize=7, ha="right", va="top", style="italic")
hl = [Rectangle((0, 0), 1, 1, fc=c_) for c_ in FC] + [Line2D([], [], color="#555555", lw=0.8, ls="--"), Line2D([], [], ls="", marker="o", ms=3.4, mfc="black", mec="black")]
lb = ["Frame 0", "Frame 1", "Frame 2", "1/3 (no frame preference)", "passes host criteria: offset high, n ≥ 5,000, F0 ≥ 0.5, top/second frame ≥ 1.5"]
fig.legend(hl, lb, loc="lower left", bbox_to_anchor=(LM / W, 0.0), ncol=2, frameon=False, fontsize=7, handlelength=1.4, columnspacing=1.6)
for ext, kw in (("pdf", {}), ("svg", {}), ("png", dict(dpi=300))):
    try: fig.savefig(f"{FIG}/Fig4v2_frame_periodicity_grouped_bars.{ext}", facecolor="white", **kw)
    except PermissionError: print("LOCKED", ext)
plt.close(fig)
print(D[D.host_criteria_pass].groupby("sample").read_length.apply(list).to_dict())
