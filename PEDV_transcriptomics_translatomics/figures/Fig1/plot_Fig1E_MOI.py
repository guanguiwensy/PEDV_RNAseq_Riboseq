#!/usr/bin/env python3
"""Figure 1E (3'/5' coverage ratio) 按 MOI 重绘。
复用 plot_Fig1_MOI.py 的配色/版式(与新 D 同版式)与 plot_Fig1DEF.py 的 E 口径(3p_main / 5p_main, 未去重 raw, MAPQ=255)。
数据: data/FigE_ratio_per_sample_MOI.csv(逐样本, 已含 MOI 列)。零重复: 不画误差棒, 不跨 MOI 平均。"""
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
import matplotlib.pyplot as plt
import plot_Fig1_MOI as P
from plot_Fig1_MOI import FIG, MM, CHX_C, HAR_C, MOI_LAB, TIMES, clean, save

E = pd.read_csv(f"{P.OUT}/data/FigE_ratio_per_sample_MOI.csv")
E = E.sort_values(["MOI", "hpi", "treatment"]).reset_index(drop=True)
assert len(E) == 20 and (E.groupby(["MOI", "treatment", "hpi"]).size() == 1).all()
WIN = "3′ window = S start–N end (20,610–27,625) minus ORF3 deletion flanks (24,901–25,100); " \
      "5′ window = ORF1a+ORF1b from 1,101 to 20,613 minus ORF1a deletion flanks (3,298–3,499); mean depth ratio"


def kin(ax, series, col, ylim, yt):
    for lab, d, c, m, filled, ls, off in series:
        d = d.sort_values("hpi")
        ax.plot(d.hpi + off, d[col], color=c, lw=1.4, ls=ls, marker=m, ms=5,
                mfc=c if filled else "white", mec=c, mew=1.2, label=lab, zorder=3)
    ax.set_yscale("log"); ax.set_ylim(*ylim)
    ax.set_yticks(yt); ax.set_yticklabels([f"{v:g}" for v in yt]); ax.minorticks_off()
    ax.set_xticks(TIMES); ax.set_xlim(2, 17)
    ax.set_xlabel("Hours post infection (hpi)", fontsize=10)
    ax.set_ylabel("3′/5′ coverage ratio", fontsize=10)
    clean(ax)


def series_both(col):
    ser = []
    for moi, ls, filled, moff in ((1.0, "-", True, 0.12), (0.1, "--", False, -0.12)):
        for tr, c, m, toff in (("CHX", CHX_C, "o", -0.2), ("HAR", HAR_C, "s", 0.2)):
            ser.append((f"{tr} · MOI {MOI_LAB[moi]}", E[(E.MOI == moi) & (E.treatment == tr)], c, m, filled, ls, toff + moff))
    return ser


def figE():
    ycols = {"E_ratio_raw": "3′/5′ ratio (non-deduplicated, uniquely mapped)",
             "E_ratio_dedup": "3′/5′ ratio (coordinate-deduplicated; sensitivity)"}
    cols = ["sample_id", "MOI", "treatment", "hpi", "E_ratio_raw", "E_ratio_dedup"]
    ET, ETD = [50, 100, 200, 400], [20, 50, 100, 200]

    # ---- main: MOI = 1, raw ----
    d1 = E[E.MOI == 1.0]
    fig, ax = plt.subplots(figsize=(88 * MM, 72 * MM))
    kin(ax, [("CHX", d1[d1.treatment == "CHX"], CHX_C, "o", True, "-", -0.2),
             ("HAR", d1[d1.treatment == "HAR"], HAR_C, "s", True, "-", 0.2)], "E_ratio_raw", (40, 500), ET)
    ax.legend(frameon=False, fontsize=8, loc="upper right", handlelength=1.4)
    ax.text(0.98, 0.04, "MOI = 1", transform=ax.transAxes, ha="right", fontsize=8, fontweight="bold")
    ax.set_title("3′/5′ coverage ratio", fontsize=10, fontweight="bold", pad=8)
    fig.tight_layout()
    save(fig, "Fig1E_main_MOI1", 88, 72, "log y 40-500 (data 100-277); one sample per point; no error bars; raw")
    d1[cols].assign(window_note=WIN).to_csv(f"{FIG}/Fig1E_main_MOI1_source_data.csv", index=False)

    # ---- supplement: both MOI, raw ----
    fig, ax = plt.subplots(figsize=(95 * MM, 75 * MM))
    kin(ax, series_both("E_ratio_raw"), "E_ratio_raw", (40, 500), ET)
    ax.legend(frameon=False, fontsize=7.5, loc="upper right", handlelength=2.0, ncol=1)
    ax.set_title("3′/5′ coverage ratio, by MOI", fontsize=10, fontweight="bold", pad=8)
    fig.tight_layout()
    save(fig, "Fig1E_suppl_bothMOI", 95, 75, "4 series; filled solid = MOI 1, open dashed = MOI 0.1; no error bars; raw")
    E[cols].assign(window_note=WIN).to_csv(f"{FIG}/Fig1E_suppl_bothMOI_source_data.csv", index=False)

    # ---- supplement: dedup sensitivity ----
    fig, ax = plt.subplots(figsize=(95 * MM, 75 * MM))
    kin(ax, series_both("E_ratio_dedup"), "E_ratio_dedup", (15, 300), [20, 50, 100, 200])
    ax.legend(frameon=False, fontsize=7.5, loc="upper right", handlelength=2.0, ncol=1)
    ax.set_title("3′/5′ coverage ratio, deduplicated", fontsize=10, fontweight="bold", pad=8)
    fig.tight_layout()
    save(fig, "Fig1E_suppl_dedup_bothMOI", 95, 75, "sensitivity: coordinate dedup (no UMI); 4 series; no error bars")
    E[cols].assign(window_note=WIN).to_csv(f"{FIG}/Fig1E_suppl_dedup_bothMOI_source_data.csv", index=False)


if __name__ == "__main__":
    figE()
    print("E done; ranges raw", E.E_ratio_raw.min().round(1), E.E_ratio_raw.max().round(1),
          "dedup", E.E_ratio_dedup.min().round(1), E.E_ratio_dedup.max().round(1))
