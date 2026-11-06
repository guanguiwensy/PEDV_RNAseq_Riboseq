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
p = PEDV_P+"Fig2_Supp_junction_dynamics/scripts/s4_plot_supp_figure.py"
s = open(p, encoding="utf-8").read()


def rep(a, b):
    global s
    assert a in s, a
    s = s.replace(a, b)


rep('(38.0, "Fraction of noncanonical\\njunction reads")', '(35.0, "Fraction of noncanonical\\njunction reads")')
rep("ax = ax_mm(GX[moi], ytop, GW, 19)", "ax = ax_mm(GX[moi], ytop, GW, 17)")
rep("-0.30, str(t)", "-0.27, str(t)")
rep('ax.text(-0.01, -0.30, "hpi"', 'ax.text(-0.01, -0.27, "hpi"')
rep("txt(GX[0.1] + GW + 4, 68.5,", "txt(GX[0.1] + GW + 4, 61.0,")
rep("yB = 74", "yB = 66")
rep("ax = ax_mm(GX[moi], yB + 8, GW, 26)", "ax = ax_mm(GX[moi], yB + 8, GW, 22)")
rep("yC = 116", "yC = 106")
rep("ax = ax_mm(16 + ci * 44, yC + 12, 31, 50)", "ax = ax_mm(16 + ci * 44, yC + 11, 31, 42)")
rep("bbox_to_anchor=(12 / W, 1 - (yC + 70) / H), ncol=2", "bbox_to_anchor=(12 / W, 1 - (yC + 66) / H), ncol=2")
rep("hx, hy, hw, hh = 118, yC + 12, 60, 62", "hx, hy, hw, hh = 122, yC + 11, 56, 54")
rep("axD.tick_params(length=0, pad=1)", "axD.tick_params(length=0, pad=9)")
rep("side = ax_mm(hx - 2.6, hy, 1.8, hh)", "side = ax_mm(hx - 2.4, hy, 1.8, hh)")
rep("cb = fig.add_axes([(hx + 6) / W, 1 - (hy + hh + 15.0) / H, 40 / W, 2.2 / H])",
    "cb = fig.add_axes([(hx + 2) / W, 1 - (hy + hh + 13.0) / H, 34 / W, 2.2 / H])")
rep('cb.set_xlabel("Row z-score of log2(reads per M viral fragments + 1)", fontsize=5.5, labelpad=1)',
    'cb.set_xlabel("Row z-score of log2(reads per M viral fragments + 1)", fontsize=5.3, labelpad=1)')
rep('txt(hx + 49, hy + hh + 13.4, "hatched: not\\ninformative (low depth)", fontsize=5.0, va="top", ha="left", color="#444")',
    'txt(hx + 38, hy + hh + 12.2, "hatched: not informative\\n(low depth)", fontsize=5.0, va="top", ha="left", color="#444")')
rep("bbox_to_anchor=((hx - 26) / W, 1 - (hy + hh + 22.5) / H), ncol=3",
    "bbox_to_anchor=((hx - 28) / W, 1 - (hy + hh + 20.0) / H), ncol=3")
rep('txt(hx - 26, hy + hh + 29.5, "† authenticity under audit (Fig. 2A, hotspot d).", fontsize=5.2, color="#444", va="top")',
    'txt(hx - 28, hy + hh + 25.0, "† authenticity under audit (Fig. 2A, hotspot d).", fontsize=5.2, color="#444", va="top")')
rep("yE = 213", "yE = 205")
rep("ax = ax_mm(xE[ci], yE + 10 + ri * 22, wE, 18)", "ax = ax_mm(xE[ci], yE + 10 + ri * 21, wE, 16)")
rep("if len(picked) == 3:", "if len(picked) == 2:")
rep("W, H = 180, 266", "W, H = 180, 262")
open(p, "w", encoding="utf-8").write(s)
print("patched")
