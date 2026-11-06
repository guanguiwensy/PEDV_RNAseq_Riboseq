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


rep('layers = [(9.5, "Fraction of all\\nclassified reads"), (31.5, "Fraction of\\nnoncanonical")]',
    'layers = [(9.5, "Fraction of all\\nclassified reads"), (33.0, "Fraction of\\nnoncanonical")]')
rep("ax.text((xpos[2 * i] + xpos[2 * i + 1]) / 2, -0.62, str(t)", "ax.text((xpos[2 * i] + xpos[2 * i + 1]) / 2, -0.33, str(t)")
rep('ax.text(-0.012, -0.62, "hpi"', 'ax.text(-0.012, -0.33, "hpi"')
rep("yB = 52", "yB = 55")
rep("yC = 92", "yC = 96")
rep("ax = ax_mm(16 + ci * 40, yC + 6, 28, 44)", "ax = ax_mm(19 + ci * 42, yC + 6, 25, 42)")
rep("bbox_to_anchor=(2 / W, 1 - (yC + 56) / H), ncol=1, frameon=False, fontsize=FS,\n           handlelength=1.8, labelspacing=0.2, borderaxespad=0)",
    "bbox_to_anchor=(2 / W, 1 - (yC + 55) / H), ncol=2, frameon=False, fontsize=FS,\n           handlelength=1.6, labelspacing=0.2, columnspacing=1.0, borderaxespad=0)")
rep('''fig.legend(handles=hD, loc="upper left", bbox_to_anchor=(2 / W, 1 - (yC + 73) / H), ncol=1, frameon=False, fontsize=FS,
           handlelength=1.0, labelspacing=0.2, borderaxespad=0, title="Panel D", title_fontsize=FS, alignment="left")
txt(2, yC + 73 + 21.5, "Hatched: not informative (low depth)", va="top")
txt(2, yC + 73 + 25.5, "† authenticity under audit", va="top")''',
    '''fig.legend(handles=[Patch(color=TR_COL["CHX"], label="CHX (column strip)"), Patch(color=TR_COL["HAR"], label="HAR (column strip)")],
           loc="upper left", bbox_to_anchor=(2 / W, 1 - (yC + 70) / H), ncol=2, frameon=False, fontsize=FS,
           handlelength=1.0, columnspacing=1.0, borderaxespad=0)
txt(2, yC + 75, "Panel D: row colours as in C; hatched cells, not informative (low depth); † authenticity under audit.", va="top")''')
open(p, "w", encoding="utf-8").write(s)
print("ok")
