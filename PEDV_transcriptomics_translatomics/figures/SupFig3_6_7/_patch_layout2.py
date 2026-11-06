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


rep('txt(hx - 28, hy + hh + 25.0, "† authenticity under audit (Fig. 2A, hotspot d).", fontsize=5.2, color="#444", va="top")',
    'txt(hx - 28, hy + hh + 28.5, "† authenticity under audit (Fig. 2A, hotspot d).", fontsize=5.2, color="#444", va="top")')
rep("bbox_to_anchor=((hx - 28) / W, 1 - (hy + hh + 20.0) / H), ncol=3",
    "bbox_to_anchor=((hx - 28) / W, 1 - (hy + hh + 19.0) / H), ncol=2")
rep('''        for p in picked:
            ax.annotate(f"{int(p - RB/2 + 1)}", (p, 800), fontsize=4.8, ha="center", va="bottom", color="#333")''',
    '''        for k_, p in enumerate(picked):
            ax.annotate(f"{int(p - RB/2 + 1)}", (p, 800), fontsize=4.8, ha="right" if k_ == 0 else "left", va="bottom", color="#333")''')
open(p, "w", encoding="utf-8").write(s)
print("ok")
