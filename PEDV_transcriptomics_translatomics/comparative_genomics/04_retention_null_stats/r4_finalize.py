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
import hashlib,glob,os,markdown,base64
O=PEDV_K+r"Supplementary_Figure_L/"
css="body{font-family:'Microsoft YaHei',Arial,sans-serif;max-width:980px;margin:2em auto;padding:0 1em;line-height:1.6;color:#111;background:#fff}table{border-collapse:collapse;margin:1em 0}td,th{border:1px solid #bbb;padding:4px 8px;font-size:14px}th{background:#f2f2f2}code{background:#f4f4f4;padding:1px 4px}img{max-width:100%}"
t=open(O+"REVISION_AUDIT_REPORT.md",encoding="utf-8").read()
t=t.replace("168 条中 154 条属 GII-c（旧命名中的 S-INDEL），只有 7 条同时是 GI-b","168 条中 156 条属 GII-c（旧命名中的 S-INDEL）、7 条为 GI-b、4 条为 GI-a、1 条 unresolved")
open(O+"REVISION_AUDIT_REPORT.md","w",encoding="utf-8").write(t)
imgs="".join("<p><img src='data:image/png;base64,%s'></p>"%base64.b64encode(open(O+"work/rev_page%d.png"%k,'rb').read()).decode() for k in (1,2))
r=open(O+"RESULTS_en.md",encoding="utf-8").read()
open(O+"REVISION_AUDIT_REPORT.html","w",encoding="utf-8").write("<!doctype html><meta charset='utf-8'><title>Supplementary Figure L revised audit</title><style>%s</style>%s<h2>Results (EN)</h2>%s<h2>Figure preview</h2>%s"%(css,markdown.markdown(t,extensions=['tables']),markdown.markdown(r),imgs))
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
out=["%s  %s"%(sha(p),os.path.basename(p)) for p in sorted(glob.glob(O+"*.tsv")+glob.glob(O+"*.csv")+glob.glob(O+"Supplementary_Figure_L_REVISED*"))]
open(O+"PROVENANCE_SHA256_REVISED_outputs.txt","w").write("\n".join(out))
open(O+"README.md","a",encoding="utf-8").write("\n## REVISED (2026-10-12)\nStart with `REVISION_AUDIT_REPORT.md`. Revised figure: `Supplementary_Figure_L_REVISED*`. Audit tables: `SupL_*.tsv`. The first-version figure and tables are kept; where they differ, REVISED supersedes (see `SupL_old_vs_revised.tsv`). Re-run order: p3a, p3b, p5a, p6, p7, r2_audit_tables, r3_old_vs_revised, f1_make_figure, f2_svg_check, r4_finalize.\n")
print("ok")
