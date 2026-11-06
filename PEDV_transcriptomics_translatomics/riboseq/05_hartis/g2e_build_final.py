# -*- coding: utf-8 -*-
"""G2 定稿 HTML：G2_FINAL_HARTIS_v3.md + 预注册修订 3 + 关键图（复用各轮图）。"""
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
import os, re, base64, markdown
import g2b_common as G
F = [('G2c_freeze_checks/figures/G2c_per_combo_calibration.png', '图 1 · EMP_Zs 逐库×读长宿主校准（支持范围内；按基因 bootstrap 95%）'),
     ('G2d_round5/figures/G2d_host_AUG_profiles.png', '图 2 · 宿主 +1 框 AUG 的起始相关占据：集中在 AUG 的 P-site，HAR ≫ CHX'),
     ('G2d_round5/figures/G2d_transfer_overview.png', '图 3 · PEDV 非候选区：周期性与宿主支持（迁移现状）')]


def nest(md):
    out = []
    for line in md.split(chr(10)):
        m = re.match(r'^( +)([-*]|\d+\.) ', line)
        if m:
            n = len(m.group(1)); line = ' ' * (4 * (1 if n <= 4 else 2)) + line[n:]
        out.append(line)
    return chr(10).join(out)


md = lambda p: markdown.markdown(nest(open(p, encoding='utf-8').read()), extensions=['tables'])
img = ''.join(f'<figure><img alt="{c}" src="data:image/png;base64,{base64.b64encode(open(os.path.join(G.OUT, f), "rb").read()).decode()}"><figcaption>{c}</figcaption></figure>' for f, c in F)
CSS = """:root{--surface:#fcfcfb;--panel:#f4f3f0;--ink:#0b0b0b;--ink2:#52514e;--line:#e3e2dd;--accent:#2a78d6}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--surface:#1a1a19;--panel:#232321;--ink:#f0efec;--ink2:#c3c2b7;--line:#3a3a37;--accent:#3987e5}}
:root[data-theme="dark"]{--surface:#1a1a19;--panel:#232321;--ink:#f0efec;--ink2:#c3c2b7;--line:#3a3a37;--accent:#3987e5}
body{background:var(--surface);color:var(--ink);font:15px/1.65 -apple-system,"Segoe UI","Microsoft YaHei",sans-serif;margin:0;padding:24px 16px}
main{max-width:1150px;margin:0 auto}h1{font-size:24px}h2{font-size:19px;border-bottom:1px solid var(--line);padding-bottom:4px;margin-top:32px}
nav{background:var(--panel);border-radius:8px;padding:10px 14px}nav a{color:var(--accent);margin-right:18px}
table{border-collapse:collapse;font-size:12.5px;margin:8px 0;display:block;overflow-x:auto}th,td{border-bottom:1px solid var(--line);padding:4px 7px;vertical-align:top}th{color:var(--ink2);text-align:left}
code{background:var(--panel);padding:1px 4px;border-radius:3px;font-size:12.5px}figure{margin:16px 0}figure img{width:100%;background:#fff;border-radius:6px}figcaption{color:var(--ink2);font-size:13px}
.part{border-top:3px solid var(--accent);margin-top:40px;padding-top:6px}"""
html = f"""<!doctype html><html lang="zh"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>HARTIS_v3 G2 定稿</title><style>{CSS}</style></head><body><main>
<nav><a href="#a">定稿说明</a><a href="#f">图</a><a href="#p">预注册修订 3</a></nav>
<section id="a" class="part">{md(os.path.join(G.V3, 'G2_FINAL_HARTIS_v3.md'))}</section>
<section id="f" class="part"><h2>图</h2>{img}</section>
<section id="p" class="part">{md(os.path.join(G.V3, 'PREREGISTRATION_HARTIS_v3_rev3_G2_FINAL.md'))}</section>
</main></body></html>"""
out = os.path.join(G.V3, 'report_G2_FINAL_HARTIS_v3.html'); open(out, 'w', encoding='utf-8').write(html); print(out, len(html))
