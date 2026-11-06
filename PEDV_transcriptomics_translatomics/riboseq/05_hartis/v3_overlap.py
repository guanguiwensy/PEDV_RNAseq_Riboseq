"""CDS 片段相交判断（v3_02 重叠剔除用）。

GTF 坐标为 1-based 闭区间 [start, end]。两个闭区间 [x1, y1] 与 [x2, y2] 相交，当且仅当 x1 <= y2 且 x2 <= y1
（共享端点，例如 y1 == x2，算相交；相邻但不共享碱基，例如 x2 == y1 + 1，不算相交）。

历史：v3_02 v1 写成了 `x1 <= y2 and x2 >= y1`（第二个不等号方向错误）。该式在两个区间"一前一后、互不相交"
时（x2 > y1）为真，在真正部分重叠时（x2 < y1）多半为假，因此把"基因跨度交错、CDS 片段不相交"的转录本判成重叠。
旧谓词保留为 `segments_intersect_v1_buggy`，只用于对账与单元测试。
"""
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
import collections


def segments_intersect(a, b):
    (x1, y1), (x2, y2) = a, b
    return x1 <= y2 and x2 <= y1


def segments_intersect_v1_buggy(a, b):
    (x1, y1), (x2, y2) = a, b
    return x1 <= y2 and x2 >= y1


def cds_overlap_set(cds_iv, tchrom, pred=segments_intersect):
    """cds_iv: {transcript: [(start, end), ...]}；tchrom: {transcript: chrom}。
    返回与另一条转录本（任一链）有 CDS 片段相交的转录本集合。按跨度排序扫描，只比较跨度相交的转录本对。"""
    span = {t: (tchrom[t], min(a for a, b in v), max(b for a, b in v)) for t, v in cds_iv.items()}
    bych = collections.defaultdict(list)
    for t, (c, a, b) in span.items(): bych[c].append((a, b, t))
    overlap = set()
    for c, L in bych.items():
        L.sort()
        for i in range(len(L)):
            a, b, t = L[i]; j = i + 1
            while j < len(L) and L[j][0] <= b:
                t2 = L[j][2]
                if any(pred(s1, s2) for s1 in cds_iv[t] for s2 in cds_iv[t2]):
                    overlap.add(t); overlap.add(t2)
                j += 1
    return overlap
