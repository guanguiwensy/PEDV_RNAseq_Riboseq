import os as _os, sys as _sys
_d = _os.path.dirname(_os.path.abspath(__file__))
while _d != _os.path.dirname(_d) and not _os.path.isdir(_os.path.join(_d, "config")):
    _d = _os.path.dirname(_d)
_sys.path.insert(0, _os.path.join(_d, "config"))
from pedv_paths import PEDV_P, PEDV_K, PEDV_R, PEDV_COV, PEDV_SERVER_ROOT, REPO_ROOT, SCRIPT_DIRS, ORIGIN
_sys.path.extend(_x for _x in SCRIPT_DIRS if _x not in _sys.path)
__REPO_FILE__ = globals().get("__REPO_FILE__") or _os.path.abspath(__file__)
_o = ORIGIN.get(_os.path.relpath(_os.path.abspath(__file__), REPO_ROOT).replace(chr(92), "/"))
if _o:
    __file__ = _o
import collections


def segments_intersect(a, b):
    (x1, y1), (x2, y2) = a, b
    return x1 <= y2 and x2 <= y1


def segments_intersect_v1_buggy(a, b):
    (x1, y1), (x2, y2) = a, b
    return x1 <= y2 and x2 >= y1


def cds_overlap_set(cds_iv, tchrom, pred=segments_intersect):
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
