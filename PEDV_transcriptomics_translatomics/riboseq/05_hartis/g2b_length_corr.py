# -*- coding: utf-8 -*-
"""同一文库不同读长在同一宿主窗口上的残差相关（只读 per_window 主口径），用于判断各读长能否视为独立证据。
注意：per_window 只含各读长 m ≥ 10 的留出窗口；两读长都 ≥ 10 的窗口才进入。统计量：z_cf（上下文标准化残差）、r − μ。"""
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
import os, glob, numpy as np, pandas as pd
from scipy.stats import spearmanr
import g2b_common as G
PW = os.path.join(G.OUT, 'tables', 'per_window'); rows = []
for sid, d in G.SPEC.items():
    Ls = sorted(d)
    if len(Ls) < 2: continue
    P = {L: pd.read_csv(os.path.join(PW, f'calib_{sid}_{L}_main.tsv.gz'), sep='\t', usecols=['tx_idx', 'blk', 'w', 'm', 'n1', 'mu', 'z_cf']) for L in Ls}
    for i in range(len(Ls)):
        for j in range(i + 1, len(Ls)):
            a, b = Ls[i], Ls[j]
            for w in (5, 50):
                M = P[a][P[a].w == w].merge(P[b][P[b].w == w], on=['tx_idx', 'blk'], suffixes=('_a', '_b'))
                if len(M) < 30: continue
                rz = spearmanr(M.z_cf_a, M.z_cf_b).correlation
                rr = spearmanr(M.n1_a / M.m_a - M.mu_a, M.n1_b / M.m_b - M.mu_b).correlation
                rows.append(dict(sid=sid, label=G.LABEL[sid], L_a=a, L_b=b, w=w, n_windows=len(M), n_genes=M.tx_idx.nunique(),
                                 spearman_z=rz, spearman_resid_share=rr))
pd.DataFrame(rows).to_csv(os.path.join(G.OUT, 'tables', 'G2b_cross_length_residual_corr.tsv'), sep='\t', index=False, float_format='%.4g')
print(pd.DataFrame(rows).round(3).to_string())
