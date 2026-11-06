# -*- coding: utf-8 -*-
"""J_POOL（多读长联合）的单独校准范围：与 EMP_Zs 相同的单元级支持规则与深度档认证（宿主按基因 5 折留出窗口）。
支持规则：K ≥ 100 且半径 ≤ 3 格（±0.75 log2）且 φ̂ 未外推且有效基因数 ≥ 30 且合并 m < 300。
输出 G2d_round5/JPOOL_calibration_by_depth.tsv、JPOOL_support_map.tsv。"""
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
from scipy.stats import beta as betad
import g2b_common as G
PW = os.path.join(G.OUT, 'tables', 'per_window'); OUT = os.path.join(G.OUT, 'G2d_round5'); os.makedirs(OUT, exist_ok=True)
MB = [('10-29', 10, 30), ('30-99', 30, 100), ('100-299', 100, 300), ('>=300', 300, np.inf)]
rng = np.random.default_rng(G.SEED)


def cp(k, n):
    if n == 0: return np.nan, np.nan
    return (betad.ppf(.025, k, n - k + 1) if k > 0 else 0.0), (betad.ppf(.975, k + 1, n - k) if k < n else 1.0)


def gb(p, g, a):
    ug, inv = np.unique(g, return_inverse=True)
    hit = np.bincount(inv, weights=(p <= a).astype(float)); tot = np.bincount(inv).astype(float)
    W = rng.multinomial(len(ug), np.full(len(ug), 1 / len(ug)), size=1000); est = (W @ hit) / np.maximum(W @ tot, 1)
    return np.percentile(est, 2.5), np.percentile(est, 97.5)


rows = []
for f in sorted(glob.glob(os.path.join(PW, 'joint_*_w*.tsv.gz'))):
    T = pd.read_csv(f, sep='\t'); sid, w, Ls = T.sid.iloc[0], int(T.w.iloc[0]), T.lengths.iloc[0]
    T['in_support'] = (~T.p_J_POOL.isna() & (T.J_radius <= 3) & (T.phi_extrapolated == 0) & (T.J_eff_genes >= 30) & (T.m < 300))
    for scope, col in (('all_assessable', None), ('in_support', 'in_support')):
        for lab, lo, hi in [('all', 10, np.inf)] + MB:
            e = T[(T.m >= lo) & (T.m < hi)]
            if len(e) == 0: continue
            use = ~e.p_J_POOL.isna() if col is None else e[col]
            q = e.p_J_POOL.values[use]; qs = e.pswap_J_POOL.values[use]; g = e.tx_idx.values[use]; n = len(q)
            r = dict(scope=scope, sid=sid, label=G.LABEL[sid], lengths=Ls, w=w, depth_bin=lab, n_eval=len(e), n_assessable=n,
                     frac_out_or_NA=1 - n / len(e), K_median=np.nanmedian(e.J_K.values[use]) if n else np.nan,
                     eff_genes_median=np.nanmedian(e.J_eff_genes.values[use]) if n else np.nan)
            for a, nm in ((0.05, '05'), (0.01, '01'), (0.001, '001')):
                k = int((q <= a).sum()); r[f'le{nm}'] = k / n if n else np.nan; r[f'le{nm}_cp_lo'], r[f'le{nm}_cp_hi'] = cp(k, n)
                r[f'swap_le{nm}'] = (qs <= a).sum() / n if n else np.nan
            r['le05_gb_lo'], r['le05_gb_hi'] = gb(q, g, .05) if n >= 20 else (np.nan, np.nan)
            r['le01_gb_lo'], r['le01_gb_hi'] = gb(q, g, .01) if n >= 20 else (np.nan, np.nan)
            ev = n >= 200
            lib = ev and (r['le05'] > 0.08 or r['le01'] > 0.02 or r['swap_le05'] > 0.08 or r['swap_le01'] > 0.02)
            r['liberal_supported'] = bool(ev and (r['le05_gb_lo'] > 0.05 or r['le01_gb_lo'] > 0.01))
            r['verdict'] = 'not_evaluable(n<200)' if not ev else 'liberal' if lib else 'conservative' if r['le05'] < 0.02 else 'calibrated'
            rows.append(r)
C = pd.DataFrame(rows); C.to_csv(os.path.join(OUT, 'JPOOL_calibration_by_depth.tsv'), sep='\t', index=False, float_format='%.5g')
sup = []
for (sid, w), d in C[(C.scope == 'in_support') & (C.depth_bin != 'all')].groupby(['sid', 'w']):
    contig = []
    for b, lo, hi in MB[:-1]:
        if ((d.depth_bin == b) & (d.verdict == 'calibrated') & ~d.liberal_supported).any(): contig.append(b)
        else: break
    sup.append(dict(sid=sid, label=G.LABEL[sid], lengths=d.lengths.iloc[0], w=w, contiguous_certified=';'.join(contig) or 'none',
                    certified_pooled_m_range={'10-29': '10-29', '30-99': '10-99', '100-299': '10-299'}[contig[-1]] if contig else 'none'))
S = pd.DataFrame(sup); S.to_csv(os.path.join(OUT, 'JPOOL_support_map.tsv'), sep='\t', index=False)
pd.set_option('display.width', 250)
print(S.to_string(index=False))
c = C[(C.scope == 'in_support') & (C.depth_bin == 'all')]
print(c[['label', 'lengths', 'w', 'n_assessable', 'frac_out_or_NA', 'le05', 'le05_gb_lo', 'le05_gb_hi', 'le01', 'le01_gb_hi', 'le001_cp_hi', 'swap_le05', 'swap_le01', 'verdict']].round(4).to_string(index=False))
