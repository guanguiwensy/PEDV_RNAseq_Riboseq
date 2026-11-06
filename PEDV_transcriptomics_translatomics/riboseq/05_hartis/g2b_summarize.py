# -*- coding: utf-8 -*-
"""G2b 汇总：按深度分档的校准、对调/置换检验、深度描述、深度效应回归（按基因整群 bootstrap）。"""
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
from scipy.stats import beta as betad, kstest
import g2b_common as G
T = os.path.join(G.OUT, 'tables'); PW = os.path.join(T, 'per_window')
MODELS = ['BIN_GLOBAL', 'BIN_CTX', 'BB_RHOm', 'BB_CTX_RHOm', 'BB_CTX_cf', 'EMP_R', 'EMP_R_mid', 'EMP_Z', 'EMP_Z_cf', 'EMP_Z_cf_mid', 'EMP_Zs', 'EMP_Zs_tight', 'EMP_Zs_noctx', 'EMP_Zs_gcap']
SWAP = ['BB_CTX_RHOm', 'BB_CTX_cf', 'EMP_R', 'EMP_R_mid', 'EMP_Z', 'EMP_Z_cf', 'EMP_Z_cf_mid', 'EMP_Zs', 'EMP_Zs_tight', 'EMP_Zs_noctx', 'EMP_Zs_gcap']
MB = [10, 30, 100, 300, np.inf]; MBL = ['10-29', '30-99', '100-299', '>=300']


def cp(k, n):
    lo = betad.ppf(.025, k, n - k + 1) if k > 0 else 0.0; hi = betad.ppf(.975, k + 1, n - k) if k < n else 1.0
    return lo, hi


def summ(p, genes):
    ok = ~np.isnan(p); q = p[ok]; n = len(q); r = dict(n_eval=len(p), n_assessable=n, frac_not_assessable=1 - n / max(len(p), 1),
                                                          n_genes=len(set(genes[ok])))
    for a in (0.05, 0.01, 0.001):
        k = int((q <= a).sum()); lo, hi = cp(k, n) if n else (np.nan, np.nan)
        r[f'le{a}'] = k / n if n else np.nan; r[f'le{a}_lo'] = lo; r[f'le{a}_hi'] = hi
    r['median_p'] = float(np.median(q)) if n else np.nan
    r['KS_p'] = float(kstest(q, 'uniform').pvalue) if n >= 20 else np.nan
    return r


def verdict(r):
    if r['n_assessable'] < 200: return 'not_evaluable(n<200)'
    ok05 = 0.02 <= r['le0.05'] <= 0.08; ok01 = r['le0.01'] <= 0.02
    if ok05 and ok01: return 'calibrated'
    if r['le0.05'] > 0.08 or r['le0.01'] > 0.02: return 'liberal'
    return 'conservative'


rows = []; desc = []
for f in sorted(glob.glob(os.path.join(PW, 'calib_*.tsv.gz'))):
    P = pd.read_csv(f, sep='\t')
    sid, L, var = P.sid.iloc[0], int(P.L.iloc[0]), P.variant.iloc[0]
    for w, d in P.groupby('w'):
        q99 = d.m.quantile(0.99)
        bins = [('all', d)] + [(lab, d[(d.m >= lo) & (d.m < hi)]) for lab, lo, hi in zip(MBL, MB[:-1], MB[1:])] + [('top1%m', d[d.m >= q99])]
        for lab, e in bins:
            if len(e) == 0: continue
            for mdl in MODELS:
                r = summ(e[f'p_{mdl}'].values, e.tx_idx.values); r.update(sid=sid, label=G.LABEL[sid], L=L, variant=var, w=w, depth_bin=lab,
                                                                     m_min=int(e.m.min()), m_max=int(e.m.max()), model=mdl, test='F1')
                r['verdict'] = verdict(r); rows.append(r)
            for mdl in SWAP:
                r = summ(e[f'pswap_{mdl}'].values, e.tx_idx.values); r.update(sid=sid, label=G.LABEL[sid], L=L, variant=var, w=w, depth_bin=lab,
                                                                         m_min=int(e.m.min()), m_max=int(e.m.max()), model=mdl, test='F2_swap')
                r['verdict'] = verdict(r); rows.append(r)
R = pd.DataFrame(rows)
R['qualified_length'] = [(s, L) in G.QUALIFIED for s, L in zip(R.sid, R.L)]
cols = ['sid', 'label', 'L', 'qualified_length', 'variant', 'w', 'depth_bin', 'm_min', 'm_max', 'model', 'test', 'n_eval', 'n_assessable',
        'frac_not_assessable', 'n_genes', 'le0.05', 'le0.05_lo', 'le0.05_hi', 'le0.01', 'le0.01_lo', 'le0.01_hi', 'le0.001', 'le0.001_lo',
        'le0.001_hi', 'median_p', 'KS_p', 'verdict']
R[cols].to_csv(os.path.join(T, 'G2b_calibration_by_depth.tsv'), sep='\t', index=False, float_format='%.5g')

D = pd.concat([pd.read_csv(f, sep='\t') for f in glob.glob(os.path.join(PW, 'desc_*.tsv'))])
D['label'] = D.sid.map(G.LABEL)
D.to_csv(os.path.join(T, 'G2b_depth_description.tsv'), sep='\t', index=False, float_format='%.5g')


# ---------- 深度效应回归：logit P(F1) = offset(上下文) + b0 + b1·log2 基因密度 + b2·log2(窗口密度/基因密度) ----------
def irls(y, m, X, off, it=30):
    b = np.zeros(X.shape[1])
    for _ in range(it):
        eta = off + X @ b; mu = 1 / (1 + np.exp(-eta)); W = m * mu * (1 - mu)
        z = eta - off + (y - m * mu) / np.maximum(W, 1e-9)
        b_new = np.linalg.solve(X.T @ (X * W[:, None]) + 1e-8 * np.eye(X.shape[1]), X.T @ (W * z))
        if np.max(np.abs(b_new - b)) < 1e-7: b = b_new; break
        b = b_new
    return b


reg = []; rng = np.random.default_rng(G.SEED)
for f in sorted(glob.glob(os.path.join(PW, 'calib_*_main.tsv.gz'))):
    P = pd.read_csv(f, sep='\t'); sid, L = P.sid.iloc[0], int(P.L.iloc[0])
    for w in (5, 50):
        d = P[(P.w == w)].copy()
        if len(d) < 300: continue
        eps = 1.0 / (3 * w)
        d['lg'] = np.log2(d.gene_dens + eps); d['ll'] = np.log2((d.dens + eps) / (d.gene_dens + eps)); d['lm'] = np.log2(d.m)
        off = np.log(d.mu.clip(1e-4, 1 - 1e-4) / (1 - d.mu.clip(1e-4, 1 - 1e-4))).values
        for spec, cols_ in (('gene+local', ['lg', 'll']), ('m_only', ['lm'])):
            X = np.column_stack([np.ones(len(d))] + [d[c].values for c in cols_])
            b = irls(d.n1.values.astype(float), d.m.values.astype(float), X, off)
            genes = d.tx_idx.values; ug = np.unique(genes); idx = {g: np.where(genes == g)[0] for g in ug}
            bs = []
            for _ in range(200):
                s = np.concatenate([idx[g] for g in rng.choice(ug, len(ug), replace=True)])
                bs.append(irls(d.n1.values[s].astype(float), d.m.values[s].astype(float), X[s], off[s]))
            bs = np.array(bs)
            for i, nm in enumerate(['intercept'] + cols_):
                reg.append(dict(sid=sid, label=G.LABEL[sid], L=L, w=w, spec=spec, term=nm, coef=b[i], ci_lo=np.percentile(bs[:, i], 2.5),
                                ci_hi=np.percentile(bs[:, i], 97.5), n_windows=len(d), n_genes=len(ug),
                                note='系数为 log-odds/每 log2 单位；offset = 上下文期望'))
pd.DataFrame(reg).to_csv(os.path.join(T, 'G2b_depth_regression.tsv'), sep='\t', index=False, float_format='%.5g')
print('done', len(R), len(D), len(reg))
