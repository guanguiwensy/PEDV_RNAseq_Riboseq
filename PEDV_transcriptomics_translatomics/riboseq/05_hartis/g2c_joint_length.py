# -*- coding: utf-8 -*-
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
import os, sys, numpy as np, pandas as pd
from scipy.stats import chi2
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import g2b_common as G
from g2b_calibrate import emp_p, phi_fit, match_diag, SHR, TRAIN_MMIN, EVAL_MMIN
OUT = os.path.join(G.OUT, 'G2c_freeze_checks'); PW = os.path.join(G.OUT, 'tables', 'per_window')
WINS = (5, 50)


def length_windows(C, f, w):
    trc = C[(C.fold != f) & (C.m > 0)]; p0 = trc.n1.sum() / trc.m.sum()
    g = trc.groupby('ctx')[['n1', 'm']].sum(); pj = ((g.n1 + SHR * p0) / (g.m + SHR)).reindex(C.ctx).fillna(p0).values
    pj_cf = pj.copy()
    for g2 in range(G.NFOLD):
        if g2 == f: continue
        sub = C[(C.fold != f) & (C.fold != g2) & (C.m > 0)]; q0 = sub.n1.sum() / sub.m.sum()
        h = sub.groupby('ctx')[['n1', 'm']].sum(); qc = (h.n1 + SHR * q0) / (h.m + SHR)
        sel = (C.fold == g2).values; pj_cf[sel] = qc.reindex(C.ctx[sel]).fillna(q0).values
    W = G.windows(C, w, pj); Wc = G.windows(C, w, pj_cf)[['tx_idx', 'blk', 'E1', 'V']].rename(columns={'E1': 'E1cf', 'V': 'Vcf'})
    W = W.merge(Wc, on=['tx_idx', 'blk'])
    tr = (W.fold != f).values
    W['E1u'] = np.where(tr, W.E1cf, W.E1); W['Vu'] = np.where(tr, W.Vcf, W.V)
    W['z'] = (W.n1 - W.E1u) / np.sqrt(W.Vu.clip(lower=1e-9))
    t = W[tr & (W.m >= TRAIN_MMIN)]; ph = phi_fit(np.log2(t.m.values), t.z.values)
    W['zs'] = np.where(W.m > 0, W.z / np.sqrt(ph(np.log2(W.m.clip(lower=1)))), 0.0)
    return W[['tx_idx', 'blk', 'fold', 'n0', 'n1', 'n2', 'm', 'E1u', 'Vu', 'zs']]


rows, rowsF = [], []
for sid, spec in G.SPEC.items():
    Ls = [L for L in sorted(spec) if (sid, L) in G.QUALIFIED]
    if len(Ls) < 2: continue
    Cs = {L: G.codon_table(sid, L) for L in Ls}
    for w in WINS:
        recs = []
        for f in range(G.NFOLD):
            parts = []
            for L in Ls:
                Wl = length_windows(Cs[L], f, w).set_index(['tx_idx', 'blk'])
                Wl.columns = [f'{c}_{L}' if c != 'fold' else f'fold_{L}' for c in Wl.columns]; parts.append(Wl)
            J = pd.concat(parts, axis=1)
            J['fold'] = J[[f'fold_{L}' for L in Ls]].bfill(axis=1).iloc[:, 0]
            for c in ('n0', 'n1', 'n2', 'm', 'E1u', 'Vu', 'zs'):
                for L in Ls: J[f'{c}_{L}'] = J[f'{c}_{L}'].fillna(0)
            J['m'] = sum(J[f'm_{L}'] for L in Ls); J['n1'] = sum(J[f'n1_{L}'] for L in Ls); J['n2'] = sum(J[f'n2_{L}'] for L in Ls)
            J['n0'] = sum(J[f'n0_{L}'] for L in Ls); J['E1'] = sum(J[f'E1u_{L}'] for L in Ls); J['V'] = sum(J[f'Vu_{L}'] for L in Ls)
            J = J[J.m > 0].reset_index()
            J['z'] = (J.n1 - J.E1) / np.sqrt(J.V.clip(lower=1e-9))
            use = np.stack([(J[f'm_{L}'] >= TRAIN_MMIN).values for L in Ls]); k = use.sum(0)
            J['S'] = np.where(k > 0, sum(J[f'zs_{L}'].values * use[i] for i, L in enumerate(Ls)) / np.sqrt(np.maximum(k, 1)), 0.0)
            J['lm'] = np.log2(J.m); J['ld'] = np.log2(J.n0 / (3 * w) + 1 / (3 * w))
            tr = J[(J.fold != f) & (J.m >= TRAIN_MMIN)]; te = J[(J.fold == f) & (J.m >= EVAL_MMIN)].copy()
            ph = phi_fit(tr.lm.values, tr.z.values)
            zs_tr = tr.z.values / np.sqrt(ph(tr.lm.values)); zs_te = te.z.values / np.sqrt(ph(te.lm.values))

            lbv = np.floor(tr.lm.values / 0.25).astype(int); vc = pd.Series(lbv).value_counts(); okb = vc[vc >= 30].index.values
            xmin, xmax = (okb.min() + 0.5) * 0.25, (okb.max() + 0.5) * 0.25
            te['phi_extrapolated'] = ((te.lm.values < xmin) | (te.lm.values > xmax)).astype(int)
            dg = match_diag(tr.lm.values, tr.ld.values, tr.tx_idx.values, te.lm.values, te.ld.values, (2, 3, 5))
            for k_, v_ in dg.items(): te[f'J_{k_}'] = v_
            te['p_J_POOL'], _, _ = emp_p(zs_tr, tr.lm.values, tr.ld.values, zs_te, te.lm.values, te.ld.values, True)
            te['pswap_J_POOL'], _, _ = emp_p(zs_tr, tr.lm.values, tr.ld.values, zs_te, te.lm.values, te.ld.values, False)
            te['p_J_SUMZ'], _, _ = emp_p(tr.S.values, tr.lm.values, tr.ld.values, te.S.values, te.lm.values, te.ld.values, True)
            te['pswap_J_SUMZ'], _, _ = emp_p(tr.S.values, tr.lm.values, tr.ld.values, te.S.values, te.lm.values, te.ld.values, False)
            recs.append(te[['tx_idx', 'blk', 'm', 'p_J_POOL', 'pswap_J_POOL', 'p_J_SUMZ', 'pswap_J_SUMZ', 'phi_extrapolated',
                            'J_K', 'J_n_genes', 'J_max_gene_frac', 'J_eff_genes', 'J_radius']])
        T = pd.concat(recs)
        T.assign(sid=sid, w=w, lengths='/'.join(map(str, Ls))).to_csv(os.path.join(G.OUT, 'tables', 'per_window', f'joint_{sid}_w{w}.tsv.gz'), sep='	', index=False, float_format='%.6g')

        Pl = []
        for L in Ls:
            x = pd.read_csv(os.path.join(PW, f'calib_{sid}_{L}_main.tsv.gz'), sep='\t', usecols=['tx_idx', 'blk', 'w', 'p_EMP_Zs', 'pswap_EMP_Zs'])
            Pl.append(x[x.w == w].drop(columns='w').set_index(['tx_idx', 'blk']).add_suffix(f'_{L}'))
        F = pd.concat(Pl, axis=1, join='inner').dropna()
        for side, col in (('F1', 'p_EMP_Zs'), ('F2_swap', 'pswap_EMP_Zs')):
            X = -2 * sum(np.log(F[f'{col}_{L}']) for L in Ls); pf = chi2.sf(X, 2 * len(Ls))
            rowsF.append(dict(sid=sid, label=G.LABEL[sid], lengths='/'.join(map(str, Ls)), w=w, test=side, method='FISHER_naive',
                              n=len(pf), le05=(pf <= .05).mean(), le01=(pf <= .01).mean()))
        Tf = T.set_index(['tx_idx', 'blk']).loc[T.set_index(['tx_idx', 'blk']).index.intersection(F.index)]
        for meth in ('J_POOL', 'J_SUMZ'):
            for side, col in (('F1', f'p_{meth}'), ('F2_swap', f'pswap_{meth}')):
                for subset, D in (('all_m>=10', T), ('fisher_subset', Tf)):
                    q = D[col].dropna().values
                    rows.append(dict(sid=sid, label=G.LABEL[sid], lengths='/'.join(map(str, Ls)), w=w, test=side, method=meth, subset=subset,
                                     n=len(q), not_assessable=D[col].isna().mean(), le05=(q <= .05).mean(), le01=(q <= .01).mean()))
        print(sid, w, 'done', flush=True)
R = pd.concat([pd.DataFrame(rows), pd.DataFrame(rowsF).assign(subset='fisher_subset', not_assessable=0.0)])
R.to_csv(os.path.join(OUT, 'joint_length_calibration.tsv'), sep='\t', index=False, float_format='%.5g')
pd.set_option('display.width', 250)
print(R.pivot_table(index=['label', 'lengths', 'w'], columns=['method', 'subset', 'test'], values='le05').round(3).to_string())
