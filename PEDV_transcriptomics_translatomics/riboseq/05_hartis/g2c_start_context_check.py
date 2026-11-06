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
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import g2b_common as G
from g2b_calibrate import emp_p, phi_fit, SHR, TRAIN_MMIN, EVAL_MMIN
OUT = os.path.join(G.OUT, 'G2c_freeze_checks'); os.makedirs(OUT, exist_ok=True)
NEAR = {'CTG', 'GTG', 'TTG', 'AAG', 'ACG', 'AGG', 'ATA', 'ATC', 'ATT'}
_B = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
CODE_AUG = 0 * 16 + 3 * 4 + 2
CODE_NEAR = np.array([_B[c[0]] * 16 + _B[c[1]] * 4 + _B[c[2]] for c in NEAR])
FRACS = (0.25, 0.5, 1.0); WINS = (1, 5); MAX_OTHER = 30000


def ctx_keys(C, mode):
    c = C.ctx.values
    if mode == 'ctx256': return c
    if mode == 'ctx3p16': return np.where(c == 256, 16, c % 16)
    return np.zeros(len(c), int)                                           # noctx


def rate_table(C, keys, train_mask):
    tr = train_mask & (C.m.values > 0)
    p0 = C.n1.values[tr].sum() / C.m.values[tr].sum()
    df = pd.DataFrame(dict(k=keys[tr], n1=C.n1.values[tr], m=C.m.values[tr])).groupby('k').sum()
    return (df.n1 + SHR * p0) / (df.m + SHR), p0


def rates(C, keys, train_mask):
    pc, p0 = rate_table(C, keys, train_mask)
    return pd.Series(keys).map(pc).fillna(p0).values, p0


def mode_key(b, mode):
    if mode == 'ctx256': return np.where((b == 4).any(0), 256, b[0] * 64 + b[1] * 16 + b[2] * 4 + b[3])
    if mode == 'ctx3p16': return np.where((b[2:] == 4).any(0), 16, b[2] * 4 + b[3])
    return np.zeros(b.shape[1], int)


def run(sid, L):
    tx = G.tx_table(); SEQ, ST = G._c['seq']; t = tx.set_index('tx_idx')
    C = G.codon_table(sid, L)

    keyN = t.cds_len_nt // 3
    rows_mu, rows_cal, rows_pow = [], [], []
    rng = np.random.default_rng(G.SEED)
    for f in range(G.NFOLD):
        trm = (C.fold.values != f)
        P = {}
        for mode in ('ctx256', 'ctx3p16', 'noctx'):
            keys = ctx_keys(C, mode)
            pj, p0 = rates(C, keys, trm)
            pj_cf = pj.copy()
            for g2 in range(G.NFOLD):
                if g2 == f: continue
                q, _ = rates(C, keys, (C.fold.values != f) & (C.fold.values != g2)); sel = C.fold.values == g2; pj_cf[sel] = q[sel]
            P[mode] = (pj, pj_cf, p0, rate_table(C, keys, trm)[0])
        for w in WINS:

            null = {}
            for mode, (pj, pj_cf, p0, _) in P.items():
                W = G.windows(C, w, pj_cf); W = W[(W.fold != f) & (W.m >= TRAIN_MMIN)].copy()
                z = (W.n1 - W.E1) / np.sqrt(W.V.clip(lower=1e-9)); lm = np.log2(W.m.values); ld = np.log2(W.dens.values + 1 / (3 * w))
                ph = phi_fit(lm, z.values)
                null[mode] = (z.values / np.sqrt(ph(lm)), lm, ld, ph)

            te = C[C.fold.values == f]; off = G.SPEC[sid][L]
            recs = []
            for k_tx, d in te.groupby('tx_idx'):
                N = int(keyN[k_tx])
                if N - 21 < w: continue
                arr = {c: np.zeros(N + 1) for c in ('n0', 'n1', 'n2')}
                for c in arr: arr[c][d.j.values] = d[c].values
                cs = int(t.loc[k_tx, 'cds_start_tx']); base = ST[k_tx] + cs
                jj = np.arange(N + 1); pos = base + 3 * jj + 1 - off
                bb = np.stack([SEQ[pos], SEQ[pos + 1], SEQ[pos + L - 1], SEQ[pos + L]])
                pjm = {mode: pd.Series(mode_key(bb, mode)).map(P[mode][3]).fillna(P[mode][2]).values for mode in P}
                tri_b = np.stack([SEQ[base + 3 * jj + 1], SEQ[base + 3 * jj + 2], SEQ[base + 3 * jj + 3]])
                js = np.arange(10, N - 10 - w + 1)
                if len(js) == 0: continue
                wsum = lambda a: (lambda c: c[js + w] - c[js])(np.concatenate([[0.0], np.cumsum(a)]))
                code = np.where((tri_b == 4).any(0), -1, tri_b[0] * 16 + tri_b[1] * 4 + tri_b[2])[js]
                cat = np.where(code == CODE_AUG, 'AUG', np.where(np.isin(code, CODE_NEAR), 'near_cognate', np.where(code < 0, 'N', 'other')))
                mj = arr['n1'] + arr['n2']
                r = dict(tx_idx=np.full(len(js), k_tx), j=js, cat=cat, n0=wsum(arr['n0']), n1=wsum(arr['n1']), n2=wsum(arr['n2']))
                r['m'] = r['n1'] + r['n2']
                for mode in P:
                    pv = pjm[mode]
                    r[f'E1_{mode}'] = wsum(mj * pv); r[f'V_{mode}'] = wsum(mj * pv * (1 - pv))
                    r[f'pjs_{mode}'] = pv[js]; r[f'mu_unw_{mode}'] = wsum(pv) / w
                recs.append(pd.DataFrame(r))
            if not recs: continue
            A = pd.concat(recs, ignore_index=True); A = A[A.cat != 'N']
            E = A[A.m >= EVAL_MMIN].copy()
            oth = A[A.cat == 'other']
            if len(oth) > MAX_OTHER: A = pd.concat([A[A.cat != 'other'], oth.sample(MAX_OTHER, random_state=G.SEED + f)])

            for cat, d in A.groupby('cat'):
                for mode in P:
                    rows_mu.append(dict(sid=sid, L=L, fold=f, w=w, cat=cat, mode=mode, n=len(d),
                                        mean_mu_minus_p0=(d[f'mu_unw_{mode}'] - P[mode][2]).mean(),
                                        mean_anchor_codon_mu_minus_p0=(d[f'pjs_{mode}'] - P[mode][2]).mean()))

            for mode in P:
                zs_tr, lm_tr, ld_tr, ph = null[mode]
                for frac in (0.0,) + FRACS:
                    D = (np.ceil(frac * E.m.values) if frac > 0 else np.zeros(len(E)))
                    m2 = E.m.values + D; n12 = E.n1.values + D
                    pj0 = E[f'pjs_{mode}'].values
                    E1 = E[f'E1_{mode}'].values + D * pj0; V = E[f'V_{mode}'].values + D * pj0 * (1 - pj0)
                    lm = np.log2(m2); ld = np.log2(E.n0.values / (3 * w) + 1 / (3 * w))
                    zs = (n12 - E1) / np.sqrt(np.maximum(V, 1e-9)) / np.sqrt(ph(lm))
                    p, K, _ = emp_p(zs_tr, lm_tr, ld_tr, zs, lm, ld, True)
                    for cat in ('AUG', 'near_cognate', 'other'):
                        s_ = (E.cat.values == cat) & ~np.isnan(p)
                        rec = dict(sid=sid, L=L, fold=f, w=w, mode=mode, spike_frac=frac, cat=cat, n=int(s_.sum()),
                                   k05=int((p[s_] <= 0.05).sum()), k01=int((p[s_] <= 0.01).sum()))
                        (rows_cal if frac == 0 else rows_pow).append(rec)
    return pd.DataFrame(rows_mu), pd.DataFrame(rows_cal), pd.DataFrame(rows_pow)


if __name__ == '__main__':
    sid, L = sys.argv[1], int(sys.argv[2])
    a, b, c = run(sid, L)
    for nm, d in (('mu', a), ('cal', b), ('pow', c)):
        d.to_csv(os.path.join(OUT, f'_start_ctx_{nm}_{sid}_{L}.tsv'), sep='\t', index=False)
    print(sid, L, len(a), len(b), len(c))
