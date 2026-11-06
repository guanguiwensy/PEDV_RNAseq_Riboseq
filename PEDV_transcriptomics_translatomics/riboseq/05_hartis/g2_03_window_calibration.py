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
from scipy.stats import binom, betabinom, kstest, beta as betad
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import g2_hostcodon as G
T = G.T; SHR = 20.0


def cp(k, n):
    return (betad.ppf(.025, k, n - k + 1) if k > 0 else 0.0), (betad.ppf(.975, k + 1, n - k) if k < n else 1.0)


def add_p(D, tr):
    p0 = tr.loc[tr.m > 0, 'n1'].sum() / tr.loc[tr.m > 0, 'm'].sum()
    g = tr[tr.m > 0].groupby(['b5F1', 'b5F2', 'b3F1', 'b3F2'])[['n1', 'm']].sum()
    g['pj'] = (g.n1 + SHR * p0) / (g.m + SHR)
    D = D.join(g.pj, on=['b5F1', 'b5F2', 'b3F1', 'b3F2']); D['pj'] = D.pj.fillna(p0)
    return D, p0


def windows(D, wtype):
    D = D.sort_values(['tx_idx', 'j'])
    if wtype == 'single': key = D.index.values
    elif wtype == 'cod5': key = D.tx_idx.astype(np.int64) * 100000 + (D.j - 10) // 5
    elif wtype == 'cod50': key = D.tx_idx.astype(np.int64) * 100000 + (D.j - 10) // 50
    else: key = D.tx_idx.values
    D = D.assign(E1=D.m * D.pj, V=D.m * D.pj * (1 - D.pj), key=key)
    W = D.groupby('key').agg(tx_idx=('tx_idx', 'first'), n0=('n0', 'sum'), n1=('n1', 'sum'), n2=('n2', 'sum'), m=('m', 'sum'),
                             E1=('E1', 'sum'), V=('V', 'sum'), nc=('j', 'size'))
    if wtype == 'cod5': W = W[W.nc == 5]
    if wtype == 'cod50': W = W[W.nc == 50]
    W = W[W.m > 0].copy(); W['z'] = (W.n1 - W.E1) / np.sqrt(np.maximum(W.V, 1e-9)); W['ph'] = W.E1 / W.m
    W['q'] = W.m / (W.m + W.n0)
    return W


summ, detail = [], []
for lib in G.LIBS:
    D = G.codon_table(lib)
    D = D[(D.b5F1 < 4) & (D.b5F2 < 4) & (D.b3F1 < 4) & (D.b3F2 < 4)]
    D, p0 = add_p(D, D[~D.heldout])
    for wt in ['single', 'cod5', 'cod50', 'gene']:
        Wtr = windows(D[~D.heldout], wt); Wte = windows(D[D.heldout], wt)
        edges = np.unique(np.r_[np.quantile(Wtr.m, np.linspace(0, 1, 11))[:-1], np.quantile(Wtr.m, .99), np.inf])
        Wtr['mb'] = np.searchsorted(edges, Wtr.m, side='right') - 1; Wte['mb'] = np.searchsorted(edges, Wte.m, side='right') - 1
        eq = np.unique(np.r_[np.quantile(Wtr.n0, np.linspace(0, 1, 11))[:-1], np.inf])
        Wtr['nb'] = np.searchsorted(eq, Wtr.n0, side='right') - 1; Wte['nb'] = np.searchsorted(eq, Wte.n0, side='right') - 1
        ztr = {b: np.sort(d.z.values) for b, d in Wtr.groupby('mb')}
        rho = {}
        for b, d in Wtr.groupby('mb'):
            mb_ = d.m.mean(); rho[b] = max(0.0, (d.z.var() - 1) / max(mb_ - 1, 1)) if len(d) > 30 else 0.0
        qtr = {b: np.sort(d.q.values) for b, d in Wtr.groupby('nb')}
        pe, pb, pbin, pbin0, pq = [], [], [], [], []
        for r in Wte.itertuples():
            zz = ztr.get(r.mb); K = len(zz) if zz is not None else 0
            pe.append((1 + K - np.searchsorted(zz, r.z, side='left')) / (1 + K) if K >= 100 else np.nan)
            n1 = int(r.n1); m = int(r.m); ph = min(max(r.ph, 1e-6), 1 - 1e-6); rh = rho.get(r.mb, 0.0)
            if rh > 1e-6:
                s_ = (1 - rh) / rh; pb.append(betabinom.sf(n1 - 1, m, ph * s_, (1 - ph) * s_))
            else: pb.append(binom.sf(n1 - 1, m, ph))
            pbin.append(binom.sf(n1 - 1, m, ph)); pbin0.append(binom.sf(n1 - 1, m, p0))
            qq = qtr.get(r.nb); Kq = len(qq) if qq is not None else 0
            pq.append((1 + Kq - np.searchsorted(qq, r.q, side='left')) / (1 + Kq) if Kq >= 100 else np.nan)
        Wte['p_emp'] = pe; Wte['p_bb'] = pb; Wte['p_bin'] = pbin; Wte['p_bin0'] = pbin0; Wte['p_anomaly'] = pq
        Wte['library'] = lib; Wte['wtype'] = wt; detail.append(Wte[['library', 'wtype', 'tx_idx', 'n0', 'n1', 'n2', 'm', 'ph', 'z', 'mb', 'p_emp', 'p_bb', 'p_bin', 'p_bin0', 'p_anomaly']])
        groups = [('all', Wte)] + [(f'mbin{b}', d) for b, d in Wte.groupby('mb')]
        for gname, d in groups:
            for meth in ['p_emp', 'p_bb', 'p_bin', 'p_bin0', 'p_anomaly']:
                x = d[meth].dropna(); n = len(x)
                if n == 0: continue
                row = dict(library=lib, wtype=wt, stratum=gname, m_min=int(d.m.min()), m_max=int(d.m.max()), method=meth, n=n,
                           ks_p=kstest(x, 'uniform').pvalue if n > 20 else np.nan, median_p=x.median(),
                           train_K=len(ztr.get(int(gname[4:]), [])) if gname.startswith('mbin') else np.nan)
                for a in (.05, .01, .001):
                    k = int((x <= a).sum()); lo, hi = cp(k, n)
                    row[f'f_{a}'] = k / n; row[f'lo_{a}'] = lo; row[f'hi_{a}'] = hi
                summ.append(row)
        print(lib, wt, 'train', len(Wtr), 'test', len(Wte), 'm range', int(Wtr.m.min()), int(Wtr.m.max()), flush=True)
S = pd.DataFrame(summ); S.to_csv(os.path.join(T, 'g2_window_calibration_summary.tsv'), sep='\t', index=False, float_format='%.5g')
pd.concat(detail).to_csv(os.path.join(T, 'g2_window_calibration_detail.tsv.gz'), sep='\t', index=False, float_format='%.5g')
pd.set_option('display.width', 260)
a = S[S.stratum == 'all']
for f in ['f_0.05', 'f_0.01', 'f_0.001']:
    print(f); print(a.pivot_table(index=['library', 'wtype'], columns='method', values=f).round(4).to_string())
