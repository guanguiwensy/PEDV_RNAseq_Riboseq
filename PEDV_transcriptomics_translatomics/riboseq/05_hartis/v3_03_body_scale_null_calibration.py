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
import os, numpy as np, pandas as pd
from scipy.stats import kstest, binom, betabinom, beta as betad
from scipy.optimize import minimize
H = os.path.dirname(os.path.abspath(__file__)); T = os.path.join(H, '..', 'tables')
SEED = 20261027; MMIN = 10; KMIN = 100; CAL = (0.5, 0.5, 0.03)
G = pd.read_csv(os.path.join(T, 'host_gene_level_29nt_frames.tsv.gz'), sep='\t')
G = G[G.train_eligible & (G.m >= MMIN) & G.GC.notna()].copy()
G['lm'] = np.log2(G.m); G['ld'] = np.log2(G.n0_density + 1e-4)


def cp(k, n):
    lo = betad.ppf(.025, k, n - k + 1) if k > 0 else 0.0; hi = betad.ppf(.975, k + 1, n - k) if k < n else 1.0
    return lo, hi


def fit_bb(k, n):
    def nll(t):
        a, b = np.exp(t); return -betabinom.logpmf(k, n, a, b).sum()
    p = k.sum() / n.sum(); r = minimize(nll, np.log([p * 20, (1 - p) * 20]), method='Nelder-Mead', options=dict(maxiter=400, xatol=1e-4, fatol=1e-4))
    return np.exp(r.x)


def evaluate(train, test, which='n1'):
    other = 'n2' if which == 'n1' else 'n1'
    tr_r = (train[which] / train.m).values; p0 = train[which].sum() / train.m.sum()
    X = train[['lm', 'ld', 'GC']].values
    out = []
    for _, t in test.iterrows():
        x = np.array([t.lm, t.ld, t.GC]); sel = None; lev = 0
        for lev in range(3):
            c = np.array(CAL) * (1.5 ** lev)
            k = (np.abs(X - x) <= c).all(1)
            if k.sum() >= KMIN: sel = np.where(k)[0]; break
        n1 = int(t[which]); m = int(t.m); r = n1 / m
        if sel is None:
            out.append(dict(K=int(k.sum()), level='insufficient', p_emp=np.nan, p_bb=np.nan, p_bin=binom.sf(n1 - 1, m, p0))); continue
        rg = tr_r[sel]; pe = (1 + (rg >= r).sum()) / (1 + len(sel))
        a, b = fit_bb(train[which].values[sel], train.m.values[sel])
        out.append(dict(K=len(sel), level=lev, p_emp=pe, p_bb=betabinom.sf(n1 - 1, m, a, b), p_bin=binom.sf(n1 - 1, m, p0),
                        bb_alpha=a, bb_beta=b))
    return pd.DataFrame(out, index=test.index)


rng = np.random.default_rng(SEED)
res = []
for lib, d in G.groupby('library'):
    d = d.copy(); idx = rng.permutation(d.index.values); half = len(idx) // 2
    tr, te = d.loc[idx[:half]], d.loc[idx[half:]]
    te = te.sample(min(len(te), 1500), random_state=SEED)
    for test_name, which, trn, tst in [('heldout_F1', 'n1', tr, te), ('swap_F2_as_candidate', 'n2', tr, te)]:
        e = evaluate(trn, tst, which); e['library'] = lib; e['test'] = test_name
        e = e.join(tst[['m', 'n0_density', 'transcript_id']]); res.append(e)

    flip_tr = rng.random(len(tr)) < .5; flip_te = rng.random(len(te)) < .5
    trp = tr.copy(); tep = te.copy()
    trp.loc[flip_tr, ['n1', 'n2']] = trp.loc[flip_tr, ['n2', 'n1']].values; tep.loc[flip_te, ['n1', 'n2']] = tep.loc[flip_te, ['n2', 'n1']].values
    e = evaluate(trp, tep, 'n1'); e['library'] = lib; e['test'] = 'label_permutation'; e = e.join(tep[['m', 'n0_density', 'transcript_id']]); res.append(e)
    print(lib, 'done', flush=True)
Rr = pd.concat(res); Rr.to_csv(os.path.join(T, 'calibration_body_gene_level_pvalues.tsv.gz'), sep='\t', index=False, float_format='%.5g')
summ = []
for (lib, test), d in Rr.groupby(['library', 'test']):
    subsets = {'all': d, 'top10pct_m': d[d.m >= d.m.quantile(.9)], 'top10pct_n0density': d[d.n0_density >= d.n0_density.quantile(.9)], 'm_lt_30': d[d.m < 30]}
    for sn, s in subsets.items():
        for meth in ['p_emp', 'p_bb', 'p_bin']:
            x = s[meth].dropna(); n = len(x)
            if n == 0: continue
            row = dict(library=lib, test=test, subset=sn, method=meth, n=n, frac_insufficient=float(s[meth].isna().mean()),
                       ks_uniform_p=kstest(x, 'uniform').pvalue if n > 10 else np.nan, median_p=x.median())
            for a in (.05, .01, .001):
                k = int((x <= a).sum()); lo, hi = cp(k, n)
                row[f'frac_le_{a}'] = k / n; row[f'ci_lo_{a}'] = lo; row[f'ci_hi_{a}'] = hi
            summ.append(row)
S = pd.DataFrame(summ); S.to_csv(os.path.join(T, 'calibration_body_gene_level_summary.tsv'), sep='\t', index=False, float_format='%.4g')
pd.set_option('display.width', 260)
print(S[S.subset == 'all'][['library', 'test', 'method', 'n', 'frac_le_0.05', 'frac_le_0.01', 'frac_le_0.001', 'ks_uniform_p', 'median_p']].to_string())
print(S[(S.subset != 'all') & (S.test == 'heldout_F1')][['library', 'subset', 'method', 'n', 'frac_le_0.05', 'frac_le_0.01', 'ks_uniform_p']].to_string())
