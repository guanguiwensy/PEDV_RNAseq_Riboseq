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
import os, glob, numpy as np, pandas as pd
from scipy.stats import beta as betad
import g2b_common as G
PW = os.path.join(G.OUT, 'tables', 'per_window'); OUT = os.path.join(G.OUT, 'G2c_freeze_checks'); os.makedirs(OUT, exist_ok=True)
MB = [('10-29', 10, 30), ('30-99', 30, 100), ('100-299', 100, 300), ('>=300', 300, np.inf)]
NBOOT = 1000; rng = np.random.default_rng(G.SEED)


def cp(k, n):
    if n == 0: return np.nan, np.nan
    return (betad.ppf(.025, k, n - k + 1) if k > 0 else 0.0), (betad.ppf(.975, k + 1, n - k) if k < n else 1.0)


def gboot(p, genes, a):
    ug, inv = np.unique(genes, return_inverse=True)
    hit = np.bincount(inv, weights=(p <= a).astype(float), minlength=len(ug)); tot = np.bincount(inv, minlength=len(ug)).astype(float)
    W = rng.multinomial(len(ug), np.full(len(ug), 1 / len(ug)), size=NBOOT)
    est = (W @ hit) / np.maximum(W @ tot, 1)
    return np.percentile(est, 2.5), np.percentile(est, 97.5)


def calib_row(e, use, meta):
    q = e.p_EMP_Zs.values[use]; qs = e.pswap_EMP_Zs.values[use]; g = e.tx_idx.values[use]; n = len(q)
    r = dict(**meta, n_eval=len(e), n_assessable=n, n_genes=len(np.unique(g)), frac_out_or_NA=1 - n / max(len(e), 1))
    for a, nm in ((0.05, '05'), (0.01, '01'), (0.001, '001')):
        k = int((q <= a).sum()); r[f'le{nm}'] = k / n if n else np.nan
        r[f'le{nm}_cp_lo'], r[f'le{nm}_cp_hi'] = cp(k, n)
        r[f'le{nm}_gb_lo'], r[f'le{nm}_gb_hi'] = gboot(q, g, a) if n >= 20 else (np.nan, np.nan)
        r[f'swap_le{nm}'] = (qs <= a).sum() / n if n else np.nan
    ev = n >= 200
    r['liberal_point'] = bool(ev and (r['le05'] > 0.08 or r['le01'] > 0.02))
    r['liberal_supported'] = bool(ev and (r['le05_gb_lo'] > 0.05 or r['le01_gb_lo'] > 0.01))
    r['swap_liberal_point'] = bool(ev and (r['swap_le05'] > 0.08 or r['swap_le01'] > 0.02))
    r['verdict'] = ('not_evaluable(n<200)' if not ev else 'liberal' if (r['liberal_point'] or r['swap_liberal_point'])
                    else 'conservative' if r['le05'] < 0.02 else 'calibrated')
    return r


cal, eff, mc = [], [], []
for f in sorted(glob.glob(os.path.join(PW, 'calib_*_main.tsv.gz'))):
    P = pd.read_csv(f, sep='\t'); sid, L = P.sid.iloc[0], int(P.L.iloc[0])
    if (sid, L) not in G.QUALIFIED: continue
    P['assessable'] = ~P.p_EMP_Zs.isna()
    P['in_support'] = P.assessable & (P.EMP_Zs_radius <= 3) & (P.phi_extrapolated == 0) & (P.EMP_Zs_eff_genes >= 30) & (P.m < 300)
    for w, d in P.groupby('w'):
        for scope in ('all_assessable', 'in_support'):
            col = 'assessable' if scope == 'all_assessable' else 'in_support'
            for lab, lo, hi in [('all', 10, np.inf)] + MB:
                e = d[(d.m >= lo) & (d.m < hi)]
                if len(e) == 0: continue
                use = e[col].values
                meta = dict(scope=scope, sid=sid, label=G.LABEL[sid], L=L, w=w, depth_bin=lab)
                cal.append(calib_row(e, use, meta))
                x = e[use]
                if len(x):
                    eff.append(dict(**meta, n_units=len(x), frac_of_eval=len(x) / len(e),
                                    K_median=x.EMP_Zs_K.median(), K_p05=x.EMP_Zs_K.quantile(.05), K_min=x.EMP_Zs_K.min(),
                                    genes_in_matched_median=x.EMP_Zs_n_genes.median(), genes_in_matched_min=x.EMP_Zs_n_genes.min(),
                                    eff_genes_median=x.EMP_Zs_eff_genes.median(), eff_genes_p05=x.EMP_Zs_eff_genes.quantile(.05),
                                    eff_genes_min=x.EMP_Zs_eff_genes.min(), max_gene_frac_median=x.EMP_Zs_max_gene_frac.median(),
                                    max_gene_frac_max=x.EMP_Zs_max_gene_frac.max(), frac_radius_widened=(x.EMP_Zs_radius > 2).mean(),
                                    frac_phi_extrapolated=x.phi_extrapolated.mean(), min_attainable_p_median=(1 / (x.EMP_Zs_K + 1)).median(),
                                    frac_reach_BH69=(1 / (x.EMP_Zs_K + 1) <= 0.05 / 69).mean(),
                                    frac_reach_BH567=(1 / (x.EMP_Zs_K + 1) <= 0.05 / 567).mean()))
            e = d[d[col]]
            if len(e) == 0: continue
            pv = e.groupby('tx_idx').sample(1, random_state=G.SEED).p_EMP_Zs.values
            for N in (20, 69, 200, 567):
                if len(pv) < 2 * N: continue
                any_bh = any_bon = 0
                for _ in range(2000):
                    s = np.sort(rng.choice(pv, N, replace=False))
                    any_bh += (s <= 0.05 * np.arange(1, N + 1) / N).any(); any_bon += (s[0] <= 0.05 / N)
                mc.append(dict(scope=scope, sid=sid, label=G.LABEL[sid], L=L, w=w, N_tests=N, pool_windows=len(pv),
                               P_any_BH_q05=any_bh / 2000, P_any_Bonferroni=any_bon / 2000, min_p_in_pool=pv.min(),
                               n_pool_le_BH1=int((pv <= 0.05 / N).sum()), nominal_upper=0.05))

C = pd.DataFrame(cal); C.to_csv(os.path.join(OUT, 'EMP_Zs_calibration_per_library_length.tsv'), sep='\t', index=False, float_format='%.5g')
E = pd.DataFrame(eff); E.to_csv(os.path.join(OUT, 'EMP_Zs_effective_background_size.tsv'), sep='\t', index=False, float_format='%.5g')
M = pd.DataFrame(mc); M.to_csv(os.path.join(OUT, 'EMP_Zs_multiplicity_null_montecarlo.tsv'), sep='\t', index=False, float_format='%.5g')

sup = []
for (scope, sid, L, w), d in C[C.depth_bin != 'all'].groupby(['scope', 'sid', 'L', 'w']):
    cert = [b for b, lo, hi in MB[:-1] if ((d.depth_bin == b) & (d.verdict == 'calibrated') & ~d.liberal_supported).any()]
    contig = []
    for b, lo, hi in MB[:-1]:
        if b in cert: contig.append(b)
        else: break
    mmax = {'10-29': 29, '30-99': 99, '100-299': 299}[contig[-1]] if contig else None
    sup.append(dict(scope=scope, sid=sid, label=G.LABEL[sid], L=L, w=w, certified_bins=';'.join(cert) or 'none',
                    contiguous_certified=';'.join(contig) or 'none', certified_m_range=f'10-{mmax}' if mmax else 'none'))
S = pd.DataFrame(sup); S.to_csv(os.path.join(OUT, 'EMP_Zs_host_support_map.tsv'), sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 300)
print(S[S.scope == 'in_support'].pivot_table(index=['label', 'L'], columns='w', values='certified_m_range', aggfunc='first').to_string())
c = C[(C.scope == 'in_support')]
print(c[c.depth_bin == 'all'][['label', 'L', 'w', 'n_assessable', 'frac_out_or_NA', 'le05', 'le05_gb_lo', 'le05_gb_hi', 'le01', 'le01_gb_hi', 'le001', 'le001_cp_hi', 'swap_le05', 'verdict']].round(4).to_string(index=False))
print(c[(c.depth_bin != 'all') & (c.verdict == 'liberal')][['label', 'L', 'w', 'depth_bin', 'n_assessable', 'le05', 'le01', 'swap_le05', 'swap_le01', 'liberal_supported']])
print(M.groupby(['scope', 'w', 'N_tests'])[['P_any_BH_q05', 'P_any_Bonferroni']].agg(['min', 'max']).round(4))
