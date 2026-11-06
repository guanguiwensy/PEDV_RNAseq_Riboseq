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
H = os.path.dirname(os.path.abspath(__file__)); V3 = os.path.join(H, '..'); R = os.path.join(V3, '..'); T = os.path.join(V3, 'tables')
SEED = 20261027; rng = np.random.default_rng(SEED)
WIN = {'W3': (-1, 1), 'W6': (-1, 4), 'W9': (-3, 5), 'W15': (-3, 11), 'W30': (-3, 26)}
KS = [0.05, 0.1, 0.25, 0.5]
PL = pd.read_csv(os.path.join(R, 'host_start_proximal_HAR', 'results', 'host_position_level.tsv.gz'), sep='\t', comment='#')
PL = PL[PL.read_length == 29]
REL = np.arange(-30, 91)
rows, shp_rows = [], []
for lib, d in PL.groupby('library'):
    piv = d.pivot_table(index='gene_id', columns='rel_pos', values='psite_count', aggfunc='sum', fill_value=0).reindex(columns=REL, fill_value=0)
    M = piv.values.astype(float); ix = lambda r: r + 30
    agg = M.sum(0); shape = agg[ix(-3):ix(26) + 1]; shape = shape / shape.sum()      # rel −3..+26
    for k_, v in zip(range(-3, 27), shape): shp_rows.append((lib, k_, v))
    for w, (a, b) in WIN.items():
        xs = [x for x in range(33, 91) if x % 3 == 1 and x + b <= 90]
        if not xs: continue
        rel = np.arange(a, b + 1); fr = rel % 3
        n1, n2, n0, inj = [], [], [], {k: [] for k in KS}
        for x in xs:
            seg = M[:, ix(x + a):ix(x + b) + 1]
            c1 = seg[:, fr == 0].sum(1); c2 = seg[:, fr == 1].sum(1); c0 = seg[:, fr == 2].sum(1)
            n1.append(c1); n2.append(c2); n0.append(c0)
            sh = np.zeros(len(rel)); o = np.isin(rel, np.arange(-3, 27)); sh[o] = shape[rel[o] + 3]
            for k in KS:
                lam = (k * np.maximum(c0, 1))[:, None] * sh[None, :]
                add = rng.poisson(lam)
                inj[k].append((add[:, fr == 0].sum(1), add[:, fr == 1].sum(1)))
        n1 = np.concatenate(n1); n2 = np.concatenate(n2); n0 = np.concatenate(n0); m = n1 + n2
        ok = m >= 10
        lm = np.log2(np.maximum(m, 1)); r = np.where(m > 0, n1 / np.maximum(m, 1), np.nan)

        order = np.argsort(lm[ok]); lm_ok = lm[ok][order]; r_ok = r[ok][order]

        def pemp(lmt, rt):
            lo_, hi_ = np.searchsorted(lm_ok, lmt - .5), np.searchsorted(lm_ok, lmt + .5)
            seg = r_ok[lo_:hi_]
            if len(seg) < 100: return np.nan
            return (1 + (seg >= rt).sum()) / (1 + len(seg))
        neg_p = np.array([pemp(l_, r_) for l_, r_ in zip(lm[ok], r[ok])])
        rec = dict(library=lib, window=w, n_negative_sites=int(ok.sum()), neg_frac_p_le_0_05=np.nanmean(neg_p <= .05),
                   neg_frac_p_le_0_01=np.nanmean(neg_p <= .01), median_m=float(np.median(m[ok])))
        for k in KS:
            a1 = np.concatenate([q[0] for q in inj[k]]); a2 = np.concatenate([q[1] for q in inj[k]])
            m2 = m + a1 + a2; r2 = np.where(m2 > 0, (n1 + a1) / np.maximum(m2, 1), np.nan); ok2 = m2 >= 10
            pp = np.array([pemp(np.log2(mm), rr) for mm, rr in zip(m2[ok2], r2[ok2])])
            rec[f'sens_p0.05_k{k}'] = np.nanmean(pp <= .05); rec[f'sens_p0.01_k{k}'] = np.nanmean(pp <= .01)
        rows.append(rec)
    print(lib, 'done', flush=True)
Rr = pd.DataFrame(rows); Rr.to_csv(os.path.join(T, 'start_window_spikein_host.tsv'), sep='\t', index=False, float_format='%.4g')
pd.DataFrame(shp_rows, columns=['library', 'rel', 'shape_weight']).to_csv(os.path.join(T, 'host_start_profile_shape_29nt.tsv'), sep='\t', index=False, float_format='%.5g')
pd.set_option('display.width', 260)
print(Rr.to_string())
