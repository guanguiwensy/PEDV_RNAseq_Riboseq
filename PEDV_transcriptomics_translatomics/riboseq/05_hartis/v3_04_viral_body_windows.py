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
import sample_registry as SR
import os, sys, numpy as np, pandas as pd
from scipy.stats import binom
H = os.path.dirname(os.path.abspath(__file__)); V3 = os.path.join(H, '..'); R = os.path.join(V3, '..'); T = os.path.join(V3, 'tables')
RG = os.path.join(R, 'riboseq_GIb_mapping_psite')
GL = 27959
ORFS = [('ORF1a', 293, 12622), ('ORF1b', 12592, 20613), ('S', 20610, 24758), ('ORF3', 24758, 25033),
        ('E', 25361, 25600), ('M', 25608, 26288), ('N', 26300, 27625)]
CSTART = {'ORF1a': 293, 'S': 20610, 'ORF3': 24758, 'E': 25361, 'M': 25608, 'N': 26300}
CSTOP = {'ORF1a': 12620, 'ORF1b': 20611, 'S': 24756, 'ORF3': 25031, 'E': 25598, 'M': 26286, 'N': 27623}
LIBS = {'HAR5': SR.lib("RIBO_VERO_HAR_1MOI_5h"), 'HAR8': SR.lib("RIBO_VERO_HAR_1MOI_8h"), 'HAR16': SR.lib("RIBO_VERO_HAR_1MOI_16h"), 'CHX5': SR.lib("RIBO_VERO_CHX_1MOI_5h"), 'CHX8': SR.lib("RIBO_VERO_CHX_1MOI_8h"),
        'CHX16': SR.lib("RIBO_VERO_CHX_1MOI_16h")}
HOSTLIB = {'HAR5': 'HAR5', 'HAR8': 'HAR8', 'HAR16': 'HAR16', 'CHX5': 'CHX5', 'CHX8': 'CHX8', 'CHX16': 'CHX16'}
SEQ = ' ' + ''.join(l.strip() for l in open(os.path.join(R, 'refbuild_GIb', 'v1.1', 'virus', 'PEDV_GIb.fa')) if not l.startswith('>')).upper()


PS = pd.read_csv(os.path.join(RG, 'psite_counts_PEDV_GIb.tsv.gz'), sep='\t', comment='#', usecols=['sample_id', 'strand', 'position', 'read_length', 'psite_count'])
PS = PS[(PS.strand == '+') & (PS.read_length == 29)]
A = {}
for lab, s in LIBS.items():
    a = np.zeros(GL + 2)
    if s == SR.lib("RIBO_VERO_CHX_1MOI_16h"):
        m = pd.read_csv(os.path.join(RG, '_work', 'server', 'stats', s, 'pedv_5p_matrix.tsv.gz'), sep='\t')
        m = m[(m.strand == '+') & (m.L == 29)]; g = m.groupby(m.pos5 + 12).n.sum(); g = g[g.index <= GL]
    else:
        g = PS[PS.sample_id == s].groupby('position').psite_count.sum()
    a[g.index.values] = g.values; A[lab] = a
pd.DataFrame([dict(library=k, sample_id=LIBS[k], plus_29nt_psites=int(v.sum())) for k, v in A.items()]).to_csv(os.path.join(T, 'viral_29nt_totals.tsv'), sep='\t', index=False)


parent = np.array([''] * (GL + 2), dtype=object); pst = np.zeros(GL + 2, int); cover = np.zeros(GL + 2, int)
for n, a, b in ORFS:
    cover[a:b + 1] += 1
    idx = np.arange(a, b + 1); free = parent[idx] == ''; parent[idx[free]] = n; pst[idx[free]] = a
excl = np.zeros(GL + 2, bool)
excl[12586:12623] = True
for s in CSTART.values(): excl[max(1, s - 12):s + 31] = True
for s in CSTOP.values(): excl[max(1, s - 30):min(GL, s + 30) + 1] = True
excl[26702 - 15:26702 + 16] = True
excl[cover != 1] = True
pos = np.arange(GL + 2)
frame = np.where(pst > 0, (pos - pst) % 3, -1)


def window_counts(lo, hi, par, a):
    p = np.arange(lo, hi + 1)
    p = p[(~excl[p]) & (parent[p] == par)]
    f = frame[p]; c = a[p]
    return len(p), c[f == 0].sum(), c[f == 1].sum(), c[f == 2].sum(), p


V = pd.read_csv(os.path.join(R, 'HARTIS_design', 'tables', 'viral_plus1_start_architecture.tsv'), sep='\t')
FAM = V.groupby('family_id').agg(parent=('parent_canonical_ORF', 'first'), first_start=('start_position', 'min'),
                                 stop=('stop_position', 'first'), n_starts=('start_position', 'size')).reset_index()
nm = pd.read_csv(os.path.join(R, 'global_perturbation_screen', 'plus1_orf_capacity', 'panelC_family_name_mapping.tsv'), sep='\t')
FAM['family_display'] = FAM.family_id.map(dict(zip(nm.original_family_id, nm.displayed_name)))


HG = pd.read_csv(os.path.join(T, 'host_gene_level_29nt_frames.tsv.gz'), sep='\t')
HG = HG[HG.train_eligible & (HG.m >= 10) & HG.GC.notna()].copy()
HG['lm'] = np.log2(HG.m); HG['ld'] = np.log2(HG.n0_density + 1e-4); HG['r'] = HG.n1 / HG.m
CAL = np.array([0.5, 0.5, 0.03]); KMIN = 100


def host_match(lib, m, dens, gc):
    if lib is None or m < 1: return None, 0, 'no_host_library'
    d = HG[HG.library == lib]; X = d[['lm', 'ld', 'GC']].values; x = np.array([np.log2(m), np.log2(dens + 1e-4), gc])
    for lev in range(3):
        k = (np.abs(X - x) <= CAL * 1.5 ** lev).all(1)
        if k.sum() >= KMIN: return d[k], int(k.sum()), f'level{lev}'

    out = []
    for j, nmv in enumerate(['m', 'n0_density', 'GC']):
        if x[j] > X[:, j].max() or x[j] < X[:, j].min(): out.append(nmv + '_out_of_host_range')
    return None, int(k.sum()), 'insufficient' + ('|' + ','.join(out) if out else '')


rows = []
for _, f in FAM.iterrows():
    lo, hi = f.first_start + 30, f.stop - 30
    for lab, a in A.items():
        nt, n0, n1, n2, p = window_counts(lo, hi, f.parent, a)
        m = int(n1 + n2); dens = n0 / nt if nt else 0
        gc = np.mean([SEQ[i] in 'GC' for i in p]) if nt else np.nan
        sub, K, lev = host_match(HOSTLIB[lab], m, dens, gc) if nt else (None, 0, 'empty_window')
        if sub is not None:
            rq = np.quantile(sub.r.values, .95); n1min = int(np.ceil(rq * m)); med = np.median(sub.r.values)
        else:
            rq = n1min = med = np.nan
        p0 = HG[HG.library == HOSTLIB[lab]].pipe(lambda d: d.n1.sum() / d.m.sum()) if HOSTLIB[lab] else np.nan
        rows.append(dict(scale='body', library=lab, sample_id=LIBS[lab], family_id=f.family_id, family_display=f.family_display, parent_ORF=f.parent,
                         window_lo=lo, window_hi=hi, usable_nt=nt, n0=int(n0), m=m, n1='withheld_until_freeze', n2='withheld_until_freeze',
                         n0_density=dens, window_GC=gc, matched_host_sites=K, match_level=lev,
                         null_median_F1_share=med, null_q95_F1_share=rq, MDE_F1_share_delta=(rq - med) if sub is not None else np.nan,
                         MDE_n1_min_at_p0_05=n1min, fixed_p0_library=p0,
                         m_gate=('pass' if m >= 30 else ('low' if m >= 10 else 'fail')),
                         detectability=('assessable' if (m >= 30 and sub is not None) else ('low_depth' if (m >= 10 and sub is not None) else 'not_assessable'))))

for _, v in V.iterrows():
    s0 = int(v.start_position)
    for w, (a_, b_) in {'W3': (-1, 1), 'W9': (-3, 5)}.items():
        for lab, a in A.items():
            p = np.arange(s0 + a_, s0 + b_ + 1); p = p[(parent[p] == v.parent_canonical_ORF) & (cover[p] == 1)]
            ex = bool(excl[s0])
            f = frame[p]; c = a[p]
            n0 = c[f == 0].sum(); m = c[f == 1].sum() + c[f == 2].sum()
            rows.append(dict(scale=f'start_{w}', library=lab, sample_id=LIBS[lab], family_id=v.family_id,
                             family_display=FAM.set_index('family_id').family_display.get(v.family_id), parent_ORF=v.parent_canonical_ORF,
                             start_position=s0, start_codon=v.start_codon, window_lo=s0 + a_, window_hi=s0 + b_, usable_nt=len(p),
                             n0=int(n0), m=int(m), n1='withheld_until_freeze', n2='withheld_until_freeze', in_excluded_interval=ex,
                             matched_host_sites=0, match_level='host_start_scale_training_unavailable_locally',
                             m_gate=('pass' if m >= 30 else ('low' if m >= 10 else 'fail')),
                             detectability=('not_assessable' if (ex or m < 10) else 'not_assessable_no_host_null')))
B = pd.DataFrame(rows); B.to_csv(os.path.join(T, 'viral_body_windows_candidates_blinded.tsv'), sep='\t', index=False, float_format='%.4g')
pd.set_option('display.width', 250)
print(B.groupby(['library', 'detectability']).size().unstack(fill_value=0))
print(B.groupby(['library', 'match_level']).size().unstack(fill_value=0).to_string())
print(B[B.family_display.isin(['N-x1', 'N-x2', 'ORF1a-x22'])][['library', 'family_display', 'usable_nt', 'n0', 'm', 'n0_density', 'matched_host_sites', 'match_level', 'detectability']].to_string())


occ = np.zeros(GL + 2, bool)
for _, f in FAM.iterrows(): occ[max(1, f.first_start - 60):min(GL, f.stop + 60) + 1] = True
for s in V.start_position.values: occ[max(1, s - 60):min(GL, s + 60) + 1] = True
nrows = []
for n, a0, b0 in ORFS:
    st = a0 + 30
    while st + 150 <= b0 - 30:
        lo, hi = st, st + 149; st += 150
        if occ[lo:hi + 1].any(): continue
        for lab, a in A.items():
            nt, n0, n1, n2, p = window_counts(lo, hi, n, a)
            if nt < 90: continue
            m = int(n1 + n2); dens = n0 / nt; gc = np.mean([SEQ[i] in 'GC' for i in p])
            sub, K, lev = host_match(HOSTLIB[lab], m, dens, gc) if m >= 10 else (None, 0, 'm_lt_10')
            row = dict(library=lab, parent=n, lo=lo, hi=hi, usable_nt=nt, n0=int(n0), n1=int(n1), n2=int(n2), m=m, n0_density=dens, GC=gc,
                       matched_host_sites=K, match_level=lev)
            if sub is not None:
                r1 = n1 / m; r2 = n2 / m
                row['p_emp_F1'] = (1 + (sub.r.values >= r1).sum()) / (1 + K)
                row['p_emp_F2_swap'] = (1 + ((1 - sub.r.values) >= r2).sum()) / (1 + K)
                p0 = HG[HG.library == HOSTLIB[lab]].pipe(lambda d: d.n1.sum() / d.m.sum())
                row['p_bin_F1'] = binom.sf(n1 - 1, m, p0); row['p_bin_F2'] = binom.sf(n2 - 1, m, 1 - p0)
                row['host_null_median'] = np.median(sub.r.values)
            nrows.append(row)
NC = pd.DataFrame(nrows); NC.to_csv(os.path.join(T, 'viral_noncandidate_body_windows.tsv'), sep='\t', index=False, float_format='%.5g')
print('noncandidate windows per lib:', NC.groupby('library').size().to_dict())
from scipy.stats import kstest
sm = []
for (lab, reg), d in NC.assign(region=np.where(NC.parent.isin(['ORF1a', 'ORF1b']), 'ORF1ab', 'structural')).groupby(['library', 'region']):
    x1 = d.p_emp_F1.dropna(); x2 = d.p_emp_F2_swap.dropna()
    sm.append(dict(library=lab, region=reg, n_windows=len(d), n_with_host_null=len(x1),
                   viral_F1_share_median=(d.n1 / d.m.replace(0, np.nan)).median(), host_null_median=d.host_null_median.median(),
                   F1_frac_p_le_0_05=(x1 <= .05).mean() if len(x1) else np.nan, F1_ks_p=kstest(x1, 'uniform').pvalue if len(x1) > 10 else np.nan,
                   F2swap_frac_p_le_0_05=(x2 <= .05).mean() if len(x2) else np.nan, F2swap_ks_p=kstest(x2, 'uniform').pvalue if len(x2) > 10 else np.nan,
                   binF1_frac_p_le_0_05=(d.p_bin_F1.dropna() <= .05).mean() if len(x1) else np.nan))
SM = pd.DataFrame(sm); SM.to_csv(os.path.join(T, 'viral_noncandidate_body_calibration_summary.tsv'), sep='\t', index=False, float_format='%.4g')
print(SM.to_string())
