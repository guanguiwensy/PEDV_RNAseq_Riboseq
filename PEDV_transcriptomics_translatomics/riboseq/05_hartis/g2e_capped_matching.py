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
import g2d_transfer_feasibility as TF
OUT = os.path.join(G.OUT, 'G2e_final'); os.makedirs(OUT, exist_ok=True)
SUP = pd.read_csv(os.path.join(G.OUT, 'G2c_freeze_checks', 'EMP_Zs_host_support_map.tsv'), sep='\t')
SUP = SUP[SUP.scope == 'in_support'].set_index(['sid', 'L', 'w']).certified_m_range
BIN = {'10-29': (10, 30), '10-99': (30, 100), '10-299': (100, 300)}
WINS = (1, 5)
val, vir, virsum = [], [], []
for sid, L in sorted(G.QUALIFIED):
    lab = G.LABEL[sid]; meta = dict(sid=sid, label=lab, drug='HAR' if 'HAR' in lab else 'CHX', L=L)
    C, pj, pc_all, p0 = TF.host_reference(sid, L)
    D = TF.viral_codons(sid, L)
    D['mu'] = pd.Series(D.ctx).map(pc_all).fillna(p0).values; D['E1'] = D.m * D.mu; D['V'] = D.m * D.mu * (1 - D.mu)
    for w in WINS:
        HW = G.windows(C, w, pj); HW = HW[HW.m >= 10].copy(); HW['d'] = (HW.n1 - HW.E1) / HW.m
        cert = SUP.get((sid, L, w), 'none')

        if cert in ('10-99', '10-299'):
            tr = HW[HW.fold.isin([0, 1, 2]) & (HW.m >= 30) & (HW.m < 100)].d.values
            for nm, lo, hi in (('same_bin_30-99', 30, 100), ('one_bin_up_100-299', 100, 300)):
                te = HW[HW.fold.isin([3, 4]) & (HW.m >= lo) & (HW.m < hi)].d.values
                if len(te) < 30 or len(tr) < 100: continue
                srt = np.sort(tr); p = (1 + len(srt) - np.searchsorted(srt, te - 1e-12, side='left')) / (1 + len(srt))
                ps = (1 + np.searchsorted(srt, te + 1e-12, side='right')) / (1 + len(srt))
                val.append(dict(**meta, w=w, test=nm, K_ref=len(tr), n_test=len(te), le05=(p <= .05).mean(), le01=(p <= .01).mean(),
                                swap_le05=(ps <= .05).mean(), sd_ref=tr.std(), sd_test=te.std()))

        if cert == 'none':
            virsum.append(dict(**meta, w=w, certified_m_range='none', status='not_assessable(no_certified_host_bin)')); continue
        lo, hi = BIN[cert]; ref = HW[(HW.m >= lo) & (HW.m < hi)].d.values; mref, sref = ref.mean(), ref.std(); srt = np.sort(ref)
        VW = TF.viral_windows(D, w); VW = VW[VW.m >= hi].copy()
        if len(VW) == 0:
            virsum.append(dict(**meta, w=w, certified_m_range=cert, status='no_viral_window_above_cap')); continue
        VW['d'] = (VW.n1 - VW.E1) / VW.m; VW['Z_cap'] = (VW.d - mref) / sref
        VW['p_cap_upper'] = (1 + len(srt) - np.searchsorted(srt, VW.d.values - 1e-12, side='left')) / (1 + len(srt))
        VW['p_cap_lower'] = (1 + np.searchsorted(srt, VW.d.values + 1e-12, side='right')) / (1 + len(srt))
        VW['flag'] = 'capped_reference'; VW['ref_bin'] = f'{lo}-{hi - 1}'; VW['K_ref'] = len(ref)
        for k, v in meta.items(): VW[k] = v
        VW['w'] = w; vir.append(VW)
        for orf, d in VW.groupby('ORF', sort=False):
            virsum.append(dict(**meta, w=w, certified_m_range=cert, ref_bin=f'{lo}-{hi - 1}', K_ref=len(ref), ORF=orf, n_windows=len(d),
                               m_median=d.m.median(), m_max=d.m.max(), Z_cap_median=d.Z_cap.median(), Z_cap_q10=d.Z_cap.quantile(.1), Z_cap_q90=d.Z_cap.quantile(.9),
                               upper_tail_le05=(d.p_cap_upper <= .05).mean(), lower_tail_le05=(d.p_cap_lower <= .05).mean(),
                               min_attainable_p=1 / (len(ref) + 1), flag='capped_reference', status='exploratory'))
    print(sid, L, flush=True)
pd.DataFrame(val).to_csv(os.path.join(OUT, 'capped_host_validation.tsv'), sep='\t', index=False, float_format='%.5g')
pd.concat(vir).to_csv(os.path.join(OUT, 'capped_pedv_noncandidate_windows.tsv.gz'), sep='\t', index=False, float_format='%.5g')
pd.DataFrame(virsum).to_csv(os.path.join(OUT, 'capped_pedv_noncandidate_summary.tsv'), sep='\t', index=False, float_format='%.5g')
print('done')
