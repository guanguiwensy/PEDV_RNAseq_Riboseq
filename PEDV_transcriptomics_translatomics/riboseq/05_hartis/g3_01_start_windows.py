"""G3-01 起始近端窗口验证（单点窗 W1 与五密码子窗 W5 = s, s+3, s+6, s+9, s+12）。
记号：候选起点 s 在 F1（s = 3j+1），单点窗 = 密码子单元 j（n0=c[3j], n1=c[s], n2=c[s+1]）；W5 = 单元 j..j+4 之和。
(a) 宿主：模型（Mb53）只用训练基因 CDS 主体（j ≥ 10）估计；零分布 = 训练基因主体同类型窗口同 m 分层的 z；
    在留出基因上按"距 annotated 起始的密码子数 j"（0..40）评估 p ≤ 0.05 比例与 z 均值 —— 检查起始附近的特殊占据。
(b) 机制：F0 读段 5′ 端碱基（canonical 5′ 端，3j−12）对该密码子非经典框比例 q = m/(m+n0) 的影响（主体、训练基因）。
(c) 病毒候选起点：W1/W5 的 n0、m（n1/n2 盲藏）、是否在宿主支持范围内、距 canonical 起始；超出支持范围记 out_of_host_support，不算 p。
"""
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
import sample_registry as SR
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import g2_hostcodon as G
T = G.T; SHR = 20.0; GL = 27959
R = os.path.join(G.V3, '..'); RG = os.path.join(R, 'riboseq_GIb_mapping_psite')
res, mech, viral = [], [], []
V = pd.read_csv(os.path.join(R, 'HARTIS_design', 'tables', 'viral_plus1_start_architecture.tsv'), sep='\t')
nmap = pd.read_csv(os.path.join(R, 'global_perturbation_screen', 'plus1_orf_capacity', 'panelC_family_name_mapping.tsv'), sep='\t')
disp = dict(zip(nmap.original_family_id, nmap.displayed_name))
PSTART = {'ORF1a': 293, 'ORF1b': 12592, 'S': 20610, 'ORF3': 24758, 'E': 25361, 'M': 25608, 'N': 26300}
CANON = [293, 20610, 24758, 25361, 25608, 26300]
PS = pd.read_csv(os.path.join(RG, 'psite_counts_PEDV_GIb.tsv.gz'), sep='\t', comment='#', usecols=['sample_id', 'strand', 'position', 'read_length', 'psite_count'])
PS = PS[(PS.strand == '+') & (PS.read_length == 29)]
for lib in G.LIBS:
    D = G.codon_table(lib, include_start_region=True)
    D = D[(D.b5F1 < 4) & (D.b5F2 < 4) & (D.b3F1 < 4) & (D.b3F2 < 4)]
    body_tr = D[(~D.heldout) & (D.j >= 10) & (D.j <= D.N - 11)]
    trm = body_tr[body_tr.m > 0]; p0 = trm.n1.sum() / trm.m.sum()
    g = trm.groupby(['b5F1', 'b5F2', 'b3F1', 'b3F2'])[['n1', 'm']].sum(); g['pj'] = (g.n1 + SHR * p0) / (g.m + SHR)
    D = D.join(g.pj, on=['b5F1', 'b5F2', 'b3F1', 'b3F2']); D['pj'] = D.pj.fillna(p0)
    D = D.sort_values(['tx_idx', 'j']); D['E1'] = D.m * D.pj; D['Vv'] = D.m * D.pj * (1 - D.pj)
    # 构造 W1 / W5（W5 以 j 为起点的 5 个连续单元，滚动求和）
    for wt, w in [('W1', 1), ('W5', 5)]:
        agg = D.groupby('tx_idx')[['n0', 'n1', 'm', 'E1', 'Vv']].transform(lambda x: x[::-1].rolling(w, min_periods=w).sum()[::-1]) if w > 1 else D[['n0', 'n1', 'm', 'E1', 'Vv']]
        W = pd.concat([D[['tx_idx', 'j', 'N', 'heldout']], agg], axis=1).dropna()
        W = W[(W.m > 0) & (W.j + w - 1 <= W.N - 11)].copy()
        W['z'] = (W.n1 - W.E1) / np.sqrt(np.maximum(W.Vv, 1e-9))
        Wtr = W[(~W.heldout) & (W.j >= 10)]
        edges = np.unique(np.r_[np.quantile(Wtr.m, np.linspace(0, 1, 11))[:-1], np.quantile(Wtr.m, .99), np.inf])
        zz = {b: np.sort(d.z.values) for b, d in Wtr.assign(mb=np.searchsorted(edges, Wtr.m, side='right') - 1).groupby('mb')}
        Wte = W[W.heldout & (W.j <= 40)].copy(); Wte['mb'] = np.searchsorted(edges, Wte.m, side='right') - 1
        pe = []
        for r in Wte.itertuples():
            z_ = zz.get(r.mb); K = len(z_) if z_ is not None else 0
            pe.append((1 + K - np.searchsorted(z_, r.z, side='left')) / (1 + K) if K >= 100 else np.nan)
        Wte['p_emp'] = pe
        for jb, d in Wte.groupby(pd.cut(Wte.j, [-1, 0, 1, 2, 4, 6, 9, 14, 20, 30, 40])):
            x = d.p_emp.dropna()
            res.append(dict(library=lib, window=wt, codon_bin=str(jb), n=len(x), frac_p_le_0_05=(x <= .05).mean() if len(x) else np.nan,
                            frac_p_le_0_01=(x <= .01).mean() if len(x) else np.nan, mean_z=d.z.mean(),
                            F1share_obs=d.n1.sum() / max(d.m.sum(), 1), F1share_pred=d.E1.sum() / max(d.m.sum(), 1)))
        host_m_hi = np.quantile(Wtr.m, .995); host_n0_hi = np.quantile(Wtr.n0, .995)
        # (c) 病毒候选起点（盲藏 n1/n2）
        s_ = G.LIBS[lib]; a = np.zeros(GL + 8)
        if s_ == SR.lib("RIBO_VERO_CHX_1MOI_16h"):
            mm = pd.read_csv(os.path.join(RG, '_work', 'server', 'stats', s_, 'pedv_5p_matrix.tsv.gz'), sep='\t')
            mm = mm[(mm.strand == '+') & (mm.L == 29)]; gg = mm.groupby(mm.pos5 + 12).n.sum(); gg = gg[gg.index <= GL]
        else:
            gg = PS[PS.sample_id == s_].groupby('position').psite_count.sum()
        a[gg.index.values] = gg.values
        for v in V.itertuples():
            s0 = int(v.start_position); n0 = sum(a[s0 - 1 + 3 * k] for k in range(w)); m = sum(a[s0 + 3 * k] + a[s0 + 1 + 3 * k] for k in range(w))
            dc = min(abs(s0 - c) for c in CANON)
            sup = (m <= host_m_hi) and (n0 <= host_n0_hi)
            viral.append(dict(library=lib, window=wt, family_id=v.family_id, family_display=disp.get(v.family_id, ''), start_position=s0,
                              start_codon=v.start_codon, parent_ORF=v.parent_canonical_ORF, dist_to_canonical_start=dc, n0=int(n0), m=int(m),
                              n1='withheld_until_freeze', n2='withheld_until_freeze', host_m_p995=host_m_hi, host_n0_p995=host_n0_hi,
                              m_gate=('pass' if m >= 30 else ('low' if m >= 10 else 'fail')),
                              support=('in_host_support' if sup else 'out_of_host_support'),
                              in_canonical_start_shoulder=bool(any(-12 <= s0 - c <= 30 for c in CANON)),
                              detectability=('not_assessable' if m < 10 else ('out_of_host_support' if not sup else ('assessable' if m >= 30 else 'low_depth')))))
    # (b) 机制
    B = D[(~D.heldout) & (D.j >= 10) & (D.b5F0 < 4)]
    for base, nm in enumerate('ACGT'):
        d = B[B.b5F0 == base]
        mech.append(dict(library=lib, canonical_5p_base=nm, codons=len(d), offframe_frac=d.m.sum() / max((d.m + d.n0).sum(), 1),
                         F1share=d.n1.sum() / max(d.m.sum(), 1)))
    print(lib, 'done', flush=True)
Rr = pd.DataFrame(res); Rr.to_csv(os.path.join(T, 'g3_host_start_proximal_calibration.tsv'), sep='\t', index=False, float_format='%.4g')
M = pd.DataFrame(mech); M.to_csv(os.path.join(T, 'g3_mechanism_canonical_5p_base.tsv'), sep='\t', index=False, float_format='%.4g')
VV = pd.DataFrame(viral); VV.to_csv(os.path.join(T, 'g3_viral_candidate_start_windows_blinded.tsv'), sep='\t', index=False, float_format='%.4g')
pd.set_option('display.width', 260)
print(Rr.pivot_table(index=['window', 'codon_bin'], columns='library', values='frac_p_le_0_05', sort=False).round(3).to_string())
print(Rr.pivot_table(index=['window', 'codon_bin'], columns='library', values='mean_z', sort=False).round(2).to_string())
print(M.pivot_table(index='canonical_5p_base', columns='library', values='offframe_frac').round(3).to_string())
print(VV.groupby(['library', 'window', 'detectability']).size().unstack(fill_value=0).to_string())
