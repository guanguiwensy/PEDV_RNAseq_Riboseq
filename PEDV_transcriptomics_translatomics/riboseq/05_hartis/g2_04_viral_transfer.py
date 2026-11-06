"""G2-04 宿主逐位点模型向 PEDV 的迁移诊断（只用非候选区段；不作为纯技术泄漏的假定，只作外部诊断）。
病毒密码子：父 ORF 读框；j = (p − ORF_start)/3；29 nt P-site（+ 链，offset 12）；813 由 5′ 矩阵重建。
病毒非候选区段：去掉所有 +1 家族 [最 5′ 起点 − 30, 终止 + 30]、全部候选起点 ±30、排除区间（12,586–12,622；canonical 起始 −12..+30 即 j < 10；终止前 30 nt；26,702 ± 15；ORF 重叠区）。
诊断：
 (1) 组成：端点碱基组合 (b5F1,b5F2,b3F1,b3F2) 的分布（宿主训练密码子 vs 病毒非候选密码子），Jensen–Shannon 距离；GC。
 (2) 读框偏差：F0 份额、F1/(F1+F2)（观测）vs 宿主模型预测。
 (3) 覆盖结构：窗口 m、n0 分布；前 1 个密码子占窗口 m 的比例（尖峰度）；是否在宿主支持范围内。
 (4) 迁移校准：在宿主支持范围内的病毒非候选窗口上，用宿主同类型、同 m 分层的 z 经验零分布计算 p；超出支持范围 → out_of_host_support，不计算 p。
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
from scipy.stats import kstest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import g2_hostcodon as G
T = G.T; V3 = G.V3; R = os.path.join(V3, '..'); RG = os.path.join(R, 'riboseq_GIb_mapping_psite'); GL = 27959; SHR = 20.0
ORFS = [('ORF1a', 293, 12622), ('ORF1b', 12592, 20613), ('S', 20610, 24758), ('ORF3', 24758, 25033), ('E', 25361, 25600), ('M', 25608, 26288), ('N', 26300, 27625)]
SEQ = ''.join(l.strip() for l in open(os.path.join(R, 'refbuild_GIb', 'v1.1', 'virus', 'PEDV_GIb.fa')) if not l.startswith('>')).upper()
LUT = {c: i for i, c in enumerate('ACGT')}
sa = np.array([LUT.get(c, 4) for c in SEQ], np.int8)        # 0-based：位置 p（1-based）→ sa[p−1]
PS = pd.read_csv(os.path.join(RG, 'psite_counts_PEDV_GIb.tsv.gz'), sep='\t', comment='#', usecols=['sample_id', 'strand', 'position', 'read_length', 'psite_count'])
PS = PS[(PS.strand == '+') & (PS.read_length == 29)]
cover = np.zeros(GL + 2, int)
for _, a, b in ORFS: cover[a:b + 1] += 1
V = pd.read_csv(os.path.join(R, 'HARTIS_design', 'tables', 'viral_plus1_start_architecture.tsv'), sep='\t')
cand = np.zeros(GL + 2, bool)
for f, d in V.groupby('family_id'): cand[max(1, d.start_position.min() - 30):min(GL, d.stop_position.iloc[0] + 30) + 1] = True
for s in V.start_position.values: cand[max(1, s - 30):s + 31] = True
excl = np.zeros(GL + 2, bool); excl[12586:12623] = True; excl[26702 - 15:26702 + 16] = True; excl[cover != 1] = True


def viral_array(lib):
    s = G.LIBS[lib]; a = np.zeros(GL + 8)
    if s == SR.lib("RIBO_VERO_CHX_1MOI_16h"):
        m = pd.read_csv(os.path.join(RG, '_work', 'server', 'stats', s, 'pedv_5p_matrix.tsv.gz'), sep='\t')
        m = m[(m.strand == '+') & (m.L == 29)]; g = m.groupby(m.pos5 + 12).n.sum(); g = g[g.index <= GL]
    else:
        g = PS[PS.sample_id == s].groupby('position').psite_count.sum()
    a[g.index.values] = g.values
    return a


def viral_codons(lib):
    a = viral_array(lib); rows = []
    for n, o0, o1 in ORFS:
        N = (o1 - o0 + 1) // 3
        js = np.arange(10, N - 11)
        p = o0 + 3 * js                                       # F0 P-site（1-based）
        ok = ~(cand[p] | cand[p + 1] | cand[p + 2] | excl[p] | excl[p + 1] | excl[p + 2] | excl[p + 3])
        js, p = js[ok], p[ok]
        b = lambda off: sa[p + off - 1]
        rows.append(pd.DataFrame(dict(ORF=n, j=js, pos=p, n0=a[p], n1=a[p + 1], n2=a[p + 2], n0n=a[p + 3],
                                      b5F1=b(1 - 12), b5F2=b(2 - 12), b3F1=b(1 + 16), b3F2=b(2 + 16))))
    D = pd.concat(rows, ignore_index=True); D['m'] = D.n1 + D.n2
    return D


def seg_windows(D, wt):
    """病毒窗口只在连续的非候选密码子段内平铺。"""
    D = D.sort_values(['ORF', 'j']).copy()
    D['seg'] = ((D.j.diff() != 1) | (D.ORF != D.ORF.shift())).cumsum()
    if wt == 'single': D['key'] = np.arange(len(D))
    else:
        w = 5 if wt == 'cod5' else 50
        D['rank'] = D.groupby('seg').cumcount(); D['key'] = D.seg * 100000 + D['rank'] // w
    D = D.assign(E1=D.m * D.pj, V=D.m * D.pj * (1 - D.pj))
    W = D.groupby('key').agg(ORF=('ORF', 'first'), pos_lo=('pos', 'min'), n0=('n0', 'sum'), n1=('n1', 'sum'), n2=('n2', 'sum'), m=('m', 'sum'),
                             E1=('E1', 'sum'), V=('V', 'sum'), nc=('j', 'size'), top1=('m', 'max'))
    if wt == 'cod5': W = W[W.nc == 5]
    if wt == 'cod50': W = W[W.nc == 50]
    W = W[W.m > 0].copy(); W['z'] = (W.n1 - W.E1) / np.sqrt(np.maximum(W.V, 1e-9)); W['ph'] = W.E1 / W.m
    W['top1_share'] = W.top1 / W.m
    return W


comp, cal, orfs, supp = [], [], [], []
for lib in G.LIBS:
    H = G.codon_table(lib); H = H[(H.b5F1 < 4) & (H.b5F2 < 4) & (H.b3F1 < 4) & (H.b3F2 < 4)]
    tr = H[~H.heldout]; trm = tr[tr.m > 0]; p0 = trm.n1.sum() / trm.m.sum()
    g = trm.groupby(['b5F1', 'b5F2', 'b3F1', 'b3F2'])[['n1', 'm']].sum(); g['pj'] = (g.n1 + SHR * p0) / (g.m + SHR)
    Vd = viral_codons(lib); Vd = Vd[(Vd.b5F1 < 4) & (Vd.b5F2 < 4) & (Vd.b3F1 < 4) & (Vd.b3F2 < 4)]
    Vd = Vd.join(g.pj, on=['b5F1', 'b5F2', 'b3F1', 'b3F2']); Vd['pj'] = Vd.pj.fillna(p0)
    # (1) 组成（按 m 加权的端点碱基组合分布）
    hk = trm.groupby(['b5F1', 'b5F2', 'b3F1', 'b3F2']).m.sum(); vk = Vd[Vd.m > 0].groupby(['b5F1', 'b5F2', 'b3F1', 'b3F2']).m.sum()
    idx = hk.index.union(vk.index); P_ = hk.reindex(idx, fill_value=0).values + 1e-9; Q_ = vk.reindex(idx, fill_value=0).values + 1e-9
    P_, Q_ = P_ / P_.sum(), Q_ / Q_.sum(); M_ = (P_ + Q_) / 2
    js = np.sqrt(0.5 * (P_ * np.log2(P_ / M_)).sum() + 0.5 * (Q_ * np.log2(Q_ / M_)).sum())
    hb = trm[['b5F1', 'b5F2']].values; vb = Vd[Vd.m > 0][['b5F1', 'b5F2']].values
    comp.append(dict(library=lib, JS_distance_endbase256_mweighted=js, host_codons_m_gt0=len(trm), viral_noncand_codons_m_gt0=int((Vd.m > 0).sum()),
                     host_GC_at_5p_F1=np.isin(hb[:, 0], [1, 2]).mean(), viral_GC_at_5p_F1=np.isin(vb[:, 0], [1, 2]).mean() if len(vb) else np.nan,
                     host_F0_share=tr[['n0']].sum().iloc[0] / tr[['n0', 'n1', 'n2']].sum().sum(),
                     viral_noncand_F0_share=Vd.n0.sum() / max(Vd[['n0', 'n1', 'n2']].sum().sum(), 1),
                     host_F1share=p0, viral_noncand_F1share_obs=Vd.n1.sum() / max(Vd.m.sum(), 1),
                     viral_noncand_F1share_pred=(Vd.m * Vd.pj).sum() / max(Vd.m.sum(), 1)))
    for orf, d in Vd.groupby('ORF'):
        orfs.append(dict(library=lib, ORF=orf, codons=len(d), m=int(d.m.sum()), n0=int(d.n0.sum()), F0_share=d.n0.sum() / max(d[['n0', 'n1', 'n2']].sum().sum(), 1),
                         F1share_obs=d.n1.sum() / max(d.m.sum(), 1), F1share_pred_host_model=(d.m * d.pj).sum() / max(d.m.sum(), 1)))
    # (4) 迁移校准
    for wt in ['single', 'cod5', 'cod50']:
        # 宿主训练窗口（与 g2_03 同口径）
        Hd = tr.join(g.pj, on=['b5F1', 'b5F2', 'b3F1', 'b3F2']); Hd['pj'] = Hd.pj.fillna(p0)
        Hd = Hd.sort_values(['tx_idx', 'j'])
        if wt == 'single': key = np.arange(len(Hd))
        else:
            w = 5 if wt == 'cod5' else 50; key = Hd.tx_idx.astype(np.int64) * 100000 + (Hd.j - 10) // w
        Hd = Hd.assign(E1=Hd.m * Hd.pj, V=Hd.m * Hd.pj * (1 - Hd.pj), key=key)
        HW = Hd.groupby('key').agg(n0=('n0', 'sum'), n1=('n1', 'sum'), m=('m', 'sum'), E1=('E1', 'sum'), V=('V', 'sum'), nc=('j', 'size'), top1=('m', 'max'))
        if wt != 'single': HW = HW[HW.nc == (5 if wt == 'cod5' else 50)]
        HW = HW[HW.m > 0]; HW['z'] = (HW.n1 - HW.E1) / np.sqrt(np.maximum(HW.V, 1e-9))
        edges = np.unique(np.r_[np.quantile(HW.m, np.linspace(0, 1, 11))[:-1], np.quantile(HW.m, .99), np.inf])
        HW['mb'] = np.searchsorted(edges, HW.m, side='right') - 1
        zz = {b: np.sort(d.z.values) for b, d in HW.groupby('mb')}
        m_hi = np.quantile(HW.m, .995); n0_hi = np.quantile(HW.n0, .995)
        VW = seg_windows(Vd, wt)
        VW['in_host_support'] = (VW.m <= m_hi) & (VW.n0 <= n0_hi)
        VW['mb'] = np.searchsorted(edges, VW.m, side='right') - 1
        pe = []
        for r in VW.itertuples():
            if not r.in_host_support: pe.append(np.nan); continue
            z_ = zz.get(r.mb); K = len(z_)
            pe.append((1 + K - np.searchsorted(z_, r.z, side='left')) / (1 + K) if K >= 100 else np.nan)
        VW['p_emp'] = pe; VW['library'] = lib; VW['wtype'] = wt
        supp.append(VW)
        x = VW.p_emp.dropna()
        cal.append(dict(library=lib, wtype=wt, viral_windows=len(VW), in_host_support=int(VW.in_host_support.sum()),
                        out_of_host_support=int((~VW.in_host_support).sum()), host_m_p995=m_hi, host_n0_p995=n0_hi,
                        viral_m_median=VW.m.median(), viral_m_max=VW.m.max(), n_with_p=len(x),
                        frac_p_le_0_05=(x <= .05).mean() if len(x) else np.nan, frac_p_le_0_01=(x <= .01).mean() if len(x) else np.nan,
                        frac_p_ge_0_95=(x >= .95).mean() if len(x) else np.nan,
                        ks_p=kstest(x, 'uniform').pvalue if len(x) > 20 else np.nan, mean_z=VW.loc[VW.in_host_support, 'z'].mean(),
                        host_top1_share_median=(HW.top1 / HW.m).median(), viral_top1_share_median=VW.top1_share.median()))
    print(lib, 'done', flush=True)
C = pd.DataFrame(comp); C.to_csv(os.path.join(T, 'g2_transfer_composition.tsv'), sep='\t', index=False, float_format='%.4g')
O = pd.DataFrame(orfs); O.to_csv(os.path.join(T, 'g2_transfer_by_ORF_noncandidate.tsv'), sep='\t', index=False, float_format='%.4g')
CA = pd.DataFrame(cal); CA.to_csv(os.path.join(T, 'g2_transfer_calibration.tsv'), sep='\t', index=False, float_format='%.4g')
pd.concat(supp).to_csv(os.path.join(T, 'g2_transfer_viral_windows.tsv.gz'), sep='\t', index=False, float_format='%.5g')
pd.set_option('display.width', 260)
print(C.round(3).to_string()); print(O.pivot_table(index='ORF', columns='library', values=['F1share_obs', 'F1share_pred_host_model']).round(3).to_string())
print(CA.round(3).to_string())
