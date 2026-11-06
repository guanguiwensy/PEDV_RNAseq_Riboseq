# -*- coding: utf-8 -*-
"""G2d：宿主背景（EMP_Zs 的组成部分）向 PEDV 迁移的可行性——只做描述与诊断，不计算任何病毒 p 值，不碰候选。
病毒单位：7 个 canonical ORF 的父读框；只用**非候选**密码子（遮蔽：69 个 +1 家族 [最 5′ 起点 − 30, 终止 + 30]、全部 567 个候选起点 ±30、
12,586–12,622、26,702 ± 15、ORF 重叠区、canonical 起始后 30 nt 与终止前 30 nt）。PEDV 非候选区**不假定为纯阴性**。
病毒计数：+ 链 P-site；060/061/811/062/065/064 来自 psite_counts_PEDV_GIb（RGP_v1.3，含剪接 read，无法分离）；
          813 由 pedv_5p_matrix（只含 5′ 端 22 nt 内无剪接的 read）按 offset 12 重建。
宿主参照：同库×读长、主口径、tier_main，按基因 5 折的样本外上下文期望；φ̂_host(m) 由宿主窗口估计。
诊断（逐库×读长 × ORF）：
  Q4 周期性：非候选密码子 F0/F1/F2 份额（与宿主同口径 F0 比较）
  Q1 读框偏好：观测 F1/(F1+F2) vs 宿主上下文期望（O − E）；上下文效应可比性 = 16 组 5′ 端二碱基上宿主与病毒 F1 份额的加权相关
  Q2 离散：病毒窗口 z（宿主上下文期望）方差 vs 宿主 φ̂(m)，在同一 m 档内比较（比值 = 病毒 Var z / 宿主 Var z）
  Q3 深度/支持：病毒窗口 m 分布；按 EMP_Zs 单元级规则（匹配到宿主窗口 K ≥ 100 且半径 ≤ ±0.75、宿主 φ̂ 不外推、有效基因数 ≥ 30、m < 300、
     深度档已认证）判定 in_support / out_of_host_support 的比例
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
import g2b_common as G
from g2b_calibrate import match_diag, SHR, TRAIN_MMIN
OUT = os.path.join(G.OUT, 'G2d_round5'); os.makedirs(OUT, exist_ok=True)
R = os.path.abspath(os.path.join(G.V3, '..')); RG = os.path.join(R, 'riboseq_GIb_mapping_psite'); GL = 27959
ORFS = [('ORF1a', 293, 12622), ('ORF1b', 12592, 20613), ('S', 20610, 24758), ('ORF3', 24758, 25033), ('E', 25361, 25600), ('M', 25608, 26288), ('N', 26300, 27625)]
SEQ = ''.join(l.strip() for l in open(os.path.join(R, 'refbuild_GIb', 'v1.1', 'virus', 'PEDV_GIb.fa')) if not l.startswith('>')).upper()
sa = np.array([{'A': 0, 'C': 1, 'G': 2, 'T': 3}.get(c, 4) for c in SEQ] + [4] * 80, np.int16)   # 位置 p（1-based）→ sa[p−1]
cover = np.zeros(GL + 80, int)
for _, a, b in ORFS: cover[a:b + 1] += 1
V = pd.read_csv(os.path.join(R, 'HARTIS_design', 'tables', 'viral_plus1_start_architecture.tsv'), sep='\t', usecols=['family_id', 'start_position', 'stop_position'])
mask = np.zeros(GL + 80, bool)
for f, d in V.groupby('family_id'): mask[max(1, d.start_position.min() - 30):min(GL, d.stop_position.max() + 30) + 1] = True
for s in V.start_position.values: mask[max(1, s - 30):s + 31] = True
mask[12586:12623] = True; mask[26702 - 15:26702 + 16] = True; mask[cover != 1] = True
SUP = pd.read_csv(os.path.join(G.OUT, 'G2c_freeze_checks', 'EMP_Zs_host_support_map.tsv'), sep='\t')
SUP = SUP[SUP.scope == 'in_support'].set_index(['sid', 'L', 'w']).certified_m_range
WINS = (1, 5, 50); MB = [5, 10, 30, 100, 300, 1000, 1e9]
_ps = {}


def viral_counts(sid, L):
    a = np.zeros(GL + 80)
    if sid == SR.lib("RIBO_VERO_CHX_1MOI_16h"):
        m = pd.read_csv(os.path.join(RG, '_work', 'server', 'stats', sid, 'pedv_5p_matrix.tsv.gz'), sep='\t')
        m = m[(m.strand == '+') & (m.L == L)]; g = m.groupby(m.pos5 + G.SPEC[sid][L]).n.sum()
    else:
        if 'P' not in _ps:
            _ps['P'] = pd.read_csv(os.path.join(RG, 'psite_counts_PEDV_GIb.tsv.gz'), sep='\t', comment='#',
                                   usecols=['sample_id', 'strand', 'position', 'read_length', 'offset', 'psite_count'])
        P = _ps['P']; P = P[(P.sample_id == sid) & (P.strand == '+') & (P.read_length == L)]
        assert (P.offset.astype(int) == G.SPEC[sid][L]).all(), (sid, L, P.offset.unique())
        g = P.groupby('position').psite_count.sum()
    g = g[(g.index >= 1) & (g.index <= GL)]; a[g.index.values] = g.values
    return a


def viral_codons(sid, L):
    a = viral_counts(sid, L); off = G.SPEC[sid][L]; rows = []
    for nm, o0, o1 in ORFS:
        N = (o1 - o0 + 1) // 3; js = np.arange(10, N - 10); p = o0 + 3 * js
        ok = ~(mask[p] | mask[p + 1] | mask[p + 2] | mask[p + 3])
        js, p = js[ok], p[ok]
        t = p + 1 - off                                                   # F1 读段 5′ 端（1-based）
        b = np.stack([sa[t - 1], sa[t], sa[t + L - 2], sa[t + L - 1]])
        ctx = np.where((b == 4).any(0), 256, b[0] * 64 + b[1] * 16 + b[2] * 4 + b[3])
        rows.append(pd.DataFrame(dict(ORF=nm, j=js, pos=p, n0=a[p], n1=a[p + 1], n2=a[p + 2], ctx=ctx)))
    D = pd.concat(rows, ignore_index=True); D['m'] = D.n1 + D.n2
    return D


def viral_windows(D, w):
    D = D.sort_values(['ORF', 'j']).copy()
    D['seg'] = ((D.j.diff() != 1) | (D.ORF != D.ORF.shift())).cumsum()
    D['k'] = D.groupby('seg').cumcount() // w
    W = D.groupby(['ORF', 'seg', 'k']).agg(n0=('n0', 'sum'), n1=('n1', 'sum'), n2=('n2', 'sum'), m=('m', 'sum'), E1=('E1', 'sum'),
                                          V=('V', 'sum'), nc=('j', 'size'), start=('pos', 'min')).reset_index()
    W = W[W.nc == w]; W['dens'] = W.n0 / (3 * w)
    return W


def host_reference(sid, L):
    C = G.codon_table(sid, L)
    tr = C[C.m > 0]; p0 = tr.n1.sum() / tr.m.sum()
    g = tr.groupby('ctx')[['n1', 'm']].sum(); pc_all = (g.n1 + SHR * p0) / (g.m + SHR)
    pj = np.zeros(len(C))
    for f in range(G.NFOLD):                                         # 宿主窗口用样本外期望
        s = C[(C.fold != f) & (C.m > 0)]; q0 = s.n1.sum() / s.m.sum(); h = s.groupby('ctx')[['n1', 'm']].sum()
        qc = (h.n1 + SHR * q0) / (h.m + SHR); sel = (C.fold == f).values; pj[sel] = qc.reindex(C.ctx[sel]).fillna(q0).values
    return C, pj, pc_all, p0


def main():
    rows_q4, rows_q1, rows_ctx, rows_q2, rows_q3 = [], [], [], [], []
    for sid, L in sorted(G.QUALIFIED):
        lab = G.LABEL[sid]; meta = dict(sid=sid, label=lab, drug='HAR' if 'HAR' in lab else 'CHX', L=L)
        C, pj, pc_all, p0 = host_reference(sid, L)
        D = viral_codons(sid, L)
        D['mu'] = pd.Series(D.ctx).map(pc_all).fillna(p0).values; D['E1'] = D.m * D.mu; D['V'] = D.m * D.mu * (1 - D.mu)
        hF = C[['n0', 'n1', 'n2']].sum(); hF0 = hF.n0 / hF.sum(); hshare = hF.n1 / (hF.n1 + hF.n2)
        for orf, d in D.groupby('ORF', sort=False):
            t = d[['n0', 'n1', 'n2']].sum(); tot = t.sum()
            rows_q4.append(dict(**meta, ORF=orf, noncand_codons=len(d), reads=int(tot), F0=t.n0 / tot if tot else np.nan, F1=t.n1 / tot if tot else np.nan,
                                F2=t.n2 / tot if tot else np.nan, host_F0=hF0, F0_minus_host=(t.n0 / tot - hF0) if tot else np.nan))
            mm = d.m.sum()
            rows_q1.append(dict(**meta, ORF=orf, m=int(mm), F1share_obs=d.n1.sum() / mm if mm else np.nan, F1share_hostctx_expected=d.E1.sum() / mm if mm else np.nan,
                                O_minus_E=(d.n1.sum() - d.E1.sum()) / mm if mm else np.nan, host_pooled_F1share=hshare))
        # 上下文效应可比性：16 组 5′ 端二碱基（b5F1, b5F2）
        hc = C.assign(g5=C.ctx // 16).groupby('g5')[['n1', 'm']].sum(); vc = D[D.ctx < 256].assign(g5=lambda x: x.ctx // 16).groupby('g5')[['n1', 'm']].sum()
        J = hc.join(vc, lsuffix='_h', rsuffix='_v').dropna(); J = J[J.m_v >= 50]
        if len(J) >= 5:
            xh = J.n1_h / J.m_h; xv = J.n1_v / J.m_v; wv = J.m_v
            mh, mv = np.average(xh, weights=wv), np.average(xv, weights=wv)
            r = np.sum(wv * (xh - mh) * (xv - mv)) / np.sqrt(np.sum(wv * (xh - mh) ** 2) * np.sum(wv * (xv - mv) ** 2))
            rows_ctx.append(dict(**meta, n_groups=len(J), weighted_corr_host_vs_viral=r, host_range=f'{xh.min():.2f}-{xh.max():.2f}',
                                 viral_range=f'{xv.min():.2f}-{xv.max():.2f}', viral_minus_host_mean=mv - mh))
        for w in WINS:
            HW = G.windows(C, w, pj); HW = HW[HW.m >= TRAIN_MMIN].copy()
            HW['z'] = (HW.n1 - HW.E1) / np.sqrt(HW.V.clip(lower=1e-9)); HW['lm'] = np.log2(HW.m); HW['ld'] = np.log2(HW.dens + 1 / (3 * w))
            lbv = np.floor(HW.lm / 0.25).astype(int); vc_ = lbv.value_counts(); okb = vc_[vc_ >= 30].index.values
            xmin, xmax = (okb.min() + .5) * .25, (okb.max() + .5) * .25
            VW = viral_windows(D, w); VW = VW[VW.m >= 1].copy()
            VW['z'] = (VW.n1 - VW.E1) / np.sqrt(VW.V.clip(lower=1e-9)); VW['lm'] = np.log2(VW.m); VW['ld'] = np.log2(VW.dens + 1 / (3 * w))
            # Q2：同 m 档的 z 方差
            for lo, hi in zip(MB[:-1], MB[1:]):
                h = HW[(HW.m >= lo) & (HW.m < hi)]; v = VW[(VW.m >= lo) & (VW.m < hi)]
                if len(v) < 5: continue
                rows_q2.append(dict(**meta, w=w, m_bin=f'{lo}-{int(hi) - 1 if hi < 1e9 else "inf"}', n_host=len(h), n_viral=len(v),
                                    host_var_z=h.z.var() if len(h) > 1 else np.nan, viral_var_z=v.z.var(), var_ratio_viral_host=v.z.var() / h.z.var() if len(h) > 30 else np.nan,
                                    host_mean_z=h.z.mean() if len(h) else np.nan, viral_mean_z=v.z.mean(),
                                    viral_ORFs=';'.join(f'{o}:{c}' for o, c in v.ORF.value_counts().items())))
            # Q3：支持判定（与 EMP_Zs 冻结规则相同；不计算 p）
            E = VW[VW.m >= 10].copy()
            if len(E):
                dg = match_diag(HW.lm.values, HW.ld.values, HW.tx_idx.values, E.lm.values, E.ld.values, (2, 3, 5))
                cert = SUP.get((sid, L, w), 'none'); mmax = {'none': 0, '10-29': 30, '10-99': 100, '10-299': 300}[cert]
                E['K'] = dg['K']; E['radius'] = dg['radius']; E['eff'] = dg['eff_genes']
                E['phi_ext'] = (E.lm < xmin) | (E.lm > xmax)
                E['reason'] = np.select([E.m >= 300, E.m >= mmax, E.phi_ext, E.K.isna() | (E.radius > 3), E.eff < 30],
                                        ['m>=300', 'depth_bin_not_certified', 'phi_extrapolated', 'K<100_or_radius>0.75', 'eff_genes<30'], 'in_support')
                for orf, d in E.groupby('ORF', sort=False):
                    rec = dict(**meta, w=w, ORF=orf, viral_windows_m10=len(d), m_median=d.m.median(), m_max=d.m.max(), certified_m_range=cert,
                               in_support=(d.reason == 'in_support').mean())
                    for rs in ('m>=300', 'depth_bin_not_certified', 'phi_extrapolated', 'K<100_or_radius>0.75', 'eff_genes<30'): rec[rs] = (d.reason == rs).mean()
                    ins = d[d.reason == 'in_support']
                    rec['in_support_mean_z'] = ins.z.mean() if len(ins) else np.nan
                    rows_q3.append(rec)
        print(sid, L, 'done', flush=True)

    for nm, rows in (('Q4_periodicity', rows_q4), ('Q1_frame_preference', rows_q1), ('Q1_context_comparability', rows_ctx),
                     ('Q2_dispersion', rows_q2), ('Q3_support', rows_q3)):
        pd.DataFrame(rows).to_csv(os.path.join(OUT, f'transfer_{nm}.tsv'), sep='\t', index=False, float_format='%.5g')
    print('all done')


if __name__ == '__main__':
    main()
