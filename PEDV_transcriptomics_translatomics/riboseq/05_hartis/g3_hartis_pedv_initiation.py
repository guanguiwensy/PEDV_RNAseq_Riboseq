# -*- coding: utf-8 -*-
"""HARTIS_v3 · PEDV +1 / +2 框候选起点的起始相关占据（initiation-associated occupancy）。
规则全部按 PREREGISTRATION_HARTIS_v3_rev3_G2_FINAL.md 冻结版执行，本脚本不调任何参数：
  候选：每个父 canonical ORF 内，候选框（+1 或 +2，相对父读框）上的 AUG / near-cognate{CUG,GUG,UUG,AUA,AUU,AUC,ACG}，
        延伸到该框第一个终止，ORF > 40 aa（与 plus1_orf_capacity/s15 相同规则；+1 框应复现 567 起点 / 69 家族）；
        家族 = 同父 ORF、同候选框、同终止。
  尺度：W1 = 起始密码子所在密码子单元 j（候选框计数为 n1(+1) 或 n2(+2)）；W5 = j..j+4（候选框位点 s, s+3, …, s+12）。
  统计：EMP_Zs（宿主端点上下文期望、φ̂(m) 标准化、按 log2 m 与 log2 F0 密度局部匹配，K ≥ 100），各库×读长独立；
        +1 框取 z 的上尾，+2 框取下尾（= F2 富集，宿主对调检验已认证）；Z_cand 为符号调整后的 z_s（正 = 候选框富集）。
  支持：K ≥ 100 且半径 ≤ 3 格、φ̂ 不外推、有效基因数 ≥ 30、m < 300、深度档已认证 → in_support（给 p_emp、Z_cand）；
        m ≥ 认证上限 → capped_reference（探索性 Z_cap、p_cap）；其余 → out_of_host_support；m < 10 → not_assessable_low_m；
        窗口触及排除区（12,586–12,622；26,702 ± 15；ORF 重叠区；父 ORF 起始后 30 nt；终止前 30 nt）→ not_assessable_excluded。
  多读长：各读长独立；W5 另给 J_POOL（读长计数合并，各读长自己的上下文期望；J_POOL 宿主认证区间单独）。
  迁移标记（不阻断计算，只标记；迁移闸标准未定）：periodicity_ok = 父 ORF 非候选区 F0 ≥ 宿主 F0 − 0.10；context_ok = 该库×读长端点上下文相关 ≥ 0.5。
  证据类别：initiation_associated_occupancy。**不使用、不报告任何下游延伸证据。**
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
import os, sys, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); V3 = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(V3, 'G2b_host_background', 'scripts'))
import g2b_common as G
import g2d_transfer_feasibility as TF
from g2b_calibrate import emp_p, phi_fit, match_diag
OUT = os.path.abspath(os.path.join(HERE, '..')); TB = os.path.join(OUT, 'tables'); os.makedirs(TB, exist_ok=True)
GL = TF.GL; SEQ = TF.SEQ; ORFS = TF.ORFS
START = {'ATG', 'CTG', 'GTG', 'TTG', 'ATA', 'ATT', 'ATC', 'ACG'}; STOP = {'TAA', 'TAG', 'TGA'}
cod = lambda p: SEQ[p - 1:p + 2] if 1 <= p and p + 2 <= GL else 'NA'
G2B = os.path.join(V3, 'G2b_host_background')
SUP1 = pd.read_csv(os.path.join(G2B, 'G2c_freeze_checks', 'EMP_Zs_host_support_map.tsv'), sep='\t')
SUP1 = SUP1[SUP1.scope == 'in_support'].set_index(['sid', 'L', 'w']).certified_m_range
SUPJ = pd.read_csv(os.path.join(G2B, 'G2d_round5', 'JPOOL_support_map.tsv'), sep='\t').set_index(['sid', 'w']).certified_pooled_m_range
CAP = {'none': 0, '10-29': 30, '10-99': 100, '10-299': 300}; REFBIN = {'10-29': (10, 30), '10-99': (30, 100), '10-299': (100, 300)}
CTXC = pd.read_csv(os.path.join(G2B, 'G2d_round5', 'transfer_Q1_context_comparability.tsv'), sep='\t').set_index(['sid', 'L']).weighted_corr_host_vs_viral
PER = pd.read_csv(os.path.join(G2B, 'G2d_round5', 'transfer_Q4_periodicity.tsv'), sep='\t').set_index(['sid', 'L', 'ORF'])

# ---------------- 1. 候选枚举 ----------------
frames = {}
for nm, s, e in ORFS:
    gfp = (s - 1) % 3; frames[nm] = (gfp, s, e)
rows = []
for nm, (gfp, ps, pe) in frames.items():
    for rel in (1, 2):
        c = (gfp + rel) % 3
        for p in range(ps, pe + 1):
            if (p - 1) % 3 != c or cod(p) not in START: continue
            q = p
            while q + 2 <= GL and cod(q) not in STOP: q += 3
            if q + 2 > GL: continue
            aa = (q - p) // 3
            if aa <= 40: continue
            same = [o for o, a, b in ORFS if o != nm and (a - 1) % 3 == c and a <= q + 2 and b >= p]
            rows.append(dict(frame=f'+{rel}', parent=nm, start_position=p, start_codon=cod(p).replace('T', 'U'), stop_position=q,
                             ORF_length_aa=aa, family_id=f'P{rel}_{nm}_stop{q}',
                             candidate_frame_equals_canonical=';'.join(same) or ''))
CAND = pd.DataFrame(rows).sort_values(['frame', 'parent', 'start_position']).reset_index(drop=True)
A = pd.read_csv(os.path.join(V3, '..', 'HARTIS_design', 'tables', 'viral_plus1_start_architecture.tsv'), sep='\t')
p1 = CAND[CAND.frame == '+1']
chk = dict(plus1_starts_new=len(p1), plus1_starts_ref=len(A), plus1_families_new=p1.family_id.nunique(), plus1_families_ref=A.family_id.nunique(),
           plus1_set_identical=set(zip(p1.parent, p1.start_position, p1.stop_position)) == set(zip(A.parent_canonical_ORF, A.start_position, A.stop_position)))
print('候选枚举复核', chk, '| +2 起点', (CAND.frame == '+2').sum(), '家族', CAND[CAND.frame == '+2'].family_id.nunique(), flush=True)
CAND.to_csv(os.path.join(TB, 'candidates_plus1_plus2.tsv'), sep='\t', index=False)
pd.Series(chk).to_csv(os.path.join(TB, 'candidate_enumeration_check.tsv'), sep='\t', header=['value'])

# 排除区（与宿主训练同口径的起止近端 + 设计排除区）
cover = np.zeros(GL + 80, int)
for _, a, b in ORFS: cover[a:b + 1] += 1
excl = np.zeros(GL + 80, bool); excl[12586:12623] = True; excl[26702 - 15:26702 + 16] = True; excl[cover != 1] = True
CAND['o0'] = CAND.parent.map(lambda n: frames[n][1]); CAND['N'] = CAND.parent.map(lambda n: (frames[n][2] - frames[n][1] + 1) // 3)
CAND['j'] = (CAND.start_position - CAND.o0 - CAND.frame.str[1].astype(int)) // 3


def viral_codon_table(sid, L):
    """父 ORF 全部密码子（不遮蔽候选），含端点上下文；与 TF.viral_codons 相同的计数与上下文定义。"""
    a = TF.viral_counts(sid, L); off = G.SPEC[sid][L]; out = {}
    for nm, o0, o1 in ORFS:
        N = (o1 - o0 + 1) // 3; js = np.arange(0, N); p = o0 + 3 * js; t = p + 1 - off
        b = np.stack([TF.sa[t - 1], TF.sa[t], TF.sa[t + L - 2], TF.sa[t + L - 1]])
        ctx = np.where((b == 4).any(0), 256, b[0] * 64 + b[1] * 16 + b[2] * 4 + b[3])
        out[nm] = pd.DataFrame(dict(j=js, pos=p, n0=a[p], n1=a[p + 1], n2=a[p + 2], ctx=ctx))
    return out


def host_null(sid, L, w):
    """宿主训练窗口（全部基因；样本外上下文期望）→ φ̂、匹配所需数组、capped 参照。"""
    C, pj, pc_all, p0 = TF.host_reference(sid, L)
    W = G.windows(C, w, pj); W = W[W.m >= 5].copy()
    W['z'] = (W.n1 - W.E1) / np.sqrt(W.V.clip(lower=1e-9)); W['lm'] = np.log2(W.m); W['ld'] = np.log2(W.dens + 1 / (3 * w))
    ph = phi_fit(W.lm.values, W.z.values)
    lbv = np.floor(W.lm / 0.25).astype(int); vc = lbv.value_counts(); okb = vc[vc >= 30].index.values
    rng_phi = ((okb.min() + .5) * .25, (okb.max() + .5) * .25)
    W['zs'] = W.z / np.sqrt(ph(W.lm.values)); W['d'] = (W.n1 - W.E1) / W.m
    return dict(W=W, ph=ph, rng=rng_phi, pc=pc_all, p0=p0)


def window_stats(VC, cand, w, mu_fn):
    """对每个候选：窗口内 n0/n1/n2、E1、V（候选框计数在 n1 或 n2）。"""
    rec = []
    for r in cand.itertuples():
        d = VC[r.parent]; jj = np.arange(r.j, r.j + w)
        ok_range = (r.j >= 10) and (r.j + w - 1 <= r.N - 11)
        pos = r.o0 + 3 * jj
        exc = (not ok_range) or bool(excl[np.concatenate([pos, pos + 1, pos + 2])].any())
        x = d.iloc[jj[jj < len(d)]]
        mu = mu_fn(x.ctx.values); m_j = (x.n1 + x.n2).values
        rec.append(dict(n0=x.n0.sum(), n1=x.n1.sum(), n2=x.n2.sum(), E1=(m_j * mu).sum(), V=(m_j * mu * (1 - mu)).sum(), excluded=exc))
    return pd.DataFrame(rec, index=cand.index)


def evaluate(S, cand, null, cert, w, method):
    """S：窗口统计；返回逐候选的状态与统计量（+2 框用下尾）。"""
    W = null['W']; S = S.copy(); S['m'] = S.n1 + S.n2
    S['sign'] = np.where(cand.frame == '+1', 1.0, -1.0)
    S['n_cand'] = np.where(cand.frame == '+1', S.n1, S.n2)
    S['z'] = (S.n1 - S.E1) / np.sqrt(S.V.clip(lower=1e-9)); S['lm'] = np.log2(S.m.clip(lower=1)); S['ld'] = np.log2(S.n0 / (3 * w) + 1 / (3 * w))
    S['phi_extrapolated'] = (S.lm < null['rng'][0]) | (S.lm > null['rng'][1])
    S['Z_cand'] = S.sign * S.z / np.sqrt(null['ph'](S.lm.values))
    cap = CAP[cert]; S['certified_m_range'] = cert
    S['status'] = np.select([S.excluded, S.m < 10, (cap > 0) & (S.m >= cap), cap == 0],
                            ['not_assessable_excluded', 'not_assessable_low_m', 'capped_reference', 'out_of_host_support(no_certified_bin)'], 'pending')
    pend = S.status == 'pending'
    for c in ('K', 'eff_genes', 'max_gene_frac', 'radius', 'p_emp'): S[c] = np.nan
    if pend.any():
        e = S[pend]; zs = e.z.values / np.sqrt(null['ph'](e.lm.values))
        dg = match_diag(W.lm.values, W.ld.values, W.tx_idx.values, e.lm.values, e.ld.values, (2, 3, 5))
        pu, _, _ = emp_p(W.zs.values, W.lm.values, W.ld.values, zs, e.lm.values, e.ld.values, True)
        pl, _, _ = emp_p(W.zs.values, W.lm.values, W.ld.values, zs, e.lm.values, e.ld.values, False)
        S.loc[pend, 'K'] = dg['K']; S.loc[pend, 'eff_genes'] = dg['eff_genes']; S.loc[pend, 'max_gene_frac'] = dg['max_gene_frac']; S.loc[pend, 'radius'] = dg['radius']
        S.loc[pend, 'p_emp'] = np.where(e.sign.values > 0, pu, pl)
        e = S[pend]
        ins = (e.K >= 100) & (e.radius <= 3) & (e.eff_genes >= 30) & ~e.phi_extrapolated & (e.m < 300)
        S.loc[pend, 'status'] = np.where(ins, 'in_support', 'out_of_host_support')
        S.loc[S.status != 'in_support', 'p_emp'] = np.nan
    S.loc[S.status != 'in_support', 'Z_cand'] = np.nan          # 支持范围外不报告 Z_s（capped 另给 Z_cap）
    # capped_reference：份额残差 d，相对最高认证档
    capm = S.status == 'capped_reference'
    for c in ('Z_cap', 'p_cap', 'K_ref'): S[c] = np.nan
    if capm.any():
        lo, hi = REFBIN[cert]; ref = W[(W.m >= lo) & (W.m < hi)].d.values; srt = np.sort(ref)
        e = S[capm]; d = (e.n1 - e.E1) / e.m
        S.loc[capm, 'Z_cap'] = e.sign * (d - ref.mean()) / ref.std()
        up = (1 + len(srt) - np.searchsorted(srt, d.values - 1e-12, side='left')) / (1 + len(srt))
        lw = (1 + np.searchsorted(srt, d.values + 1e-12, side='right')) / (1 + len(srt))
        S.loc[capm, 'p_cap'] = np.where(e.sign > 0, up, lw); S.loc[capm, 'K_ref'] = len(ref); S.loc[capm, 'ref_bin'] = f'{lo}-{hi - 1}'
    S['method'] = method; S['window'] = f'W{w}'
    return S


res = []
for sid, Ls in G.SPEC.items():
    lab = G.LABEL[sid]; drug = 'HAR' if 'HAR' in lab else 'CHX'
    qual = [L for L in sorted(Ls) if (sid, L) in G.QUALIFIED]
    pooled = {}
    for L in qual:
        VC = viral_codon_table(sid, L)
        for w in (1, 5):
            null = host_null(sid, L, w)
            mu_fn = lambda ctx, pc=null['pc'], p0=null['p0']: pd.Series(ctx).map(pc).fillna(p0).values
            S = window_stats(VC, CAND, w, mu_fn)
            if w == 5: pooled[L] = (S, null)
            ev = evaluate(S, CAND, null, SUP1.get((sid, L, w), 'none'), w, f'EMP_Zs_L{L}')
            ev['sid'] = sid; ev['label'] = lab; ev['drug'] = drug; ev['read_length'] = L
            ev['periodicity_ok'] = [bool(PER.F0[(sid, L, o)] >= PER.host_F0[(sid, L, o)] - 0.10) if (sid, L, o) in PER.index else False for o in CAND.parent]
            ev['context_ok'] = bool(CTXC.get((sid, L), np.nan) >= 0.5)
            res.append(pd.concat([CAND, ev], axis=1))
        print(sid, L, flush=True)
    # J_POOL（W5，≥ 2 个合格读长）
    if len(qual) >= 2:
        Ws = []
        for L in qual:
            C, pj, pc_all, p0 = TF.host_reference(sid, L); W = G.windows(C, 5, pj); Ws.append(W.set_index(['tx_idx', 'blk'])[['n0', 'n1', 'n2', 'm', 'E1', 'V']])
        idx = Ws[0].index
        for w_ in Ws[1:]: idx = idx.union(w_.index)
        J = sum(w_.reindex(idx).fillna(0) for w_ in Ws).reset_index(); J = J[J.m >= 5].copy()
        J['dens'] = J.n0 / 15; J['z'] = (J.n1 - J.E1) / np.sqrt(J.V.clip(lower=1e-9)); J['lm'] = np.log2(J.m); J['ld'] = np.log2(J.dens + 1 / 15)
        ph = phi_fit(J.lm.values, J.z.values); lbv = np.floor(J.lm / 0.25).astype(int); vc = lbv.value_counts(); okb = vc[vc >= 30].index.values
        J['zs'] = J.z / np.sqrt(ph(J.lm.values)); J['d'] = (J.n1 - J.E1) / J.m
        nullJ = dict(W=J, ph=ph, rng=((okb.min() + .5) * .25, (okb.max() + .5) * .25))
        SJ = sum(pooled[L][0][['n0', 'n1', 'n2', 'E1', 'V']] for L in qual); SJ['excluded'] = pooled[qual[0]][0].excluded
        ev = evaluate(SJ, CAND, nullJ, SUPJ.get((sid, 5), 'none'), 5, 'J_POOL_' + '/'.join(map(str, qual)))
        ev['sid'] = sid; ev['label'] = lab; ev['drug'] = drug; ev['read_length'] = '/'.join(map(str, qual))
        ev['periodicity_ok'] = [all(bool(PER.F0[(sid, L, o)] >= PER.host_F0[(sid, L, o)] - 0.10) if (sid, L, o) in PER.index else False for L in qual) for o in CAND.parent]
        ev['context_ok'] = all(bool(CTXC.get((sid, L), np.nan) >= 0.5) for L in qual)
        res.append(pd.concat([CAND, ev], axis=1))
        print(sid, 'J_POOL', flush=True)
R = pd.concat(res, ignore_index=True)
R['transfer_flag'] = np.where(R.periodicity_ok & R.context_ok, 'checks_passed', 'transfer_unverified')
R['evidence_class'] = 'initiation_associated_occupancy'
R['tail_flag'] = 'tail_uncertified_below_0.01'
R.to_csv(os.path.join(TB, 'HARTIS_v3_PEDV_start_level.tsv.gz'), sep='\t', index=False, float_format='%.5g')
print('rows', len(R)); print(R.groupby(['frame', 'window', 'status']).size().unstack(fill_value=0).to_string())
