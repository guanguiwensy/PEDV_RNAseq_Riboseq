# -*- coding: utf-8 -*-
"""G2d 探索性分析：宿主 CDS 内 +1 框（F1）AUG 处的候选框富集——技术偏差还是与翻译起始一致的特征？
不修改 EMP_Zs；不碰病毒；每个合格库×读长独立；主口径（未剪接 read、tier_main、CDS 内部 j ∈ [10, N−11]、n_sense_tx = 1）。
单位 = 单个密码子 j 的 F1 位（P-site 在 3j+1 = F1 框三联体的第 1 位）。类别：F1 框三联体 = ATG（AUG）/ 9 种 near-cognate / 其他。
输出（G2d_round5/aug_*.tsv）：
  A1 富集分布：各类别 z_s（EMP_Zs 标准化残差，按基因 5 折样本外）分位数、名义 p ≤ 0.05/0.01 比例（m ≥ 10）
  A2 匹配比较：同库×读长内按（端点上下文 256 组 × log2 F0 档 × log2 m 档）分层，AUG − 其他、NC − 其他 的 z_s 均值差与
     p≤0.05 比例差（分层加权），按基因整群 bootstrap 95% 区间
  A3 回归：n1 ~ Binomial(m, logit = 上下文 offset + β_AUG + β_NC + F0/m/GC±15/CDS 位置分档)，按基因整群 bootstrap；
     另含 AUG × Kozak（−3 = A/G 且 +4 = G）与 AUG × +1 框 ORF 长度（< 10、10–49、≥ 50 密码子）交互
  B  位置分布：以 F1 框三联体第 1 位为 0，−30..+60 nt 的逐位点 P-site 计数，按锚点局部 F0 均值归一化后取平均；
     另按"+1 框终止密码子位于 +15..+45 nt"的锚点比较终止前后的 F1 框信号
  C  重现性：同一位点（tx, j）在两个库×读长之间 z_s 的 Spearman 相关与 p ≤ 0.05 共现富集（AUG vs 其他）
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
from scipy.stats import spearmanr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import g2b_common as G
from g2b_calibrate import emp_p, phi_fit, SHR, TRAIN_MMIN, EVAL_MMIN
OUT = os.path.join(G.OUT, 'G2d_round5'); os.makedirs(OUT, exist_ok=True)
B4 = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
NEAR = {'CTG', 'GTG', 'TTG', 'AAG', 'ACG', 'AGG', 'ATA', 'ATC', 'ATT'}
code = lambda s: B4[s[0]] * 16 + B4[s[1]] * 4 + B4[s[2]]
AUG = code('ATG'); NC = np.array([code(c) for c in NEAR]); STOPS = np.array([code(c) for c in ('TAA', 'TAG', 'TGA')])
OFFS = np.arange(-30, 61); rng = np.random.default_rng(G.SEED)


def irls(y, m, X, off, it=40):
    b = np.zeros(X.shape[1])
    for _ in range(it):
        eta = off + X @ b; mu = 1 / (1 + np.exp(-eta)); W = m * mu * (1 - mu)
        z = eta - off + (y - m * mu) / np.maximum(W, 1e-9)
        bn = np.linalg.solve(X.T @ (X * W[:, None]) + 1e-6 * np.eye(X.shape[1]), X.T @ (W * z))
        if np.max(np.abs(bn - b)) < 1e-7: return bn
        b = bn
    return b


def site_table(sid, L):
    """返回逐密码子表（含零计数的内部密码子不需要：只保留有计数者），附类别、样本外 μ、z_s、p、序列特征。"""
    tx = G.tx_table(); SEQ, ST = G._c['seq']; t = tx.set_index('tx_idx')
    C = G.codon_table(sid, L).copy()
    cs = t.cds_start_tx.reindex(C.tx_idx).values; base = ST[C.tx_idx.values] + cs
    p1 = base + 3 * C.j.values + 1                                    # F1 框三联体第 1 位（序列数组下标）
    tri = np.stack([SEQ[p1], SEQ[p1 + 1], SEQ[p1 + 2]]); bad = (tri == 4).any(0)
    tc = np.where(bad, -1, tri[0] * 16 + tri[1] * 4 + tri[2])
    C['cat'] = np.where(tc == AUG, 'AUG', np.where(np.isin(tc, NC), 'near_cognate', np.where(bad, 'N', 'other')))
    C['kozak_strong'] = ((np.isin(SEQ[p1 - 3], [0, 2])) & (SEQ[p1 + 3] == 2)).astype(int)
    win = np.stack([SEQ[p1 + k] for k in range(-15, 16)]); C['gc15'] = ((win == 1) | (win == 2)).mean(0)
    # +1 框 ORF 长度：从该三联体起沿 F1 框找第一个终止密码子（最多 300 个密码子）
    orfl = np.full(len(C), 300)
    alive = np.ones(len(C), bool)
    for k in range(1, 301):
        q = p1 + 3 * k; b3 = np.stack([SEQ[q], SEQ[q + 1], SEQ[q + 2]])
        st = alive & ~(b3 == 4).any(0) & np.isin(b3[0] * 16 + b3[1] * 4 + b3[2], STOPS)
        orfl[st] = k; alive &= ~st
        if not alive.any(): break
    C['plus1_orf_codons'] = orfl
    C['mu'] = np.nan; C['zs'] = np.nan; C['p'] = np.nan
    for f in range(G.NFOLD):
        trm = (C.fold.values != f) & (C.m.values > 0)
        p0 = C.n1.values[trm].sum() / C.m.values[trm].sum()
        g = pd.DataFrame(dict(k=C.ctx.values[trm], n1=C.n1.values[trm], m=C.m.values[trm])).groupby('k').sum()
        pc = (g.n1 + SHR * p0) / (g.m + SHR)
        # 训练单密码子窗的交叉拟合期望
        mu_tr = np.full(len(C), np.nan)
        for g2 in range(G.NFOLD):
            if g2 == f: continue
            sub = (C.fold.values != f) & (C.fold.values != g2) & (C.m.values > 0); q0 = C.n1.values[sub].sum() / C.m.values[sub].sum()
            h = pd.DataFrame(dict(k=C.ctx.values[sub], n1=C.n1.values[sub], m=C.m.values[sub])).groupby('k').sum()
            qc = (h.n1 + SHR * q0) / (h.m + SHR); sel = C.fold.values == g2
            mu_tr[sel] = pd.Series(C.ctx.values[sel]).map(qc).fillna(q0).values
        te = C.fold.values == f
        mu_te = pd.Series(C.ctx.values[te]).map(pc).fillna(p0).values
        C.loc[te, 'mu'] = mu_te
        trw = (C.fold.values != f) & (C.m.values >= TRAIN_MMIN)
        z_tr = (C.n1.values[trw] - C.m.values[trw] * mu_tr[trw]) / np.sqrt(C.m.values[trw] * mu_tr[trw] * (1 - mu_tr[trw]))
        lm_tr = np.log2(C.m.values[trw]); ld_tr = np.log2(C.n0.values[trw] / 3 + 1 / 3)
        ph = phi_fit(lm_tr, z_tr); zs_tr = z_tr / np.sqrt(ph(lm_tr))
        tm = te & (C.m.values >= 1)
        z_te = (C.n1.values[tm] - C.m.values[tm] * C.mu.values[tm]) / np.sqrt(C.m.values[tm] * C.mu.values[tm] * (1 - C.mu.values[tm]))
        lm_te = np.log2(C.m.values[tm]); zs_te = z_te / np.sqrt(ph(lm_te)); C.loc[tm, 'zs'] = zs_te
        tp = te & (C.m.values >= EVAL_MMIN)
        z = (C.n1.values[tp] - C.m.values[tp] * C.mu.values[tp]) / np.sqrt(C.m.values[tp] * C.mu.values[tp] * (1 - C.mu.values[tp]))
        lm = np.log2(C.m.values[tp]); ld = np.log2(C.n0.values[tp] / 3 + 1 / 3)
        p, K, _ = emp_p(zs_tr, lm_tr, ld_tr, z / np.sqrt(ph(lm)), lm, ld, True); C.loc[tp, 'p'] = p
    return C[C.cat != 'N']


def gboot_mean_diff(a, ga, b, gb_, nb=300):
    """按基因整群 bootstrap 的均值差（a − b）95% 区间；a、b 为数值，ga/gb 为基因标签。"""
    da = pd.DataFrame(dict(g=ga, v=a)).groupby('g').v.agg(['sum', 'size']); db = pd.DataFrame(dict(g=gb_, v=b)).groupby('g').v.agg(['sum', 'size'])
    est = []
    for _ in range(nb):
        ia = rng.integers(0, len(da), len(da)); ib = rng.integers(0, len(db), len(db))
        est.append(da['sum'].values[ia].sum() / da['size'].values[ia].sum() - db['sum'].values[ib].sum() / db['size'].values[ib].sum())
    return np.percentile(est, 2.5), np.percentile(est, 97.5)


def analyse(sid, L):
    C = site_table(sid, L); lab = G.LABEL[sid]; drug = 'HAR' if 'HAR' in lab else 'CHX'
    meta = dict(sid=sid, label=lab, drug=drug, L=L)
    A1 = []
    for cat, d in C[C.m >= EVAL_MMIN].groupby('cat'):
        z = d.zs.dropna(); p = d.p.dropna()
        A1.append(dict(**meta, cat=cat, n=len(z), n_genes=d.tx_idx.nunique(), zs_mean=z.mean(), zs_q50=z.median(), zs_q90=z.quantile(.9),
                       zs_q99=z.quantile(.99), le05=(p <= .05).mean(), le01=(p <= .01).mean(), F1share_pooled=d.n1.sum() / d.m.sum(),
                       mu_pooled=(d.mu * d.m).sum() / d.m.sum()))
    # A2 分层匹配
    D = C[(C.m >= EVAL_MMIN) & C.zs.notna() & C.p.notna()].copy()
    D['s'] = D.ctx.astype(str) + '|' + np.floor(np.log2(D.n0 + 1) / 0.5).astype(int).astype(str) + '|' + np.floor(np.log2(D.m) / 0.5).astype(int).astype(str)
    A2 = []
    for cat in ('AUG', 'near_cognate'):
        x = D[D.cat == cat]; o = D[D.cat == 'other']
        so = o.groupby('s').agg(oz=('zs', 'mean'), op=('p', lambda v: (v <= .05).mean()), on=('zs', 'size'))
        xm = x.join(so, on='s'); xm = xm[xm.on >= 5]
        if len(xm) < 20: continue
        dz = (xm.zs - xm.oz).values; dp = ((xm.p <= .05).astype(float) - xm.op).values
        # 整群 bootstrap：按基因重抽 cat 位点（对照用分层均值，视为固定）
        gg = xm.tx_idx.values; ug, inv = np.unique(gg, return_inverse=True)
        sz = np.bincount(inv, weights=dz); sp = np.bincount(inv, weights=dp); cnt = np.bincount(inv).astype(float)
        W = rng.multinomial(len(ug), np.full(len(ug), 1 / len(ug)), size=1000)
        ez = (W @ sz) / (W @ cnt); ep = (W @ sp) / (W @ cnt)
        A2.append(dict(**meta, cat=cat, n_matched=len(xm), n_genes=len(ug), dz_mean=dz.mean(), dz_lo=np.percentile(ez, 2.5), dz_hi=np.percentile(ez, 97.5),
                       d_le05=dp.mean(), d_le05_lo=np.percentile(ep, 2.5), d_le05_hi=np.percentile(ep, 97.5),
                       cat_le05=(xm.p <= .05).mean(), matched_other_le05=xm.op.mean()))
    # A3 回归（m ≥ 5，样本外上下文 offset）
    R = C[(C.m >= 5) & C.mu.notna()].copy()
    off = np.log(R.mu.clip(1e-4, 1 - 1e-4) / (1 - R.mu.clip(1e-4, 1 - 1e-4))).values
    def dummies(v, nm):
        q = pd.qcut(v, 5, labels=False, duplicates='drop'); return pd.get_dummies(q, prefix=nm, drop_first=True).astype(float)
    base_cols = pd.concat([dummies(np.log2(R.n0 + 1), 'n0'), dummies(np.log2(R.m), 'm'), dummies(R.gc15, 'gc'),
                           dummies(R.j / R.N, 'pos')], axis=1).reset_index(drop=True)
    isA = (R.cat == 'AUG').astype(float).values; isN = (R.cat == 'near_cognate').astype(float).values
    kz = R.kozak_strong.values; ol = R.plus1_orf_codons.values
    Xd = pd.DataFrame(dict(AUG=isA, NC=isN, AUGxKozakStrong=isA * kz, AUGxORF10_49=isA * ((ol >= 10) & (ol < 50)), AUGxORF50plus=isA * (ol >= 50)))
    X = np.column_stack([np.ones(len(R)), Xd.values, base_cols.values]); names = ['intercept'] + list(Xd.columns) + list(base_cols.columns)
    y = R.n1.values.astype(float); mm = R.m.values.astype(float)
    b = irls(y, mm, X, off)
    genes = R.tx_idx.values; ug, inv = np.unique(genes, return_inverse=True); order = np.argsort(inv); bounds = np.searchsorted(inv[order], np.arange(len(ug) + 1))
    bs = []
    for _ in range(100):
        pick = rng.integers(0, len(ug), len(ug)); idx = np.concatenate([order[bounds[g]:bounds[g + 1]] for g in pick])
        bs.append(irls(y[idx], mm[idx], X[idx], off[idx]))
    bs = np.array(bs)
    A3 = [dict(**meta, term=nm, coef=b[i], ci_lo=np.percentile(bs[:, i], 2.5), ci_hi=np.percentile(bs[:, i], 97.5), n_codons=len(R),
               n_AUG=int(isA.sum()), n_AUG_kozak_strong=int((isA * kz).sum())) for i, nm in enumerate(names) if nm in Xd.columns]
    # B 位置分布
    tx = G.tx_table().set_index('tx_idx'); Pp = G._c[('pos', sid)]
    Pp = Pp[(Pp.read_length == L) & (Pp.n_sense_tx == 1) & (Pp.splice_class == 0) & Pp.tx_idx.isin(C.tx_idx.unique())]
    dense = {k: d for k, d in Pp.groupby('tx_idx')}
    prof = {c: np.zeros(len(OFFS)) for c in ('AUG', 'near_cognate', 'other')}; nanch = {c: 0 for c in prof}
    stopseg = []
    allA = C.copy()
    SEQ, ST = G._c['seq']
    for k, d in allA.groupby('tx_idx'):
        if k not in dense: continue
        Lc = int(tx.loc[k, 'cds_len_nt']); arr = np.zeros(Lc + 200); q = dense[k]; q = q[(q.rel >= 0) & (q.rel < Lc)]
        np.add.at(arr, q.rel.values + 60, q['count'].values)
        for cat in prof:
            dd = d[d.cat == cat]
            if cat == 'other' and len(dd) > 3: dd = dd.sample(3, random_state=int(k))   # 其他：每基因至多 3 个，控制规模
            if len(dd) == 0: continue
            a0 = 3 * dd.j.values + 1                                              # 锚点 rel
            ok = (a0 - 30 >= 30) & (a0 + 60 <= Lc - 31)
            a0 = a0[ok]
            if len(a0) == 0: continue
            M = arr[(a0[:, None] + OFFS[None, :]) + 60]
            f0 = M[:, ((OFFS + 1) % 3) == 0]                                       # CDS F0 框：相对锚点 −1、+2、…
            loc = f0.mean(1); keep = loc > 0
            if not keep.any(): continue
            prof[cat] += (M[keep] / loc[keep, None]).sum(0); nanch[cat] += int(keep.sum())
            if cat in ('AUG', 'near_cognate'):
                ol = dd.plus1_orf_codons.values[ok][keep]; sel = (ol >= 5) & (ol <= 15)
                for row, olen in zip(M[keep][sel] / loc[keep][sel, None], ol[sel]):
                    stop_nt = 3 * olen                                             # +1 框终止密码子第 1 位相对锚点
                    pre = row[(OFFS >= 3) & (OFFS < stop_nt) & (OFFS % 3 == 0)].mean()
                    post = row[(OFFS >= stop_nt + 3) & (OFFS <= stop_nt + 30) & (OFFS % 3 == 0)].mean()
                    stopseg.append(dict(cat=cat, pre=pre, post=post))
    Bp = []
    for cat in prof:
        if nanch[cat] == 0: continue
        v = prof[cat] / nanch[cat]
        for o_, x in zip(OFFS, v):
            Bp.append(dict(**meta, cat=cat, offset_nt=int(o_), frame_rel_anchor=('candidate(F1)' if o_ % 3 == 0 else 'CDS(F0)' if (o_ + 1) % 3 == 0 else 'F2'),
                           norm_psite=x, n_anchors=nanch[cat]))
    SS = pd.DataFrame(stopseg)
    Bs = [dict(**meta, cat=c, n_anchors=len(d), F1_before_plus1_stop=d.pre.mean(), F1_after_plus1_stop=d.post.mean(),
               ratio_after_before=d.post.mean() / max(d.pre.mean(), 1e-9)) for c, d in SS.groupby('cat')] if len(SS) else []
    sites = C.loc[C.m >= EVAL_MMIN, ['tx_idx', 'j', 'cat', 'zs', 'p', 'm']].assign(sid=sid, L=L)
    return A1, A2, A3, Bp, Bs, sites


if __name__ == '__main__':
    out = {k: [] for k in ('A1', 'A2', 'A3', 'B', 'Bstop')}; S = []
    for sid, L in sorted(G.QUALIFIED):
        r = analyse(sid, L)
        for k, v in zip(('A1', 'A2', 'A3', 'B', 'Bstop'), r[:5]): out[k] += v
        S.append(r[5]); print(sid, L, 'done', flush=True)
    for k, v in out.items(): pd.DataFrame(v).to_csv(os.path.join(OUT, f'aug_{k}.tsv'), sep='\t', index=False, float_format='%.5g')
    S = pd.concat(S); S.to_csv(os.path.join(OUT, 'aug_sites_m10.tsv.gz'), sep='\t', index=False, float_format='%.5g')
    # C 重现性：同一位点在两个库×读长之间
    S['key'] = S.sid + ':' + S.L.astype(str); keys = sorted(S.key.unique()); rows = []
    piv = {k: S[S.key == k].set_index(['tx_idx', 'j']) for k in keys}
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            a, b = piv[keys[i]], piv[keys[j]]; M = a.join(b, how='inner', lsuffix='_a', rsuffix='_b')
            for cat in ('AUG', 'near_cognate', 'other'):
                x = M[M.cat_a == cat]
                if len(x) < 30: continue
                ha = x.p_a <= .05; hb = x.p_b <= .05
                rows.append(dict(pair=f'{keys[i]}~{keys[j]}', same_library=keys[i].split(':')[0] == keys[j].split(':')[0], cat=cat, n_sites=len(x),
                                 spearman_zs=spearmanr(x.zs_a, x.zs_b).correlation, P_hitB_given_hitA=hb[ha].mean() if ha.sum() else np.nan,
                                 P_hitB=hb.mean(), n_hitA=int(ha.sum())))
    pd.DataFrame(rows).to_csv(os.path.join(OUT, 'aug_C_reproducibility.tsv'), sep='\t', index=False, float_format='%.5g')
    print('all done')
