# -*- coding: utf-8 -*-
"""G2b：宿主背景 P(F1 | F1+F2, F0, 读长, 文库) 的模型比较与按基因 5 折交叉验证校准（只用宿主；不碰病毒候选、不碰 PRICE）。

单元：CDS 内部不重叠 w 密码子窗（w = 1, 5, 50）。m = n1 + n2，dens = n0/(3w)。
每折：用其余 4 折基因拟合，在本折基因上评估（训练窗口 m ≥ 5；评估窗口 m ≥ 10）。
模型（单侧 F1 富集 p 值）：
  BIN_GLOBAL   二项(m, p0)，p0 = 训练集合并 F1 份额（已知偏松的基准线）
  BIN_CTX      二项(m, μ)，μ = Σ m_j p_ctx(j) / m（256 组端点上下文，向 p0 收缩，伪计数 20）
  BB_RHOm      beta-binomial(m, p0, ρ(m))：ρ 按训练窗口 log2 m 分档（宽 0.5）MLE，无上下文
  BB_CTX_RHOm  beta-binomial(m, μ, ρ(m))：同上，均值用上下文期望
  EMP_R        条件经验零分布（预注册主方法，无上下文）：r = n1/m；训练窗口按 |Δlog2 m| ≤ 0.5 且 |Δlog2 dens| ≤ 0.5 匹配
               （0.25 网格 ±2 格；K < 100 时放宽到 ±3、±5 格，仍不足记 not_assessable），p = (1 + #{r_g ≥ r})/(1 + K)
  EMP_Z        同上匹配，统计量为上下文标准化残差 z = (n1 − E1)/√V
  —— 第二轮新增（修正样本内偏差与离散性）——
  EMP_Z_cf     同 EMP_Z，但训练窗口的上下文期望用交叉拟合（训练折 g 的密码子用排除 f、g 两折估计的上下文份额），
               评估窗口用全部训练折估计，二者都对各自数据"样本外"
  BB_CTX_cf    同 BB_CTX_RHOm，ρ(m) 用交叉拟合的训练窗口均值估计
  EMP_R_mid / EMP_Z_cf_mid  并列按 mid-p 处理：p = (#{>} + 0.5·(#{=} + 1))/(K + 1)（检验单元自身计作一个并列）
  —— 第三轮新增（深度相关离散）——
  诊断：读段在位点上成团，z 的方差随 m 近似线性增大；匹配集在 ±0.5 log2 m 内被更多的小 m 窗口主导，
        导致 EMP_Z 两尾偏松、EMP_R 偏保守（同一机制，方向相反）。
  EMP_Zs       z_s = z_cf / √φ̂(m)，φ̂(m) = 训练窗口 z_cf 方差随 log2 m 的分档（宽 0.25，每档 ≥ 30 个窗口）线性插值；
               之后按 EMP_Z 同样方式匹配（±0.5 log2）
  EMP_Zs_tight 同上，匹配半径 ±0.25 log2（不足 K 时放宽到 ±0.5、±0.75）
对调（F2 富集）检验：EMP_R、EMP_Z、BB_CTX_RHOm 及上述新增模型。
输出：每个 (库, 读长, 口径) 一个 tables/calib_<sid>_<L>_<variant>.tsv.gz（逐留出窗口 p 值）与 desc_<...>.tsv（深度描述）。
用法: g2b_calibrate.py <sid> <L> <variant: main|allsplice|tiersens|perm> [windows=1,5,50]
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
import os, sys, time, numpy as np, pandas as pd
from scipy.stats import binom, betabinom
from scipy.special import betaln
from scipy.optimize import minimize_scalar
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import g2b_common as G

SHR = 20.0; GRID = 0.25; KMIN = 100; TRAIN_MMIN = 5; EVAL_MMIN = 10


def bb_rho(n1, m, mu):
    """给定逐单元均值 μ，ρ 的 MLE（bounded，log ρ ∈ [−12, −0.05]）。"""
    n2 = m - n1; mu = np.clip(mu, 1e-4, 1 - 1e-4)
    def nll(lr):
        r = np.exp(lr); s = (1 - r) / r; a = mu * s; b = (1 - mu) * s
        return -(betaln(n1 + a, n2 + b) - betaln(a, b)).sum()
    res = minimize_scalar(nll, bounds=(-12, -0.05), method='bounded', options={'xatol': 1e-3})
    return float(np.exp(res.x))


def rho_by_mbin(tr, mu_tr, rng):
    lb = np.floor(np.log2(tr.m.values) / 0.5).astype(int)
    out = {}
    for b in np.unique(lb):
        idx = np.where(lb == b)[0]
        if len(idx) < 50: continue
        if len(idx) > 30000: idx = rng.choice(idx, 30000, replace=False)
        out[b] = bb_rho(tr.n1.values[idx], tr.m.values[idx], mu_tr[idx])
    return out


def lookup_rho(rho, m):
    lb = np.floor(np.log2(m) / 0.5).astype(int); ks = np.array(sorted(rho)); vs = np.array([rho[k] for k in ks])
    i = np.clip(np.searchsorted(ks, lb), 0, len(ks) - 1); i0 = np.clip(i - 1, 0, len(ks) - 1)
    pick = np.where(np.abs(ks[i0] - lb) < np.abs(ks[i] - lb), i0, i)
    return vs[pick], (lb > ks.max()) | (lb < ks.min())


def bb_sf(n1, m, mu, rho):
    mu = np.clip(mu, 1e-4, 1 - 1e-4); s = (1 - rho) / rho
    return betabinom.sf(n1 - 1, m, mu * s, (1 - mu) * s)


def emp_p(tr_stat, tr_lm, tr_ld, te_stat, te_lm, te_ld, upper=True, mid=False, radii=(2, 3, 5)):
    """局部匹配经验 p：训练单元落在 (log2 m, log2 dens) 网格；评估单元取 ±2/±3/±5 格内的训练单元。"""
    gi = np.floor(tr_lm / GRID).astype(int); gj = np.floor(tr_ld / GRID).astype(int)
    tdf = pd.DataFrame(dict(i=gi, j=gj, s=tr_stat))
    cells = {k: np.sort(v.s.values) for k, v in tdf.groupby(['i', 'j'])}
    ti = np.floor(te_lm / GRID).astype(int); tj = np.floor(te_ld / GRID).astype(int)
    p = np.full(len(te_stat), np.nan); K = np.zeros(len(te_stat), int); rad = np.zeros(len(te_stat), int)
    order = pd.DataFrame(dict(i=ti, j=tj)).groupby(['i', 'j']).indices
    ckeys = np.array(list(cells.keys())) if cells else np.zeros((0, 2), int)
    for (a, b), idx in order.items():
        for r in radii:
            sel = (np.abs(ckeys[:, 0] - a) <= r) & (np.abs(ckeys[:, 1] - b) <= r)
            if not sel.any(): continue
            arr = np.sort(np.concatenate([cells[tuple(k)] for k in ckeys[sel]]))
            if len(arr) >= KMIN:
                x = te_stat[idx]
                ge = len(arr) - np.searchsorted(arr, x - 1e-12, side='left'); gt = len(arr) - np.searchsorted(arr, x + 1e-12, side='right')
                le = np.searchsorted(arr, x + 1e-12, side='right'); lt = np.searchsorted(arr, x - 1e-12, side='left')
                if mid:
                    a_, b_ = (gt, ge - gt) if upper else (lt, le - lt)
                    p[idx] = (a_ + 0.5 * (b_ + 1)) / (1 + len(arr))
                else:
                    p[idx] = (1 + (ge if upper else le)) / (1 + len(arr))
                K[idx] = len(arr); rad[idx] = r
                break
    return p, K, rad


GCAP = 5


def phi_fit(lm, z):
    """φ̂(m)：z 方差随 log2 m 分档（宽 0.25，每档 ≥ 30）的线性插值函数；档不足 2 个时返回 None。"""
    lb = np.floor(lm / 0.25).astype(int); vv = pd.Series(z).groupby(lb).agg(['var', 'size']); vv = vv[vv['size'] >= 30]
    if len(vv) < 2: return None
    xc = (vv.index.values + 0.5) * 0.25; yc = np.maximum(vv['var'].values, 1e-3)
    return lambda x: np.interp(x, xc, yc)


def cap_per_gene(genes, cap, rng):
    s = pd.Series(np.arange(len(genes)), index=genes)
    keep = s.groupby(level=0).apply(lambda x: x.values if len(x) <= cap else rng.choice(x.values, cap, replace=False))
    return np.sort(np.concatenate(keep.values))


def match_diag(tr_lm, tr_ld, tr_gene, te_lm, te_ld, radii):
    """与 emp_p 相同的匹配规则下，每个评估单元的匹配集诊断。"""
    gi = np.floor(tr_lm / GRID).astype(int); gj = np.floor(tr_ld / GRID).astype(int)
    cells = {k: v.g.values for k, v in pd.DataFrame(dict(i=gi, j=gj, g=tr_gene)).groupby(['i', 'j'])}
    ckeys = np.array(list(cells.keys())) if cells else np.zeros((0, 2), int)
    ti = np.floor(te_lm / GRID).astype(int); tj = np.floor(te_ld / GRID).astype(int)
    n = len(te_lm); out = {k: np.full(n, np.nan) for k in ('K', 'n_genes', 'max_gene_frac', 'eff_genes', 'radius')}
    for (a, b), idx in pd.DataFrame(dict(i=ti, j=tj)).groupby(['i', 'j']).indices.items():
        for r in radii:
            sel = (np.abs(ckeys[:, 0] - a) <= r) & (np.abs(ckeys[:, 1] - b) <= r)
            if not sel.any(): continue
            g = np.concatenate([cells[tuple(k)] for k in ckeys[sel]])
            if len(g) >= KMIN:
                _, c = np.unique(g, return_counts=True); f = c / c.sum()
                out['K'][idx] = len(g); out['n_genes'][idx] = len(c); out['max_gene_frac'][idx] = f.max()
                out['eff_genes'][idx] = 1.0 / (f ** 2).sum(); out['radius'][idx] = r
                break
    return out


def run(sid, L, variant, wins):
    t0 = time.time(); rng = np.random.default_rng(G.SEED)
    splice = 'all' if variant == 'allsplice' else 'main'; tier = 'tier_sens' if variant == 'tiersens' else 'tier_main'
    C = G.codon_table(sid, L, splice=splice, tier=tier)
    if variant == 'perm':                       # 标签置换：逐密码子随机交换 n1/n2，之后全部流程不变
        sw = rng.random(len(C)) < 0.5
        n1 = np.where(sw, C.n2, C.n1); n2 = np.where(sw, C.n1, C.n2); C = C.assign(n1=n1, n2=n2)
    gd = G.gene_density(C)
    out_p, out_d = [], []
    for f in range(G.NFOLD):
        trc = C[(C.fold != f) & (C.m > 0)]
        p0 = trc.n1.sum() / trc.m.sum()
        g = trc.groupby('ctx')[['n1', 'm']].sum(); pc = (g.n1 + SHR * p0) / (g.m + SHR)
        pj = pc.reindex(C.ctx).fillna(p0).values
        pj_cf = pj.copy()                       # 交叉拟合：训练折 g 的密码子用排除 f、g 两折的上下文份额
        for g2 in range(G.NFOLD):
            if g2 == f: continue
            sub = C[(C.fold != f) & (C.fold != g2) & (C.m > 0)]; q0 = sub.n1.sum() / sub.m.sum()
            h = sub.groupby('ctx')[['n1', 'm']].sum(); qc = (h.n1 + SHR * q0) / (h.m + SHR)
            sel = (C.fold == g2).values; pj_cf[sel] = qc.reindex(C.ctx[sel]).fillna(q0).values
        for w in wins:
            W = G.windows(C, w, pj)
            Wc = G.windows(C, w, pj_cf)[['tx_idx', 'blk', 'E1', 'V']].rename(columns={'E1': 'E1cf', 'V': 'Vcf'})
            W = W.merge(Wc, on=['tx_idx', 'blk'])
            W = W[W.m > 0].copy()
            # 训练窗口用交叉拟合期望，评估窗口用全部训练折期望（都是样本外）
            trmask = (W.fold != f).values
            W['mu_cf'] = np.where(trmask, W.E1cf / W.m, W.E1 / W.m)
            W['z_cf'] = np.where(trmask, (W.n1 - W.E1cf) / np.sqrt(W.Vcf.clip(lower=1e-9)), (W.n1 - W.E1) / np.sqrt(W.V.clip(lower=1e-9)))
            W['mu'] = W.E1 / W.m; W['z'] = (W.n1 - W.E1) / np.sqrt(W.V.clip(lower=1e-9)); W['r'] = W.n1 / W.m
            W['lm'] = np.log2(W.m); W['ld'] = np.log2(W.dens + 1.0 / (3 * w))
            tr = W[(W.fold != f) & (W.m >= TRAIN_MMIN)]; te = W[(W.fold == f) & (W.m >= EVAL_MMIN)].copy()
            if len(te) == 0 or len(tr) < KMIN: continue
            te['p_BIN_GLOBAL'] = binom.sf(te.n1 - 1, te.m, p0)
            te['p_BIN_CTX'] = binom.sf(te.n1 - 1, te.m, te.mu)
            rho0 = rho_by_mbin(tr, np.full(len(tr), p0), rng); rhoc = rho_by_mbin(tr, tr.mu.values, rng)
            r0, ext0 = lookup_rho(rho0, te.m.values); rc, extc = lookup_rho(rhoc, te.m.values)
            te['rho_noctx'] = r0; te['rho_ctx'] = rc; te['rho_extrapolated'] = extc.astype(int)
            te['p_BB_RHOm'] = bb_sf(te.n1.values, te.m.values, np.full(len(te), p0), r0)
            te['p_BB_CTX_RHOm'] = bb_sf(te.n1.values, te.m.values, te.mu.values, rc)
            te['pswap_BB_CTX_RHOm'] = bb_sf(te.n2.values, te.m.values, 1 - te.mu.values, rc)
            rcf = rho_by_mbin(tr, tr.mu_cf.values, rng); rcf_, _ = lookup_rho(rcf, te.m.values); te['rho_ctx_cf'] = rcf_
            te['p_BB_CTX_cf'] = bb_sf(te.n1.values, te.m.values, te.mu.values, rcf_)
            te['pswap_BB_CTX_cf'] = bb_sf(te.n2.values, te.m.values, 1 - te.mu.values, rcf_)
            for nm, st in (('EMP_R', 'r'), ('EMP_Z', 'z')):
                p, K, rad = emp_p(tr[st].values, tr.lm.values, tr.ld.values, te[st].values, te.lm.values, te.ld.values, True)
                ps, _, _ = emp_p(tr[st].values, tr.lm.values, tr.ld.values, te[st].values, te.lm.values, te.ld.values, False)
                te[f'p_{nm}'] = p; te[f'pswap_{nm}'] = ps; te[f'K_{nm}'] = K; te[f'rad_{nm}'] = rad
            for nm, st, mid in (('EMP_R_mid', 'r', True), ('EMP_Z_cf', 'z_cf', False), ('EMP_Z_cf_mid', 'z_cf', True)):
                p, K, rad = emp_p(tr[st].values, tr.lm.values, tr.ld.values, te[st].values, te.lm.values, te.ld.values, True, mid)
                ps, _, _ = emp_p(tr[st].values, tr.lm.values, tr.ld.values, te[st].values, te.lm.values, te.ld.values, False, mid)
                te[f'p_{nm}'] = p; te[f'pswap_{nm}'] = ps
            # 深度相关离散 φ(m)：训练窗口 z_cf 方差随 log2 m 分档
            lb = np.floor(tr.lm.values / 0.25).astype(int); vv = pd.Series(tr.z_cf.values).groupby(lb).agg(['var', 'size'])
            vv = vv[vv['size'] >= 30]
            if len(vv) >= 2:
                xc = (vv.index.values + 0.5) * 0.25; yc = np.maximum(vv['var'].values, 1e-3)
                phi_tr = np.interp(tr.lm.values, xc, yc); phi_te = np.interp(te.lm.values, xc, yc)
                te['phi_m'] = phi_te; te['phi_extrapolated'] = ((te.lm.values > xc.max()) | (te.lm.values < xc.min())).astype(int)
                zs_tr = tr.z_cf.values / np.sqrt(phi_tr); zs_te = te.z_cf.values / np.sqrt(phi_te)
                for nm, radii in (('EMP_Zs', (2, 3, 5)), ('EMP_Zs_tight', (1, 2, 3))):
                    p, K, rad = emp_p(zs_tr, tr.lm.values, tr.ld.values, zs_te, te.lm.values, te.ld.values, True, False, radii)
                    ps, _, _ = emp_p(zs_tr, tr.lm.values, tr.ld.values, zs_te, te.lm.values, te.ld.values, False, False, radii)
                    te[f'p_{nm}'] = p; te[f'pswap_{nm}'] = ps; te[f'K_{nm}'] = K
                # —— 第四轮（冻结前检查）——
                # (a) EMP_Zs 匹配集诊断：基因数、单基因最大占比、有效基因数（逆 Simpson）、匹配半径
                dg = match_diag(tr.lm.values, tr.ld.values, tr.tx_idx.values, te.lm.values, te.ld.values, (2, 3, 5))
                for k_, v_ in dg.items(): te[f'EMP_Zs_{k_}'] = v_
                # (b) 每基因至多 GCAP 个背景窗口（随机抽取），检验少数高表达基因是否主导
                cap_idx = cap_per_gene(tr.tx_idx.values, GCAP, rng)
                p, K, _ = emp_p(zs_tr[cap_idx], tr.lm.values[cap_idx], tr.ld.values[cap_idx], zs_te, te.lm.values, te.ld.values, True)
                ps, _, _ = emp_p(zs_tr[cap_idx], tr.lm.values[cap_idx], tr.ld.values[cap_idx], zs_te, te.lm.values, te.ld.values, False)
                te['p_EMP_Zs_gcap'] = p; te['pswap_EMP_Zs_gcap'] = ps; te['K_EMP_Zs_gcap'] = K
            else:
                for nm in ('EMP_Zs', 'EMP_Zs_tight', 'EMP_Zs_gcap'): te[f'p_{nm}'] = np.nan; te[f'pswap_{nm}'] = np.nan
            # (c) 去掉端点上下文的 EMP_Zs（μ = 训练折库常数 p0），用于检查上下文是否吸收起始相关信号
            z0_tr = (tr.n1.values - tr.m.values * p0) / np.sqrt(tr.m.values * p0 * (1 - p0))
            z0_te = (te.n1.values - te.m.values * p0) / np.sqrt(te.m.values * p0 * (1 - p0))
            ph = phi_fit(tr.lm.values, z0_tr)
            if ph is not None:
                zs0_tr = z0_tr / np.sqrt(ph(tr.lm.values)); zs0_te = z0_te / np.sqrt(ph(te.lm.values))
                p, K, _ = emp_p(zs0_tr, tr.lm.values, tr.ld.values, zs0_te, te.lm.values, te.ld.values, True)
                ps, _, _ = emp_p(zs0_tr, tr.lm.values, tr.ld.values, zs0_te, te.lm.values, te.ld.values, False)
                te['p_EMP_Zs_noctx'] = p; te['pswap_EMP_Zs_noctx'] = ps
            else:
                te['p_EMP_Zs_noctx'] = np.nan; te['pswap_EMP_Zs_noctx'] = np.nan
            te['w'] = w; te['p0'] = p0
            te['gene_dens'] = gd.reindex(te.tx_idx).values
            out_p.append(te.drop(columns=['lm', 'ld', 'E1', 'V', 'E1cf', 'Vcf']))
            if f == 0:   # 深度描述只在第 0 折的全部窗口上做一次（不涉及 p 值）
                A = W[W.m >= TRAIN_MMIN].copy(); A['gene_dens'] = gd.reindex(A.tx_idx).values
                A['mbin'] = pd.cut(A.m, [5, 10, 30, 100, 300, 1000, np.inf], right=False)
                A['dbin'] = pd.qcut(A.dens.rank(method='first'), 5, labels=False)
                for (mb, db), d in A.groupby(['mbin', 'dbin'], observed=True):
                    zz = d.z.values
                    out_d.append(dict(w=w, mbin=str(mb), dens_quintile=int(db), n=len(d), n_genes=d.tx_idx.nunique(),
                                      m_mean=d.m.mean(), dens_median=d.dens.median(), F1share_pooled=d.n1.sum() / d.m.sum(),
                                      F1share_ctx_expected=d.E1.sum() / d.m.sum(), z_mean=zz.mean(), z_var=zz.var(ddof=1),
                                      rho_mom=max(0.0, (zz.var(ddof=1) - 1) / max(d.m.mean() - 1, 1e-9))))
    P = pd.concat(out_p); P.insert(0, 'L', L); P.insert(0, 'sid', sid); P.insert(2, 'variant', variant)
    os.makedirs(os.path.join(G.OUT, 'tables', 'per_window'), exist_ok=True)
    P.to_csv(os.path.join(G.OUT, 'tables', 'per_window', f'calib_{sid}_{L}_{variant}.tsv.gz'), sep='\t', index=False, float_format='%.6g')
    if out_d:
        D = pd.DataFrame(out_d); D.insert(0, 'L', L); D.insert(0, 'sid', sid); D.insert(2, 'variant', variant)
        D.to_csv(os.path.join(G.OUT, 'tables', 'per_window', f'desc_{sid}_{L}_{variant}.tsv'), sep='\t', index=False)
    print(sid, L, variant, 'codons', len(C), 'eval windows', len(P), 'elapsed', round(time.time() - t0), flush=True)


if __name__ == '__main__':
    sid, L, variant = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    wins = [int(x) for x in (sys.argv[4] if len(sys.argv) > 4 else '1,5,50').split(',')]
    run(sid, L, variant, wins)
