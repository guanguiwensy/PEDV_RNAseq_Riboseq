"""G2-02 逐位点宿主模型：各因素对 F1/F2 分配的解释能力（按基因留出，训练基因拟合，留出基因评估）。
评估量：留出密码子（m ≥ 1）上的二项对数似然 LL = Σ[n1 log p + n2 log(1−p)]；
  McFadden 伪 R² = 1 − LL(模型)/LL(M0)；每 read 的 ΔLL。
模型：
  M0     常数 p
  Mb5    F1/F2 读段 5′ 端碱基对 (b5F1, b5F2)（16 组）
  Mb53   5′ 与 3′ 端碱基 (b5F1, b5F2, b3F1, b3F2)（256 组，向 M0 收缩，伪计数 20 reads）
  Mcod   P 位密码子 × A 位密码子（4,096 组，收缩）
  Mocc   占据：λ1 = k·(n0+ε)，λ2 = (n0n+ε)，k 在训练集上网格搜索
  Mocc_alt 占据（对照）：λ1 = k·(n0n+ε)，λ2 = (n0+ε)（检验 F1/F2 各自跟随哪个相邻密码子）
  Moccb  占据 × 端点碱基速率：λ1 = a(b5F1,b3F1)·(n0+ε)，λ2 = b(b5F2,b3F2)·(n0n+ε)
  Mfull  Moccb × (CDS 相对位置 5 档 × 局部深度 5 档) 的赔率校正
ε = 0.5。所有分组参数只用训练基因。
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
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import g2_hostcodon as G
T = G.T; EPS = 0.5; SHR = 20.0


def ll(n1, n2, p):
    p = np.clip(p, 1e-6, 1 - 1e-6); return float((n1 * np.log(p) + n2 * np.log(1 - p)).sum())


def grp_p(tr, te, keys, p0):
    g = tr.groupby(keys)[['n1', 'm']].sum()
    g['p'] = (g.n1 + SHR * p0) / (g.m + SHR)
    return te.join(g.p, on=keys).p.fillna(p0).values


def rates(tr, keyF1, keyF2):
    a = tr.groupby(keyF1).apply(lambda d: (d.n1.sum() + 1) / (d.n0.sum() + EPS * len(d) + 1), include_groups=False)
    b = tr.groupby(keyF2).apply(lambda d: (d.n2.sum() + 1) / (d.n0n.sum() + EPS * len(d) + 1), include_groups=False)
    return a, b


rows, mech = [], []
for lib in G.LIBS:
    D = G.codon_table(lib)
    D = D[(D.b5F1 < 4) & (D.b5F2 < 4) & (D.b3F1 < 4) & (D.b3F2 < 4)]
    tr, te = D[~D.heldout], D[D.heldout]
    tr_m, te = tr[tr.m > 0], te[te.m > 0].copy()
    p0 = tr_m.n1.sum() / tr_m.m.sum()
    preds = {'M0': np.full(len(te), p0)}
    preds['Mb5'] = grp_p(tr_m, te, ['b5F1', 'b5F2'], p0)
    preds['Mb53'] = grp_p(tr_m, te, ['b5F1', 'b5F2', 'b3F1', 'b3F2'], p0)
    preds['Mcod'] = grp_p(tr_m, te, ['Pcodon', 'Acodon'], p0)
    best = None
    for k in np.exp(np.linspace(-2.5, 2.5, 101)):
        v = ll(tr_m.n1, tr_m.n2, k * (tr_m.n0 + EPS) / (k * (tr_m.n0 + EPS) + tr_m.n0n + EPS))
        if best is None or v > best[0]: best = (v, k)
    k = best[1]; preds['Mocc'] = k * (te.n0 + EPS) / (k * (te.n0 + EPS) + te.n0n + EPS)
    best2 = None
    for k2 in np.exp(np.linspace(-2.5, 2.5, 101)):
        v = ll(tr_m.n1, tr_m.n2, k2 * (tr_m.n0n + EPS) / (k2 * (tr_m.n0n + EPS) + tr_m.n0 + EPS))
        if best2 is None or v > best2[0]: best2 = (v, k2)
    preds['Mocc_alt'] = best2[1] * (te.n0n + EPS) / (best2[1] * (te.n0n + EPS) + te.n0 + EPS)
    trc = tr[(tr.n0 + tr.n0n) > 0]
    a, b = rates(trc, ['b5F1', 'b3F1'], ['b5F2', 'b3F2'])
    la = te.join(a.rename('a'), on=['b5F1', 'b3F1']).a.values * (te.n0.values + EPS)
    lb = te.join(b.rename('b'), on=['b5F2', 'b3F2']).b.values * (te.n0n.values + EPS)
    preds['Moccb'] = la / (la + lb)
    # Mfull：位置 × 深度 赔率校正（训练集上估计）
    def bins(d):
        return pd.cut(d.relpos, [0, .2, .4, .6, .8, 1.01], labels=False), pd.qcut(np.log1p(d.loc_n0), 5, labels=False, duplicates='drop')
    la_tr = tr_m.join(a.rename('a'), on=['b5F1', 'b3F1']).a.values * (tr_m.n0.values + EPS)
    lb_tr = tr_m.join(b.rename('b'), on=['b5F2', 'b3F2']).b.values * (tr_m.n0n.values + EPS)
    ptr = la_tr / (la_tr + lb_tr)
    rb, db = bins(tr_m); tmp = pd.DataFrame(dict(rb=rb.values, db=db.values, n1=tr_m.n1.values, e1=tr_m.m.values * ptr, m=tr_m.m.values, ep=ptr))
    cal = tmp.groupby(['rb', 'db']).apply(lambda d: (d.n1.sum() / max(d.m.sum() - d.n1.sum(), 1)) / (d.e1.sum() / max(d.m.sum() - d.e1.sum(), 1)), include_groups=False)
    rbt, dbt = bins(te)
    q = np.quantile(np.log1p(tr_m.loc_n0), [.2, .4, .6, .8])
    dbt = np.searchsorted(q, np.log1p(te.loc_n0.values))
    adj = np.array([cal.get((r_, d_), 1.0) for r_, d_ in zip(rbt.values, dbt)])
    od = preds['Moccb'] / (1 - preds['Moccb']) * adj; preds['Mfull'] = od / (1 + od)
    L0 = ll(te.n1, te.n2, preds['M0']); nreads = te.m.sum()
    for mname, p in preds.items():
        L = ll(te.n1, te.n2, p)
        rows.append(dict(library=lib, model=mname, test_codons=len(te), test_reads=int(nreads), LL=L, pseudoR2_vs_M0=1 - L / L0,
                         dLL_per_read=(L - L0) / nreads, k_occ=k if mname == 'Mocc' else (best2[1] if mname == 'Mocc_alt' else np.nan)))
    # 机制：端点碱基的速率（F1 用 b5F1，F2 用 b5F2；F0 读段比例按 b5F0）
    for base, nm in enumerate('ACGT'):
        d1 = trc[trc.b5F1 == base]; d2 = trc[trc.b5F2 == base]
        d3 = trc[trc.b3F1 == base]; d4 = trc[trc.b3F2 == base]
        mech.append(dict(library=lib, base=nm, rate_F1_by_5p=d1.n1.sum() / max(d1.n0.sum(), 1), rate_F2_by_5p=d2.n2.sum() / max(d2.n0n.sum(), 1),
                         rate_F1_by_3p=d3.n1.sum() / max(d3.n0.sum(), 1), rate_F2_by_3p=d4.n2.sum() / max(d4.n0n.sum(), 1)))
    print(lib, 'done', flush=True)
R = pd.DataFrame(rows); R.to_csv(os.path.join(T, 'g2_model_comparison_heldout.tsv'), sep='\t', index=False, float_format='%.5g')
M = pd.DataFrame(mech); M.to_csv(os.path.join(T, 'g2_end_base_rates.tsv'), sep='\t', index=False, float_format='%.5g')
pd.set_option('display.width', 250)
print(R.pivot_table(index='library', columns='model', values='pseudoR2_vs_M0').round(4)[['Mb5', 'Mb53', 'Mcod', 'Mocc', 'Mocc_alt', 'Moccb', 'Mfull']].to_string())
print(M.pivot_table(index='base', columns='library', values=['rate_F1_by_5p', 'rate_F2_by_5p']).round(3).to_string())
