# -*- coding: utf-8 -*-
"""G2d AUG 探索的补充（修正两处设计问题）：
B'：位置分布改用**由序列选定的无偏锚点**（全部内部密码子中 F1 框三联体为 AUG / near-cognate / 其他，不以该处是否有 read 为条件；
    第一版 B 只取有计数的密码子，会抬高锚点附近的信号）。每个锚点按其 −30..+60 nt 内 CDS 框（F0）的均值归一化（均值 = 0 的锚点不计）；
    "其他"每基因随机至多 20 个。另比较 +1 框终止前后的 F1 框信号（终止位于 +15..+45 nt 的锚点），并以"其他"锚点为对照。
A2'：匹配比较改为按（log2 F0 档 × log2 m 档）分层（z_s 已含 256 组上下文校正），另一变体再加 16 组 5′ 端二碱基（b5F1, b5F2）。
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
import g2b_common as G
OUT = os.path.join(G.OUT, 'G2d_round5'); rng = np.random.default_rng(G.SEED)
B4 = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
code = lambda s: B4[s[0]] * 16 + B4[s[1]] * 4 + B4[s[2]]
AUG = code('ATG'); NC = np.array([code(c) for c in ('CTG', 'GTG', 'TTG', 'AAG', 'ACG', 'AGG', 'ATA', 'ATC', 'ATT')])
STOPS = np.array([code(c) for c in ('TAA', 'TAG', 'TGA')]); OFFS = np.arange(-30, 61)
S = pd.read_csv(os.path.join(OUT, 'aug_sites_m10.tsv.gz'), sep='\t')
rowsB, rowsS, rowsA = [], [], []
tx = G.tx_table(); t = tx.set_index('tx_idx')
for sid, L in sorted(G.QUALIFIED):
    lab = G.LABEL[sid]; meta = dict(sid=sid, label=lab, drug='HAR' if 'HAR' in lab else 'CHX', L=L)
    SEQ, ST = G._c['seq']
    C = G.codon_table(sid, L)
    # —— A2' ——
    s = S[(S.sid == sid) & (S.L == L)].merge(C[['tx_idx', 'j', 'n0', 'ctx']], on=['tx_idx', 'j'])
    s['b0'] = np.floor(np.log2(s.n0 + 1) / 0.5).astype(int); s['bm'] = np.floor(np.log2(s.m) / 0.5).astype(int); s['c16'] = s.ctx // 16
    for variant, keys in (('F0xm', ['b0', 'bm']), ('F0xm_x_5p16', ['b0', 'bm', 'c16'])):
        st = s.groupby(keys)
        for cat in ('AUG', 'near_cognate'):
            o = s[s.cat == 'other'].groupby(keys).agg(oz=('zs', 'mean'), op=('p', lambda v: (v <= .05).mean()), on=('zs', 'size'))
            x = s[s.cat == cat].join(o, on=keys); x = x[x.on >= 10]
            if len(x) < 20: continue
            dz = (x.zs - x.oz).values; dp = ((x.p <= .05).astype(float) - x.op).values
            ug, inv = np.unique(x.tx_idx.values, return_inverse=True); cnt = np.bincount(inv).astype(float)
            W = rng.multinomial(len(ug), np.full(len(ug), 1 / len(ug)), size=1000)
            ez = (W @ np.bincount(inv, weights=dz)) / (W @ cnt); ep = (W @ np.bincount(inv, weights=dp)) / (W @ cnt)
            rowsA.append(dict(**meta, variant=variant, cat=cat, n_matched=len(x), n_genes=len(ug), dz=dz.mean(), dz_lo=np.percentile(ez, 2.5),
                              dz_hi=np.percentile(ez, 97.5), cat_le05=(x.p <= .05).mean(), matched_other_le05=x.op.mean(), d_le05=dp.mean(),
                              d_le05_lo=np.percentile(ep, 2.5), d_le05_hi=np.percentile(ep, 97.5)))
    # —— B' ——
    P = G._c[('pos', sid)]; P = P[(P.read_length == L) & (P.n_sense_tx == 1) & (P.splice_class == 0)]
    keep_tx = set(C.tx_idx.unique()); prof = {c: np.zeros(len(OFFS)) for c in ('AUG', 'near_cognate', 'other')}; n = {c: 0 for c in prof}
    seg = {c: [] for c in prof}
    for k, d in P[P.tx_idx.isin(keep_tx)].groupby('tx_idx'):
        Lc = int(t.loc[k, 'cds_len_nt']); arr = np.zeros(Lc + 200); d = d[(d.rel >= 0) & (d.rel < Lc)]
        np.add.at(arr, d.rel.values + 60, d['count'].values)
        base = ST[k] + int(t.loc[k, 'cds_start_tx']); N = Lc // 3
        js = np.arange(20, N - 31)                                            # 锚点 −30..+60 nt 完全落在 CDS 内部（去起止近端）
        if len(js) == 0: continue
        a0 = 3 * js + 1; tri = np.stack([SEQ[base + a0], SEQ[base + a0 + 1], SEQ[base + a0 + 2]])
        tc = np.where((tri == 4).any(0), -1, tri[0] * 16 + tri[1] * 4 + tri[2])
        cats = np.where(tc == AUG, 'AUG', np.where(np.isin(tc, NC), 'near_cognate', np.where(tc < 0, 'N', 'other')))
        for c in prof:
            sel = np.where(cats == c)[0]
            if c == 'other' and len(sel) > 20: sel = rng.choice(sel, 20, replace=False)
            if len(sel) == 0: continue
            aa = a0[sel]; M = arr[aa[:, None] + OFFS[None, :] + 60]
            loc = M[:, ((OFFS + 1) % 3) == 0].mean(1); ok = loc > 0
            if not ok.any(): continue
            Mn = M[ok] / loc[ok, None]; prof[c] += Mn.sum(0); n[c] += int(ok.sum())
            # +1 框终止：沿 F1 框找第一个终止（5..15 个密码子内）
            aok = aa[ok]; olen = np.full(len(aok), 999)
            for kk in range(1, 16):
                q = base + aok + 3 * kk; b3 = np.stack([SEQ[q], SEQ[q + 1], SEQ[q + 2]])
                hit = (olen == 999) & ~(b3 == 4).any(0) & np.isin(b3[0] * 16 + b3[1] * 4 + b3[2], STOPS); olen[hit] = kk
            sel2 = (olen >= 5) & (olen <= 15)
            for row, ol in zip(Mn[sel2], olen[sel2]):
                sn = 3 * ol
                seg[c].append((row[(OFFS >= 3) & (OFFS < sn) & (OFFS % 3 == 0)].mean(), row[(OFFS >= sn + 3) & (OFFS <= sn + 30) & (OFFS % 3 == 0)].mean(),
                               row[(OFFS >= 3) & (OFFS < sn) & ((OFFS + 1) % 3 == 0)].mean(), row[(OFFS >= sn + 3) & (OFFS <= sn + 30) & ((OFFS + 1) % 3 == 0)].mean()))
    for c in prof:
        if n[c] == 0: continue
        v = prof[c] / n[c]
        for o_, x in zip(OFFS, v):
            rowsB.append(dict(**meta, cat=c, offset_nt=int(o_), frame=('candidate(F1)' if o_ % 3 == 0 else 'CDS(F0)' if (o_ + 1) % 3 == 0 else 'F2'), norm_psite=x, n_anchors=n[c]))
        if seg[c]:
            a = np.array(seg[c])
            rowsS.append(dict(**meta, cat=c, n_anchors=len(a), F1_before=a[:, 0].mean(), F1_after=a[:, 1].mean(), F1_after_over_before=a[:, 1].mean() / a[:, 0].mean(),
                              F0_before=a[:, 2].mean(), F0_after=a[:, 3].mean(), F0_after_over_before=a[:, 3].mean() / a[:, 2].mean()))
    print(sid, L, n, flush=True)
pd.DataFrame(rowsA).to_csv(os.path.join(OUT, 'aug_A2b_matched.tsv'), sep='\t', index=False, float_format='%.5g')
pd.DataFrame(rowsB).to_csv(os.path.join(OUT, 'aug_B2_profile_unbiased.tsv'), sep='\t', index=False, float_format='%.5g')
pd.DataFrame(rowsS).to_csv(os.path.join(OUT, 'aug_B2_plus1stop.tsv'), sep='\t', index=False, float_format='%.5g')
print('done')
