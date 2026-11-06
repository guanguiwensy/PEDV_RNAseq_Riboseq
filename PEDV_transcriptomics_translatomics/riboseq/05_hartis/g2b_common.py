# -*- coding: utf-8 -*-
"""G2b 公共模块：宿主逐密码子表（逐库 × 逐读长）、基因分折、窗口聚合。只读 _server 数据，不写任何既有文件。

密码子 j（rel 3j = F0 P-site）：n0 = c[3j]，n1 = c[3j+1]（F1），n2 = c[3j+2]（F2）；只用 n_sense_tx = 1 的位点；
任何一个位点存在 n_sense_tx > 1 记录的密码子整体剔除。CDS 内部：j ∈ [10, N−11]（去起始后 30 nt、终止前 30 nt）。
splice 口径：main = 只用 splice_class 0（未剪接，与病毒 read 一致）；all = 0/1/2 合并（敏感性）。
端点上下文（与 HARTIS_v3 g2_hostcodon 的 Mb53 一致，按读长推广）：读长 L、offset o 时，F1 读段 5′ 端 t = cs + 3j + 1 − o，
  b5F1 = seq[t]，b5F2 = seq[t+1]，b3F1 = seq[t+L−1]，b3F2 = seq[t+L]；ctx = 64·b5F1 + 16·b5F2 + 4·b3F1 + b3F2（含 N 记 256）。
基因 5 折：按 gene_id 随机排列（seed 20261027）后 rank mod 5。
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
import os, gzip, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
V3 = os.path.abspath(os.path.join(HERE, '..', '..'))
SV = os.path.join(V3, '_server'); GB = os.path.join(SV, 'G1b_G2'); AUD = os.path.join(SV, 'G1b_audit')
OUT = os.path.abspath(os.path.join(HERE, '..'))
LIBS = [(SR.lib("RIBO_VERO_CHX_1MOI_5h"), 'CHX5h'), (SR.lib("RIBO_VERO_HAR_1MOI_5h"), 'HAR5h'), (SR.lib("RIBO_VERO_CHX_1MOI_8h"), 'CHX8h'), (SR.lib("RIBO_VERO_HAR_1MOI_8h"), 'HAR8h'),
        (SR.lib("RIBO_VERO_CHX_1MOI_16h"), 'CHX16h'), (SR.lib("RIBO_VERO_HAR_1MOI_16h"), 'HAR16h')]
LABEL = dict(LIBS)
# Module 1 M1_v1.0 expanded 集（读长:offset）
SPEC = {SR.lib("RIBO_VERO_CHX_1MOI_5h"): {29: 12, 30: 12}, SR.lib("RIBO_VERO_HAR_1MOI_5h"): {29: 12}, SR.lib("RIBO_VERO_CHX_1MOI_8h"): {29: 12, 30: 12, 31: 12},
        SR.lib("RIBO_VERO_HAR_1MOI_8h"): {29: 12, 30: 12}, SR.lib("RIBO_VERO_CHX_1MOI_16h"): {29: 12, 30: 12, 31: 12},
        SR.lib("RIBO_VERO_HAR_1MOI_16h"): {29: 12, 30: 12}}
QUALIFIED = {(s, L) for s, d in SPEC.items() for L in d}
FS = {'OAZ1', 'OAZ2', 'OAZ3', 'AZIN1', 'PEG10', 'PNCK', 'MAP3K14', 'CCR5', 'EDN1', 'PLEKHM2', 'SGIP1', 'SFRS9', 'SRSF9'}
SEED = 20261027; NFOLD = 5
_c = {}


def tx_table():
    if 'tx' in _c: return _c['tx']
    tx = pd.read_csv(os.path.join(SV, SR.lib("RIBO_VERO_CHX_1MOI_8h"), 'transcripts.tsv'), sep='\t')
    fl = pd.read_csv(os.path.join(AUD, 'G2_training_flags.tsv'), sep='\t', usecols=['tx_idx', 'tier_main', 'split_class', 'fs_list'])
    tx = tx.merge(fl, on='tx_idx')
    tx['tier_sens'] = (~tx.gene_name.fillna('').isin(FS) & (tx.fs_list == 0)).astype(int)
    genes = np.array(sorted(tx.gene_id.unique()))
    rng = np.random.default_rng(SEED); perm = rng.permutation(len(genes))
    fold = dict(zip(genes[perm], np.arange(len(genes)) % NFOLD))
    tx['fold'] = tx.gene_id.map(fold)
    seqs = {}
    with gzip.open(os.path.join(GB, 'annot', 'tx_seq.fa.gz'), 'rt') as f:
        k = None
        for line in f:
            if line[0] == '>': k = int(line[1:].split()[0])
            else: seqs[k] = line.strip()
    LUT = np.full(256, 4, np.int16)
    for ch, v in {'A': 0, 'C': 1, 'G': 2, 'T': 3}.items(): LUT[ord(ch)] = v
    PAD = 100; starts = np.zeros(len(tx), np.int64); parts = []; cur = 0
    for k in tx.tx_idx:
        s = LUT[np.frombuffer(seqs[k].encode(), np.uint8)]
        parts.append(np.full(PAD, 4, np.int16)); parts.append(s); parts.append(np.full(PAD, 4, np.int16))
        starts[k] = cur + PAD; cur += len(s) + 2 * PAD
    _c['tx'] = tx; _c['seq'] = (np.concatenate(parts), starts)
    return tx


def codon_table(sid, L, splice='main', tier='tier_main'):
    """返回 CDS 内部逐密码子 DataFrame：tx_idx, fold, j, N, n0, n1, n2, m, ctx（只保留 n0+n1+n2 > 0 的密码子，窗口聚合时零密码子不贡献计数）。"""
    tx = tx_table(); SEQ, ST = _c['seq']; off = SPEC[sid][L]
    key = ('pos', sid)
    if key not in _c:
        _c[key] = pd.read_csv(os.path.join(GB, sid, f'pos_multiL_{sid}.tsv.gz'), sep='\t')
    P = _c[key]; P = P[P.read_length == L]
    t = tx.set_index('tx_idx')
    ok = t.index[t[tier] == 1]
    P = P[P.tx_idx.isin(ok) & (P.rel >= 0)]
    Nmap = (t.cds_len_nt // 3)
    P = P.assign(j=P.rel // 3, fr=P.rel % 3)
    P = P[(P.j >= 10) & (P.j <= Nmap.reindex(P.tx_idx).values - 11)]
    bad = set(zip(P.tx_idx[P.n_sense_tx > 1], P.j[P.n_sense_tx > 1]))
    P = P[P.n_sense_tx == 1]
    if splice == 'main': P = P[P.splice_class == 0]
    g = P.groupby(['tx_idx', 'j', 'fr'])['count'].sum().unstack('fr', fill_value=0).reindex(columns=[0, 1, 2], fill_value=0)
    g.columns = ['n0', 'n1', 'n2']; g = g.reset_index()
    if bad:
        bk = pd.MultiIndex.from_tuples(list(bad), names=['tx_idx', 'j'])
        g = g[~pd.MultiIndex.from_frame(g[['tx_idx', 'j']]).isin(bk)]
    g['m'] = g.n1 + g.n2
    g['N'] = Nmap.reindex(g.tx_idx).values; g['fold'] = t.fold.reindex(g.tx_idx).values
    g['gene_id'] = t.gene_id.reindex(g.tx_idx).values
    pos = ST[g.tx_idx.values] + t.cds_start_tx.reindex(g.tx_idx).values + 3 * g.j.values + 1 - off
    b = np.stack([SEQ[pos], SEQ[pos + 1], SEQ[pos + L - 1], SEQ[pos + L]])
    g['ctx'] = np.where((b == 4).any(0), 256, b[0] * 64 + b[1] * 16 + b[2] * 4 + b[3]).astype(np.int16)
    return g.reset_index(drop=True)


def gene_density(C):
    """基因内部 F0 密度（P-site/nt），用于区分基因表达水平与窗口局部峰。"""
    s = C.groupby('tx_idx').agg(n0=('n0', 'sum'), N=('N', 'first'))
    return (s.n0 / (3 * (s.N - 20).clip(lower=1))).rename('gene_dens')


def windows(C, w, pj):
    """把逐密码子聚合为不重叠的 w 密码子窗（只取完整窗）。pj = 每个密码子的上下文期望 F1 份额。"""
    blk = (C.j.values - 10) // w
    nblk = (C.N.values - 21 + 1) // w        # 内部密码子数 = N − 20
    keep = blk < nblk
    D = pd.DataFrame(dict(tx_idx=C.tx_idx.values[keep], fold=C.fold.values[keep], blk=blk[keep],
                          n0=C.n0.values[keep], n1=C.n1.values[keep], n2=C.n2.values[keep], m=C.m.values[keep],
                          E1=(C.m.values * pj)[keep], V=(C.m.values * pj * (1 - pj))[keep]))
    W = D.groupby(['tx_idx', 'blk'], sort=False).agg(fold=('fold', 'first'), n0=('n0', 'sum'), n1=('n1', 'sum'), n2=('n2', 'sum'),
                                                     m=('m', 'sum'), E1=('E1', 'sum'), V=('V', 'sum')).reset_index()
    W['dens'] = W.n0 / (3 * w)
    return W
