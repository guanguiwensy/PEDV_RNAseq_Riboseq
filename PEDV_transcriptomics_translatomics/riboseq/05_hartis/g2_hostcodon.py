"""G2 公共模块：宿主逐密码子表（逐库）。
密码子单元 j（相对 CDS 起点，F0 P-site 在 rel 3j）：n0 = c[3j]，n1 = c[3j+1]（F1），n2 = c[3j+2]（F2），n0n = c[3j+3]，n0p = c[3j−3]。
29 nt 读段：P-site 在 p 时 5′ 端在 p−12、3′ 端在 p+16（转录本坐标，剪接 reads 也成立）。
上下文：b5F1 = seq[3j+1−12]，b5F2 = seq[3j+2−12]，b3F1 = seq[3j+1+16]，b3F2 = seq[3j+2+16]；P 位密码子 = seq[3j:3j+3]，A 位 = seq[3j+3:3j+6]。
训练集剔除：split-CDS（≤10 nt 缺口）；起始密码子非 ATG；CDS < 300 nt；GTF 层 CDS 重叠（v3_02）；写死的移码基因名单；
含任何 n_sense_tx > 1 位点的转录本；密码子 j < 10（rel < +30）或 j > N − 11（终止前 30 nt 内）。
基因留出：按 gene_id 的稳定哈希二分（seed 20261027），训练 / 留出各半。"""
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
import os, gzip, hashlib, numpy as np, pandas as pd
H = os.path.dirname(os.path.abspath(__file__)); V3 = os.path.join(H, '..'); T = os.path.join(V3, 'tables'); SV = os.path.join(V3, '_server')
LIBS = {'CHX5': SR.lib("RIBO_VERO_CHX_1MOI_5h"), 'HAR5': SR.lib("RIBO_VERO_HAR_1MOI_5h"), 'CHX8': SR.lib("RIBO_VERO_CHX_1MOI_8h"), 'HAR8': SR.lib("RIBO_VERO_HAR_1MOI_8h"), 'CHX16': SR.lib("RIBO_VERO_CHX_1MOI_16h"), 'HAR16': SR.lib("RIBO_VERO_HAR_1MOI_16h")}
FS_GENES = {'OAZ1', 'OAZ2', 'OAZ3', 'AZIN1', 'PEG10', 'PNCK', 'MAP3K14', 'CCR5', 'EDN1', 'PLEKHM2', 'SGIP1', 'SFRS9'}
BASE = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
FL = 60
_c = {}


def tx_table():
    if 'tx' in _c: return _c['tx']
    TX = pd.read_csv(os.path.join(T, 'g2_tx_seq.tsv.gz'), sep='\t')
    G = pd.read_csv(os.path.join(T, 'host_gene_level_29nt_frames.tsv.gz'), sep='\t', usecols=['transcript_id', 'gene_id', 'excl_overlap']).drop_duplicates('transcript_id')
    TX = TX.merge(G, on='transcript_id', how='left')
    TX['excl_overlap'] = TX.excl_overlap.fillna(False).astype(bool)
    TX['heldout'] = [int(hashlib.md5(f'20261027|{g}'.encode()).hexdigest(), 16) % 2 == 1 for g in TX.gene_id.fillna(TX.transcript_id)]
    _c['tx'] = TX
    return TX


def load_pos(lib):
    P = pd.read_csv(os.path.join(SV, LIBS[lib], f'pos_{LIBS[lib]}.tsv.gz'), sep='\t')
    multi = set(P.loc[P.n_sense_tx > 1, 'tx_idx'])
    P = P.groupby(['tx_idx', 'rel'], as_index=False)['count'].sum()       # 合并 sc5（不丢弃）
    return P, multi


def codon_table(lib, include_start_region=False):
    """返回逐密码子 DataFrame（训练资格基因；include_start_region=True 时含 j ≥ 0 的起始近端密码子，用于 G3）。"""
    key = ('cod', lib, include_start_region)
    if key in _c: return _c[key]
    TX = tx_table(); P, multi = load_pos(lib)
    ok = (TX.split_cds_gap_le10 == 0) & (TX.start_codon == 'ATG') & (TX.cds_len_nt >= 300) & (~TX.excl_overlap) & \
         (~TX.gene_name.fillna('').isin(FS_GENES)) & (~TX.tx_idx.isin(multi))
    TXo = TX[ok].set_index('tx_idx')
    P = P[P.tx_idx.isin(TXo.index)]
    rows = []
    grp = {k: v for k, v in P.groupby('tx_idx')}
    LUT = np.full(256, 4, np.int8)
    for ch, v in BASE.items(): LUT[ord(ch)] = v
    for k, r in TXo.iterrows():
        if k not in grp: continue                     # 无任何 29 nt 计数的转录本不贡献信息
        L = int(r.cds_len_nt); cs = int(r.cds_start_tx)
        sa = LUT[np.frombuffer(r.seq.encode(), np.uint8)]
        sa = np.concatenate([np.full(FL, 4, np.int8), sa, np.full(FL, 4, np.int8)])   # 两端补 N，下标 + FL
        a = np.zeros(L + 2 * FL + 6, np.int32)
        d = grp[k]; a[d.rel.values + FL] = d['count'].values
        N = L // 3
        j0 = 0 if include_start_region else 10
        js = np.arange(j0, N - 11)
        if len(js) == 0: continue
        p0 = 3 * js + FL; t0 = cs + 3 * js + FL
        g = lambda off: sa[t0 + off]
        def codon(off):
            b = np.stack([g(off), g(off + 1), g(off + 2)])
            c = b[0].astype(np.int16) * 16 + b[1] * 4 + b[2]
            return np.where((b == 4).any(0), 64, c).astype(np.int16)
        rows.append(pd.DataFrame(dict(tx_idx=k, gene_id=r.gene_id, heldout=r.heldout, j=js, N=N,
                                      n0=a[p0], n1=a[p0 + 1], n2=a[p0 + 2], n0n=a[p0 + 3], n0p=a[p0 - 3],
                                      b5F1=g(1 - 12), b5F2=g(2 - 12), b5F0=g(-12), b3F1=g(1 + 16), b3F2=g(2 + 16),
                                      Pcodon=codon(0), Acodon=codon(3))))
    D = pd.concat(rows, ignore_index=True)
    D['m'] = D.n1 + D.n2
    D['relpos'] = D.j / D.N
    # 局部深度：±5 个密码子的 F0 均值（同一转录本内）
    D['loc_n0'] = D.groupby('tx_idx').n0.transform(lambda x: x.rolling(11, center=True, min_periods=1).mean())
    _c[key] = D
    return D
