"""v3-02 宿主 29 nt F1/F2 泄漏的经验关系（非参数分箱，不假定函数形式）。
两个本地可用的数据层：
 G 基因层（body 尺度）：host_gene_cds_psite.tsv（逐转录本 × 读长 × 读框，整段 CDS）。协变量：文库、n0 密度（n0/CDS 长）、m=n1+n2、CDS 长、GC（cds.all.fa）。
 P 起始近端位点层：host_position_level.tsv.gz（rel −30..+90，逐位点 × 读长；只含 depth_0_60≥10 的基因、计数>0 的行）。
   协变量：距起始（密码子）、局部 n0、5′ 端碱基（F1 读段 5′ 端 = P−12）。065 在该表中只有 30 nt 行，**29 nt 无数据**。
训练集剔除规则（本地可执行部分）：CDS < 300 nt；与另一 canonical CDS 在基因组上重叠（同链或反链）；起始近端 rel < +30（P 层）；
  已知移码基因 OAZ1/OAZ2/AZIN1/PEG10 等（名单见 EXCLUDED_GENES）；uORF 区：注释中无 uORF，**无法剔除**（如实报告）；
  终止近端：本地无逐位点数据，**基因层无法剔除终止近端 ±30**（如实报告）。
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
import os, re, gzip, collections, numpy as np, pandas as pd
H = os.path.dirname(os.path.abspath(__file__)); V3 = os.path.join(H, '..'); R = os.path.join(V3, '..'); T = os.path.join(V3, 'tables')
RG = os.path.join(R, 'riboseq_GIb_mapping_psite'); ROOT = os.path.join(R, '..')
LIB = {SR.lib("RIBO_VERO_CHX_1MOI_5h"): 'CHX5', SR.lib("RIBO_VERO_CHX_1MOI_8h"): 'CHX8', SR.lib("RIBO_VERO_CHX_1MOI_16h"): 'CHX16', SR.lib("RIBO_VERO_HAR_1MOI_5h"): 'HAR5', SR.lib("RIBO_VERO_HAR_1MOI_8h"): 'HAR8', SR.lib("RIBO_VERO_HAR_1MOI_16h"): 'HAR16'}
EXCLUDED_GENES = {'OAZ1', 'OAZ2', 'OAZ3', 'AZIN1', 'PEG10', 'PNCK', 'MAP3K14', 'CCR5', 'EDN1', 'PLEKHM2', 'SGIP1', 'SFRS9'}  # 已知/报道过的程序性移码或移码候选（名单写死）
MINCDS = 300

# ---- 注释：canonical CDS 基因组区间，用于重叠剔除；GC 来自 cds.all.fa ----
gtf = os.path.join(ROOT, 'refdata_vero', 'Chlorocebus_sabaeus.ChlSab1.1.116.gtf.gz')
cds_iv = collections.defaultdict(list); tchrom = {}
for line in gzip.open(gtf, 'rt'):
    if line.startswith('#'): continue
    f = line.split('\t')
    if f[2] != 'CDS' or 'Ensembl_canonical' not in f[8] or 'transcript_biotype "protein_coding"' not in f[8]: continue
    t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
    cds_iv[t].append((int(f[3]), int(f[4]))); tchrom[t] = f[0]
# 2026-10-26 修复：v1 的片段相交谓词写成 `x1 <= y2 and x2 >= y1`（错误，旧版见 v3_02_leakage_descriptive_v1_buggy_overlap.py）。
# 现改用 v3_overlap.segments_intersect（闭区间 x1 <= y2 且 x2 <= y1），单元测试见 test_v3_overlap.py。
import sys; sys.path.insert(0, H)
from v3_overlap import cds_overlap_set
overlap = cds_overlap_set(cds_iv, tchrom)
gc = {}
seq = None; tid = None
for line in gzip.open(os.path.join(ROOT, 'refdata_vero_deliver', 'Chlorocebus_sabaeus.ChlSab1.1.cds.all.fa.gz'), 'rt'):
    if line.startswith('>'):
        if tid and seq: s = ''.join(seq).upper(); gc[tid] = (s.count('G') + s.count('C')) / max(len(s), 1)
        tid = line[1:].split()[0].split('.')[0]; seq = []
    else: seq.append(line.strip())
if tid and seq: s = ''.join(seq).upper(); gc[tid] = (s.count('G') + s.count('C')) / max(len(s), 1)

# ---- G 基因层 ----
rows = []
for s, lab in LIB.items():
    g = pd.read_csv(os.path.join(RG, '_work', 'server', 'host_psite', s, 'host_gene_cds_psite.tsv'), sep='\t')
    g = g[g.L == 29].copy()
    g['library'] = lab; g['sample_id'] = s
    rows.append(g)
G = pd.concat(rows)
G['n0'] = G.f0; G['n1'] = G.f1; G['n2'] = G.f2; G['m'] = G.n1 + G.n2
G['n0_density'] = G.n0 / G.cds_len_nt; G['GC'] = G.transcript_id.map(gc)
G['excl_short'] = G.cds_len_nt < MINCDS; G['excl_overlap'] = G.transcript_id.isin(overlap); G['excl_frameshift_gene'] = G.gene_name.isin(EXCLUDED_GENES)
G['train_eligible'] = ~(G.excl_short | G.excl_overlap | G.excl_frameshift_gene)
G.to_csv(os.path.join(T, 'host_gene_level_29nt_frames.tsv.gz'), sep='\t', index=False, float_format='%.5g')
cnt = G.groupby('library').agg(genes_with_29nt=('transcript_id', 'size'), excl_short=('excl_short', 'sum'), excl_overlap=('excl_overlap', 'sum'),
                               excl_frameshift=('excl_frameshift_gene', 'sum'), eligible=('train_eligible', 'sum'),
                               eligible_m_ge10=('m', lambda x: 0)).reset_index()
cnt['eligible_m_ge10'] = [int(((G.library == l) & G.train_eligible & (G.m >= 10)).sum()) for l in cnt.library]
cnt['eligible_m_ge30'] = [int(((G.library == l) & G.train_eligible & (G.m >= 30)).sum()) for l in cnt.library]
cnt.to_csv(os.path.join(T, 'leakage_training_set_counts_gene_level.tsv'), sep='\t', index=False)
print(cnt.to_string())

E = G[G.train_eligible & (G.m >= 10)].copy()
out = []
for lab, d in E.groupby('library'):
    pooled = d.n1.sum() / d.m.sum()
    for var, nb in [('n0_density', 10), ('m', 10), ('cds_len_nt', 5), ('GC', 5)]:
        q = pd.qcut(d[var], nb, duplicates='drop')
        for b, dd in d.groupby(q, observed=True):
            f = dd.n1 / dd.m
            # 过度离散（矩估计 ρ）：Var(f) 相对二项方差的超出
            pbar = dd.n1.sum() / dd.m.sum(); mbar = dd.m.mean()
            vb = (pbar * (1 - pbar) / dd.m).mean(); vo = f.var()
            rho = max(0.0, (vo - vb) / (pbar * (1 - pbar) - vb)) if pbar * (1 - pbar) - vb > 0 else np.nan
            out.append(dict(library=lab, covariate=var, bin=str(b), n_genes=len(dd), pooled_F1_share=pbar,
                            median_gene_F1_share=f.median(), q10=f.quantile(.1), q90=f.quantile(.9), rho_overdispersion=rho,
                            library_pooled=pooled))
O = pd.DataFrame(out); O.to_csv(os.path.join(T, 'leakage_gene_level_bins.tsv'), sep='\t', index=False, float_format='%.4g')
pd.set_option('display.width', 250)
print(O[O.covariate == 'n0_density'][['library', 'bin', 'n_genes', 'pooled_F1_share', 'median_gene_F1_share', 'q10', 'q90', 'rho_overdispersion']].to_string())

# ---- P 起始近端位点层 ----
PL = pd.read_csv(os.path.join(R, 'host_start_proximal_HAR', 'results', 'host_position_level.tsv.gz'), sep='\t', comment='#')
PL = PL[PL.read_length == 29]
SEQ = pd.read_csv(os.path.join(R, 'HARTIS_design', 'tables', 'host_start_seq_m30_p90.tsv'), sep='\t').drop_duplicates('gene_id').set_index('gene_id')
libmap = {'CHX_5h': 'CHX5', 'CHX_8h': 'CHX8', 'CHX_16h': 'CHX16', 'HAR_5h': 'HAR5', 'HAR_8h': 'HAR8', 'HAR_16h': 'HAR16'}
po = []
for lib, d in PL.groupby('library'):
    lab = libmap[lib]
    piv = d.pivot_table(index='gene_id', columns='rel_pos', values='psite_count', aggfunc='sum', fill_value=0).reindex(columns=range(-30, 91), fill_value=0)
    piv = piv[piv.index.isin(SEQ.index)]
    a = piv.values
    # 逐密码子 j（P-site 3j 为 F0）：n0=a[3j]，n1=a[3j+1]，n2=a[3j+2]（同一密码子内的 +1、+2 位）
    for j in range(-9, 29):
        c0 = 3 * j + 30
        if c0 + 2 > 120: break
        n0, n1, n2 = a[:, c0].sum(), a[:, c0 + 1].sum(), a[:, c0 + 2].sum()
        po.append(dict(library=lab, codon_index=j, rel_nt=3 * j, n0=int(n0), n1=int(n1), n2=int(n2),
                       F1_share=n1 / (n1 + n2) if n1 + n2 else np.nan, F0_share=n0 / (n0 + n1 + n2) if n0 + n1 + n2 else np.nan))
    # 5′ 端碱基上下文（codon_index ≥ 10，即 rel ≥ +30 的 CDS 内部）
    seqs = SEQ.loc[piv.index, 'seq_m30_p90'].values
    ctx = collections.defaultdict(lambda: np.zeros(2))
    for j in range(10, 29):
        c0 = 3 * j + 30
        if c0 + 2 > 120: break
        for gi, sq in enumerate(seqs):
            # F1 读段 P=3j+1 → 5′ 端 = 3j−11；F2 读段 P=3j+2 → 5′ 端 = 3j−10（同一密码子的 +2 位）
            b1 = sq[c0 + 1 - 12]; b2 = sq[c0 + 2 - 12]
            ctx[(b1, b2)] += [a[gi, c0 + 1], a[gi, c0 + 2]]
    for (b1, b2), v in ctx.items():
        if v.sum() > 0:
            po.append(dict(library=lab, codon_index=f'ctx_{b1}{b2}', rel_nt=np.nan, n0=np.nan, n1=int(v[0]), n2=int(v[1]),
                           F1_share=v[0] / v.sum(), F0_share=np.nan))
PO = pd.DataFrame(po); PO.to_csv(os.path.join(T, 'leakage_start_proximal_by_codon_and_context.tsv'), sep='\t', index=False, float_format='%.4g')
print(PO[PO.codon_index.astype(str).str.match(r'^-?\d+$') & PO.codon_index.astype(str).isin([str(x) for x in [-3, -1, 0, 1, 2, 3, 5, 8, 10, 15, 20, 25]])].pivot_table(
    index='codon_index', columns='library', values='F1_share').round(3).to_string())
print(PO[PO.codon_index.astype(str).str.startswith('ctx')].pivot_table(index='codon_index', columns='library', values='F1_share').round(3).to_string())
