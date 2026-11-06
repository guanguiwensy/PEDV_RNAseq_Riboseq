# -*- coding: utf-8 -*-
"""量化 v3_02 CDS 重叠判断笔误的影响：旧输出（tables/_archive_v1_buggy_overlap_20261026/）vs 修复后重跑（tables/）。
只比较宿主侧、已重跑的 v3_02 / v3_03 / g2_02 / g2_03 输出；涉及病毒窗口的 v3_04、g2_04、g3_01 未重跑（按用户指示不碰病毒候选）。"""
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
import os, gzip, re, collections, numpy as np, pandas as pd
import g2b_common as G
V3 = G.V3; T = os.path.join(V3, 'tables'); A = os.path.join(T, '_archive_v1_buggy_overlap_20261026')
OUT = os.path.join(G.OUT, 'G2c_freeze_checks'); rows = []
o = pd.read_csv(os.path.join(A, 'host_gene_level_29nt_frames.tsv.gz'), sep='\t'); n = pd.read_csv(os.path.join(T, 'host_gene_level_29nt_frames.tsv.gz'), sep='\t')
ot = set(o.loc[o.excl_overlap, 'transcript_id']); nt = set(n.loc[n.excl_overlap, 'transcript_id'])
# 全注释层（不限有 29 nt 计数）的重叠集合
import sys; sys.path.insert(0, os.path.join(V3, 'scripts'))
from v3_overlap import cds_overlap_set, segments_intersect_v1_buggy
cds_iv = collections.defaultdict(list); tchrom = {}
for line in gzip.open(os.path.join(V3, '..', '..', 'refdata_vero', 'Chlorocebus_sabaeus.ChlSab1.1.116.gtf.gz'), 'rt'):
    if line.startswith('#'): continue
    f = line.split('\t')
    if f[2] != 'CDS' or 'Ensembl_canonical' not in f[8] or 'transcript_biotype "protein_coding"' not in f[8]: continue
    t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1); cds_iv[t].append((int(f[3]), int(f[4]))); tchrom[t] = f[0]
new_all = cds_overlap_set(cds_iv, tchrom); old_all = cds_overlap_set(cds_iv, tchrom, pred=segments_intersect_v1_buggy)
summ = dict(annotation_transcripts=len(cds_iv), overlap_old_bugged=len(old_all), overlap_fixed=len(new_all),
            false_positive_old_only=len(old_all - new_all), missed_by_old_new_only=len(new_all - old_all), both=len(old_all & new_all))
pd.Series(summ).to_csv(os.path.join(OUT, 'v3_02_overlap_set_comparison.tsv'), sep='\t', header=['n'])
pd.DataFrame(dict(transcript_id=sorted(new_all - old_all))).to_csv(os.path.join(OUT, 'v3_02_true_overlaps_missed_by_old.tsv'), sep='\t', index=False)
# 训练集计数
co = pd.read_csv(os.path.join(A, 'leakage_training_set_counts_gene_level.tsv'), sep='\t').set_index('library')
cn = pd.read_csv(os.path.join(T, 'leakage_training_set_counts_gene_level.tsv'), sep='\t').set_index('library')
C = co.join(cn, lsuffix='_old', rsuffix='_new'); C.to_csv(os.path.join(OUT, 'v3_02_training_counts_old_vs_new.tsv'), sep='\t')
# 泄漏分箱：库合并 F1 份额、ρ 范围
bo = pd.read_csv(os.path.join(A, 'leakage_gene_level_bins.tsv'), sep='\t'); bn = pd.read_csv(os.path.join(T, 'leakage_gene_level_bins.tsv'), sep='\t')
B = []
for lib in sorted(bo.library.unique()):
    x, y = bo[bo.library == lib], bn[bn.library == lib]
    B.append(dict(library=lib, pooled_F1_old=x.library_pooled.iloc[0], pooled_F1_new=y.library_pooled.iloc[0],
                  rho_range_old=f'{x.rho_overdispersion.min():.3f}-{x.rho_overdispersion.max():.3f}',
                  rho_range_new=f'{y.rho_overdispersion.min():.3f}-{y.rho_overdispersion.max():.3f}',
                  dens_bin_F1_range_old=f'{x[x.covariate=="n0_density"].pooled_F1_share.min():.3f}-{x[x.covariate=="n0_density"].pooled_F1_share.max():.3f}',
                  dens_bin_F1_range_new=f'{y[y.covariate=="n0_density"].pooled_F1_share.min():.3f}-{y[y.covariate=="n0_density"].pooled_F1_share.max():.3f}'))
pd.DataFrame(B).to_csv(os.path.join(OUT, 'v3_02_leakage_bins_old_vs_new.tsv'), sep='\t', index=False)
# v3_03 基因层校准
k = ['library', 'test', 'subset', 'method']
so = pd.read_csv(os.path.join(A, 'calibration_body_gene_level_summary.tsv'), sep='\t'); sn = pd.read_csv(os.path.join(T, 'calibration_body_gene_level_summary.tsv'), sep='\t')
S = so.merge(sn, on=k, suffixes=('_old', '_new'))
S = S[['library', 'test', 'subset', 'method', 'n_old', 'n_new', 'frac_le_0.05_old', 'frac_le_0.05_new', 'frac_le_0.01_old', 'frac_le_0.01_new', 'ks_uniform_p_old', 'ks_uniform_p_new']]
S['d_le05'] = S['frac_le_0.05_new'] - S['frac_le_0.05_old']; S.to_csv(os.path.join(OUT, 'v3_03_calibration_old_vs_new.tsv'), sep='\t', index=False, float_format='%.5g')
# g2_02 模型比较
mo = pd.read_csv(os.path.join(A, 'g2_model_comparison_heldout.tsv'), sep='\t'); mn = pd.read_csv(os.path.join(T, 'g2_model_comparison_heldout.tsv'), sep='\t')
Mx = mo.merge(mn, on=['library', 'model'], suffixes=('_old', '_new'))[['library', 'model', 'test_codons_old', 'test_codons_new', 'pseudoR2_vs_M0_old', 'pseudoR2_vs_M0_new']]
Mx['d_R2'] = Mx.pseudoR2_vs_M0_new - Mx.pseudoR2_vs_M0_old; Mx.to_csv(os.path.join(OUT, 'g2_02_model_comparison_old_vs_new.tsv'), sep='\t', index=False, float_format='%.5g')
# g2_03 窗口校准
wo = pd.read_csv(os.path.join(A, 'g2_window_calibration_summary.tsv'), sep='\t'); wn = pd.read_csv(os.path.join(T, 'g2_window_calibration_summary.tsv'), sep='\t')
Wx = wo.merge(wn, on=['library', 'wtype', 'stratum', 'method'], suffixes=('_old', '_new'))[['library', 'wtype', 'stratum', 'method', 'n_old', 'n_new', 'f_0.05_old', 'f_0.05_new', 'f_0.01_old', 'f_0.01_new']]
Wx['d_f05'] = Wx['f_0.05_new'] - Wx['f_0.05_old']; Wx.to_csv(os.path.join(OUT, 'g2_03_window_calibration_old_vs_new.tsv'), sep='\t', index=False, float_format='%.5g')
pd.set_option('display.width', 250)
print(summ); print(C[['excl_overlap_old', 'excl_overlap_new', 'eligible_old', 'eligible_new', 'eligible_m_ge10_old', 'eligible_m_ge10_new']])
print(pd.DataFrame(B).to_string(index=False))
print('v3_03 |d le05| max', S.d_le05.abs().max().round(4), S.groupby('method').d_le05.agg(['min', 'max']).round(4))
print('g2_02 d_R2', Mx.d_R2.abs().max().round(4), Mx.pivot_table(index='library', columns='model', values='pseudoR2_vs_M0_new').round(3).to_string())
e = Wx[Wx.method == 'p_emp']; print('g2_03 p_emp d_f05', e.d_f05.abs().max().round(4)); print(e[e.stratum.str.contains('top|1%|99', regex=True)].sort_values('d_f05').round(4).to_string(index=False))
