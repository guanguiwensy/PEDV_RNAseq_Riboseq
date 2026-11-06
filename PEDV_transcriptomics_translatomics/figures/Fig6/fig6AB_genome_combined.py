# -*- coding: utf-8 -*-
"""Figure 6A+B 合并版：上 = HAR HARTIS W1（+1 框候选起点），中 = CHX 29 nt PRICE（+1 框家族），下 = PEDV 经典 ORF 轨道，三者共用基因组横轴。
简化展示（按用户要求）：不标注文库数（k/n）与 q；不画 E-x1（PRICE 记录为经典 M 读框）；HAR 点大小不编码 Z，只看 P。
统计量与 Fig6_HAR_genome / Fig6_PRICE_CHX_genome 完全相同（读同一张汇总表），旧版两图保留不动。

数据：../HARTIS_PRICE_start_summary.tsv（HAR，逐起点）、../HARTIS_PRICE_family_summary.tsv（PRICE 与标签，逐家族）。
HAR：y = −log10(该起点在全部 HAR 文库 × 合格读长中的最小名义 P)；实心圆 = in_support（EMP_Zs），空心三角 = capped_reference（探索性，未经病毒高深度校准），分开画、不合并。
PRICE：每个家族一个点（ORF 级检验，不复制到各起点）；x = 给出最小 P 的 CHX 库所报告的 ORF 起点；y = −log10(3 个 CHX 库中最小原始 P)。
无有效 P 的起点 / PRICE 未调出的家族不绘制（不赋 P = 1）。
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
import os, numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import fig6_style as S
plt.rcParams.update({'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic', 'mathtext.bf': 'Arial:bold'})

NAME = 'Fig6_HAR_PRICE_genome_combined'
D = pd.read_csv(os.path.join(S.OUT, 'HARTIS_PRICE_start_summary.tsv'), sep='\t', dtype={'frame': str}); D = D[D.frame == '+1']
F = pd.read_csv(os.path.join(S.OUT, 'HARTIS_PRICE_family_summary.tsv'), sep='\t', dtype={'frame': str}); F = F[F.frame == '+1']
ins = D[D.HAR_W1_min_nominal_p.notna()]; cap = D[D.HAR_W1_capped_min_p.notna()]
PR = F[F.PRICE_called.fillna(False) & F.PRICE_best_p.notna() & ~F.PRICE_frame_mismatch.astype(bool)].copy()   # 去掉 E-x1（经典读框）
assert (PR.PRICE_best_p > 0).all() and 'E-x1' not in set(PR.family_name)
PR['y'] = S.nlp(PR.PRICE_best_p)

fig = plt.figure(figsize=(7.4, 3.05))
L, W = 0.075, 0.905
axh = fig.add_axes([L, 0.60, W, 0.33]); axp = fig.add_axes([L, 0.27, W, 0.29], sharex=axh); axg = fig.add_axes([L, 0.115, W, 0.10], sharex=axh)

# ---------- 上：HAR ----------
axh.axhline(S.P05, color='0.6', lw=0.5, ls='--', zorder=0); axh.text(250, S.P05 + 0.05, 'P = 0.05', fontsize=5, color='0.45')
axh.scatter(ins.start_position, S.nlp(ins.HAR_W1_min_nominal_p), s=9, c=ins.parent.map(S.REGION_COLOR), edgecolors='white', linewidths=0.25, zorder=3)
axh.scatter(cap.start_position, S.nlp(cap.HAR_W1_capped_min_p), s=10, marker='^', facecolors='none', edgecolors=cap.parent.map(S.REGION_COLOR), linewidths=0.55, zorder=2)
yh = max(S.nlp(ins.HAR_W1_min_nominal_p).max(), S.nlp(cap.HAR_W1_capped_min_p).max()) * 1.25
axh.set_ylim(0, yh); axh.set_ylabel('−log$_{10}$ $P_{\\mathrm{HAR}}$\n(HARTIS W1)', fontsize=6.5); axh.tick_params(labelbottom=False, labelsize=6)
# 标签（只写家族名）：N-x1、N-x2、S-x5 + HAR in_support P ≤ 0.05 的家族；位置规则同 Fig6_HAR_genome（in_support 名义支持 → in_support 最佳点；
#   否则 capped 在 ≥ 2 个库 p ≤ 0.05 → capped 最佳点；否则 in_support 最佳点）
pts, txt = [], []
for r in F[F.family_name.isin(['N-x1', 'N-x2', 'S-x5']) | (F.HAR_W1_best_p <= 0.05)].itertuples():
    ncap = int(r.HAR_W1_capped_n_libs_pcap05) if pd.notna(r.HAR_W1_capped_best_p) else 0
    if (pd.notna(r.HAR_W1_best_p) and r.HAR_W1_best_p <= 0.05) or (pd.notna(r.HAR_W1_best_p) and ncap < 2):
        pts.append((r.HAR_W1_best_start, S.nlp(r.HAR_W1_best_p)))
    elif pd.notna(r.HAR_W1_capped_best_p):
        pts.append((r.HAR_W1_capped_best_start, S.nlp(r.HAR_W1_capped_best_p)))
    else:
        continue
    txt.append(r.family_name)
allh = np.r_[np.c_[ins.start_position, S.nlp(ins.HAR_W1_min_nominal_p)], np.c_[cap.start_position, S.nlp(cap.HAR_W1_capped_min_p)]]
S.place_labels(axh, pts, txt, allh, fontsize=5.4, r_pt=(6, 10, 15, 21, 28, 36))
labs_h = txt

# ---------- 中：PRICE ----------
axp.axhline(S.P05, color='0.6', lw=0.5, ls='--', zorder=0); axp.text(250, S.P05 + 0.08, 'P = 0.05', fontsize=5, color='0.45')
axp.scatter(PR.PRICE_best_start, PR.y, s=11, c=PR.parent_region.map(S.REGION_COLOR), edgecolors='white', linewidths=0.25, zorder=3)
axp.set_ylim(0, PR.y.max() * 1.3); axp.set_ylabel('−log$_{10}$ $P_{\\mathrm{PRICE}}$\n(CHX, 29 nt)', fontsize=6.5)
axp.tick_params(labelbottom=False, labelsize=6)
lab = PR[PR.family_name.isin(['N-x1', 'N-x2', 'S-x5']) | (PR.PRICE_best_p < 0.05)]
S.place_labels(axp, list(zip(lab.PRICE_best_start, lab.y)), list(lab.family_name), PR[['PRICE_best_start', 'y']].values, fontsize=5.4,
               r_pt=(6, 10, 15, 21, 28, 36))
labs_p = list(lab.family_name)

# ---------- 下：基因组轨道 ----------
S.genome_track(axg); axg.tick_params(labelsize=6); axg.xaxis.label.set_size(7)

# 图例（顶部一行）
reg = [Line2D([], [], marker='o', ls='', color=S.REGION_COLOR[k], ms=3.5, label=k) for k in S.REGION_COLOR]
shp = [Line2D([], [], marker='o', ls='', color='0.35', ms=3.5, label='HAR in_support'),
       Line2D([], [], marker='^', ls='', mfc='none', mec='0.35', ms=4, label='HAR capped_reference (exploratory)')]
fig.legend(handles=reg + shp, loc='upper center', bbox_to_anchor=(L + W / 2, 1.0), ncol=9, fontsize=5.8, handletextpad=0.2, columnspacing=0.9)
fig.savefig(os.path.join(S.OUT, f'{NAME}.pdf')); fig.savefig(os.path.join(S.OUT, f'{NAME}.png'), dpi=300)
print('HAR in_support starts', len(ins), 'capped starts', len(cap), 'labels', labs_h)
print('PRICE families plotted', len(PR), 'labels', labs_p)
