# -*- coding: utf-8 -*-
"""Figure 6 联合散点图 v2：PRICE（CHX 29 nt，ORF 级）× HARTIS W1（HAR，起始相关占据），+1 框家族，每个家族一个点。
只改作图与展示逻辑；统计量全部读取已有输出，不重跑 PRICE / HARTIS，不改任何原始表。旧版 Fig6_HAR_PRICE_joint.* 保留不动。

输入（只读）：
  ../HARTIS_PRICE_family_summary.tsv            家族名、PRICE 最佳 P（3 个 CHX 库最小原始 P）及库/起点/q、PRICE 状态
  ../../G3_pedv_initiation/tables/HARTIS_v3_PEDV_start_level.tsv.gz   HARTIS 逐记录原始结果（起点 × 文库 × 读长）
  orf_calling_GIb/_price_final/annotated_29nt_named/<lib>.plus.orfs.annotated.tsv   PRICE 逐库原始输出（仅用于核对横坐标）
输出：
  ../Fig6_HAR_PRICE_joint_v2.pdf / .png
  ../Fig6_HAR_PRICE_joint_v2_plotdata.tsv       每个 +1 家族一行（含未绘制家族及原因）

展示规则（预先固定，不因个别家族破例）：
  1. 纵轴 HAR 记录选择：在该家族全部起点 × HAR 文库 × 合格读长的 W1 记录中，
     a) 若存在 status = in_support 且 p_emp 为有效数值（有限、0 < p ≤ 1）的记录 → 取 p_emp 最小的那一条（实心圆）；
     b) 否则若存在 status = capped_reference 且 p_cap 为有效数值的记录 → 取 p_cap 最小的那一条（空心三角，探索性）；
     c) 否则不绘制（不赋 P = 1），计入"无有效 HAR P"。
     并列最小 P 时按 (start_position, label, read_length) 升序取第一条（本数据无并列，脚本中核对）。
     两类 P 值均写入 plotdata 表，不删除。最小名义 P 含"多选一"选择效应，不是校准的家族级 P。
  2. 横轴：PRICE_best_p（3 个 CHX 库中最小原始 P，名义）；按 PRICE_best_lib 回到 PRICE 原始输出核对。
     PRICE 未调出 / 记录为经典读框（E-x1）的家族不绘制。
  3. 点大小：所选 HAR 记录本身的 W1 富集效应量 log2(O/E) = log2(n_cand / E1)，
     n_cand = 该起点 W1 中候选框（+1 = F1）读数，E1 = 按 256 上下文模型给出的期望 F1 读数（同一行）；n_cand = 0 时 log2 O/E = −inf，取最小点。
     选它而非 Z：Z（Z_cand / Z_cap）与 −log10 P 高度单调相关（重复编码纵轴），且 in_support 的 Z_cand 与 capped 的 Z_cap 是不同统计量；
     log2(O/E) 对两类记录定义相同，是效应大小而非证据强度。映射：size = 12 + 26·clip(log2 O/E, 0, 3)（pt²）；log2 O/E ≤ 0 取最小点。
  4. 不画任何连接不同数据点的线；只保留标签到其自身散点的短引导线，并检查引导线不经过其他散点。
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
import os, numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import spearmanr
import fig6_style as S
plt.rcParams.update({'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic', 'mathtext.bf': 'Arial:bold'})   # 数学文本也用 Arial（仅本图）

B = os.path.abspath(os.path.join(S.OUT, '..', '..'))
F = pd.read_csv(os.path.join(S.OUT, 'HARTIS_PRICE_family_summary.tsv'), sep='\t', dtype={'frame': str})
F = F[F.frame == '+1'].copy()
R = pd.read_csv(os.path.join(B, 'HARTIS_v3', 'G3_pedv_initiation', 'tables', 'HARTIS_v3_PEDV_start_level.tsv.gz'), sep='\t',
                dtype={'frame': str, 'read_length': str})
R = R[R['sid'].isin(SR.WHITELIST_LIBRARIES)].copy(); SR.assert_whitelisted(R['sid'].unique(), 'HARTIS start-level table')   # final sample table only; HARTIS fits are per library (independent)
R = R[(R.frame == '+1') & (R.drug == 'HAR') & (R.window == 'W1')].copy()
PL = {'CHX5': SR.lib("RIBO_VERO_CHX_1MOI_5h"), 'CHX8': SR.lib("RIBO_VERO_CHX_1MOI_8h"), 'CHX16': SR.lib("RIBO_VERO_CHX_1MOI_16h")}
PRAW = {k: pd.read_csv(os.path.join(B, 'orf_calling_GIb', '_price_final', 'annotated_29nt_named', f'{v}.plus.orfs.annotated.tsv'), sep='\t')
        for k, v in PL.items()}
LABEL_ALWAYS = ['N-x1', 'N-x2', 'S-x5', 'N-x3', 'S-x1', 'ORF1a-x19', 'ORF1b-x4', 'ORF1b-x14']
EFF = lambda e: 12 + 26 * np.clip(np.nan_to_num(e, nan=0.0), 0, 3)


def valid_p(s):
    s = pd.to_numeric(s, errors='coerce')
    return s.notna() & np.isfinite(s) & (s > 0) & (s <= 1)


def pick(x, pcol):
    x = x[valid_p(x[pcol])]
    if not len(x): return None, 0
    pmin = x[pcol].min(); t = x[x[pcol] == pmin].sort_values(['start_position', 'label', 'read_length'])
    return t.iloc[0], len(t)


rows = []
for f in F.itertuples():
    x = R[R.family_id == f.family_id]
    bi, ti = pick(x[x.status == 'in_support'], 'p_emp')
    bc, tc = pick(x[x.status == 'capped_reference'], 'p_cap')
    r = dict(family_id=f.family_id, family_name=f.family_name, parent_region=f.parent_region, stop_position=f.stop_position,
             n_start_variants=f.n_start_variants,
             PRICE_status=f.PRICE_status, PRICE_best_p=f.PRICE_best_p, PRICE_best_lib=f.PRICE_best_lib, PRICE_best_start=f.PRICE_best_start,
             PRICE_best_q_BH_within_file=f.PRICE_best_q_BH_within_file, PRICE_n_libs_p05=f.PRICE_n_libs_p05, PRICE_n_libs_called=f.PRICE_n_libs_called,
             PRICE_frame_mismatch=f.PRICE_frame_mismatch,
             HAR_in_support_best_p=np.nan, HAR_in_support_best_start=np.nan, HAR_in_support_best_lib=None, HAR_in_support_best_read_length=None,
             HAR_in_support_Z_cand=np.nan, HAR_in_support_n_tied_min=ti,
             HAR_capped_best_p=np.nan, HAR_capped_best_start=np.nan, HAR_capped_best_lib=None, HAR_capped_best_read_length=None,
             HAR_capped_Z_cap=np.nan, HAR_capped_n_tied_min=tc,
             HAR_n_libs_valid_in_support=x[(x.status == 'in_support') & valid_p(x.p_emp)].label.nunique(),
             HAR_n_libs_in_support_p05=x[(x.status == 'in_support') & valid_p(x.p_emp) & (x.p_emp <= 0.05)].label.nunique(),
             HAR_n_libs_capped_pcap05=x[(x.status == 'capped_reference') & valid_p(x.p_cap) & (x.p_cap <= 0.05)].label.nunique())
    if bi is not None:
        r.update(HAR_in_support_best_p=bi.p_emp, HAR_in_support_best_start=bi.start_position, HAR_in_support_best_lib=bi.label,
                 HAR_in_support_best_read_length=bi.read_length, HAR_in_support_Z_cand=bi.Z_cand)
    if bc is not None:
        r.update(HAR_capped_best_p=bc.p_cap, HAR_capped_best_start=bc.start_position, HAR_capped_best_lib=bc.label,
                 HAR_capped_best_read_length=bc.read_length, HAR_capped_Z_cap=bc.Z_cap)
    sel = bi if bi is not None else bc
    r['plot_HAR_source'] = 'in_support' if bi is not None else ('capped_reference' if bc is not None else 'none')
    if sel is not None:
        r.update(plot_HAR_p=sel.p_emp if bi is not None else sel.p_cap, plot_HAR_start=sel.start_position, plot_HAR_start_codon=sel.start_codon,
                 plot_HAR_lib=sel.label, plot_HAR_read_length=sel.read_length,
                 plot_HAR_Z=sel.Z_cand if bi is not None else sel.Z_cap, plot_HAR_Z_type='Z_cand' if bi is not None else 'Z_cap',
                 plot_HAR_n_cand=sel.n_cand, plot_HAR_E1=sel.E1, plot_HAR_m=sel.m,
                 plot_HAR_OE=sel.n_cand / sel.E1, plot_HAR_log2OE=np.log2(sel.n_cand / sel.E1) if sel.n_cand > 0 else -np.inf)   # n_cand = 0 → O/E = 0 → −inf（最小点）
    # 是否绘制及原因（不赋 P = 1）
    if f.PRICE_frame_mismatch: reason = 'not_plotted: PRICE record is canonical frame (not +1 evidence)'
    elif f.PRICE_status.startswith('not_called'): reason = 'not_plotted: PRICE not called (no PRICE P)'
    elif sel is None: reason = 'not_plotted: no valid HAR W1 P (in_support or capped)'
    else: reason = 'plotted'
    r['plot_status'] = reason
    rows.append(r)
D = pd.DataFrame(rows)
P = D[D.plot_status == 'plotted'].copy()
P['x'] = S.nlp(P.PRICE_best_p); P['y'] = S.nlp(P.plot_HAR_p); P['size_pt2'] = EFF(P.plot_HAR_log2OE)
D = D.merge(P[['family_id', 'x', 'y', 'size_pt2']], on='family_id', how='left')

# ---------- 核对 ----------
chk = []
chk.append(('每个家族一行', D.family_id.is_unique and len(D) == len(F)))
chk.append(('家族名全部来自 family_summary 且非空/非 unlabeled', D.family_name.notna().all() and not D.family_name.str.contains('unlabel', case=False).any()))
chk.append(('in_support 最小 P 无并列', int(D.HAR_in_support_n_tied_min.max()) <= 1))
chk.append(('capped 最小 P 无并列', int(D.HAR_capped_n_tied_min.max()) <= 1))
m = D.merge(F[['family_id', 'HAR_W1_best_p', 'HAR_W1_best_start', 'HAR_W1_Z_at_best', 'HAR_W1_capped_best_p', 'HAR_W1_capped_best_start']], on='family_id')
eq = lambda a, b: bool((((a - b).abs() < 1e-12) | (a.isna() & b.isna())).all())
chk.append(('in_support 最佳 P/起点/Z 与 family_summary 一致', eq(m.HAR_in_support_best_p, m.HAR_W1_best_p) and eq(m.HAR_in_support_best_start, m.HAR_W1_best_start)
            and eq(m.HAR_in_support_Z_cand, m.HAR_W1_Z_at_best)))
chk.append(('capped 最佳 P/起点与 family_summary 一致', eq(m.HAR_capped_best_p, m.HAR_W1_capped_best_p) and eq(m.HAR_capped_best_start, m.HAR_W1_capped_best_start)))
ok = True
for r in P.itertuples():   # 横坐标回到 PRICE 原始输出核对
    t = PRAW[r.PRICE_best_lib]; t = t[(t.our_start_1based == int(r.PRICE_best_start)) & (t.our_stop_end_1based == int(r.stop_position) + 2)]
    ok &= len(t) == 1 and abs(float(t['p value'].iloc[0]) - r.PRICE_best_p) < 1e-12
chk.append(('横坐标 P 与 PRICE 原始逐库输出一致（库/起点/终止）', ok))
ok = True
for r in P.itertuples():   # 纵坐标与点大小回到 HARTIS 原始逐记录核对：同一行
    t = R[(R.family_id == r.family_id) & (R.start_position == r.plot_HAR_start) & (R.label == r.plot_HAR_lib) & (R.read_length == r.plot_HAR_read_length)]
    pc, zc = ('p_emp', 'Z_cand') if r.plot_HAR_source == 'in_support' else ('p_cap', 'Z_cap')
    ok &= len(t) == 1 and t.status.iloc[0] == r.plot_HAR_source and abs(t[pc].iloc[0] - r.plot_HAR_p) < 1e-12 and abs(t[zc].iloc[0] - r.plot_HAR_Z) < 1e-9 \
        and abs(t.n_cand.iloc[0] / t.E1.iloc[0] - r.plot_HAR_OE) < 1e-12 and t.E1.iloc[0] > 0
chk.append(('纵坐标 P、Z、log2(O/E) 来自同一条 HARTIS 原始记录', ok))
chk.append(('in_support 有效时一律用 in_support', bool((P[P.HAR_in_support_best_p.notna()].plot_HAR_source == 'in_support').all())))
chk.append(('绘制点 P 均为 (0,1] 数值，未见 P = 1 填补', bool(((P.plot_HAR_p > 0) & (P.plot_HAR_p <= 1) & (P.PRICE_best_p > 0)).all())))
chk.append(('必标家族全部在绘图集内', set(LABEL_ALWAYS) <= set(P.family_name)))
rho_z = spearmanr(P.plot_HAR_Z, P.y).correlation; rho_e = spearmanr(P.plot_HAR_OE, P.y).correlation

# ---------- 作图 ----------
fig = plt.figure(figsize=(5.8, 4.3)); ax = fig.add_axes([0.10, 0.12, 0.58, 0.82])
xmax = max(P.x.max(), 2) * 1.12; ymax = max(P.y.max(), 2) * 1.15
ax.axvline(S.P05, color='0.6', lw=0.6, ls='--', zorder=0); ax.axhline(S.P05, color='0.6', lw=0.6, ls='--', zorder=0)
ax.text(S.P05 + 0.04, ymax * 0.995, 'P = 0.05', fontsize=5.5, color='0.45', va='top', rotation=90)
ax.text(xmax * 0.995, S.P05 + 0.03, 'P = 0.05', fontsize=5.5, color='0.45', ha='right', va='bottom')
col = P.parent_region.map(S.REGION_COLOR)
d = P.plot_HAR_source == 'in_support'
ax.scatter(P.x[d], P.y[d], s=P.size_pt2[d], c=col[d], edgecolors='white', linewidths=0.4, zorder=3)
ax.scatter(P.x[~d], P.y[~d], s=P.size_pt2[~d] * 1.3, marker='^', facecolors='none', edgecolors=col[~d], linewidths=0.8, zorder=3)
ax.set_xlim(0, xmax); ax.set_ylim(0, ymax)
ax.set_xlabel('−log$_{10}$ $P_{\\mathrm{PRICE}}$  (CHX, 29 nt; min. raw P over 3 libraries)')
ax.set_ylabel('−log$_{10}$ $P_{\\mathrm{HAR}}$  (HARTIS W1; min. nominal P)')

# 标签：必标 8 个 + 任一轴名义 P < 0.05 的家族（不另设特例）
lab = P[P.family_name.isin(LABEL_ALWAYS) | (P.PRICE_best_p < 0.05) | (P.plot_HAR_p <= 0.05)]
fig.canvas.draw(); T = ax.transData; dpi = fig.dpi; bb = ax.get_window_extent()
XY = T.transform(P[['x', 'y']].values); RAD = np.sqrt(P.size_pt2.values * np.where(d, 1, 1.3)) / 2 * dpi / 72   # 每个点的像素半径
idx = {fid: i for i, fid in enumerate(P.family_id)}
REF = T.transform((S.P05, S.P05))   # 参考线像素位置


def seg_dist(p, a, b):
    ab = b - a; t = np.clip(np.dot(p - a, ab) / max(np.dot(ab, ab), 1e-9), 0, 1); return np.hypot(*(a + t * ab - p))


def nearest_point(box):
    """离标签框（像素）边缘最近的散点下标（按点边缘计）。"""
    dx = np.maximum(np.maximum(box[0] - XY[:, 0], XY[:, 0] - box[2]), 0); dy = np.maximum(np.maximum(box[1] - XY[:, 1], XY[:, 1] - box[3]), 0)
    return int(np.argmin(np.hypot(dx, dy) - RAD))


def seg_hits_box(a, b, box, pad=1.5):
    """线段 a→b 是否穿过矩形 box（像素，四周外扩 pad）：在线段上密采样判断。"""
    t = np.linspace(0, 1, 60)[:, None]; q = a + t * (b - a)
    return bool(np.any((q[:, 0] > box[0] - pad) & (q[:, 0] < box[2] + pad) & (q[:, 1] > box[1] - pad) & (q[:, 1] < box[3] + pad)))


placed, leaders = [], []; n_leader = 0
dirs = [(0, 1), (1, 0), (-1, 0), (0.7, 0.7), (-0.7, 0.7), (0.7, -0.7), (-0.7, -0.7), (0, -1), (0.9, 0.4), (-0.9, 0.4), (0.9, -0.4), (-0.9, -0.4)]
fs = 6.0
for r in lab.sort_values('y', ascending=False).itertuples():
    i = idx[r.family_id]; px, py = XY[i]; w = 0.55 * fs * len(r.family_name) * dpi / 72; h = 1.3 * fs * dpi / 72
    best = None
    for strict in (True, False):   # 先要求引导线不跨 P = 0.05 参考线；找不到再放宽（参考线不是数据点）
        for rad in ((2, 4, 6, 10, 15, 21) if strict else (2, 4, 6, 10, 15, 21, 28, 36, 46, 58)):
            for dx, dy in dirs:
                off = RAD[i] + rad * dpi / 72
                cx, cy = px + dx * off, py + dy * off
                x0 = cx - w / 2 if abs(dx) < 0.3 else (cx if dx > 0 else cx - w); y0 = cy - h / 2 if abs(dy) < 0.3 else (cy if dy > 0 else cy - h)
                box = (x0, y0, x0 + w, y0 + h)
                if box[0] < bb.x0 + 2 or box[2] > bb.x1 - 2 or box[1] < bb.y0 + 2 or box[3] > bb.y1 - 2: continue
                if any(not (box[2] < b[0] - 1 or box[0] > b[2] + 1 or box[3] < b[1] - 1 or box[1] > b[3] + 1) for b in placed): continue
                if box[0] - 2 < REF[0] < box[2] + 2 or box[1] - 2 < REF[1] < box[3] + 2: continue   # 标签框不压 P = 0.05 参考线
                # 标签框不压任何散点
                if np.any((XY[:, 0] + RAD > box[0] - 1.5) & (XY[:, 0] - RAD < box[2] + 1.5) & (XY[:, 1] + RAD > box[1] - 1.5) & (XY[:, 1] - RAD < box[3] + 1.5)): continue
                # 引导线（框上离点最近处 → 点）不经过其他散点、不穿过已放标签
                a = np.array([np.clip(px, box[0], box[2]), np.clip(py, box[1], box[3])]); pnt = np.array([px, py])
                if any(seg_dist(XY[j], a, pnt) < RAD[j] + 2 for j in range(len(XY)) if j != i): continue
                if any(seg_hits_box(a, pnt, b) for b in placed): continue
                if nearest_point(box) != i: continue   # 标签框最近的散点必须是它自己的家族（避免指向歧义）
                if strict and ((a[0] - REF[0]) * (px - REF[0]) < 0 or (a[1] - REF[1]) * (py - REF[1]) < 0): continue
                best = (box, a); break
            if best: break
        if best: break
    if best is None: raise SystemExit(f'无法为 {r.family_name} 找到不压点/不穿点的标签位置')
    box, a = best; placed.append(box); leaders.append((r.family_name, i, a))
    tx, ty = T.inverted().transform(((box[0] + box[2]) / 2, (box[1] + box[3]) / 2))
    ax.text(tx, ty, r.family_name, fontsize=fs, ha='center', va='center', color='0.12', zorder=6)
    if np.hypot(*(a - np.array([px, py]))) > RAD[i] + 2.5:     # 只有标签离点较远时才画引导线，线止于该点边缘
        ex, ey = np.array([px, py]) + (a - np.array([px, py])) / np.hypot(*(a - np.array([px, py]))) * (RAD[i] + 0.8)
        (lx0, ly0), (lx1, ly1) = T.inverted().transform([a, (ex, ey)])
        ax.plot([lx0, lx1], [ly0, ly1], color='0.45', lw=0.4, zorder=5, solid_capstyle='butt'); n_leader += 1
chk.append(('每个标签框最近的散点是其自身家族', all(nearest_point(bx) == i for bx, (_, i, _) in zip(placed, leaders))))
chk.append(('被标注散点与其他散点不重叠（引导线终点无歧义）', all(np.sort(np.hypot(*(XY - XY[i]).T))[1] > RAD[i] for _, i, _ in leaders)))
chk.append(('图中线条 = 2 条 P=0.05 参考线 + 引导线（无点间连线）', len(ax.lines) == 2 + n_leader))

# 图例
reg = [Line2D([], [], marker='o', ls='', color=S.REGION_COLOR[k], ms=4, label=k) for k in S.REGION_COLOR if k in set(P.parent_region)]
shp = [Line2D([], [], marker='o', ls='', color='0.35', ms=4.5, label='in_support\n(EMP_Zs, formal)'),
       Line2D([], [], marker='^', ls='', mfc='none', mec='0.35', ms=5, label='capped_reference only\n(exploratory, not calibrated)')]
szs = [plt.scatter([], [], s=EFF(e), color='0.6', label=l) for e, l in ((0, '≤ 0'), (1, '1'), (2, '2'), (3, '≥ 3'))]
fig.legend(handles=reg, loc='upper left', bbox_to_anchor=(0.70, 0.95), title='Parent ORF', title_fontsize=6.5, alignment='left')
fig.legend(handles=shp, loc='upper left', bbox_to_anchor=(0.70, 0.61), title='HAR P source', title_fontsize=6.5, alignment='left')
fig.legend(handles=szs, loc='upper left', bbox_to_anchor=(0.70, 0.40), title='HAR W1 log$_2$(O/E)\n(same record as y)', title_fontsize=6.5,
           ncol=2, columnspacing=0.8, alignment='left')
nn = D.plot_status.value_counts()
fig.text(0.705, 0.13, f'+1 ORF families: {len(D)}; plotted: {len(P)}\n'
         f'not plotted (no P assigned):\n  PRICE not called: {nn.get("not_plotted: PRICE not called (no PRICE P)", 0)}\n'
         f'  no valid HAR P: {nn.get("not_plotted: no valid HAR W1 P (in_support or capped)", 0)}\n'
         f'  PRICE record canonical frame: {nn.get("not_plotted: PRICE record is canonical frame (not +1 evidence)", 0)}',
         fontsize=5.6, color='0.3', va='bottom')
fig.savefig(os.path.join(S.OUT, 'Fig6_HAR_PRICE_joint_v2.pdf')); fig.savefig(os.path.join(S.OUT, 'Fig6_HAR_PRICE_joint_v2.png'), dpi=300)

D['labelled'] = D.family_name.isin([n for n, _, _ in leaders])
D['p_value_note'] = ('x = min raw PRICE P over 3 CHX libraries; y = min nominal HARTIS W1 P over starts x HAR libraries x read lengths '
                     '(in_support preferred, capped_reference only if no valid in_support P); both nominal, not joint, not calibrated for selection; '
                     'size = log2(n_cand/E1) of the same HARTIS record')
D.to_csv(os.path.join(S.OUT, 'Fig6_HAR_PRICE_joint_v2_plotdata.tsv'), sep='\t', index=False)
for k, v in chk: print(('PASS ' if v else 'FAIL ') + k)
print(f'spearman(Z, -log10P_HAR) = {rho_z:.3f}; spearman(log2OE, -log10P_HAR) = {rho_e:.3f}; log2OE range {P.plot_HAR_log2OE.min():.2f}..{P.plot_HAR_log2OE.max():.2f}')
print(P[['family_name', 'plot_HAR_source', 'x', 'y', 'plot_HAR_Z', 'plot_HAR_log2OE', 'plot_HAR_n_cand', 'plot_HAR_E1']].sort_values('y').to_string())
print('labels:', [n for n, _, _ in leaders])
