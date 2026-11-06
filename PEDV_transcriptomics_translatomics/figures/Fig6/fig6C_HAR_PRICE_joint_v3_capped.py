# -*- coding: utf-8 -*-
"""Figure 6 联合散点图 v3（capped 纵轴）：PRICE（CHX 29 nt，ORF 级）× HARTIS HAR W1 capped_reference，+1 框家族，每个家族一个点。
只改作图口径；统计量全部读取已有输出，不重跑 PRICE / HARTIS，不改任何原始表。旧版 v1（Fig6_HAR_PRICE_joint.*）与 v2（*_v2.*）原样保留。

输入（只读）：
  ../HARTIS_PRICE_family_summary.tsv            家族名、父 ORF、PRICE 最佳 P（3 个 CHX 库最小原始 P）及库/起点/q、经典读框标记
  ../../G3_pedv_initiation/tables/HARTIS_v3_PEDV_start_level.tsv.gz   HARTIS 逐记录原始结果（起点 × 文库 × 读长）
  orf_calling_GIb/_price_final/annotated_29nt_named/<lib>.plus.orfs.annotated.tsv   PRICE 逐库原始输出（仅用于核对横坐标）
输出：
  ../Fig6_HAR_PRICE_joint_v3_capped.pdf / .png
  ../Fig6_HAR_PRICE_joint_v3_capped_plotdata.tsv   每个 +1 家族一行（含未绘制家族及原因；in_support p_emp 仅作对照列保留）

展示规则（预先固定）：
  1. 只用 +1 框家族；排除经典读框 ORF：candidate_frame_equals_canonical 非空（候选读框即某经典 ORF 的读框）或 PRICE 记录为经典读框
     （本数据仅 E-x1：其 +1 框即 M 读框，PRICE 记录为经典 M ORF）。
  2. 纵轴：HAR、W1、status = capped_reference、p_cap 为有效数值（有限、0 < p ≤ 1）的全部记录（起点 × HAR 文库 × 合格读长）中取 p_cap 最小的一条，
     绘制 −log10(p_cap)。不使用 in_support 的 p_emp（即使该家族有 in_support 记录）。并列最小 P 时按 (start, library, read length) 升序取第一条。
     无 capped_reference 记录的家族不绘制、不赋 P = 1，在输出表与图中计数。
  3. 横轴：PRICE_best_p（3 个 CHX 库中最小原始 P）；按 PRICE_best_lib 回到 PRICE 原始输出核对。PRICE 未调出的家族不绘制。
  4. 点大小：所选 capped 记录本身的 W1 富集效应量 log2(O/E) = log2(n_cand / E1)（同一行；n_cand = 候选框 F1 读数，E1 = 256 上下文模型期望 F1 读数）；
     size = 12 + 26·clip(log2 O/E, 0, 3) pt²，≤ 0（含 n_cand = 0 → −inf）取最小点。选 O/E 而非 Z_cap：与 v2 同一定义便于对照，且为效应大小而非证据强度；
     Z_cap 保留在表中。两者与 −log10 p_cap 的 Spearman 相关运行时打印（本数据 Z_cap 0.82、O/E 0.85，点大小与纵轴部分冗余）。
  5. 不画任何连接不同数据点的线；只保留标签到其自身散点的短引导线，并检查引导线不经过其他散点、标签最近的点就是自身。
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

NAME = 'Fig6_HAR_PRICE_joint_v3_capped'
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
LABEL_FOCUS = ['N-x1', 'N-x2', 'S-x1', 'S-x5', 'ORF1b-x13', 'ORF1b-x14', 'ORF1a-x19', 'M-x1']
EFF = lambda e: 12 + 26 * np.clip(np.nan_to_num(e, nan=0.0, neginf=-1.0), 0, 3)


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
    bc, tc = pick(x[x.status == 'capped_reference'], 'p_cap')
    bi, _ = pick(x[x.status == 'in_support'], 'p_emp')
    canonical = pd.notna(f.candidate_frame_equals_canonical) or bool(f.PRICE_frame_mismatch)
    r = dict(family_id=f.family_id, family_name=f.family_name, parent_region=f.parent_region, stop_position=f.stop_position,
             candidate_frame_equals_canonical=f.candidate_frame_equals_canonical, excluded_canonical_frame=canonical,
             PRICE_status=f.PRICE_status, PRICE_best_p=f.PRICE_best_p, PRICE_best_lib=f.PRICE_best_lib, PRICE_best_start=f.PRICE_best_start,
             PRICE_best_q_BH_within_file=f.PRICE_best_q_BH_within_file, PRICE_n_libs_p05=f.PRICE_n_libs_p05, PRICE_n_libs_called=f.PRICE_n_libs_called,
             HAR_capped_n_records=int((x.status == 'capped_reference').sum()),
             HAR_capped_n_libs=x[(x.status == 'capped_reference') & valid_p(x.p_cap)].label.nunique(),
             HAR_capped_n_libs_pcap05=x[(x.status == 'capped_reference') & valid_p(x.p_cap) & (x.p_cap <= 0.05)].label.nunique(),
             HAR_capped_n_tied_min=tc,
             HAR_in_support_best_p_reference_only=bi.p_emp if bi is not None else np.nan,
             HAR_in_support_n_libs_valid_reference_only=x[(x.status == 'in_support') & valid_p(x.p_emp)].label.nunique())
    if bc is not None:
        r.update(plot_HAR_status='capped_reference', plot_HAR_p_cap=bc.p_cap, plot_HAR_Z_cap=bc.Z_cap, plot_HAR_start=bc.start_position,
                 plot_HAR_start_codon=bc.start_codon, plot_HAR_lib=bc.label, plot_HAR_sid=bc.sid, plot_HAR_read_length=bc.read_length,
                 plot_HAR_K_ref=bc.K_ref, plot_HAR_ref_bin=bc.ref_bin, plot_HAR_m=bc.m, plot_HAR_n_cand=bc.n_cand, plot_HAR_E1=bc.E1,
                 plot_HAR_OE=bc.n_cand / bc.E1, plot_HAR_log2OE=np.log2(bc.n_cand / bc.E1) if bc.n_cand > 0 else -np.inf)
    if canonical: reason = 'not_plotted: canonical reading frame (excluded)'
    elif f.PRICE_status.startswith('not_called'): reason = 'not_plotted: PRICE not called (no PRICE P)'
    elif bc is None: reason = 'not_plotted: no HAR W1 capped_reference record'
    else: reason = 'plotted'
    r['plot_status'] = reason
    rows.append(r)
D = pd.DataFrame(rows)
P = D[D.plot_status == 'plotted'].copy()
P['x'] = S.nlp(P.PRICE_best_p); P['y'] = S.nlp(P.plot_HAR_p_cap); P['size_pt2'] = EFF(P.plot_HAR_log2OE)
D = D.merge(P[['family_id', 'x', 'y', 'size_pt2']], on='family_id', how='left')

# ---------- 核对 ----------
chk = []
chk.append(('每个家族一行', D.family_id.is_unique and len(D) == len(F)))
chk.append(('家族名全部来自 family_summary 且非空/非 unlabeled', D.family_name.notna().all() and not D.family_name.str.contains('unlabel', case=False).any()))
chk.append(('绘图集中无经典读框家族', not P.excluded_canonical_frame.any() and 'E-x1' not in set(P.family_name)))
chk.append(('capped 最小 P 无并列', int(D.HAR_capped_n_tied_min.max()) <= 1))
m = D.merge(F[['family_id', 'HAR_W1_capped_best_p', 'HAR_W1_capped_best_start', 'HAR_W1_capped_Zcap']], on='family_id')
eq = lambda a, b: bool((((a - b).abs() < 1e-9) | (a.isna() & b.isna())).all())
chk.append(('capped 最佳 p_cap/起点/Z_cap 与 family_summary 一致', eq(m.plot_HAR_p_cap, m.HAR_W1_capped_best_p) and eq(m.plot_HAR_start, m.HAR_W1_capped_best_start)
            and eq(m.plot_HAR_Z_cap, m.HAR_W1_capped_Zcap)))
ok = True
for r in P.itertuples():   # 横坐标回到 PRICE 原始输出核对
    t = PRAW[r.PRICE_best_lib]; t = t[(t.our_start_1based == int(r.PRICE_best_start)) & (t.our_stop_end_1based == int(r.stop_position) + 2)]
    ok &= len(t) == 1 and abs(float(t['p value'].iloc[0]) - r.PRICE_best_p) < 1e-12
chk.append(('横坐标 P 与 PRICE 原始逐库输出一致（库/起点/终止）', ok))
ok = True
for r in P.itertuples():   # 纵坐标、Z_cap、O/E 回到 HARTIS 原始逐记录：唯一一行、status = capped_reference
    t = R[(R.family_id == r.family_id) & (R.start_position == r.plot_HAR_start) & (R.label == r.plot_HAR_lib) & (R.read_length == r.plot_HAR_read_length)]
    ok &= len(t) == 1 and t.status.iloc[0] == 'capped_reference' and abs(t.p_cap.iloc[0] - r.plot_HAR_p_cap) < 1e-12 \
        and abs(t.Z_cap.iloc[0] - r.plot_HAR_Z_cap) < 1e-9 and t.E1.iloc[0] > 0 and abs(t.n_cand.iloc[0] / t.E1.iloc[0] - r.plot_HAR_OE) < 1e-12
    ok &= abs(r.y - S.nlp(r.plot_HAR_p_cap)) < 1e-12
chk.append(('纵坐标 = −log10(p_cap)，p_cap、Z_cap、O/E 来自同一条 capped_reference 原始记录', ok))
nn = P.set_index('family_name')
ok = all(abs(nn.loc[k, 'y'] - S.nlp(nn.loc[k, 'plot_HAR_p_cap'])) < 1e-12 and abs(nn.loc[k, 'plot_HAR_p_cap'] - nn.loc[k, 'HAR_in_support_best_p_reference_only']) > 1e-6
         for k in ('N-x1', 'N-x2'))
chk.append(('N-x1/N-x2 纵坐标为 capped p_cap，而非 in_support p_emp', ok))
chk.append(('绘制点 P 均为 (0,1] 数值，未见 P = 1 填补', bool(((P.plot_HAR_p_cap > 0) & (P.plot_HAR_p_cap <= 1) & (P.PRICE_best_p > 0)).all())))
rho_z = spearmanr(P.plot_HAR_Z_cap, P.y).correlation; rho_e = spearmanr(P.plot_HAR_OE, P.y).correlation

# ---------- 作图 ----------
fig = plt.figure(figsize=(5.8, 4.3)); ax = fig.add_axes([0.10, 0.12, 0.58, 0.82])
xmax = max(P.x.max(), 2) * 1.12; ymax = max(P.y.max(), 2) * 1.12
ax.axvline(S.P05, color='0.6', lw=0.6, ls='--', zorder=0); ax.axhline(S.P05, color='0.6', lw=0.6, ls='--', zorder=0)
ax.text(S.P05 + 0.04, ymax * 0.995, 'P = 0.05', fontsize=5.5, color='0.45', va='top', rotation=90)
ax.text(xmax * 0.995, S.P05 + 0.03, 'P = 0.05', fontsize=5.5, color='0.45', ha='right', va='bottom')
col = P.parent_region.map(S.REGION_COLOR)
ax.scatter(P.x, P.y, s=P.size_pt2, c=col, edgecolors='white', linewidths=0.4, zorder=3)
ax.set_xlim(-0.03 * xmax, xmax); ax.set_ylim(-0.03 * ymax, ymax)   # 轻微负边距，避免贴轴的点被裁切
ax.set_xlabel('−log$_{10}$ $P_{\\mathrm{PRICE}}$  (CHX, 29 nt; min. raw P over 3 libraries)')
ax.set_ylabel('−log$_{10}$ $p_{\\mathrm{cap}}$  (HAR W1 capped_reference; min. over family)')

# 标签：重点 8 个（在绘图集内者）+ 任一轴名义 P < 0.05 的家族
lab = P[P.family_name.isin(LABEL_FOCUS) | (P.PRICE_best_p < 0.05) | (P.plot_HAR_p_cap <= 0.05)]
fig.canvas.draw(); T = ax.transData; dpi = fig.dpi; bb = ax.get_window_extent()
XY = T.transform(P[['x', 'y']].values); RAD = np.sqrt(P.size_pt2.values) / 2 * dpi / 72
idx = {fid: i for i, fid in enumerate(P.family_id)}
REF = T.transform((S.P05, S.P05))


def seg_dist(p, a, b):
    ab = b - a; t = np.clip(np.dot(p - a, ab) / max(np.dot(ab, ab), 1e-9), 0, 1); return np.hypot(*(a + t * ab - p))


def nearest_point(box):
    """离标签框（像素）边缘最近的散点下标（按点边缘计）。"""
    dx = np.maximum(np.maximum(box[0] - XY[:, 0], XY[:, 0] - box[2]), 0); dy = np.maximum(np.maximum(box[1] - XY[:, 1], XY[:, 1] - box[3]), 0)
    return int(np.argmin(np.hypot(dx, dy) - RAD))


def seg_hits_box(a, b, box, pad=1.5):
    t = np.linspace(0, 1, 60)[:, None]; q = a + t * (b - a)
    return bool(np.any((q[:, 0] > box[0] - pad) & (q[:, 0] < box[2] + pad) & (q[:, 1] > box[1] - pad) & (q[:, 1] < box[3] + pad)))


placed, leaders = [], []; n_leader = 0
dirs = [(1, 0), (0, 1), (-1, 0), (0.7, 0.7), (-0.7, 0.7), (0.7, -0.7), (-0.7, -0.7), (0, -1), (0.9, 0.4), (-0.9, 0.4), (0.9, -0.4), (-0.9, -0.4)]
fs = 6.0
lab = lab.assign(focus=lab.family_name.isin(LABEL_FOCUS)).sort_values(['focus', 'y'], ascending=[False, False])   # 重点家族先放
skipped = []
for r in lab.itertuples():
    i = idx[r.family_id]; px, py = XY[i]; w = 0.55 * fs * len(r.family_name) * dpi / 72; h = 1.3 * fs * dpi / 72
    best = None
    for lvl in (0, 1, 2):   # 0：标签紧贴、引导线不跨参考线；1：放宽距离；2：再允许标签框压参考线（参考线不是数据点）
        strict = lvl == 0
        for rad in ((1.5, 3, 5, 8, 12) if strict else (1.5, 3, 5, 8, 12, 17, 23, 30, 40, 52)):
            for dx, dy in dirs:
                off = RAD[i] + rad * dpi / 72
                cx, cy = px + dx * off, py + dy * off
                x0 = cx - w / 2 if abs(dx) < 0.3 else (cx if dx > 0 else cx - w); y0 = cy - h / 2 if abs(dy) < 0.3 else (cy if dy > 0 else cy - h)
                box = (x0, y0, x0 + w, y0 + h)
                if box[0] < bb.x0 + 2 or box[2] > bb.x1 - 2 or box[1] < bb.y0 + 2 or box[3] > bb.y1 - 2: continue
                if any(not (box[2] < b[0] - 1 or box[0] > b[2] + 1 or box[3] < b[1] - 1 or box[1] > b[3] + 1) for b in placed): continue
                if lvl < 2 and (box[0] - 2 < REF[0] < box[2] + 2 or box[1] - 2 < REF[1] < box[3] + 2): continue
                if np.any((XY[:, 0] + RAD > box[0] - 1.5) & (XY[:, 0] - RAD < box[2] + 1.5) & (XY[:, 1] + RAD > box[1] - 1.5) & (XY[:, 1] - RAD < box[3] + 1.5)): continue
                a = np.array([np.clip(px, box[0], box[2]), np.clip(py, box[1], box[3])]); pnt = np.array([px, py])
                if np.hypot(*(a - pnt)) > RAD[i] + 2.5:   # 只检查实际画出的引导线段（标签框 → 自身点边缘）
                    e = pnt + (a - pnt) / np.hypot(*(a - pnt)) * (RAD[i] + 0.8)
                    if any(seg_dist(XY[j], a, e) < RAD[j] + 2 for j in range(len(XY)) if j != i): continue
                    if any(seg_hits_box(a, e, b) for b in placed): continue
                nb = nearest_point(box)
                if nb != i:
                    # 坐标几乎重合的点对（两点圆面重叠）：允许标签放在远离邻点的一侧（仅放宽轮），引导线止于自身点远侧边缘
                    ov = np.hypot(*(XY[nb] - pnt)) < RAD[i] + RAD[nb]
                    cvec = np.array([(box[0] + box[2]) / 2, (box[1] + box[3]) / 2]) - pnt
                    if strict or not ov or np.dot(cvec, XY[nb] - pnt) >= 0: continue
                if strict and ((a[0] - REF[0]) * (px - REF[0]) < 0 or (a[1] - REF[1]) * (py - REF[1]) < 0): continue
                best = (box, a); break
            if best: break
        if best: break
    if best is None:
        if r.focus: raise SystemExit(f'无法为重点家族 {r.family_name} 找到不压点/不穿点的标签位置')
        skipped.append(r.family_name); continue   # 非重点的附加标签：放不下就不标（不强行压点），记录在表中
    box, a = best; placed.append(box); leaders.append((r.family_name, i, a))
    tx, ty = T.inverted().transform(((box[0] + box[2]) / 2, (box[1] + box[3]) / 2))
    ax.text(tx, ty, r.family_name, fontsize=fs, ha='center', va='center', color='0.12', zorder=6)
    if np.hypot(*(a - np.array([px, py]))) > RAD[i] + 2.5:     # 标签离点较远时才画引导线，线止于该点边缘
        ex, ey = np.array([px, py]) + (a - np.array([px, py])) / np.hypot(*(a - np.array([px, py]))) * (RAD[i] + 0.8)
        (lx0, ly0), (lx1, ly1) = T.inverted().transform([a, (ex, ey)])
        ax.plot([lx0, lx1], [ly0, ly1], color='0.45', lw=0.4, zorder=5, solid_capstyle='butt'); n_leader += 1
OVL = [(P.family_name.iloc[i], P.family_name.iloc[j]) for i in range(len(XY)) for j in range(i + 1, len(XY)) if np.hypot(*(XY[i] - XY[j])) < RAD[i] + RAD[j]]
amb = [n for bx, (n, i, _) in zip(placed, leaders) if nearest_point(bx) != i]
chk.append(('标签框最近散点为自身；例外仅限坐标几乎重合的点对且标签在远离邻点一侧', all(any(n in pr for pr in OVL) for n in amb)))
chk.append(('图中线条 = 2 条 P=0.05 参考线 + 引导线（无点间连线）', len(ax.lines) == 2 + n_leader))
focus_in = [k for k in LABEL_FOCUS if k in set(P.family_name)]
chk.append(('绘图集内的重点家族全部已标注', set(focus_in) <= {n for n, _, _ in leaders}))

# 图例与计数
reg = [Line2D([], [], marker='o', ls='', color=S.REGION_COLOR[k], ms=4, label=k) for k in S.REGION_COLOR if k in set(P.parent_region)]
szs = [plt.scatter([], [], s=EFF(e), color='0.6', label=l) for e, l in ((0, '≤ 0'), (1, '1'), (2, '2'), (3, '≥ 3'))]
fig.legend(handles=reg, loc='upper left', bbox_to_anchor=(0.70, 0.95), title='Parent ORF', title_fontsize=6.5, alignment='left')
fig.legend(handles=szs, loc='upper left', bbox_to_anchor=(0.70, 0.62), title='HAR W1 log$_2$(O/E)\n(same record as y)', title_fontsize=6.5,
           ncol=2, columnspacing=0.8, alignment='left')
nn = D.plot_status.value_counts(); miss_focus = [k for k in LABEL_FOCUS if k not in set(P.family_name)]
note = (f'y: HAR capped_reference p_cap only\n(exploratory; not calibrated for\nhigh viral depth). Not a joint P.\n\n'
        f'+1 ORF families: {len(D)}; plotted: {len(P)}\nnot plotted (no P assigned):\n'
        f'  no capped_reference record: {nn.get("not_plotted: no HAR W1 capped_reference record", 0)}\n'
        f'  PRICE not called: {nn.get("not_plotted: PRICE not called (no PRICE P)", 0)}\n'
        f'  canonical reading frame: {nn.get("not_plotted: canonical reading frame (excluded)", 0)}')
if OVL:
    note += '\n' + '\n'.join(f'  {u} and {v} nearly coincide;\n  labelled on opposite sides' for u, v in OVL)
if miss_focus:
    st = D.set_index('family_name').plot_status
    note += '\n' + '\n'.join(f'  {k} (focus): no capped record' if 'capped' in st[k] else f'  {k} (focus): {st[k].replace("not_plotted: ", "")}' for k in miss_focus)
fig.text(0.705, 0.11, note, fontsize=5.6, color='0.3', va='bottom')
fig.savefig(os.path.join(S.OUT, f'{NAME}.pdf')); fig.savefig(os.path.join(S.OUT, f'{NAME}.png'), dpi=300)

D['labelled'] = D.family_name.isin([n for n, _, _ in leaders]); D['label_skipped_no_clear_position'] = D.family_name.isin(skipped)
D['p_value_note'] = ('x = min raw PRICE P over 3 CHX libraries (29 nt); y = min HARTIS W1 capped_reference p_cap over starts x HAR libraries x read lengths '
                     '(exploratory, not calibrated for high viral depth; in_support p_emp NOT used, kept only as reference column); '
                     'both nominal, not joint, not calibrated for family-level selection; size = log2(n_cand/E1) of the same capped record')
D.to_csv(os.path.join(S.OUT, f'{NAME}_plotdata.tsv'), sep='\t', index=False)
for k, v in chk: print(('PASS ' if v else 'FAIL ') + k)
print(f'spearman(Z_cap, -log10 p_cap) = {rho_z:.3f}; spearman(O/E, -log10 p_cap) = {rho_e:.3f}')
print(D.plot_status.value_counts().to_dict()); print('focus not plotted:', {k: D.set_index('family_name').loc[k, 'plot_status'] for k in miss_focus})
print(P[['family_name', 'x', 'y', 'plot_HAR_p_cap', 'HAR_in_support_best_p_reference_only', 'plot_HAR_Z_cap', 'plot_HAR_log2OE', 'plot_HAR_start', 'plot_HAR_lib',
         'plot_HAR_read_length', 'HAR_capped_n_libs_pcap05']].sort_values('y').to_string())
print('labels:', [n for n, _, _ in leaders], 'skipped extra labels:', skipped)
print('overlapping point pairs (pixel discs overlap):', OVL, '; labels whose nearest point is the overlapping neighbour:', amb)
