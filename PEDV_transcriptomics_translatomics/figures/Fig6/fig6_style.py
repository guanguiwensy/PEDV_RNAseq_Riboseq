# -*- coding: utf-8 -*-
"""Figure 6 共用样式：区域配色、字体、PEDV 经典 ORF 示意轨道、标签避让、P 值下限处理。"""
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
import os, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
plt.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
                     'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'axes.spines.top': False, 'axes.spines.right': False,
                     'pdf.fonttype': 42, 'ps.fonttype': 42, 'legend.frameon': False, 'legend.fontsize': 6.5})
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.abspath(os.path.join(HERE, '..'))
GENOME_LEN = 27959
CANON = [('ORF1a', 293, 12622), ('ORF1b', 12592, 20613), ('S', 20610, 24758), ('ORF3', 24758, 25033), ('E', 25361, 25600), ('M', 25608, 26288), ('N', 26300, 27625)]
REGION_COLOR = {'ORF1a': '#2a78d6', 'ORF1b': '#1baf7a', 'S': '#eb6834', 'ORF3': '#eda100', 'E': '#e87ba4', 'M': '#6250d6', 'N': '#e34948'}
P05 = -np.log10(0.05)


def nlp(p):
    return -np.log10(p)


def genome_track(ax):
    """经典 ORF 示意：坐标为实际基因组位置；ORF1a/ORF1b 重叠与 S/ORF3 相接时分两行画。"""
    rows = {'ORF1a': 1, 'ORF1b': 0, 'S': 1, 'ORF3': 0, 'E': 1, 'M': 0, 'N': 1}
    for nm, a, b in CANON:
        y = rows[nm] * 0.55
        ax.add_patch(Rectangle((a, y), b - a + 1, 0.42, facecolor=REGION_COLOR[nm], edgecolor='none', alpha=0.9))
        w = b - a + 1
        if w > 900:      # 长 ORF：名字写在框内
            ax.text((a + b) / 2, y + 0.21, nm, ha='center', va='center', fontsize=5.5, color='white', fontweight='bold', clip_on=False)
        else:            # 短 ORF：名字写在框外（上行写在上方，下行写在下方）
            ax.text((a + b) / 2, y + 0.55 if rows[nm] else -0.18, nm, ha='center', va='center', fontsize=5, color=REGION_COLOR[nm], clip_on=False)
    ax.set_ylim(-0.35, 1.15); ax.set_yticks([]); ax.spines['left'].set_visible(False)
    ax.set_xlim(0, GENOME_LEN); ax.set_xlabel('PEDV genome position (nt)')


def repel_labels(ax, pts, texts, x_span, y_span, min_dx=0.06, min_dy=0.07, start_dy=0.06, color='0.15', fontsize=6, x_lim=None):
    """简单的贪心避让：标签按 x 排序，依次放在点上方，若与已放标签重叠则继续上移；用引导线连回原点。坐标单位为数据坐标。"""
    placed = []
    order = np.argsort([p[0] for p in pts])
    for i in order:
        x, y = pts[i]; ty = y + start_dy * y_span; tx = x
        if x_lim is not None: tx = min(max(tx, x_lim[0]), x_lim[1])          # 标签不越出坐标轴
        for _ in range(60):
            clash = [q for q in placed if abs(q[0] - tx) < min_dx * x_span and abs(q[1] - ty) < min_dy * y_span]
            if not clash: break
            ty += 0.5 * min_dy * y_span
        placed.append((tx, ty))
        ax.annotate(texts[i], xy=(x, y), xytext=(tx, ty), fontsize=fontsize, ha='center', va='bottom', color=color,
                    arrowprops=dict(arrowstyle='-', lw=0.4, color='0.4', shrinkA=0, shrinkB=1.5))
    return placed


def save(fig, name):
    fig.savefig(os.path.join(OUT, f'{name}.pdf')); fig.savefig(os.path.join(OUT, f'{name}.png'), dpi=300)


def place_labels(ax, pts, texts, all_xy, fontsize=6, color='0.15', r_pt=(10, 16, 24, 34, 46, 60)):
    """按像素避让：对每个标签依次尝试若干半径 × 方位的偏移（单位：点），选第一个与已放标签不重叠、
    且不覆盖任何数据点（all_xy，数据坐标）的位置；标签框与坐标轴边界保持在轴内；用引导线连回原点。"""
    fig = ax.figure; fig.canvas.draw(); T = ax.transData; dpi = fig.dpi
    ax_bb = ax.get_window_extent()
    P = T.transform(np.asarray(all_xy, float)) if len(all_xy) else np.zeros((0, 2))
    placed = []
    dirs = [(0, 1), (0.7, 0.7), (-0.7, 0.7), (1, 0), (-1, 0), (0.9, 0.4), (-0.9, 0.4), (0.7, -0.7), (-0.7, -0.7), (0, -1)]
    order = np.argsort([-p[1] for p in pts])
    for i in order:
        px, py = T.transform(pts[i]); w = 0.52 * fontsize * len(texts[i]) * dpi / 72; h = 1.25 * fontsize * dpi / 72
        best = None
        for r in r_pt:
            for dx, dy in dirs:
                cx = px + dx * r * dpi / 72; cy = py + dy * r * dpi / 72
                bx0 = cx - w / 2 if dx == 0 else (cx if dx > 0 else cx - w); by0 = cy if dy >= 0 else cy - h
                box = (bx0, by0, bx0 + w, by0 + h)
                if box[0] < ax_bb.x0 or box[2] > ax_bb.x1 or box[1] < ax_bb.y0 or box[3] > ax_bb.y1: continue
                if any(not (box[2] < b[0] or box[0] > b[2] or box[3] < b[1] or box[1] > b[3]) for b in placed): continue
                if len(P) and np.any((P[:, 0] > box[0] - 2) & (P[:, 0] < box[2] + 2) & (P[:, 1] > box[1] - 2) & (P[:, 1] < box[3] + 2)): continue
                best = box; break
            if best: break
        if best is None:      # 退路：放在点正上方最远半径处
            r = r_pt[-1]; best = (px - w / 2, py + r * dpi / 72, px + w / 2, py + r * dpi / 72 + h)
        placed.append(best)
        tx, ty = T.inverted().transform(((best[0] + best[2]) / 2, (best[1] + best[3]) / 2))
        ax.annotate(texts[i], xy=pts[i], xytext=(tx, ty), fontsize=fontsize, ha='center', va='center', color=color,
                    arrowprops=dict(arrowstyle='-', lw=0.4, color='0.45', shrinkA=1, shrinkB=1.5), zorder=6)
