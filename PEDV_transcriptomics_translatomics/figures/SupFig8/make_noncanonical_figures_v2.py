#!/usr/bin/env python3
"""Two figures under the reclassified scheme (Noncanonical_reclass/build_new_tables.py).

Counted classes: truncated (canonical frame, N-terminally truncated), frameshifted, out-of-frame ORF.
Excluded from all panels and statistics: translates_canonical (123 junctions, kept as noncanonical sgRNAs).
Panels use only the counted classes; the excluded class is shown once, in the workflow and class summary, as not counted.

Inputs: ../../Noncanonical_reclass/tables/ (T1_new_ORF_ranking, T3_new_time_course, T4_new_class_summary, T4_new_distinct_ORFs)
Outputs: noncanonical_v2_group_figure.* and noncanonical_v2_genome_top20.* in this folder; source_data/.
"""
import os as _os, sys as _sys
_d = _os.path.dirname(_os.path.abspath(__file__))
while _d != _os.path.dirname(_d) and not _os.path.isdir(_os.path.join(_d, "config")):
    _d = _os.path.dirname(_d)
_sys.path.insert(0, _os.path.join(_d, "config"))
from pedv_paths import PEDV_P, PEDV_K, PEDV_R, PEDV_COV, PEDV_SERVER_ROOT, REPO_ROOT, SCRIPT_DIRS, ORIGIN
_sys.path.extend(_x for _x in SCRIPT_DIRS if _x not in _sys.path)
__REPO_FILE__ = globals().get("__REPO_FILE__") or _os.path.abspath(__file__)
_o = ORIGIN.get(_os.path.relpath(_os.path.abspath(__file__), REPO_ROOT).replace(chr(92), "/"))
if _o:
    __file__ = _o
from pathlib import Path
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Arc
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, LogFormatterMathtext

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
TAB = BASE / 'Noncanonical_reclass' / 'tables'
SRC = HERE / 'source_data'
SRC.mkdir(exist_ok=True)
MM = 1 / 25.4
BLUE = '#2878B5'; RED = '#D94841'; ORANGE = '#E69F00'; GRAY = '#8C8C8C'; LGRAY = '#BDBDBD'
DARK = '#222222'; GREYTXT = '#444444'; FRAME = '#9A9A9A'
TP = [3, 5, 8, 12, 16]
CLS = {'truncated': ('Truncated canonical frame', BLUE), 'frameshift': ('Frameshifted', RED),
       'out-of-frame ORF': ('Out-of-frame ORF', ORANGE)}
FS = 7.2; FST = 8.2; FSP = 12
plt.rcParams.update({'font.family': ['Arial', 'DejaVu Sans'], 'font.size': FS, 'axes.labelsize': FS,
                     'xtick.labelsize': FS, 'ytick.labelsize': FS, 'pdf.fonttype': 42, 'svg.fonttype': 'none',
                     'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
                     'xtick.major.size': 2, 'ytick.major.size': 2})
CANON = {'ORF1a': 293, 'S': 20597, 'ORF3': 24718, 'E': 25342, 'M': 25597, 'N': 26285}
GENES = [('ORF1a', 293, 12622, '#9ECAE1'), ('ORF1b', 12592, 20613, '#4292C6'), ('S', 20610, 24758, '#FDAE6B'),
         ('ORF3', 24758, 25033, '#A1D99B'), ('E', 25361, 25600, '#FDD0A2'), ('M', 25608, 26288, '#74C476'),
         ('N', 26300, 27625, '#FCAE91')]

# ---------------------------------------------------------------- inputs
T1 = pd.read_csv(TAB / 'T1_new_ORF_ranking.csv')
T3 = pd.read_csv(TAB / 'T3_new_time_course.csv')
T4 = pd.read_csv(TAB / 'T4_new_class_summary.csv')
TD = pd.read_csv(TAB / 'T4_new_distinct_ORFs.csv').set_index('cls').distinct_ORFs.to_dict()
EV = pd.read_csv(TAB / 'events_with_new_ORF_assignment.csv')
PS = pd.read_csv(TAB / 'per_library_support_counted.csv')
T1['label'] = T1.start_gene.astype(str).str.replace('intergenic/UTR', "3'UTR") + ':' + T1.orf_start.astype(int).astype(str)
for c in CLS:
    assert c in set(T1.cls)

# canonical reference values: Fig. 2B gene totals (20 libraries); ORF1a = dedup reads spanning 293 (20 libraries)
CANON_READS = pd.read_csv(BASE / 'GIb_rnaseq_reanalysis/figures/Fig2B_composition_source_data.csv').set_index('segment').reads.to_dict()
_d = BASE / 'GIb_rnaseq_reanalysis' / 'figdata' / 'fig1_depth'
_v = []
for f in sorted(_d.glob('R2608219*.dedup.depth.tsv')) + sorted(_d.glob('R2608220*.dedup.depth.tsv')):
    x = pd.read_csv(f, sep='\t', header=None, names=['chr', 'pos', 'depth'])
    _v.append(int(x.loc[x.pos == 293, 'depth'].iloc[0]))
assert len(_v) == 20
CANON_READS['ORF1a'] = int(sum(_v))
CANON_MAX = float(max(CANON_READS[g] for g in ['S', 'ORF3', 'E', 'M', 'N']))


def top_n(cls, n):
    d = T1[T1.cls == cls].sort_values('rank_within_class').head(n).copy()
    return d


def time_sum(orf_start, cls):
    s = T3[(T3.cls == cls) & (T3.orf_start == orf_start)]
    return [int(s[s.hpi == h].reads.sum()) for h in TP]


def frame_fig(fig, W, H, x, y, w, h):
    fig.add_artist(Rectangle((x / W, 1 - (y + h) / H), w / W, h / H, transform=fig.transFigure,
                             fill=False, ec=FRAME, lw=0.6))


# =====================================================================================  group figure
def fig_group():
    W, H = 180, 248
    fig = plt.figure(figsize=(W * MM, H * MM))

    def ax_mm(x, y, w, h):
        return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])

    def t_mm(x, y, s, **kw):
        kw.setdefault('fontsize', FS)
        return fig.text(x / W, 1 - y / H, s, **kw)

    def head(letter, x, y, text):
        t_mm(x, y, letter, fontsize=FSP, fontweight='bold', va='top', ha='left')
        t_mm(x + 9, y + 1.2, text, fontsize=FST, fontweight='bold', va='top', ha='left')

    t_mm(W / 2, 2, 'Noncanonical sgRNA-associated candidate ORFs (reclassified)', fontsize=9.5,
         fontweight='bold', ha='center', va='top')

    # ---------------- A: workflow
    head('A', 1, 7, 'Analysis workflow and classification')
    ax = ax_mm(0, 12, W, 55)
    ax.set_xlim(0, W); ax.set_ylim(55, 0); ax.axis('off')

    def box(x, y, w, h, face, edge, ls='-'):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0,rounding_size=2.2', fc=face, ec=edge,
                                    lw=0.8, linestyle=ls))

    def arrow(x1, y1, x2, y2):
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle='-|>', color=GRAY, lw=0.7, mutation_scale=7, shrinkA=0, shrinkB=0))

    n_ev = {r.new_class: (int(r.events), int(r.reads)) for r in T4.itertuples()}
    box(56, 0, 70, 12, '#F4F4F4', '#888888')
    ax.text(91, 3.2, 'Noncanonical sgRNA events (TRS-L dependent)', ha='center', va='center', fontsize=FST, fontweight='bold')
    ax.add_patch(Rectangle((66, 8.4), 12, 2.4, fc=BLUE, ec='none'))
    ax.add_patch(Rectangle((104, 8.4), 12, 2.4, fc=RED, ec='none'))
    ax.add_patch(Arc((91, 8.4), 22, 7, theta1=180, theta2=360, color=GRAY, lw=0.7))
    ax.text(72, 11.2, "5' donor", ha='center', va='center', fontsize=FS - 0.4, color=GREYTXT)
    ax.text(110, 11.2, "3' acceptor", ha='center', va='center', fontsize=FS - 0.4, color=GREYTXT)
    arrow(91, 12, 91, 15)
    box(56, 15, 70, 12, '#F4F4F4', '#888888')
    ax.text(91, 18.2, 'ORF assignment (sgRNA ATG…stop)', ha='center', va='center', fontsize=FST, fontweight='bold')
    ax.add_patch(Rectangle((66, 23.2), 30, 2.4, fc=BLUE, ec='none'))
    ax.add_patch(Rectangle((100, 23.2), 22, 2.4, fc=RED, ec='none'))
    ax.text(91, 20.6, 'in canonical frame / out of frame', ha='center', va='center', fontsize=FS - 0.4, color=GREYTXT)
    ax.plot([91, 91], [27, 30], color=GRAY, lw=0.7)
    centers = [25, 69, 113, 157]
    ax.plot([centers[0], centers[-1]], [30, 30], color=GRAY, lw=0.7)
    for cx in centers:
        arrow(cx, 30, cx, 32)
    tr = n_ev['truncated (counted)']; fs = n_ev['frameshift (counted)']; oo = n_ev['out-of-frame ORF (counted)']
    ex = n_ev['translates_canonical (excluded)']
    cls_boxes = [
        (6, '#DDECF7', BLUE, '-', [('Truncated', True), ('canonical frame', False),
                                   (f'{TD["truncated"]} ORFs · {tr[0]} events', False)]),
        (50, '#FBE4E1', RED, '-', [('Frameshifted', True), ('different frame from CDS', False),
                                   (f'{TD["frameshift"]} ORFs · {fs[0]} events', False)]),
        (94, '#FDF1DC', ORANGE, '-', [('Out-of-frame ORF', True), ('ATG…stop, not canonical frame', False),
                                      (f'{TD["out-of-frame ORF"]} ORFs · {oo[0]} events', False)]),
        (138, '#F2F2F2', GRAY, '--', [('Not counted', True), ('translates canonical ORF', False),
                                      (f'{ex[0]} events (kept as noncanonical)', False)]),
    ]
    for x, face, edge, ls, lines in cls_boxes:
        box(x, 32, 38, 13, face, edge, ls)
        for i, (s, bold) in enumerate(lines):
            ax.text(x + 19, 32 + 13 * (i + 0.5) / len(lines), s, ha='center', va='center', fontsize=FS - 0.2,
                    fontweight='bold' if bold else 'normal', color=DARK)
    box(36, 47, 108, 7, '#EAF4EC', '#4D7F5A')
    ax.text(90, 50.5, 'Rank high-frequency candidates and evaluate temporal dynamics (counted classes only)',
            ha='center', va='center', fontsize=FS - 0.2, color=DARK)
    arrow(25, 45, 52, 47); arrow(69, 45, 80, 47); arrow(113, 45, 90, 47)

    # ---------------- B: class summary (counted classes; excluded shown as not counted)
    head('B', 1, 69, 'Junction events, supporting reads and distinct ORFs by class')
    order = ['truncated (counted)', 'frameshift (counted)', 'out-of-frame ORF (counted)']
    tick = ['Truncated', 'Frameshift', 'Out-of-\nframe']
    cols = [BLUE, RED, ORANGE]
    specs = [(0, 'Junction events'), (1, 'Supporting reads'), (2, 'Distinct ORFs')]
    B_src = []
    for k, title in specs:
        x = [12, 72, 132][k]
        axb = ax_mm(x, 80, 44, 24)
        if k == 2:
            vals = [TD['truncated'], TD['frameshift'], TD['out-of-frame ORF']]
        else:
            vals = [int(T4.set_index('new_class').loc[c, 'events' if k == 0 else 'reads']) for c in order]
        axb.bar(range(3), vals, color=cols, width=0.62, edgecolor='none')
        for i, v in enumerate(vals):
            axb.text(i, v + max(vals) * 0.03, f'{v:,}', ha='center', va='bottom', fontsize=FS)
        axb.set_xticks(range(3)); axb.set_xticklabels(tick, fontsize=FS)
        axb.set_ylim(0, max(vals) * 1.2)
        axb.set_title(title, fontsize=FS, fontweight='bold', loc='left', pad=2)
        axb.spines['top'].set_visible(False); axb.spines['right'].set_visible(False)
        axb.tick_params(axis='y', labelsize=FS - 0.4)
        B_src += [dict(panel='B', metric=title, cls=c, value=int(v)) for c, v in zip(order, vals)]
    ex_txt = f"Not counted: {ex[0]} events / {ex[1]:,} reads translate a canonical ORF (kept as noncanonical sgRNAs)."
    t_mm(12, 108, ex_txt, fontsize=FS - 0.4, color=GREYTXT, va='top')
    pd.DataFrame(B_src).to_csv(SRC / 'v2_panelB_class_summary.csv', index=False)

    # ---------------- C: top 8 truncated and top 8 frameshifted (counted)
    head('C', 1, 114, 'Top candidates ranked by supporting junctions')
    for x, cls, col, title in [(28, 'truncated', BLUE, 'Truncated canonical-frame'), (118, 'frameshift', RED, 'Frameshifted')]:
        d = top_n(cls, 8)
        axc = ax_mm(x, 126, 58, 32)
        yy = np.arange(len(d))[::-1]
        mx = d.n_supporting_junctions.max()
        axc.barh(yy, d.n_supporting_junctions, color=col, height=0.66)
        for y_, (_, r) in zip(yy, d.iterrows()):
            axc.text(r.n_supporting_junctions + mx * 0.02, y_, f'{int(r.n_supporting_junctions)}', va='center', fontsize=FS)
            axc.text(mx * 1.32, y_, f'reads={int(r.reads):,}, libs={int(r.n_libraries_detected)}', va='center',
                     fontsize=FS - 0.3, color=GREYTXT)
        axc.set_yticks(yy); axc.set_yticklabels(d.label, fontsize=FS)
        axc.set_xlim(0, mx * 2.35)
        axc.set_xticks([0, round(mx / 2), mx])
        axc.set_title(title + ' (junctions)', fontsize=FS, fontweight='bold', loc='left', pad=2)
        axc.spines['top'].set_visible(False); axc.spines['right'].set_visible(False)
    # ---------------- D/E: time heatmaps (top 8 each)
    head('D', 1, 166, 'Truncated candidates across time')
    head('E', 96, 166, 'Frameshifted candidates across time')
    heat_rows = []
    for x_ax, x_cb, cls, lbl in [(27, 88, 'truncated', 'D'), (118, 170, 'frameshift', 'E')]:
        d = top_n(cls, 8)
        mat = np.array([time_sum(int(r.orf_start), cls) for _, r in d.iterrows()], dtype=float)
        logm = np.log2(mat + 1)
        axh = ax_mm(x_ax, 174, 58 if lbl == 'D' else 48, 28)
        im = axh.imshow(logm, aspect='auto', interpolation='nearest', cmap='viridis')
        vmax = logm.max()
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                axh.text(j, i, f'{int(mat[i, j]):,}', ha='center', va='center', fontsize=FS - 0.4,
                         color='white' if logm[i, j] < 0.55 * vmax else DARK)
        axh.set_yticks(range(len(d))); axh.set_yticklabels(d.label, fontsize=FS)
        axh.set_xticks(range(len(TP))); axh.set_xticklabels([str(t) for t in TP], fontsize=FS)
        axh.set_xlabel('Hours post infection (hpi)', fontsize=FS)
        axh.tick_params(length=0)
        for s in axh.spines.values():
            s.set_visible(False)
        cb = fig.colorbar(im, cax=ax_mm(x_cb, 174, 2.6, 28))
        cb.set_label('log$_2$(reads + 1)', fontsize=FS - 0.4)
        cb.ax.tick_params(labelsize=FS - 0.4, length=1.5)
        for (_, r), row in zip(d.iterrows(), mat):
            heat_rows.append(dict(panel=lbl, label=r.label, cls=cls, **{f'reads_{t}h': v for t, v in zip(TP, row)}))
    pd.DataFrame(heat_rows).to_csv(SRC / 'v2_panelDE_heatmap_sums.csv', index=False)
    axc.set_xticks  # noqa (kept for layout)
    # ---------------- F: class trend + distinct ORFs
    head('F', 1, 216, 'Class-level temporal trend')
    axf = ax_mm(30, 224, 72, 18)
    F_src = []
    for cls in ['truncated', 'frameshift', 'out-of-frame ORF']:
        lab, col = CLS[cls]
        sums = [int(T3[(T3.cls == cls) & (T3.hpi == h)].reads.sum()) for h in TP]
        axf.plot(TP, sums, marker='o', ms=3.4, lw=1.4, color=col, label=lab)
        for t, v in zip(TP, sums):
            axf.text(t, v + max(sums + [1]) * 0.04, f'{v:,}', ha='center', va='bottom', fontsize=FS - 0.4, color=col)
            F_src.append(dict(cls=cls, hpi=t, summed_reads=v))
    axf.set_xticks(TP); axf.set_xlabel('Hours post infection (hpi)', fontsize=FS)
    axf.set_ylabel('Summed supporting reads', fontsize=FS)
    axf.set_ylim(0, T3.groupby(['cls', 'hpi']).reads.sum().max() * 1.3)
    axf.legend(frameon=False, fontsize=FS - 0.4, loc='upper left')
    axf.spines['top'].set_visible(False); axf.spines['right'].set_visible(False)
    pd.DataFrame(F_src).to_csv(SRC / 'v2_panelF_class_trend.csv', index=False)
    rows = []
    for cls in ['truncated', 'frameshift', 'out-of-frame ORF']:
        s = T3[T3.cls == cls]
        rows.append([CLS[cls][0]] + [str(s[(s.hpi == t) & (s.reads > 0)][['start_gene', 'orf_start']].drop_duplicates().shape[0]) for t in TP])
    t_mm(112, 222, 'Distinct ORFs with ≥1 read', fontsize=FS, fontweight='bold', va='top', ha='left')
    axt = ax_mm(110, 230, 68, 16)
    axt.axis('off')
    tb = axt.table(cellText=rows, colLabels=['class'] + [f'{t} hpi' for t in TP], loc='center', cellLoc='center',
                   colWidths=[0.34] + [0.13] * 5)
    tb.set_fontsize(FS - 0.4); tb.scale(1, 1.3)
    for cell in tb.get_celld().values():
        cell.set_linewidth(0.4)

    for (fy, fh) in [(4, 63), (67, 44), (113, 52), (165, 47), (215, 31)]:
        frame_fig(fig, W, H, 1, fy, 178, fh)
    for ext, kw in [('png', dict(dpi=300)), ('pdf', {}), ('svg', {})]:
        fig.savefig(HERE / f'noncanonical_v2_group_figure.{ext}', facecolor='white', **kw)
    plt.close(fig)


# =====================================================================================  genome figure
def fig_genome():
    W, H = 180, 160
    fig = plt.figure(figsize=(W * MM, H * MM))

    def ax_mm(x, y, w, h):
        return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])

    def t_mm(x, y, s, **kw):
        kw.setdefault('fontsize', FS)
        return fig.text(x / W, 1 - y / H, s, **kw)

    def head(letter, y, text):
        t_mm(1, y, letter, fontsize=FSP, fontweight='bold', va='top', ha='left')
        t_mm(10, y + 1.2, text, fontsize=FST, fontweight='bold', va='top', ha='left')

    t_mm(W / 2, 1.5, 'Genome positions and support of candidate ORFs (reclassified, counted classes)',
         fontsize=9.5, fontweight='bold', ha='center', va='top')
    XL = (0, 28000)

    def track(y, h):
        ax = ax_mm(28, y, 150, h)
        ax.set_xlim(*XL); ax.set_ylim(-1.8, 1.0)
        for name, s, e, c in GENES:
            ax.add_patch(Rectangle((s, 0.1), e - s, 0.8, fc=c, ec='white', lw=0.3))
        for name, xm in {'ORF1a': 7000, 'ORF1b': 16300, 'S': 22700}.items():
            ax.text(xm, 0.5, name, ha='center', va='center', fontsize=FS, color=DARK)
        for name, (xm, yb) in {'ORF3': (24900, -0.6), 'E': (25480, -1.45), 'M': (25960, -0.6), 'N': (27050, -1.45)}.items():
            ax.text(xm, yb, name, ha='center', va='center', fontsize=FS - 0.4, color=DARK)
        ax.set_yticks([])
        for sp in ('left', 'top', 'right'):
            ax.spines[sp].set_visible(False)
        ax.set_xticks(np.arange(0, 28001, 5000))
        ax.set_xticklabels([f'{v:,}' for v in np.arange(0, 28001, 5000)], fontsize=FS - 0.4)

    # ---------------- A: canonical reference (reads, log)
    head('A', 4, 'Canonical sgRNA reference (reads; ORF1a = reads spanning position 293)')
    axA = ax_mm(28, 12, 150, 24)
    axA.set_xlim(*XL); axA.set_yscale('log'); axA.set_ylim(1, CANON_MAX * 1.6)
    axA.yaxis.set_major_locator(LogLocator(base=10, numticks=6))
    axA.yaxis.set_major_formatter(LogFormatterMathtext(base=10))
    axA.set_ylabel('Reads\n(log scale)', fontsize=FS)
    axA.set_xticks([])
    for s_ in ('top', 'right', 'bottom'):
        axA.spines[s_].set_visible(False)
    for i, (name, pos) in enumerate(sorted(CANON.items(), key=lambda kv: kv[1])):
        r = CANON_READS[name]
        axA.vlines(pos, 1, r, color='#333333', lw=0.9, zorder=2)
        axA.scatter([pos], [r], s=20, color='#333333', zorder=3)
        axA.annotate(f'{name} {r:,}', xy=(pos, r), xytext=(pos + 200, r * [1.6, 3.0, 6.0, 12.0, 24.0][i % 5]),
                     textcoords='data', rotation=50, ha='left', va='bottom', fontsize=FS - 0.2, color=DARK,
                     arrowprops=dict(arrowstyle='-', color='#333333', lw=0.4, alpha=0.6, shrinkA=0, shrinkB=0),
                     annotation_clip=False)
    track(36, 9)

    # ---------------- B/C: candidates
    src_rows = []

    def lolli(h, df, colours, letter, title, yhead, legend_colour):
        head(letter, yhead, title)
        ax = ax_mm(28, yhead + 7, 150, h)
        x = df.orf_start.values.astype(float)
        y_ = np.maximum(df.reads.values.astype(float), 1)
        j = df.n_supporting_junctions.values.astype(float)
        ax.set_yscale('log')
        ax.set_ylim(1, CANON_MAX * 1.6)
        ax.axhline(CANON_MAX, color=GRAY, linestyle=(0, (4, 2)), lw=0.8, zorder=1)
        ax.text(XL[1] - 150, CANON_MAX * 1.12, f'canonical maximum (N): {CANON_MAX:,.0f} reads', ha='right',
                va='bottom', fontsize=FS - 0.4, color=GREYTXT)
        ax.set_xlim(*XL)
        ax.vlines(x, 1, y_, color=colours, lw=0.9, alpha=0.8, zorder=2)
        ax.scatter(x, y_, s=8 + 1.6 * j, color=colours, edgecolors='white', linewidths=0.4, alpha=0.9, zorder=3)
        factors = [1.6, 3.0, 6.0, 12.0, 24.0]
        for k, idx in enumerate(np.argsort(x)):
            ax.annotate(df.label.iloc[idx], xy=(x[idx], y_[idx]),
                        xytext=(x[idx] + 200, y_[idx] * factors[k % 5]), textcoords='data', rotation=50,
                        ha='left', va='bottom', fontsize=FS - 0.2, color=DARK,
                        arrowprops=dict(arrowstyle='-', color=legend_colour, lw=0.4, alpha=0.6, shrinkA=0, shrinkB=0),
                        annotation_clip=False)
        for name, pos in CANON.items():
            r = CANON_READS[name]
            ax.vlines(pos, 1, r, color='#333333', lw=0.9, zorder=2)
            ax.scatter([pos], [r], s=20, color='#333333', zorder=3)
        ax.yaxis.set_major_locator(LogLocator(base=10, numticks=6))
        ax.yaxis.set_major_formatter(LogFormatterMathtext(base=10))
        ax.set_ylabel('Supporting reads\n(log scale)', fontsize=FS)
        ax.set_xticks([])
        for sp in ('top', 'right', 'bottom'):
            ax.spines[sp].set_visible(False)
        handles = [Line2D([], [], marker='o', ls='none', color=legend_colour, markersize=np.sqrt(8 + 1.6 * v) * 0.9,
                          label=str(v)) for v in (10, 30, 60, 90)]
        ax.legend(handles=handles, title='Supporting junctions', frameon=True, framealpha=1, edgecolor='#888888',
                  fontsize=FS - 0.4, title_fontsize=FS - 0.4, loc='upper left', ncol=4, handletextpad=0.2,
                  columnspacing=0.9, borderaxespad=0.3)
        track(yhead + 7 + h + 0.5, 9)
        for _, r in df.iterrows():
            src_rows.append(dict(panel=letter, cls=r.cls, rank=int(r.rank_within_class), label=r.label,
                                 orf_start=int(r.orf_start), start_gene=r.start_gene, reads=int(r.reads),
                                 supporting_junctions=int(r.n_supporting_junctions)))

    tr20 = top_n('truncated', 20)
    fs20 = top_n('frameshift', 20)
    oo = T1[T1.cls == 'out-of-frame ORF'].sort_values('rank_within_class')
    lolli(24, tr20, BLUE, 'B', 'Top 20 truncated canonical-frame candidates', 47, BLUE)
    cfs = pd.concat([fs20, oo]).reset_index(drop=True)
    cfs['colour'] = cfs.cls.map({'frameshift': RED, 'out-of-frame ORF': ORANGE})
    lolli(34, cfs.drop(columns=['colour']), RED, 'C',
          'Top 20 frameshifted candidates, and out-of-frame ORFs (orange)', 91, RED)
    t_mm(W / 2, 151, 'Point x = candidate start ATG; size = supporting junctions; y = supporting reads. '
                     'Excluded: junctions translating a canonical ORF. Dashed line: canonical maximum (N).',
         fontsize=FS - 0.6, color=GREYTXT, ha='center', va='top')
    pd.DataFrame(src_rows).to_csv(SRC / 'v2_genome_candidates.csv', index=False)
    for (fy, fh) in [(3, 43), (45, 44), (89, 55)]:
        frame_fig(fig, W, H, 1, fy, 178, fh)
    for ext, kw in [('png', dict(dpi=300)), ('pdf', {}), ('svg', {})]:
        fig.savefig(HERE / f'noncanonical_v2_genome_top20.{ext}', facecolor='white', **kw)
    plt.close(fig)


if __name__ == '__main__':
    fig_group()
    fig_genome()
    print('done')
