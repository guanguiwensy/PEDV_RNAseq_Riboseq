#!/usr/bin/env python3
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
# Figure 2 plotting (PEDV junction landscape). Input: figures/Fig2_work/events_main.tsv (from fig2_classify.py)
import numpy as np, pandas as pd, json, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.colors import Normalize
from matplotlib import cm

B = PEDV_P+'GIb_rnaseq_reanalysis/'
W = B + 'figures/Fig2_work/'
OUT = B + 'figures/'
plt.rcParams.update({'font.family': 'Arial', 'pdf.fonttype': 42, 'svg.fonttype': 'none',
                     'axes.linewidth': 0.6, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
                     'xtick.major.size': 2, 'ytick.major.size': 2, 'font.size': 7})
FS_T, FS_L, FS_K = 8, 8, 7
MM = 1 / 25.4
GL = 27959
E = pd.read_csv(W + 'events_main.tsv', sep='\t')
TOT = int(E.reads.sum())
GENES = [("5\u2032UTR", 1, 292), ('ORF1a', 293, 12622), ('ORF1b', 12592, 20613), ('S', 20610, 24758),
         ('ORF3', 24758, 25033), ('E', 25361, 25600), ('M', 25608, 26288), ('N', 26300, 27625),
         ("3\u2032UTR", 27626, GL)]
GENE_FILL = {"5\u2032UTR": '#EDEDED', 'ORF1a': '#F3E3EC', 'ORF1b': '#E1EBF5', 'S': '#E6F1E3', 'ORF3': '#FBF0D9',
             'E': '#EAE3F3', 'M': '#F8E3E3', 'N': '#E0F0F0', "3\u2032UTR": '#EDEDED'}
CANON_GENES = ['S', 'ORF3', 'E', 'M', 'N']
CANON_COL = {'S': '#9ECAE1', 'ORF3': '#6BAED6', 'E': '#4292C6', 'M': '#2171B5', 'N': '#08306B'}
COL = {'dep': '#B8860B', 'distant': '#4DAF4A', 'local': '#984EA3',
       'E_in': '#00838F', 'E_out': '#E66100', 'E_na': '#9E9E9E', 'E_nc': '#CFCFCF',
       'F_in': '#1B7837', 'F_out': '#A6DBA0', 'F_nc': '#BDBDBD',
       'G_in': '#762A83', 'G_out': '#C2A5CF', 'G_nc': '#BDBDBD'}
HOT = (27620, 27682)


def axmm(fig, x, y, w, h, size, **kw):
    return fig.add_axes([x / size[0], y / size[1], w / size[0], h / size[1]], **kw)


def new_fig(w, h):
    f = plt.figure(figsize=(w * MM, h * MM), dpi=300)
    f.patch.set_facecolor('white')
    return f


def letter(fig, x, y, s, size):
    fig.text(x / size[0], y / size[1], s, fontsize=13, fontweight='bold', va='top', ha='left')


def gene_bar(ax, y0, y1, xscale=1000.0, orient='h'):
    for n, s, e in GENES:
        if orient == 'h':
            ax.add_patch(Rectangle((s / xscale, y0), (e - s + 1) / xscale, y1 - y0, fc=GENE_FILL[n], ec='#9A9A9A', lw=0.3))
        else:
            ax.add_patch(Rectangle((y0, s / xscale), y1 - y0, (e - s + 1) / xscale, fc=GENE_FILL[n], ec='#9A9A9A', lw=0.3))


# ---------------------------------------------------------------- data products
def hot_list():
    dep = E[E.cls == 'TRSL_dep_noncanonical'].sort_values('reads', ascending=False).head(3)
    ind = E[E.cls.str.contains('indep')].sort_values('reads', ascending=False)
    ind = pd.concat([ind[(ind.donor == HOT[0]) & (ind.acceptor == HOT[1])], ind[~((ind.donor == HOT[0]) & (ind.acceptor == HOT[1]))].head(2)])
    h = pd.concat([dep, ind]).reset_index(drop=True)
    h['label'] = list('abcdef')
    h['note'] = ['TRS-L-dep noncanonical'] * 3 + ['TRS-L-indep, authenticity pending (audit in progress)'] + ['TRS-L-indep'] * 2
    return h[['label', 'donor', 'acceptor', 'reads', 'n_samples', 'cls', 'category', 'note']]


HOTS = hot_list()
BIN = 100
NB = int(np.ceil(GL / BIN))


def matrices():
    M = np.zeros((NB, NB))
    for d, a, r in zip(E.donor, E.acceptor, E.reads):
        M[(a - 1) // BIN, (d - 1) // BIN] += r  # row=acceptor, col=donor
    return M


MAT = matrices()


# ---------------------------------------------------------------- Panel A
def panel_A(fig, ox, oy, size):
    """local frame 0..113 x 0..123 ; (ox,oy) offset in fig mm"""
    def A(x, y, w, h, **kw): return axmm(fig, ox + x, oy + y, w, h, size, **kw)
    hm = A(25, 31, 78, 76)
    Z = np.where(MAT > 0, np.log2(np.maximum(MAT, 1)), np.nan)
    cmap = plt.get_cmap('viridis').copy(); cmap.set_bad('white')
    hm.imshow(Z, origin='lower', extent=[0, NB * BIN / 1000, 0, NB * BIN / 1000], cmap=cmap, vmin=0, vmax=16,
              interpolation='nearest', aspect='auto')
    hm.set_xlim(0, GL / 1000); hm.set_ylim(0, GL / 1000)
    hm.plot([0, GL / 1000], [0, GL / 1000], color='#BDBDBD', lw=0.5, ls='--')
    hm.set_xticks([0, 5, 10, 15, 20, 25]); hm.set_yticks([])
    hm.tick_params(labelbottom=False, labelsize=FS_K)
    for sp in hm.spines.values(): sp.set_linewidth(0.6)
    # markers
    for r in HOTS.itertuples():
        pend = (r.donor == HOT[0] and r.acceptor == HOT[1])
        x, y = r.donor / 1000, r.acceptor / 1000
        hm.plot(x, y, marker='D' if pend else 'o', ms=5 if pend else 4.5, mfc='none', mec='#D55E00', mew=0.8)
        hm.text(x + 0.9 if not pend else x - 0.5, y + (0.15 if r.label != 'f' else -0.7) if not pend else y + 1.0, r.label, fontsize=FS_K, color='#D55E00',
                ha='left' if not pend else 'right', va='center')
    # colorbar (inset, empty lower-right triangle)
    cax = A(60, 47, 38, 2.2)
    cb = fig.colorbar(cm.ScalarMappable(Normalize(0, 16), cmap=cmap), cax=cax, orientation='horizontal', extend='max')
    cb.set_ticks([0, 4, 8, 12, 16]); cb.ax.tick_params(labelsize=FS_K, length=2, width=0.5)
    cb.outline.set_linewidth(0.4)
    fig.text((ox + 79) / size[0], (oy + 53) / size[1], 'Junction-spanning reads per bin\n(log$_2$; 100-nt bins; clipped at 2$^{16}$)',
             fontsize=FS_K, ha='center', va='bottom', linespacing=1.15)
    fig.text((ox + 60) / size[0], (oy + 43.0) / size[1], 'Orange markers a\u2013f:\nhotspots (source data);\ndiamond d: authenticity\npending',
             fontsize=FS_K, ha='left', va='top', color='#D55E00', linespacing=1.15)
    # left marginal: acceptor frequency
    ml = A(13, 31, 10, 76)
    acc = MAT.sum(1); ml.fill_betweenx((np.arange(NB) + 0.5) * BIN / 1000, acc + 1e-9, 1e-9, step=None, color='#555555', lw=0)
    ml.set_xscale('log'); ml.set_xlim(3e6, 0.8); ml.set_ylim(0, GL / 1000)
    ml.set_xticks([1e0, 1e6]); ml.set_xticklabels(['$10^0$', '$10^6$']); ml.get_xticklabels()[0].set_ha('right'); ml.set_xlabel('reads', fontsize=FS_K, labelpad=1)
    ml.set_yticks([0, 5, 10, 15, 20, 25]); ml.tick_params(labelsize=FS_K)
    ml.set_ylabel('3\u2032 breakpoint (kb)', fontsize=FS_L, labelpad=1)
    ml.spines['top'].set_visible(False); ml.spines['left'].set_visible(False); ml.yaxis.tick_right(); ml.yaxis.set_label_position('left')
    ml.tick_params(axis='y', labelleft=True, labelright=False, left=True, right=False)
    # bottom gene bar
    gb = A(25, 26.5, 78, 3.5); gb.set_xlim(0, GL / 1000); gb.set_ylim(0, 1); gb.axis('off')
    gene_bar(gb, 0, 1)
    for n, xc, inside in [('ORF1a', 6.4, 1), ('ORF1b', 16.5, 1), ('S', 22.7, 1), ('M', 26.0, 1), ('N', 26.98, 1)]:
        gb.text(xc, 0.5, n, fontsize=FS_K, ha='center', va='center')
    # vertical gene bar (right of heatmap)
    gv = A(103.6, 31, 2.6, 76); gv.set_ylim(0, GL / 1000); gv.set_xlim(0, 1); gv.axis('off')
    gene_bar(gv, 0, 1, orient='v')
    ys = {'ORF1a': 6.4, 'ORF1b': 16.5, 'S': 22.7}
    for n, y in ys.items(): gv.text(0.5, y, n, fontsize=FS_K, rotation=90, ha='center', va='center')
    # small genes: spread labels to the right with leader lines (y in kb -> mm)
    tgt = [('ORF3', 24.9), ('E', 25.48), ('M', 25.95), ('N', 26.96)]
    k = 76 / (GL / 1000)
    pos = [t[1] * k for t in tgt]
    lab = list(pos)
    for i in range(1, len(lab)):
        if lab[i] - lab[i - 1] < 3.2: lab[i] = lab[i - 1] + 3.2
    # shift block down so top label stays inside
    over = lab[-1] - 76 + 2
    if over > 0: lab = [l - over for l in lab]
    for (n, y), p, l in zip(tgt, pos, lab):
        fig.add_artist(matplotlib.lines.Line2D([(ox + 106.2) / size[0], (ox + 107.4) / size[0]],
                                               [(oy + 31 + p) / size[1], (oy + 31 + l) / size[1]], color='#777777', lw=0.4))
        fig.text((ox + 107.6) / size[0], (oy + 31 + l) / size[1], n, fontsize=FS_K, va='center', ha='left')
    # bottom-gene-bar labels for small genes with tiers
    for n, xk, dx, tier in [("5\u2032UTR", 0.15, 3.0, 0), ('ORF3', 24.9, -7.5, 1), ('E', 25.48, -1.5, 0), ("3\u2032UTR", 27.8, 2.5, 0)]:
        xm = 25 + xk / (GL / 1000) * 78
        yb = 26.5 - 2.6 - tier * 3.0
        fig.add_artist(matplotlib.lines.Line2D([(ox + xm) / size[0], (ox + xm + dx * 0.4) / size[0]],
                                               [(oy + 26.5) / size[1], (oy + yb + 1.0) / size[1]], color='#777777', lw=0.4))
        fig.text((ox + xm + dx * 0.4) / size[0], (oy + yb) / size[1], n, fontsize=FS_K, ha='center', va='top')
    # donor marginal
    mb = A(25, 8, 78, 10.5)
    don = MAT.sum(0); mb.fill_between((np.arange(NB) + 0.5) * BIN / 1000, don + 1e-9, 1e-9, color='#555555', lw=0)
    mb.set_yscale('log'); mb.set_ylim(0.8, 3e6); mb.set_xlim(0, GL / 1000)
    mb.set_yticks([1, 1e6]); mb.set_yticklabels(['$10^0$', '$10^6$'])
    mb.set_xticks([0, 5, 10, 15, 20, 25]); mb.tick_params(labelsize=FS_K)
    mb.set_xlabel('5\u2032 breakpoint (kb)', fontsize=FS_L, labelpad=1)
    mb.spines['top'].set_visible(False); mb.spines['right'].set_visible(False)
    mb.set_ylabel('reads', fontsize=FS_K, labelpad=1)
    # schematic (method schematic, no numbers)
    sc = A(2, 111, 106, 11); sc.set_xlim(0, 106); sc.set_ylim(0, 11); sc.axis('off')
    sc.plot([4, 100], [3.2, 3.2], color='#666666', lw=1.0)
    sc.add_patch(Rectangle((4, 2.4), 14, 1.6, fc='#B8860B', ec='none'))
    sc.text(11, 5.6, 'leader', fontsize=FS_K, ha='center', va='bottom', color='#7A5A08')
    sc.add_patch(Rectangle((60, 2.4), 40, 1.6, fc='#BDBDBD', ec='none'))
    sc.text(80, 5.6, 'body', fontsize=FS_K, ha='center', va='bottom', color='#555555')
    sc.plot([18, 60], [3.2, 3.2], color='#666666', lw=0.6, ls=':')
    sc.text(39, 4.3, 'skipped genomic sequence', fontsize=FS_K, ha='center', va='bottom', color='#555555')
    sc.plot([14.5, 14.5], [1.0, 2.4], color='#D55E00', lw=0.8); sc.plot([60, 60], [1.0, 2.4], color='#D55E00', lw=0.8)
    sc.text(14.5, 0.5, "5\u2032 breakpoint", fontsize=FS_K, ha='center', va='top', color='#D55E00')
    sc.text(60, 0.5, "3\u2032 breakpoint", fontsize=FS_K, ha='center', va='top', color='#D55E00')
    sc.add_patch(Rectangle((8, 8.0), 9.0, 1.2, fc='#56B4E9', ec='none')); sc.add_patch(Rectangle((60, 8.0), 14, 1.2, fc='#56B4E9', ec='none'))
    sc.plot([17, 60], [8.6, 8.6], color='#56B4E9', lw=0.8)
    sc.text(40, 10.2, 'junction-spanning read', fontsize=FS_K, ha='center', va='center', color='#2B7BB0')
    letter(fig, ox + 1, oy + 123, 'A', size)


# ---------------------------------------------------------------- Panel B
def class_pct():
    g = E.groupby('cls').reads.sum()
    can = E[E.cls == 'canonical'].groupby('canonical_gene').reads.sum()
    return g, can


def panel_B(fig, ox, oy, size):
    """local frame 115 x 38"""
    g, can = class_pct()
    segs = [('S', can['S'], CANON_COL['S']), ('ORF3', can['ORF3'], CANON_COL['ORF3']), ('E', can['E'], CANON_COL['E']),
            ('M', can['M'], CANON_COL['M']), ('N', can['N'], CANON_COL['N']),
            ('dep', g['TRSL_dep_noncanonical'], COL['dep']), ('distant', g['TRSL_indep_distant'], COL['distant']),
            ('local', g['TRSL_indep_local'], COL['local'])]
    ax = axmm(fig, ox + 13, oy + 15, 92, 6, size)
    x = 0; rows = []
    for n, v, c in segs:
        p = v / TOT * 100
        ax.barh(0, p, left=x, color=c, edgecolor='white', linewidth=0.3, height=1)
        rows.append(dict(segment=n, reads=int(v), pct_of_junction_spanning_reads=p, start_pct=x, end_pct=x + p)); x += p
    ax.set_xlim(0, 100); ax.set_ylim(-0.5, 0.5); ax.set_yticks([])
    ax.set_xticks([0, 20, 40, 60, 80, 100]); ax.tick_params(labelsize=FS_K)
    ax.set_xlabel('Junction-spanning reads (%%; n = %s)' % format(TOT, ','), fontsize=FS_L, labelpad=1)
    for sp in ['top', 'right', 'left']: ax.spines[sp].set_visible(False)
    can_end = rows[4]['end_pct']
    ax.plot([0, 0, can_end, can_end], [0.6, 0.7, 0.7, 0.6], color='k', lw=0.6, clip_on=False)
    ax.text(0, 0.8, 'D  canonical: %.1f%%' % can_end, fontsize=FS_K, ha='left', va='bottom')
    for i, (L, txt) in enumerate([('E', 'TRS-L-dep noncanonical'), ('F', 'TRS-L-indep distant'), ('G', 'TRS-L-indep local')]):
        r = rows[5 + i]; xm = (r['start_pct'] + r['end_pct']) / 2
        lev = [1.0, 1.7, 2.4][i]
        ax.plot([xm, xm], [0.5, lev], color='#777777', lw=0.4, clip_on=False)
        ax.text(xm - 0.8, lev, '%s  %s: %.1f%%' % (L, txt, r['pct_of_junction_spanning_reads']), fontsize=FS_K, ha='right', va='center')
    fig.text((ox + 13) / size[0], (oy + 2.0) / size[1],
             'Canonical, %% of all reads: S %.2f | ORF3 %.2f | E %.2f | M %.2f | N %.2f' % tuple(can[k] / TOT * 100 for k in CANON_GENES),
             fontsize=FS_K, ha='left', va='bottom')
    letter(fig, ox + 1, oy + 38, 'B', size)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- Panel C
def species():
    n = E[E.cls == 'TRSL_dep_noncanonical'].sort_values('acceptor').copy()
    n['grp'] = (n.acceptor.diff() > 5).cumsum()
    rows = []
    for gid, x in n.groupby('grp'):
        top = x.sort_values('reads', ascending=False).iloc[0]
        rows.append(dict(species=f'TRS-L-dep noncanonical @{int(top.acceptor)}', kind='noncanonical', rep_acceptor=int(top.acceptor),
                         acceptor_min=int(x.acceptor.min()), acceptor_max=int(x.acceptor.max()), reads=int(x.reads.sum()),
                         n_events=len(x), frame=top.frame, orf_consequence=top.orf_consequence))
    nc = pd.DataFrame(rows).sort_values('reads', ascending=False)
    can = E[E.cls == 'canonical'].groupby('canonical_gene').reads.sum().sort_values(ascending=False)
    c = pd.DataFrame([dict(species=f'canonical {g}', kind='canonical', rep_acceptor=None, acceptor_min=None, acceptor_max=None,
                           reads=int(v), n_events=int((E.canonical_gene == g).sum()), frame='NA', orf_consequence='NA', gene=g) for g, v in can.items()])
    return pd.concat([c, nc], ignore_index=True)


SP = species()


def panel_C(fig, ox, oy, size):
    """local frame 115 x 52"""
    can = SP[SP.kind == 'canonical']; nc = SP[SP.kind == 'noncanonical']
    a1 = axmm(fig, ox + 13, oy + 9, 24, 40, size)
    a2 = axmm(fig, ox + 40, oy + 9, 72, 40, size, sharey=a1)
    a1.bar(np.arange(5), can.reads.values, width=0.8, color=[CANON_COL[g] for g in can.gene], linewidth=0)
    a1.set_xticks(np.arange(5)); a1.set_xticklabels(list(can.gene), fontsize=FS_K)
    a1.set_yscale('log'); a1.set_ylim(1, 5e6)
    a1.set_yticks([1, 10, 1e2, 1e3, 1e4, 1e5, 1e6]); a1.set_yticklabels(['$10^0$', '$10^1$', '$10^2$', '$10^3$', '$10^4$', '$10^5$', '$10^6$'])
    a1.tick_params(labelsize=FS_K); a1.set_xlim(-0.6, 4.6)
    a1.set_ylabel('Reads (log$_{10}$ scale)', fontsize=FS_L, labelpad=1)
    a1.set_xlabel('Canonical', fontsize=FS_L, labelpad=1)
    cols = [{'in': '#8C8C8C', 'out': '#000000'}.get(f, '#DADADA') for f in nc.frame]
    a2.bar(np.arange(len(nc)) + 1, nc.reads.values, width=0.9, color=cols, linewidth=0)
    a2.set_xlim(0, len(nc) + 1); a2.tick_params(labelsize=FS_K, labelleft=False)
    a2.set_xticks([1, 100, 200, 300, len(nc)])
    a2.set_xlabel('TRS-L-dependent noncanonical species (n = %d), ranked' % len(nc), fontsize=FS_L, labelpad=1)
    for a in (a1, a2):
        for sp in ['top', 'right']: a.spines[sp].set_visible(False)
    hs = [matplotlib.patches.Patch(fc='#8C8C8C', label='in-frame (N-terminally truncated ORF, internal AUG)'),
          matplotlib.patches.Patch(fc='#000000', label='out-of-frame (frameshifted ORF)'),
          matplotlib.patches.Patch(fc='#DADADA', label='no AUG in acceptor CDS / non-coding acceptor')]
    a2.legend(handles=hs, fontsize=FS_K, frameon=False, loc='upper right', handlelength=1.0, borderaxespad=0.1, labelspacing=0.25)
    letter(fig, ox + 1, oy + 52, 'C', size)


# ---------------------------------------------------------------- Panels D-G
def arcs(ax, sub, color_fn):
    """draw all events; lw/alpha from reads. Returns number drawn"""
    sub = sub.sort_values('reads')
    t = np.linspace(0, np.pi, 28)
    segs, lws, cols = [], [], []
    for d, a, r, c in zip(sub.donor, sub.acceptor, sub.reads, [color_fn(x) for x in sub.itertuples()]):
        xs = d / 1000 + (a - d) / 1000 * (1 - np.cos(t)) / 2
        h = 0.10 + 0.72 * np.log10(max(a - d, 50) / 50) / np.log10(GL / 50)
        ys = h * np.sin(t)
        segs.append(np.column_stack([xs, ys]))
        u = np.log10(r) / 6.3
        lws.append(0.12 + 1.6 * u)
        rgba = matplotlib.colors.to_rgba(c, 0.12 + 0.75 * u)
        cols.append(rgba)
    ax.add_collection(LineCollection(segs, linewidths=lws, colors=cols, capstyle='round'))


AH = 22.0   # arc axes height (mm)
BH = 46     # standard block height
EXG = 13


def block_h(n2, rows, extra=0):
    return 7.5 + AH + 4.0 * rows + 0.6 + 3.2 * n2 + 4.6 + extra


def dg_panel(fig, ox, oy, size, title1, title2, sub, color_fn, letter_s, n2, rows, extra=0, legend=None, labels=False, hot=False, ncol=3):
    H = block_h(n2, rows, extra)
    ax = axmm(fig, ox + 2, oy + 7.5 + extra, 60, AH, size)
    ax.set_xlim(0, GL / 1000); ax.set_ylim(0, 1.03)
    arcs(ax, sub, color_fn)
    ax.set_yticks([]); ax.set_xticks([0, 5, 10, 15, 20, 25]); ax.tick_params(labelbottom=False, labelsize=FS_K)
    for s_ in ['top', 'right', 'left']: ax.spines[s_].set_visible(False)
    if hot:
        ax.plot(HOT[0] / 1000, 0.13, marker='D', ms=4, mfc='none', mec='#D55E00', mew=0.8)
        ax.annotate('d (authenticity pending)', (HOT[0] / 1000, 0.13), xytext=(HOT[0] / 1000 - 4.0, 0.75), fontsize=FS_K, color='#D55E00',
                    ha='right', va='bottom', arrowprops=dict(arrowstyle='-', color='#D55E00', lw=0.5))
    tr = axmm(fig, ox + 2, oy + 4.0 + extra, 60, 3.5, size); tr.set_xlim(0, GL / 1000); tr.set_ylim(0, 1); tr.axis('off')
    gene_bar(tr, 0, 1)
    if labels:
        for n, xc in [('ORF1a', 6.4), ('ORF1b', 16.5), ('S', 22.7), ('M', 26.0), ('N', 26.98)]:
            tr.text(xc, 0.5, n, fontsize=FS_K, ha='center', va='center')
        for n, xk, dx, tier, ha in [("5′UTR", 0.15, 3.0, 0, 'center'), ('ORF3', 24.9, -9.0, 1, 'center'), ('E', 25.48, -1.0, 0, 'center'), ("3′UTR", 27.8, -1.0, 1, 'right')]:
            xm = 2 + xk / (GL / 1000) * 60
            xt = xm + dx * 0.4 if ha == 'center' else 62.8
            yb = 4.0 - 1.8 - tier * 3.0 + extra
            fig.add_artist(matplotlib.lines.Line2D([(ox + xm) / size[0], (ox + xt - (3.5 if ha == 'right' else 0)) / size[0]],
                                                   [(oy + 4.0 + extra) / size[1], (oy + yb + 0.9) / size[1]], color='#777777', lw=0.4))
            fig.text((ox + xt) / size[0], (oy + yb) / size[1], n, fontsize=FS_K, ha=ha, va='top')
        fig.text((ox + 2) / size[0], (oy + 0.3) / size[1], 'Genome position (kb); each arc = one junction;\nall events drawn (no read threshold).',
                 fontsize=FS_K, ha='left', va='bottom', linespacing=1.15)
    fig.text((ox + 9) / size[0], (oy + H - 0.3) / size[1], title1, fontsize=FS_T, fontweight='bold', ha='left', va='top')
    fig.text((ox + 9) / size[0], (oy + H - 4.0) / size[1], title2, fontsize=FS_K, ha='left', va='top', linespacing=1.2)
    letter(fig, ox, oy + H, letter_s, size)
    if legend:
        hs = [matplotlib.lines.Line2D([0], [0], color=c, lw=1.8, label=l) for l, c in legend]
        ax.legend(handles=hs, fontsize=FS_K, frameon=False, loc='lower left', bbox_to_anchor=(-0.02, 1.0), handlelength=0.9, labelspacing=0.15,
                  columnspacing=0.8, handletextpad=0.35, borderaxespad=0, ncol=ncol)
    return H


def scale_key(fig, ox, oy, size):
    fig.text((ox + 2) / size[0], (oy + 5.4) / size[1], 'Arc width/opacity scale with log$_{10}$(reads):', fontsize=FS_K, ha='left', va='bottom')
    for i, r in enumerate([1, 10, 100, 1000, 10000, 100000, 1000000]):
        u = np.log10(r) / 6.3
        x0 = 2 + i * 8.6
        fig.add_artist(matplotlib.lines.Line2D([(ox + x0) / size[0], (ox + x0 + 3.2) / size[0]], [(oy + 3.0) / size[1]] * 2,
                                               color=matplotlib.colors.to_rgba('#333333', 0.12 + 0.75 * u), lw=0.12 + 1.6 * u, solid_capstyle='round'))
        fig.text((ox + x0 + 3.8) / size[0], (oy + 3.0) / size[1], ('$10^{%d}$' % np.log10(r)) if r > 1 else '1', fontsize=FS_K, va='center', ha='left')


def stats(cls):
    x = E[E.cls == cls]; return len(x), int(x.reads.sum())


def pct(x): return x / TOT * 100


def panel_DG(fig, ox, oy, size):
    """column frame 63 wide"""
    can = E[E.cls == 'canonical']
    cl = lambda r: CANON_COL.get(r.canonical_gene, '#08306B')
    hD = block_h(1, 1); hE = block_h(2, 2); hF = block_h(2, 1); hG = block_h(2, 1, EXG)
    yG = 10; yF = yG + hG; yE = yF + hF; yD = yE + hE
    nD, rD = stats('canonical')
    dg_panel(fig, ox, oy + yD, size, 'TRS-L-dependent canonical: %.1f%%' % pct(rD), '%s reads; %d events' % (format(rD, ','), nD),
             can, cl, 'D', 1, 1, legend=[(g, CANON_COL[g]) for g in CANON_GENES], ncol=5)
    dep = E[E.cls == 'TRSL_dep_noncanonical']
    ef = lambda r: {'in': COL['E_in'], 'out': COL['E_out'], 'no_AUG_in_acceptor_CDS': COL['E_na'], 'noncoding': COL['E_nc']}[r.frame]
    fr = dep.groupby('frame').reads.sum(); nE, rE = stats('TRSL_dep_noncanonical')
    dg_panel(fig, ox, oy + yE, size, 'TRS-L-dep. noncanonical: %.1f%%' % pct(rE),
             '%s reads; %s events\n(legend values: %% of all reads)' % (format(rE, ','), format(nE, ',')), dep, ef, 'E', 2, 2,
             legend=[('in-frame %.2f' % pct(fr['in']), COL['E_in']), ('out-of-frame %.2f' % pct(fr['out']), COL['E_out']),
                     ('no AUG in acceptor CDS %.2f' % pct(fr['no_AUG_in_acceptor_CDS']), COL['E_na']), ('non-coding acceptor %.2f' % pct(fr['noncoding']), COL['E_nc'])], ncol=2)
    dis = E[E.cls == 'TRSL_indep_distant']
    ff = lambda r: {'in': COL['F_in'], 'out': COL['F_out'], 'noncoding': COL['F_nc']}[r.donor_frame_continuation]
    fr = dis.groupby('donor_frame_continuation').reads.sum(); nF, rF = stats('TRSL_indep_distant')
    dg_panel(fig, ox, oy + yF, size, 'TRS-L-indep. distant: %.1f%%' % pct(rF),
             '%s reads; %s events; span ≥ 1 kb\n(legend values: %% of all reads)' % (format(rF, ','), format(nF, ',')), dis, ff, 'F', 2, 1,
             legend=[('in-frame %.2f' % pct(fr['in']), COL['F_in']), ('out %.2f' % pct(fr['out']), COL['F_out']), ('non-coding end %.2f' % pct(fr['noncoding']), COL['F_nc'])])
    loc = E[E.cls == 'TRSL_indep_local']
    gf = lambda r: {'in': COL['G_in'], 'out': COL['G_out'], 'noncoding': COL['G_nc']}[r.donor_frame_continuation]
    fr = loc.groupby('donor_frame_continuation').reads.sum(); nG, rG = stats('TRSL_indep_local')
    dg_panel(fig, ox, oy + yG, size, 'TRS-L-indep. local: %.1f%%' % pct(rG),
             '%s reads; %s events; span < 1 kb\n(legend values: %% of all reads)' % (format(rG, ','), format(nG, ',')), loc, gf, 'G', 2, 1, extra=EXG,
             legend=[('in-frame %.2f' % pct(fr['in']), COL['G_in']), ('out %.2f' % pct(fr['out']), COL['G_out']), ('non-coding end %.2f' % pct(fr['noncoding']), COL['G_nc'])],
             labels=True, hot=True)
    scale_key(fig, ox, oy, size)
    return yD + hD


# ---------------------------------------------------------------- assemble
def save(fig, name):
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(OUT + f'{name}.{ext}', dpi=300, facecolor='white')
    plt.close(fig)


def main():
    CW, CH = 180, 215
    Bdf = None
    # combined
    size = (CW, CH); f = new_fig(*size)
    panel_A(f, 0, 92, size); Bdf = panel_B(f, 0, 53, size); panel_C(f, 0, 0, size); panel_DG(f, 117, 2, size)
    save(f, 'Fig2_combined')
    # single panels
    s = (115, 123); f = new_fig(*s); panel_A(f, 0, 0, s); save(f, 'Fig2A_junction_landscape')
    s = (115, 38); f = new_fig(*s); panel_B(f, 0, 0, s); save(f, 'Fig2B_composition')
    s = (115, 52); f = new_fig(*s); panel_C(f, 0, 0, s); save(f, 'Fig2C_canonical_abundance')
    s = (63, 207); f = new_fig(*s); panel_DG(f, 0, 0, s); save(f, 'Fig2D-G_arcs')
    # source data
    zz = np.argwhere(MAT > 0)
    pd.DataFrame([dict(donor_bin_start_nt=int(c * BIN + 1), acceptor_bin_start_nt=int(r * BIN + 1), bin_size_nt=BIN, reads=int(MAT[r, c]),
                       log2_reads=float(np.log2(MAT[r, c])), plotted_log2_clipped_at_16=float(min(16, np.log2(MAT[r, c])))) for r, c in zz]
                 ).to_csv(OUT + 'Fig2A_junction_landscape_source_data.csv', index=False)
    pd.DataFrame(dict(bin_start_nt=np.arange(NB) * BIN + 1, donor_reads=MAT.sum(0).astype(int), acceptor_reads=MAT.sum(1).astype(int))
                 ).to_csv(OUT + 'Fig2A_marginals_source_data.csv', index=False)
    HOTS.to_csv(OUT + 'Fig2A_hotspots_source_data.csv', index=False)
    Bdf.to_csv(OUT + 'Fig2B_composition_source_data.csv', index=False)
    SP.to_csv(OUT + 'Fig2C_canonical_abundance_source_data.csv', index=False)
    ev = E.copy()
    ev['panel'] = ev.cls.map({'canonical': 'D', 'TRSL_dep_noncanonical': 'E', 'TRSL_indep_distant': 'F', 'TRSL_indep_local': 'G'})
    ev.rename(columns={'frame': 'frame_acceptor_first_AUG_rule', 'cls': 'class', 'trsl_dependent': 'leader_region_donor'}, inplace=True)
    ev['drawn_in_arcs'] = True
    ev.drop(columns=['category']).to_csv(OUT + 'Fig2D-G_arcs_source_data.csv', index=False)
    ps = pd.read_csv(W + 'per_sample_events.tsv.gz', sep='\t')
    ps = ps[ps['sample'].isin(SR.LIB_RNA_VERO)]; SR.assert_exactly(set(ps['sample']), SR.RNA_VERO_ALIASES, 'Fig2 per-sample events')
    ps.to_csv(OUT + 'Fig2_per_sample_event_counts.csv.gz', index=False)
    print('done', len(SP), TOT)


if __name__ == '__main__':
    main()
