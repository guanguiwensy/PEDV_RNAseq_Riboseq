#!/usr/bin/env python3
"""Single panel: the top 20 counted candidate ORFs (truncated, frameshifted, out-of-frame) in one genome track.
Ranked by supporting junctions across the three counted classes; coloured by class.
Canonical maximum (N) dashed line as in the v2 genome figure. Excluded junctions are not plotted.
Inputs: ../../Noncanonical_reclass/tables/T1_new_ORF_ranking.csv, Fig2B and fig1_depth (ORF1a reference).
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
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, LogFormatterMathtext

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
TAB = BASE / 'Noncanonical_reclass' / 'tables'
SRC = HERE / 'source_data'
SRC.mkdir(exist_ok=True)
MM = 1 / 25.4
BLUE = '#2878B5'; RED = '#D94841'; ORANGE = '#E69F00'; GRAY = '#8C8C8C'; DARK = '#222222'; GREYTXT = '#444444'
CLS = {'truncated': ('Truncated canonical frame', BLUE), 'frameshift': ('Frameshifted', RED),
       'out-of-frame ORF': ('Out-of-frame ORF', ORANGE)}
FS = 7.2
plt.rcParams.update({'font.family': ['Arial', 'DejaVu Sans'], 'font.size': FS, 'pdf.fonttype': 42,
                     'svg.fonttype': 'none', 'axes.linewidth': 0.6})
GENES = [('ORF1a', 293, 12622, '#9ECAE1'), ('ORF1b', 12592, 20613, '#4292C6'), ('S', 20610, 24758, '#FDAE6B'),
         ('ORF3', 24758, 25033, '#A1D99B'), ('E', 25361, 25600, '#FDD0A2'), ('M', 25608, 26288, '#74C476'),
         ('N', 26300, 27625, '#FCAE91')]
CANON = {'ORF1a': 293, 'S': 20597, 'ORF3': 24718, 'E': 25342, 'M': 25597, 'N': 26285}

T1 = pd.read_csv(TAB / 'T1_new_ORF_ranking.csv')
T1['label'] = T1.start_gene.astype(str).str.replace('intergenic/UTR', "3'UTR") + ':' + T1.orf_start.astype(int).astype(str)
top = T1.sort_values(['n_supporting_junctions', 'reads'], ascending=False).head(20).reset_index(drop=True)
top['rank'] = np.arange(1, len(top) + 1)
top[['rank', 'cls', 'label', 'orf_start', 'n_supporting_junctions', 'reads', 'n_libraries_detected']].to_csv(
    SRC / 'top20_single_panel.csv', index=False)

canon = pd.read_csv(BASE / 'GIb_rnaseq_reanalysis/figures/Fig2B_composition_source_data.csv').set_index('segment').reads.to_dict()
depth = BASE / 'GIb_rnaseq_reanalysis' / 'figdata' / 'fig1_depth'
vals = []
for f in sorted(depth.glob('R2608219*.dedup.depth.tsv')) + sorted(depth.glob('R2608220*.dedup.depth.tsv')):
    x = pd.read_csv(f, sep='\t', header=None, names=['chr', 'pos', 'depth'])
    vals.append(int(x.loc[x.pos == 293, 'depth'].iloc[0]))
canon['ORF1a'] = int(sum(vals))
CMAX = float(max(canon[g] for g in ['S', 'ORF3', 'E', 'M', 'N']))

W, H = 180, 92
fig = plt.figure(figsize=(W * MM, H * MM))
ax = fig.add_axes([0.10, 0.25, 0.86, 0.62])
x = top.orf_start.values.astype(float)
y = np.maximum(top.reads.values.astype(float), 1)
j = top.n_supporting_junctions.values.astype(float)
cols = [CLS[c][1] for c in top.cls]
ax.set_yscale('log'); ax.set_xlim(0, 28000); ax.set_ylim(1, CMAX * 1.6)
ax.axhline(CMAX, color=GRAY, ls=(0, (4, 2)), lw=0.8, zorder=1)
ax.text(27850, CMAX * 1.1, f'canonical maximum (N): {CMAX:,.0f} reads', ha='right', va='bottom', fontsize=FS - 0.4, color=GREYTXT)
ax.vlines(x, 1, y, color=cols, lw=0.9, alpha=0.8, zorder=2)
ax.scatter(x, y, s=8 + 1.6 * j, c=cols, edgecolors='white', linewidths=0.4, alpha=0.9, zorder=3)
factors = [1.6, 3.0, 6.0, 12.0, 24.0]
for k, idx in enumerate(np.argsort(x)):
    ax.annotate(top.label.iloc[idx], xy=(x[idx], y[idx]), xytext=(x[idx] + 200, y[idx] * factors[k % 5]),
                textcoords='data', rotation=50, ha='left', va='bottom', fontsize=FS - 0.2, color=DARK,
                arrowprops=dict(arrowstyle='-', color=cols[idx], lw=0.4, alpha=0.6, shrinkA=0, shrinkB=0),
                annotation_clip=False)
for name, pos in CANON.items():
    r = canon[name]
    ax.vlines(pos, 1, r, color='#333333', lw=0.9, zorder=2)
    ax.scatter([pos], [r], s=20, color='#333333', zorder=3)
ax.yaxis.set_major_locator(LogLocator(base=10, numticks=6))
ax.yaxis.set_major_formatter(LogFormatterMathtext(base=10))
ax.set_ylabel('Supporting reads (log scale)', fontsize=FS)
ax.set_xticks([])
for s in ('top', 'right', 'bottom'):
    ax.spines[s].set_visible(False)
handles = [Line2D([], [], marker='o', ls='none', color=GRAY, markersize=np.sqrt(8 + 1.6 * v) * 0.9, label=str(v))
           for v in (10, 30, 60, 90)]
leg1 = ax.legend(handles=handles, title='Supporting junctions', frameon=True, framealpha=1, edgecolor='#888888',
                 fontsize=FS - 0.4, title_fontsize=FS - 0.4, loc='upper left', ncol=4, borderaxespad=0.3)
ax.add_artist(leg1)
cls_handles = [Line2D([], [], marker='o', ls='none', color=c, markersize=4, label=CLS[k][0]) for k, c in
               [('truncated', BLUE), ('frameshift', RED), ('out-of-frame ORF', ORANGE)]]
ax.legend(handles=cls_handles, frameon=False, fontsize=FS - 0.4, loc='lower right', borderaxespad=0.3)
# gene track
tr = fig.add_axes([0.10, 0.10, 0.86, 0.09])
tr.set_xlim(0, 28000); tr.set_ylim(-1.8, 1.0)
for name, s, e, c in GENES:
    tr.add_patch(Rectangle((s, 0.1), e - s, 0.8, fc=c, ec='white', lw=0.3))
for name, xm in {'ORF1a': 7000, 'ORF1b': 16300, 'S': 22700}.items():
    tr.text(xm, 0.5, name, ha='center', va='center', fontsize=FS, color=DARK)
for name, (xm, yb) in {'ORF3': (24900, -0.6), 'E': (25480, -1.45), 'M': (25960, -0.6), 'N': (27050, -1.45)}.items():
    tr.text(xm, yb, name, ha='center', va='center', fontsize=FS - 0.4, color=DARK)
tr.set_yticks([])
for sp in ('left', 'top', 'right'):
    tr.spines[sp].set_visible(False)
tr.set_xticks(np.arange(0, 28001, 5000))
tr.set_xticklabels([f'{v:,}' for v in np.arange(0, 28001, 5000)], fontsize=FS - 0.4)
tr.set_xlabel('Genome position (nt)', fontsize=FS)
fig.text(0.5, 0.97, 'Top 20 counted candidate ORFs (truncated, frameshifted, out-of-frame)', ha='center', va='top',
         fontsize=9, fontweight='bold')
fig.text(0.5, 0.01, 'Point x = candidate start ATG; size = supporting junctions; y = supporting reads. '
                    'Excluded junctions (translating a canonical ORF) are not shown.', ha='center', fontsize=FS - 0.4,
         color=GREYTXT)
for ext, kw in [('png', dict(dpi=300)), ('pdf', {}), ('svg', {})]:
    fig.savefig(HERE / f'top20_counted_single_panel.{ext}', facecolor='white', **kw)
print(top[['rank', 'cls', 'label', 'n_supporting_junctions', 'reads']].to_string(index=False))
