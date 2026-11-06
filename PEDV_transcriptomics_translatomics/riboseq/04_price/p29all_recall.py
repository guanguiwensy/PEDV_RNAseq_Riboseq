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
import sample_registry as SR
# Host positive control, same rule as pilot p29_recall.py / price_recall.py, for the four 29-nt scopes (+ old 16-50 nt for 813 and merged_CHX where an old host run exists in price_813)
import os, re
import numpy as np, pandas as pd
from statsmodels.stats.multitest import multipletests
O = os.environ['O']; W = os.environ['W']; P = f'{O}/price_29nt/all4'; P8 = f'{O}/price_813'
stop, gene_of, cds = {}, {}, {}
for ln in open(f'{O}/host4/chr1_reduced.gtf'):
    p = ln.rstrip('\n').split('\t')
    if len(p) < 9: continue
    m = re.search(r'transcript_id "([^"]+)"', p[8])
    if not m: continue
    t = m.group(1)
    if p[2] == 'stop_codon': stop[t] = (p[6], int(p[4]) if p[6] == '+' else int(p[3]))
    if p[2] == 'CDS': cds.setdefault(t, []).append((int(p[3]), int(p[4])))
    g = re.search(r'gene_id "([^"]+)"', p[8])
    if g: gene_of[t] = g.group(1)
TRUTH = [t for t in stop if t in cds]
gene_stops = {}
for ln in open(f'{O}/host4/chr1_all.gtf'):
    p = ln.rstrip('\n').split('\t')
    if len(p) > 8 and p[2] == 'stop_codon':
        g = re.search(r'gene_id "([^"]+)"', p[8]).group(1)
        gene_stops.setdefault(g, set()).add((p[6], int(p[4]) if p[6] == '+' else int(p[3])))
def bh(p): return multipletests(np.nan_to_num(np.asarray(p, float), nan=1.0), method='fdr_bh')[1]
def parse_price(path):
    d = pd.read_csv(path, sep='\t'); d.columns = list(d.columns[:9]) + ['reads', 'Total']
    rows = []
    for loc, pv in zip(d.Location, d['p value']):
        m = re.match(r'(.+?)([+-]):(.*)', loc)
        if m.group(1) != '1': continue
        bl = [tuple(map(int, b.split('-'))) for b in m.group(3).split('|')]
        st = m.group(2); se = max(b[1] for b in bl) if st == '+' else min(b[0] for b in bl) + 1
        rows.append((st, se, pv))
    x = pd.DataFrame(rows, columns=['st', 'se', 'p']); x['q'] = bh(x.p.values)
    return {(a, b): c for a, b, c in x.groupby(['st', 'se']).q.min().reset_index().itertuples(index=False)}
def called(path):
    q = parse_price(path)
    return {t for t in TRUTH if any(q.get(k, 1) < 0.1 for k in gene_stops.get(gene_of.get(t), [stop[t]]))}
def fc(samples):
    return pd.concat([pd.read_csv(f'{W}/host_psite/{s}/host_gene_cds_psite.tsv', sep='\t') for s in samples])
def depth(f, lens): return f[f.L.isin(lens)].groupby('transcript_id')[['f0', 'f1', 'f2']].sum().sum(axis=1)
SC = {SR.lib("RIBO_VERO_CHX_1MOI_5h"): [SR.lib("RIBO_VERO_CHX_1MOI_5h")], SR.lib("RIBO_VERO_CHX_1MOI_8h"): [SR.lib("RIBO_VERO_CHX_1MOI_8h")], SR.lib("RIBO_VERO_CHX_1MOI_16h"): [SR.lib("RIBO_VERO_CHX_1MOI_16h")], 'merged_CHX': [SR.lib("RIBO_VERO_CHX_1MOI_5h"), SR.lib("RIBO_VERO_CHX_1MOI_8h"), SR.lib("RIBO_VERO_CHX_1MOI_16h")]}
rows = []
for sc, sams in SC.items():
    f = fc(sams); DEP = {'strict_core 29+30': depth(f, (29, 30)), '29nt only': depth(f, (29,))}
    runs = [('29nt (all4)', f'{P}/host/price/{sc}.orfs.tsv')]
    if sc in (SR.lib("RIBO_VERO_CHX_1MOI_16h"), 'merged_CHX'): runs.append(('16-50nt (price_813)', f'{P8}/host/price/{sc}.orfs.tsv'))
    for lab, fn in runs:
        if not os.path.exists(fn) or os.path.getsize(fn) == 0: rows.append(dict(scope=sc, run=lab, note='no output')); continue
        c = called(fn)
        for dn, dep in DEP.items():
            for k in (10, 30, 100):
                T = [x for x in TRUTH if dep.get(x, 0) >= k]; hit = [x for x in T if x in c]
                rows.append(dict(scope=sc, run=lab, depth_proxy=dn, stratum=f'CDS P-sites>={k}', n_truth=len(T), n_called=len(hit), recall=round(len(hit) / len(T), 4) if T else np.nan))
out = pd.DataFrame(rows); out.to_csv(f'{P}/price_host_recall_29nt_all4.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); print(out.to_string())
