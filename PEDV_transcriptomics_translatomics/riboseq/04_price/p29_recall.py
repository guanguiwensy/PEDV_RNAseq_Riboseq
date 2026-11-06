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
# Host positive control for the 29-nt pilot (same rule as price_recall.py: chr1 truth, BH q<0.1 at annotated stop, strata by strict_core CDS P-sites f0+f1+f2)
# + PEDV ORF summary of 29-nt vs 16-50 nt (price_813) runs for RIBO_VERO_CHX_1MOI_16h.
import os, re, json
import numpy as np, pandas as pd
from statsmodels.stats.multitest import multipletests
O = os.environ['O']; W = os.environ['W']; P = f'{O}/price_29nt'; P8 = f'{O}/price_813'
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
fc = pd.read_csv(f'{W}/host_psite/{SR.lib("RIBO_VERO_CHX_1MOI_16h")}/host_gene_cds_psite.tsv', sep='\t')
def depth(lens): return fc[fc.L.isin(lens)].groupby('transcript_id')[['f0', 'f1', 'f2']].sum().sum(axis=1)
# depth proxies: (A) the library's strict_core 29+30 (same as the earlier table), (B) 29 nt only (the actual PRICE input)
DEP = {'strict_core 29+30': depth((29, 30)), '29nt only': depth((29,))}
rows = []
for lab, f in [('16-50nt (price_813)', f'{P8}/host/price/{SR.lib("RIBO_VERO_CHX_1MOI_16h")}.orfs.tsv'), ('29nt pilot (no -filter, MAIN)', f'{P}/host/price/{SR.lib("RIBO_VERO_CHX_1MOI_16h")}.orfs.tsv'), ('29nt BAM + -filter 28:30 (diag_b)', f'{P}/diag_b_filter28_30/{SR.lib("RIBO_VERO_CHX_1MOI_16h")}.orfs.tsv')]:
    if not os.path.exists(f) or os.path.getsize(f) == 0: rows.append(dict(run=lab, note='no output')); continue
    c = called(f)
    for dn, dep in DEP.items():
        for k in (10, 30, 100):
            T = [x for x in TRUTH if dep.get(x, 0) >= k]; hit = [x for x in T if x in c]
            rows.append(dict(run=lab, depth_proxy=dn, stratum=f'CDS P-sites>={k}', n_truth=len(T), n_called=len(hit), recall=round(len(hit) / len(T), 4) if T else np.nan))
out = pd.DataFrame(rows); out.to_csv(f'{P}/price_host_recall_29nt_vs_1650.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); print(out.to_string())
# ---- PEDV ORF summary
ANN = {'ORF1a', 'ORF1b', 'S', 'ORF3', 'E', 'M', 'N'}
def pedv(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0: return None
    d = pd.read_csv(path, sep='\t'); d.columns = list(d.columns[:9]) + ['reads', 'Total']; return d
summ = []
for lab, base in [('16-50nt', f'{P8}/virus/real/{SR.lib("RIBO_VERO_CHX_1MOI_16h")}'), ('29nt', f'{P}/virus/real/{SR.lib("RIBO_VERO_CHX_1MOI_16h")}_29nt')]:
    for st in ('plus', 'minus'):
        d = pedv(f'{base}.{st}.orfs.tsv')
        if d is None: summ.append(dict(run=lab, strand=st, note='no output')); continue
        pv = d['p value'].astype(float)
        summ.append(dict(run=lab, strand=st, n_orfs=len(d), p_lt_0_05=int((pv < 0.05).sum()), q_BH_lt_0_1=int((bh(pv) < 0.1).sum()),
                         types=';'.join(f'{k}:{v}' for k, v in d.Type.value_counts().items())))
s = pd.DataFrame(summ); s.to_csv(f'{P}/pedv_orf_summary_29nt_vs_1650.tsv', sep='\t', index=False); print(s.to_string())
