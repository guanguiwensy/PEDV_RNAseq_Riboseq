#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
import sys, os, csv, gzip, json, collections
root, per_path = sys.argv[1:3]; sids = sys.argv[3:]

ref = {}
for line in open(per_path, encoding='utf-8'):
    if line.startswith('#') or line.startswith('sample_id'): continue
    f = line.rstrip('\n').split('\t')
    if f[3] == 'host' and f[4] == '29':
        ref[f[0]] = (int(f[7]), float(f[8]), float(f[9]), float(f[10]))

KEYS = ['primary_mapped', 'NH_gt1_skipped', 'L29_reads', 'psite_NA', 'sc5_reads', 'sense_assigned', 'sense_multi_tx',
        'antisense_only', 'unassigned', 'sense_outside_cds_flank', 'elapsed_s', 'bam_size', 'n_transcripts', 'n_split_cds_tx']
rows = []
for sid in sids:
    d = os.path.join(root, sid)
    meta = json.load(open(f'{d}/meta_{sid}.json'))
    tx = {}
    with open(f'{d}/transcripts.tsv') as f:
        for r in csv.DictReader(f, delimiter='\t'): tx[int(r['tx_idx'])] = r
    clen = {k: int(r['cds_len_nt']) for k, r in tx.items()}; strand = {k: r['strand'] for k, r in tx.items()}
    body = [0, 0, 0]; body5 = [0, 0, 0]; bys = {'+': [0, 0, 0], '-': [0, 0, 0]}
    tot = 0; tot_sc5 = 0; tot_multi = 0; tot_incds = 0; tx_nonzero = set()
    with gzip.open(f'{d}/pos_{sid}.tsv.gz', 'rt') as f:
        next(f)
        for line in f:
            k, rel, sc5, ns, n = map(int, line.split('\t'))
            tot += n; tot_sc5 += n * sc5; tot_multi += n * (ns > 1)
            L = clen[k]
            if 0 <= rel <= L - 1: tot_incds += n; tx_nonzero.add(k)
            if ns != 1: continue
            if 45 <= rel <= L - 46: body[rel % 3] += n; bys[strand[k]][rel % 3] += n
            if 57 <= rel <= L - 34: body5[rel % 3] += n
    anti = 0
    with open(f'{d}/antisense_{sid}.tsv') as f:
        next(f)
        for line in f: anti += int(line.split('\t')[1])
    gen = collections.Counter()
    with gzip.open(f'{d}/genomic_splitcds_{sid}.tsv.gz', 'rt') as f:
        next(f)
        for line in f:
            x = line.rstrip('\n').split('\t'); gen[int(x[0])] += int(x[4])
    oaz = {}
    for gname in ('OAZ1', 'OAZ2'):
        ks = [k for k, r in tx.items() if r['gene_name'] == gname]
        oaz[gname] = ';'.join(f"{tx[k]['transcript_id']}|split={tx[k]['split_cds_gap_le10']}|genomic_rows_counts={gen.get(k, 0)}" for k in ks) or 'NOT_IN_ANNOTATION'
    fr = lambda v: [round(x / max(sum(v), 1), 4) for x in v]
    r = collections.OrderedDict(sid=sid)
    for k in KEYS: r[k] = int(meta.get(k, 0))
    r['pos_total_count'] = tot; r['pos_incds_count'] = tot_incds; r['tx_with_cds_counts'] = len(tx_nonzero)
    r['body_n'] = sum(body); r['body_F0'], r['body_F1'], r['body_F2'] = fr(body)
    r['body5_n'] = sum(body5); r['body5_F0'], r['body5_F1'], r['body5_F2'] = fr(body5)
    rf = ref.get(sid)
    if rf:
        r['ref_n'], r['ref_F0'], r['ref_F1'], r['ref_F2'] = rf
        r['maxabs_diff_body_vs_ref'] = round(max(abs(a - b) for a, b in zip(fr(body), rf[1:])), 4)
        r['maxabs_diff_body5_vs_ref'] = round(max(abs(a - b) for a, b in zip(fr(body5), rf[1:])), 4)
        r['ref_consistent_0.01'] = int(r['maxabs_diff_body_vs_ref'] <= 0.01)
    else:
        r['ref_n'] = r['ref_F0'] = r['ref_F1'] = r['ref_F2'] = 'NA'
        r['maxabs_diff_body_vs_ref'] = r['maxabs_diff_body5_vs_ref'] = r['ref_consistent_0.01'] = 'NA(不在 periodicity_host)'
    for s in '+-':
        nm = 'plus' if s == '+' else 'minus'
        r[f'{nm}_body_n'] = sum(bys[s]); r[f'{nm}_F0'], r[f'{nm}_F1'], r[f'{nm}_F2'] = fr(bys[s])
    r['maxabs_diff_plus_vs_minus'] = round(max(abs(a - b) for a, b in zip(fr(bys['+']), fr(bys['-']))), 4)
    sa = r['sense_assigned']; ao = r['antisense_only']
    r['antisense_only_frac_of_assigned'] = round(ao / max(sa + ao, 1), 5)
    r['sense_to_antisense_ratio'] = round(sa / max(ao, 1), 1)
    r['antisense_tx_counts_total'] = anti
    r['sc5_frac_L29_reads'] = round(r['sc5_reads'] / max(r['L29_reads'], 1), 5)
    r['sc5_frac_pos_counts'] = round(tot_sc5 / max(tot, 1), 5)
    r['multi_tx_frac_sense_reads'] = round(r['sense_multi_tx'] / max(sa, 1), 5)
    r['multi_tx_frac_pos_counts'] = round(tot_multi / max(tot, 1), 5)
    r['OAZ1'] = oaz['OAZ1']; r['OAZ2'] = oaz['OAZ2']
    rows.append(r)
with open(os.path.join(root, 'QC_summary.tsv'), 'w') as f:
    f.write('\t'.join(rows[0].keys()) + '\n')
    for r in rows: f.write('\t'.join(map(str, r.values())) + '\n')
for r in rows: print(json.dumps(r, ensure_ascii=False))
