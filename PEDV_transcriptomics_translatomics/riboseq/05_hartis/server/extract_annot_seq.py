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
import sys, os, re, gzip, collections, pysam
gtf, fa, txtsv, out = sys.argv[1:5]; os.makedirs(out, exist_ok=True)
rx = re.compile(r'(\S+) "([^"]*)"'); rt = re.compile(r'tag "([^"]*)"')
FEAT = ('exon', 'CDS', 'start_codon', 'stop_codon', 'five_prime_utr', 'three_prime_utr')
ex = collections.defaultdict(list); n = 0
with gzip.open(os.path.join(out, 'canonical_pc_features.tsv.gz'), 'wt') as fo:
    fo.write('chrom\tfeature\tstart\tend\tstrand\tphase\ttranscript_id\tgene_id\tgene_name\texon_number\ttags\n')
    for line in open(gtf):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if f[2] not in FEAT or 'tag "Ensembl_canonical"' not in f[8]: continue
        d = dict(rx.findall(f[8]))
        if d.get('transcript_biotype') != 'protein_coding': continue
        fo.write('\t'.join([f[0], f[2], f[3], f[4], f[6], f[7], d['transcript_id'], d['gene_id'], d.get('gene_name', ''),
                            d.get('exon_number', ''), ','.join(rt.findall(f[8]))]) + '\n'); n += 1
        if f[2] == 'exon': ex[d['transcript_id']].append((f[0], f[6], int(f[3]), int(f[4])))
F = pysam.FastaFile(fa); comp = str.maketrans('ACGTNacgtn', 'TGCANtgcan'); m = 0
with gzip.open(os.path.join(out, 'tx_seq.fa.gz'), 'wt') as fo:
    for line in open(txtsv):
        if line.startswith('tx_idx'): continue
        k, t = line.split('\t')[:2]
        el = ex[t]; strand = el[0][1]
        el = sorted(el, key=lambda x: x[2], reverse=(strand == '-'))
        s = ''.join(F.fetch(c, b - 1, e) if strand == '+' else F.fetch(c, b - 1, e).translate(comp)[::-1] for c, st, b, e in el)
        fo.write(f'>{k} {t}\n{s.upper()}\n'); m += 1
print('features', n, 'seqs', m)
