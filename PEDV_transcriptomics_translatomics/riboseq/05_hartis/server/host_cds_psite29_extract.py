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
import sys, os, re, json, time, gzip, collections
import numpy as np, pysam

sid, bam_path, gtf, outroot = sys.argv[1:5]
OUT = os.path.join(outroot, sid); os.makedirs(OUT, exist_ok=True)
REFC = (0, 2, 7, 8); L_KEEP = 29; OFF = 12; FLANK = 60; GAP_MAX = 10
t0 = time.time()


rx = re.compile(r'(\S+) "([^"]*)"')
ex = collections.defaultdict(list); cds = collections.defaultdict(list); tinfo = {}
for line in open(gtf):
    if line.startswith('#'): continue
    f = line.rstrip('\n').split('\t')
    if f[2] not in ('exon', 'CDS') or 'tag "Ensembl_canonical"' not in f[8]: continue
    d = dict(rx.findall(f[8]))
    if d.get('transcript_biotype') != 'protein_coding': continue
    t = d['transcript_id']; tinfo[t] = (f[0], f[6], d['gene_id'], d.get('gene_name', ''))
    (ex if f[2] == 'exon' else cds)[t].append((int(f[3]), int(f[4])))

TID, CS, CE, TL, STR, CHR, SPLIT = [], [], [], [], [], [], []
blocks = collections.defaultdict(list)
for t, el in ex.items():
    if t not in cds: continue
    chrom, strand, g, gn = tinfo[t]
    el = sorted(el, key=lambda x: x[0], reverse=(strand == '-'))
    def tx(pos):
        c = 0
        for b, e in el:
            if b <= pos <= e: return c + (pos - b if strand == '+' else e - pos)
            c += e - b + 1
        return None
    cl = cds[t]
    cs_g = min(b for b, e in cl) if strand == '+' else max(e for b, e in cl)
    ce_g = max(e for b, e in cl) if strand == '+' else min(b for b, e in cl)
    cs, ce = tx(cs_g), tx(ce_g)
    if cs is None or ce is None: continue

    sc = sorted(cl); gaps = [sc[i + 1][0] - sc[i][1] - 1 for i in range(len(sc) - 1)]
    k = len(TID); TID.append((t, g, gn)); CS.append(cs); CE.append(ce); STR.append(strand); CHR.append(chrom)
    TL.append(sum(e - b + 1 for b, e in el)); SPLIT.append(int(any(0 < x <= GAP_MAX for x in gaps)))
    c = 0
    for b, e in el:
        blocks[(chrom, strand)].append((b, e, c, k)); c += e - b + 1
HB = {kk: np.array(sorted(v), dtype=np.int64) for kk, v in blocks.items()}
CMX = {kk: np.maximum.accumulate(a[:, 1]) for kk, a in HB.items()}
NT = len(TID)


def hits(chrom, strand, p):
    a = HB.get((chrom, strand))
    if a is None: return []
    i = np.searchsorted(a[:, 0], p, side='right'); cm = CMX[(chrom, strand)]
    out = []
    j = i - 1
    while j >= 0 and cm[j] >= p:
        gs, ge, cb, k = a[j]
        if gs <= p <= ge:
            out.append((int(k), int(cb + (p - gs) if strand == '+' else cb + (ge - p))))
        j -= 1
    return out


def psite(r, rev, cig):
    if any(op == 3 for op, n in cig):
        coords = []; x = r.reference_start
        for op, n in cig:
            if op in REFC: coords.extend(range(x + 1, x + n + 1)); x += n
            elif op == 3: x += n
        if rev: coords.reverse()
        return coords[OFF] if OFF < len(coords) else None
    return (r.reference_end - OFF) if rev else (r.reference_start + 1 + OFF)


SENSE = collections.Counter(); ANTI = collections.Counter(); GEN = collections.Counter()
qc = collections.Counter()
split_tx = {k for k in range(NT) if SPLIT[k]}
split_span = collections.defaultdict(list)
for k in split_tx:
    t = TID[k][0]; cl = cds[t]; split_span[CHR[k]].append((min(b for b, e in cl), max(e for b, e in cl), STR[k], k))
bam = pysam.AlignmentFile(bam_path)
for chrom in bam.references:
    if chrom.startswith('PEDV'): continue
    for r in bam.fetch(chrom):
        if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
        qc['primary_mapped'] += 1
        nh = r.get_tag('NH') if r.has_tag('NH') else 1
        if nh > 1: qc['NH_gt1_skipped'] += 1; continue
        cig = r.cigartuples
        L0 = sum(n for op, n in cig if op in REFC)
        if L0 != L_KEEP: continue
        qc['L29_reads'] += 1
        rev = r.is_reverse; st = '-' if rev else '+'
        p = psite(r, rev, cig)
        if p is None: qc['psite_NA'] += 1; continue
        sc5 = int((cig[-1][0] == 4) if rev else (cig[0][0] == 4))
        qc['sc5_reads'] += sc5
        hs = hits(chrom, st, p); ha = hits(chrom, '+' if rev else '-', p)
        if hs:
            ns = len(hs); qc['sense_assigned'] += 1; qc['sense_multi_tx'] += (ns > 1)
            for k, c in hs:
                rel = c - CS[k]
                if -FLANK <= rel <= (CE[k] - CS[k]) + FLANK:
                    SENSE[(k, rel, sc5, ns)] += 1
                else: qc['sense_outside_cds_flank'] += 1
        elif ha:
            qc['antisense_only'] += 1
            for k, c in ha: ANTI[k] += 1
        else:
            qc['unassigned'] += 1
        for gs, ge, s2, k in split_span.get(chrom, []):
            if gs - FLANK <= p <= ge + FLANK and s2 == st: GEN[(k, chrom, p, st)] += 1
bam.close()

with gzip.open(os.path.join(OUT, f'pos_{sid}.tsv.gz'), 'wt') as f:
    f.write('tx_idx\trel\tsc5\tn_sense_tx\tcount\n')
    for (k, rel, sc5, ns), n in sorted(SENSE.items()): f.write(f'{k}\t{rel}\t{sc5}\t{ns}\t{n}\n')
with open(os.path.join(OUT, f'antisense_{sid}.tsv'), 'w') as f:
    f.write('tx_idx\tantisense_L29\n')
    for k, n in sorted(ANTI.items()): f.write(f'{k}\t{n}\n')
with gzip.open(os.path.join(OUT, f'genomic_splitcds_{sid}.tsv.gz'), 'wt') as f:
    f.write('tx_idx\tchrom\tgenomic_psite\tstrand\tcount\n')
    for (k, c, p, s2), n in sorted(GEN.items()): f.write(f'{k}\t{c}\t{p}\t{s2}\t{n}\n')
with open(os.path.join(OUT, 'transcripts.tsv'), 'w') as f:
    f.write('tx_idx\ttranscript_id\tgene_id\tgene_name\tchrom\tstrand\ttx_len\tcds_start_tx\tcds_end_tx\tcds_len_nt\tsplit_cds_gap_le10\n')
    for k, (t, g, gn) in enumerate(TID):
        f.write(f'{k}\t{t}\t{g}\t{gn}\t{CHR[k]}\t{STR[k]}\t{TL[k]}\t{CS[k]}\t{CE[k]}\t{CE[k]-CS[k]+1}\t{SPLIT[k]}\n')
meta = dict(sid=sid, bam=bam_path, bam_size=os.path.getsize(bam_path), gtf=gtf, read_length=L_KEEP, offset=OFF,
            n_transcripts=NT, n_split_cds_tx=len(split_tx), elapsed_s=round(time.time() - t0), **{k: int(v) for k, v in qc.items()})
json.dump(meta, open(os.path.join(OUT, f'meta_{sid}.json'), 'w'), indent=1, ensure_ascii=False)
print(json.dumps(meta, ensure_ascii=False))
