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

sid, bam_path, gtf, outroot, spec = sys.argv[1:6]
OUT = os.path.join(outroot, sid); os.makedirs(OUT, exist_ok=True)
REFC = (0, 2, 7, 8); FLANK = 60; GAP_MAX = 10; TRIM = 45; CLEAN_NT = 22; BODY = 45
OFFS = {int(a): int(b) for a, b in (x.split(':') for x in spec.split(','))}
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
NT = len(TID); T2K = {t: k for k, (t, g, gn) in enumerate(TID)}


oblocks = collections.defaultdict(list)
for t, bl in cds.items():
    chrom, strand = tinfo[t][0], tinfo[t][1]
    bl = sorted(bl, key=lambda x: x[0], reverse=(strand == '-'))
    tlen = sum(e - b + 1 for b, e in bl); c = 0
    for b, e in bl:
        oblocks[(chrom, strand)].append((b, e, c, tlen, T2K.get(t, -1))); c += e - b + 1
OHB = {kk: np.array(sorted(v), dtype=np.int64) for kk, v in oblocks.items()}


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


def old_coord(chrom, strand, p5):
    a = OHB.get((chrom, strand))
    if a is None: return None
    i = int(np.searchsorted(a[:, 0], p5, side='right')) - 1
    if i < 0: return None
    gs, ge, cb, tlen, k = a[i]
    if p5 > ge: return None
    return (int(cb + (p5 - gs) if strand == '+' else cb + (ge - p5)), int(tlen), int(k))


def walk(r, rev, cig):
    spliced = any(op == 3 for op, n in cig)
    L = 0; first_n = None
    for op, n in (reversed(cig) if rev else cig):
        if op in REFC: L += n
        elif op == 3 and first_n is None: first_n = L
    return spliced, first_n


def psite(r, rev, cig, off):
    if any(op == 3 for op, n in cig):
        coords = []; x = r.reference_start
        for op, n in cig:
            if op in REFC: coords.extend(range(x + 1, x + n + 1)); x += n
            elif op == 3: x += n
        if rev: coords.reverse()
        return coords[off] if off < len(coords) else None
    return (r.reference_end - off) if rev else (r.reference_start + 1 + off)


SENSE = collections.Counter(); ANTI = collections.Counter(); OLD = collections.Counter()
CT = collections.Counter(); DUP = collections.Counter(); qc = collections.Counter()
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
        off = OFFS.get(L0)
        if off is None: continue
        qc[f'L{L0}_reads'] += 1
        rev = r.is_reverse; st = '-' if rev else '+'
        spliced, first_n = walk(r, rev, cig)
        cls = 0 if not spliced else (1 if first_n >= CLEAN_NT else 2)
        p = psite(r, rev, cig, off)
        if p is None: qc[f'L{L0}_psite_NA'] += 1; continue
        sc5 = int((cig[-1][0] == 4) if rev else (cig[0][0] == 4))
        hs = hits(chrom, st, p)
        ns = len(hs); nrec = 0; new_body = 0; new_tx = -1; new_fr = -1
        for k, c in hs:
            rel = c - CS[k]; Lc = CE[k] - CS[k] + 1
            if -FLANK <= rel <= Lc - 1 + FLANK:
                SENSE[(k, L0, rel, sc5, ns, cls)] += 1; nrec += 1
            if ns == 1 and BODY <= rel <= Lc - 1 - BODY:
                new_body = 1; new_tx = k; new_fr = rel % 3
        if ns == 0:
            ha = hits(chrom, '+' if rev else '-', p)
            if ha:
                qc[f'L{L0}_antisense_only'] += 1
                for k, c in ha: ANTI[(k, L0)] += 1
            else: qc[f'L{L0}_unassigned'] += 1
        DUP[(L0, ns, nrec)] += 1

        p5 = r.reference_end if rev else r.reference_start + 1
        oc = old_coord(chrom, st, p5)
        old_body = 0; old_tx = -1; old_fr = -1
        if oc is not None and cls != 2:
            c5, tlen, ko = oc
            if TRIM <= c5 < tlen - TRIM:
                old_body = 1; old_tx = ko; old_fr = (c5 + off) % 3
                OLD[(ko, L0, old_fr)] += 1
        same_tx = int(old_body and new_body and old_tx == new_tx)
        same_fr = int(same_tx and old_fr == new_fr)
        CT[(L0, old_body, new_body, cls, min(ns, 2), same_tx, same_fr)] += 1
bam.close()

with gzip.open(os.path.join(OUT, f'pos_multiL_{sid}.tsv.gz'), 'wt') as f:
    f.write('tx_idx\tread_length\trel\tsc5\tn_sense_tx\tsplice_class\tcount\n')
    for (k, L, rel, sc5, ns, cls), n in sorted(SENSE.items()): f.write(f'{k}\t{L}\t{rel}\t{sc5}\t{ns}\t{cls}\t{n}\n')
with open(os.path.join(OUT, f'antisense_multiL_{sid}.tsv'), 'w') as f:
    f.write('tx_idx\tread_length\tantisense_count\n')
    for (k, L), n in sorted(ANTI.items()): f.write(f'{k}\t{L}\t{n}\n')
with open(os.path.join(OUT, f'oldframe_{sid}.tsv'), 'w') as f:
    f.write('tx_idx\tread_length\tf0\tf1\tf2\n')
    keys = sorted({(k, L) for k, L, fr in OLD})
    for k, L in keys: f.write(f'{k}\t{L}\t{OLD[(k, L, 0)]}\t{OLD[(k, L, 1)]}\t{OLD[(k, L, 2)]}\n')
with open(os.path.join(OUT, f'crosstab_{sid}.tsv'), 'w') as f:
    f.write('read_length\told_body\tnew_body_ns1\tsplice_class\tn_sense_tx_cap2\tsame_tx\tsame_frame\treads\n')
    for key, n in sorted(CT.items()): f.write('\t'.join(map(str, key)) + f'\t{n}\n')
with open(os.path.join(OUT, f'dup_{sid}.tsv'), 'w') as f:
    f.write('read_length\tn_sense_tx\tn_records_written\treads\n')
    for key, n in sorted(DUP.items()): f.write('\t'.join(map(str, key)) + f'\t{n}\n')
with open(os.path.join(OUT, 'transcripts.tsv'), 'w') as f:
    f.write('tx_idx\ttranscript_id\tgene_id\tgene_name\tchrom\tstrand\ttx_len\tcds_start_tx\tcds_end_tx\tcds_len_nt\tsplit_cds_gap_le10\n')
    for k, (t, g, gn) in enumerate(TID):
        f.write(f'{k}\t{t}\t{g}\t{gn}\t{CHR[k]}\t{STR[k]}\t{TL[k]}\t{CS[k]}\t{CE[k]}\t{CE[k]-CS[k]+1}\t{SPLIT[k]}\n')
meta = dict(sid=sid, bam=bam_path, bam_size=os.path.getsize(bam_path), gtf=gtf, lengths_offsets=OFFS,
            n_transcripts=NT, elapsed_s=round(time.time() - t0), NH_gt1_skipped=0, primary_mapped=0)
meta.update({k: int(v) for k, v in qc.items()})
json.dump(meta, open(os.path.join(OUT, f'meta_multiL_{sid}.json'), 'w'), indent=1, ensure_ascii=False)
print(json.dumps(meta, ensure_ascii=False))
