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
import sys, os, re, collections
import numpy as np, pysam
sid, bam_path, outroot = sys.argv[1:4]; OFFS = 12; TRIM = 45; REFC = (0, 2, 7, 8)
GTF = PEDV_SERVER_ROOT+"/refdata_dual_host/vero_ChlSab1.1/Chlorocebus_sabaeus.ChlSab1.1.116.gtf"
rx = re.compile(r'(\S+) "([^"]*)"'); cds = collections.defaultdict(list)
for line in open(GTF):
    if line.startswith("#"): continue
    f = line.split("\t")
    if f[2] != "CDS" or 'tag "Ensembl_canonical"' not in f[8]: continue
    d = dict(rx.findall(f[8]))
    if d.get("transcript_biotype") != "protein_coding": continue
    cds[d["transcript_id"]].append((f[0], f[6], int(f[3]), int(f[4])))
blocks = collections.defaultdict(list)
for t, bl in cds.items():
    chrom, strand = bl[0][0], bl[0][1]; bl = sorted(bl, key=lambda x: x[2], reverse=(strand == "-"))
    tlen = sum(e - b + 1 for _, _, b, e in bl); c = 0
    for _, _, b, e in bl: blocks[(chrom, strand)].append((b, e, c, tlen)); c += e - b + 1
HB = {k: np.array(sorted(v), dtype=np.int64) for k, v in blocks.items()}
F = np.zeros(3, np.int64); FR = collections.Counter(); n_used = 0
bam = pysam.AlignmentFile(bam_path)
for chrom in bam.references:
    if chrom.startswith("PEDV"): continue
    P = {"+": [], "-": []}
    for r in bam.fetch(chrom):
        if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
        if r.has_tag("NH") and r.get_tag("NH") != 1: continue
        cig = r.cigartuples; rev = r.is_reverse
        if any(op == 3 for op, n in cig):
            co = []; x = r.reference_start
            for op, n in cig:
                if op in REFC: co.extend(range(x + 1, x + n + 1)); x += n
                elif op == 3: x += n
            if rev: co.reverse()
            if OFFS >= len(co): continue
            p = co[OFFS]
        else:
            p = (r.reference_end - OFFS) if rev else (r.reference_start + 1 + OFFS)
        P["-" if rev else "+"].append(p)
    for st in "+-":
        a = HB.get((chrom, st))
        if a is None or not P[st]: continue
        pos = np.asarray(P[st], np.int64)
        i = np.searchsorted(a[:, 0], pos, side="right") - 1; ok = i >= 0; ii = np.where(ok, i, 0)
        gs, ge, cb, tl = a[ii, 0], a[ii, 1], a[ii, 2], a[ii, 3]; ok &= pos <= ge
        c = np.where(st == "+", cb + (pos - gs), cb + (ge - pos))
        ok &= (c >= TRIM) & (c < tl - TRIM)
        np.add.at(F, c[ok] % 3, 1); n_used += int(ok.sum())
bam.close()
os.makedirs(outroot, exist_ok=True)
with open(f"{outroot}/{sid}.rna_pseudo_psite_frame.tsv", "w") as f:
    f.write("sample_id\tpseudo_offset\tn_reads_in_CDS_body\tf0\tf1\tf2\n")
    f.write(f"{sid}\t{OFFS}\t{n_used}\t{F[0]}\t{F[1]}\t{F[2]}\n")
print(sid, "DONE", F)
