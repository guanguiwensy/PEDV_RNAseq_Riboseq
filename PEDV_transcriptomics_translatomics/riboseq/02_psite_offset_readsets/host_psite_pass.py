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
import sys, os, re, json, collections, time
import numpy as np, pysam
sid, bam_path, offj, outroot = sys.argv[1:5]
OUT = f"{outroot}/{sid}"; os.makedirs(OUT, exist_ok=True)
GTF = PEDV_SERVER_ROOT+"/refdata_dual_host/vero_ChlSab1.1/Chlorocebus_sabaeus.ChlSab1.1.116.gtf"
REFC = (0, 2, 7, 8); W = 60; LMAX = 150
OFF = {int(k): v for k, v in json.load(open(offj))[sid].items() if v is not None}
t0 = time.time()

rx = re.compile(r'(\S+) "([^"]*)"')
ex = collections.defaultdict(list); cds = collections.defaultdict(list); tinfo = {}
for line in open(GTF):
    if line.startswith("#"): continue
    f = line.split("\t")
    if f[2] not in ("exon", "CDS") or 'tag "Ensembl_canonical"' not in f[8]: continue
    d = dict(rx.findall(f[8]))
    if d.get("transcript_biotype") != "protein_coding": continue
    t = d["transcript_id"]; tinfo[t] = (f[0], f[6], d["gene_id"], d.get("gene_name", ""))
    (ex if f[2] == "exon" else cds)[t].append((int(f[3]), int(f[4])))

blocks = collections.defaultdict(list)   # (chrom,strand) -> [gs, ge, tx_base, tid_index]
TID = []; CS = []; CE = []; TL = []
for t, el in ex.items():
    if t not in cds: continue
    chrom, strand, g, gn = tinfo[t]
    el = sorted(el, key=lambda x: x[0], reverse=(strand == "-"))
    def tx(pos):
        c = 0
        for b, e in el:
            if b <= pos <= e: return c + (pos - b if strand == "+" else e - pos)
            c += e - b + 1
        return None
    cl = cds[t]
    cstart_g = min(b for b, e in cl) if strand == "+" else max(e for b, e in cl)
    cend_g = max(e for b, e in cl) if strand == "+" else min(b for b, e in cl)
    cs, ce = tx(cstart_g), tx(cend_g)
    if cs is None or ce is None: continue
    k = len(TID); TID.append((t, g, gn)); CS.append(cs); CE.append(ce); TL.append(sum(e - b + 1 for b, e in el))
    c = 0
    for b, e in el:
        blocks[(chrom, strand)].append((b, e, c, k)); c += e - b + 1
HB = {}
for kk, v in blocks.items():
    v.sort(); HB[kk] = np.array(v, dtype=np.int64)
CS = np.array(CS); CE = np.array(CE)
GC = np.zeros((len(TID), LMAX + 1, 3), np.int64)

def to_tx(chrom, strand, pos):
    a = HB.get((chrom, strand))
    if a is None or not len(pos): return np.full(len(pos), -1), np.full(len(pos), -1)
    i = np.searchsorted(a[:, 0], pos, side="right") - 1
    ok = i >= 0; ii = np.where(ok, i, 0)
    gs, ge, cb, k = a[ii, 0], a[ii, 1], a[ii, 2], a[ii, 3]
    ok &= pos <= ge
    c = np.where(strand == "+", cb + (pos - gs), cb + (ge - pos))
    return np.where(ok, c, -1), np.where(ok, k, -1)

MS = np.zeros((LMAX + 1, 2 * W + 1), np.int64); ME = np.zeros((LMAX + 1, 2 * W + 1), np.int64)

nNA = 0; nused = 0
bam = pysam.AlignmentFile(bam_path)
for chrom in bam.references:
    if chrom.startswith("PEDV"): continue
    P = {"+": [], "-": []}; LL = {"+": [], "-": []}
    for r in bam.fetch(chrom):
        if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
        cig = r.cigartuples; rev = r.is_reverse
        L0 = sum(n for op, n in cig if op in REFC); L = min(L0, LMAX)
        o = OFF.get(L)
        if o is None: nNA += 1; continue
        if any(op == 3 for op, n in cig):
            coords = []; x = r.reference_start
            for op, n in cig:
                if op in REFC: coords.extend(range(x + 1, x + n + 1)); x += n
                elif op == 3: x += n
            if rev: coords.reverse()
            if o >= len(coords): nNA += 1; continue
            p = coords[o]
        else:
            p = (r.reference_end - o) if rev else (r.reference_start + 1 + o)
        st = "-" if rev else "+"; P[st].append(p); LL[st].append(L); nused += 1
    for st in "+-":
        if not P[st]: continue
        p = np.asarray(P[st], np.int64); ln = np.asarray(LL[st], np.int64)
        c, k = to_tx(chrom, st, p)
        ok = k >= 0
        c, k, ln = c[ok], k[ok], ln[ok]
        rs = c - CS[k]; re_ = c - (CE[k] + 1)
        m = (rs >= -W) & (rs <= W); np.add.at(MS, (ln[m], rs[m] + W), 1)
        m = (re_ >= -W) & (re_ <= W); np.add.at(ME, (ln[m], re_[m] + W), 1)
        incds = (rs >= 0) & (c <= CE[k])
        np.add.at(GC, (k[incds], ln[incds], rs[incds] % 3), 1)
bam.close()
with open(f"{OUT}/host_psite_metagene.tsv", "w") as f:
    f.write("anchor\tL\trel\tn\n")
    for an, M in (("start", MS), ("stop", ME)):
        for L in range(LMAX + 1):
            for i in np.nonzero(M[L])[0]: f.write(f"{an}\t{L}\t{int(i)-W}\t{M[L,i]}\n")
with open(f"{OUT}/host_gene_cds_psite.tsv", "w") as f:
    f.write("transcript_id\tgene_id\tgene_name\tcds_len_nt\tL\tf0\tf1\tf2\n")
    for kk in np.nonzero(GC.sum((1, 2)))[0]:
        M = GC[kk]; t, g, gn = TID[kk]
        for L in np.nonzero(M.sum(1))[0]:
            f.write(f"{t}\t{g}\t{gn}\t{CE[kk]-CS[kk]+1}\t{L}\t{M[L,0]}\t{M[L,1]}\t{M[L,2]}\n")
with open(f"{OUT}/host_psite_qc.tsv", "w") as f:
    f.write(f"bam\t{bam_path}\nreads_with_psite\t{nused}\nreads_offset_NA_or_short\t{nNA}\ncanonical_pc_transcripts\t{len(TID)}\nelapsed_s\t{round(time.time()-t0)}\n")
print(sid, "DONE", nused, nNA, round(time.time() - t0))
