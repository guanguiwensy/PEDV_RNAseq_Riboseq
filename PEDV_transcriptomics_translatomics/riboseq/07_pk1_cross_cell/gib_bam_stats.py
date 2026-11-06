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
import sys, os, re, gzip, collections, time
import numpy as np
import pysam

sid, bam_path, vchr, outroot, mode = sys.argv[1:6]
OUT = f"{outroot}/{sid}"; os.makedirs(OUT, exist_ok=True)
GTF = PEDV_SERVER_ROOT+"/refdata_dual_host/pig_Sscrofa11.1/Sus_scrofa.Sscrofa11.1.110.gtf"

START_BED = sys.argv[6] if len(sys.argv) > 6 else PEDV_SERVER_ROOT+"/analysis_PK1/psite/annotation/start_codon_pig.bed"
ORF_TSV = PEDV_SERVER_ROOT+"/refdata_GIb_v1.1/virus/PEDV_GIb.orf_check.tsv"
TRIM = 45; LMAX = 150; CLEAN_NT = 22; NFIRST = 25; META_W = 300
REFC = (0, 2, 7, 8)   # M D = X
FULL = (mode == "full")
t0 = time.time()

# ---------- PEDV ORF ----------
with pysam.AlignmentFile(bam_path) as f:
    vlen = f.get_reference_length(vchr)
orfs = []
if vchr == "PEDV_GIb":
    for r in open(ORF_TSV):
        x = r.rstrip("\n").split("\t")
        if x[0] == "orf": continue
        orfs.append((x[0], int(x[1]), int(x[2])))
ONAMES = [o[0] for o in orfs]
cover = np.zeros(vlen + 2, np.int16)
for _, a, b in orfs: cover[a:b + 1] += 1
vframe = np.full(vlen + 2, -1, np.int8); vorf = np.full(vlen + 2, -1, np.int8)
for k, (_, a, b) in enumerate(orfs):
    p = np.arange(a + TRIM, b - TRIM + 1); p = p[cover[p] == 1]
    vframe[p] = (p - a) % 3; vorf[p] = k
vstarts = [(n, a) for n, a, b in orfs if n != "ORF1b"]


HB = {}
if FULL:
    rx = re.compile(r'(\S+) "([^"]*)"')
    cds = collections.defaultdict(list)
    for line in open(GTF):
        if line.startswith("#"): continue
        f = line.split("\t")
        if f[2] != "CDS" or 'tag "Ensembl_canonical"' not in f[8]: continue
        d = dict(rx.findall(f[8]))
        if d.get("transcript_biotype") != "protein_coding": continue
        cds[d["transcript_id"]].append((f[0], f[6], int(f[3]), int(f[4])))
    blocks = collections.defaultdict(list)
    for t, bl in cds.items():
        chrom, strand = bl[0][0], bl[0][1]
        bl = sorted(bl, key=lambda x: x[2], reverse=(strand == "-"))
        tlen = sum(e - b + 1 for _, _, b, e in bl); c = 0
        for _, _, b, e in bl:
            blocks[(chrom, strand)].append((b, e, c, tlen)); c += e - b + 1
    for k, v in blocks.items():
        v.sort(); HB[k] = np.array(v, dtype=np.int64)

def host_cds_coord(chrom, strand, pos):
    a = HB.get((chrom, strand))
    out = np.full(pos.shape, -1, np.int64); tl = np.zeros(pos.shape, np.int64)
    if a is None or not len(pos): return out, tl
    i = np.searchsorted(a[:, 0], pos, side="right") - 1
    ok = i >= 0; ii = np.where(ok, i, 0)
    gs, ge, cb, tlen = a[ii, 0], a[ii, 1], a[ii, 2], a[ii, 3]
    ok &= pos <= ge
    c = np.where(strand == "+", cb + (pos - gs), cb + (ge - pos))
    out[ok] = c[ok]; tl[ok] = tlen[ok]
    return out, tl


rl_host = np.zeros(LMAX + 1, np.int64); rl_vp = np.zeros(LMAX + 1, np.int64)
rl_vm = np.zeros(LMAX + 1, np.int64); rl_vpc = np.zeros(LMAX + 1, np.int64)
F5 = collections.defaultdict(lambda: np.zeros((LMAX + 1, 3), np.int64))
META5 = {"start": np.zeros((LMAX + 1, META_W), np.int64), "stop": np.zeros((LMAX + 1, META_W), np.int64)}
qc = collections.Counter()
AM = collections.defaultdict(collections.Counter)    # domain -> metric -> value
DUP = collections.Counter()
M5 = collections.Counter(); SP = collections.Counter()

def cig_info(cig, rev):
    L = 0; first_n = None; ins = dele = 0; spliced = False
    for op, n in (reversed(cig) if rev else cig):
        if op in REFC:
            L += n
            if op == 2: dele += 1
        elif op == 3:
            spliced = True
            if first_n is None: first_n = L
        elif op == 1: ins += 1
    clean = (first_n is None) or first_n >= CLEAN_NT
    sc_l = cig[0][1] if cig[0][0] == 4 else 0; sc_r = cig[-1][1] if cig[-1][0] == 4 else 0
    sc5 = sc_r if rev else sc_l
    mbases = sum(n for op, n in cig if op in (0, 7, 8))
    return L, clean, sc5, sc_l + sc_r, (ins + dele) > 0, spliced, mbases

def add_am(dom, r, L, sct, ind, spl, mb, sc5):
    a = AM[dom]; nm = r.get_tag("nM") if r.has_tag("nM") else 0
    a["reads"] += 1; a["aligned_M_bases"] += mb; a["nM_sum"] += nm; a["reads_with_mismatch"] += (nm > 0)
    a["softclip_reads"] += (sct > 0); a["softclip_bases"] += sct; a["softclip5_reads"] += (sc5 > 0)
    a["indel_reads"] += ind; a["spliced_reads"] += spl
    nh = r.get_tag("NH") if r.has_tag("NH") else 1
    a["NH_gt1_reads"] += (nh > 1)

bam = pysam.AlignmentFile(bam_path)
refs = bam.references if FULL else (vchr,)
for chrom in refs:
    pos = {"+": [], "-": []}; lens = {"+": [], "-": []}
    keys = set()
    isv = chrom == vchr
    for r in bam.fetch(chrom):
        if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
        rev = r.is_reverse; cig = r.cigartuples
        L0, clean, sc5, sct, ind, spl, mb = cig_info(cig, rev)
        L = min(L0, LMAX)
        p5 = r.reference_end if rev else r.reference_start + 1
        st = "-" if rev else "+"
        keys.add((st, p5, L0))
        if isv:
            add_am("pedv_all", r, L0, sct, ind, spl, mb, sc5)
            add_am("pedv_minus" if rev else "pedv_plus", r, L0, sct, ind, spl, mb, sc5)

            if not clean:
                coords = []; x = r.reference_start
                for op, n in cig:
                    if op in REFC: coords.extend(range(x + 1, x + n + 1)); x += n
                    elif op == 3: x += n
                if rev: coords.reverse()
                SP[(st, L0, min(sc5, 30), ",".join(map(str, coords[:NFIRST])))] += 1
            else:
                M5[(st, L0, p5, min(sc5, 30))] += 1
            if rev: rl_vm[L] += 1; continue
            rl_vp[L] += 1
            if not clean: qc["pedv_plus_unclean_5p_junction"] += 1; continue
            rl_vpc[L] += 1
            fr = vframe[p5] if p5 <= vlen else -1
            if fr >= 0:
                F5["PEDV_all"][L, fr] += 1; F5["PEDV_" + ONAMES[vorf[p5]]][L, fr] += 1
            else: qc["pedv_plus_5p_outside_unique_orf_body"] += 1
            continue
        add_am("host", r, L0, sct, ind, spl, mb, sc5)
        rl_host[L] += 1
        if not clean: qc["host_unclean"] += 1; continue
        pos[st].append(p5); lens[st].append(L)
    if isv:
        DUP[("pedv_all", "distinct")] += len(keys)
        DUP[("pedv_plus", "distinct")] += sum(1 for k in keys if k[0] == "+")
        DUP[("pedv_minus", "distinct")] += sum(1 for k in keys if k[0] == "-")
    else:
        DUP[("host", "distinct")] += len(keys)
    for s_ in "+-":
        if not pos[s_]: continue
        p = np.asarray(pos[s_], np.int64); ln = np.asarray(lens[s_], np.int64)
        c, tl = host_cds_coord(chrom, s_, p)
        ok = (c >= TRIM) & (c < tl - TRIM)
        np.add.at(F5["host"], (ln[ok], c[ok] % 3), 1)
        okS = (c >= 0) & (c < META_W); np.add.at(META5["start"], (ln[okS], c[okS]), 1)
        rel = c - tl; okE = (c >= 0) & (rel >= -META_W); np.add.at(META5["stop"], (ln[okE], rel[okE] + META_W), 1)
bam.close()
print(sid, "pass1 done", round(time.time() - t0), "s", flush=True)

def pos_in_read(r, A):
    coords = []; x = r.reference_start
    for op, n in r.cigartuples:
        if op in REFC: coords.extend(range(x + 1, x + n + 1)); x += n
        elif op == 3: x += n
    if r.is_reverse: coords.reverse()
    try: return len(coords), coords.index(A) + 1
    except ValueError: return len(coords), None

HS = collections.Counter(); HH = collections.Counter(); VS = collections.Counter()
bam = pysam.AlignmentFile(bam_path)
nstart = 0
if FULL:
    ss = set(); seen = []
    for line in open(START_BED):
        c, s0, e, t, _, st = line.rstrip("\n").split("\t"); A = int(e)
        if (c, A, st) in ss: continue
        ss.add((c, A, st)); seen.append((c, A, st))
    nstart = len(seen)
    for k, (c, A, st) in enumerate(seen):
        h = k % 2
        for r in bam.fetch(c, A - 1, A):
            if r.is_secondary or r.is_supplementary or (r.is_reverse != (st == "-")): continue
            L, j = pos_in_read(r, A)
            if j: HS[(min(L, LMAX), j)] += 1; HH[(h, min(L, LMAX), j)] += 1
for name, A in vstarts:
    for r in bam.fetch(vchr, A - 1, A):
        if r.is_secondary or r.is_supplementary or r.is_reverse: continue
        L, j = pos_in_read(r, A)
        if j: VS[(name, min(L, LMAX), j)] += 1
bam.close()


def w(name, header, rows):
    op = gzip.open if name.endswith(".gz") else open
    with op(f"{OUT}/{name}", "wt") as f:
        f.write(header + "\n")
        for r in rows: f.write("\t".join(map(str, r)) + "\n")
w("readlen_aln.tsv", "L\thost_n\tpedv_plus_n\tpedv_minus_n\tpedv_plus_clean_n",
  [(L, rl_host[L], rl_vp[L], rl_vm[L], rl_vpc[L]) for L in range(LMAX + 1) if rl_host[L] + rl_vp[L] + rl_vm[L]])
if FULL:
    w("host_startpos.tsv", "L\tj\tn", [(L, j, n) for (L, j), n in sorted(HS.items())])
    w("host_startpos_split.tsv", "half\tL\tj\tn", [(h, L, j, n) for (h, L, j), n in sorted(HH.items())])
    rows = []
    for an, M in META5.items():
        for L in range(LMAX + 1):
            for i in np.nonzero(M[L])[0]:
                rows.append((an, L, int(i) if an == "start" else int(i) - META_W, M[L, i]))
    w("host_cds_meta5.tsv", "anchor\tL\trel\tn", rows)
w("pedv_startpos.tsv", "orf\tL\tj\tn", [(o, L, j, n) for (o, L, j), n in sorted(VS.items())])
w("frame5.tsv", "domain\tL\tf0\tf1\tf2",
  [(d, L, M[L, 0], M[L, 1], M[L, 2]) for d, M in sorted(F5.items()) for L in range(LMAX + 1) if M[L].sum()])
w("pedv_5p_matrix.tsv.gz", "strand\tL\tpos5\tsc5\tn", [(*k, n) for k, n in sorted(M5.items())])
w("pedv_spliced5p.tsv.gz", "strand\tL\tsc5\tcoords\tn", [(*k, n) for k, n in sorted(SP.items())])
MK = ["reads", "aligned_M_bases", "nM_sum", "reads_with_mismatch", "softclip_reads", "softclip_bases",
      "softclip5_reads", "indel_reads", "spliced_reads", "NH_gt1_reads"]
w("aln_metrics.tsv", "domain\t" + "\t".join(MK), [(d, *[AM[d][k] for k in MK]) for d in sorted(AM)])
w("dup.tsv", "domain\treads\tdistinct_5p_keys",
  [(d, AM[d]["reads"], DUP[(d, "distinct")]) for d in sorted(AM)])
with open(f"{OUT}/qc.tsv", "w") as f:
    f.write(f"bam\t{bam_path}\nbam_size\t{os.path.getsize(bam_path)}\nmode\t{mode}\nvirus_contig\t{vchr}\tlen={vlen}\n")
    f.write(f"host_total\t{rl_host.sum()}\npedv_plus_total\t{rl_vp.sum()}\npedv_minus_total\t{rl_vm.sum()}\n")
    f.write(f"host_start_codons_used\t{nstart}\npedv_orfs\t{';'.join(f'{n}:{a}-{b}' for n,a,b in orfs)}\n")
    for k, v in qc.items(): f.write(f"{k}\t{v}\n")
    f.write(f"elapsed_s\t{round(time.time()-t0)}\n")
print(sid, "DONE", round(time.time() - t0), "s host", rl_host.sum(), "pedv+", rl_vp.sum(), "pedv-", rl_vm.sum(), dict(qc), flush=True)
