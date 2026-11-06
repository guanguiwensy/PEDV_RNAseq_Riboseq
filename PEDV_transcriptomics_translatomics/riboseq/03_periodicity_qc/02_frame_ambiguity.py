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
import os, gzip, re, collections
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.dirname(HERE); TASK = os.path.dirname(OUT)
GTF = os.path.join(os.path.dirname(TASK), "refdata_vero", "Chlorocebus_sabaeus.ChlSab1.1.116.gtf.gz")
TAB = f"{OUT}/tables"
rx = re.compile(r'(\S+) "([^"]*)"')
feats = []   # chrom, strand, start, end, phase, transcript, canonical
canon = {}
n_pc = 0
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if f[2] != "CDS": continue
        d = dict(rx.findall(f[8]))
        if d.get("transcript_biotype", d.get("gene_biotype")) != "protein_coding" and d.get("gene_biotype") != "protein_coding": continue
        is_c = 'tag "Ensembl_canonical"' in f[8]
        feats.append((f[0], f[6], int(f[3]), int(f[4]), int(f[7]) if f[7] != "." else 0, d["transcript_id"], d["gene_id"], is_c))
F = pd.DataFrame(feats, columns=["chrom", "strand", "start", "end", "phase", "tid", "gid", "canon"])
print("CDS features:", len(F), "transcripts:", F.tid.nunique(), "canonical transcripts:", F[F.canon].tid.nunique())

def cp_at(pos, r):
    i = (pos - r.start) if r.strand == "+" else (r.end - pos)
    return (i - r.phase) % 3
amb = collections.Counter(); tot = collections.Counter()
for t, g in F[F.canon].groupby("tid"): tot[t] = int((g.end - g.start + 1).sum())
can = F[F.canon]; oth = F[~F.canon]
for (chrom, strand), cg in can.groupby(["chrom", "strand"]):
    og = F[(F.chrom == chrom) & (F.strand == strand)].sort_values("start").reset_index(drop=True)
    starts = og.start.values; ends = og.end.values
    cummax_end = np.maximum.accumulate(ends)
    for r in cg.itertuples():

        hi = np.searchsorted(starts, r.end, side="right")

        lo = np.searchsorted(cummax_end, r.start, side="left")
        if lo >= hi: continue
        mask = np.zeros(r.end - r.start + 1, dtype=bool)
        for j in range(lo, hi):
            o = og.iloc[j]
            if o.tid == r.tid or ends[j] < r.start: continue
            a, b = max(r.start, o.start), min(r.end, o.end)
            if a > b: continue
            if (cp_at(a, r) != cp_at(a, o)):
                mask[a - r.start:b - r.start + 1] = True
        amb[r.tid] += int(mask.sum())
rows = [dict(transcript_id=t, cds_nt_annotated=tot[t], ambiguous_nt=amb.get(t, 0), ambiguous_fraction=amb.get(t, 0) / tot[t]) for t in tot]
A = pd.DataFrame(rows)
A.to_csv(f"{TAB}/canonical_cds_frame_ambiguity.csv", index=False)
print(A.ambiguous_fraction.describe())
print("transcripts with ambiguous_fraction>0:", int((A.ambiguous_fraction > 0).sum()), "of", len(A), "; >=0.10:", int((A.ambiguous_fraction >= 0.10).sum()))
