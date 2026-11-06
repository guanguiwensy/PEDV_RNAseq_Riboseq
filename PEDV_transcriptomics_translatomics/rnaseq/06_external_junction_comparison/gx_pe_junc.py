#!/usr/bin/env python3
"""Junction extraction for paired-end data mapped as one fastq (both mates, names suffixed /1 /2) with minimap2 splice. Same rules as phase 2/3 (N-gap>=50, MAPQ>=20, anchors 15/15, primary only, left-aligned donor/acceptor)
plus PE handling: a fragment crossing the same (donor, acceptor) in both mates is counted once; viral fragments = fragments with >=1 read >=50% aligned. No PCR de-duplication (as original).
usage: gx_pe_junc.py run reflabel (派生自 phase5/scripts/p5_pe_junc.py，仅改路径；分类/过滤规则零改动)"""
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
import sys, os, json, gzip, collections
import numpy as np, pysam
sys.path.insert(0, REPO_ROOT + '/rnaseq/06_external_junction_comparison')
from p3_lib import cigar_map, read_fa, CANON, LEADER_WIN
W = PEDV_SERVER_ROOT+"/public_replication/ESV/gII_cross_strain_validation_v1"; PH5 = PEDV_SERVER_ROOT+"/public_replication/ESV/phase5"
run, lab = sys.argv[1:3]
ref = f"{W}/_work/ref/{lab}.fa"; gib = f"{PH5}/_work/ref/PEDV_GIb.fa"
mp = f"{W}/_work/map_{lab}.npz"
if not os.path.exists(mp):
    q2g, g2q = cigar_map(ref, gib); np.savez(mp, q2g=q2g, g2q=g2q)
z = np.load(mp); q2g = z["q2g"]; L = len(read_fa(ref))
bam = pysam.AlignmentFile(f"{W}/_work/bam/{run}.{lab}.bam", "rb")
depth = np.zeros(L + 2, dtype=np.int32)
vfrag = set(); n_gap = 0; seen = set(); rows = []
def cls(dg, ag):
    if dg == 0 or ag == 0: return "nonhomologous_end"
    indep = not (LEADER_WIN[0] <= dg <= LEADER_WIN[1])
    can = [g for g, p in CANON.items() if abs(ag - p) <= 5]
    if not indep: return "TRSL_dep_canonical_" + can[0] if can else "TRSL_dep_noncanonical"
    return "TRSL_indep_distant" if abs(ag - dg) >= 1000 else "TRSL_indep_local"
for r in bam.fetch(until_eof=True):
    if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
    base = r.query_name.rsplit("/", 1)[0]
    al = sum(l for o, l in r.cigartuples if o in (0, 7, 8))
    if al >= 0.5 * r.query_length: vfrag.add(base)
    for s, e in r.get_blocks(): depth[s + 1:e + 1] += 1
    ct = r.cigartuples; pos = r.reference_start; la = 0
    for i, (o, l) in enumerate(ct):
        if o in (0, 7, 8): pos += l; la += l
        elif o == 2: pos += l
        elif o == 3:
            if l >= 50:
                n_gap += 1; d = pos; a = pos + l + 1; ba = 0
                for o2, l2 in ct[i + 1:]:
                    if o2 in (0, 7, 8): ba += l2
                    elif o2 == 3: break
                if r.mapping_quality >= 20 and la >= 15 and ba >= 15 and (base, d, a) not in seen:
                    seen.add((base, d, a))
                    rows.append((base, d, a, l, la, ba, r.reference_start + 1, int(q2g[d]) if d < len(q2g) else 0, int(q2g[a]) if a < len(q2g) else 0))
            pos += l
os.makedirs(f"{W}/junc", exist_ok=True); os.makedirs(f"{W}/cov", exist_ok=True)
np.save(f"{W}/cov/{run}.{lab}.depth.npy", depth)
with gzip.open(f"{W}/junc/{run}.{lab}.reads.tsv.gz", "wt") as f:
    f.write("read\tdonor\tacceptor\tgap\tleader_anchor\tbody_anchor\tblock_start\tdonor_GIb\tacceptor_GIb\tclass\n")
    for x in rows: f.write("\t".join(map(str, x)) + "\t" + cls(x[7], x[8]) + "\n")
json.dump(dict(run=run, ref=lab, primary_reads_in_prefilter=len(vfrag), reads_aln_any=len(vfrag), reads_aln_ge50pct=len(vfrag), gap_rows_ge50=n_gap, qc_pass_junction_reads=len(rows)), open(f"{W}/junc/{run}.{lab}.summary.json", "w"))
print(run, lab, len(vfrag), n_gap, len(rows))
