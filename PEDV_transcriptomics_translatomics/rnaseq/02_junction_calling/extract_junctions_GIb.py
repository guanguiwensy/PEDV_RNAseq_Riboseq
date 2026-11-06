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
import sys
import os
import glob
import pysam
import csv

MIN_GAP = 50


def get_tag(read, tag, default=None):
    try:
        return read.get_tag(tag)
    except KeyError:
        return default


def extract_one(bam_path, sample, out_unfiltered):
    bam = pysam.AlignmentFile(bam_path, "rb")
    unfiltered_rows = []
    n_total_reads = 0
    n_with_n_op = 0

    for read in bam.fetch(until_eof=True):
        n_total_reads += 1
        if read.is_unmapped or read.cigartuples is None:
            continue
        has_n = any(op == 3 for op, _ in read.cigartuples)
        if not has_n:
            continue
        n_with_n_op += 1

        nh = get_tag(read, "NH", 1)
        is_dup = read.is_duplicate
        is_sec = read.is_secondary
        is_sup = read.is_supplementary
        mapq = read.mapping_quality
        chrom = read.reference_name
        seq = read.query_sequence if read.query_sequence else ""

        cigartuples = read.cigartuples
        ref_pos = read.reference_start
        block_start_1based = read.reference_start + 1
        gap_index_in_read = 0

        for cig_idx, (op, length) in enumerate(cigartuples):
            if op in (0, 7, 8):
                ref_pos += length
            elif op == 1:
                pass
            elif op == 2:
                ref_pos += length
            elif op == 3:
                donor_0based = ref_pos
                gap_len = length
                acceptor_0based = ref_pos + length
                if gap_len >= MIN_GAP:
                    gap_index_in_read += 1
                    leader_end_1based = donor_0based
                    body_start_1based = acceptor_0based + 1
                    effective_leader_start_1based = block_start_1based

                    ref_len_after = 0
                    for op2, length2 in cigartuples[cig_idx + 1:]:
                        if op2 in (0, 2, 7, 8):
                            ref_len_after += length2
                        elif op2 in (3, 4, 5):
                            break
                    body_end_1based = body_start_1based + ref_len_after - 1 if ref_len_after else body_start_1based

                    unfiltered_rows.append({
                        "sample": sample, "chrom": chrom,
                        "read_name": read.query_name,
                        "leader_start": effective_leader_start_1based,
                        "leader_end": leader_end_1based,
                        "body_start": body_start_1based,
                        "body_end": body_end_1based,
                        "gap_len": gap_len,
                        "gap_index_in_read": gap_index_in_read,
                        "cigar": read.cigarstring,
                        "mapq": mapq, "NH": nh,
                        "is_duplicate": is_dup, "is_secondary": is_sec, "is_supplementary": is_sup,
                        "is_read1": read.is_read1, "is_read2": read.is_read2,
                        "junction_sequence": seq,
                    })
                ref_pos = acceptor_0based
                block_start_1based = acceptor_0based + 1
            elif op in (4, 5):
                pass

    bam.close()

    fieldnames = ["sample", "chrom", "read_name", "leader_start", "leader_end",
                  "body_start", "body_end", "gap_len", "gap_index_in_read", "cigar", "mapq", "NH",
                  "is_duplicate", "is_secondary", "is_supplementary", "is_read1", "is_read2",
                  "junction_sequence"]
    with open(out_unfiltered, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for row in unfiltered_rows:
            w.writerow(row)

    return n_total_reads, n_with_n_op, len(unfiltered_rows)


def main():
    mm2dir = sys.argv[1]
    outdir = sys.argv[2]
    samples_file = sys.argv[3]
    os.makedirs(outdir, exist_ok=True)
    with open(samples_file) as f:
        samples = [l.strip() for l in f.readlines()[1:] if l.strip()]

    for sample in samples:
        bam_path = os.path.join(mm2dir, f"{sample}_mm2.sorted.bam")
        out_fp = os.path.join(outdir, f"{sample}_mm2_unfiltered.tsv")
        stats = extract_one(bam_path, sample, out_fp)
        print(f"{sample}\ttotal_reads={stats[0]}\treads_with_N={stats[1]}\tunfiltered_rows={stats[2]}")
    print("GIB_GAPFIX_EXTRACT_ALL_DONE")


if __name__ == "__main__":
    main()
