#!/usr/bin/env python3
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
import sample_registry as SR
import csv
import os
import re
import sys

import pysam

BASE = PEDV_SERVER_ROOT+"/analysis_GIb/RNAseq"
VIRUS_CONTIG = "PEDV_GIb"
VIRUS_LEN = 27959


def parse_log_final(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            if "|" not in line:
                continue
            k, v = line.split("|", 1)
            d[k.strip()] = v.strip()
    return d


def virus_host_counts_and_mismatch(bam_path):
    virus_frag = 0
    host_frag = 0
    virus_nm_sum = 0
    virus_alen_sum = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            try:
                nh = read.get_tag("NH")
            except KeyError:
                nh = 1
            if nh != 1:
                continue
            is_virus = read.reference_name == VIRUS_CONTIG
            if read.is_read1 or (not read.is_paired):
                if is_virus:
                    virus_frag += 1
                else:
                    host_frag += 1
            if is_virus:
                try:
                    nm = read.get_tag("NM")
                except KeyError:
                    nm = 0
                virus_nm_sum += nm
                virus_alen_sum += read.query_alignment_length or 0
    virus_mismatch_pct = (100.0 * virus_nm_sum / virus_alen_sum) if virus_alen_sum > 0 else float("nan")
    return virus_frag, host_frag, virus_mismatch_pct


def coverage_stats(depth_tsv):
    if not os.path.exists(depth_tsv):
        return {"mean_depth": float("nan"), "breadth_1x_pct": float("nan"), "breadth_10x_pct": float("nan"), "n_positions": 0}
    n = 0
    total = 0
    n_ge1 = 0
    n_ge10 = 0
    with open(depth_tsv) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            d = int(parts[2])
            n += 1
            total += d
            if d >= 1:
                n_ge1 += 1
            if d >= 10:
                n_ge10 += 1
    if n == 0:
        return {"mean_depth": float("nan"), "breadth_1x_pct": float("nan"), "breadth_10x_pct": float("nan"), "n_positions": 0}
    return {
        "mean_depth": total / n,
        "breadth_1x_pct": 100.0 * n_ge1 / n,
        "breadth_10x_pct": 100.0 * n_ge10 / n,
        "n_positions": n,
    }


def main():
    manifest = sys.argv[1] if len(sys.argv) > 1 else os.path.join(BASE, "scripts", "samples.tsv")
    out_tsv = sys.argv[2] if len(sys.argv) > 2 else os.path.join(BASE, "alignment_summary.tsv")

    with open(manifest) as fh:
        samples = [l.strip() for l in fh.readlines()[1:] if l.strip()]
    samples = [s for s in samples if s in SR.WHITELIST_LIBRARIES]   # final sample table only

    rows = []
    for sample in samples:
        outdir = os.path.join(BASE, "align", sample)
        log_final = os.path.join(outdir, "Log.final.out")
        bam = os.path.join(outdir, "Aligned.sorted.bam")
        depth_tsv = os.path.join(BASE, "coverage", f"{sample}.depth.tsv")

        row = {"sample_id": sample}
        if not os.path.exists(log_final):
            row["status"] = "MISSING_LOG"
            rows.append(row)
            print(f"[WARN] {sample}: Log.final.out missing", file=sys.stderr)
            continue
        log = parse_log_final(log_final)
        total_reads = int(log.get("Number of input reads", "0"))
        uniq_pct_overall = log.get("Uniquely mapped reads %", "NA").rstrip("%")
        mismatch_pct_overall = log.get("Mismatch rate per base, %", "NA").rstrip("%")
        multi_pct = log.get("% of reads mapped to multiple loci", "NA").rstrip("%")
        too_many_loci_pct = log.get("% of reads mapped to too many loci", "NA").rstrip("%")
        unmapped_short_pct = log.get("% of reads unmapped: too short", "NA").rstrip("%")

        if os.path.exists(bam):
            virus_frag, host_frag, virus_mismatch_pct = virus_host_counts_and_mismatch(bam)
            denom = virus_frag + host_frag
            virus_fraction_pct = (100.0 * virus_frag / denom) if denom > 0 else float("nan")
        else:
            virus_frag = host_frag = virus_fraction_pct = virus_mismatch_pct = float("nan")
            print(f"[WARN] {sample}: Aligned.sorted.bam missing", file=sys.stderr)

        cov = coverage_stats(depth_tsv)

        row.update({
            "status": "OK",
            "total_input_read_pairs": total_reads,
            "host_unique_fragments": host_frag,
            "virus_unique_fragments": virus_frag,
            "virus_fraction_pct": round(virus_fraction_pct, 4) if virus_fraction_pct == virus_fraction_pct else "NA",
            "uniquely_mapped_pct_overall": uniq_pct_overall,
            "multi_mapped_pct_overall": multi_pct,
            "too_many_loci_pct_overall": too_many_loci_pct,
            "unmapped_too_short_pct_overall": unmapped_short_pct,
            "mismatch_rate_pct_overall": mismatch_pct_overall,
            "mismatch_rate_pct_virus_only": round(virus_mismatch_pct, 4) if virus_mismatch_pct == virus_mismatch_pct else "NA",
            "virus_genome_len_nt": VIRUS_LEN,
            "virus_mean_depth_dedup": round(cov["mean_depth"], 2) if cov["mean_depth"] == cov["mean_depth"] else "NA",
            "virus_breadth_1x_pct_dedup": round(cov["breadth_1x_pct"], 2) if cov["breadth_1x_pct"] == cov["breadth_1x_pct"] else "NA",
            "virus_breadth_10x_pct_dedup": round(cov["breadth_10x_pct"], 2) if cov["breadth_10x_pct"] == cov["breadth_10x_pct"] else "NA",
        })
        rows.append(row)
        print(f"[OK] {sample}: total={total_reads} virus_frag={virus_frag} host_frag={host_frag} virus%={row['virus_fraction_pct']} uniq%={uniq_pct_overall} mismatch%={mismatch_pct_overall}")

    fieldnames = [
        "sample_id", "status", "total_input_read_pairs", "host_unique_fragments", "virus_unique_fragments",
        "virus_fraction_pct", "uniquely_mapped_pct_overall", "multi_mapped_pct_overall",
        "too_many_loci_pct_overall", "unmapped_too_short_pct_overall",
        "mismatch_rate_pct_overall", "mismatch_rate_pct_virus_only",
        "virus_genome_len_nt", "virus_mean_depth_dedup", "virus_breadth_1x_pct_dedup", "virus_breadth_10x_pct_dedup",
    ]
    with open(out_tsv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    print(f"\nWrote {out_tsv}")


if __name__ == "__main__":
    main()
