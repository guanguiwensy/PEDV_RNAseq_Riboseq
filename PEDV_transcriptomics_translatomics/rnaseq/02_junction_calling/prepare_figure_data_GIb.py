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
import sample_registry as SR
import csv
import os
import re
import json
from collections import defaultdict, Counter
import pysam

BASE = PEDV_SERVER_ROOT+"/analysis_GIb/RNAseq"
GAPFIX_DIR = BASE + "/junction_atlas/gapfix_extraction"
ATLAS_DIR = BASE + "/junction_atlas"
FIGDATA = BASE + "/figdata"
MANIFEST = BASE + "/sample_manifest_GIb.tsv"
FASTA_PATH = PEDV_SERVER_ROOT+"/refdata_GIb_v1.1/virus/PEDV_GIb.fa"
CONTIG = "PEDV_GIb"

LEADER_MAX = 100
MAPQ_MIN = 20
ANCHOR_MIN = 15
CANONICAL_SITES = {20597: "S", 24718: "ORF3", 25342: "E", 25597: "M", 26285: "N"}
CANON_WINDOW = 5
GENE_COORDS = {
    "ORF1a": (293, 12622), "ORF1b": (12592, 20613), "S": (20610, 24758),
    "ORF3": (24758, 25033), "E": (25361, 25600), "M": (25608, 26288), "N": (26300, 27625),
}

CIGAR_RE = re.compile(r'(\d+)([MIDNSH=X])')


def anchors_for_read(cigar_str, gap_len):
    ops = [(int(l), op) for l, op in CIGAR_RE.findall(cigar_str)]
    n_indices = [i for i, (l, op) in enumerate(ops) if op == 'N' and l == gap_len]
    if not n_indices:
        n_indices = [i for i, (l, op) in enumerate(ops) if op == 'N']
        if not n_indices:
            return None, None
    idx = n_indices[0]
    leader_anchor = sum(l for l, op in ops[:idx] if op in ('M', '=', 'X'))
    body_anchor = 0
    for l, op in ops[idx + 1:]:
        if op in ('M', '=', 'X'):
            body_anchor += l
        elif op in ('N', 'S', 'H'):
            break
    return leader_anchor, body_anchor


def median(lst):
    if not lst:
        return None
    s = sorted(lst)
    n = len(s)
    return s[n // 2] if n % 2 == 1 else (s[n // 2 - 1] + s[n // 2]) / 2


def nearest_gene(pos):
    best_gene, best_dist = None, None
    for gene, (start, end) in GENE_COORDS.items():
        dist = pos - start
        if best_dist is None or abs(dist) < abs(best_dist):
            best_gene, best_dist = gene, dist
    return best_gene, best_dist


def load_manifest():
    meta = {}
    order = []
    with open(MANIFEST) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["sample_id"] not in SR.WHITELIST_LIBRARIES:
                continue          # final sample table only
            meta[row["sample_id"]] = row
            order.append(row["sample_id"])
    SR.assert_exactly(order, SR.RNA_VERO_ALIASES, "prepare_figure_data manifest")
    return meta, order


def load_cluster_map():
    """body_position -> (type, cluster_id_or_gene, confidence)"""
    pos_map = {}
    for site, gene in CANONICAL_SITES.items():
        for p in range(site - CANON_WINDOW, site + CANON_WINDOW + 1):
            pos_map[p] = ("canonical", gene, "canonical")
    with open(ATLAS_DIR + "/noncanonical_clusters_w5_primary.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            members = [int(x) for x in row["member_breakpoints"].split(";")]
            for m in members:
                if m in pos_map:
                    continue  # canonical takes priority (shouldn't overlap given clustering excludes canon-assigned positions)
                pos_map[m] = ("noncanonical", row["cluster_id"], row["confidence_level"])
    return pos_map


def main():
    os.makedirs(FIGDATA, exist_ok=True)
    manifest, sample_order = load_manifest()
    pos_map = load_cluster_map()
    fasta = pysam.FastaFile(FASTA_PATH)

    def get_flank(pos_1based, upstream=25, downstream=15):
        lo = max(0, pos_1based - 1 - upstream)
        hi = pos_1based - 1 + downstream
        try:
            return fasta.fetch(CONTIG, lo, hi).upper()
        except Exception:
            return ""

    print("=== re-pooling reliable classA/classB reads per (sample,donor,body) ===")
    # per-sample x body -> reads (classA only, reliable), also keep donor for 2D map / donor distribution
    sample_body_reads = defaultdict(lambda: defaultdict(int))  # sample -> body -> reads
    event_reads = defaultdict(int)   # (donor,body,cls) -> total reads pooled (all samples)
    event_reads_by_sample = defaultdict(lambda: defaultdict(int))  # (donor,body,cls) -> sample -> reads
    event_qc = {}  # (donor,body,cls) -> dict(mapq_list, la_list, ba_list) pooled, for QC pass determination consistent with stage
    body_donor_reads = defaultdict(lambda: defaultdict(int))  # body -> donor -> reads (classA only)

    for sample in sample_order:
        fp = os.path.join(GAPFIX_DIR, f"{sample}_mm2_unfiltered.tsv")
        with open(fp) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if row["is_duplicate"] == "True" or row["is_secondary"] == "True" or row["is_supplementary"] == "True" or row["NH"] != "1":
                    continue
                ls = int(row["leader_start"])
                donor = int(row["leader_end"])
                body = int(row["body_start"])
                cls = "A" if ls <= LEADER_MAX else "B"
                event_reads_by_sample[(donor, body, cls)][sample] += 1
                event_reads[(donor, body, cls)] += 1
                if cls == "A":
                    sample_body_reads[sample][body] += 1
                    body_donor_reads[body][donor] += 1
        print(f"  pooled {sample}")

    # ---------- 1. per-sample x type read totals (for RNA/J01/SGA03/CAN timecourses) ----------
    print("=== building per-sample type totals ===")
    rows = []
    for sample in sample_order:
        m = manifest[sample]
        totals = Counter()
        gene_reads = Counter()
        for body, reads in sample_body_reads[sample].items():
            info = pos_map.get(body)
            if info is None:
                totals["unclustered_or_QCfail"] += reads
                continue
            typ, ident, conf = info
            if typ == "canonical":
                totals["canonical"] += reads
                gene_reads[ident] += reads
            else:
                totals[conf] += reads
        row = {
            "sample_id": sample, "batch": m["batch"], "treatment": m["treatment"], "time_h": m["time_h"],
            "display_label": m["display_label"],
            "canonical_reads": totals["canonical"],
            "High_confidence_reads": totals["High-confidence"],
            "Probable_reads": totals["Probable"],
            "Low_confidence_reads": totals["Low-confidence-uncertain"],
            "unclustered_or_QCfail_reads": totals["unclustered_or_QCfail"],
            "S_reads": gene_reads["S"], "ORF3_reads": gene_reads["ORF3"], "E_reads": gene_reads["E"],
            "M_reads": gene_reads["M"], "N_reads": gene_reads["N"],
        }
        row["total_classA_reads"] = sum(sample_body_reads[sample].values())
        row["total_leader_body_atlas_reads"] = row["canonical_reads"] + row["High_confidence_reads"] + row["Probable_reads"]
        rows.append(row)
    fields = ["sample_id", "batch", "treatment", "time_h", "display_label", "total_classA_reads",
              "total_leader_body_atlas_reads", "canonical_reads", "High_confidence_reads", "Probable_reads",
              "Low_confidence_reads", "unclustered_or_QCfail_reads", "S_reads", "ORF3_reads", "E_reads", "M_reads", "N_reads"]
    with open(FIGDATA + "/per_sample_junction_type_totals.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"per_sample_junction_type_totals.tsv written: {len(rows)} rows")

    # ---------- 2. event-level table (donor,body,cls) pooled, with classification, for J02/J04/J05/J06 ----------
    print("=== building event-level classified table ===")
    event_rows = []
    for (donor, body, cls), reads in event_reads.items():
        if cls == "A":
            info = pos_map.get(body)
            if info is None:
                etype = "unclustered_or_QCfail"
                gene = nearest_gene(body)[0]
                conf = "NA"
            else:
                typ, ident, conf = info
                etype = "canonical" if typ == "canonical" else conf
                gene = ident if typ == "canonical" else nearest_gene(body)[0]
        else:
            etype = "internal_genomic_donor(classB)"
            gene = nearest_gene(body)[0]
            conf = "NA"
        positive_samples = len(event_reads_by_sample[(donor, body, cls)])
        event_rows.append({
            "donor_position": donor, "body_position": body, "read_class": cls,
            "event_type": etype, "genomic_region": gene, "total_reads": reads,
            "positive_samples": positive_samples,
        })
    with open(FIGDATA + "/event_level_classified.tsv", "w", newline="") as f:
        fields = ["donor_position", "body_position", "read_class", "event_type", "genomic_region",
                   "total_reads", "positive_samples"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in event_rows:
            w.writerow(r)
    print(f"event_level_classified.tsv written: {len(event_rows)} rows")

    # ---------- 3. canonical per-sample microheterogeneity (for CAN_breakpoint_microheterogeneity, CAN_composition_over_time) ----------
    print("=== building canonical per-sample microheterogeneity table ===")
    canon_rows = []
    for sample in sample_order:
        m = manifest[sample]
        for body, reads in sample_body_reads[sample].items():
            info = pos_map.get(body)
            if info is None or info[0] != "canonical":
                continue
            gene = info[1]
            canon_rows.append({"sample_id": sample, "batch": m["batch"], "treatment": m["treatment"],
                                "time_h": m["time_h"], "gene": gene, "body_position": body, "reads": reads})
    with open(FIGDATA + "/canonical_microhet_per_sample.tsv", "w", newline="") as f:
        fields = ["sample_id", "batch", "treatment", "time_h", "gene", "body_position", "reads"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in canon_rows:
            w.writerow(r)
    print(f"canonical_microhet_per_sample.tsv written: {len(canon_rows)} rows")

    # ---------- 4. per-sample x cluster matrix (atlas entries only) for SGA02 heatmap ----------
    print("=== building sample x atlas-cluster matrix ===")
    atlas_rows = list(csv.DictReader(open(ATLAS_DIR + "/sgRNA_atlas_GIb.tsv"), delimiter="\t"))
    cluster_body_map = {}  # sgRNA_id -> list of member body positions
    for r in atlas_rows:
        members = [int(x) for x in r["member_breakpoints"].split(";")] if r["member_breakpoints"] else\
                  [int(r["representative_breakpoint"])]
        cluster_body_map[r["sgRNA_id"]] = members
    matrix_rows = []
    for sample in sample_order:
        m = manifest[sample]
        row = {"sample_id": sample, "batch": m["batch"], "treatment": m["treatment"], "time_h": m["time_h"]}
        for cid, members in cluster_body_map.items():
            row[cid] = sum(sample_body_reads[sample].get(bp, 0) for bp in members)
        matrix_rows.append(row)
    fieldnames = ["sample_id", "batch", "treatment", "time_h"] + list(cluster_body_map.keys())
    with open(FIGDATA + "/sample_x_atlas_cluster_matrix.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in matrix_rows:
            w.writerow(r)
    print(f"sample_x_atlas_cluster_matrix.tsv written: {len(matrix_rows)} samples x {len(cluster_body_map)} clusters")

    # ---------- 5. sequence flanks for canonical (logo + TRS alignment) ----------
    print("=== extracting canonical junction flanking sequences ===")
    with open(FIGDATA + "/canonical_body_flanks.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["gene", "body_position", "flank_upstream25_downstream15"])
        for site, gene in CANONICAL_SITES.items():
            w.writerow([gene, site, get_flank(site, 25, 15)])
    # leader reference flank (mode donor from atlas build step)
    with open(ATLAS_DIR + "/pipeline_summary_GIb.json") as f:
        summ = json.load(f)
    leader_donor_mode = summ["leader_donor_mode"]
    with open(FIGDATA + "/leader_reference_flank.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["leader_donor_mode", "flank_upstream25_downstream5"])
        w.writerow([leader_donor_mode, get_flank(leader_donor_mode, 25, 5)])

    print("ALL_DONE_FIGDATA_PREP")


if __name__ == "__main__":
    main()
