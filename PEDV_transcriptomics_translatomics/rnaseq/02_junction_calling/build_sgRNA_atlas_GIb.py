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
import csv
import glob
import os
import re
import json
from collections import defaultdict, Counter
import pysam

BASE = PEDV_SERVER_ROOT+"/analysis_GIb/RNAseq"
GAPFIX_DIR = BASE + "/junction_atlas/gapfix_extraction"
OUT = BASE + "/junction_atlas"
MANIFEST = BASE + "/sample_manifest_GIb.tsv"
FASTA_PATH = PEDV_SERVER_ROOT+"/refdata_GIb_v1.1/virus/PEDV_GIb.fa"
CONTIG = "PEDV_GIb"
ALIGN_DIR = BASE + "/align"

LEADER_MAX = 100
MAPQ_MIN = 20
ANCHOR_MIN = 15
WINDOWS = [3, 5]
PRIMARY_WINDOW = 5
WIDTH_THRESHOLD = 12
VALLEY_RATIO_THRESHOLD = 0.5
MIN_PEAK_READS = 2
MIN_PEAK_SAMPLES = 1


HIGH_MIN_POSITIVE_SAMPLES = 6
PROBABLE_MULTI_MIN_SAMPLES = 3
PROBABLE_MULTI_MAX_SAMPLES = 5
PROBABLE_MULTI_MIN_READS = 10
PROBABLE_SINGLE_MIN_READS = 20

GENE_COORDS = {
    "ORF1a": (293, 12622), "ORF1b": (12592, 20613), "S": (20610, 24758),
    "ORF3": (24758, 25033), "E": (25361, 25600), "M": (25608, 26288), "N": (26300, 27625),
}

CANONICAL_SITES = {20597: "S", 24718: "ORF3", 25342: "E", 25597: "M", 26285: "N"}

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
    with open(MANIFEST) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            meta[row["sample_id"]] = row
    return meta


def cluster_1d_fixed(positions_sorted, window):
    clusters = []
    cur = [positions_sorted[0]]
    for p in positions_sorted[1:]:
        if p - cur[-1] <= window:
            cur.append(p)
        else:
            clusters.append(cur)
            cur = [p]
    clusters.append(cur)
    return clusters


def main():
    os.makedirs(OUT, exist_ok=True)
    manifest = load_manifest()
    samples = list(manifest.keys())
    print(f"samples: {len(samples)}")

    fasta = pysam.FastaFile(FASTA_PATH)

    def get_flank(pos_1based, upstream=25, downstream=15):
        lo = max(0, pos_1based - 1 - upstream)
        hi = pos_1based - 1 + downstream
        try:
            return fasta.fetch(CONTIG, lo, hi).upper()
        except Exception:
            return ""

    def best_alignment_similarity(seq_a, seq_b, max_shift=10):
        best = (0.0, 0, 0)
        la, lb = len(seq_a), len(seq_b)
        if la == 0 or lb == 0:
            return best
        for shift in range(-max_shift, max_shift + 1):
            matches, total = 0, 0
            for i in range(la):
                j = i + shift
                if 0 <= j < lb:
                    total += 1
                    if seq_a[i] == seq_b[j]:
                        matches += 1
            if total >= 8:
                frac = matches / total
                if frac > best[0]:
                    best = (frac, total, shift)
        return best

    def microhomology_len(donor_end, body_start, max_check=10):
        try:
            up_leader = fasta.fetch(CONTIG, max(0, donor_end - max_check), donor_end).upper()
            up_body = fasta.fetch(CONTIG, max(0, body_start - 1 - max_check), body_start - 1).upper()
        except Exception:
            return 0
        mh = 0
        for i in range(1, min(len(up_leader), len(up_body)) + 1):
            if up_leader[-i] == up_body[-i]:
                mh += 1
            else:
                break
        return mh

    # ---------- STEP 1: read gapfix extraction files, reliable-filter, classify A/B ----------
    print("=== STEP1: pooling reliable reads across 23 samples, classify A(leader-derived)/B(internal) ===")
    pair_groups = defaultdict(lambda: {"read_names": set(), "mapqs": [], "leader_anchors": [],
                                        "body_anchors": [], "leader_starts": [], "cigars": Counter()})
    n_total_rows = 0
    n_reliable = 0
    n_classA = 0
    n_classB = 0
    raw_events_all = []

    for sample in samples:
        fp = os.path.join(GAPFIX_DIR, f"{sample}_mm2_unfiltered.tsv")
        with open(fp) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                n_total_rows += 1
                raw_events_all.append(row)
                is_dup = row["is_duplicate"] == "True"
                is_sec = row["is_secondary"] == "True"
                is_sup = row["is_supplementary"] == "True"
                nh = row["NH"]
                if is_dup or is_sec or is_sup or nh != "1":
                    continue
                n_reliable += 1
                ls = int(row["leader_start"])
                donor = int(row["leader_end"])
                body = int(row["body_start"])
                gap_len = int(row["gap_len"])
                mapq = int(row["mapq"])
                cigar = row["cigar"]
                la, ba = anchors_for_read(cigar, gap_len)
                cls = "A_true_leader_derived" if ls <= LEADER_MAX else "B_internal_genomic_donor"
                if cls == "A_true_leader_derived":
                    n_classA += 1
                else:
                    n_classB += 1
                key = (sample, donor, body, cls)
                g = pair_groups[key]
                g["read_names"].add(row["read_name"])
                g["mapqs"].append(mapq)
                if la is not None:
                    g["leader_anchors"].append(la)
                if ba is not None:
                    g["body_anchors"].append(ba)
                g["leader_starts"].append(ls)
                g["cigars"][cigar] += 1

    print(f"total candidate rows(all samples,unfiltered)={n_total_rows}, reliable={n_reliable}, "
          f"classA(true_leader_derived)={n_classA}, classB(internal_genomic_donor)={n_classB}")


    raw_fields = ["sample", "chrom", "read_name", "leader_start", "leader_end", "body_start", "body_end",
                  "gap_len", "gap_index_in_read", "cigar", "mapq", "NH", "is_duplicate", "is_secondary",
                  "is_supplementary", "is_read1", "is_read2"]
    with open(OUT + "/junction_events_raw.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=raw_fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in raw_events_all:
            w.writerow(row)
    print(f"junction_events_raw.tsv written: {len(raw_events_all)} rows (未过滤,全部候选)")

    # ---------- STEP 2: pool to junction-pair (donor,body,class) level ----------
    print("=== STEP2: pooling to junction-pair level ===")
    jp = defaultdict(lambda: {"reads_by_sample": defaultdict(int), "mapqs": [], "leader_anchors": [],
                               "body_anchors": [], "timepoints": set(), "batches": set()})
    for (sample, donor, body, cls), g in pair_groups.items():
        key = (donor, body, cls)
        d = jp[key]
        d["reads_by_sample"][sample] += len(g["read_names"])
        d["mapqs"].extend(g["mapqs"])
        d["leader_anchors"].extend(g["leader_anchors"])
        d["body_anchors"].extend(g["body_anchors"])
        m = manifest.get(sample, {})
        th = m.get("time_h", "NA")
        if th not in ("NA", None, ""):
            d["timepoints"].add(th)
        batch = m.get("batch", "NA")
        d["batches"].add(batch)

    print(f"unique (donor,body,class) junction-pairs = {len(jp)}")


    donor_counter = Counter()
    for (donor, body, cls), d in jp.items():
        if cls == "A_true_leader_derived":
            donor_counter[donor] += sum(d["reads_by_sample"].values())
    leader_donor_mode = donor_counter.most_common(1)[0][0] if donor_counter else 66
    print(f"leader donor经验众数(全部classA pooled reads) = {leader_donor_mode}")
    leader_seq_ref = get_flank(leader_donor_mode, upstream=25, downstream=5)

    # ---------- STEP 3: position-level (body) aggregates, class A only ----------
    print("=== STEP3: position-level aggregates (class A / leader-derived only) ===")
    posA = defaultdict(lambda: {"reads_by_sample": defaultdict(int), "mapqs": [], "leader_anchors": [],
                                 "body_anchors": [], "donor_reads": Counter(), "timepoints": set(),
                                 "batch2_samples": set(), "batch3_samples": set()})
    for (donor, body, cls), d in jp.items():
        if cls != "A_true_leader_derived":
            continue
        p = posA[body]
        for smp, cnt in d["reads_by_sample"].items():
            p["reads_by_sample"][smp] += cnt
            b = manifest.get(smp, {}).get("batch", "NA")
            if b == "batch2":
                p["batch2_samples"].add(smp)
            elif b == "batch3":
                p["batch3_samples"].add(smp)
        p["mapqs"].extend(d["mapqs"])
        p["leader_anchors"].extend(d["leader_anchors"])
        p["body_anchors"].extend(d["body_anchors"])
        p["donor_reads"][donor] += sum(d["reads_by_sample"].values())
        p["timepoints"] |= d["timepoints"]

    def pos_qc(p):
        total_reads = sum(p["reads_by_sample"].values())
        positive_samples = len(p["reads_by_sample"])
        mq_med = median(p["mapqs"])
        la_med = median(p["leader_anchors"])
        ba_med = median(p["body_anchors"])
        min_anchor = min([x for x in (la_med, ba_med) if x is not None], default=None)
        low_mapq = (mq_med is None or mq_med < MAPQ_MIN)
        short_anchor = (min_anchor is not None and min_anchor < ANCHOR_MIN)
        low_support = (positive_samples == 1 and total_reads <= 2)
        artifact_hard_flag = bool(low_support and (short_anchor or low_mapq))
        technical_QC_pass = (not low_mapq) and (not short_anchor)
        rep_donor = p["donor_reads"].most_common(1)[0][0] if p["donor_reads"] else None
        donor_vals = sorted(p["donor_reads"].keys())
        return dict(total_reads=total_reads, positive_samples=positive_samples,
                    positive_timepoints=len(p["timepoints"]), MAPQ_median=mq_med,
                    leader_anchor_median=la_med, body_anchor_median=ba_med,
                    technical_QC_pass=technical_QC_pass, artifact_hard_flag=artifact_hard_flag,
                    representative_donor=rep_donor,
                    donor_range=f"{donor_vals[0]}-{donor_vals[-1]}" if donor_vals else "NA",
                    batch2_positive_samples=len(p["batch2_samples"]), batch3_positive_samples=len(p["batch3_samples"]))

    qcA = {body: pos_qc(p) for body, p in posA.items()}
    n_pass = sum(1 for v in qcA.values() if v["technical_QC_pass"])
    print(f"class A positions total={len(qcA)}, technical_QC_pass={n_pass}")

    # ---------- STEP 4: canonical assignment (fixed +-5nt window) ----------
    print("=== STEP4: canonical assignment (+-5nt window around 5 known GI-b sites) ===")
    canon_window = 5
    canon_clusters = {}
    assigned = set()
    for site, gene in CANONICAL_SITES.items():
        members = sorted([bp for bp in qcA if abs(bp - site) <= canon_window])
        canon_clusters[gene] = {"main_breakpoint": site, "members": members}
        assigned.update(members)
    n_canon_members = sum(len(c["members"]) for c in canon_clusters.values())
    print(f"canonical members total = {n_canon_members}")
    for gene, c in canon_clusters.items():
        print(f"  {gene}: main={c['main_breakpoint']} members={c['members']}")

    # ---------- STEP 5: STAR independent support (from Phase A SJ.out.tab, leader region only) ----------
    print("=== STEP5: building STAR SJ.out.tab cross-validation set ===")
    star_body_positions = defaultdict(int)
    for sample in samples:
        sj_path = os.path.join(ALIGN_DIR, sample, "SJ.out.tab")
        if not os.path.exists(sj_path):
            continue
        with open(sj_path) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 7 or parts[0] != CONTIG:
                    continue
                intron_start = int(parts[1])  # 1-based, first intron base
                intron_end = int(parts[2])    # 1-based, last intron base
                n_uniq = int(parts[6])
                donor = intron_start - 1
                if donor > LEADER_MAX:
                    continue
                body_start = intron_end + 1
                star_body_positions[body_start] += n_uniq
    print(f"STAR leader-derived SJ body positions (pooled 23 samples) = {len(star_body_positions)}")

    def star_support(body_pos):
        return any(abs(body_pos - sp) <= 2 for sp in star_body_positions)

    # ---------- STEP 6: noncanonical chain-linkage clustering ----------
    print("=== STEP6: noncanonical chain-linkage clustering (windows 3nt, 5nt; FIXED rule) ===")

    def build_noncanonical(window):
        pool = sorted([bp for bp, r in qcA.items() if bp not in assigned and r["technical_QC_pass"]])
        clusters = []
        if pool:
            for cl in cluster_1d_fixed(pool, window):
                clusters.append({"window": window, "members": sorted(cl)})
        return clusters, len(pool)

    noncanon5, pool5 = build_noncanonical(5)
    noncanon3, pool3 = build_noncanonical(3)
    print(f"window=5: candidate pool={pool5}, raw clusters(pre-peak-aware)={len(noncanon5)}")
    print(f"window=3: candidate pool={pool3}, raw clusters={len(noncanon3)}")

    def verify_no_gap_violation(clusters, window):
        v = 0
        for c in clusters:
            m = sorted(c["members"])
            for i in range(1, len(m)):
                if m[i] - m[i - 1] > window:
                    v += 1
        return v

    assert verify_no_gap_violation(noncanon5, 5) == 0
    assert verify_no_gap_violation(noncanon3, 3) == 0

    # ---------- STEP 7: peak-aware splitting (primary window=5 only, width>12nt) ----------
    print("=== STEP7: peak-aware splitting for window=5 clusters with width>12nt ===")

    def peak_aware_split(members):
        reads = [qcA[m]["total_reads"] for m in members]
        n = len(members)
        if n < 3:
            return [members], []
        peak_idx = []
        for i in range(n):
            left_ok = (i == 0) or (reads[i] >= reads[i - 1])
            right_ok = (i == n - 1) or (reads[i] >= reads[i + 1])
            strict = (i == 0 and reads[i] > reads[i + 1]) or (i == n - 1 and reads[i] > reads[i - 1]) or\
                     (0 < i < n - 1 and (reads[i] > reads[i - 1] or reads[i] > reads[i + 1]))
            if left_ok and right_ok and strict:
                peak_idx.append(i)
        if len(peak_idx) < 2:
            return [members], []
        cut_points = []
        diag = []
        for a, b in zip(peak_idx[:-1], peak_idx[1:]):
            seg = reads[a + 1:b]
            valley = min(seg) if seg else min(reads[a], reads[b])
            lower_peak = min(reads[a], reads[b])
            ratio = valley / lower_peak if lower_peak else 1.0
            pa_samples = qcA[members[a]]["positive_samples"]
            pb_samples = qcA[members[b]]["positive_samples"]
            qualifies = (ratio <= VALLEY_RATIO_THRESHOLD and reads[a] >= MIN_PEAK_READS and reads[b] >= MIN_PEAK_READS
                         and pa_samples >= MIN_PEAK_SAMPLES and pb_samples >= MIN_PEAK_SAMPLES)
            diag.append({"peak_a_pos": members[a], "peak_a_reads": reads[a], "peak_b_pos": members[b],
                         "peak_b_reads": reads[b], "valley_reads": valley, "valley_ratio": round(ratio, 3),
                         "qualifies_for_split": qualifies})
            if qualifies:
                valley_local_idx = a + 1 + seg.index(valley) if seg else a
                cut_points.append(valley_local_idx)
        if not cut_points:
            return [members], diag
        result = []
        start = 0
        for cp in sorted(cut_points):
            result.append(members[start:cp + 1])
            start = cp + 1
        result.append(members[start:])
        return [g for g in result if g], diag

    peak_review_log = []
    final_clusters5 = []
    for c in noncanon5:
        width = max(c["members"]) - min(c["members"]) + 1
        if width > WIDTH_THRESHOLD:
            subgroups, diag = peak_aware_split(c["members"])
            verdict = "SPLIT" if len(subgroups) > 1 else "KEPT_SINGLE_CLUSTER"
            peak_review_log.append({"orig_members": ";".join(map(str, c["members"])), "width": width,
                                     "verdict": verdict, "n_subgroups": len(subgroups),
                                     "diagnostics": json.dumps(diag, ensure_ascii=False)})
            for sg in subgroups:
                final_clusters5.append({"members": sg})
        else:
            final_clusters5.append(c)
    n_split = sum(1 for r in peak_review_log if r["verdict"] == "SPLIT")
    print(f"peak-aware review: n_reviewed={len(peak_review_log)}, SPLIT={n_split}, "
          f"final w5 cluster count={len(final_clusters5)} (pre-split raw={len(noncanon5)})")
    with open(OUT + "/peak_aware_review_log.tsv", "w", newline="") as f:
        fields = ["orig_members", "width", "verdict", "n_subgroups", "diagnostics"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in peak_review_log:
            w.writerow(r)

    # ---------- STEP 8: enrich clusters + confidence tiering (recalibrated for n=23) ----------
    print("=== STEP8: enrich clusters, confidence tiering (thresholds recalibrated for n=23 pooled) ===")

    def enrich(members, window, cid):
        member_reads = {m: qcA[m]["total_reads"] for m in members}
        cluster_total = sum(member_reads.values())
        dominant = max(member_reads, key=member_reads.get)
        secondary = sorted([m for m in members if m != dominant], key=lambda m: -member_reads[m])
        samples_all = set()
        b2_all = set()
        b3_all = set()
        for m in members:
            for smp in posA[m]["reads_by_sample"]:
                samples_all.add(smp)
                b = manifest.get(smp, {}).get("batch", "NA")
                if b == "batch2":
                    b2_all.add(smp)
                elif b == "batch3":
                    b3_all.add(smp)
        positive_samples = len(samples_all)
        qc_pass_reads = sum(member_reads[m] for m in members if qcA[m]["technical_QC_pass"])
        pct_qc_pass = qc_pass_reads / cluster_total if cluster_total else 0
        mapqs_all = [qcA[m]["MAPQ_median"] for m in members if qcA[m]["MAPQ_median"] is not None]
        la_all = [qcA[m]["leader_anchor_median"] for m in members if qcA[m]["leader_anchor_median"] is not None]
        ba_all = [qcA[m]["body_anchor_median"] for m in members if qcA[m]["body_anchor_median"] is not None]
        avg_mapq = sum(mapqs_all) / len(mapqs_all) if mapqs_all else None
        min_anchor = min(la_all + ba_all) if (la_all or ba_all) else None
        max_frac = max(member_reads.values()) / cluster_total if cluster_total else 1.0

        gene, dist_to_gene = nearest_gene(dominant)
        body_seq = get_flank(dominant, upstream=25, downstream=15)
        sim_frac, aligned_len, shift = best_alignment_similarity(leader_seq_ref, body_seq)
        trs_class = "TRS-like(>=0.5)" if sim_frac >= 0.5 else ("weak(0.3-0.5)" if sim_frac >= 0.3 else "none(<0.3)")
        rep_donor = qcA[dominant]["representative_donor"]
        mh = microhomology_len(rep_donor, dominant) if rep_donor else 0
        star_supported = star_support(dominant)

        cluster_width = max(members) - min(members) + 1
        stable_cluster = cluster_width <= (window * 2 + 1) * 3
        hard_artifact_combo = (positive_samples == 1 and cluster_total <= 2 and
                                ((min_anchor is not None and min_anchor < ANCHOR_MIN) or
                                 (avg_mapq is None or avg_mapq < MAPQ_MIN)))

        reasons = []
        if hard_artifact_combo:
            confidence = "Low-confidence-uncertain"
            reasons.append("单样本+reads<=2+锚定短或MAPQ低,符合artifact硬特征组合")
        elif positive_samples >= HIGH_MIN_POSITIVE_SAMPLES and pct_qc_pass >= 0.5 and stable_cluster:
            confidence = "High-confidence"
            reasons.append(f"positive_samples={positive_samples}>={HIGH_MIN_POSITIVE_SAMPLES}(n=23重校准阈值); QC占比{pct_qc_pass:.2f}; cluster稳定")
        elif (PROBABLE_MULTI_MIN_SAMPLES <= positive_samples <= PROBABLE_MULTI_MAX_SAMPLES and
              pct_qc_pass >= 0.5 and stable_cluster and cluster_total >= PROBABLE_MULTI_MIN_READS):
            confidence = "Probable"
            reasons.append(f"positive_samples={positive_samples}(n=23重校准区间[{PROBABLE_MULTI_MIN_SAMPLES},{PROBABLE_MULTI_MAX_SAMPLES}]); reads={cluster_total}; QC占比{pct_qc_pass:.2f}")
        elif positive_samples == 1 and cluster_total >= PROBABLE_SINGLE_MIN_READS and sim_frac >= 0.5 and star_supported:
            confidence = "Probable"
            reasons.append(f"单样本但reads={cluster_total}>={PROBABLE_SINGLE_MIN_READS}(n=23重校准); TRS={sim_frac:.2f}; STAR独立支持")
        else:
            confidence = "Low-confidence-uncertain"
            reasons.append(f"未达标准(samples={positive_samples},reads={cluster_total},QC占比{pct_qc_pass:.2f},TRS={sim_frac:.2f},STAR={star_supported})")

        return {
            "cluster_id": cid, "window": window,
            "representative_breakpoint": dominant, "secondary_breakpoints": ";".join(map(str, secondary)),
            "member_breakpoints": ";".join(map(str, sorted(members))),
            "member_breakpoint_reads": ";".join(str(member_reads[m]) for m in sorted(members)),
            "cluster_width": cluster_width, "number_of_member_positions": len(members),
            "body_cluster_start": min(members), "body_cluster_end": max(members),
            "total_reads": cluster_total, "positive_samples": positive_samples,
            "batch2_positive_samples": len(b2_all), "batch3_positive_samples": len(b3_all),
            "max_sample_fraction": round(max_frac, 4),
            "MAPQ_median": round(avg_mapq, 2) if avg_mapq is not None else None,
            "leader_anchor_median": median(la_all), "body_anchor_median": median(ba_all),
            "technical_QC_pass_fraction": round(pct_qc_pass, 4),
            "leader_donor_repr": rep_donor, "genomic_region": gene, "canonical_distance": dist_to_gene,
            "TRS_similarity": round(sim_frac, 3), "TRS_class": trs_class, "microhomology": mh,
            "STAR_independent_support": star_supported,
            "confidence_level": confidence, "QC_notes": "; ".join(reasons),
        }

    enriched5 = [enrich(c["members"], 5, f"ncsgRNA_GIb_{min(c['members'])}_{max(c['members'])}_w5")
                 for c in final_clusters5]
    enriched3 = [enrich(c["members"], 3, f"ncsgRNA_GIb_{min(c['members'])}_{max(c['members'])}_w3")
                 for c in noncanon3]

    conf5 = Counter(c["confidence_level"] for c in enriched5)
    conf3 = Counter(c["confidence_level"] for c in enriched3)
    print(f"w5(主口径) clusters={len(enriched5)} confidence={dict(conf5)}")
    print(f"w3(敏感性) clusters={len(enriched3)} confidence={dict(conf3)}")

    out_fields = ["cluster_id", "window", "representative_breakpoint", "secondary_breakpoints",
                  "member_breakpoints", "member_breakpoint_reads", "cluster_width", "number_of_member_positions",
                  "body_cluster_start", "body_cluster_end", "total_reads", "positive_samples",
                  "batch2_positive_samples", "batch3_positive_samples", "max_sample_fraction", "MAPQ_median",
                  "leader_anchor_median", "body_anchor_median", "technical_QC_pass_fraction",
                  "leader_donor_repr", "genomic_region", "canonical_distance", "TRS_similarity", "TRS_class",
                  "microhomology", "STAR_independent_support", "confidence_level", "QC_notes"]
    with open(OUT + "/noncanonical_clusters_w5_primary.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_fields, delimiter="\t")
        w.writeheader()
        for c in enriched5:
            w.writerow(c)
    with open(OUT + "/noncanonical_clusters_w3_sensitivity.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_fields, delimiter="\t")
        w.writeheader()
        for c in enriched3:
            w.writerow(c)

    # ---------- STEP 9: canonical_5_summary.tsv ----------
    print("=== STEP9: canonical_5_summary.tsv ===")
    canon_rows = []
    for gene, c in canon_clusters.items():
        members = c["members"]
        total_reads = sum(qcA[m]["total_reads"] for m in members)
        samples_all = set()
        for m in members:
            samples_all |= set(posA[m]["reads_by_sample"].keys())
        rep = max(members, key=lambda m: qcA[m]["total_reads"]) if members else c["main_breakpoint"]
        microhet_members = ";".join(f"{m}:{qcA[m]['total_reads']}reads" for m in sorted(members, key=lambda x: -qcA[x]["total_reads"]))
        canon_rows.append({
            "gene": gene, "main_breakpoint_used_for_window": c["main_breakpoint"],
            "representative_body_site_observed": rep, "total_reads": total_reads,
            "positive_samples": len(samples_all), "n_microheterogeneity_members": len(members),
            "microheterogeneity_members_detail": microhet_members,
        })
    with open(OUT + "/canonical_5_summary.tsv", "w", newline="") as f:
        fields = ["gene", "main_breakpoint_used_for_window", "representative_body_site_observed", "total_reads",
                   "positive_samples", "n_microheterogeneity_members", "microheterogeneity_members_detail"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in canon_rows:
            w.writerow(r)
    for r in canon_rows:
        print(f"  {r['gene']}: 主位点={r['representative_body_site_observed']} reads={r['total_reads']} "
              f"samples={r['positive_samples']} microhet_members={r['n_microheterogeneity_members']}")

    # ---------- STEP 10: atlas (canonical + High + Probable) ----------
    print("=== STEP10: sgRNA_atlas_GIb.tsv (canonical + High-confidence + Probable, window=5主口径) ===")
    atlas_fields = ["sgRNA_id", "type", "genomic_region", "representative_breakpoint", "body_cluster_start",
                     "body_cluster_end", "member_breakpoints", "total_reads", "positive_samples",
                     "batch2_positive_samples", "batch3_positive_samples", "confidence_level",
                     "TRS_similarity", "STAR_independent_support", "QC_notes"]
    atlas_rows = []
    for gene, c in canon_clusters.items():
        total_reads = sum(qcA[m]["total_reads"] for m in c["members"])
        samples_all = set()
        b2, b3 = set(), set()
        for m in c["members"]:
            samples_all |= set(posA[m]["reads_by_sample"].keys())
            b2 |= posA[m]["batch2_samples"]
            b3 |= posA[m]["batch3_samples"]
        atlas_rows.append({
            "sgRNA_id": f"canonical_{gene}", "type": "canonical", "genomic_region": gene,
            "representative_breakpoint": c["main_breakpoint"], "body_cluster_start": min(c["members"]) if c["members"] else c["main_breakpoint"],
            "body_cluster_end": max(c["members"]) if c["members"] else c["main_breakpoint"],
            "member_breakpoints": ";".join(map(str, c["members"])), "total_reads": total_reads,
            "positive_samples": len(samples_all), "batch2_positive_samples": len(b2), "batch3_positive_samples": len(b3),
            "confidence_level": "canonical(not_tiered)", "TRS_similarity": "", "STAR_independent_support": "",
            "QC_notes": "canonical主位点(+-5nt窗口,基于refbuild实测观测body起点)",
        })
    for c in enriched5:
        if c["confidence_level"] in ("High-confidence", "Probable"):
            atlas_rows.append({
                "sgRNA_id": c["cluster_id"], "type": "noncanonical", "genomic_region": c["genomic_region"],
                "representative_breakpoint": c["representative_breakpoint"],
                "body_cluster_start": c["body_cluster_start"], "body_cluster_end": c["body_cluster_end"],
                "member_breakpoints": c["member_breakpoints"], "total_reads": c["total_reads"],
                "positive_samples": c["positive_samples"], "batch2_positive_samples": c["batch2_positive_samples"],
                "batch3_positive_samples": c["batch3_positive_samples"], "confidence_level": c["confidence_level"],
                "TRS_similarity": c["TRS_similarity"], "STAR_independent_support": c["STAR_independent_support"],
                "QC_notes": c["QC_notes"],
            })
    with open(OUT + "/sgRNA_atlas_GIb.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=atlas_fields, delimiter="\t")
        w.writeheader()
        for r in atlas_rows:
            w.writerow(r)

    n_hc = sum(1 for c in enriched5 if c["confidence_level"] == "High-confidence")
    n_prob = sum(1 for c in enriched5 if c["confidence_level"] == "Probable")
    n_low = sum(1 for c in enriched5 if c["confidence_level"] == "Low-confidence-uncertain")
    summary = {
        "n_samples": len(samples), "n_canon_members": n_canon_members,
        "n_noncanon_pool_w5": pool5, "n_noncanon_pool_w3": pool3,
        "n_noncanon_clusters_raw_w5": len(noncanon5), "n_noncanon_clusters_w5_after_peakaware": len(enriched5),
        "n_noncanon_clusters_w3": len(enriched3),
        "confidence_w5": dict(conf5), "confidence_w3": dict(conf3),
        "peak_aware_reviewed": len(peak_review_log), "peak_aware_split": n_split,
        "atlas_rows": len(atlas_rows), "atlas_canonical": len(canon_clusters),
        "atlas_high": n_hc, "atlas_probable": n_prob, "n_low_confidence_excluded_from_atlas": n_low,
        "leader_donor_mode": leader_donor_mode,
    }
    with open(OUT + "/pipeline_summary_GIb.json", "w") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2, default=str)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print("ALL_DONE_GIB_ATLAS_BUILD")


if __name__ == "__main__":
    main()
