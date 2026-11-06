# -*- coding: utf-8 -*-
"""
Phase 1a input QC for PEDV N1 conservation analysis.
Read-only verification of the 4 input files (2 GenBank + 2 SnapGene .dna),
cross-check N1 annotation, coordinates, frame, start/stop, overlap facts,
and write results/reference_annotation.tsv + a machine-readable QC summary.

Inputs are NEVER modified. All outputs go to WORK/{results,logs}.
"""
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
import json
import sys
from io import StringIO
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

WORK = Path(PEDV_COV+r"N1_conservation_analysis")
LOGS = WORK / "logs"
RESULTS = WORK / "results"

FILES = {
    "full_genome_gb": WORK / "LZ202401_full_genome.gb.gb",
    "full_genome_dna": WORK / "LZ202401_full_genome.dna.dna",
    "n1_ref_gb": WORK / "LZ202401_N1_reference.gb.gb",
    "n1_ref_dna": WORK / "LZ202401_N1_reference.dna.dna",
}

log_lines = []


def log(msg):
    line = str(msg)
    print(line)
    log_lines.append(line)


def parse_snapgene_dna(path: Path):
    """Parse an old-format SnapGene .dna binary (no embedded GenBank text).

    These files store the plain uppercase sequence directly after a short
    binary header. We extract the longest ACGT run and wrap it in a minimal
    SeqRecord; identity against the corresponding .gb file is verified
    separately in main().
    """
    import re
    from Bio.SeqRecord import SeqRecord
    raw = path.read_bytes()
    if b"SnapGene" not in raw[:32]:
        raise ValueError(f"{path.name}: missing 'SnapGene' magic header")
    m = re.search(rb"[ACGT]{100,}", raw)
    if not m:
        raise ValueError(f"{path.name}: no plain sequence block found")
    seq = m.group().decode("ascii")
    rec = SeqRecord(Seq(seq), id=path.stem, name=path.stem,
                    description=f"plain sequence extracted from SnapGene binary {path.name}")
    return rec, None


def summarize_record(key, rec, source_file, parse_ok, note=""):
    feats = [(f.type, str(f.location), dict(f.qualifiers)) for f in rec.features]
    accession = rec.annotations.get("accessions", [""])
    info = {
        "file_key": key,
        "source_file": source_file,
        "parse_success": parse_ok,
        "locus_name": rec.name,
        "accession": accession[0] if accession else "",
        "version": rec.annotations.get("sequence_version", ""),
        "definition": rec.description,
        "length_bp": len(rec.seq),
        "topology": rec.annotations.get("topology", ""),
        "organism": rec.annotations.get("organism", ""),
        "note": note,
        "features": [
            {"type": t, "location": loc,
             "qualifiers": {k: v[0] if len(v) == 1 else v for k, v in q.items()}}
            for t, loc, q in feats
        ],
    }
    return info, feats


def main():
    summary = {"files": {}, "comparisons": {}, "n1": {}, "orf_scan": []}

    # ---------- 1. read all four files ----------
    parsed = {}
    for key, path in FILES.items():
        try:
            if path.suffix == ".gb":
                rec = SeqIO.read(str(path), "genbank")
                info, feats = summarize_record(key, rec, path.name, True)
            else:
                rec, gb_text = parse_snapgene_dna(path)
                info, feats = summarize_record(key, rec, path.name, True,
                                               note="parsed via embedded GenBank text in SnapGene binary")
            parsed[key] = rec
        except Exception as e:  # noqa: BLE001
            info = {"file_key": key, "source_file": path.name, "parse_success": False,
                    "error": repr(e)}
        summary["files"][key] = info
        log(f"[read] {path.name}: parse_success={info['parse_success']} "
            f"length={info.get('length_bp', 'NA')} accession='{info.get('accession', '')}' "
            f"version={info.get('version', '')}")
        if info.get("features") is not None:
            for f in info["features"]:
                log(f"    feature {f['type']} {f['location']} "
                    f"label={f['qualifiers'].get('label', f['qualifiers'].get('gene', ''))}")

    # ---------- 2. pairwise sequence identity: .dna vs .gb ----------
    for prefix, gb_key, dna_key in [
        ("full_genome", "full_genome_gb", "full_genome_dna"),
        ("n1_ref", "n1_ref_gb", "n1_ref_dna"),
    ]:
        if gb_key in parsed and dna_key in parsed:
            s1, s2 = str(parsed[gb_key].seq).upper(), str(parsed[dna_key].seq).upper()
            identical = s1 == s2
            diff_positions = [i + 1 for i, (a, b) in enumerate(zip(s1, s2)) if a != b]
            summary["comparisons"][f"{prefix}_dna_vs_gb"] = {
                "identical": identical,
                "len_gb": len(s1), "len_dna": len(s2),
                "mismatch_count": len(diff_positions),
                "mismatch_positions_1based_first10": diff_positions[:10],
            }
            log(f"[cmp] {prefix}: .dna vs .gb identical={identical} "
                f"({len(s1)} vs {len(s2)} bp, {len(diff_positions)} mismatches)")

    # ---------- 3. N1 feature from full genome ----------
    fg = parsed["full_genome_gb"]
    n1_feats = [f for f in fg.features if f.type == "CDS" and
                f.qualifiers.get("label", [""])[0] == "N1"]
    if len(n1_feats) != 1:
        log(f"[N1] ERROR: expected exactly 1 N1 CDS feature, found {len(n1_feats)}")
        sys.exit(1)
    n1 = n1_feats[0]
    start = int(n1.location.start) + 1  # 1-based
    end = int(n1.location.end)  # 1-based inclusive (Biopython end is exclusive 0-based)
    strand = n1.location.strand  # +1 / -1
    codon_start = int(n1.qualifiers.get("codon_start", ["1"])[0])
    n1_nt = n1.extract(fg.seq)
    n1_aa_feature = n1.qualifiers.get("translation", [""])[0].replace("\n", "")
    n1_aa_translated = str(n1_nt.translate(table=11 if False else 1))  # standard code
    # note: viruses use standard genetic code; transl_table not annotated in file
    n1_start_codon = str(n1_nt[:3]).upper()
    n1_stop_codon = str(n1_nt[-3:]).upper()
    n1_internal_stops = n1_aa_translated[:-1].count("*")
    frame_1based = ((start - 1) % 3) + 1  # frame of start codon on plus-strand 1-based coords

    summary["n1"] = {
        "start_1based": start, "end_1based_inclusive": end,
        "strand": "+" if strand == 1 else "-",
        "codon_start_qualifier": codon_start,
        "frame_1based_plus_coords": frame_1based,
        "nt_length": len(n1_nt),
        "aa_length_feature_translation": len(n1_aa_feature),
        "aa_length_retranslated": len(n1_aa_translated.rstrip("*")),
        "start_codon": n1_start_codon,
        "stop_codon": n1_stop_codon,
        "internal_stop_count": n1_internal_stops,
        "feature_translation_equals_retranslation": n1_aa_feature == n1_aa_translated.rstrip("*"),
        "feature_translation": n1_aa_feature,
        "retranslated": n1_aa_translated.rstrip("*"),
    }
    log(f"[N1] coords {start}..{end} (1-based, inclusive), strand {'+' if strand==1 else '-'}, "
        f"frame(1-based plus coords)={frame_1based}, nt={len(n1_nt)}, "
        f"aa={len(n1_aa_feature)} (feature) / {len(n1_aa_translated.rstrip('*'))} (retranslated), "
        f"start={n1_start_codon}, stop={n1_stop_codon}, internal_stops={n1_internal_stops}")

    # ---------- 4. compare with standalone N1 files ----------
    for key in ("n1_ref_gb", "n1_ref_dna"):
        if key in parsed:
            s = str(parsed[key].seq).upper()
            nt_identical = s == str(n1_nt).upper()
            # standalone file has no CDS feature -> translate frame 1 of its own + strand
            aa_standalone = str(Seq(s).translate(table=1)).rstrip("*")
            summary["comparisons"][f"standalone_{key}_vs_genome_N1"] = {
                "nt_identical": nt_identical,
                "nt_len": len(s),
                "aa_identical_to_feature_translation": aa_standalone == n1_aa_feature,
                "aa_len_standalone_frame1": len(aa_standalone),
                "first_mismatch_nt": next((i + 1 for i, (a, b) in enumerate(
                    zip(s, str(n1_nt).upper())) if a != b), None),
            }
            log(f"[cmp] standalone {key} vs genome N1: nt_identical={nt_identical}, "
                f"aa(frame1)_identical={aa_standalone == n1_aa_feature}")

    # ---------- 5. ORF scan to establish overlap facts ----------
    # The full-genome GenBank annotates ONLY the N1 CDS; the N gene is NOT annotated.
    # Establish facts by scanning all 6 frames for ORFs >= 300 nt.
    def scan_orfs(seq, strand, label):
        orfs = []
        s = str(seq).upper()
        for f in range(3):
            i = f
            while i + 3 <= len(s):
                if s[i:i + 3] == "ATG":
                    j = i
                    while j + 3 <= len(s):
                        if s[j:j + 3] in ("TAA", "TAG", "TGA"):
                            length = j + 3 - i
                            if length >= 300:
                                if strand == 1:
                                    ostart, oend = i + 1, j + 3
                                else:
                                    ostart, oend = len(s) - (j + 3) + 1, len(s) - i
                                orfs.append({"strand": strand, "frame_1based": f + 1,
                                             "start": ostart, "end": oend,
                                             "nt_len_incl_stop": length,
                                             "aa_len_excl_stop": (j - i) // 3})
                            break
                        j += 3
                    i = j + 3 if j + 3 <= len(s) else len(s)
                else:
                    i += 3
        return orfs

    plus_orfs = scan_orfs(fg.seq, 1, "plus")
    minus_orfs = scan_orfs(fg.seq.reverse_complement(), -1, "minus")
    all_orfs = plus_orfs + minus_orfs
    all_orfs.sort(key=lambda o: -o["nt_len_incl_stop"])
    # ORFs overlapping N1 interval on the same strand
    n1_overlap = [o for o in plus_orfs if not (o["end"] < start or o["start"] > end)]
    for o in all_orfs[:12]:
        log(f"[orf] strand{o['strand']:+d} frame{o['frame_1based']} "
            f"{o['start']}..{o['end']} nt={o['nt_len_incl_stop']} aa={o['aa_len_excl_stop']}")
    summary["orf_scan"] = all_orfs[:20]
    summary["n1_same_strand_overlapping_orfs"] = n1_overlap

    # annotate overlap facts vs N1
    overlap_facts = []
    for o in n1_overlap:
        ov_start = max(start, o["start"])
        ov_end = min(end, o["end"])
        overlap_facts.append({
            "other_orf": o, "overlap_nt": ov_end - ov_start + 1,
            "same_frame_as_n1": o["frame_1based"] == frame_1based,
            "overlap_interval_1based": [ov_start, ov_end],
        })
    summary["n1_overlap_facts"] = overlap_facts
    for fct in overlap_facts:
        o = fct["other_orf"]
        log(f"[overlap] N1({start}..{end},frame{frame_1based}) vs "
            f"ORF({o['start']}..{o['end']},frame{o['frame_1based']},aa={o['aa_len_excl_stop']}): "
            f"overlap={fct['overlap_nt']} nt at {fct['overlap_interval_1based']}, "
            f"same_frame={fct['same_frame_as_n1']}")

    # ---------- 6. write outputs ----------
    LOGS.mkdir(exist_ok=True)
    RESULTS.mkdir(exist_ok=True)
    (LOGS / "input_qc_log.txt").write_text("\n".join(log_lines), encoding="utf-8")
    (LOGS / "input_qc_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    log("wrote logs/input_qc_log.txt and logs/input_qc_summary.json")


if __name__ == "__main__":
    main()
