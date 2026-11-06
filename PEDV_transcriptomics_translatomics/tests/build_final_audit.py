#!/usr/bin/env python3
"""Builds tests/FINAL_SAMPLE_WHITELIST_AUDIT.tsv (figure/panel level) from
   * the registry (aliases, SRR),
   * tests/chain_results.jsonl (what each script read/wrote + exit codes of the guarded run),
   * tests/FINAL_SAMPLE_WHITELIST_AUDIT_FILES.tsv (content scan of every input file, from audit_whitelist.py),
   * the manual figure map below (which stages feed which panel, what was recomputed, differences found, unfinished steps).
Usage: python tests/build_final_audit.py tests/chain_results.jsonl tests/FINAL_SAMPLE_WHITELIST_AUDIT_FILES.tsv tests/FINAL_SAMPLE_WHITELIST_AUDIT.tsv
"""
import os, sys, json, csv, collections
_d = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(_d), "config"))
import sample_registry as SR

G = {
    "RIBO6": list(SR.RIBO_VERO_ALIASES),
    "RNA20": list(SR.RNA_VERO_ALIASES),
    "RNA12_time_pooled": [a for a in SR.RNA_VERO_ALIASES if int(SR.SAMPLES[a]["hpi"]) in (5, 8, 16)],
    "RNA6_partners": [SR.RNA_PARTNER[a] for a in SR.RIBO_VERO_ALIASES],
    "RIBO_CHX3": [SR.ribo_vero("CHX", h) for h in (5, 8, 16)],
    "RIBO_HAR3": [SR.ribo_vero("HAR", h) for h in (5, 8, 16)],
    "PK1_16h": ["RIBO_PK1_CHX_1MOI_16h", "RIBO_PK1_HAR_1MOI_16h", "RNA_PK1_CHX_1MOI_16h", "RNA_PK1_HAR_1MOI_16h"],
    "VERO16_PAIR": ["RIBO_VERO_CHX_1MOI_16h", "RIBO_VERO_HAR_1MOI_16h", "RNA_VERO_CHX_1MOI_16h", "RNA_VERO_HAR_1MOI_16h"],
    "EXT4": ["SRR27751926", "SRR27751927", "SRR27751928", "SRR31828424"],
    "NONE": [],
}
# figure, panels, stages, sample groups, recomputed, old_vs_new, consistent_with_submitted, unfinished
FIG = [
    ("Fig 1", "B-F", "F1_2_S2_7", ["RNA20"], "per-sample tables rebuilt from the stored counts/depths of the 20 libraries (assert_exactly)", "none: all regenerated data files identical to the original", "yes (data); Illustrator layout not compared", "server counting (STAR, samtools depth) not re-run"),
    ("Fig 2", "A-G", "F1_2_S2_7", ["RNA20"], "events recomputed by fig2_classify from junction_events_raw (rows of the 20 libraries selected explicitly + assert); arcs/plots rerun", "none (events_main identical); arcs table of the original has 2 extra columns not written by the archived script", "yes (data)", "minimap2/junction extraction not re-run; stored per-fragment table also holds 3 earlier libraries (filtered, not used)"),
    ("Fig 3", "A-H", "F3_perm", ["RNA20"], "full span-stratified permutation (seed 42, 10,000 draws) recomputed from the filtered events; see tests/chain_results", "see FINAL run log (identical tables expected and checked)", "data yes; eight-panel assembly only in Illustrator", "composite figure script does not exist"),
    ("Fig 4", "A-E", "C0_RGP_tables,C1_periodicity,C2_shared,F4", ["RIBO6", "RNA20"], "whole upstream chain rebuilt from the per-library server statistics of the six whitelisted libraries (RGP tables, offsets, length sets, frame tables, shared metrics) then make_fig4_v2", "all 8 source tables identical to the submitted ones; the 5 libraries shared with the old run are identical, the CHX 16 h rows now come from the whitelisted library; bootstrap CI columns of the frame table change (random stream)", "yes (values identical)", "BAM-level steps (server) not re-run"),
    ("Fig 5", "A,C", "C0_RGP_tables,F5_data,S13,F5_plot", ["RIBO6", "RNA12_time_pooled"], "P-site CPM tracks recomputed from the whitelisted RGP table; RNA depth = mean of the 20 whitelisted libraries; candidate evidence recomputed (step0-4, g1, g2)", "CHX 16 h library replaced: 29 of 176 columns of the 5' evidence table and 29 of 184 of the terminal-ORF table change (all 16 h columns, e.g. HAR/CHX ratio, CHX body reads, occupancy); 124 of 252 terminal candidates change the 16 h CHX body-signal call", "PARTIAL: track shapes/peaks agree (22, 113, 271, ~293); submitted y-scale (sqrt CPM up to ~500) not reproducible from archived data", "original plotting script withheld; independent reimplementation drawn; submitted normalisation unknown"),
    ("Fig 5", "B,D", "F5_data,S13,F5_plot", ["RIBO6", "RNA12_time_pooled"], "candidate catalogues recomputed; red/blue (signal) calls follow the explicit rules in the independent script", "red 5' ORFs: 005,008,009,010,011 (= submitted); the submitted red ORF 023 (UUG@271) is blue under the stated rule; 3' red set differs from the submitted figure", "NOT confirmed: submitted signal rule unrecoverable", "authors must confirm the signal rule"),
    ("Fig 6", "A-D", "F6_AC", ["NONE"], "sequence-only (no experimental libraries)", "none (20 tables identical)", "yes for A,B; C values identical (post-hoc plot); D post-hoc plot", "original C/D plotting programs not found"),
    ("Fig 6", "E,F", "F6_EF", ["RIBO_CHX3", "RIBO_HAR3"], "start-level HARTIS rows of the six libraries selected + asserted (per-library fits are independent); PRICE merged_CHX = three whitelisted CHX libraries; tables and plots rerun", "plotted HAR/PRICE values identical; CHX columns of the family/start summary tables change because the removed library no longer contributes (19 of 74 and 24 of 61 columns; 2 exploratory flags change)", "yes (plot data identical)", "HARTIS G1/G2 and PRICE runs themselves are server outputs and were not re-run"),
    ("Sup 2", "A-E", "F1_2_S2_7", ["RNA20"], "QC plots rerun with the 20-library selection", "none", "yes (data)", "-"),
    ("Sup 3,6,7", "all", "F1_2_S2_7", ["RNA20"], "s1/s3/s5 recomputed from the filtered events; plots rerun", "none (35 files identical)", "yes (data)", "server denominator (viral_fragments_mm2) stored input"),
    ("Sup 5", "A,B", "F1_2_S2_7", ["RNA20"], "replot from the Fig 2D-G arcs table", "none", "yes", "-"),
    ("Sup 4,8,9", "all", "-", ["EXT4", "RNA20"], "NOT re-run", "-", "not assessed", "scripts that also process unapproved datasets (HX, GDS01) are withheld (withheld/); the 4 approved external runs + RNA20 need a re-implementation and server-processed reads"),
    ("Sup 10", "A-E", "C0_RGP_tables,C2_shared,S10_12", ["RIBO6"], "shared metrics rebuilt from the whitelisted tables; panel tables rebuilt", "CHX 16 h rows of the shared metrics change (as expected); the Sup 10 source tables written by the plotting scripts are identical to the submitted ones", "yes (data)", "-"),
    ("Sup 11", "A-F", "C1_periodicity,S10_12", ["RIBO6"], "frame tables rebuilt for the six libraries", "identical after sorting (row order only)", "yes (data)", "-"),
    ("Sup 12", "A-C", "C2_shared,S10_12", ["RIBO6"], "viral periodicity tables rebuilt from per-library stats (all six)", "identical", "yes (data)", "-"),
    ("Sup 13", "A,B", "S13", ["RIBO6", "RNA6_partners"], "3' UTR audit rerun with the registry RNA partners (CHX 8 h / CHX 16 h partners replaced by whitelisted libraries)", "region-level Ribo tables unchanged; RNA-depth based columns change (RNA_depth_fold_of_library_mean in 516 of 704 matched rows, RPF/RNA in 20; RNA-only time course 161 -> 140 rows)", "plots rerun; submitted-figure comparison by data only", "-"),
    ("Sup 14", "A-G", "-", ["NONE"], "NOT re-run (sequence-only; external SynPlot2 workflow missing)", "-", "not assessed", "make_pairs.py and synplot2 executable not found"),
    ("Sup 15", "A,B", "F6_AC", ["RIBO_CHX3", "RIBO_HAR3", "RNA20"], "CHX/HAR totals reproduced exactly from the whitelisted P-sites (CHX 060+811+813 and HAR 061+062+064; ORF1a CHX 181,189; HAR 594,200 ...); RNA depth from the 20 libraries", "plotted numbers identical; the sensitivity column 'RNA_rel_depth_all23_median' is now 'all20'", "yes (numbers)", "Fig6_canonical_background_audit.tsv is a stored table (verified by recomputation, not regenerated)"),
    ("Sup 16", "A,B", "S16", ["NONE"], "sequence-only; featured families fixed", "none (28 tables identical)", "yes", "-"),
    ("Sup 17", "A-D", "-", ["RIBO6"], "NOT re-run", "-", "not assessed", "PRICE locus inputs are server outputs"),
    ("Sup 18", "A,B", "-", ["PK1_16h", "VERO16_PAIR"], "NOT re-run: code restricted to the whitelisted 16 h PK1 libraries (8 h libraries removed); stored PK1 tables not recomputed", "-", "not assessed", "PK1 P-site tables are server-derived; 8 h PK1 libraries still present in the stored upstream tables"),
]


def main(chain, files, out):
    stages = {}
    for l in open(chain, encoding="utf-8"):
        r = json.loads(l); stages[r["test"]] = r
    hits = collections.defaultdict(list)
    for r in csv.DictReader(open(files, encoding="utf-8"), delimiter="\t"):
        if r["content_status"] == "non_whitelist_content":
            hits[r["stage"]].append(os.path.basename(r["input_file(root-relative)"].replace("\\", "/")))
    rows = []
    for fig, panel, st, groups, recomputed, diff, consistent, unfinished in FIG:
        al = [a for g in groups for a in G[g]]
        srr = [SR.srr(a) if a in SR.SAMPLES else a for a in al]
        ids = [s for s in st.split(",") if s and s != "-"]
        run_ok = all(all(x[1] == 0 for x in stages[s]["scripts"]) for s in ids if s in stages) if ids else None
        ran = "not_run" if not ids else ("yes" if run_ok else "FAILED")
        forb = sum(len(d["forbidden"]) for s in ids if s in stages for d in stages[s]["detail"])
        inputs = sorted({"/".join(p.replace("\\", "/").strip("/").split("/")[:2]) for s in ids if s in stages for d in stages[s]["detail"] for p in d["reads"]})
        inputs = [i for i in inputs if not i.lower().endswith(".py")]
        nonwl = sorted({f for s in ids for f in hits.get(s, [])})
        if not ids:
            wl = "n/a (not run)" if al else "n/a"
        elif not al:
            wl = "n/a (no experimental libraries)"
        elif forb:
            wl = "NO (guard fired)"
        elif nonwl:
            wl = "yes after explicit whitelist filter at load (stored multi-library inputs: %s)" % ",".join(nonwl[:6])
        else:
            wl = "yes (no non-whitelisted id in any input)"
        rows.append([fig, panel, ";".join(al), ";".join(srr), "; ".join(inputs[:14]), wl, recomputed, ran, diff, consistent, unfinished])
    with open(out, "w", encoding="utf-8", newline="") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n")
        w.writerow(["figure", "panel", "sample_aliases_used", "SRR", "input_data_roots(root-relative dirs actually read)", "all_inputs_whitelisted", "recomputed",
                    "plot_script_ran", "old_vs_new_statistics", "consistent_with_submitted_figure", "unfinished_steps"])
        w.writerows(rows)
    print("rows", len(rows))


if __name__ == "__main__":
    main(*sys.argv[1:4])
