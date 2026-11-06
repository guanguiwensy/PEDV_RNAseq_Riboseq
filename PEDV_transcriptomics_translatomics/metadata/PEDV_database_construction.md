# PEDV genome set used in Fig 6A–D and Sup. Fig. 14–16 (1,083 genomes)

**List**: `PEDV_1083_genomes_accession_list.tsv` — 1,083 rows, one per GenBank accession.version (all public NCBI Nucleotide records; no sequences are redistributed here).
Columns: strain / country / date (from the record), genome length, the two provenance flags below, NCBI status at the access date, QC status, end completeness, number of non-ACGT bases (IUPAC codes are *kept*), S-gene lineage assigned in this study, near-duplicate cluster information, analysis decision (all `included_main_analysis`).

## How the set was built (from the project's sequence-database scripts, `comparative_genomics/01_sequence_database/`)
* Source: NCBI Nucleotide via E-utilities (esearch/efetch); access date **2026-09-24** (`ACCESS_DATE` in `phase2_step2_ncbi_download.py`). E-utilities expose no database version; the date is the version stamp.
* Inclusion: complete genome/sequence records of PEDV with an unambiguous virus name; excluded: partial, synthetic, patent, vector, clone, unknown-origin records; length 25–34 kb; N fraction ≤ 1 %; exact-duplicate sequences removed; recombinant strains kept and only flagged (`recombination_flag`).
* Final analysis set: **1,083 genomes** (27,238–28,381 nt), all with `qc_status = pass`. Genomes whose family status could not be assessed (end-incomplete or ambiguous bases, 294 genomes with ≥ 30 reference nt not covered) are *not assessable*, never counted as loss.
* Lineage: S-gene nearest-panel assignment (SGP v1.0, Guo et al. 2019 anchors): GII-a 752, GII-c 158, GII-b 96, GI-b 30, GI-a 23, unresolved 24.

## Reconciliation of the different counts found in the project
| number | meaning |
|---|---|
| **1,083** | the analysis set (this list; identical accession set in `PEDV_genome_database_QC.csv` and `PEDV_target_manifest.tsv`) |
| 1,041 | records that were already in the earlier local file `PEDV.fasta` (flag `in_original_PEDV_fasta_1041`: 1,041 yes / 42 no) |
| 1,040 + 43 | the manifest's own source label (`sequence_source_in_manifest`: 1,040 `cov_existing`, 43 `ncbi_download`). The two labels are *not* the same partition (they differ by one record); the accession set itself is identical. Not resolved further — both flags are given. |
| 1,084 | `wc -l` of a 1,083-row table including its header — not a genome count |
| 42 | records not returned by the NCBI query at the access date (`not_in_current_ncbi_query`); they are still public GenBank accessions |
| 5 | accessions for which a newer NCBI version existed at the access date (`ncbi_current_version_at_access` = 2) — the listed version is the one analysed |

## Reference genome (not part of the 1,083)
`PEDV_GIb` = refdata_GIb_v1.1, 27,959 nt: a reference-guided consensus of the authors' own CV777-group RNA-seq data (sample RNA_VERO_HAR_0.1MOI_16h; pipeline "single-sample assembly SSA_v1.0"), with 5′ positions 1–7 filled from CV777 and position 8 corrected A→G from read evidence (v1 → v1.1).
**It is colinear with the public record KP162057.1 (same length, no indel) and differs at exactly 8 positions.** The analysis reference has SHA-256 `48c36027…8412b92`.
