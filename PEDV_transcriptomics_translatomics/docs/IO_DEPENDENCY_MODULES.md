# Input / output dependencies by module

A per-script heuristic table (`IO_DEPENDENCY_TABLE.tsv`, auto-extracted from path literals and `open/read_csv/to_csv/savefig` calls) accompanies this file. This file gives the curated module-level view. "Server" = produced on the HPC and stored in the project tree; "Local" = produced by scripts in this repository.

| Module | Main inputs | Main outputs | Feeds figure(s) | Upstream dependency |
|---|---|---|---|---|
| `rnaseq/01_alignment` | trimmed RNA-seq FASTQ (23 libraries), merged STAR index oh149 | BAM, minimap2 PAF/BAM | Fig 1-3 | Server (STAR 2.7.11b, minimap2 2.30) |
| `rnaseq/02_junction_calling` | minimap2 alignments, `PEDV_GIb` FASTA | per-fragment junction table | Fig 2, 3, Sup 3-9 | 01 |
| `rnaseq/03_classification_abundance` | junction table, annotation, depth counts | `events_main.tsv`, per-sample class counts, Fig 1 source data | Fig 1, 2, Sup 2 | 02; Fig 1 counts from server (`fig1DEF_counts.sh`) |
| `rnaseq/04_junction_dynamics` | events table, server denominator | s1-s5 tables | Sup 3, 6, 7 | 03 |
| `rnaseq/05_trs_microhomology_pairing` | events table, FASTA (v1), ViennaRNA | microhomology / MFE statistics, permutation nulls | Fig 3 | 03; ViennaRNA 2.7.2 |
| `rnaseq/06_external_junction_comparison` | public RNA-seq (PRJNA1069035, PRJNA1203018) mapped to strain assemblies | cross-dataset junction tables | Sup 4, 8, 9 | Server (public reads) |
| `riboseq/01_alignment` | Ribo-seq FASTQ, merged index oh40 | BAMs | all Ribo figs | Server (cutadapt, STAR, `--outFilterMultimapNmax 1`) |
| `riboseq/02_psite_offset_readsets` | BAM-derived read-length x frame tables | P-site offsets (method B), read sets strict_core / fixed_29_30 / expanded / short_RPF / host_criteria | Fig 4-6, Sup 10-13 | 01 |
| `riboseq/03_periodicity_qc` | stored offset tables, one-pass stats | shared metrics, frame fractions, Sup 10-12 source data | Fig 4, Sup 10-12 | 02 |
| `riboseq/04_price` | per-library BAMs, host annotation | PRICE ORF calls (`orfs.tsv`) | Fig 6E/F, Sup 17 | Server (PRICE 1.0.5) |
| `riboseq/05_hartis` | host per-position data (G1), PRICE outputs, P-site tables | start-level HARTIS table, HARTIS x PRICE summary | Fig 6E/F | 02, 04; G1/G2 on server |
| `riboseq/06_gRNA_terminal_translation` | P-site tables, stored audit tables | 5'/3' terminal windows | Fig 5, Sup 13 | 02 |
| `riboseq/07_pk1_cross_cell` | PK1 and Vero P-site profiles | correlation tables | Sup 18 | Server (PK1 mapping) |
| `comparative_genomics/01_*` | NCBI PEDV genomes | filtered, de-duplicated FASTA | Fig 6, Sup 14-16 | Network access (not run here) |
| `comparative_genomics/02_*` | FASTA | MAFFT alignment, lineage tree/labels | Sup 14, 16 | MAFFT 7.525 |
| `comparative_genomics/03_orf_scan` | alignments, reference | stop-anchored ORF families, start-codon classes, homology mapping | Fig 6, Sup 16 | 02 |
| `comparative_genomics/04_retention_null_stats` | ORF tables, Ribo audit table | G1-G4 tables, composition and synonymous nulls (seed 20261019), Sup 15/16 tables | Fig 6A-C, Sup 15, 16 | 03; Ribo audit table from 02/03 of riboseq |
| `comparative_genomics/05_sup14_synplot` | per-sequence alignments | three-frame heat maps | Sup 14 | FRASynPlot (external) |
| `figures/*` | tables above | PDF/SVG/PNG panels + source-data CSVs | see `FIGURE_CODE_INDEX.tsv` | corresponding analysis module |
