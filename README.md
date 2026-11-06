# Integrated Transcriptomic and Translatomic Profiling Reveals Noncanonical RNA Junctions and Hidden Coding Potential in Porcine Epidemic Diarrhea Virus

> **Status:** first public code archive (v0.1). The analysis code, configuration, sample registry and test records are uploaded; see [README_CODE_PACKAGE.md](README_CODE_PACKAGE.md) for the full guide. Raw sequencing data are registered under NCBI BioProject PRJNA1554149 but **not yet public**. Reproduction level reached: recomputation from stored analysis source tables and replotting from archived intermediate results; **full reproduction from raw reads (FASTQ) has not been done**.

## Overview

Porcine epidemic diarrhea virus (PEDV) is an alphacoronavirus and a major enteric pathogen of pigs. This study combines time-course RNA-seq, Ribo-seq with cycloheximide (CHX) and harringtonine (HAR) treatment, and comparative genomics to profile the PEDV transcriptome and translatome. The analyses address coronavirus RNA transcription, noncanonical RNA junctions, ribosome occupancy on viral RNA, and the hidden coding potential of the PEDV genome. This directory will hold the code and reproduction notes for the computational analyses reported in the manuscript.

## Repository scope

Only in-house analysis and plotting code for the three modules below is in scope.

### RNA-seq analysis

- Read preprocessing, quality control and alignment
- Quantification of viral and host reads
- PEDV genome coverage
- Identification and classification of canonical and noncanonical RNA junctions
- Abundance, distribution and temporal dynamics of sgRNAs
- TRS-L similarity, microhomology and related junction analyses
- Plotting code for the corresponding manuscript figures

### Ribo-seq analysis

- Read preprocessing and alignment to reference genomes
- RPF length distribution and P-site offset estimation
- Reading-frame periodicity and ribosome occupancy
- CHX- and HAR-based translation signals
- Ribosome occupancy in the 5′ and 3′ UTRs
- Prediction and filtering of noncanonical ORFs
- PRICE and HARTIS (in-house method) analyses and plots

RPF length selection and library choice differ between figures. They will be documented for each analysis as the code is added.

### Comparative genomics

- Collection and processing of coronavirus and PEDV reference genomes
- Multiple sequence alignment and redundancy reduction
- Genome-wide scanning for overlapping ORFs in PEDV
- Distribution of potential ORFs across genomic regions
- Sequence conservation and homology of overlapping ORFs
- Comparative genomics and evolution plots

Third-party tools (for example MAFFT, CD-HIT and SynPlot2) are not redistributed. Only our own wrapper code, configuration and workflow steps will be included.

### Out of scope

Western blot, RT-qPCR, PRM/mass spectrometry, GraphPad Prism analyses, schematic illustrations, and exploratory analyses that are not part of the reported results.

## Organization

```
PEDV_transcriptomics_translatomics/
├── README.md                 this overview
├── README_CODE_PACKAGE.md    detailed guide (levels of reproduction, sample registry, set-up, limitations)
├── LICENSE                   MIT
├── config/                   path configuration and the sample registry (config/sample_registry.py)
├── metadata/                 sample.txt (30 libraries), reference/genome notes, 1,083-genome accession list, figure-to-sample map
├── rnaseq/                   RNA-seq modules (alignment, junctions, classification, dynamics, TRS/microhomology)
├── riboseq/                  Ribo-seq modules (alignment, P-site offsets, periodicity, PRICE, HARTIS, gRNA terminal translation, PK1)
├── comparative_genomics/     sequence database, alignment/lineage, ORF scan, retention and null models
├── figures/                  final plotting scripts, one directory per figure
├── environment/              software versions and environment drafts
├── docs/                     per-figure code index, run order, test report, issue register, release blockers
└── tests/                    reproduction/audit scripts, logs and the final sample-whitelist audit
```

## Data availability

* Authors' sequencing data (20 Vero RNA-seq, 6 Vero Ribo-seq, 2 PK-1 RNA-seq, 2 PK-1 Ribo-seq libraries): NCBI BioProject **PRJNA1554149**, registered, not yet released; library and run accessions are listed in `metadata/sample.txt`.
* Public external data keep their own accessions (`metadata/public_data_sources.tsv`).
* 1,083 PEDV genomes: `metadata/PEDV_1083_genomes_accession_list.tsv` and `metadata/PEDV_database_construction.md`.
* Reference `PEDV_GIb` (27,959 nt) is rebuilt exactly from GenBank KP162057.1 plus 8 substitutions with `metadata/reconstruct_reference_v1_1.py` (SHA-256 verified).

## Software requirements

See `environment/` (`software_versions.tsv`, `requirements.txt`, `environment.yml` draft). Versions marked TO_CONFIRM in `software_versions.tsv` still have to be checked.

## Figure reproduction

The manuscript contains Figures 1–6 and Supplementary Figures 1–19. Some panels come from wet-lab experiments, mass spectrometry, GraphPad Prism or schematics and are not covered here. The per-panel table is `docs/FIGURE_CODE_INDEX.tsv`; run order in `docs/RUN_ORDER.md`; what was reproduced, what differs and what is unfinished: `docs/REPRODUCTION_TEST_REPORT.md`, `docs/SAMPLE_WHITELIST_DIFF_REPORT.md`, `docs/UNRESOLVED_ISSUES.md` and `tests/FINAL_SAMPLE_WHITELIST_AUDIT.tsv`. Not re-run: Sup. Figs 4, 8, 9, 14, 17, 18; Fig. 5 is only partly reproduced (independent plotting script).

## Citation

*Integrated Transcriptomic and Translatomic Profiling Reveals Noncanonical RNA Junctions and Hidden Coding Potential in Porcine Epidemic Diarrhea Virus.* Manuscript submitted.

Full citation details will be added after publication.
