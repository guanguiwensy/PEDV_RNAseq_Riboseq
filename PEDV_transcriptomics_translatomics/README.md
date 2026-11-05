# Integrated Transcriptomic and Translatomic Profiling Reveals Noncanonical RNA Junctions and Hidden Coding Potential in Porcine Epidemic Diarrhea Virus

> **Status:** This directory is being prepared. The analysis code has not been uploaded yet. Scripts, documentation and reproduction instructions will be added step by step and checked against the final manuscript.

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

## Planned organization

The structure below is a plan. The subdirectories do not exist yet, and final directory and file names will follow the code once it is uploaded.

```
PEDV_transcriptomics_translatomics/
├── README.md
├── RNAseq/
├── Riboseq/
├── Comparative_genomics/
└── Documentation/
```

## Data availability

How to obtain the raw sequencing data, reference sequences and required intermediate files will be described here once the data are deposited. Accession numbers are not yet available.

## Software requirements

To be added. This section will list the computing environment, software versions, Python and R packages, and external tools needed for each module.

## Figure reproduction

The manuscript contains Figures 1–6 and Supplementary Figures 1–19. Some panels come from wet-lab experiments, mass spectrometry, GraphPad Prism or schematics and are not covered by this repository. A table linking each computational figure panel to its script, input files and run instructions will be added as the code is released. Figures cannot be reproduced from this directory at present.

## Citation

*Integrated Transcriptomic and Translatomic Profiling Reveals Noncanonical RNA Junctions and Hidden Coding Potential in Porcine Epidemic Diarrhea Virus.* Manuscript submitted.

Full citation details will be added after publication.
