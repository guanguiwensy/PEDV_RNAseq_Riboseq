# PEDV translatome / transcriptome analysis code (release candidate)

> **Status: first code archive (v0.1) of the PEDV transcriptomics/translatomics analyses. Raw sequencing data are not yet downloadable: BioProject PRJNA1554149 is registered but not yet public, so data-driven re-runs are currently limited.**
> **Reproduction status in one line: class B (recompute from stored analysis source tables) and class A (replot from archived intermediate results) were done and checked; class C (full reproduction from raw sequencing data, FASTQ -> counts) has NOT been done for any figure.**
> Scope: RNA-seq, Ribo-seq and comparative-genomics analyses behind Figures 1–6 and Supplementary Figures 2–18.
> Out of scope (deliberately not included): Western blots, RT-qPCR, PRM/mass-spectrometry, GraphPad Prism files, schematics.

This repository is a *curated subset* (223 scripts incl. 3 added after submission: two post-hoc Fig 6C/6D scripts and the low-memory Fig 3 runner) of a much larger working directory. It keeps the code that produced the
numbers and plots of the final figures; discarded versions, debugging scripts and exploratory directories were left out.
Nothing in the original project was modified, and no result of the submitted paper was recomputed with changed parameters.

## 1. What is where

```
config/                 path configuration (pedv_paths.py), example config, script->original-location map
environment/            software versions actually used, pip/conda drafts
metadata/               samples, references, public data accessions, per-figure sample / RPF-length usage
rnaseq/                 01 alignment -> 02 junction calling -> 03 classification/abundance -> 04 junction dynamics
                        -> 05 TRS similarity / microhomology / pairing -> 06 external-dataset junction comparison
riboseq/                01 alignment -> 02 P-site offset + read-length sets -> 03 periodicity/QC -> 04 PRICE
                        -> 05 HARTIS -> 06 gRNA terminal translation -> 07 PK1 (pig cell) cross-cell comparison
comparative_genomics/   01 sequence database -> 02 alignment/homology/lineage -> 03 ORF scan
                        -> 04 retention + null models + statistics -> 05 Sup. Fig. 14 (SynPlot2 visualisation)
figures/                final plotting scripts, one directory per figure (kept separate from the analysis code)
docs/                   code index per figure, input/output dependencies, run order, test report, open issues
tests/                  commands used for the reproduction tests, first-round ledger and the authoritative isolated re-run results
docs/RELEASE_BLOCKERS.md  items that must be fixed before publication vs. documented known limitations
```

## 2. Three levels of "reproduced" (used throughout the docs)

| Level | Meaning |
|---|---|
| **L1 code provenance** | the script that produced the panel has been identified (string / docstring / output-file evidence), not guessed from appearance |
| **L2 numeric check** | the script was run here on stored inputs and its output tables were compared byte-for-byte (or value-for-value) with the stored/submitted ones |
| **L3 figure reproduction** | the plotting script was run and the regenerated figure was compared with the submitted figure (layout, labels, values) |

A panel is only called *reproduced* when the level is stated. Many submitted figures were assembled or edited in Illustrator;
those edits are **not** in the code and are listed in `docs/UNRESOLVED_ISSUES.md`.
Independently of L1-L3, every test carries a **reproduction class**:

| Class | Meaning | Status |
|---|---|---|
| **A** | replot from archived intermediate results (plot script run on stored source/plot-data tables) | done for the figures listed in the index |
| **B** | recompute from stored analysis source tables (junction event tables, P-site count tables, HARTIS/PRICE start-level tables, ORF-family tables, ...) | **done** (numeric tables byte-identical to the stored ones; see `tests/isolated_run_results.tsv`) |
| **C** | full reproduction from raw sequencing data (trimming, alignment, counting, PRICE, HARTIS host extraction) | **not done**: raw reads are not in the repository, the server was not used |

See `docs/FIGURE_CODE_INDEX.tsv` for the status of every panel and `docs/REPRODUCTION_TEST_REPORT.md` for what was actually run.

## 2b. Sample registry and whitelist (final sample set)

* **Only the 30 experimental libraries of `metadata/sample.txt` may enter RNA-seq / Ribo-seq results** (BioProject PRJNA1554149: 20 Vero RNA-seq, 6 Vero Ribo-seq, 2 PK-1 RNA-seq, 2 PK-1 Ribo-seq). Four external public runs (3 piglet intestine runs of PRJNA1069035 and one LLC-PK run of PRJNA1203018) are approved for the comparative junction analysis; no other external RNA-seq/Ribo-seq data are used.
* Library numbers, SRR accessions and FASTQ names live **only** in `metadata/sample.txt`. Code uses stable aliases such as `RNA_VERO_CHX_0.1MOI_16h`, `RIBO_VERO_CHX_1MOI_16h`, `RIBO_PK1_HAR_1MOI_16h` and resolves them through `config/sample_registry.py` (`SR.lib(alias)`); FASTQ/BAM files are not renamed.
* **Guard (always on, no fallback):** `sample_registry` installs a runtime hook that stops the program if a file whose path contains a library number outside the table is opened, and the loaders assert that the sample set they use is exactly the expected whitelisted set. `tests/audit_whitelist.py` additionally scans the *contents* of every input file that a guarded run opened; stored multi-library tables are filtered explicitly right after loading and re-asserted.
* Libraries that had been processed earlier but are not in the final table are not used anywhere. Results that depended on them were recomputed from the whitelisted libraries' stored per-library data; see `docs/SAMPLE_WHITELIST_DIFF_REPORT.md` (old vs new numbers) and `tests/FINAL_SAMPLE_WHITELIST_AUDIT.tsv` (per figure/panel: aliases, SRR, inputs, recomputation, plot run, differences, unfinished steps).

## 3. Setting up

1. Python 3.11 with `environment/requirements.txt` (matplotlib 3.10, numpy 2.x, pandas 2.3, scipy 1.17, biopython 1.87).
2. Copy `config/paths.example.yaml` to `config/paths.local.yaml` and point `PEDV_P`, `PEDV_K`, `PEDV_R` to **working copies**
   of the data trees (see `docs/IO_DEPENDENCY_MODULES.md`). Many scripts write their outputs next to their inputs, so never point
   these variables at the original data.
3. Server-side scripts additionally need `config/env.example.sh` (copy, edit, `source`). **Never use the root partition for
   caches or temporary files.**
4. Scripts were written for a particular directory layout. Each Python script therefore sets `__file__` to its *original*
   location inside the working tree (see the "path shim" at the top of every script and `config/script_origin.tsv`), so that
   the scripts' own `dirname(__file__)` logic keeps working unchanged. Run them directly: `python rnaseq/03_classification_abundance/fig2_classify.py 55 72 1000 main`.

## 4. Important differences between figures (do not unify)

* **Read lengths** differ by figure: Fig. 4 / Sup. 10–12 use per-library *host-criteria* lengths; Fig. 6E–F and Sup. 17–18 use 29-nt
  P-sites; Sup. 15 uses strict_core 29/30 nt. See `metadata/figure_sample_usage.tsv`.
* **CHX 16 h**: every figure uses `RIBO_VERO_CHX_1MOI_16h` (the CHX 16 h library of the final sample table). Earlier tables that were built with a different CHX 16 h library were rebuilt (Fig. 5, Sup. 13, shared metrics; see the diff report).
* **Fig. 6F** plots the *capped_reference* HARTIS p-value (exploratory, not calibrated for high viral depth) — see `docs/UNRESOLVED_ISSUES.md`.
* RNA-seq design: 2 MOI (0.1, 1) x CHX/HAR x 3/5/8/12/16 hpi, one library per cell, **no biological replicates**; Ribo-seq is MOI 1.

## 5. Post-hoc additions, exclusions and curation patches

* **Fig. 6C / 6D scripts are post-hoc reproduction scripts**, not the original plotting programs (those were not found). Differences from the submitted figures: 6C — none in the plotted values (six printed numbers checked); 6D — lane placement and red outlines inferred, N-x1/N-x2 lane differs.
* **Fig. 3**: the permutation data are reproduced (full run, seed 42, 10,000 draws); the eight-panel final figure was assembled in Illustrator and no script draws it.
* **Fig. 5**: the original plotting script is not in the public package (`docs/EXCLUDED_FROM_PUBLIC_PACKAGE.md`); `figures/Fig5/fig5_independent_reimplementation.py` is an independent reimplementation (explicit signal rules, whitelisted data) that does not reproduce the submitted y-scale or the submitted red/blue calls.
* `figures/Fig6/fig6C_reproduction_from_source_data.py`, `fig6D_reproduction_from_source_data.py`: **written after submission**, not the original plotting
  programs (those were not found). They redraw Fig. 6C/6D from stored source tables, recompute no statistic, and check the plotted values against the submitted figure.
  Lane placement and red outlines of Fig. 6D were inferred from the submitted figure (cosmetic only).
* `rnaseq/05_trs_microhomology_pairing/run_s4_new_lowmem.py`: single-process runner for the full Fig. 3 permutation analysis (seed 42, 10,000 draws, span-stratified background);
  it only replaces the original `Pool(26)`.
* `config/sup16_featured_families.json`: the 9 featured families of Sup. 16 as a fixed input; scripts abort if the working file disagrees.
* All other edits to archived scripts (path literals, Windows separators, a few guards) are listed in `docs/PATH_REWRITE_LOG.tsv` and `docs/CURATION_PATCHES.tsv`.

## 6. Data availability

* **Authors' own sequencing data** (RNA-seq, Ribo-seq of Vero cells; PK1 pig-cell Ribo-seq): NCBI BioProject **PRJNA1554149** — registered, **not yet released to the public**; run-level accessions will be listed after release. (The PK1 libraries are expected to belong to this BioProject; to be confirmed by the authors before release.)
* **Public external data** (Sup. 4/9: PRJNA1069035, PRJNA1203018; others in `metadata/public_data_sources.tsv`) keep their own accessions and are not redistributed.
* **1,083 PEDV genomes** (Fig. 6A–D, Sup. 14–16): accession list in `metadata/PEDV_1083_genomes_accession_list.tsv`, construction rules and count reconciliation in `metadata/PEDV_database_construction.md` (NCBI access date 2026-09-24).
* **Reference `PEDV_GIb` (refdata_GIb_v1.1, 27,959 nt)**: not stored here; it is exactly KP162057.1 plus 8 substitutions — rebuild and verify with `metadata/reconstruct_reference_v1_1.py` (SHA-256 recorded). Host: Vero ChlSab1.1 (Ensembl 116), PK1: Sscrofa11.1.110 (`metadata/reference_genomes.md`).

## 7. Known limitations (state in the README of the public release)

* Figures were assembled / edited in Illustrator; regenerated panels match data and layout, not pixels.
* Server-side steps (trimming, STAR, minimap2 fragment mapping, PRICE, HARTIS host-position extraction, PK1 mapping) are not runnable from this repository;
  their stored outputs are the inputs of the tests (`docs/SERVER_SCRIPTS_REQUIRED.md`).
* Sup. 4, 8, 9, 14, 17, 18 were not re-run; Fig. 5 tables were not recomputed and its plotting script is withheld; Fig. 6C/6D scripts are post-hoc; Fig. 3 has no single 8-panel plotting script.
* The RNA-seq design has no biological replicates.

## 8. Licence

MIT (see `LICENSE`). The copyright line names "the authors of the PEDV transcriptomics/translatomics analyses" generically — replace it by the legal names/institution if required. The licence covers the code in this repository only; the sequencing data (BioProject PRJNA1554149) and external public data have their own terms. One script of unclear provenance was left out (`docs/EXCLUDED_FROM_PUBLIC_PACKAGE.md`).
