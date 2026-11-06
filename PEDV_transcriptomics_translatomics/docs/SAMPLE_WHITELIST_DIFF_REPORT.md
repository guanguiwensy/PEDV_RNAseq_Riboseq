# Sample whitelist correction: what was wrong, what was recomputed, what changed

Scope: RNA-seq and Ribo-seq results of the 30 libraries of `metadata/sample.txt` only. Comparative-genomics inputs (public virus genomes, reference) are not experimental libraries and are not restricted. The four approved external RNA-seq runs are listed in the README.

## 1. Libraries found outside the final table (and where)

| library (not in the table) | role in the old analyses | where it entered |
|---|---|---|
| the earlier CHX 16 h Ribo-seq library (`XN…065`, now excluded) | CHX 16 h of the old six-library set | Fig 5 package/chain, shared metrics (Sup 10), RGP tables, HARTIS CHX summaries, G2 calibration list, periodicity tables, `s15` screen |
| two further Ribo-seq libraries (`XN…063`, `XN…066`) and one 12 h CHX library (`XN…812`) | processed on the server, never in a figure | host P-site pass loop, module-1 tables (stored), sample lists |
| three earlier RNA-seq libraries (`EXCLUDED_LIBRARY/28/29`) | batch-2 RNA | RNA depth for Fig 5 (time-point pooling of RNA depth), Sup 13 RNA partners (CHX 8 h, CHX 16 h), Sup 15 sensitivity median ("all 23"), RNA manifests, legacy RNA plots |
| PK-1 8 h Ribo-seq and RNA-seq libraries | PK1 pipeline | PK1 P-site/correlation scripts (not used in the submitted Sup 18 panels, which show 16 h only) |
| IPEC-J2 external runs (PRJNA730044, 18 runs: HX, GDS01) | external junction comparison | scripts of the Sup 4/8/9 module (moved to `withheld/`) |

The path guard added in this round found two dependencies that earlier checks had missed: `step2_evidence.py` and `g1_gRNA5p_analysis.py` read the depth file of an earlier RNA library (time-point pooling over all manifest entries).

## 2. What was recomputed (class B: from stored per-library data of the whitelisted libraries; class C was not possible)
* RGP tables (alignment summary, offsets, length sets, P-site tables, host metagenes) rebuilt for the six Vero Ribo-seq libraries -> frame/periodicity tables -> shared metrics -> Fig 4, Sup 10-12 data/plots.
* Fig 5 chain: junction re-measurement, catalogue, evidence (step0-4), 5' evidence (g1/g2), 3' audit (audit3utr, b, c), tracks; independent plot written.
* Sup 15: CHX/HAR totals of the stored background audit table were **reproduced exactly** from the whitelisted P-sites (e.g. ORF1a CHX 181,189; HAR 594,200; N CHX 2,851,226) and the RNA-depth medians recomputed from the 20 libraries.
* Fig 6E/F: HARTIS start-level rows restricted to the six libraries (per-library fits are independent), PRICE merged CHX = three whitelisted libraries; tables and plots rerun.
* Fig 1-3, Sup 2-7: inputs restricted/asserted to the 20 Vero RNA-seq libraries; rerun (plus the full Fig 3 permutation, see the test log).

## 3. Old versus new numbers

| figure | table | rows matched | columns changed | main changes |
|---|---|---|---|---|
| Fig 4 / Sup 11 | `frame_fractions_by_length_v2.csv` | 66/66 | 0 of 13 | identical (row order differs) |
| Fig 4 | 8 source tables + numbers | all | 0 | identical |
| Sup 10 | `sample_funnel.tsv`, offsets | matched libraries | 0 | the five libraries shared with the old set are identical; the CHX 16 h row now belongs to the whitelisted library |
| Sup 10 | `length_summary.tsv` | per library | 6 of 8 | only the CHX 16 h library's rows change |
| Fig 5 (5') | `PEDV_gRNA_5prime_candidate_evidence.tsv` | 24/24 | 29 of 176 | all 16 h CHX columns: e.g. CHX body reads at 16 h of the first candidate 123 -> 531, HAR/CHX ratio 1.50 -> 0.32, occupied codons 5 -> 12; hotspot attribution of 2 candidates shifts by one nt; hotspot windows of 10 candidates gain one position |
| Fig 5 (3') | `PEDV_terminal_ORF_candidates.tsv` | 252/252 | 29 of 184 | all 16 h CHX columns; `CHX_body_signal_16h` changes from False to True for 124 candidates and `n_CHX_supported_timepoints` increases for 124 (more in-frame CHX body support with the whitelisted library); the HAR start signal at 16 h (which is judged against the CHX 16 h library) disappears for 12 candidates |
| Fig 5 | junction re-measurement, 5' catalogue | all | 0 | identical |
| Sup 13 | 3'UTR time course (region level) | 704 matched rows | RNA-depth columns | RNA-only time course 161 -> 140 rows; `RNA_depth_fold_of_library_mean` changes in 516 rows, `RPF_over_RNA_descriptive` in 20; region-level Ribo columns (P-sites) unchanged |
| Sup 13 | candidate CHX body evidence, HAR responses | 159/159 | 0 | identical |
| Fig 6E/F | plotted HAR/PRICE values (`plotdata`) | all | 0 | identical (HAR and PRICE columns are unaffected) |
| Fig 6E/F | HARTIS x PRICE family summary | 74/74 | 19 of 74 | CHX-only summary columns (best start, p, number of valid libraries) change because the removed library no longer contributes; `exploratory_candidate` changes for 2 families |
| Fig 6E/F | HARTIS x PRICE start summary | 578/578 | 24 of 61 | CHX-only columns |
| Sup 15 | G4 activity tables | 7 regions | 0 (numbers) | `RNA_rel_depth_all23_median` -> `..._all20_median` (sensitivity column only); all plotted values identical |
| Sup 16, Fig 6A-D | sequence-only tables | all | 0 | identical |
| periodicity | `read_length_frame_counts_and_fractions.csv` | 63 rows x 6 libraries | bootstrap CI only | `F0_CI95_transcript_bootstrap_*` of three libraries (random-stream order); counts and fractions identical |

## 4. Consequences for the submitted figures
* **Not affected numerically (recomputed and identical):** Fig 1-4, Fig 6, Sup 2-7, 10-12, 15, 16 (plotted values).
* **Affected:** Fig 5 (all CHX 16 h-dependent quantities and the signal calls of B/D) and Sup 13 (RNA-based columns). These need new panels or an explicit statement; the independent Fig 5 script gives one version, but its colour rule and y-scale are not those of the submitted figure.
* **Not re-run (unfinished):** Sup 4/8/9, 14, 17, 18.
