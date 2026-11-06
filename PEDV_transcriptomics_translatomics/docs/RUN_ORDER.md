# Run order

Set the paths first (see `config/paths.example.yaml`, `config/env.example.sh`):

```bash
export PEDV_P=/path/to/project_tree_P      # main PEDV project tree
export PEDV_K=/path/to/project_tree_K      # comparative-genomics tree
export PEDV_R=/path/to/reproducibility_package
export PEDV_SERVER_ROOT=/path/to/server_outputs
```

Scripts keep their original relative layout; `config/script_origin.tsv` maps each repository script to its original location so `__file__`-relative paths still resolve.

## Module order

1. **RNA-seq** — `01_alignment` → `02_junction_calling` → `03_classification_abundance` → `04_junction_dynamics` → `05_trs_microhomology_pairing`; `06_external_junction_comparison` is independent (public data).
2. **Ribo-seq** — `01_alignment` → `02_psite_offset_readsets` → `03_periodicity_qc` → `04_price` → `05_hartis` → `06_gRNA_terminal_translation`; `07_pk1_cross_cell` is independent.
3. **Comparative genomics** — `01` (download/filter/dedup) → `02` (alignment, lineage) → `03_orf_scan` → `04_retention_null_stats` → `05_sup14_synplot`.
4. **Figures** — run the script(s) listed per figure after its module has produced the tables.

Steps marked *server* in `SERVER_SCRIPTS_REQUIRED.md` are not runnable from this repository alone; their outputs are the stored inputs used by the local tests.

## Per figure (see `FIGURE_CODE_INDEX.tsv` for panel-level status)

| Figure | Commands (from repo root, after setting paths) |
|---|---|
| Fig 1, Sup 2 | `python rnaseq/03_classification_abundance/fig1_moi_data.py`; `python figures/Fig1/plot_Fig1_MOI.py`; `python figures/Fig1/plot_Fig1E_MOI.py`; `python figures/SupFig2/plot_SuppFig1_QC_MOI.py` |
| Fig 2 | `python rnaseq/03_classification_abundance/fig2_classify.py 55 72 1000 main` (**arguments required**: the script default leader-donor window is 55-71, the submitted tables use 55-72; running without arguments changes 41 of 43,753 event rows); `python figures/Fig2/fig2_plot.py` |
| Fig 3 (data) | `python rnaseq/05_trs_microhomology_pairing/run_s4_new_lowmem.py --workers 1` (single process, ~18 min, same seed 42 / 10,000 permutations / span-stratified swap chain as `s4_new.py`; replaces the original `Pool(26)`). Reproduces the permutation **data** only; the eight-panel final assembly was done in Illustrator. |
| Sup 3, 6, 7 | `s1_core_counts.py`, `s3_panel_stats.py`, `s5_hotspots.py`, then `figures/SupFig3_6_7/*.py` |
| Sup 5 | `python figures/SupFig5/plot_canonical_variants.py` (needs the Fig 2D-G arcs table) |
| Fig 4 | `python figures/Fig4/make_fig4_v2.py` |
| Sup 10-12 | `build_shared_metrics.py`, `prepare_panels.py`, then `figures/SupFig10_12/*.py` |
| Fig 6C (**post-hoc reproduction script**, not the original) | `python figures/Fig6/fig6C_reproduction_from_source_data.py OUTDIR` (needs `R6-E_panel_data_ORF3_included.csv`, produced by `figures/SupFig15/fig_E_onesided.py`) |
| Fig 6D (**post-hoc reproduction script**, not the original; lanes/outlines inferred) | `python figures/Fig6/fig6D_reproduction_from_source_data.py OUTDIR` (needs Sup 16 `SourceData_PanelA_families_gt40.csv`) |
| Fig 6A,B, Sup 15 | `python comparative_genomics/04_retention_null_stats/{g1_g2_g4_summary,g1_conservation,g2_length,g3_g4_tables,g4_activity,stats_abc}.py`; `python figures/Fig6/figs_abc.py`; `python figures/SupFig15/fig_E_onesided.py` |
| Fig 6E, F | `python riboseq/05_hartis/fig6_build_tables.py`; `python figures/Fig6/fig6AB_genome_combined.py`; `python figures/Fig6/fig6C_HAR_PRICE_joint_v3_capped.py` |
| Sup 16 | **order matters**: `04_retention_null_stats/p6_aggregate.py` -> `s1_C12_featured_and_data.py` -> `p7_panels_and_qc.py` -> `r2_audit_tables.py` -> `r3_old_vs_revised.py` -> `r4_finalize.py`; then `python figures/SupFig16/f3_panelC.py`, `f4_make_figure_v3.py`. The 9 featured families are the fixed input `config/sup16_featured_families.json`; scripts abort if `work/featured.json` disagrees. |
| Fig 5 | analysis scripts in `riboseq/06_gRNA_terminal_translation/` (g1 -> g2 -> step* ...; CHX 16 h = EXCLUDED_LIBRARY); the plotting script is **not released** (`EXCLUDED_FROM_PUBLIC_PACKAGE.md`); tables not recomputed |
| Sup 13 | `python figures/SupFig13/landscape_5p_style.py`; `landscape_3p_style.py` (load `riboseq/06_gRNA_terminal_translation/audit3utr.py`; need `SupM_alt_ORF_evidence/cache/psite_by_length.pkl` and `Fig5_3UTR_translation_audit/_work/`, stored server-derived intermediates) |
| Sup 4, 8, 9, 14, 17, 18 | scripts present but **not re-run**; see `UNRESOLVED_ISSUES.md` L6 |

**RPF length sets are not unified**: each figure uses the set stated in `metadata/figure_sample_usage.tsv`. CHX 16 h = RIBO_VERO_CHX_1MOI_16h for Fig 4, Fig 6, Sup 10-13, 15, 17, 18; **Fig 5 input tables use EXCLUDED_LIBRARY** (per the Fig 5 package README and its scripts) — each figure legend/methods must name its library.


All commands were exercised in `tests/isolated_run_results.tsv` (classes A/B only); nothing here reproduces results from raw sequencing data (class C).


## Whitelist rebuild chain (the order used for the final audit; see tests/chain_results.jsonl)
Set `PEDV_P/K` to data-only working copies (no script files). Then, with the registry in place (`metadata/sample.txt`):
1. `riboseq/02_psite_offset_readsets/build_A_alignment_qc.py` -> `build_B_offset_sets.py` -> `build_C_psite_tables.py` -> `build_D_host_tables.py` (RGP tables of the six whitelisted Ribo-seq libraries)
2. `riboseq/03_periodicity_qc/01_frame_counts.py`, `02_frame_ambiguity.py`, `03_tables_audit_figures.py`, `04_strand_check.py`, then `build_shared_metrics.py`, `prepare_panels.py`, `phase1_audit.py`
3. Fig 4: `figures/Fig4/make_fig4_v2.py`; Sup 10-12: `figures/SupFig10_12/*.py`
4. Fig 5: `riboseq/06_gRNA_terminal_translation/step0_junction_remeasure.py`, `step1_build_catalog.py`, `step2_evidence.py`, `step3_spike_frame.py`, `step4_tables.py`, `g1_gRNA5p_analysis.py`, `g2_gRNA5p_tables.py`; Sup 13 / 3' audit: `audit3utr.py`, `audit3utr_b.py`, `audit3utr_c.py`, `figures/SupFig13/landscape_*_style.py`; plot: `figures/Fig5/fig5_independent_reimplementation.py`
5. Fig 6 / Sup 15 / Sup 16 / Fig 1-3 / Sup 2-7 as in the per-figure table above.
Each script stops with a message if a non-whitelisted library appears (`config/sample_registry.py`). `tests/audit_whitelist.py` + `tests/build_final_audit.py` rebuild `tests/FINAL_SAMPLE_WHITELIST_AUDIT*.tsv` from a guarded run log.
