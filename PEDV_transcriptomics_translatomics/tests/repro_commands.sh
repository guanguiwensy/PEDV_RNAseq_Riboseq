#!/usr/bin/env bash
# Reproduction entry points. Run with populated PEDV_P, PEDV_K, and PEDV_R roots.
set -euo pipefail
: "${PEDV_P:?}" "${PEDV_K:?}" "${PEDV_R:?}"
cd "$(dirname "$0")/.."
run() { echo ">>> $*"; python "$@"; }
# RNA-seq
run rnaseq/03_classification_abundance/fig2_classify.py 55 72 1000 main   # T01 (arguments are mandatory)
run figures/Fig2/fig2_plot.py                                      # T02,T03
for s in s1_core_counts s3_panel_stats s5_hotspots; do run rnaseq/04_junction_dynamics/$s.py; done   # T04
for s in s9_panelE_standalone s8_panelC_standalone s4_plot_supp_figure; do run figures/SupFig3_6_7/$s.py; done  # T05
run rnaseq/03_classification_abundance/fig1_moi_data.py            # T17
run figures/Fig1/plot_Fig1_MOI.py; run figures/Fig1/plot_Fig1E_MOI.py; run figures/SupFig2/plot_SuppFig1_QC_MOI.py  # T17,T18
run rnaseq/05_trs_microhomology_pairing/run_s4_new_lowmem.py --workers 1   # T27 (single process, ~15-20 min)
# Ribo-seq
run figures/Fig4/make_fig4_v2.py                                   # T13,T14
run riboseq/03_periodicity_qc/build_shared_metrics.py; run riboseq/03_periodicity_qc/prepare_panels.py  # T15
run figures/SupFig10_12/make_supfig12_v2_ABC.py; run figures/SupFig10_12/make_grouped_frame_bars_v2.py; run figures/SupFig10_12/plot_supplementary_figure_13.py  # T15,T16
run riboseq/05_hartis/fig6_build_tables.py                         # T10
run figures/Fig6/fig6AB_genome_combined.py; run figures/Fig6/fig6C_HAR_PRICE_joint_v3_capped.py  # T11,T12
# Comparative genomics
for s in g1_g2_g4_summary g1_conservation g2_length g3_g4_tables g4_activity stats_abc; do run comparative_genomics/04_retention_null_stats/$s.py; done  # T06
run figures/Fig6/figs_abc.py; run figures/SupFig15/fig_E_onesided.py   # T07,T08
for s in p6_aggregate s1_C12_featured_and_data p7_panels_and_qc; do run comparative_genomics/04_retention_null_stats/$s.py; done
run figures/SupFig16/f3_panelC.py; run figures/SupFig16/f4_make_figure_v3.py   # T21
run figures/SupFig5/plot_canonical_variants.py --input "$PEDV_P/GIb_rnaseq_reanalysis/figures/Fig2D-G_arcs_source_data.csv" --output out_sup5   # T22
run figures/SupFig13/landscape_5p_style.py; run figures/SupFig13/landscape_3p_style.py   # T29
