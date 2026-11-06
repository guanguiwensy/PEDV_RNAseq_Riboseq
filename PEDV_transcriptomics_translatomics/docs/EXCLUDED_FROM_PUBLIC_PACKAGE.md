# Files deliberately left out of the public package

| File | Used for | Reason | How to include it later |
|---|---|---|---|
| `figures/Fig5/plot_Fig5_recolored.py` (original: `Fig5_recolor_Fig4palette/plot_Fig5_recolored.py`) | Plotting of Fig 5 A–D and Sup. Fig. S5 panels from the stored gRNA 5′ package tables | The script is a recoloured and fixed version of an externally generated script (the project folder `gpt_plot_script_verified/plot_Fig5_stylematched_portable_FIXED.py`). Authorship, terms of use and licence are not recorded anywhere in the project, and the authors have not stated that it may be published. | Authors confirm authorship/terms and choose a licence, or re-implement the plotting independently; then restore from the project tree. |

Consequences: the original Fig 5 plotting script is not part of the package. In the sample-whitelist round the Fig 5 analysis chain was recomputed with the whitelisted libraries and `figures/Fig5/fig5_independent_reimplementation.py` was written from scratch (it does not reuse the withheld script); it does not reproduce the submitted y-scale or colour calls (see UNRESOLVED_ISSUES L5).

## Withheld or removed in the sample-whitelist round
| file(s) | action | reason |
|---|---|---|
| `bv_common.py`, `s02_four_class.py`, `s03_highfreq.py`, `s14_nc_sgrna_tables_coding.py`, `s16_nc_correlation.py` (rnaseq/06), `s06_figs.py`, `s17_nc_corr_figs.py` (figures/SupFig4_9) | moved to `withheld/` (not in the public ZIP) | they process external datasets that are not approved (IPEC-J2 runs of PRJNA730044); Sup 4/8/9 need a re-implementation restricted to the 4 approved external runs |
| `plot_figures_GIb.py`, `step1_join_check.py` | removed | legacy / diagnostic scripts built on a 23-library RNA set that includes libraries outside the final table |
| the comparison script of the CHX 16 h library against the earlier one | removed | it compares against a library that is not in the final sample table |
