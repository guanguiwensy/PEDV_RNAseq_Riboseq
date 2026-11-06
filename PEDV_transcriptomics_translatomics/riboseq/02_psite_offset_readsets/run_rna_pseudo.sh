#!/bin/bash
R=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/analysis_GIb/Riboseq_GIb_psite
A=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/analysis_GIb/RNAseq/align
PY=${PYTHON_BIN:-python}
cd $R; mkdir -p rna_pseudo_psite
REG="${PYTHON_BIN:-python3} ${PEDV_REPO:?set PEDV_REPO to the repository root}/config/sample_registry.py"
for id in $(for a in RNA_VERO_CHX_0.1MOI_5h RNA_VERO_CHX_1MOI_5h RNA_VERO_CHX_0.1MOI_12h RNA_VERO_CHX_1MOI_12h RNA_VERO_CHX_0.1MOI_16h RNA_VERO_CHX_1MOI_16h RNA_VERO_HAR_0.1MOI_5h RNA_VERO_HAR_1MOI_5h RNA_VERO_HAR_0.1MOI_8h RNA_VERO_HAR_1MOI_8h RNA_VERO_HAR_0.1MOI_16h RNA_VERO_HAR_1MOI_16h; do $REG lib $a; done); do echo $id; done | \
  xargs -P 12 -I{} bash -c "$PY $R/scripts/rna_pseudo_psite.py {} $A/{}/Aligned.sorted.bam $R/rna_pseudo_psite > $R/logs/rna_pseudo_{}.log 2>&1; echo {} rc=\$?"
echo RNA_PSEUDO_DONE $(date)
