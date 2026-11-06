#!/bin/bash
R=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/analysis_GIb/Riboseq_GIb_psite
PY=${PYTHON_BIN:-python}
cd $R; mkdir -p host_psite
for id in $(${PYTHON_BIN:-python3} "${PEDV_REPO:?set PEDV_REPO to the repository root}/config/sample_registry.py" libs ribo_vero); do echo $id; done | \
  xargs -P 7 -I{} bash -c "$PY $R/scripts/host_psite_pass.py {} $R/align/{}/{}.Aligned.sortedByCoord.out.bam $R/scripts/offset_table_GIb.json $R/host_psite > $R/logs/host_psite_{}.log 2>&1; echo {} rc=\$?"
echo HOST_PSITE_DONE $(date)
