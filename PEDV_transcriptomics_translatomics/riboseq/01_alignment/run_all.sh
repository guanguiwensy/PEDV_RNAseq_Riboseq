#!/bin/bash
R=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/analysis_GIb/Riboseq_GIb_psite
cd $R
echo "ALL_START $(date)" >> logs/run_all.log
tail -n +2 scripts/samples_GIb.tsv | cut -f1 | xargs -P 3 -I{} bash scripts/run_star_one.sh {} >> logs/run_all.log 2>&1
echo "ALL_FINISHED $(date)" >> logs/run_all.log
