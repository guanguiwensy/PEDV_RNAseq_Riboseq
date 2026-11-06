#!/bin/bash
# 主比对完成后依次：BAM 统计（6 并行）与多重比对诊断（3 并行）
R=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/analysis_GIb/Riboseq_GIb_psite
cd $R
while ! grep -q ALL_FINISHED logs/run_all.log; do sleep 30; done
bash scripts/run_stats.sh &
tail -n +2 scripts/samples_GIb.tsv | cut -f1 | xargs -P 3 -I{} bash scripts/run_multimap_diag.sh {} >> logs/multimap_diag.log 2>&1
wait
echo AFTER_STAR_ALL_DONE $(date) >> logs/multimap_diag.log
