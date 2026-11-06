#!/bin/bash
R=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/analysis_GIb/Riboseq_GIb_psite
PY=${PYTHON_BIN:-python}
Q=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/analysis_batch23/Riboseq_qc_mapping/bam_out
cd $R; mkdir -p stats stats_oldref
{
tail -n +2 scripts/samples_GIb.tsv | cut -f1 | while read id; do echo "$id $R/align/$id/$id.Aligned.sortedByCoord.out.bam PEDV_GIb $R/stats full"; done
for b in $Q/*/*.Aligned.sortedByCoord.out.bam; do n=$(basename $b .Aligned.sortedByCoord.out.bam); v=$(basename $(dirname $b)); echo "${n#*_} $b $v $R/stats_oldref virus_only"; done
} > logs/stats_jobs.txt
xargs -P 6 -L 1 bash -c "$PY $R/scripts/gib_bam_stats.py \$0 \$1 \$2 \$3 \$4 > $R/logs/stats_\$(basename \$3)_\$0.log 2>&1; echo \"\$0 \$4 rc=\$?\"" < logs/stats_jobs.txt >> logs/stats_all.log
echo STATS_ALL_DONE $(date) >> logs/stats_all.log
