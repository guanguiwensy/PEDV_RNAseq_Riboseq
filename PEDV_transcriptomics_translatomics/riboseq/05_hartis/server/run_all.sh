#!/bin/bash
# HARTIS_v3 G1 第 4 步：7 库全量抽取（并行，每库 1 核 1 日志；meta 已存在的库跳过 = 断点续跑）
H=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/analysis_GIb/HARTIS_v3_hostpos
PY=${PYTHON_BIN:-python}
GTF=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/refdata_dual_host/vero_ChlSab1.1/Chlorocebus_sabaeus.ChlSab1.1.116.gtf
A=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/analysis_GIb/Riboseq_GIb_psite/align
export TMPDIR=$H/_tmp
cd $H
echo "=== $(date '+%F %T') STEP4 全量抽取开始 ===" >> logs/steps.log; df -h / >> logs/steps.log
REG="${PYTHON_BIN:-python3} ${PEDV_REPO:?set PEDV_REPO to the repository root}/config/sample_registry.py"
for s in $($REG libs ribo_vero); do
  if [ -s $s/meta_$s.json ]; then echo "$s 已完成，跳过" >> logs/steps.log; continue; fi
  ( echo "$(date '+%F %T') start" > logs/extract_$s.log
    $PY scripts/host_cds_psite29_extract.py $s $A/$s/$s.Aligned.sortedByCoord.out.bam $GTF $H > logs/extract_$s.out 2> logs/extract_$s.err
    echo "$(date '+%F %T') exit=$?" >> logs/extract_$s.log ) &
done
wait
echo "=== $(date '+%F %T') STEP4 全量抽取结束 ===" >> logs/steps.log; df -h / >> logs/steps.log
grep -H exit logs/extract_*.log >> logs/steps.log
touch logs/STEP4.DONE
