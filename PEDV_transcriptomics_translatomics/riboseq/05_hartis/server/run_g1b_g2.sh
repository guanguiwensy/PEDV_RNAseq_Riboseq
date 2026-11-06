#!/bin/bash
# HARTIS_v3 G1b/G2：多读长 + 新旧口径审计（7 库并行）与注释/序列抽取；已完成者跳过
H=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/analysis_GIb/HARTIS_v3_hostpos/G1b_G2
PY=${PYTHON_BIN:-python}
REF=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/refdata_dual_host/vero_ChlSab1.1
GTF=$REF/Chlorocebus_sabaeus.ChlSab1.1.116.gtf
FA=$REF/Chlorocebus_sabaeus.ChlSab1.1.dna.toplevel.fa
A=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/analysis_GIb/Riboseq_GIb_psite/align
export TMPDIR=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/analysis_GIb/HARTIS_v3_hostpos/_tmp
cd $H
REG="${PYTHON_BIN:-python3} ${PEDV_REPO:?set PEDV_REPO to the repository root}/config/sample_registry.py"
declare -A SPEC
SPEC[$($REG lib RIBO_VERO_CHX_1MOI_5h)]=29:12,30:12; SPEC[$($REG lib RIBO_VERO_HAR_1MOI_5h)]=29:12; SPEC[$($REG lib RIBO_VERO_CHX_1MOI_8h)]=29:12,30:12,31:12; SPEC[$($REG lib RIBO_VERO_HAR_1MOI_8h)]=29:12,30:12;
SPEC[$($REG lib RIBO_VERO_HAR_1MOI_16h)]=29:12,30:12; SPEC[$($REG lib RIBO_VERO_CHX_1MOI_16h)]=29:12,30:12,31:12
echo "=== $(date '+%F %T') G1b/G2 开始 ===" >> logs/steps.log; df -h / >> logs/steps.log
for s in "${!SPEC[@]}"; do
  if [ -s $s/meta_multiL_$s.json ]; then echo "$s 已完成，跳过" >> logs/steps.log; continue; fi
  ( echo "$(date '+%F %T') start spec=${SPEC[$s]}" > logs/multiL_$s.log
    $PY scripts/host_cds_psite_multiL_audit.py $s $A/$s/$s.Aligned.sortedByCoord.out.bam $GTF $H ${SPEC[$s]} > logs/multiL_$s.out 2> logs/multiL_$s.err
    echo "$(date '+%F %T') exit=$?" >> logs/multiL_$s.log ) &
done
if [ ! -s annot/tx_seq.fa.gz ]; then
  ( $PY scripts/extract_annot_seq.py $GTF $FA ../$($REG lib RIBO_VERO_CHX_1MOI_8h)/transcripts.tsv annot > logs/annot.out 2> logs/annot.err; echo "annot exit=$?" >> logs/annot.out ) &
fi
wait
echo "=== $(date '+%F %T') G1b/G2 结束 ===" >> logs/steps.log; df -h / >> logs/steps.log
grep -H exit logs/multiL_*.log >> logs/steps.log; tail -1 logs/annot.out >> logs/steps.log
touch logs/G1b_G2.DONE
