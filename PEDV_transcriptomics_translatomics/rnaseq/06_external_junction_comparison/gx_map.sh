#!/bin/bash


G=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/public_replication/ESV/gII_cross_strain_validation_v1
P5=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/public_replication/ESV/phase5
MM2=${MINIMAP2_BIN:-minimap2}; ST=${SAMTOOLS_BIN:-samtools}
export TMPDIR=$G/_work/tmp
run=$1; ref=$2
nice -n 10 $MM2 -t 4 -ax splice -un -G 30000 -k 13 -w 5 --secondary=no --sam-hit-only $G/_work/ref/$ref.fa $P5/_work/fq2/$run.vir.fq.gz 2>$G/logs/map_${run}_$ref.log | $ST view -u -F 0x4 - | $ST sort -T $TMPDIR/s_${run}_$ref -o $G/_work/bam/$run.$ref.bam -
$ST index $G/_work/bam/$run.$ref.bam
echo "$(date) map $run $ref done; df: $(df -h / | tail -1)" >> $G/PROGRESS.md
