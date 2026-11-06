#!/bin/bash
# GI-b新参考RNA-seq重比对(23样本,批2+批3合并聚类,不分株)
# STAR参数沿用analysis_sgRNA_junction/scripts/run_qc_align.sh的junction友好参数,索引改为merged host+virus GI-b oh149
set -uo pipefail

BASE=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/analysis_GIb/RNAseq
REF=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/refdata_GIb_v1.1
STAR=${STAR_BIN:-STAR}
SAMTOOLS=${SAMTOOLS_BIN:-samtools}
MM2=${MINIMAP2_BIN:-minimap2}
QC=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/analysis_sgRNA_junction/qc
MANIFEST=$1
MAXPAR=${2:-4}
THREADS=10

mkdir -p $BASE/align $BASE/mm2 $BASE/coverage $BASE/logs

run_one() {
  sample=$1
  logf=$BASE/logs/${sample}.align.log
  echo "[$(date)] START $sample" > $logf
  outdir=$BASE/align/${sample}
  mkdir -p $outdir

  $STAR --runMode alignReads \
    --genomeDir $REF/merged/star_ChlSab1.1_PEDV_GIb_oh149 \
    --readFilesIn $QC/${sample}_1.trim.fq.gz $QC/${sample}_2.trim.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix ${outdir}/ \
    --runThreadN $THREADS \
    --alignIntronMax 30000 \
    --alignMatesGapMax 30000 \
    --scoreGapNoncan 0 --scoreGapGCAG 0 --scoreGapATAC 0 \
    --outFilterMismatchNmax 6 \
    --outFilterMismatchNoverLmax 0.08 \
    --outFilterMultimapNmax 10 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outSJfilterOverhangMin 8 8 8 8 \
    --outSJfilterCountUniqueMin 1 1 1 1 \
    --outSJfilterCountTotalMin 1 1 1 1 \
    --outSJfilterDistToOtherSJmin 0 0 0 0 \
    --outFilterType Normal \
    --outSAMattributes NH HI AS nM NM \
    >> $logf 2>&1
  echo "[$(date)] STAR_DONE $sample" >> $logf

  # 排序+索引(供病毒占比/错配率统计用,STAR原始比对,未去重)
  $SAMTOOLS sort -@ 6 -o $outdir/Aligned.sorted.bam $outdir/Aligned.out.bam >> $logf 2>&1
  $SAMTOOLS index $outdir/Aligned.sorted.bam >> $logf 2>&1

  # 去重流程(坐标法,非UMI建库),去重后BAM供覆盖度使用
  $SAMTOOLS sort -n -@ 6 -o $outdir/nsort.bam $outdir/Aligned.out.bam >> $logf 2>&1
  $SAMTOOLS fixmate -m -@ 6 $outdir/nsort.bam $outdir/fixmate.bam >> $logf 2>&1
  $SAMTOOLS sort -@ 6 -o $outdir/csort.bam $outdir/fixmate.bam >> $logf 2>&1
  $SAMTOOLS markdup -@ 6 -s $outdir/csort.bam $outdir/dedup.bam >> $logf 2>&1
  $SAMTOOLS index $outdir/dedup.bam >> $logf 2>&1
  rm -f $outdir/nsort.bam $outdir/fixmate.bam $outdir/csort.bam $outdir/Aligned.out.bam
  echo "[$(date)] DEDUP_DONE $sample" >> $logf

  # 病毒全基因组逐碱基覆盖度(去重后)
  $SAMTOOLS depth -a -r PEDV_GIb $outdir/dedup.bam > $BASE/coverage/${sample}.depth.tsv 2>> $logf

  # minimap2 splice比对(独立交叉验证,病毒单独索引,未去重原始reads)
  $MM2 -t 8 -ax splice -un -G 30000 -k 13 -w 5 --secondary=no \
    $REF/virus/minimap2/PEDV_GIb.splice.mmi $QC/${sample}_1.trim.fq.gz $QC/${sample}_2.trim.fq.gz 2>> $logf \
    | $SAMTOOLS sort -@ 4 -o $BASE/mm2/${sample}_mm2.sorted.bam -
  $SAMTOOLS index $BASE/mm2/${sample}_mm2.sorted.bam >> $logf 2>&1

  echo "[$(date)] ALL_DONE $sample" >> $logf
}

while read -r sample; do
  [ -z "$sample" ] && continue
  while [ "$(jobs -pr | wc -l)" -ge "$MAXPAR" ]; do sleep 5; done
  echo "[$(date)] LAUNCHING $sample (当前并发数: $(jobs -pr | wc -l))"
  run_one "$sample" &
done < <(tail -n +2 "$MANIFEST")
wait
echo "[$(date)] PIPELINE_ALL_SAMPLES_DONE"
