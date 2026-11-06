#!/bin/bash
set -euo pipefail
id=$1
R=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/analysis_GIb/Riboseq_GIb_psite
STAR=${STAR_AVX2_BIN:-STAR-avx2}
SAM=${SAMTOOLS_BIN:-samtools}
PIGZ=$(command -v pigz)
IDX=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/refdata_GIb_v1.1/merged/star_ChlSab1.1_PEDV_GIb_oh40
fq=$(awk -F"\t" -v s=$id '$1==s{print $4}' $R/scripts/samples_GIb.tsv)
[[ -s $fq ]] || { echo "[ERR] missing $fq"; exit 2; }
out=$R/align/$id; mkdir -p $out; pre=$out/$id.
if [[ -s ${pre}Log.final.out ]]; then echo "[$id] skip (done)"; exit 0; fi
rm -rf $R/tmp/$id
echo "[$id] START $(date)"
$STAR --runThreadN 8 --genomeDir $IDX --readFilesIn $fq --readFilesCommand "pigz -dc" \
  --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $pre --limitBAMsortRAM 20000000000 \
  --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --quantMode GeneCounts \
  --outReadsUnmapped Fastx --outTmpDir $R/tmp/$id \
  > $R/logs/$id.STAR.stdout.log 2> $R/logs/$id.STAR.stderr.log
$PIGZ -p 8 -f ${pre}Unmapped.out.mate1
$SAM index -@4 ${pre}Aligned.sortedByCoord.out.bam
$SAM idxstats ${pre}Aligned.sortedByCoord.out.bam > ${pre}idxstats.txt
$SAM flagstat ${pre}Aligned.sortedByCoord.out.bam > ${pre}flagstat.txt
rm -rf $R/tmp/$id
echo "[$id] ALL_DONE $(date)"
