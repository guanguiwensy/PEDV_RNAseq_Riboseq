#!/bin/bash
set -eu
R=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/analysis_GIb/Riboseq_GIb_psite
BT=${BOWTIE2_BIN_DIR:?set BOWTIE2_BIN_DIR}
SAM=${SAMTOOLS_BIN:-samtools}
W=$R/qc_contamination; mkdir -p $W; cd $W
$BT/bowtie2 --version 2>&1 | head -n 2 > bowtie2_version.txt || true
$SAM --version 2>&1 | head -n 1 > samtools_version.txt || true
[[ -s myco_idx.1.bt2 ]] || $BT/bowtie2-build --threads 8 $R/refdata_myco/mycoplasma_panel.fa myco_idx > build.log 2>&1
tail -n +2 $R/scripts/samples_GIb.tsv | while IFS=$'\t' read id th dr fq; do
  [[ -s $id.q20.idxstats.tsv ]] && continue
  $BT/bowtie2 -x myco_idx -U $fq -p 12 --no-unal 2> $id.bt2.log | $SAM sort -@ 4 -o $id.myco.bam -
  $SAM index $id.myco.bam
  $SAM idxstats $id.myco.bam > $id.idxstats.tsv
  $SAM view -b -q 20 $id.myco.bam > $id.q20.bam; $SAM index $id.q20.bam
  $SAM idxstats $id.q20.bam > $id.q20.idxstats.tsv
  echo "$id done $(date)"
done
echo MYCO_ALL_DONE $(date)
