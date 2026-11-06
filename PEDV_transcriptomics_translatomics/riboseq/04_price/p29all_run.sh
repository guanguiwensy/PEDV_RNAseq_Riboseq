#!/bin/bash
# PRICE 29-nt, 4 scopes (060, 811, 813, merged_CHX). New tree: $O/price_29nt/all4. Recipe = pilot p29_run.sh (host training WITHOUT -filter; virus -filter 29:29; host-trained models reused for plus/minus).
source "${PEDV_ENV_ORF:?set PEDV_ENV_ORF (environment file that was sourced on the original server; see docs)}"
set -o pipefail
export PATH=${GEDI_ENV:?set GEDI_ENV (conda env with Gedi/PRICE)}/bin:${GEDI_HOME:?set GEDI_HOME}/bin:$PATH JAVA_HOME=${GEDI_ENV:?set GEDI_ENV (conda env with Gedi/PRICE)} JAVA_TOOL_OPTIONS=-Djava.io.tmpdir=${PEDV_TMPDIR:?set PEDV_TMPDIR (must NOT be on the root partition)}
SAM=${SAMTOOLS_BIN:-samtools}; PY=${PYTHON_BIN:-python}
P=$O/price_29nt/all4; PP=$O/price_29nt
mkdir -p $P/virus/real $P/host/pricebam_29nt $P/host/price $P/logs
L=$P/logs/stage.log; echo "start $(date)" > $L
# 813: reuse pilot BAMs/models (identical recipe); copy (not move) into the new tree
REG="${PYTHON_BIN:-python3} ${PEDV_REPO:?set PEDV_REPO to the repository root}/config/sample_registry.py"
S_CHX5=$($REG lib RIBO_VERO_CHX_1MOI_5h); S_CHX8=$($REG lib RIBO_VERO_CHX_1MOI_8h); S_CHX16=$($REG lib RIBO_VERO_CHX_1MOI_16h)
cp -n $PP/virus/${S_CHX16}.PEDV.29nt.bam $PP/virus/${S_CHX16}.PEDV.29nt.bam.bai $P/virus/
cp -n $PP/host/pricebam_29nt/${S_CHX16}.all.29nt.bam $PP/host/pricebam_29nt/${S_CHX16}.all.29nt.bam.bai $P/host/pricebam_29nt/
cp -n $PP/logs/filter_counts_29nt_${S_CHX16}.json $P/logs/
for ext in model start.model noise.model orfs.tsv; do cp -n $PP/host/price/${S_CHX16}.$ext $P/host/price/; done
for f in $PP/virus/real/${S_CHX16}_29nt.*; do cp -n $f $P/virus/real/; done
cp -n $PP/logs/price_host_${S_CHX16}_29nt.log $PP/logs/price_virus_${S_CHX16}_29nt.*.log $P/logs/
echo "813 reused from pilot $(date)" >> $L
( $PY $PP/p29all_filter.py ${S_CHX5} > $P/logs/filter_29nt_${S_CHX5}.log 2>&1 ) &
( $PY $PP/p29all_filter.py ${S_CHX8} > $P/logs/filter_29nt_${S_CHX8}.log 2>&1 ) &
wait
echo "filters done $(date)" >> $L
for s in ${S_CHX5} ${S_CHX8}; do [ -s $P/host/pricebam_29nt/$s.all.29nt.bam.bai ] || { echo "filter failed $s" >> $L; exit 1; }; done
$SAM merge -f -@4 $P/virus/CHX.all.PEDV.29nt.bam $P/virus/${S_CHX5}.PEDV.29nt.bam $P/virus/${S_CHX8}.PEDV.29nt.bam $P/virus/${S_CHX16}.PEDV.29nt.bam && $SAM index $P/virus/CHX.all.PEDV.29nt.bam
$SAM merge -f -@8 $P/host/pricebam_29nt/merged_CHX.all.29nt.bam $P/host/pricebam_29nt/${S_CHX5}.all.29nt.bam $P/host/pricebam_29nt/${S_CHX8}.all.29nt.bam $P/host/pricebam_29nt/${S_CHX16}.all.29nt.bam && $SAM index -@4 $P/host/pricebam_29nt/merged_CHX.all.29nt.bam
echo "merge done $(date)" >> $L
cd $P/host/price
hostrun() { s=$1; bam=$2; T0=$(date +%s)
  gedi -e Price -keep -reads $bam -genomic hostfull -prefix $s -fdr 1 -nthreads 8 > $P/logs/price_host_${s}_29nt.log 2>&1
  echo "host $s exit=$? secs=$(( $(date +%s)-T0 )) $(date)" >> $L; }
hostrun ${S_CHX5} $P/host/pricebam_29nt/${S_CHX5}.all.29nt.bam &
hostrun ${S_CHX8} $P/host/pricebam_29nt/${S_CHX8}.all.29nt.bam &
wait
hostrun merged_CHX $P/host/pricebam_29nt/merged_CHX.all.29nt.bam
run() { sc=$1; bam=$2; strand=$3; pre=${sc}_29nt.$strand; T0=$(date +%s)
  cd $P/virus/real
  for ext in model start.model noise.model; do cp $P/host/price/$sc.$ext $P/virus/real/$pre.$ext; done
  gedi -e Price -filter 29:29 -reads $bam -genomic PEDVgib_$strand -prefix $pre -fdr 1 -nthreads 3 > $P/logs/price_virus_$pre.log 2>&1
  echo "virus $sc $strand exit=$? secs=$(( $(date +%s)-T0 ))" >> $L
  rm -f $P/virus/real/$pre.model $P/virus/real/$pre.start.model $P/virus/real/$pre.noise.model $P/virus/real/$pre.estimateData $P/virus/real/$pre.clusters.csv $P/virus/real/$pre.maxPos
}
for sc in ${S_CHX5} ${S_CHX8}; do run $sc $P/virus/$sc.PEDV.29nt.bam plus & run $sc $P/virus/$sc.PEDV.29nt.bam minus & done
run merged_CHX $P/virus/CHX.all.PEDV.29nt.bam plus & run merged_CHX $P/virus/CHX.all.PEDV.29nt.bam minus & wait
echo "ALL DONE $(date)" >> $L
