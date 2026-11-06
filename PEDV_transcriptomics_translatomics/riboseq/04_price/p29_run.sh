#!/bin/bash
# PRICE 29-nt pilot for RIBO_VERO_CHX_1MOI_16h. New tree only: $O/price_29nt. Nothing existing is modified.
source "${PEDV_ENV_ORF:?set PEDV_ENV_ORF (environment file that was sourced on the original server; see docs)}"
set -o pipefail
export PATH=${GEDI_ENV:?set GEDI_ENV (conda env with Gedi/PRICE)}/bin:${GEDI_HOME:?set GEDI_HOME}/bin:$PATH JAVA_HOME=${GEDI_ENV:?set GEDI_ENV (conda env with Gedi/PRICE)} JAVA_TOOL_OPTIONS=-Djava.io.tmpdir=${PEDV_TMPDIR:?set PEDV_TMPDIR (must NOT be on the root partition)}
PY=${PYTHON_BIN:-python}
REG="${PYTHON_BIN:-python3} ${PEDV_REPO:?set PEDV_REPO to the repository root}/config/sample_registry.py"
P=$O/price_29nt; S=$($REG lib RIBO_VERO_CHX_1MOI_16h)
mkdir -p $P/virus/real $P/host/pricebam_29nt $P/host/price $P/logs
L=$P/logs/stage.log; chk; echo "start $(date)" > $L
[ -s $P/host/pricebam_29nt/$S.all.29nt.bam.bai ] || $PY $P/p29_filter.py > $P/logs/filter_29nt.log 2>&1 || { echo "filter failed" >> $L; exit 1; }
echo "filter done $(date)" >> $L
cd $P/host/price
T0=$(date +%s)
# NOTE: attempt 1 with "-filter 29:29" FAILED (see attempt1_filter29_failed); host training therefore uses the 29nt-only BAM with NO -filter
gedi -e Price -keep -reads $P/host/pricebam_29nt/$S.all.29nt.bam -genomic hostfull -prefix $S -fdr 1 -nthreads 8 > $P/logs/price_host_${S}_29nt.log 2>&1
echo "host exit=$? secs=$(( $(date +%s)-T0 )) $(date)" >> $L
run() { strand=$1; pre=${S}_29nt.$strand
  cd $P/virus/real
  for ext in model start.model noise.model; do cp $P/host/price/$S.$ext $P/virus/real/$pre.$ext; done
  gedi -e Price -filter 29:29 -reads $P/virus/$S.PEDV.29nt.bam -genomic PEDVgib_$strand -prefix $pre -fdr 1 -nthreads 3 > $P/logs/price_virus_$pre.log 2>&1
  echo "virus $strand exit=$?" >> $L
  rm -f $P/virus/real/$pre.model $P/virus/real/$pre.start.model $P/virus/real/$pre.noise.model $P/virus/real/$pre.estimateData $P/virus/real/$pre.clusters.csv $P/virus/real/$pre.maxPos
}
run plus & run minus & wait
echo "ALL DONE $(date)" >> $L
