#!/bin/bash
B=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/analysis_GIb/RNAseq; O=$B/fig1DEF_tmp; mkdir -p $O; rm -f $O/*.counts.tsv $O/DONE
for i in $(seq 82190 82209); do echo R260$i; done | xargs -P 8 -I{} bash $O/fig1DEF_counts.sh {} > $O/run.log 2>&1
echo DONE > $O/DONE
