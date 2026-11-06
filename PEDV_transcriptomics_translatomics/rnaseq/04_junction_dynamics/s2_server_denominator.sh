#!/bin/bash


set -euo pipefail
R=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/analysis_GIb/RNAseq
O=$R/fig2supp_denominator; mkdir -p $O/tmp
SM=$(seq 26082190 26082209 | sed 's/^/R/')
printf "sample_id\tviral_frag_primary_all\tviral_frag_primary_q20\tviral_frag_primary_q20_noflagdup\tviral_aln_total_primary\n" > $O/viral_fragments_mm2.tsv
for s in $SM; do
  B=$R/mm2/${s}_mm2.sorted.bam
  a=$(samtools view -F 0x904 $B PEDV_GIb | cut -f1 | sort -u -T $O/tmp | wc -l)
  b=$(samtools view -F 0x904 -q 20 $B PEDV_GIb | cut -f1 | sort -u -T $O/tmp | wc -l)
  c=$(samtools view -F 0x904 -q 20 $B PEDV_GIb | wc -l)
  d=$(samtools view -c -F 0x904 $B PEDV_GIb)
  printf "%s\t%s\t%s\t%s\t%s\n" $s $a $b $c $d >> $O/viral_fragments_mm2.tsv
done
cat $O/viral_fragments_mm2.tsv
