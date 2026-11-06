#!/bin/bash
# 在 R1 上从 minimap2 BAM(与 junction calling 同一比对)统计每样本病毒 fragment 数, 作为 junction reads per million 的分母。
# 只读原始 BAM; 输出写到 data1 下新目录; sort 临时目录也在 data1(不写根分区)。
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
