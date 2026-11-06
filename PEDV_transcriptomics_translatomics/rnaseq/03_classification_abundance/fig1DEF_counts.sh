#!/bin/bash
# 逐样本: 唯一比对(MAPQ=255)fragment计数(宿主/病毒, 未去重+去重) + 病毒contig逐碱基depth(未去重+去重) + 链取向计数
S=$1
B=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/analysis_GIb/RNAseq
O=$B/fig1DEF_tmp; mkdir -p $O
SAM=${SAMTOOLS_BIN:-samtools}
for v in raw dedup; do
  if [ $v = raw ]; then bam=$B/align/$S/Aligned.sorted.bam; else bam=$B/align/$S/dedup.bam; fi
  # 排除 unmapped(0x4)/secondary(0x100)/supplementary(0x800)；去重版另排除 duplicate(0x400)
  if [ $v = raw ]; then F=0x904; else F=0xD04; fi
  # fragment = read1记录 + (read2且mate未比对)记录
  tot=$(( $($SAM view -c -q 255 -f 0x40 -F $F $bam) + $($SAM view -c -q 255 -f 0x88 -F $F $bam) ))
  vir=$(( $($SAM view -c -q 255 -f 0x40 -F $F $bam PEDV_GIb) + $($SAM view -c -q 255 -f 0x88 -F $F $bam PEDV_GIb) ))
  fw=$($SAM view -c -q 255 -f 0x42 -F $(printf '0x%x' $((F|0x10))) $bam PEDV_GIb)   # 病毒上proper pair的read1正向
  rv=$($SAM view -c -q 255 -f 0x52 -F $F $bam PEDV_GIb)                            # read1反向
  echo -e "$S\t$v\t$tot\t$vir\t$fw\t$rv" >> $O/${S}.counts.tsv
  $SAM depth -a -Q 255 -r PEDV_GIb $bam > $O/${S}.$v.depth.tsv
done
