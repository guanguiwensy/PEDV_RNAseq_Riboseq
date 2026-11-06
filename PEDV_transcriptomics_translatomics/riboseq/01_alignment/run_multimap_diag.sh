#!/bin/bash
set -euo pipefail
id=$1
R=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/analysis_GIb/Riboseq_GIb_psite
STAR=${STAR_AVX2_BIN:-STAR-avx2}
SAM=${SAMTOOLS_BIN:-samtools}
IDX=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/refdata_GIb_v1.1/merged/star_ChlSab1.1_PEDV_GIb_oh40
D=$R/multimap_diag/$id; mkdir -p $D; pre=$D/$id.mm100.
[[ -s ${pre}class.tsv ]] && { echo "[$id] skip"; exit 0; }
rm -rf $R/tmp/mm_$id
$STAR --runThreadN 8 --genomeDir $IDX --readFilesIn $R/align/$id/$id.Unmapped.out.mate1.gz --readFilesCommand "pigz -dc" \
  --outSAMtype BAM Unsorted --outFileNamePrefix $pre --outFilterMultimapNmax 100 --outFilterMismatchNmax 2 \
  --outSAMattributes NH --outTmpDir $R/tmp/mm_$id > ${pre}stdout.log 2>&1
$SAM view ${pre}Aligned.out.bam | awk -F'\t' -v V=PEDV_GIb '
  function flush(){ if(q!=""){ k=(hv&&hh)?"host_and_PEDV":(hv?"PEDV_only":"host_only"); c[nh>1?"NH_gt1":"NH_eq1", k]++ } }
  $1!=q{ flush(); q=$1; hv=0; hh=0; nh=0 }
  { if($3==V) hv=1; else hh=1; for(i=12;i<=NF;i++) if($i~/^NH:i:/){ nh=substr($i,6)+0 } }
  END{ flush(); print "nh_class\tref_class\treads"; for(x in c){ split(x,a,SUBSEP); print a[1]"\t"a[2]"\t"c[x] } }' > ${pre}class.tsv
rm -f ${pre}Aligned.out.bam; rm -rf $R/tmp/mm_$id
echo "[$id] MM_DONE $(date)"
