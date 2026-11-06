#!/bin/bash
# 独立 NH 复核：samtools view -F 0x904（primary mapped），统计全部 / NH>1 / 无 NH 标签 的 read 数（含 PEDV contig）
H=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/analysis_GIb/HARTIS_v3_hostpos
ST=${SAMTOOLS_BIN:-samtools}
A=${PEDV_SERVER_ROOT:?set PEDV_SERVER_ROOT}/analysis_GIb/Riboseq_GIb_psite/align
export TMPDIR=$H/_tmp
cd $H
REG="${PYTHON_BIN:-python3} ${PEDV_REPO:?set PEDV_REPO to the repository root}/config/sample_registry.py"
for s in $($REG libs ribo_vero); do
  ( $ST view -F 0x904 $A/$s/$s.Aligned.sortedByCoord.out.bam | awk -v s=$s 'BEGIN{OFS="\t"} {nh=-1; for(i=12;i<=NF;i++) if(substr($i,1,5)=="NH:i:"){nh=substr($i,6)+0}; t++; if(nh>1)c++; if(nh<0)m++; if($3!~/^PEDV/)h++} END{print s,t,h,c+0,m+0}' > logs/nh_check_$s.tsv ) &
done
wait
( printf "sid\tprimary_mapped_all\tprimary_mapped_host\tNH_gt1\tNH_missing\n"; cat logs/nh_check_X*.tsv ) > logs/nh_check.tsv
touch logs/NHCHECK.DONE
