#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""§3–§4：比对汇总、新旧参考对比、读长分布（宿主/病毒分开）、LPF、rRNA 残留、重复率、支原体污染。
只读 _work/server/ 下从 R1 同步回来的统计；不删除、不过滤任何 PEDV reads。"""
# >>> release-candidate path shim (added during repository curation; no analytical change) >>>
import os as _os, sys as _sys
_d = _os.path.dirname(_os.path.abspath(__file__))
while _d != _os.path.dirname(_d) and not _os.path.isdir(_os.path.join(_d, "config")):
    _d = _os.path.dirname(_d)
_sys.path.insert(0, _os.path.join(_d, "config"))
from pedv_paths import PEDV_P, PEDV_K, PEDV_R, PEDV_COV, PEDV_SERVER_ROOT, REPO_ROOT, SCRIPT_DIRS, ORIGIN
_sys.path.extend(_x for _x in SCRIPT_DIRS if _x not in _sys.path)  # sibling modules of other repo modules (own directory keeps priority)
__REPO_FILE__ = globals().get("__REPO_FILE__") or _os.path.abspath(__file__)  # real location of this script inside the repository (sibling-script loading must use this, not the virtual __file__)
_o = ORIGIN.get(_os.path.relpath(_os.path.abspath(__file__), REPO_ROOT).replace(chr(92), "/"))
if _o:
    __file__ = _o  # virtual ORIGINAL location inside the working tree: keeps the scripts' own dirname(__file__)-relative path logic unchanged
# <<< release-candidate path shim <<<
import os, sys, csv, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

W = WORK
OLD_SUM = {r["seq_id"]: r for r in rd(os.path.join(os.path.dirname(OUT), "batch23_riboseq_qc_mapping_output/batch23_riboseq_qc_mapping_summary.tsv"))}
VAL = [r for r in rd(f"{RB}/validation/old_vs_new_star_metrics.tsv") if r["data"] == "RIBO"]

def aln(path):
    return {r["domain"]: {k: int(v) for k, v in r.items() if k != "domain"} for r in rd(path)}
def dupt(path):
    return {r["domain"]: (int(r["reads"]), int(r["distinct_5p_keys"])) for r in rd(path)}
def qcf(path):
    d = {}
    for l in open(path, encoding="utf-8"):
        x = l.rstrip("\n").split("\t"); d[x[0]] = x[1] if len(x) > 1 else ""
    return d
def mm_class(path):
    d = collections.Counter()
    for r in rd(path): d[(r["nh_class"], r["ref_class"])] += int(r["reads"])
    return d
pct = lambda a, b: f"{100*a/b:.3f}" if b else "NA"
f4 = lambda x: f"{x:.4f}"

BIOTYPE = {r["gene_id"]: r["gene_biotype"] for r in rd(f"{W}/host_gene_biotype.tsv")}
summ, comp, rlrows, lpf, contam = [], [], [], [], []
for sid in IDS:
    m = mcols(sid)
    fin = star_final(f"{W}/align/{sid}/{sid}.Log.final.out")
    A = aln(f"{W}/stats/{sid}/aln_metrics.tsv"); D = dupt(f"{W}/stats/{sid}/dup.tsv"); Q = qcf(f"{W}/stats/{sid}/qc.tsv")
    META[sid]["bam_size"] = Q["bam_size"]
    inp = int(fin["Number of input reads"]); uniq = int(fin["Uniquely mapped reads number"])
    many = int(fin["Number of reads mapped to too many loci"]); multi = int(fin["Number of reads mapped to multiple loci"])
    short = int(fin["Number of reads unmapped: too short"]); other = int(fin["Number of reads unmapped: other"])
    mism = int(fin["Number of reads unmapped: too many mismatches"])
    h, vp, vm, va = A["host"], A["pedv_plus"], A["pedv_minus"], A["pedv_all"]
    assert h["reads"] + va["reads"] == uniq, (sid, h["reads"], va["reads"], uniq)
    MM = mm_class(f"{W}/multimap/{sid}.mm100.class.tsv")
    mm_tot = sum(v for (a, b), v in MM.items() if a == "NH_gt1")
    row = {**m,
           "input_reads(去rRNA后)": inp, "uniquely_mapped": uniq, "uniquely_mapped_pct": pct(uniq, inp),
           "multimapped_2_to_1(multiple loci)": multi, "mapped_too_many_loci(>1 locus; outFilterMultimapNmax=1)": many,
           "too_many_loci_pct": pct(many, inp), "unmapped_too_short": short, "unmapped_too_short_pct": pct(short, inp),
           "unmapped_other": other, "unmapped_too_many_mismatches": mism,
           "host_unique_reads": h["reads"], "PEDV_unique_reads": va["reads"], "PEDV_plus_reads": vp["reads"], "PEDV_minus_reads": vm["reads"],
           "PEDV_pct_of_unique": pct(va["reads"], uniq),
           "PEDV_mismatch_per_100_aligned_bases(nM)": f4(100 * va["nM_sum"] / va["aligned_M_bases"]),
           "PEDV_mean_mismatch_per_read(nM)": f4(va["nM_sum"] / va["reads"]),
           "PEDV_reads_with_mismatch_frac": f4(va["reads_with_mismatch"] / va["reads"]),
           "PEDV_softclip_read_frac": f4(va["softclip_reads"] / va["reads"]),
           "PEDV_softclip5_read_frac": f4(va["softclip5_reads"] / va["reads"]),
           "PEDV_softclip_bases_per_read": f4(va["softclip_bases"] / va["reads"]),
           "PEDV_indel_read_frac": f4(va["indel_reads"] / va["reads"]),
           "PEDV_spliced_reads": va["spliced_reads"],
           "host_mismatch_per_100_aligned_bases(nM)": f4(100 * h["nM_sum"] / h["aligned_M_bases"]),
           "host_softclip_read_frac": f4(h["softclip_reads"] / h["reads"]),
           "STAR_mismatch_rate_per_base_all": fin["Mismatch rate per base, %"],
           "NH_gt1_in_BAM(应为0)": h["NH_gt1_reads"] + va["NH_gt1_reads"],
           "diag_multimap_reads_total(NH>1, mm100)": mm_tot,
           "diag_multimap_host_only": MM[("NH_gt1", "host_only")], "diag_multimap_PEDV_only": MM[("NH_gt1", "PEDV_only")],
           "diag_multimap_host_and_PEDV": MM[("NH_gt1", "host_and_PEDV")],
           "diag_multimap_host_and_PEDV_pct_of_PEDV_unique": pct(MM[("NH_gt1", "host_and_PEDV")], va["reads"]),
           "diag_rescued_unique_from_unmapped(NH=1, mm100; 主结果不含)": sum(v for (a, b), v in MM.items() if a == "NH_eq1"),
           "dup_host_frac": f4(1 - D["host"][1] / D["host"][0]), "dup_PEDV_frac": f4(1 - D["pedv_all"][1] / D["pedv_all"][0]),
           "PEDV_distinct_5p_keys(去重敏感性)": D["pedv_all"][1], "host_distinct_5p_keys(去重敏感性)": D["host"][1],
           "bam_size_bytes": Q["bam_size"]}
    summ.append(row)

    # ---- 新旧参考对比（旧参考数字取既有交付 + 旧 BAM 只读统计；不重跑旧流程） ----
    o = OLD_SUM[sid]; ofin = star_final(f"{W}/oldref/{sid}.Log.final.out"); OA = aln(f"{W}/stats_oldref/{sid}/aln_metrics.tsv")["pedv_all"]
    old_v = int(o["pedv_mapped_reads"]); old_h = int(o["host_mapped_reads"]); old_u = int(o["star_uniquely_mapped_reads"])
    assert OA["reads"] == old_v, (sid, OA["reads"], old_v)
    c = {**m,
         "old_PEDV_unique": old_v, "GIb_PEDV_unique": va["reads"], "PEDV_diff": va["reads"] - old_v, "PEDV_ratio_GIb_over_old": f4(va["reads"] / old_v),
         "old_host_unique": old_h, "GIb_host_unique": h["reads"], "host_diff": h["reads"] - old_h, "host_ratio": f4(h["reads"] / old_h),
         "old_total_unique": old_u, "GIb_total_unique": uniq, "total_unique_ratio": f4(uniq / old_u),
         "old_PEDV_mean_mismatch_per_read": f4(OA["nM_sum"] / OA["reads"]), "GIb_PEDV_mean_mismatch_per_read": row["PEDV_mean_mismatch_per_read(nM)"],
         "old_PEDV_mismatch_per_100nt": f4(100 * OA["nM_sum"] / OA["aligned_M_bases"]), "GIb_PEDV_mismatch_per_100nt": row["PEDV_mismatch_per_100_aligned_bases(nM)"],
         "old_PEDV_reads_with_mismatch_frac": f4(OA["reads_with_mismatch"] / OA["reads"]), "GIb_PEDV_reads_with_mismatch_frac": row["PEDV_reads_with_mismatch_frac"],
         "old_PEDV_softclip_read_frac": f4(OA["softclip_reads"] / OA["reads"]), "GIb_PEDV_softclip_read_frac": row["PEDV_softclip_read_frac"],
         "old_PEDV_softclip_bases_per_read": f4(OA["softclip_bases"] / OA["reads"]), "GIb_PEDV_softclip_bases_per_read": row["PEDV_softclip_bases_per_read"],
         "old_STAR_mismatch_rate_all": ofin["Mismatch rate per base, %"], "GIb_STAR_mismatch_rate_all": fin["Mismatch rate per base, %"],
         "old_too_many_loci": int(ofin["Number of reads mapped to too many loci"]), "GIb_too_many_loci": many,
         "old_too_short": int(ofin["Number of reads unmapped: too short"]), "GIb_too_short": short,
         "old_source": "旧病毒参考（该样本当时所用）：PEDV/宿主/总唯一数取既有 Ribo QC 交付汇总表；错配/soft-clip 为旧 BAM 只读统计（同一脚本 virus_only 模式）"}
    v = [x for x in VAL if x["sample"] == sid]
    if v:
        vn = [x for x in v if x["reference"].startswith(("new", "v11"))]; vo = [x for x in v if x["reference"].startswith("old")]
        vn = [x for x in vn if x["reference"] == "new_GIb"] or vn
        if vn and vo:
            vn, vo = vn[0], vo[0]
            c["refbuild_subsample_3M_PEDV_ratio_new_over_old"] = f4(int(vn["virus_unique_fragments"]) / int(vo["virus_unique_fragments"]))
            c["refbuild_subsample_3M_mismatch_per_read_old_to_new"] = f"{vo['virus_mean_mismatch_per_read']}→{vn['virus_mean_mismatch_per_read']}"
            c["refbuild_subsample_3M_softclip_frac_old_to_new"] = f"{vo['virus_softclip_read_frac']}→{vn['virus_softclip_read_frac']}"
            c["refbuild_subsample_3M_host_ratio_new_over_old"] = f4(int(vn["host_unique_fragments"]) / int(vo["host_unique_fragments"]))
    comp.append(c)

    # ---- 读长分布（宿主 / PEDV+ / PEDV− 分开） ----
    R = rd(f"{W}/stats/{sid}/readlen_aln.tsv")
    ht = sum(int(r["host_n"]) for r in R); vt = sum(int(r["pedv_plus_n"]) for r in R); mt = sum(int(r["pedv_minus_n"]) for r in R)
    for r in R:
        L = int(r["L"])
        rlrows.append({**m, "aligned_length": L if L < 150 else "≥150", "length_zone": zone(L),
                       "host_reads": r["host_n"], "host_pct": pct(int(r["host_n"]), ht),
                       "PEDV_plus_reads": r["pedv_plus_n"], "PEDV_plus_pct": pct(int(r["pedv_plus_n"]), vt),
                       "PEDV_minus_reads": r["pedv_minus_n"], "PEDV_minus_pct": pct(int(r["pedv_minus_n"]), mt)})
    sm = lambda col, lo, hi: sum(int(r[col]) for r in R if lo <= int(r["L"]) <= hi)
    BIG = 10 ** 6
    lpf.append({**m, "host_total": ht, "host_LPF_ge35": sm("host_n", 35, BIG), "host_LPF_ge35_pct": pct(sm("host_n", 35, BIG), ht),
                "host_LPF_35_36_pct": pct(sm("host_n", 35, 36), ht), "host_LPF_ge37_pct": pct(sm("host_n", 37, BIG), ht),
                "PEDV_plus_total": vt, "PEDV_plus_LPF_ge35": sm("pedv_plus_n", 35, BIG), "PEDV_plus_LPF_ge35_pct": pct(sm("pedv_plus_n", 35, BIG), vt),
                "PEDV_plus_LPF_35_36_pct": pct(sm("pedv_plus_n", 35, 36), vt), "PEDV_plus_LPF_ge37_pct": pct(sm("pedv_plus_n", 37, BIG), vt),
                "PEDV_minus_total(单列)": mt, "PEDV_minus_LPF_ge35": sm("pedv_minus_n", 35, BIG),
                **{f"host_zone_{z}_pct": pct(sm("host_n", lo, hi), ht) for z, lo, hi in (("lt23", 0, 22), ("23_28", 23, 28), ("29_30", 29, 30), ("31_34", 31, 34))},
                **{f"PEDV_plus_zone_{z}_pct": pct(sm("pedv_plus_n", lo, hi), vt) for z, lo, hi in (("lt23", 0, 22), ("23_28", 23, 28), ("29_30", 29, 30), ("31_34", 31, 34))},
                "status": "LPF 全部保留、不删除，只统计"})

    # ---- 污染 / rRNA 残留 / 重复 ----
    bt = open(f"{W}/myco/{sid}.bt2.log", encoding="utf-8").read().split("\n")
    nin = int(bt[0].split()[0]); n0 = int(bt[2].split()[0]); nal = nin - n0
    idx = {l.split("\t")[0]: int(l.split("\t")[2]) for l in open(f"{W}/myco/{sid}.idxstats.tsv") if l.split("\t")[0] != "*"}
    q20 = {l.split("\t")[0]: int(l.split("\t")[2]) for l in open(f"{W}/myco/{sid}.q20.idxstats.tsv") if l.split("\t")[0] != "*"}
    sp = collections.Counter(); sp20 = collections.Counter()
    for k, n in idx.items(): sp[k.split("__")[0]] += n
    for k, n in q20.items(): sp20[k.split("__")[0]] += n
    q20t = sum(sp20.values())
    rp = collections.Counter(); rtot = 0
    for l in open(f"{W}/align/{sid}/{sid}.ReadsPerGene.out.tab"):
        g, u, *_ = l.rstrip("\n").split("\t")
        if g.startswith("N_"): continue
        b = BIOTYPE.get(g, "PEDV_or_unknown"); rp[b] += int(u)
    rr = rp["rRNA"] + rp["Mt_rRNA"]
    contam.append({**m, "denominator_reads(去rRNA后，与6月检测同口径)": nin,
                   "myco_aligned_reads(all)": nal, "myco_pct(all)": f"{100*nal/nin:.4f}",
                   "myco_MAPQ20_reads": q20t, "myco_MAPQ20_pct": f"{100*q20t/nin:.4f}",
                   "Mhyorhinis_share_of_MAPQ20": f"{sp20['Mhyorhinis']/q20t:.3f}" if q20t else "NA",
                   "species_MAPQ20": ";".join(f"{k}:{v}" for k, v in sp20.most_common() if v),
                   "myco_pct_of_STAR_input_unmapped_pool": f"{100*nal/(inp-uniq):.4f}",
                   "normalization_denominator_affected": "否：本任务 RPM 分母 = 该样本 STAR 唯一比对 reads（宿主+PEDV），支原体 reads 不进入；"
                                                         f"若改用去rRNA后 total reads 为分母，会被抬高 {100*nal/nin:.3f}%（all）/ {100*q20t/nin:.3f}%（MAPQ≥20）",
                   "host_unique_on_rRNA_genes(5S/5.8S)": rp["rRNA"], "host_unique_on_Mt_rRNA": rp["Mt_rRNA"],
                   "rRNA_residual_pct_of_host_unique(仅可测部分)": f"{100*rr/h['reads']:.4f}",
                   "rRNA_18S_28S_residual": "待确认：去rRNA索引与宿主GTF均不含18S/28S序列，本任务未测",
                   "dup_host_frac": row["dup_host_frac"], "dup_PEDV_frac": row["dup_PEDV_frac"],
                   "bowtie2": "2.5.4 默认参数 --no-unal（6月检测为 2.3.5.1，同参数）"})

prov = provenance("riboseq_GIb", [
    "inputs: STAR Log.final.out / ReadsPerGene / idxstats（align/<id>/）；gib_bam_stats.py 输出（stats/<id>/）；multimap_diag/<id>/*.class.tsv；qc_contamination/",
    "STAR 2.7.11b: --runThreadN 8 --readFilesCommand 'pigz -dc' --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 20000000000 "
    "--outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --quantMode GeneCounts (+ --outReadsUnmapped Fastx --outTmpDir，纯输出选项)",
    "index: refdata_GIb_v1.1/merged/star_ChlSab1.1_PEDV_GIb_oh40；input = 已去接头、去rRNA的 reads（09-25 QC 产出，未重做）",
    "mismatch = STAR nM 标签（SE 即每条 read 错配数）；aligned bases = CIGAR M/=/X；soft-clip = CIGAR S",
    "scripts: scripts/build_A_alignment_qc.py（本地）；R1 scripts/run_star_one.sh, gib_bam_stats.py, run_multimap_diag.sh, run_myco.sh"])
write(f"{OUT}/riboseq_GIb_alignment_summary.tsv", summ, prov.replace("riboseq_GIb", "riboseq_GIb_alignment_summary"))
write(f"{OUT}/riboseq_ref_comparison.tsv", comp, prov.replace("riboseq_GIb", "riboseq_ref_comparison") +
      "# 旧参考数字只作对照，不作为 GI-b 结论；refbuild_subsample_* 列 = 建参考时 300 万 reads 抽样验证（refbuild_GIb/validation）\n")
write(f"{OUT}/riboseq_readlength_distribution.tsv", rlrows, prov.replace("riboseq_GIb", "riboseq_readlength_distribution") +
      "# 读长 = 比对区块长度（CIGAR M/=/X/D，跳过 N，不含 soft-clip）；宿主与 PEDV 分开；PEDV− 单列；全部读长保留\n")
write(f"{OUT}/riboseq_LPF_summary.tsv", lpf, prov.replace("riboseq_GIb", "riboseq_LPF_summary"))
write(f"{OUT}/riboseq_qc_contamination.tsv", contam, prov.replace("riboseq_GIb", "riboseq_qc_contamination") +
      "# 支原体：8 物种 panel（mycoplasma_panel.fa md5 b20a11d9…），bowtie2 比对去 rRNA 后 reads；不因污染剔除任何样本\n")
print("A done")
for r in comp: print(r["sample_id"], r["PEDV_ratio_GIb_over_old"], r["old_PEDV_mean_mismatch_per_read"], "->", r["GIb_PEDV_mean_mismatch_per_read"],
                     r["old_PEDV_softclip_read_frac"], "->", r["GIb_PEDV_softclip_read_frac"], r["host_ratio"], r.get("refbuild_subsample_3M_PEDV_ratio_new_over_old", ""))
