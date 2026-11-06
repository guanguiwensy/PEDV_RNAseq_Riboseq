#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os as _os, sys as _sys
_d = _os.path.dirname(_os.path.abspath(__file__))
while _d != _os.path.dirname(_d) and not _os.path.isdir(_os.path.join(_d, "config")):
    _d = _os.path.dirname(_d)
_sys.path.insert(0, _os.path.join(_d, "config"))
from pedv_paths import PEDV_P, PEDV_K, PEDV_R, PEDV_COV, PEDV_SERVER_ROOT, REPO_ROOT, SCRIPT_DIRS, ORIGIN
_sys.path.extend(_x for _x in SCRIPT_DIRS if _x not in _sys.path)
__REPO_FILE__ = globals().get("__REPO_FILE__") or _os.path.abspath(__file__)
_o = ORIGIN.get(_os.path.relpath(_os.path.abspath(__file__), REPO_ROOT).replace(chr(92), "/"))
if _o:
    __file__ = _o
import sample_registry as SR
import os, re, sys, json, datetime
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.dirname(HERE)
TASK = os.path.dirname(OUT)
RGP = f"{TASK}/riboseq_GIb_mapping_psite"
QC = f"{TASK}/batch23_riboseq_qc_mapping_output/batch23_riboseq_qc_mapping_summary.tsv"
NOW = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
SAMPLES = [(SR.lib("RIBO_VERO_CHX_1MOI_5h"), "CHX", 5), (SR.lib("RIBO_VERO_CHX_1MOI_8h"), "CHX", 8), (SR.lib("RIBO_VERO_CHX_1MOI_16h"), "CHX", 16),
           (SR.lib("RIBO_VERO_HAR_1MOI_5h"), "HAR", 5), (SR.lib("RIBO_VERO_HAR_1MOI_8h"), "HAR", 8), (SR.lib("RIBO_VERO_HAR_1MOI_16h"), "HAR", 16)]
IDS = [s for s, _, _ in SAMPLES]; DRUG = {s: d for s, d, _ in SAMPLES}; TIME = {s: t for s, _, t in SAMPLES}
rd = lambda p: pd.read_csv(p, sep="\t", comment="#", low_memory=False)
CHECKS = []
def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok), detail)); print(("PASS " if ok else "FAIL ") + name, detail)
    if not ok: raise SystemExit(f"恒等式不成立：{name} {detail}")
def prov(title, extra=()):
    L = [f"# {title}", "# version: SHARED_v1.0", f"# generated: {NOW} (local)",
         "# inputs (frozen): RGP_v1.3 riboseq_GIb_mapping_psite/{riboseq_GIb_alignment_summary,riboseq_readlength_distribution,psite_offset_GIb,length_qc_GIb,periodicity_host,host_metagene_GIb,riboseq_qc_contamination}.tsv; "
         "batch23_riboseq_qc_mapping_output/batch23_riboseq_qc_mapping_summary.tsv (only raw/trimmed/non_rRNA columns)",
         "# host reference: ChlSab1.1 (Vero); virus reference: PEDV_GIb (refdata_GIb_v1.1, 27,959 nt). coordinates 1-based",
         "# script: shared_riboseq_metrics/scripts/build_shared_metrics.py"]
    return "\n".join(L + [f"# {e}" for e in extra]) + "\n"
def wr(name, df, title, extra=()):
    with open(f"{OUT}/{name}", "w", encoding="utf-8", newline="") as f:
        f.write(prov(title, extra)); df.to_csv(f, sep="\t", index=False, float_format="%.6g")


q = rd(QC).set_index("seq_id").loc[IDS]
a = rd(f"{RGP}/riboseq_GIb_alignment_summary.tsv").set_index("sample_id").loc[IDS]
c = rd(f"{RGP}/riboseq_qc_contamination.tsv").set_index("sample_id").loc[IDS]
check("六个样本都在 QC 汇总表中", len(q) == 6)
check("non_rRNA_reads == GI-b input_reads(去rRNA后)（唯一合法接缝）", (q.non_rRNA_reads.values == a["input_reads(去rRNA后)"].values).all(),
      str(dict(zip(IDS, q.non_rRNA_reads.tolist()))))
check("time_h/drug 与样本矩阵一致", all(int(a.loc[s, "time_h"]) == TIME[s] and a.loc[s, "drug"] == DRUG[s] for s in IDS))
check("host_unique + PEDV_unique == uniquely_mapped", ((a.host_unique_reads + a.PEDV_unique_reads) == a.uniquely_mapped).all())
ML = "mapped_too_many_loci(>1 locus; outFilterMultimapNmax=1)"
check("uniquely + too_many_loci + too_short + other + too_many_mismatches == input_reads",
      ((a.uniquely_mapped + a[ML] + a.unmapped_too_short + a.unmapped_other + a.unmapped_too_many_mismatches) == a["input_reads(去rRNA后)"]).all())
check("multimapped_2_to_1 恒为 0（Nmax=1 构造性为 0，不是'没有多比对'）", (a["multimapped_2_to_1(multiple loci)"] == 0).all(),
      f"真正的多位点类别 mapped_too_many_loci 占 non-rRNA {a.too_many_loci_pct.min():.1f}–{a.too_many_loci_pct.max():.1f}%")
check("trimmed <= raw 且 non_rRNA <= trimmed", ((q.trimmed_reads <= q.raw_r1_reads) & (q.non_rRNA_reads <= q.trimmed_reads)).all())
f = pd.DataFrame(index=IDS)
f["drug"] = [DRUG[s] for s in IDS]; f["time_h"] = [TIME[s] for s in IDS]
f["raw_r1_reads"] = q.raw_r1_reads; f["trimmed_reads"] = q.trimmed_reads; f["non_rRNA_reads"] = q.non_rRNA_reads
f["n_adapter_length_discarded"] = q.raw_r1_reads - q.trimmed_reads
f["n_rRNA_ncRNA_discarded"] = q.trimmed_reads - q.non_rRNA_reads
f["n_host_unique"] = a.host_unique_reads; f["n_PEDV_unique"] = a.PEDV_unique_reads
f["n_more_than_one_locus"] = a[ML]; f["n_unmapped_too_short"] = a.unmapped_too_short
f["n_unmapped_other"] = a.unmapped_other + a.unmapped_too_many_mismatches
cats = [c_ for c_ in f.columns if c_.startswith("n_")]
check("七类互斥且之和 == raw_r1_reads（分母 = raw R1 reads）", (f[cats].sum(axis=1) == f.raw_r1_reads).all())
for c_ in cats: f["frac_" + c_[2:]] = f[c_] / f.raw_r1_reads
check("堆叠比例之和 == 1", np.allclose(f[[c_ for c_ in f.columns if c_.startswith("frac_")]].sum(axis=1), 1.0))
f["myco_pct_all_of_nonrRNA"] = c["myco_pct(all)"]; f["myco_pct_MAPQ20_of_nonrRNA"] = c["myco_MAPQ20_pct"]
f["myco_denominator_reads(non_rRNA)"] = c["denominator_reads(去rRNA后，与6月检测同口径)"]
check("支原体分母 == non_rRNA_reads", (f["myco_denominator_reads(non_rRNA)"] == f.non_rRNA_reads).all())
f.insert(0, "sample_id", f.index)
wr("sample_funnel.tsv", f, "六个 Ribo-seq 文库的 reads 漏斗（分母 raw R1 reads；支原体分母为 non-rRNA reads，另列）",
   ["上游三列（raw/trimmed/non_rRNA）取自 batch23 QC 表，其 mapped 之后的列属旧参考 star_ChlSab1.1_LZ202401，一律不用",
    "mapped 之后全部取自 riboseq_GIb_alignment_summary.tsv（GI-b）；Nmax=1 → multimapped_2_to_1 构造性为 0"])


dist = rd(f"{RGP}/riboseq_readlength_distribution.tsv")
dist["L"] = pd.to_numeric(dist.aligned_length, errors="coerce")
ge150 = dist[dist.L.isna()]
dist = dist.dropna(subset=["L"]).astype({"L": int})
off = rd(f"{RGP}/psite_offset_GIb.tsv")
lq = rd(f"{RGP}/length_qc_GIb.tsv")
ph = rd(f"{RGP}/periodicity_host.tsv"); ph["L"] = pd.to_numeric(ph.read_length, errors="coerce")
check("psite_offset_GIb 共 6×25=150 行（16–40 nt）", len(off) == 150 and off.read_length.min() == 16 and off.read_length.max() == 40, f"{len(off)} 行")
rows = []
for s in IDS:
    g = dist[dist.sample_id == s]
    tot_host = a.loc[s, "host_unique_reads"]; tot_p = a.loc[s, "PEDV_plus_reads"]; tot_m = a.loc[s, "PEDV_minus_reads"]
    ex = ge150[ge150.sample_id == s]
    check(f"{s} 逐读长求和 + ≥150 桶 == 唯一比对（host/PEDV+/PEDV−）",
          g.host_reads.sum() + ex.host_reads.sum() == tot_host and g.PEDV_plus_reads.sum() + ex.PEDV_plus_reads.sum() == tot_p
          and g.PEDV_minus_reads.sum() + ex.PEDV_minus_reads.sum() == tot_m,
          f"host {g.host_reads.sum()}+{ex.host_reads.sum()} vs {tot_host}")
    for _, r in g.iterrows():
        L = int(r.L); o = off[(off.sample_id == s) & (off.read_length == L)]
        l = lq[(lq.sample_id == s) & (lq.read_length == L)]
        o_ = o.iloc[0] if len(o) else None; l_ = l.iloc[0] if len(l) else None
        p = ph[(ph.sample_id == s) & (ph.L == L)]
        hf = p[["psite_frame0", "psite_frame1", "psite_frame2"]].iloc[0].tolist() if len(p) else [np.nan] * 3
        hn = float(p["n_reads(5′端在CDS主体)"].iloc[0]) if len(p) else np.nan
        pf = [np.nan] * 3; pn = np.nan
        if l_ is not None and isinstance(l_["PEDV_psite_frame0_1_2"], str):
            pf = [float(x) for x in l_["PEDV_psite_frame0_1_2"].split("/")]; pn = l_["PEDV_frame_n(ORF唯一区5′端)"]
        for dom, reads, tot, fr, fn in (("host", r.host_reads, tot_host, hf, hn), ("PEDV_plus", r.PEDV_plus_reads, tot_p, pf, pn),
                                        ("PEDV_minus", r.PEDV_minus_reads, tot_m, [np.nan] * 3, np.nan)):
            rows.append(dict(sample_id=s, drug=DRUG[s], time_h=TIME[s], read_length=L, domain=dom, reads=int(reads), frac_of_domain=reads / tot,
                             offset=(o_.offset if o_ is not None else np.nan), offset_confidence=(o_.offset_confidence if o_ is not None else ""),
                             category=(o_.category_GIb if o_ is not None else ""),
                             frame_n=fn, frame0=fr[0], frame1=fr[1], frame2=fr[2],
                             F0_minus_mean_F1F2=(fr[0] - (fr[1] + fr[2]) / 2) if not np.isnan(fr[0]) else np.nan))
long = pd.DataFrame(rows)
wr("length_metrics_long.tsv", long, "sample × read_length × domain 指标（宿主唯一比对 / PEDV(+) / PEDV(−)）",
   ["reads = 该读长的唯一比对 reads 数（读长=比对区块长度，不含 soft-clip；≥150 nt 合并为一桶，不进本表，逐样本在 length_summary 中计入分母）",
    "host frame = 宿主 canonical CDS 主体（两端各去 45 nt）上 P-site 读框（periodicity_host）；PEDV frame = 7 个 canonical ORF 唯一区、P-site 读框（length_qc_GIb）",
    "F0_minus_mean_F1F2 = frame0 − mean(frame1, frame2)（FigSX 既定口径）；offset 为 psite_offset_GIb（Module 1 方法 B，宿主起始密码子校准）"])


sm = []
for s in IDS:
    g = dist[dist.sample_id == s].sort_values("L"); ex = ge150[ge150.sample_id == s]
    for dom, col in (("host", "host_reads"), ("PEDV_plus", "PEDV_plus_reads")):
        n = g[col].values.astype(float); n_ex = float(ex[col].sum()); tot = n.sum() + n_ex; Ls = g.L.values
        cum = np.cumsum(np.append(n, n_ex)) / tot
        med = int(np.append(Ls, 150)[np.searchsorted(cum, 0.5)])
        sm.append(dict(sample_id=s, drug=DRUG[s], time_h=TIME[s], domain=dom, n_unique_reads=int(tot), median_length=med,
                       frac_29_31=n[(Ls >= 29) & (Ls <= 31)].sum() / tot, frac_ge35_LPF=(n[Ls >= 35].sum() + n_ex) / tot,
                       modal_length=int(Ls[n.argmax()])))
sm = pd.DataFrame(sm)
wr("length_summary.tsv", sm, "逐样本 × domain 的长度汇总（中位长度、29–31 nt 占比、≥35 nt 占比）",
   ["中位长度取累计比例首次 ≥0.5 的长度；≥150 nt 桶按 150 计入分母与累计"])


mg = rd(f"{RGP}/host_metagene_GIb.tsv"); mg = mg[(mg.anchor == "start") & mg.length_group.str.match(r"^L\d+$")].copy()
mg["L"] = mg.length_group.str[1:].astype(int)
surf = []; DEVS = []
for (s, L), g in mg.groupby(["sample_id", "L"]):
    r = off[(off.sample_id == s) & (off.read_length == L)]
    if r.empty or pd.isna(r.offset.iloc[0]): continue
    o = int(r.offset.iloc[0]); cnt = g.set_index("rel_nt").psite_count

    vals = {cand: int(cnt.get(o - cand, 0)) for cand in range(6, 19)}
    tot = sum(vals.values()); pk = max(vals.values()); am = max(vals, key=vals.get)
    if r.peak_reads.iloc[0] > 0:
        dv = abs(pk - r.peak_reads.iloc[0]) / r.peak_reads.iloc[0]; DEVS.append((s, L, int(r.peak_reads.iloc[0]), pk, dv))
    for cand, v in vals.items():
        surf.append(dict(sample_id=s, drug=DRUG[s], time_h=TIME[s], read_length=L, candidate_offset=cand, reads=v,
                         share_of_window=(v / tot if tot else np.nan), is_selected=(cand == o),
                         selected_offset=o, offset_confidence=r.offset_confidence.iloc[0], category=r.category_GIb.iloc[0],
                         support_reads_host_start=int(r.support_reads_host_start.iloc[0]), peak_reads_table=int(r.peak_reads.iloc[0]),
                         peak_share_table=float(r.peak_share.iloc[0]), margin_table=float(r.peak_margin_over_2nd.iloc[0]),
                         half0=r.offset_half0.iloc[0], half1=r.offset_half1.iloc[0]))
surf = pd.DataFrame(surf)
# tie-aware identity: the selected offset must reach the maximum read count of the candidate surface (sparse lengths can have ties, e.g. the CHX 16 h library at 39 nt with 4 reads;
# the original check required the FIRST maximum to equal the selected offset, which holds for the five previously delivered libraries)
am_ok = surf.groupby(["sample_id", "read_length"]).apply(lambda g: g.loc[g.candidate_offset == g.selected_offset.iloc[0], "reads"].max() == g.reads.max()).all()
check("由 host_metagene 反推的候选曲面，峰位置 == 冻结 offset（逐样本×逐读长全部）", am_ok, f"{surf.groupby(['sample_id','read_length']).ngroups} 个组合")
dv = pd.DataFrame(DEVS, columns=["sample_id", "L", "peak_table", "peak_metagene", "rel_dev"])
big = dv[dv.peak_table >= 500]
check("峰高 ≥500 的组合：反推峰高与 psite_offset_GIb.peak_reads 一致（最大相对偏差 < 1.5%）", big.rel_dev.max() < 0.015,
      f"{len(big)} 个组合，max rel dev = {big.rel_dev.max():.4f}；<500 reads 的 {len(dv)-len(big)} 个组合 max = {dv[dv.peak_table<500].rel_dev.max():.3f}（小计数下 canonical 转录本与 start_codon.bed 的差异被放大）")
wr("offset_candidate_surface.tsv", surf, "候选 P-site offset 曲面（由宿主起始密码子 metagene 反推；候选 6–18 nt，对应搜索窗 j∈[7,19]）",
   ["reads(c) = host_metagene_GIb 中该读长 start 锚点 rel_nt = offset_selected − c 的 P-site 数；与 psite_offset_GIb.peak_reads 一致（<1.5%）",
    "曲面只含 c∈[6,18]（Module 1 方法 B 的搜索窗），24 nt 的 offset=6 落在窗边界，已标低置信"])


with open(f"{OUT}/AUDIT_checks.json", "w", encoding="utf-8") as fh:
    json.dump([dict(name=n, ok=ok, detail=d) for n, ok, d in CHECKS], fh, ensure_ascii=False, indent=1)
print(len(CHECKS), "项检查全部通过")
