#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Sup Fig 13：全部源数据（纯本地；RGP_v1.3 冻结表 + stats/ 统计 + 第二批 RIBO_VERO_CHX_1MOI_16h 的 stats）。
- 六个文库：CHX 5/8/16 h = RIBO_VERO_CHX_1MOI_5h / RIBO_VERO_CHX_1MOI_8h / RIBO_VERO_CHX_1MOI_16h，HAR 5/8/16 h = RIBO_VERO_HAR_1MOI_5h / RIBO_VERO_HAR_1MOI_8h / RIBO_VERO_HAR_1MOI_16h（每条件 1 个库，时间点非重复）。
- offset：Module 1 方法 B 判定函数在 stats/*/host_startpos*.tsv 上逐库重新推导（与已交付 psite_offset_GIb 在五个库 121/121 一致）；不改动已采用的 offset。
- 病毒读长集合（主）：host_criteria = 各库 25–35 nt 中通过宿主判据的读长；strict_core / 29–31 nt / all_accepted 作补充口径。
- ORF 级 P-site 计数用真实 P-site 位置（+ 链），ORF 唯一区 = [起点+45, 终点−45] 去掉与任何其他 ORF 重叠的碱基；frame = (P-site − ORF 起点) mod 3。"""
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
import sample_registry as SR
import os, gzip, csv, json, collections, datetime
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.dirname(HERE); TASK = os.path.dirname(OUT)
RGP = f"{TASK}/riboseq_GIb_mapping_psite"; WORK = f"{RGP}/_work/server"; SD = f"{OUT}/source_data"; os.makedirs(SD, exist_ok=True)
NOW = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
SAMPLES = [(SR.lib("RIBO_VERO_CHX_1MOI_5h"), "CHX", 5), (SR.lib("RIBO_VERO_CHX_1MOI_8h"), "CHX", 8), (SR.lib("RIBO_VERO_CHX_1MOI_16h"), "CHX", 16),
           (SR.lib("RIBO_VERO_HAR_1MOI_5h"), "HAR", 5), (SR.lib("RIBO_VERO_HAR_1MOI_8h"), "HAR", 8), (SR.lib("RIBO_VERO_HAR_1MOI_16h"), "HAR", 16)]
IDS = [s for s, _, _ in SAMPLES]; DRUG = {s: d for s, d, _ in SAMPLES}; TIME = {s: t for s, _, t in SAMPLES}; LAB = {s: f"{d} {t} h" for s, d, t in SAMPLES}
ORFS = [("ORF1a", 293, 12622), ("ORF1b", 12592, 20613), ("S", 20610, 24758), ("ORF3", 24758, 25033), ("E", 25361, 25600), ("M", 25608, 26288), ("N", 26300, 27625)]
GLEN = 27959; TRIM = 45; J_LO, J_HI = 7, 19
BASE_F0 = 0.80
LOW, MISSING = 1000, 100
rdt = lambda p: pd.read_csv(p, sep="\t", comment="#", low_memory=False)

def wr(name, df, title, extra=()):
    head = [f"# {title}", "# version: SUP13_v1.0", f"# generated: {NOW} (local)",
            f"# libraries: CHX 5 h {SR.lib('RIBO_VERO_CHX_1MOI_5h')}, CHX 8 h {SR.lib('RIBO_VERO_CHX_1MOI_8h')}, CHX 16 h {SR.lib('RIBO_VERO_CHX_1MOI_16h')}, HAR 5 h {SR.lib('RIBO_VERO_HAR_1MOI_5h')}, HAR 8 h {SR.lib('RIBO_VERO_HAR_1MOI_8h')}, HAR 16 h {SR.lib('RIBO_VERO_HAR_1MOI_16h')}; one library per condition (time points are not replicates)",
            "# virus reference PEDV_GIb (27,959 nt, refdata_GIb_v1.1, PEDV_GIb.fa sha256 48c36027c9b24ddcdbbec34537dea126ade378ad4094e5d72054f697c8412b92); host ChlSab1.1; coordinates 1-based",
            "# script: code/prepare_panels.py"] + [f"# {e}" for e in extra]
    with open(f"{SD}/{name}", "w", encoding="utf-8", newline="") as f:
        f.write("\n".join(head) + "\n"); df.to_csv(f, index=False, float_format="%.6g")

# ============================== offset 与读长集合 ==============================
def peak_stats(x, L):
    x = x[1:L + 1]; n = int(x.sum())
    if n == 0: return None, 0, 0, 0.0, 0.0
    w = x[J_LO - 1:J_HI]; k = int(np.argmax(w)); srt = np.sort(w); pk = int(x[J_LO + k - 1])
    return J_LO + k - 1, n, pk, pk / n, (srt[-1] / srt[-2] if srt[-2] > 0 else np.inf)
OFF = {}; CONF = {}
for s in IDS:
    d = f"{WORK}/stats/{s}"
    hs = collections.defaultdict(lambda: np.zeros(600, np.int64)); hh = collections.defaultdict(lambda: np.zeros(600, np.int64))
    for r in csv.DictReader(open(f"{d}/host_startpos.tsv"), delimiter="\t"): hs[int(r["L"])][int(r["j"])] += int(r["n"])
    for r in csv.DictReader(open(f"{d}/host_startpos_split.tsv"), delimiter="\t"): hh[(int(r["half"]), int(r["L"]))][int(r["j"])] += int(r["n"])
    OFF[s], CONF[s] = {}, {}
    for L in sorted(hs):
        o, n, pk, share, margin = peak_stats(hs[L], L); o0 = peak_stats(hh[(0, L)], L)[0]; o1 = peak_stats(hh[(1, L)], L)[0]
        c = o is not None and n >= 500 and share >= 0.10 and margin >= 1.2 and o not in (J_LO - 1, J_HI - 1); st = o is not None and o0 == o and o1 == o
        OFF[s][L] = o; CONF[s][L] = "high" if (c and st) else "medium" if (c or st) else "low"
po = rdt(f"{RGP}/psite_offset_GIb.tsv"); chk = [int(OFF[r.sample_id].get(int(r.read_length)) == int(r.offset)) for _, r in po.iterrows() if r.sample_id in OFF and not pd.isna(r.offset)]
assert all(chk), "offset 重新推导与 psite_offset_GIb 不一致"
cs = pd.read_csv(f"{TASK}/read_length_frame_periodicity/read_length_frame_counts_and_fractions.csv", comment="#")
SR.assert_table_whitelisted(cs, "sample", "read_length_frame_counts"); SR.assert_exactly(set(cs["sample"]), SR.RIBO_VERO_ALIASES, "read_length_frame_counts")
body = cs[cs.region == "CDS_body_45_45"].copy()
body["hc"] = (body.offset_confidence == "high") & (body.total_host_CDS_psites >= 5000) & (body.F0 >= 0.5) & (body.top_to_second_frame_ratio >= 1.5)
SETDEF = {}
for s in IDS:
    b = body[(body["sample"] == s) & body.read_length.between(25, 35)]
    SETDEF[(s, "host_criteria")] = sorted(b[b.hc].read_length.astype(int))
    SETDEF[(s, "strict_core")] = sorted(b[b.length_set_category == "strict_core"].read_length.astype(int))
    SETDEF[(s, "29-31 nt")] = [29, 30, 31]
    SETDEF[(s, "all accepted (every length with an offset)")] = None
SETNAMES = ["host_criteria", "strict_core", "29-31 nt", "all accepted (every length with an offset)"]

# ============================== P-site 表（+ 链）：非剪接 / 剪接分开 ==============================
# 非剪接 reads：P-site = 5′ 端 + offset；剪接 reads（5′ 端 22 nt 内有 leader–body 跳跃）：沿比对碱基取第 offset+1 个，跨 junction 的 P-site 归属本质歧义，
# 因此主结果 ORF 级读框只用非剪接 reads（与冻结的 frame5 口径一致），剪接 reads 的 P-site 单列并做 sensitivity。
psite = {s: collections.Counter() for s in IDS}      # 非剪接 (L, pos) -> n
psite_sp = {s: collections.Counter() for s in IDS}   # 剪接
for s in IDS:
    with gzip.open(f"{WORK}/stats/{s}/pedv_5p_matrix.tsv.gz", "rt") as f:
        next(f)
        for l in f:
            st, L, p5, sc, n = l.rstrip("\n").split("\t")
            if st != "+": continue
            L = int(L); o = OFF[s].get(min(L, 150))
            if o is None: continue
            psite[s][(L, int(p5) + o)] += int(n)
    with gzip.open(f"{WORK}/stats/{s}/pedv_spliced5p.tsv.gz", "rt") as f:
        next(f)
        for l in f:
            st, L, sc, coords, n = l.rstrip("\n").split("\t")
            if st != "+": continue
            L = int(L); o = OFF[s].get(min(L, 150)); cs_ = coords.split(",")
            if o is None or o >= len(cs_): continue
            psite_sp[s][(L, int(cs_[o]))] += int(n)
# 对账：本脚本重建的 + 链 P-site 总数 == 已交付 psite_counts_PEDV_GIb（五个库）
big = rdt(f"{RGP}/psite_counts_PEDV_GIb.tsv.gz")
for s in IDS:
    mine = sum(psite[s].values()) + sum(psite_sp[s].values())
    if True:   # all six whitelisted libraries are in the delivered table
        assert mine == int(big[(big.sample_id == s) & (big.strand == "+")].psite_count.sum()), (s, mine)
print("rebuilt + strand P-site totals equal psite_counts_PEDV_GIb for all six whitelisted libraries")

# ============================== ORF 唯一区 ==============================
region = np.zeros(GLEN + 2, np.int8)       # 0 = 不属于任何 ORF 唯一区
cover = np.zeros(GLEN + 2, np.int8)
for _, a, b in ORFS: cover[a:b + 1] += 1
for k, (nm, a, b) in enumerate(ORFS, start=1):
    for p in range(a + TRIM, b - TRIM + 1):
        if cover[p] == 1: region[p] = k
n_unique = {nm: int((region == k).sum()) for k, (nm, a, b) in enumerate(ORFS, start=1)}

def frame_counts(s, Ls, spliced=False):
    """返回 {orf: np.array(3)}：frame = (P-site − ORF 起点) mod 3；spliced=True 时只统计剪接 reads 的 P-site"""
    out = {nm: np.zeros(3, np.int64) for nm, _, _ in ORFS}
    for (L, p), n in (psite_sp[s] if spliced else psite[s]).items():
        if Ls is not None and L not in Ls: continue
        k = region[p] if 0 <= p <= GLEN else 0
        if k:
            nm, a, b = ORFS[k - 1]; out[nm][(p - a) % 3] += n
    return out
# ---- 主口径：frame5（5′ 端落在 ORF 唯一区、按 offset 换算 P-site 帧；与冻结的 periodicity_virus.tsv 同口径）----
F5CNT = {}
for s in IDS:
    F = collections.defaultdict(lambda: np.zeros(3, np.int64))
    for r in csv.DictReader(open(f"{WORK}/stats/{s}/frame5.tsv"), delimiter="\t"): F[(r["domain"], int(r["L"]))] += [int(r["f0"]), int(r["f1"]), int(r["f2"])]
    F5CNT[s] = F
def frame_counts5(s, Ls):
    out = {}
    lens = sorted({L for (d, L) in F5CNT[s] if d == "PEDV_all"}) if Ls is None else Ls
    for nm, a, b in ORFS:
        p = np.zeros(3, np.int64)
        for L in lens:
            o = OFF[s].get(L)
            if o is None: continue
            f = F5CNT[s][(f"PEDV_{nm}", L)]; p += np.array([f[(k - o) % 3] for k in range(3)])
        out[nm] = p
    return out
rowsC = []
FR = {}
for s in IDS:
    for setn in SETNAMES:
        Ls = SETDEF[(s, setn)]
        fc = frame_counts5(s, None if Ls is None else list(Ls)); FR[(s, setn)] = fc
        fs = frame_counts(s, None if Ls is None else set(Ls), spliced=True)
        tot = np.zeros(3, np.int64)
        for nm, a, b in ORFS + [("all_7ORF", 0, 0)]:
            if nm != "all_7ORF": c = fc[nm]; tot = tot + c; nsp = int(fs[nm].sum())
            else: c = tot; nsp = int(sum(v.sum() for v in fs.values()))
            n = int(c.sum())
            f = c / n if n else np.array([np.nan] * 3)
            status = "valid" if n >= LOW else ("low_depth" if n >= MISSING else "missing")
            rowsC.append(dict(sample=s, treatment=DRUG[s], time_hpi=TIME[s], read_set=setn, lengths=("all" if Ls is None else ";".join(map(str, Ls))), orf=nm,
                              unique_region_nt=(sum(n_unique.values()) if nm == "all_7ORF" else n_unique[nm]), usable_psites=n, N0=int(c[0]), N1=int(c[1]), N2=int(c[2]),
                              F0=f[0], F1=f[1], F2=f[2], spliced_psites_in_region_not_used=nsp, frame0_minus_one_third=f[0] - 1 / 3,
                              ribosome_fraction_f=(f[0] - 1 / 3) / (BASE_F0 - 1 / 3) if n else np.nan, status=status))
C = pd.DataFrame(rowsC)
assert ((C.N0 + C.N1 + C.N2) == C.usable_psites).all()
wr("panel_C_ORF_frame_enrichment.csv", C, "Panel C: ORF-level reading-frame distribution of PEDV (+) footprints (ORF-unique regions) under four read-length sets",
   ["counting unit: reads whose 5' end falls in the ORF-unique region [start+45, end-45] minus any nt covered by another annotated ORF (ORF1a/ORF1b 31 nt, ORF1b/S 4 nt, S/ORF3 1 nt excluded); frame of the P-site = (5' end + length-specific offset - ORF start) mod 3",
    "this is the frozen RGP_v1.3 periodicity_virus.tsv definition (verified against the delivered table); reads with a leader-body junction within 22 nt of the 5' end are not counted (counts of their P-sites in the region are in spliced_psites_in_region_not_used)",
    "ribosome_fraction_f = (F0 - 1/3) / (0.80 - 1/3) (two-component model; may be negative or >1); status: valid >=1,000 usable reads, low_depth 100-999, missing <100",
    "ORF1b is counted in its own frame (-1 relative to ORF1a); usable_psites = number of reads (one P-site per read)"])
Bdf = C[(C.orf == "all_7ORF")][["sample", "treatment", "time_hpi", "read_set", "lengths", "usable_psites", "N0", "N1", "N2", "F0", "F1", "F2", "ribosome_fraction_f", "status"]]
wr("panel_B_viral_frame_periodicity.csv", Bdf, "Panel B: reading-frame distribution of PEDV (+) footprints pooled over the seven ORF-unique regions")
# 与 delivered periodicity_virus.tsv 对账（strict_core，五个库）
pv = rdt(f"{RGP}/periodicity_virus.tsv"); pv = pv[pv.read_length == "pooled_strict_core"]
dd = []
for s in IDS:
    if s == SR.lib("RIBO_VERO_CHX_1MOI_16h"): continue
    for nm, a, b in ORFS:
        mine = C[(C["sample"] == s) & (C.read_set == "strict_core") & (C.orf == nm)].F0.iloc[0] * 100
        deliv = float(pv[(pv.sample_id == s) & (pv.domain == f"PEDV_{nm}")].phase0_pct.iloc[0]); dd.append(abs(mine - deliv))
assert max(dd) < 0.06, max(dd)
print("frame5-based ORF frame0 equals delivered periodicity_virus.tsv (max |diff| = %.3f pp)" % max(dd))

# ---- sensitivity：真实 P-site 位置口径（非剪接 reads，P-site 落在 ORF 唯一区）----
srows = []
for s in IDS:
    for setn in ("host_criteria", "strict_core"):
        Ls = SETDEF[(s, setn)]; fp = frame_counts(s, set(Ls))
        for nm, a, b in ORFS:
            c = fp[nm]; n = int(c.sum()); m5 = C[(C["sample"] == s) & (C.read_set == setn) & (C.orf == nm)].iloc[0]
            top = 0.0
            srows.append(dict(sample=s, read_set=setn, orf=nm, frame5_n=int(m5.usable_psites), frame5_F0=m5.F0, psite_position_n=n, psite_position_F0=c[0] / n if n else np.nan))
M = pd.DataFrame(srows); M["F0_diff_psite_position_minus_frame5"] = M.psite_position_F0 - M.frame5_F0
M.to_csv(f"{SD}/sensitivity_psite_position_vs_frame5.csv", index=False, float_format="%.5g")
print("psite-position vs frame5: median |dF0| = %.4f; rows |dF0|>0.05: %d of %d" % (M.F0_diff_psite_position_minus_frame5.abs().median(), int((M.F0_diff_psite_position_minus_frame5.abs() > 0.05).sum()), len(M)))
print(M.reindex(M.F0_diff_psite_position_minus_frame5.abs().sort_values(ascending=False).index).head(5).round(3).to_string())


# ============================== Panel A：链向 ==============================
al = rdt(f"{RGP}/riboseq_GIb_alignment_summary.tsv").set_index("sample_id")
host_ann = {r["sample_id"]: r for r in csv.DictReader(open(f"{TASK}/riboseq_qc_figs_R1_R3/_work/host_annotation/host_annotation_all.tsv"), delimiter="\t") if r["sample_id"] in SR.WHITELIST_LIBRARIES}   # final sample table only
host_ann[SR.lib("RIBO_VERO_CHX_1MOI_16h")] = list(csv.DictReader(open(f"{TASK}/riboseq_qc_figs_R1_R3/_work/host_annotation/{SR.lib('RIBO_VERO_CHX_1MOI_16h')}.host_annotation.tsv"), delimiter="\t"))[0]
rowsA = []; rowsL = []
for s in IDS:
    d = f"{WORK}/stats/{s}"
    am = {r["domain"]: r for r in csv.DictReader(open(f"{d}/aln_metrics.tsv"), delimiter="\t")}
    dup = {r["domain"]: r for r in csv.DictReader(open(f"{d}/dup.tsv"), delimiter="\t")}
    rl = pd.read_csv(f"{d}/readlen_aln.tsv", sep="\t")
    pl, mi, host = int(am["pedv_plus"]["reads"]), int(am["pedv_minus"]["reads"]), int(am["host"]["reads"])
    if s in al.index:
        assert int(al.loc[s, "PEDV_plus_reads"]) == pl and int(al.loc[s, "PEDV_minus_reads"]) == mi
    assert int(am["pedv_all"]["reads"]) == pl + mi
    ha = {c: int(host_ann[s][c]) for c in ("CDS", "5UTR", "3UTR", "other_exon", "intron", "antisense", "intergenic")}
    L = rl.L.values
    def stat(col):
        v = rl[col].values.astype(float); cum = np.cumsum(v) / v.sum()
        return int(L[np.searchsorted(cum, 0.5)]), v[(L >= 29) & (L <= 31)].sum() / v.sum()
    mp, fp = stat("pedv_plus_n"); mm, fm = stat("pedv_minus_n")
    # − 链 P-site 在 ORF 唯一区内的最大帧占比（方向无关）
    cnt = collections.Counter()
    with gzip.open(f"{WORK}/stats/{s}/pedv_5p_matrix.tsv.gz", "rt") as f:
        next(f)
        for l in f:
            st, Lr, p5, sc, n = l.rstrip("\n").split("\t")
            if st != "-": continue
            Lr = int(Lr)
            if Lr not in SETDEF[(s, "host_criteria")]: continue
            o = OFF[s].get(Lr)
            if o is None: continue
            p = int(p5) - o
            if 0 <= p <= GLEN and region[p]: cnt[p % 3] += int(n)
    tm = sum(cnt.values()); dom = max(cnt.values()) / tm if tm else np.nan
    rowsA.append(dict(sample=s, treatment=DRUG[s], time_hpi=TIME[s], plus_reads=pl, minus_reads=mi, strand_assignable_reads=pl + mi, plus_fraction=pl / (pl + mi),
                      minus_fraction=mi / (pl + mi), plus_over_minus=pl / mi, unassignable_reads=0,
                      PEDV_pct_of_unique_mapped=100 * (pl + mi) / (pl + mi + host), host_antisense_fraction_of_host_psites=ha["antisense"] / sum(ha.values()),
                      median_length_plus=mp, median_length_minus=mm, fraction_29_31_plus=fp, fraction_29_31_minus=fm,
                      softclip_read_frac_minus=int(am["pedv_minus"]["softclip_reads"]) / mi, mismatch_per100_minus=100 * int(am["pedv_minus"]["nM_sum"]) / int(am["pedv_minus"]["aligned_M_bases"]),
                      dup_frac_minus=1 - int(dup["pedv_minus"]["distinct_5p_keys"]) / mi, minus_psites_in_ORF_unique_regions_host_criteria=tm,
                      minus_dominant_frame_fraction=dom, minus_coverage=("insufficient" if tm < LOW else "sufficient")))
    for dom_, col in (("PEDV (+)", "pedv_plus_n"), ("PEDV (-)", "pedv_minus_n")):
        v = rl[col].values.astype(float); tot = v.sum()
        for Lr, x in zip(L, v):
            if 16 <= Lr <= 60: rowsL.append(dict(sample=s, treatment=DRUG[s], time_hpi=TIME[s], strand=dom_, aligned_length=int(Lr), reads=int(x), fraction_of_strand_reads=x / tot))
A = pd.DataFrame(rowsA)
wr("panel_A_strand_counts.csv", A, "Panel A: strand assignment of PEDV uniquely mapped reads (STAR, NH=1)",
   ["strand = alignment orientation; the Ribo-seq libraries retain strand (RA3 adapter ligation) but no library-prep record was found, so strand-specificity is supported only by the host antisense baseline",
    "plus + minus = PEDV_unique; unassignable_reads = 0 by construction (every alignment has a strand)",
    "host_antisense_fraction_of_host_psites: P-sites falling on the antisense strand of an annotated gene (unstranded-agnostic: includes real overlapping opposite-strand genes)"])
wr("panel_A_strand_length_profiles.csv", pd.DataFrame(rowsL), "Panel A: aligned-length profile of PEDV (+) and (-) reads (all uniquely mapped; denominator = reads of the same strand)")

# ============================== Panel D：长度分布 ==============================
rowsD = []; rowsQ = []
for s in IDS:
    d = f"{WORK}/stats/{s}"
    rl = pd.read_csv(f"{d}/readlen_aln.tsv", sep="\t").set_index("L")
    F = collections.defaultdict(lambda: np.zeros(3, np.int64))
    for r in csv.DictReader(open(f"{d}/frame5.tsv"), delimiter="\t"): F[(r["domain"], int(r["L"]))] += [int(r["f0"]), int(r["f1"]), int(r["f2"])]
    prof = {"PEDV (+) all uniquely mapped": {L: int(rl.loc[L, "pedv_plus_n"]) for L in rl.index},
            "PEDV (+) 5' end in ORF-unique region": {L: int(F[("PEDV_all", L)].sum()) for L in rl.index},
            "Host CDS-body (5' end)": {L: int(F[("host", L)].sum()) for L in rl.index}}
    den = {"PEDV (+) all uniquely mapped": "PEDV (+) uniquely mapped reads", "PEDV (+) 5' end in ORF-unique region": "PEDV footprints with 5' end in ORF-unique regions",
           "Host CDS-body (5' end)": "host CDS-body footprints"}
    for nm, dct in prof.items():
        tot = sum(dct.values()); Ls = sorted(dct)
        v = np.array([dct[L] for L in Ls], float); cum = np.cumsum(v) / tot
        med = int(Ls[np.searchsorted(cum, 0.5)])
        f2931 = v[[29 <= L <= 31 for L in Ls]].sum() / tot
        short = v[[L <= 27 for L in Ls]].sum() / tot; long_ = v[[L >= 32 for L in Ls]].sum() / tot
        for L in Ls:
            if 16 <= L <= 60: rowsD.append(dict(sample=s, treatment=DRUG[s], time_hpi=TIME[s], profile=nm, denominator=den[nm], aligned_length=L, reads=dct[L], fraction=dct[L] / tot))
        rowsQ.append(dict(sample=s, treatment=DRUG[s], time_hpi=TIME[s], profile=nm, denominator=den[nm], n_reads=tot, median_length=med, fraction_29_31=f2931,
                          fraction_le27=short, fraction_ge32=long_))
D = pd.DataFrame(rowsD); Q = pd.DataFrame(rowsQ)
# 16 hpi 专项：加 frame0（host_criteria，ORF 唯一区，P-site 位置口径）
fr0 = {s: C[(C["sample"] == s) & (C.read_set == "host_criteria") & (C.orf == "all_7ORF")].F0.iloc[0] for s in IDS}
Q["frame0_viral_ORF_unique_host_criteria"] = Q["sample"].map(fr0).where(Q.profile == "PEDV (+) 5' end in ORF-unique region")
bd = body[body.read_length.isin(range(25, 36))]
hf0 = {}
for s in IDS:
    Ls = SETDEF[(s, "host_criteria")]; b = bd[(bd["sample"] == s) & bd.read_length.isin(Ls)]
    hf0[s] = b.N0.sum() / (b.N0 + b.N1 + b.N2).sum()
Q["frame0_host_CDS_body_host_criteria"] = Q["sample"].map(hf0).where(Q.profile == "Host CDS-body (5' end)")
wr("panel_D_host_PEDV_length_profiles.csv", D, "Panel D: aligned-length profiles; denominators differ by profile (see 'denominator')",
   ["no smoothing; PEDV (+) 5' end in ORF-unique region = reads of frame5.tsv domain PEDV_all (7 canonical ORF unique interiors, junction-clean)", "host CDS-body (5' end) = frame5.tsv domain host (canonical protein-coding CDS minus 45 nt at both ends)"])
wr("panel_D_length_QC_summary.csv", Q, "Panel D: summary of footprint-length features per library and profile",
   ["fraction_le27 / fraction_ge32 = shorter / longer fragments; frame0 columns use the host_criteria length set of each library"])

# ============================== 整合 QC 表 ==============================
rows = []
Bm = C[(C.read_set == "host_criteria") & (C.orf == "all_7ORF")].set_index("sample")
Bs = C[(C.read_set == "strict_core") & (C.orf == "all_7ORF")].set_index("sample")
for s in IDS:
    a = A[A["sample"] == s].iloc[0]
    q = Q[(Q["sample"] == s)].set_index("profile")
    n_eval = int(((C["sample"] == s) & (C.read_set == "host_criteria") & (C.orf != "all_7ORF") & (C.status == "valid")).sum())
    rows.append(dict(sample=s, treatment=DRUG[s], time_hpi=TIME[s], PEDV_uniquely_mapped=int(a.strand_assignable_reads), strand_assignable_reads=int(a.strand_assignable_reads),
                     positive_strand_fraction=a.plus_fraction, viral_CDS_usable_psites_host_criteria=int(Bm.loc[s, "usable_psites"]), viral_frame0_host_criteria=Bm.loc[s, "F0"],
                     viral_ribosome_fraction_f_host_criteria=Bm.loc[s, "ribosome_fraction_f"], viral_frame0_strict_core=Bs.loc[s, "F0"],
                     host_CDS_frame0_host_criteria=hf0[s], viral_median_length_ORF_footprints=int(q.loc["PEDV (+) 5' end in ORF-unique region", "median_length"]),
                     host_median_length_CDS_footprints=int(q.loc["Host CDS-body (5' end)", "median_length"]),
                     viral_fraction_29_31_ORF_footprints=q.loc["PEDV (+) 5' end in ORF-unique region", "fraction_29_31"],
                     host_fraction_29_31_CDS_footprints=q.loc["Host CDS-body (5' end)", "fraction_29_31"], evaluable_ORFs_of_7=n_eval,
                     PEDV_pct_of_unique_mapped=a.PEDV_pct_of_unique_mapped, minus_fraction=a.minus_fraction))
I = pd.DataFrame(rows)
wr("integrated_viral_translation_QC.csv", I, "Integrated viral translation-feature QC per library",
   ["evaluable ORF = ORF with >=1,000 usable P-sites under the host_criteria length set", f"permutation-null false-pass rate (ribotricer) is NOT available for {SR.lib('RIBO_VERO_CHX_1MOI_16h')} (not run); quality tier is by viral frame0 only"])
from scipy.stats import spearmanr
r1 = spearmanr(I.PEDV_pct_of_unique_mapped, I.viral_frame0_host_criteria); r2 = spearmanr(I.PEDV_pct_of_unique_mapped, I.viral_ribosome_fraction_f_host_criteria)
r3 = spearmanr(I.PEDV_pct_of_unique_mapped, I.viral_frame0_strict_core)
print("Spearman PEDV%% unique vs frame0 host_criteria: %.2f (p=%.2f); vs f: %.2f (p=%.2f); vs frame0 strict_core: %.2f (p=%.2f)" % (r1.correlation, r1.pvalue, r2.correlation, r2.pvalue, r3.correlation, r3.pvalue))
json.dump({"spearman_pedv_pct_vs_frame0_hc": [r1.correlation, r1.pvalue], "spearman_pedv_pct_vs_f_hc": [r2.correlation, r2.pvalue], "spearman_pedv_pct_vs_frame0_sc": [r3.correlation, r3.pvalue],
           "sets": {f"{s}|{k}": v for (s, k), v in SETDEF.items() if v is not None}}, open(f"{OUT}/_numbers.json", "w"), indent=1)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
print(I.round(3).to_string())
print(A[["sample", "minus_fraction", "host_antisense_fraction_of_host_psites", "minus_dominant_frame_fraction", "minus_psites_in_ORF_unique_regions_host_criteria"]].round(4).to_string())
print(Q[Q.profile.str.startswith("PEDV (+) 5'")][["sample", "n_reads", "median_length", "fraction_29_31", "fraction_le27", "fraction_ge32", "frame0_viral_ORF_unique_host_criteria"]].round(3).to_string())
print(Q[Q.profile.str.startswith("Host")][["sample", "n_reads", "median_length", "fraction_29_31", "fraction_le27", "fraction_ge32", "frame0_host_CDS_body_host_criteria"]].round(3).to_string())
print(Q[Q.profile.str.startswith("PEDV (+) all")][["sample", "n_reads", "median_length", "fraction_29_31", "fraction_le27", "fraction_ge32"]].round(3).to_string())
