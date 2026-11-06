#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Sup Fig 13 Phase 1 可行性审计（纯本地；只读 RGP_v1.3 / stats / psite 表）。
CHX 16 h = RIBO_VERO_CHX_1MOI_16h（最终样本表中的 CHX 16 h 文库）。
输出 phase1/tables/*.csv，供 PHASE1_AUDIT_REPORT.md 引用。"""
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
import os, gzip, csv, json, collections, math
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.dirname(HERE); TASK = os.path.dirname(OUT)
RGP = f"{TASK}/riboseq_GIb_mapping_psite"; WORK = f"{RGP}/_work/server"
TAB = f"{OUT}/phase1_tables"; os.makedirs(TAB, exist_ok=True)
SAMPLES = [(SR.lib("RIBO_VERO_CHX_1MOI_5h"), "CHX", 5), (SR.lib("RIBO_VERO_CHX_1MOI_8h"), "CHX", 8), (SR.lib("RIBO_VERO_CHX_1MOI_16h"), "CHX", 16),
           (SR.lib("RIBO_VERO_HAR_1MOI_5h"), "HAR", 5), (SR.lib("RIBO_VERO_HAR_1MOI_8h"), "HAR", 8), (SR.lib("RIBO_VERO_HAR_1MOI_16h"), "HAR", 16)]
IDS = [s for s, _, _ in SAMPLES]; LAB = {s: f"{d} {t}h" for s, d, t in SAMPLES}
ORFS = [("ORF1a", 293, 12622), ("ORF1b", 12592, 20613), ("S", 20610, 24758), ("ORF3", 24758, 25033), ("E", 25361, 25600), ("M", 25608, 26288), ("N", 26300, 27625)]
J_LO, J_HI = 7, 19
rdt = lambda p: pd.read_csv(p, sep="\t", comment="#", low_memory=False)

def peak_stats(x, L):
    x = x[1:L + 1]; n = int(x.sum())
    if n == 0: return None, 0, 0, 0.0, 0.0
    w = x[J_LO - 1:J_HI]; k = int(np.argmax(w)); srt = np.sort(w); pk = int(x[J_LO + k - 1])
    return J_LO + k - 1, n, pk, pk / n, (srt[-1] / srt[-2] if srt[-2] > 0 else np.inf)
def derive_offsets(sid):
    d = f"{WORK}/stats/{sid}"
    hs = collections.defaultdict(lambda: np.zeros(600, np.int64)); hh = collections.defaultdict(lambda: np.zeros(600, np.int64))
    for r in csv.DictReader(open(f"{d}/host_startpos.tsv"), delimiter="\t"): hs[int(r["L"])][int(r["j"])] += int(r["n"])
    for r in csv.DictReader(open(f"{d}/host_startpos_split.tsv"), delimiter="\t"): hh[(int(r["half"]), int(r["L"]))][int(r["j"])] += int(r["n"])
    off, conf = {}, {}
    for L in sorted({L for L in hs}):
        o, n, pk, share, margin = peak_stats(hs[L], L); o0 = peak_stats(hh[(0, L)], L)[0]; o1 = peak_stats(hh[(1, L)], L)[0]
        c = o is not None and n >= 500 and share >= 0.10 and margin >= 1.2 and o not in (J_LO - 1, J_HI - 1)
        st = o is not None and o0 == o and o1 == o
        off[L] = o; conf[L] = ("high" if (c and st) else "medium" if (c or st) else "low")
    return off, conf
OFF, CONF = {}, {}
for s in IDS: OFF[s], CONF[s] = derive_offsets(s)
# 与已交付 psite_offset_GIb 比对（五个库）
po = rdt(f"{RGP}/psite_offset_GIb.tsv"); mism = 0; tot = 0
for _, r in po.iterrows():
    if r.sample_id in OFF and not pd.isna(r.offset):
        tot += 1; mism += int(OFF[r.sample_id].get(int(r.read_length)) != int(r.offset))
print(f"offset re-derivation vs psite_offset_GIb: {tot - mism}/{tot} identical")

# ---------- 逐 ORF 周期性（frame5 口径，与 periodicity_virus.tsv 同）----------
def frame5(sid):
    F = collections.defaultdict(lambda: np.zeros(3, np.int64))
    for r in csv.DictReader(open(f"{WORK}/stats/{sid}/frame5.tsv"), delimiter="\t"): F[(r["domain"], int(r["L"]))] += [int(r["f0"]), int(r["f1"]), int(r["f2"])]
    return F
pf = lambda f, off: np.array([f[(k - off) % 3] for k in range(3)])
# 读长集合：strict_core（冻结）、host_criteria（新）、all、29-31
cs = pd.read_csv(f"{TASK}/read_length_frame_periodicity/read_length_frame_counts_and_fractions.csv", comment="#")
SR.assert_table_whitelisted(cs, "sample", "read_length_frame_counts"); SR.assert_exactly(set(cs["sample"]), SR.RIBO_VERO_ALIASES, "read_length_frame_counts")
body = cs[cs.region == "CDS_body_45_45"].copy()
body["hc"] = (body.offset_confidence == "high") & (body.total_host_CDS_psites >= 5000) & (body.F0 >= 0.5) & (body.top_to_second_frame_ratio >= 1.5)
SETS = {}
for s in IDS:
    b = body[(body["sample"] == s) & body.read_length.between(25, 35)]
    SETS[(s, "host_criteria")] = sorted(b[b.hc].read_length.astype(int))
    SETS[(s, "strict_core")] = sorted(b[b.length_set_category == "strict_core"].read_length.astype(int))
    SETS[(s, "29-31nt")] = [29, 30, 31]
    SETS[(s, "all_25-35")] = list(range(25, 36))
rows = []
for s in IDS:
    F = frame5(s)
    for setname in ("strict_core", "host_criteria", "29-31nt", "all_25-35"):
        Ls = SETS[(s, setname)]
        for nm, a, b in ORFS + [("all_7ORF", 0, 0)]:
            dom = "PEDV_all" if nm == "all_7ORF" else f"PEDV_{nm}"
            p = np.zeros(3, np.int64)
            for L in Ls:
                o = OFF[s].get(L)
                if o is None: continue
                p += pf(F[(dom, L)], o)
            n = int(p.sum())
            f0 = p[0] / n if n else np.nan
            rows.append(dict(sample=s, label=LAB[s], read_set=setname, lengths=";".join(map(str, Ls)), orf=nm, n_psites=n, N0=int(p[0]), N1=int(p[1]), N2=int(p[2]),
                             frame0=f0, f_ribosome=(f0 - 1 / 3) / (0.80 - 1 / 3) if n else np.nan))
P = pd.DataFrame(rows); P.to_csv(f"{TAB}/orf_frame0_by_set.csv", index=False, float_format="%.5g")
# 复核任务书给出的 strict_core 表（五个库）
told = {"ORF1a": [38.9, 40.3, 67.6, 39.7, 35.9, 48.6], "ORF1b": [33.0, 32.7, 53.2, 34.6, 44.0, 37.4], "S": [64.3, 55.9, 78.2, 59.8, 48.0, 82.9],
        "ORF3": [25.2, 14.4, 77.9, 18.6, 13.8, 54.4], "E": [74.3, 57.6, 76.1, 60.1, 46.7, 89.8], "M": [70.2, 77.0, 79.3, 59.2, 48.2, 83.0],
        "N": [80.5, 77.0, 83.9, 68.8, 65.5, 86.0], "all_7ORF": [67.8, 50.3, 82.2, 56.0, 55.0, 81.3]}
told_order = [SR.lib("RIBO_VERO_CHX_1MOI_5h"), SR.lib("RIBO_VERO_HAR_1MOI_5h"), SR.lib("RIBO_VERO_CHX_1MOI_8h"), SR.lib("RIBO_VERO_HAR_1MOI_8h"), SR.lib("RIBO_VERO_CHX_1MOI_16h"), SR.lib("RIBO_VERO_HAR_1MOI_16h")]   # CHX5 HAR5 CHX8 HAR8 CHX16 HAR16
chk = []
pv = rdt(f"{RGP}/periodicity_virus.tsv"); pv = pv[pv.read_length == "pooled_strict_core"]
for orf, vals in told.items():
    for s, v in zip(told_order, vals):
        mine = P[(P["sample"] == s) & (P.read_set == "strict_core") & (P.orf == orf)].frame0.iloc[0] * 100
        deliv = pv[(pv.sample_id == s) & (pv.domain == ("PEDV_all_7ORF" if orf == "all_7ORF" else f"PEDV_{orf}"))].phase0_pct
        chk.append(dict(orf=orf, sample=s, task_text=v, delivered_table=float(deliv.iloc[0]) if len(deliv) else np.nan, recomputed=round(mine, 1)))
C = pd.DataFrame(chk); C.to_csv(f"{TAB}/periodicity_table_recheck.csv", index=False)
print("task text vs delivered periodicity_virus: max |diff| =", (C.task_text - C.delivered_table).abs().max(), "; recomputed vs delivered:", (C.recomputed - C.delivered_table).abs().max())

# ---------- 链向 ----------
al = rdt(f"{RGP}/riboseq_GIb_alignment_summary.tsv").set_index("sample_id")
strand = []
for s in IDS:
    d = f"{WORK}/stats/{s}"
    am = {r["domain"]: r for r in csv.DictReader(open(f"{d}/aln_metrics.tsv"), delimiter="\t")}
    dup = {r["domain"]: r for r in csv.DictReader(open(f"{d}/dup.tsv"), delimiter="\t")}
    rl = pd.read_csv(f"{d}/readlen_aln.tsv", sep="\t")
    pl, mi = int(am["pedv_plus"]["reads"]), int(am["pedv_minus"]["reads"])
    host = int(am["host"]["reads"])
    L = rl.L.values
    def f2931(col): v = rl[col].values.astype(float); return v[(L >= 29) & (L <= 31)].sum() / v.sum() if v.sum() else np.nan
    def med(col):
        v = rl[col].values.astype(float); cum = np.cumsum(v) / v.sum(); return int(L[np.searchsorted(cum, 0.5)])
    def pct(dom, num, den): return int(am[dom][num]) / max(int(am[dom][den]), 1)
    strand.append(dict(sample=s, label=LAB[s], plus_reads=pl, minus_reads=mi, minus_fraction=mi / (pl + mi), plus_over_minus=pl / mi,
                       PEDV_pct_of_unique=100 * (pl + mi) / (pl + mi + host),
                       median_len_plus=med("pedv_plus_n"), median_len_minus=med("pedv_minus_n"),
                       f2931_plus=f2931("pedv_plus_n"), f2931_minus=f2931("pedv_minus_n"),
                       softclip_frac_minus=int(am["pedv_minus"]["softclip_reads"]) / mi, softclip_frac_plus=int(am["pedv_plus"]["softclip_reads"]) / pl,
                       mismatch_per100_minus=100 * int(am["pedv_minus"]["nM_sum"]) / max(int(am["pedv_minus"]["aligned_M_bases"]), 1),
                       mismatch_per100_plus=100 * int(am["pedv_plus"]["nM_sum"]) / max(int(am["pedv_plus"]["aligned_M_bases"]), 1),
                       dup_frac_minus=1 - int(dup["pedv_minus"]["distinct_5p_keys"]) / mi, dup_frac_plus=1 - int(dup["pedv_plus"]["distinct_5p_keys"]) / pl,
                       spliced_frac_minus=int(am["pedv_minus"]["spliced_reads"]) / mi, spliced_frac_plus=int(am["pedv_plus"]["spliced_reads"]) / pl))
S = pd.DataFrame(strand); S.to_csv(f"{TAB}/strand_overview.csv", index=False, float_format="%.5g")
# 与 alignment_summary 对账
for s in IDS:
    if s in al.index:
        assert int(al.loc[s, "PEDV_plus_reads"]) == int(S[S["sample"] == s].plus_reads.iloc[0]) and int(al.loc[s, "PEDV_minus_reads"]) == int(S[S["sample"] == s].minus_reads.iloc[0])
print("strand counts match riboseq_GIb_alignment_summary for the 5 libraries in it")

# − 链读框：用 5′ 端矩阵 + offset，在 ORF 唯一区内统计 (P-site 坐标 mod 3) 的最大帧占比（方向无关，与 + 链同一口径比较）
def trimmed_regions():
    regs = []
    for nm, a, b in ORFS:
        others = [(a2, b2) for n2, a2, b2 in ORFS if n2 != nm]
        mask = np.zeros(b - a + 1, bool)
        for a2, b2 in others:
            lo, hi = max(a, a2), min(b, b2)
            if lo <= hi: mask[lo - a:hi - a + 1] = True
        for p in range(a + 45, b - 45 + 1):
            if not mask[p - a]: regs.append((p, nm))
    return {p: nm for p, nm in regs}
REG = trimmed_regions()
minus_frames = []
for s in IDS:
    for st in ("+", "-"):
        cnt = collections.Counter()
        with gzip.open(f"{WORK}/stats/{s}/pedv_5p_matrix.tsv.gz", "rt") as f:
            next(f)
            for l in f:
                sg, L, p5, sc, n = l.rstrip("\n").split("\t")
                if sg != st: continue
                L = int(L)
                if L not in SETS[(s, "host_criteria")]: continue
                o = OFF[s].get(L)
                if o is None: continue
                p = int(p5) + o if st == "+" else int(p5) - o
                if p in REG: cnt[p % 3] += int(n)
        t = sum(cnt.values()); fr = [cnt[i] / t if t else np.nan for i in range(3)]
        minus_frames.append(dict(sample=s, label=LAB[s], strand=st, n_psites_in_ORF_unique_regions=t, frame_mod3_0=fr[0], frame_mod3_1=fr[1], frame_mod3_2=fr[2],
                                 dominant_frame_fraction=max(fr) if t else np.nan))
MF = pd.DataFrame(minus_frames); MF.to_csv(f"{TAB}/minus_vs_plus_frame_dominance.csv", index=False, float_format="%.5g")

# 共变（n=6，Spearman）
from scipy.stats import spearmanr
cov = []
for col in ("softclip_frac_minus", "mismatch_per100_minus", "dup_frac_minus", "spliced_frac_minus", "PEDV_pct_of_unique", "median_len_minus", "f2931_minus"):
    r = spearmanr(S.minus_fraction, S[col]); cov.append(dict(variable=col, spearman_rho=r.correlation, p=r.pvalue))
pd.DataFrame(cov).to_csv(f"{TAB}/minus_fraction_covariation_n6.csv", index=False, float_format="%.4g")

# ---------- 坐标核对 ----------
oc = pd.read_csv(f"{TASK}/refbuild_GIb/v1.1/virus/PEDV_GIb.orf_check.tsv", sep="\t")
oc[["orf", "GIb_start", "GIb_end"]].to_csv(f"{TAB}/orf_coordinates_refbuild.csv", index=False)
fa = "".join(l.strip() for l in open(f"{TASK}/refbuild_GIb/v1.1/virus/PEDV_GIb.fa") if not l.startswith(">"))
print("genome length", len(fa), "| 12586-12592:", fa[12585:12592], "| 12580-12600:", fa[12579:12600])
import hashlib; print("PEDV_GIb.fa sha256:", hashlib.sha256(open(f"{TASK}/refbuild_GIb/v1.1/virus/PEDV_GIb.fa", "rb").read()).hexdigest())

pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
print(oc[["orf", "GIb_start", "GIb_end"]].to_string())
print(P[(P.read_set == "strict_core")].pivot(index="orf", columns="label", values="frame0").round(3)[[LAB[s] for s in IDS]].to_string())
print(P[(P.read_set == "host_criteria")].pivot(index="orf", columns="label", values="frame0").round(3)[[LAB[s] for s in IDS]].to_string())
print(P[(P.read_set == "host_criteria")].pivot(index="orf", columns="label", values="f_ribosome").round(2)[[LAB[s] for s in IDS]].to_string())
print(P[(P.read_set == "host_criteria")].pivot(index="orf", columns="label", values="n_psites")[[LAB[s] for s in IDS]].to_string())
print(S[["label", "plus_reads", "minus_reads", "minus_fraction", "PEDV_pct_of_unique", "median_len_minus", "f2931_plus", "f2931_minus", "softclip_frac_minus", "mismatch_per100_minus", "dup_frac_minus"]].round(4).to_string())
print(MF.round(3).to_string()); print(pd.DataFrame(cov).round(3).to_string())
print({k: v for k, v in SETS.items() if k[1] in ("strict_core", "host_criteria")})
