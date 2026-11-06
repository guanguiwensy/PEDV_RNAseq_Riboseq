#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""PK1_psite_v1 §2–§3：猪宿主 + GI-b 病毒上独立计算 P-site offset（Module 1 方法 B）、C1/C2/C3、四套读长集合。
改自 riboseq_GIb_mapping_psite/scripts/build_B_offset_sets.py（RGP_v1.3）。判定函数 peak_stats / pframe / fmetrics / c2_c3 / 集合规则逐行保持与原脚本一致；
PK1 改动：（1）不再读取 Module 1 的 Vero 表，改为 Vero 冻结值（29/30 nt = 12）三列对照「Vero 值 / PK1 值 / 是否一致」；（2）加 qc_status（pass / fail / insufficient_depth）；
（3）31–33 nt 若 offset 半数不一致或宿主 frame0 低于 0.40，标「frame分裂-待确认」；（4）输出每（样本, 读长）的 PEDV reads 数及每 ORF 深度；（5）输出名 *_PK1。
不删除任何 reads；集合只是分析口径。"""
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
import os, sys, collections, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

# ---------- 与 RGP_v1.3 / M1 完全相同的判定函数 ----------
def peak_stats(x, L):
    x = x[1:L + 1]; n = int(x.sum())
    if n == 0: return None, 0, 0, 0.0, 0.0
    w = x[J_LO - 1:J_HI]; k = int(np.argmax(w)); srt = np.sort(w)
    pk = int(x[J_LO + k - 1])
    return J_LO + k - 1, n, pk, pk / n, (srt[-1] / srt[-2] if srt[-2] > 0 else np.inf)
def pframe(f, off): return np.array([f[(k - off) % 3] for k in range(3)])
def fmetrics(p):
    t = p.sum()
    if t == 0: return None
    fr = p / t; s = np.sort(fr)
    return {"frac": fr, "dom": int(np.argmax(fr)), "domfrac": float(s[-1]),
            "ratio": float(s[-1] / s[-2]) if s[-2] > 0 else np.inf,
            "cfe": float(fr[0] / ((fr[1] + fr[2]) / 2)) if (fr[1] + fr[2]) > 0 else np.inf}
def c2_c3(f5_counts, off):
    n = int(f5_counts.sum())
    if n == 0: return False, False, 0, None
    s = np.sort(f5_counts / n); c2 = n >= 1000 and s[-1] >= 1.5 * s[-2]
    if off is None: return c2, False, n, None
    pf = pframe(f5_counts, off) / n; c3 = int(np.argmax(pf)) == 0 and pf[0] >= 0.40
    return c2, c3, n, pf
f3 = lambda x: f"{x:.3f}" if x is not None and np.isfinite(x) else ("inf" if x is not None else "NA")
f2 = lambda x: f"{x:.2f}" if x is not None and np.isfinite(x) else ("inf" if x is not None else "NA")
fr3 = lambda a: "/".join(f"{v:.3f}" for v in a) if a is not None else "NA"

def load(sid):
    d = f"{WORK}/stats/{sid}"
    hs = collections.defaultdict(lambda: np.zeros(600, np.int64))
    for r in rd(f"{d}/host_startpos.tsv"): hs[int(r["L"])][int(r["j"])] += int(r["n"])
    hh = collections.defaultdict(lambda: np.zeros(600, np.int64))
    for r in rd(f"{d}/host_startpos_split.tsv"): hh[(int(r["half"]), int(r["L"]))][int(r["j"])] += int(r["n"])
    rl = {int(r["L"]): {k: int(v) for k, v in r.items()} for r in rd(f"{d}/readlen_aln.tsv")}
    f5 = collections.defaultdict(lambda: np.zeros(3, np.int64))
    for r in rd(f"{d}/frame5.tsv"): f5[(r["domain"], int(r["L"]))] += [int(r["f0"]), int(r["f1"]), int(r["f2"])]
    return hs, hh, rl, f5

offs, src, sets_long, sets_wide, split31, vero_diff, depth_rows, orf_depth = [], [], [], [], [], [], [], []
OFFSET_TABLE = {}
for sid in IDS:
    m = mcols(sid)
    hs, hh, rl, f5 = load(sid)
    htot = sum(r["host_n"] for r in rl.values()); vtot = sum(r["pedv_plus_n"] for r in rl.values())
    info = {}; passed = {}; off = {}
    for L in sorted(rl):
        o, n, pk, share, margin = peak_stats(hs[L], L)
        o0 = peak_stats(hh[(0, L)], L)[0]; o1 = peak_stats(hh[(1, L)], L)[0]
        conf = o is not None and n >= 500 and share >= 0.10 and margin >= 1.2 and o not in (J_LO - 1, J_HI - 1)
        stab = o is not None and o0 == o and o1 == o
        nbo = [peak_stats(hs[x], x)[0] for x in (L - 1, L + 1)]
        nb = o is not None and all(v is None or abs(v - o) <= 1 for v in nbo)
        c1 = conf and stab
        level = ("NA" if o is None else "high" if c1 else "medium" if (conf or stab) else "low")
        hm = fmetrics(pframe(f5[("host", L)], o)) if o is not None else None
        h5 = fmetrics(f5[("host", L)])
        c2, c3, vn, vpf = c2_c3(f5[("PEDV_all", L)], o)
        vm = fmetrics(pframe(f5[("PEDV_all", L)], o)) if o is not None else None
        v5 = fmetrics(f5[("PEDV_all", L)])
        passed[L] = c1 and c2 and c3; off[L] = o
        host_shallow = n < MIN_HOST_START; vir_shallow = vn < MIN_PEDV_C2
        qc_status = "pass" if passed[L] else ("insufficient_depth" if (host_shallow or vir_shallow) else "fail")
        info[L] = dict(o=o, o0=o0, o1=o1, n=n, pk=pk, share=share, margin=margin, conf=conf, stab=stab, nb=nb, c1=c1,
                       level=level, hm=hm, h5=h5, c2=c2, c3=c3, vn=vn, vm=vm, v5=v5, qc=passed[L], qcs=qc_status,
                       host_shallow=host_shallow, vir_shallow=vir_shallow)
    OFFSET_TABLE[sid] = off
    sets = {"strict_core": [L for L in CORE if passed.get(L)], "fixed_29_30": list(CORE)}
    sets["expanded"] = sets["strict_core"] + [L for L in EXT if passed.get(L)]
    sets["short_RPF"] = [L for L in SHORT if passed.get(L)]
    def category(L):
        return ("strict_core" if L in sets["strict_core"] else "expanded_only" if L in sets["expanded"]
                else "short_RPF" if L in sets["short_RPF"] else "long_protected_fragment" if L >= LONG_MIN else "not_included")
    for L in sorted(rl):
        I = info[L]; r = rl[L]; o = I["o"]; Lk = L if L < 150 else "≥150"
        base = {**m, "read_length": Lk, "length_zone": zone(L)}
        vero = VERO_FROZEN.get(L)
        vero_s = vero if vero is not None else "NA(无冻结值)"
        same = (("一致" if o == vero else "不一致") if (vero is not None and o is not None) else "NA")
        hm, vm, h5, v5 = I["hm"], I["vm"], I["h5"], I["v5"]
        crit = (f"n={I['n']};peak={I['pk']};share={f3(I['share'])};margin={f2(I['margin'])};"
                f"half0={I['o0'] if I['o0'] is not None else 'NA'};half1={I['o1'] if I['o1'] is not None else 'NA'}")
        if L <= 40:
            offs.append({**base, "offset": o if o is not None else "NA", "判据数值": crit,
                         "support_reads_host_start": I["n"], "peak_reads": I["pk"], "peak_share": f3(I["share"]),
                         "peak_margin_over_2nd": f2(I["margin"]),
                         "offset_half0": I["o0"] if I["o0"] is not None else "NA", "offset_half1": I["o1"] if I["o1"] is not None else "NA",
                         "offset_peak_confidence": I["conf"], "offset_stability": I["stab"], "neighbor_consistency_flag(仅标记)": I["nb"],
                         "C1": I["c1"], "offset_confidence": I["level"],
                         "host_start_depth": ("insufficient_depth" if I["host_shallow"] else "ok"),
                         "Vero冻结值": vero_s, "PK1值": o if o is not None else "NA", "与Vero冻结值是否一致": same,
                         "category_PK1": category(L)})
        row = {**base, "host_reads": r["host_n"], "host_pct": f"{100*r['host_n']/htot:.3f}",
               "PEDV_plus_reads": r["pedv_plus_n"], "PEDV_plus_pct": f"{100*r['pedv_plus_n']/vtot:.3f}", "PEDV_minus_reads": r["pedv_minus_n"],
               "offset": o if o is not None else "NA", "offset_confidence": I["level"], "C1": I["c1"],
               "host_5p_frame_n": int(f5[("host", L)].sum()), "host_psite_frame0_1_2": fr3(hm["frac"]) if hm else "NA",
               "host_phase0_pct(CDS主体)": f"{100*hm['frac'][0]:.2f}" if hm else "NA",
               "host_dominant_frame": hm["dom"] if hm else "NA", "host_periodicity_ratio": f2(h5["ratio"]) if h5 else "NA",
               "host_correct_frame_enrichment": f2(hm["cfe"]) if hm else "NA",
               "PEDV_frame_n(ORF唯一区5′端)": I["vn"], "PEDV_psite_frame0_1_2": fr3(vm["frac"]) if vm else "NA",
               "PEDV_dominant_frame": vm["dom"] if vm else "NA", "PEDV_periodicity_ratio(C2指标)": f2(v5["ratio"]) if v5 else "NA",
               "PEDV_correct_frame_enrichment": f2(vm["cfe"]) if vm else "NA",
               "C2_periodicity": I["c2"] if L <= 45 else "NA", "C3_frame0": I["c3"] if L <= 45 else "NA", "QC_pass_C1_C2_C3": I["qc"],
               "qc_status": I["qcs"] if L <= 45 else "NA",
               "in_strict_core": L in sets["strict_core"], "in_fixed_29_30": L in sets["fixed_29_30"],
               "in_expanded": L in sets["expanded"], "in_short_RPF": L in sets["short_RPF"], "category": category(L)}
        src.append(row)
        # 每（样本, 读长）的 PEDV reads 数（含 insufficient_depth 标记）
        if L <= 60 or r["pedv_plus_n"] > 0:
            depth_rows.append({**base, "PEDV_plus_reads": r["pedv_plus_n"], "PEDV_minus_reads": r["pedv_minus_n"],
                               "PEDV_plus_clean_5p_reads": r["pedv_plus_clean_n"], "PEDV_5p_in_ORF_unique_body_reads": I["vn"],
                               "C2_depth_requirement": MIN_PEDV_C2,
                               "depth_flag": "insufficient_depth" if I["vir_shallow"] else "ok",
                               "host_start_codon_reads": I["n"],
                               "host_depth_flag": "insufficient_depth" if I["host_shallow"] else "ok"})
        if 31 <= L <= 33:
            alt = {}
            for oo in (11, 12, 13, 14):
                hmo = fmetrics(pframe(f5[("host", L)], oo)); vmo = fmetrics(pframe(f5[("PEDV_all", L)], oo))
                alt[oo] = (hmo, vmo)
            hx = hs[L][1:L + 1]
            split = (I["o0"] != I["o1"]) or (hm is not None and hm["frac"][0] < 0.40) or (o is None)
            split31.append({**base, "offset": o if o is not None else "NA", "offset_confidence": I["level"],
                            "offset_half0": I["o0"], "offset_half1": I["o1"], "support_reads": I["n"],
                            "startpos_j12_j13_j14_reads(offset11/12/13)": f"{hx[11]}/{hx[12]}/{hx[13]}",
                            "peak_margin_over_2nd": f2(I["margin"]),
                            **{f"host_frame0_if_offset_{oo}": f3(alt[oo][0]["frac"][0]) if alt[oo][0] else "NA" for oo in (12, 13)},
                            **{f"PEDV_frame0_if_offset_{oo}": f3(alt[oo][1]["frac"][0]) if alt[oo][1] else "NA" for oo in (12, 13)},
                            "host_5p_frame_f0_f1_f2": "/".join(map(str, f5[("host", L)])),
                            "C2": I["c2"], "C3": I["c3"], "category": category(L),
                            "frame_split_flag": "frame分裂-待确认（不硬给数）" if split else "无分裂迹象"})
    # ORF 逐读长深度（7 个 canonical ORF 的唯一区 5′端 reads）
    for dom in [f"PEDV_{n}" for n, a, b in ORFS]:
        for L in sorted(rl):
            nn = int(f5[(dom, L)].sum())
            if nn == 0 and not (26 <= L <= 35): continue
            orf_depth.append({**m, "ORF": dom[5:], "read_length": L, "reads_5p_in_ORF_unique_body": nn,
                              "depth_flag": "insufficient_depth" if nn < MIN_PEDV_C2 else "ok"})
    # Vero 差异清单（仅 29/30 nt 有冻结值）
    for L in CORE:
        o = off.get(L)
        vero_diff.append({**m, "read_length": L, "Vero冻结值": VERO_FROZEN[L], "PK1值": o if o is not None else "NA",
                          "与Vero冻结值是否一致": ("一致" if o == VERO_FROZEN[L] else "不一致") if o is not None else "NA",
                          "offset_confidence": info[L]["level"] if L in info else "NA",
                          "qc_status": info[L]["qcs"] if L in info else "NA"})
    # ---- 四套集合 ----
    w = {**m}
    pctv = lambda Ls: 100 * sum(rl[L]["pedv_plus_n"] for L in Ls if L in rl) / vtot
    pcth = lambda Ls: 100 * sum(rl[L]["host_n"] for L in Ls if L in rl) / htot
    for k in ("strict_core", "fixed_29_30", "expanded", "short_RPF"):
        Ls = sets[k]; w[k] = ",".join(f"{L}:{off[L]}" for L in Ls) or "NONE"
        w[f"{k}_pct_PEDV_plus"] = f"{pctv(Ls):.2f}"; w[f"{k}_pct_host"] = f"{pcth(Ls):.2f}"
    w["LPF_ge35_pct_PEDV_plus"] = f"{pctv([L for L in rl if L >= 35]):.2f}"
    sets_wide.append(w)
    for L in sorted(rl):
        if L > 40 and L < LONG_MIN: continue
        sets_long.append({**m, "read_length": L if L < 150 else "≥150", "length_zone": zone(L), "offset": off[L] if off[L] is not None else "NA",
                          "QC_pass_C1_C2_C3": passed[L], "qc_status": info[L]["qcs"] if L <= 45 else "NA",
                          "in_strict_core": L in sets["strict_core"], "in_fixed_29_30": L in sets["fixed_29_30"],
                          "in_expanded": L in sets["expanded"], "in_short_RPF": L in sets["short_RPF"], "category": category(L),
                          "PEDV_plus_reads": rl[L]["pedv_plus_n"],
                          "note": "分析口径，不是过滤器：该读长全部 reads 保留在 P-site 表中"})

prov = provenance("X", [
    "method: Module 1 M1_v1.0 方法 B（docs/METHODS_quality_gated_selection.md 原文）；判定函数逐行同 RGP_v1.3 / m1_build_tables.py",
    "host anchors: 猪 Sscrofa11.1.110 GTF 中全部单段 3 nt start_codon 去重（annotation/start_codon_pig.bed，方法 B 原定义）；canonical-only 敏感性集见 stats_canon/",
    "PEDV frame5: GI-b 7 个 canonical ORF（orf_check.tsv），ORF 两端各去 45 nt、去 ORF 重叠区；5′端 22 nt 内有剪接跳跃的 reads 不进 frame 统计（与 M1 同）",
    "Vero 冻结值（strict_core 29/30 nt = 12）只作对照列；不一致时以 PK1 实测为准",
    "scripts: psite/scripts/build_B_offset_sets.py（改自 RGP_v1.3；原版保留于 psite/scripts_orig/）"] + bam_prov())
P = lambda t: prov.replace("# X", f"# {t}")
write(f"{OUT}/psite_offset_PK1.tsv", offs, P("psite_offset_PK1（≤40 nt；≥41 nt 的 offset 见 length_qc_PK1.tsv）"))
write(f"{OUT}/length_qc_PK1.tsv", src, P("length_qc_PK1：逐样本×逐读长全部判据（source data）"))
write(f"{OUT}/readlength_sets_PK1.tsv", sets_long, P("readlength_sets_PK1：逐样本×读长成员表（分析口径，不是过滤器）"))
write(f"{OUT}/readlength_sets_PK1_wide.tsv", sets_wide, P("readlength_sets_PK1_wide：四套集合（读长:offset）"))
write(f"{OUT}/readlen31_33_offset_frame_PK1.tsv", split31, P("31–33 nt offset 与 frame 复核（frame 分裂标 待确认）"))
write(f"{OUT}/vero_vs_pk1_offset_diff.tsv", vero_diff, P("Vero 冻结值 vs PK1 实测 offset（29/30 nt）"))
write(f"{OUT}/pedv_depth_by_length_PK1.tsv", depth_rows, P("每（样本, 读长）的 PEDV reads 数与深度标记（insufficient_depth：ORF 唯一区 5′端 reads < 1000）"))
write(f"{OUT}/pedv_orf_depth_by_length_PK1.tsv", orf_depth, P("7 个 canonical ORF 唯一区逐读长 5′端 reads 数（深度不足标 insufficient_depth，不下翻译结论）"))
os.makedirs(f"{OUT}/_work", exist_ok=True)
json.dump({sid: {str(L): o for L, o in OFFSET_TABLE[sid].items()} for sid in IDS},
          open(f"{OUT}/_work/offset_table_PK1.json", "w"), indent=0)
for w in sets_wide: print(w["sample_id"], w["strict_core"], "| fixed_29_30", w["fixed_29_30"], "| exp", w["expanded"], "| short", w["short_RPF"])
for v in vero_diff: print("VERO", v["sample_id"], v["read_length"], v["Vero冻结值"], v["PK1值"], v["与Vero冻结值是否一致"], v["qc_status"])
