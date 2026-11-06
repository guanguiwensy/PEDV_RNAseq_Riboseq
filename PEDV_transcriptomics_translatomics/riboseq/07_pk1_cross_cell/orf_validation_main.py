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
import os, sys, gzip, csv, pickle, hashlib, collections, math
import numpy as np

ROOT = PEDV_P+""; HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
PREREG = HERE + "/PREREG_cross_cell_classification.md"
PREREG_SHA = hashlib.sha256(open(PREREG, "rb").read()).hexdigest()
VLEN = 27959
COND = [("CHX", 16), ("HAR", 16)]
VSID = {("CHX", 16): SR.lib("RIBO_VERO_CHX_1MOI_16h"), ("HAR", 16): SR.lib("RIBO_VERO_HAR_1MOI_16h")}; PSID = {("CHX", 16): SR.lib("RIBO_PK1_CHX_1MOI_16h"), ("HAR", 16): SR.lib("RIBO_PK1_HAR_1MOI_16h")}
HOST = {"Vero": {"CHX": (0.8454, 0.099, 0.0555), "HAR": (0.8637, 0.0815, 0.0548)}, "PK1": {"CHX": (0.9257, 0.0409, 0.0334), "HAR": (0.7454, 0.1430, 0.1116)}}
PARENT = {"ORF1a": (293, 12622), "ORF1b": (12592, 20613), "S": (20610, 24758), "ORF3": (24758, 25033), "E": (25361, 25600), "M": (25608, 26288), "N": (26300, 27625)}


cache = pickle.load(open(f"{ROOT}/SupM_alt_ORF_evidence/cache/psite_by_length.pkl", "rb"))["A"]
ARR = {"Vero": {d: np.asarray(cache[VSID[(d, 16)]][29], float)[:VLEN + 1] for d in ("CHX", "HAR")}, "PK1": {d: np.zeros(VLEN + 1) for d in ("CHX", "HAR")}}
rev = {v: k for k, v in PSID.items()}
with gzip.open(f"{ROOT}/PK1_psite/psite_counts_PEDV_GIb_PK1_29nt.tsv.gz", "rt", encoding="utf-8") as f:
    for line in f:
        if line.startswith("#") or line.startswith("sample_id"): continue
        p = line.rstrip("\n").split("\t")
        if p[3] != "+" or p[5] != "29" or p[0] not in rev: continue
        ARR["PK1"][rev[p[0]][0]][int(p[4])] += int(p[8])
TOT = {c: {d: ARR[c][d].sum() for d in ("CHX", "HAR")} for c in ARR}


EXCL = np.zeros(VLEN + 2, bool)
EXCL[12586:12623] = True; EXCL[26687:26718] = True; EXCL[20610:20614] = True; EXCL[24758:24759] = True
for n, (a, b) in PARENT.items():
    EXCL[max(1, a - 12):a + 31] = True; EXCL[b - 30:b + 31] = True
def parent_body(parent):
    a, b = PARENT[parent]; m = np.zeros(VLEN + 2, bool); m[a + 45:b - 44] = True
    ar = np.arange(VLEN + 2)
    for n2, (a2, b2) in PARENT.items():
        if n2 != parent: m &= ~((ar >= a2) & (ar <= b2))
    return m & ~EXCL
def frame_of(p, parent): return (p - PARENT[parent][0]) % 3
def f0p(arr, parent):
    m = parent_body(parent); idx = np.nonzero(m)[0]
    if len(idx) == 0 or arr[idx].sum() == 0: return float("nan")
    fr = (idx - PARENT[parent][0]) % 3
    return float(arr[idx][fr == 0].sum() / arr[idx].sum())


cl = [l for l in open(f"{ROOT}/PK1_orf_evidence/candidate_list_frozen.tsv", encoding="utf-8") if not l.startswith("#")]
CL = list(csv.DictReader(cl, delimiter="\t"))
FAM = [r for r in CL if r["source_group"] == "HARTIS_PRICE_family"]
famsum = {r["family_id"]: r for r in csv.DictReader([l for l in open(f"{ROOT}/HARTIS_v3/G4_Figure6_integration/HARTIS_PRICE_family_summary.tsv", encoding="utf-8") if not l.startswith("#")], delimiter="\t")}
PRIORITY = ["N-x1", "N-x2", "S-x1", "ORF1b-x13", "ORF1a-x22", "E-x1", "M-x1", "M-x2"]
def has_vero_evidence(r):
    e = r["vero_strongest_evidence_type"]; return e.startswith("PRICE nominal only") or e.startswith("HARTIS HAR W1 nominal only") or e.startswith("HARTIS capped_reference")
SELA = [r for r in FAM if r["alias_name"] in PRIORITY or has_vero_evidence(r)]
print("analysis A families:", len(SELA), "of", len(FAM))
def parse_fam(r):
    sv = [(int(x.split(":")[0]), x.split(":")[1]) for x in r["start_variants"].split(";") if x]
    parent = r["parent_region"]; k = 1 if r["relative_frame"] == "+1" else 2
    return dict(alias=r["alias_name"], fid=r["source_candidate_id"], parent=parent, k=k, starts=sv, stop=int(r["stop"]), ev=r["vero_strongest_evidence_type"], span=(min(s for s, _ in sv), int(r["stop"]) + 2))
FAMS = {r["alias_name"]: parse_fam(r) for r in FAM}
bad = []
for a, f in FAMS.items():
    for s, c in f["starts"]:
        if frame_of(s, f["parent"]) != f["k"]: bad.append((a, s, frame_of(s, f["parent"]), f["k"]))
print("frame inconsistencies (alias,start,frame,expected):", bad[:10], len(bad))


def region_gate(f, cell):
    if cell != "PK1": return None
    lo, hi = f["span"]
    if not (hi < 25361 or lo > 25600): return "E"
    if not (hi < 24758 or lo > 25033): return "ORF3"
    return None
def start_metrics(arr, s, k, parent):
    P = arr[s]; cand = [arr[p] for p in range(max(1, s - 30), min(VLEN, s + 30) + 1) if p != s and frame_of(p, parent) == k]
    mean = float(np.mean(cand)) if cand else 0.0
    fold = (P / mean) if mean > 0 else (float("inf") if P > 0 else 0.0)
    rank = 1 + sum(1 for x in cand if x > P)
    w5 = [0.0, 0.0, 0.0]
    for p in range(s, min(VLEN, s + 12) + 1): w5[frame_of(p, parent)] += arr[p]
    return dict(P=int(P), fold=fold, rank=rank, w5=w5, w5_total=sum(w5), m=w5[1] + w5[2])
def start_class(m, gate):
    if gate == "E" or gate == "ORF3": return "not_assessable", "区域门槛：%s 区唯一区深度 <1000" % gate
    if m["w5_total"] < 20: return "not_assessable", "W5 总 P-site %d < 20" % m["w5_total"]
    if m["P"] >= 100 and m["fold"] >= 20 and m["rank"] == 1: return "start_supported", "P_s=%d≥100, fold=%.1f≥20, rank=1" % (m["P"], m["fold"])
    if m["P"] >= 20 and m["fold"] >= 5 and m["rank"] <= 2: return "start_suggestive", "P_s=%d≥20, fold=%.1f≥5, rank=%d≤2" % (m["P"], m["fold"], m["rank"])
    return "start_not_detected", "P_s=%d, fold=%.1f, rank=%d 未达 suggestive" % (m["P"], m["fold"] if np.isfinite(m["fold"]) else 1e9, m["rank"])
def tier(m):
    mm = m["m"]
    return "below_floor(m<10)" if mm < 10 else ("low_m_limited(10-29)" if mm < 30 else "not_certifiable(m>=30)")
def body_metrics(arr, f, cell, drug):
    parent, k = f["parent"], f["k"]; lo, hi = min(s for s, _ in f["starts"]) + 30, f["stop"] - 30
    pos = np.arange(lo, hi + 1) if hi > lo else np.array([], int)
    pos = pos[~EXCL[pos]] if len(pos) else pos
    fr = (pos - PARENT[parent][0]) % 3 if len(pos) else pos
    cnt = [float(arr[pos][fr == j].sum()) if len(pos) else 0.0 for j in range(3)]; B = sum(cnt)

    m = parent_body(parent).copy(); m[lo:hi + 1] = False; idx = np.nonzero(m)[0]
    vtot = arr[idx].sum(); vfr = (idx - PARENT[parent][0]) % 3
    v_c = float(arr[idx][vfr == k].sum() / vtot) if vtot > 0 else float("nan")
    s_c = cnt[k] / B if B else float("nan"); s_o = cnt[3 - k] / B if B else float("nan")
    cp = sorted([arr[p] for p in pos if frame_of(p, parent) == k], reverse=True) if len(pos) else []
    npos = len(cp); nz = sum(1 for x in cp if x > 0)
    z_c = nz / npos if npos else float("nan")
    top3 = sum(cp[:3]); s_c3 = (cnt[k] - top3) / (B - top3) if B - top3 > 0 else float("nan")
    return dict(B=B, F0=cnt[0] / B if B else float("nan"), s_c=s_c, s_o=s_o, v_c=v_c, z_c=z_c, s_c3=s_c3, F0p=f0p(arr, parent), h_c=HOST[cell][drug][k], lo=lo, hi=hi, nz=nz, npos=npos)
def body_class(b, gate, parent):
    if gate == "E": return "not_assessable", "区域门槛：E 区 29/30 nt 唯一区深度 <1000"
    if b["B"] < 500: return "not_assessable", "体内 P-site B=%d < 500" % b["B"]
    if not (b["F0p"] >= 0.45): return "not_assessable", "父 ORF %s 的 canonical 框份额 F0p=%.2f < 0.45（背景接近 1/3，候选框富集无法解读）" % (parent, b["F0p"])
    r = []
    ok_s = (b["s_c"] >= 2.0 * b["h_c"]) and (b["s_c"] >= 1.5 * b["s_o"]) and (b["s_c"] >= 1.5 * b["v_c"]) and (b["z_c"] >= 0.5) and (b["s_c3"] >= 1.2 * b["v_c"])
    ok_g = (b["s_c"] >= 1.5 * b["h_c"]) and (b["s_c"] >= 1.2 * b["s_o"]) and (b["s_c"] >= 1.2 * b["v_c"]) and (b["z_c"] >= 0.3)
    txt = "s_c=%.3f, s_o=%.3f, 宿主 h_c=%.3f, 父 ORF v_c=%.3f, z_c=%.2f, s_c3=%.3f, F0p=%.2f" % (b["s_c"], b["s_o"], b["h_c"], b["v_c"], b["z_c"], b["s_c3"], b["F0p"])
    if ok_s: return "body_supported", txt
    if ok_g: return "body_suggestive", txt
    return "body_not_detected", txt
ORD = {"not_assessable": 0, "start_not_detected": 1, "body_not_detected": 1, "start_suggestive": 2, "body_suggestive": 2, "start_supported": 3, "body_supported": 3}
def overall(sc, bc):
    a, b = ORD[sc], ORD[bc]
    if a == 0 and b == 0: return "not_assessable"
    if a == 3 and b == 3: return "supported"
    if a >= 2 or b >= 2: return "suggestive"
    return "not_detected"

def evaluate(f, cell):
    gate = region_gate(f, cell); harr, carr = ARR[cell]["HAR"], ARR[cell]["CHX"]; variants = []
    for s, codon in f["starts"]:
        m = start_metrics(harr, s, f["k"], f["parent"]); cls, why = start_class(m, gate); variants.append(dict(s=s, codon=codon, m=m, cls=cls, why=why, tier=tier(m)))
    order = {"start_supported": 3, "start_suggestive": 2, "start_not_detected": 1, "not_assessable": 0}
    best = max(variants, key=lambda v: (order[v["cls"]], v["m"]["P"])) if variants else None
    b = body_metrics(carr, f, cell, "CHX"); bcls, bwhy = body_class(b, gate, f["parent"])
    return dict(gate=gate, variants=variants, best=best, body=b, bcls=bcls, bwhy=bwhy, overall=overall(best["cls"], bcls))
EVAL = {}
for r in FAM:
    f = FAMS[r["alias_name"]]
    EVAL[f["alias"]] = {cell: evaluate(f, cell) for cell in ("PK1", "Vero")}
fmt = lambda x, nd=3: ("%." + str(nd) + "f") % x if isinstance(x, float) and np.isfinite(x) else ("inf" if isinstance(x, float) and np.isinf(x) else str(x))
def rowfor(alias):
    f = FAMS[alias]; fs = famsum[f["fid"]]; row = dict(family=alias, family_id=f["fid"], priority="是" if alias in PRIORITY else "", parent=f["parent"], frame="+%d" % f["k"], stop_codon_first_nt=f["stop"],
                                                n_start_variants=len(f["starts"]), vero_frozen_evidence=("PRICE_status=%s; HAR_W1_status=%s; joint=%s" % (fs["PRICE_status"], fs["HAR_W1_status"], fs["joint_category_exploratory"])),
                                                prereg_sha256=PREREG_SHA[:16])
    for cell in ("PK1", "Vero"):
        e = EVAL[alias][cell]; b = e["best"]; bd = e["body"]
        pre = cell + "_"
        row.update({pre + "region_gate": e["gate"] or "", pre + "HAR_best_start": "%d %s" % (b["s"], b["codon"]), pre + "HAR_start_P": b["m"]["P"], pre + "HAR_start_fold": fmt(b["m"]["fold"], 1), pre + "HAR_start_rank": b["m"]["rank"],
                    pre + "HAR_W5_n0/n1/n2": "%d/%d/%d" % tuple(int(x) for x in b["m"]["w5"]), pre + "HAR_host_cert_tier": b["tier"], pre + "HAR_start_class": b["cls"], pre + "HAR_start_reason": b["why"],
                    pre + "CHX_body_window": "%d-%d" % (bd["lo"], bd["hi"]), pre + "CHX_body_B": int(bd["B"]), pre + "CHX_parent_F0p": fmt(bd["F0p"], 2), pre + "CHX_s_c": fmt(bd["s_c"]), pre + "CHX_s_o": fmt(bd["s_o"]), pre + "CHX_host_h_c": fmt(bd["h_c"]),
                    pre + "CHX_parent_v_c": fmt(bd["v_c"]), pre + "CHX_z_c": fmt(bd["z_c"], 2), pre + "CHX_s_c3": fmt(bd["s_c3"]), pre + "CHX_body_class": e["bcls"], pre + "CHX_body_reason": e["bwhy"], pre + "overall_class": e["overall"]})
    ceq = famsum[f["fid"]].get("candidate_frame_equals_canonical", "")
    row["candidate_frame_equals_canonical"] = ceq
    if ceq:
        row["Vero_overall_class"] = "excluded_same_frame_as_canonical(%s)" % ceq; row["Vero_note"] = "候选框与 canonical ORF %s 同框；Vero 冻结分析（PRICE frame_mismatch）已排除，本规则对 Vero 的计算值保留但不解释" % ceq
    row["PRICE_PK1"] = "not_run（无猪 gedi 索引与宿主训练模型；非 not_detected）"; row["HARTIS_pvalue_PK1"] = "not_computed（HAR 16 h 宿主背景不足，不套用 Vero 背景）"
    row["note"] = "探索性、未校正；n=1；PK1 CHX 29 nt 周期性不足；ORF1a/ORF1b 的 +1/+2 不据此宣称移码翻译" if f["parent"] in ("ORF1a", "ORF1b") else "探索性、未校正；n=1；PK1 CHX 29 nt 周期性不足"
    return row
A_ROWS = [rowfor(r["alias_name"]) for r in SELA]
def wtsv(name, rows, header):
    with open(f"{OUT}/{name}", "w", encoding="utf-8", newline="") as fh:
        fh.write("\n".join(header) + "\n"); w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(rows)
H = ["# Classification rules: PREREG_cross_cell_classification.md (sha256 %s)" % PREREG_SHA,
     f"# Samples: PK1 16 h ({SR.lib('RIBO_PK1_CHX_1MOI_16h')}, {SR.lib('RIBO_PK1_HAR_1MOI_16h')}) versus Vero 16 h ({SR.lib('RIBO_VERO_CHX_1MOI_16h')}, {SR.lib('RIBO_VERO_HAR_1MOI_16h')}); 29 nt; n=1; MOI=1 for both cell lines",
     "# Exploratory, unadjusted results. Supported/suggestive labels indicate consistency with candidate-frame enrichment, not statistical significance. Lack of detection does not imply lack of translation. PK1 CHX 29-nt periodicity is weak (frame 0 approximately 0.47)."]
wtsv("VERO_KNOWN_ORF_PK1.tsv", A_ROWS, H)

vrows = []
for r in FAM:
    f = FAMS[r["alias_name"]]
    for cell in ("PK1", "Vero"):
        e = EVAL[f["alias"]][cell]
        for v in e["variants"]:
            vrows.append({"family": f["alias"], "in_analysis_A": "是" if r in SELA else "", "cell": cell, "start": v["s"], "codon": v["codon"], "P_s": v["m"]["P"], "fold": fmt(v["m"]["fold"], 1), "rank": v["m"]["rank"], "W5_n0/n1/n2": "%d/%d/%d" % tuple(int(x) for x in v["m"]["w5"]),
                          "W5_m": int(v["m"]["m"]), "host_cert_tier": v["tier"], "class": v["cls"], "reason": v["why"]})
wtsv("start_variants_PK1_Vero.tsv", vrows, H)


selA = {r["alias_name"] for r in SELA}
C_ROWS = []
ORDV = {"supported": 3, "suggestive": 2, "not_detected": 1, "not_assessable": 0}
for r in FAM:
    a = r["alias_name"]
    if a in selA: continue
    ep, ev = EVAL[a]["PK1"], EVAL[a]["Vero"]
    cat = ""
    if ep["overall"] == "supported" and ev["overall"] not in ("supported", "suggestive"): cat = "1: PK1 supported，Vero 不是 supported/suggestive"
    elif ep["overall"] == "suggestive" and ev["overall"] in ("not_detected",): cat = "2: PK1 suggestive，Vero not_detected"
    elif ep["overall"] in ("supported", "suggestive") and ev["overall"] == "not_assessable": cat = "4: Vero not_assessable，无法判断是否 PK1 特有"
    elif ep["overall"] in ("supported", "suggestive"): cat = "同类：Vero 也是 supported/suggestive（非新）"
    C_ROWS.append((ORDV[ep["overall"]], ep["best"]["m"]["fold"] if np.isfinite(ep["best"]["m"]["fold"]) else 1e9, ep["best"]["m"]["P"], a, cat))
C_ROWS.sort(key=lambda x: (-x[0], -x[1], -x[2]))
top = [x for x in C_ROWS if x[0] >= 2][:10]
crow = []
for _, _, _, a, cat in top:
    row = rowfor(a); row["category(PREREG §5)"] = cat; row["rank_rule"] = "总类 > HAR 起始 fold > P_s；探索性、未校正，不得写成显著"; crow.append(row)
cat_counts = collections.Counter(x[4] for x in C_ROWS if x[0] >= 2)
if not crow:
    crow = [{"result": "没有任何不在分析 A 里的家族在 PK1 16 h 达到 suggestive 或 supported（按冻结规则）", "n_families_scanned": len(C_ROWS)}]
wtsv("PK1_NEW_ORF_CANDIDATES.tsv", crow, H + ["# 新候选扫描：既有 74 个家族 578 个起点；不在分析 A 的家族 %d 个；PK1 达到 suggestive/supported 的家族 %d 个；类别计数 %s" % (len(C_ROWS), sum(1 for x in C_ROWS if x[0] >= 2), dict(cat_counts)),
                                                    "# PRICE 在 PK1 上未运行；条件差异检验（类别 3）未做（n=1，无重复，待确认）"])
print("C scan: non-A families", len(C_ROWS), "PK1 >=suggestive", sum(1 for x in C_ROWS if x[0] >= 2), dict(cat_counts))


HOT = [("5UTR", 22), ("5UTR", 113), ("5UTR", 271), ("3UTR", 27707), ("3UTR", 27820), ("3UTR", 27919), ("3UTR", 27920)]
def win(arr, lo, hi): return arr[lo:hi + 1]
B_ROWS = []
for reg, h in HOT:
    for drug in ("CHX", "HAR"):
        v, p = ARR["Vero"][drug], ARR["PK1"][drug]; lo, hi = h - 10, min(VLEN, h + 10)
        cv, cp = v / TOT["Vero"][drug] * 1e6, p / TOT["PK1"][drug] * 1e6
        pv = int(np.argmax(v[lo:hi + 1])) + lo; pp = int(np.argmax(p[lo:hi + 1])) + lo
        bg = np.concatenate([p[max(1, h - 100):lo], p[hi + 1:min(VLEN, h + 100) + 1]]); bgmean = float(bg.mean()) if len(bg) else 0.0
        bgv = np.concatenate([v[max(1, h - 100):lo], v[hi + 1:min(VLEN, h + 100) + 1]]); bgvm = float(bgv.mean()) if len(bgv) else 0.0
        P = p[pp]; Pv = v[pv]; fold = P / bgmean if bgmean > 0 else float("inf")
        tot = p[lo:hi + 1].sum()
        if tot < 20: cls, why = "not_assessable", "PK1 窗口总 P-site %d < 20" % tot
        elif abs(pp - pv) <= 1 and P >= 20 and fold >= 10: cls, why = "reproduced", "峰位置差 %d nt, PK1 峰 %d≥20, 峰/背景均值 %s≥10" % (pp - pv, P, fmt(fold, 1))
        elif abs(pp - pv) <= 3 and P >= 20: cls, why = "weak", "峰位置差 %d nt（≤3），PK1 峰 %d≥20，峰/背景均值 %s" % (pp - pv, P, fmt(fold, 1))
        else: cls, why = "not_reproduced", "峰位置差 %d nt，PK1 峰 %d，峰/背景均值 %s 不满足 weak" % (pp - pv, P, fmt(fold, 1))
        B_ROWS.append({"type": "hotspot", "region": reg, "target": "%d" % h, "drug": drug, "Vero_peak_pos": pv, "Vero_peak_Psites": int(Pv), "Vero_peak_CPM29": round(float(cv[pv]), 1), "Vero_bg_mean": round(bgvm, 3),
                       "PK1_peak_pos": pp, "PK1_peak_Psites": int(P), "PK1_peak_CPM29": round(float(cp[pp]), 1), "PK1_bg_mean": round(bgmean, 3), "peak_shift_nt": pp - pv, "PK1_fold_over_bg": fmt(fold, 1), "PK1_window_total_Psites": int(tot), "class": cls, "reason": why,
                       "note": "热点是 P-site 位置，不是起始密码子；n=1；PK1 CHX 29 nt 周期性不足" if drug == "CHX" else "热点是 P-site 位置，不是起始密码子；n=1"})

hotwin = [(h - 10, h + 10) for _, h in HOT]
for reg, lo, hi in (("5UTR", 1, 292), ("3UTR", 27626, 27959)):
    for drug in ("CHX", "HAR"):
        p = ARR["PK1"][drug]; seg = [(p[x], x) for x in range(lo, hi + 1) if not any(a <= x <= b for a, b in hotwin)]
        seg.sort(reverse=True)
        for P, x in seg[:5]:
            bg = np.concatenate([p[max(1, x - 100):x - 10], p[x + 11:min(VLEN, x + 100) + 1]]); bm = float(bg.mean()) if len(bg) else 0.0
            fold = P / bm if bm > 0 else float("inf")
            B_ROWS.append({"type": "PK1_top_nonhotspot_site", "region": reg, "target": str(x), "drug": drug, "PK1_peak_pos": x, "PK1_peak_Psites": int(P), "PK1_peak_CPM29": round(float(P / TOT["PK1"][drug] * 1e6), 1), "PK1_bg_mean": round(bm, 3), "PK1_fold_over_bg": fmt(fold, 1),
                           "Vero_peak_Psites": int(ARR["Vero"][drug][x]), "class": "listed_only", "reason": "PREREG §4：只列出，不做“新热点”结论（P-site≥50 且 峰/背景≥20 者标 candidate_for_review）" + ("；candidate_for_review" if (P >= 50 and fold >= 20) else ""), "note": "n=1；探索性"})

frz = {}
for r in CL:
    if r["source_group"] in ("3UTR_audit", "gRNA_5prime_package", "terminal_ORF") and isinstance(r["start"], str) and r["start"].isdigit():
        frz.setdefault(int(r["start"]), []).append(r)
UTR_STARTS = [(27629, "AUG（手稿 3′UTR 候选）"), (27797, "AUU（手稿 3′UTR 候选）"), (27866, "AUC（手稿 3′UTR 候选）"), (27706, "UUG（最近的冻结起点，对应热点 27707）"), (27710, "CUG（最近的冻结起点，对应热点 27707）"),
              (27823, "ACG（最近的冻结起点，对应热点 27820）"), (27918, "CUG（最近的冻结起点，对应热点 27919/27920）")]
gr = sorted({int(r["start"]) for r in CL if r["source_group"] == "gRNA_5prime_package"})
UTR_STARTS += [(s, "gRNA 5′ 候选") for s in gr]
def utr_start_row(s, label):
    rr = frz.get(s, []); stopv = int(rr[0]["stop"]) if rr else None; codon = rr[0]["start_codon"] if rr else "?"; src = ",".join(sorted({x["source_group"] for x in rr})) or "冻结表中无完全相同起点"
    out = {"type": "UTR_candidate_start", "region": "3UTR" if s > 27625 else "5UTR", "target": "%d %s" % (s, codon), "label": label, "frozen_sources": src, "frozen_stop": stopv if stopv else ""}
    for cell in ("Vero", "PK1"):
        arr = ARR[cell]["HAR"]; cand = [arr[p] for p in range(max(1, s - 30), min(VLEN, s + 30) + 1) if p != s and (p - s) % 3 == 0]
        mean = float(np.mean(cand)) if cand else 0.0; P = arr[s]; fold = P / mean if mean > 0 else (float("inf") if P > 0 else 0.0); rank = 1 + sum(1 for x in cand if x > P)
        w5 = float(arr[s:s + 13].sum())
        if w5 < 20: cls, why = "not_assessable", "W5 总 P-site %d < 20" % w5
        elif P >= 100 and fold >= 20 and rank == 1: cls, why = "start_supported", "P_s=%d, fold=%s, rank=%d" % (P, fmt(fold, 1), rank)
        elif P >= 20 and fold >= 5 and rank <= 2: cls, why = "start_suggestive", "P_s=%d, fold=%s, rank=%d" % (P, fmt(fold, 1), rank)
        else: cls, why = "start_not_detected", "P_s=%d, fold=%s, rank=%d" % (P, fmt(fold, 1), rank)

        hi = stopv if stopv else s + 60; seg = np.arange(s, min(VLEN, hi) + 1); fr = (seg - s) % 3; carr = ARR[cell]["CHX"]
        cnt = [float(carr[seg][fr == j].sum()) for j in range(3)]; B = sum(cnt)
        out.update({cell + "_HAR_W1_P": int(P), cell + "_HAR_fold": fmt(fold, 1), cell + "_HAR_rank": rank, cell + "_HAR_W5_total": int(w5), cell + "_HAR_class": cls, cell + "_HAR_reason": why,
                    cell + "_CHX_inframe/other1/other2_share": ("%.3f/%.3f/%.3f" % tuple(c / B for c in cnt)) if B else "NA", cell + "_CHX_cand_region_Psites": int(B)})
    out["class"] = "PK1:%s | Vero:%s" % (out["PK1_HAR_class"], out["Vero_HAR_class"]); out["reason"] = "PREREG §4；3′UTR 内无 canonical 框，不做 CHX 体判定，只列候选框份额；n=1；探索性"
    out["note"] = "候选框 = 以该起点为 0 的读框；PK1 CHX 29 nt 周期性不足"
    return out
for s, lab in UTR_STARTS: B_ROWS.append(utr_start_row(s, lab))
keys = list(dict.fromkeys(k for r in B_ROWS for k in r))
with open(f"{OUT}/UTR_SIGNAL_VALIDATION.tsv", "w", encoding="utf-8", newline="") as fh:
    fh.write("\n".join(H + ["# 分析 B：热点窗口 ±10 nt；UTR 候选起点；热点是 P-site 位置，不是起始密码子；冻结表没有完全相同起点时并列最近的起点，不合并"]) + "\n")
    w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", restval=""); w.writeheader(); w.writerows(B_ROWS)
CEQ = {r['alias_name']: famsum[r['source_candidate_id']].get('candidate_frame_equals_canonical', '') for r in FAM}
pickle.dump(dict(CEQ=CEQ, EVAL=EVAL, FAMS=FAMS, SELA=[r["alias_name"] for r in SELA], TOP=[x[3] for x in top], ARR=ARR, TOT=TOT, UTR=B_ROWS, PARENT=PARENT), open(f"{OUT}/_orf_eval.pkl", "wb"))
cnt = collections.Counter((a, EVAL[a]["PK1"]["overall"]) for a in [r["alias_name"] for r in SELA])
print("A PK1 overall:", collections.Counter(EVAL[r["alias_name"]]["PK1"]["overall"] for r in SELA), "Vero:", collections.Counter(EVAL[r["alias_name"]]["Vero"]["overall"] for r in SELA))
print("done")
