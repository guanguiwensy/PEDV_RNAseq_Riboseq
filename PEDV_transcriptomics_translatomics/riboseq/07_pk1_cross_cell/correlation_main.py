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
import os, sys, gzip, csv, pickle, json, collections
import numpy as np
from scipy.stats import rankdata, pearsonr

ROOT = PEDV_P+""
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = HERE + "/data"
OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
VLEN = 27959; RNG = np.random.default_rng(20261103); NBOOT = 500
COND = [("CHX", 16), ("HAR", 16)]
CL = lambda c: "%s %d h" % c
VSID = {("CHX", 16): SR.lib("RIBO_VERO_CHX_1MOI_16h"), ("HAR", 16): SR.lib("RIBO_VERO_HAR_1MOI_16h")}
PSID = {("CHX", 16): SR.lib("RIBO_PK1_CHX_1MOI_16h"), ("HAR", 16): SR.lib("RIBO_PK1_HAR_1MOI_16h")}
VRNA = {("CHX", 16): SR.lib("RNA_VERO_CHX_1MOI_16h"), ("HAR", 16): SR.lib("RNA_VERO_HAR_1MOI_16h")}
PRNA = {("CHX", 16): SR.lib("RNA_PK1_CHX_1MOI_16h"), ("HAR", 16): SR.lib("RNA_PK1_HAR_1MOI_16h")}

REG = collections.OrderedDict([("whole_virus", (1, 27959)), ("5UTR", (1, 292)), ("ORF1a", (293, 12622)), ("ORF1b", (12623, 20609)), ("S", (20610, 24757)),
                               ("ORF3_E_M", (24758, 26299)), ("N", (26300, 27625)), ("3UTR", (27626, 27959))])
FC = {0: "#7B3294", 1: "#56B4E9", 2: "#999999"}


cache = pickle.load(open(f"{ROOT}/SupM_alt_ORF_evidence/cache/psite_by_length.pkl", "rb"))["A"]
VERO = {L: {c: np.asarray(cache[VSID[c]][L], float)[:VLEN + 1] for c in COND} for L in (29, 30)}
PK1 = {L: {c: np.zeros(VLEN + 1) for c in COND} for L in (29, 30)}
rev = {v: k for k, v in PSID.items()}
with gzip.open(f"{ROOT}/PK1_psite/psite_counts_PEDV_GIb_PK1.tsv.gz", "rt", encoding="utf-8") as f:
    for line in f:
        if line.startswith("#") or line.startswith("sample_id"): continue
        p = line.rstrip("\n").split("\t")
        if p[3] != "+" or p[5] not in ("29", "30"): continue
        if rev[p[0]] not in COND: continue
        PK1[int(p[5])][rev[p[0]]][int(p[4])] += int(p[8])
def load_depth(path):
    a = np.zeros(VLEN + 1); op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        for l in f:
            c = l.rstrip("\n").split("\t")
            if len(c) == 3: a[int(c[1])] = float(c[2])
            elif len(c) == 2: a[int(c[0])] = float(c[1])
    return a
RNA_V = {c: load_depth(f"{DATA}/vero_rna_depth_samecond_{VRNA[c]}.tsv.gz") for c in COND}
RNA_P = {c: load_depth(f"{DATA}/rna_depth_{PRNA[c]}.tsv.gz") for c in COND}
cpm = lambda a: a / a.sum() * 1e6
CV = {L: {c: cpm(VERO[L][c]) for c in COND} for L in (29, 30)}; CP = {L: {c: cpm(PK1[L][c]) for c in COND} for L in (29, 30)}
RV = {c: RNA_V[c] / RNA_V[c].sum() * 1e6 for c in COND}; RP = {c: RNA_P[c] / RNA_P[c].sum() * 1e6 for c in COND}


def pear(x, y):
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0: return np.nan
    return float(np.corrcoef(x, y)[0, 1])
def spear(x, y):
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0: return np.nan
    return pear(rankdata(x), rankdata(y))
def binsum(a, lo, hi, b):
    v = a[lo:hi + 1]; n = len(v) // b * b
    return v[:n].reshape(-1, b).sum(1)
def block_boot(fn, x, y, block, nboot=NBOOT):
    n = len(x); block = max(5, min(block, n // 3)); nb = int(np.ceil(n / block)); res = []
    starts_max = n - block
    for _ in range(nboot):
        st = RNG.integers(0, starts_max + 1, nb)
        idx = (st[:, None] + np.arange(block)[None, :]).ravel()[:n]
        v = fn(x[idx], y[idx])
        if not np.isnan(v): res.append(v)
    if len(res) < 50: return (np.nan, np.nan)
    return (float(np.percentile(res, 2.5)), float(np.percentile(res, 97.5)))
rows = []
def add(analysis, level, target, cond, metric, val, lo=np.nan, hi=np.nan, n_coord="", n_either="", n_both="", norm="", note=""):
    rows.append(dict(analysis=analysis, level=level, target=target, condition=cond, metric=metric, value=("" if val is None or (isinstance(val, float) and np.isnan(val)) else round(val, 4)),
                     ci95_lo=("" if np.isnan(lo) else round(lo, 4)), ci95_hi=("" if np.isnan(hi) else round(hi, 4)), n_coordinates=n_coord, n_nonzero_either=n_either,
                     frac_nonzero_either=(round(n_either / n_coord, 4) if n_coord and n_either != "" else ""), n_nonzero_both=n_both, normalization=norm, note=note))

PAIRNOTE = "MOI=1 for PK1 and Vero"
CHXNOTE = "PK1 CHX 29 nt 周期性不足（C2 fail）"
def rnote(c):
    return PAIRNOTE + ("；" + CHXNOTE if c[0] == "CHX" else "")


for L in (29, 30):
    main = (L == 29)
    for c in COND:
        x_all, y_all = CV[L][c], CP[L][c]
        for rname, (lo, hi) in REG.items():
            x = x_all[lo:hi + 1]; y = y_all[lo:hi + 1]; n = hi - lo + 1
            nz_e = int(((x > 0) | (y > 0)).sum()); nz_b = int(((x > 0) & (y > 0)).sum())
            base = dict(cond=CL(c), n_coord=n, n_either=nz_e, n_both=nz_b)
            lvl = "region_L%d" % L
            nm = "CPM%d（P-sites/10^6 个 %d nt PEDV(+) P-sites）" % (L, L)

            r_p = pear(x, y); r_s = spear(x, y)
            ci_p = block_boot(pear, x, y, 300) if main else (np.nan, np.nan); ci_s = block_boot(spear, x, y, 300) if main else (np.nan, np.nan)
            add("7.1", lvl, rname, base["cond"], "Pearson_perpos_allpos", r_p, *ci_p, n, nz_e, nz_b, nm, rnote(c))
            add("7.1", lvl, rname, base["cond"], "Spearman_perpos_allpos", r_s, *ci_s, n, nz_e, nz_b, nm, rnote(c))

            m = (x > 0) | (y > 0)
            add("7.4", lvl, rname, base["cond"], "Pearson_perpos_nonzero_union", pear(x[m], y[m]), n_coord=n, n_either=nz_e, n_both=nz_b, norm=nm, note="仅 Vero 或 PK1 非零的位点")
            add("7.4", lvl, rname, base["cond"], "Spearman_perpos_nonzero_union", spear(x[m], y[m]), n_coord=n, n_either=nz_e, n_both=nz_b, norm=nm, note="仅 Vero 或 PK1 非零的位点")
            # C log1p
            add("7.4", lvl, rname, base["cond"], "Pearson_log1p_perpos", pear(np.log1p(x), np.log1p(y)), n_coord=n, n_either=nz_e, n_both=nz_b, norm="log1p(CPM%d)" % L, note="降低高峰主导")

            for b in (3, 30):
                xb, yb = binsum(x_all, lo, hi, b), binsum(y_all, lo, hi, b)
                ci_pb = block_boot(pear, xb, yb, max(1, 300 // b)) if (main and b == 30) else (np.nan, np.nan)
                ci_sb = block_boot(spear, xb, yb, max(1, 300 // b)) if (main and b == 30) else (np.nan, np.nan)
                add("7.4", lvl, rname, base["cond"], "Pearson_bin%dnt" % b, pear(xb, yb), *ci_pb, len(xb), int(((xb > 0) | (yb > 0)).sum()), int(((xb > 0) & (yb > 0)).sum()), "每 %d nt 求和后 CPM%d" % (b, L), "分箱，不平滑；n_coordinates 为箱数")
                add("7.4", lvl, rname, base["cond"], "Spearman_bin%dnt" % b, spear(xb, yb), *ci_sb, len(xb), int(((xb > 0) | (yb > 0)).sum()), int(((xb > 0) & (yb > 0)).sum()), "每 %d nt 求和后 CPM%d" % (b, L), "分箱，不平滑；n_coordinates 为箱数")

            tx, ty = np.percentile(x[x > 0], 99.5) if (x > 0).any() else np.inf, np.percentile(y[y > 0], 99.5) if (y > 0).any() else np.inf
            mk = (x <= tx) & (y <= ty)
            add("7.4", lvl, rname, base["cond"], "Spearman_perpos_mask_top0.5pct", spear(x[mk], y[mk]), n_coord=int(mk.sum()), n_either=nz_e, n_both=nz_b, norm=nm, note="掩蔽两库各自非零位点中最高 0.5% 的位点")
            add("7.4", lvl, rname, base["cond"], "Pearson_perpos_mask_top0.5pct", pear(x[mk], y[mk]), n_coord=int(mk.sum()), n_either=nz_e, n_both=nz_b, norm=nm, note="掩蔽两库各自非零位点中最高 0.5% 的位点")

            if main:
                xb, yb = binsum(x_all, lo, hi, 30), binsum(y_all, lo, hi, 30)
                rv, rp = binsum(RV[c], lo, hi, 30), binsum(RP[c], lo, hi, 30)
                def resid(v, r):
                    X = np.vstack([np.ones_like(r), np.log1p(r)]).T
                    beta = np.linalg.lstsq(X, np.log1p(v), rcond=None)[0]; return np.log1p(v) - X @ beta
                if len(xb) > 5:
                    ex, ey = resid(xb, rv), resid(yb, rp)
                    add("7.4", lvl, rname, base["cond"], "Pearson_bin30_RNAresid", pear(ex, ey), *(block_boot(pear, ex, ey, 10)), len(xb), int(((xb > 0) | (yb > 0)).sum()), int(((xb > 0) & (yb > 0)).sum()), "log1p(30-nt CPM29) 对各库自身 log RNA 密度（per 10^6 summed viral depth）做线性残差后相关", "RNA 丰度控制后的敏感性")
                    te_x, te_y = xb / (rv + 1e-9), yb / (rp + 1e-9)
                    add("7.4", lvl, rname, base["cond"], "Spearman_bin30_P_over_RNA", spear(te_x, te_y), n_coord=len(xb), n_either=int(((xb > 0) | (yb > 0)).sum()), n_both=int(((xb > 0) & (yb > 0)).sum()), norm="30-nt P-site CPM29 / 30-nt RNA 密度", note="P-site 对 RNA 密度的比值（翻译效率样指标，不是翻译效率）")

    if main:
        for c in COND:
            x_all, y_all = CV[29][c], CP[29][c]
            keep = np.ones(VLEN + 1, bool); keep[0] = False
            for nm_, (lo, hi) in (("M", (25608, 26288)), ("N", (26300, 27625))): keep[lo:hi + 1] = False
            x, y = x_all[keep], y_all[keep]
            ci_s = block_boot(spear, x, y, 300)
            add("7.1", "whole_virus_sensitivity", "excl_M_N", CL(c), "Spearman_perpos_allpos", spear(x, y), *ci_s, int(keep.sum()), int(((x > 0) | (y > 0)).sum()), int(((x > 0) & (y > 0)).sum()), "CPM29", "去掉 M（25608–26288）与 N（26300–27625）后的全病毒；" + rnote(c))
            add("7.1", "whole_virus_sensitivity", "excl_M_N", CL(c), "Pearson_perpos_allpos", pear(x, y), n_coord=int(keep.sum()), norm="CPM29")

            xs, ys = [], []
            for rn in ("5UTR", "ORF1a", "ORF1b", "S", "ORF3_E_M", "N", "3UTR"):
                lo, hi = REG[rn]; xv, yv = x_all[lo:hi + 1], y_all[lo:hi + 1]
                xs.append(xv / xv.sum() if xv.sum() else xv); ys.append(yv / yv.sum() if yv.sum() else yv)
            xe, ye = np.concatenate(xs), np.concatenate(ys)
            add("7.1", "whole_virus_sensitivity", "region_equal_weight", CL(c), "Spearman_perpos_allpos", spear(xe, ye), *block_boot(spear, xe, ye, 300), len(xe), int(((xe > 0) | (ye > 0)).sum()), int(((xe > 0) & (ye > 0)).sum()), "每个区域内先除以区域总和，再拼接", "各区域等权，防止 N/M 高覆盖区主导；" + rnote(c))
            add("7.1", "whole_virus_sensitivity", "region_equal_weight", CL(c), "Pearson_perpos_allpos", pear(xe, ye), n_coord=len(xe), norm="每个区域内先除以区域总和，再拼接")


cl = [l for l in open(f"{ROOT}/PK1_orf_evidence/candidate_list_frozen.tsv", encoding="utf-8") if not l.startswith("#")]
CAND = {}
for r in csv.DictReader(cl, delimiter="\t"):
    if r["alias_name"] in ("N-x1", "N-x2"):
        CAND[r["alias_name"]] = ([int(x.split(":")[0]) for x in r["start_variants"].split(";") if x], int(r["stop"]))
frame = lambda p: (p - 26300) % 3
KEY = []
for name, (starts, stop) in CAND.items():
    s0 = min(starts); lo, hi = s0 + 30, stop - 30
    pos = np.arange(lo, hi + 1)
    for c in COND:
        x_all, y_all = CV[29][c], CP[29][c]
        x, y = x_all[lo:hi + 1], y_all[lo:hi + 1]
        base = dict(n_coord=len(pos), n_either=int(((x > 0) | (y > 0)).sum()), n_both=int(((x > 0) & (y > 0)).sum()))
        nm = "CPM29"
        add("7.2", "key_ORF_body", name, CL(c), "Spearman_allframes_perpos", spear(x, y), *block_boot(spear, x, y, 60), norm=nm, note="ORF 体窗口 %d-%d；块长 60 nt（窗口较短）；%s" % (lo, hi, rnote(c)), **base)
        add("7.2", "key_ORF_body", name, CL(c), "Pearson_allframes_perpos", pear(x, y), *block_boot(pear, x, y, 60), norm=nm, note="同上", **base)
        for fr in (0, 1, 2):
            m = np.array([frame(p) == fr for p in pos]); xf, yf = x[m], y[m]
            add("7.2", "key_ORF_body", name, CL(c), "Spearman_frame%d_perpos" % fr, spear(xf, yf), *block_boot(spear, xf, yf, 20), n_coord=int(m.sum()), n_either=int(((xf > 0) | (yf > 0)).sum()), n_both=int(((xf > 0) & (yf > 0)).sum()), norm=nm,
                note=("候选读框(+1)" if fr == 1 else "frame%d" % fr) + "；frame=(pos-26300) mod 3；" + rnote(c))

        for cellname, a in (("Vero", x), ("PK1", y)):
            tot = sum(a[i] for i in range(len(a)))
            fs = [sum(a[i] for i in range(len(a)) if frame(pos[i]) == k) / tot if tot else float("nan") for k in range(3)]
            add("7.2", "key_ORF_body", name, CL(c), "frame_shares_%s(F0/F1/F2)" % cellname, float("nan"), norm=nm, note="%s: %.3f / %.3f / %.3f" % ((cellname,) + tuple(fs)))

        xb, yb = binsum(x_all, lo, hi, 3), binsum(y_all, lo, hi, 3)
        add("7.2", "key_ORF_body", name, CL(c), "Spearman_codon_bin3nt", spear(xb, yb), *block_boot(spear, xb, yb, 20), len(xb), int(((xb > 0) | (yb > 0)).sum()), int(((xb > 0) & (yb > 0)).sum()), "3-nt 箱求和 CPM29", "ORF 体 CHX/HAR 分布相似性；" + rnote(c))

    zlo, zhi = s0 - 15, s0 + 45
    for c in [cc for cc in COND if cc[0] == "HAR"]:
        x, y = CV[29][c][zlo:zhi + 1], CP[29][c][zlo:zhi + 1]
        pv, pp = int(np.argmax(x)) + zlo, int(np.argmax(y)) + zlo
        f1v = [CV[29][c][p] if frame(p) == 1 else 0 for p in range(zlo, zhi + 1)]; f1p = [CP[29][c][p] if frame(p) == 1 else 0 for p in range(zlo, zhi + 1)]
        tv, tp = int(np.argmax(f1v)) + zlo, int(np.argmax(f1p)) + zlo
        add("7.2", "key_ORF_HAR_start_window", name, CL(c), "Spearman_perpos_window", spear(x, y), norm="CPM29", n_coord=len(x), n_either=int(((x > 0) | (y > 0)).sum()), n_both=int(((x > 0) & (y > 0)).sum()),
            note="窗口 %d-%d（首个起点 −15 到 +45；n 小，无 bootstrap）；最大峰位置 Vero %d / PK1 %d（差 %d nt）；候选框(+1)最大峰位置 Vero %d / PK1 %d（差 %d nt）；%s" % (zlo, zhi, pv, pp, pp - pv, tv, tp, tp - tv, rnote(c)))
        for s in sorted(CAND[name][0]):
            wlo, whi = s, s + 12
            xv, yv = CV[29][c][wlo:whi + 1], CP[29][c][wlo:whi + 1]
KEY_STARTS = {"N-x1": [26325], "N-x2": [26703, 26760]}
for name, ss in KEY_STARTS.items():
    for s in ss:
        for c in [cc for cc in COND if cc[0] == "HAR"]:
            xv, yv = VERO[29][c], PK1[29][c]
            loc_v = np.mean([CV[29][c][p] for p in range(s - 30, s + 31) if frame(p) == 1 and p != s]); loc_p = np.mean([CP[29][c][p] for p in range(s - 30, s + 31) if frame(p) == 1 and p != s])
            add("7.2", "key_ORF_HAR_start_site", "%s@%d" % (name, s), CL(c), "start_peak_over_local_F1_mean", float("nan"), norm="CPM29",
                n_coord=61, note="Vero: P-sites=%d, 倍数=%.1f；PK1: P-sites=%d, 倍数=%.1f；%s" % (xv[s], CV[29][c][s] / loc_v if loc_v else float('nan'), yv[s], CP[29][c][s] / loc_p if loc_p else float('nan'), rnote(c)))


HOT = [("5UTR", 22), ("5UTR", 113), ("5UTR", 271), ("3UTR", 27707), ("3UTR", 27820), ("3UTR", 27919), ("3UTR", 27920)]
for reg, h in HOT:
    for c in COND:
        x_all, y_all = CV[29][c], CP[29][c]
        wlo, whi = h - 10, min(VLEN, h + 10)
        x, y = x_all[wlo:whi + 1], y_all[wlo:whi + 1]
        pv, pp = int(np.argmax(x)) + wlo, int(np.argmax(y)) + wlo
        bg_v = float(np.median(np.concatenate([x_all[h - 100:wlo], x_all[whi + 1:h + 101]]))); bg_p = float(np.median(np.concatenate([y_all[h - 100:wlo], y_all[whi + 1:h + 101]])))
        add("7.3", "UTR_hotspot_window_pm10", "%s@%d" % (reg, h), CL(c), "Spearman_perpos_window", spear(x, y), norm="CPM29", n_coord=len(x), n_either=int(((x > 0) | (y > 0)).sum()), n_both=int(((x > 0) & (y > 0)).sum()),
            note="窗口 ±10 nt（n=21，无 CI）；最大峰位置 Vero %d / PK1 %d（差 %d nt）；峰值 CPM Vero %.0f / PK1 %.0f；±100 nt 背景中位数 Vero %.1f / PK1 %.1f；%s" % (pv, pp, pp - pv, x.max(), y.max(), bg_v, bg_p, rnote(c)))
with open(f"{OUT}/PK1_VERO_CORRELATION.tsv", "w", encoding="utf-8", newline="") as f:
    f.write("# PK1_VERO_CORRELATION v1.0；只含 PK1 某条件 vs Vero 同条件的相关；PEDV_GIb（refdata_GIb_v1.1）坐标；29 nt 主口径，30 nt 敏感性；区域为非重叠分区（先到先得）\n")
    f.write("# 区域：5UTR 1-292；ORF1a 293-12622；ORF1b 12623-20609；S 20610-24757；ORF3_E_M 24758-26299；N 26300-27625；3UTR 27626-27959\n")
    f.write("# ci95：移动块 bootstrap（块长 300 nt，500 次；窗口较短的项见 note）；n_coordinates 为位点数（分箱项为箱数）；frac_nonzero_either = 两库至少一个非零的位点占比\n")
    f.write("# Includes the two 16-h pairs only; 8-h samples are excluded because they are not independent of the 16-h samples. MOI=1 for PK1 and Vero. PK1 CHX 29-nt periodicity is weak.\n")
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(rows)
print(len(rows), "rows")

np.savez(f"{OUT}/_plotdata.npz", **{"V%d_%s" % (L, "_".join(map(str, c))): CV[L][c] for L in (29,) for c in COND}, **{"P%d_%s" % (L, "_".join(map(str, c))): CP[L][c] for L in (29,) for c in COND})
