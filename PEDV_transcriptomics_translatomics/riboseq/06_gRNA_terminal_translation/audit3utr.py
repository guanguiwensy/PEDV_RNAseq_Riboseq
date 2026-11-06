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
import os, sys, gzip, pickle, json, itertools, collections
import numpy as np, pandas as pd
ROOT = PEDV_P+r""
OUT = ROOT + r"/Fig5_3UTR_translation_audit"; T = OUT + r"/tables"; W = OUT + r"/_work"
L = 27959
SEQ = "".join(l.strip() for l in open(ROOT + r"/refbuild_GIb/virus/PEDV_GIb.fa") if not l.startswith(">")).upper().replace("U", "T")
assert len(SEQ) == L
cache = pickle.load(open(ROOT + r"/SupM_alt_ORF_evidence/cache/psite_by_length.pkl", "rb"))
A = cache["A"]
LIBS = {("CHX", 5): SR.lib("RIBO_VERO_CHX_1MOI_5h"), ("CHX", 8): SR.lib("RIBO_VERO_CHX_1MOI_8h"), ("CHX", 16): SR.lib("RIBO_VERO_CHX_1MOI_16h"), ("HAR", 5): SR.lib("RIBO_VERO_HAR_1MOI_5h"), ("HAR", 8): SR.lib("RIBO_VERO_HAR_1MOI_8h"), ("HAR", 16): SR.lib("RIBO_VERO_HAR_1MOI_16h")}
LAB = {v: f"{k[0]} {k[1]}h" for k, v in LIBS.items()}
HOST_SET = {SR.lib("RIBO_VERO_CHX_1MOI_5h"): [28, 29, 30], SR.lib("RIBO_VERO_CHX_1MOI_8h"): [26, 27, 28, 29, 30, 31], SR.lib("RIBO_VERO_CHX_1MOI_16h"): [26, 28, 29, 30, 31], SR.lib("RIBO_VERO_HAR_1MOI_5h"): [29, 30], SR.lib("RIBO_VERO_HAR_1MOI_8h"): [29, 30, 32], SR.lib("RIBO_VERO_HAR_1MOI_16h"): [28, 29, 30, 31]}
STRICT = [29, 30]
RNA_PAIR = {(t, h): SR.lib(SR.RNA_PARTNER[SR.ribo_vero(t, h)]) for t in ("CHX", "HAR") for h in (5, 8, 16)}   # RNA-seq partner (Vero, MOI 1, same drug and time) from the registry
TIMES = (5, 8, 16)
NSTART, NSTOP0, NSTOP1 = 26300, 27623, 27625
UTR3 = (27626, 27959)
REG = {"N_CDS_body": (26345, 27562), "N_last300": (27326, 27625), "N_stop_pm30": (27593, 27655), "N_stop_pm60": (27563, 27685), "N_stop_pm120": (27503, 27745),
       "UTR3_P(27626-27685)": (27626, 27685), "UTR3_Z1(27686-27776)": (27686, 27776), "UTR3_Z2(27777-27867)": (27777, 27867), "UTR3_Z3(27868-27959)": (27868, 27959),
       "UTR3_all": UTR3, "M_body": (25653, 26243), "S_body": (20655, 24713)}
UTR_ZONES = ["UTR3_P(27626-27685)", "UTR3_Z1(27686-27776)", "UTR3_Z2(27777-27867)", "UTR3_Z3(27868-27959)"]
CAND_REGS = ["N_CDS_body", "N_stop_pm30", "N_stop_pm60"] + UTR_ZONES + ["UTR3_all", "M_body", "S_body"]
OWN_FRAME = {"M_body": 25608, "S_body": 20610}

def arr_set(sid, lens):
    a = np.zeros(L + 2)
    for l in lens:
        if l in A[sid]: a += A[sid][l]
    return a
S_FORMAL = {sid: arr_set(sid, HOST_SET[sid]) for sid in HOST_SET}
S_STRICT = {sid: arr_set(sid, STRICT) for sid in HOST_SET}
def tot(a): return a[1:L + 1].sum()
MEAN_DEN = {sid: tot(S_FORMAL[sid]) / L for sid in S_FORMAL}
CPM_DEN = {sid: tot(S_FORMAL[sid]) / 1e6 for sid in S_FORMAL}

def rna_depth(r, kind="raw"):
    p = ROOT + (r"/GIb_rnaseq_reanalysis/figdata/fig1_depth" if kind == "raw" else r"/GIb_rnaseq_reanalysis/coverage")
    f = f"{p}/{r}.raw.depth.tsv" if kind == "raw" else f"{p}/{r}.depth.tsv"
    if kind == "raw" and not os.path.exists(f): f = ROOT + rf"/GIb_rnaseq_reanalysis/figdata/fig1_depth_extra/{r}.raw.depth.tsv"
    d = pd.read_csv(f, sep="\t", header=None, names=["c", "p", "d"]); a = np.zeros(L + 2); a[d.p.values] = d.d.values; return a
RNA = {k: rna_depth(r) for k, r in RNA_PAIR.items()}
RNA_MEAN = {k: a[1:L + 1].mean() for k, a in RNA.items()}

def rng(a, lo, hi): return a[max(lo, 1):min(hi, L) + 1]
def depth_class(n): return "Not assessable" if n < 30 else ("Low depth" if n < 100 else "Measurable")
def jsd(p, q):
    p = np.asarray(p, float); q = np.asarray(q, float); p = p / p.sum(); q = q / q.sum(); m = (p + q) / 2
    f = lambda a, b: np.sum(np.where(a > 0, a * np.log2(a / np.where(b > 0, b, 1)), 0))
    return 0.5 * f(p, m) + 0.5 * f(q, m)


LO, HI = 27326, 27959
rows = []
for (dr, t), sid in LIBS.items():
    fa, sa = S_FORMAL[sid], S_STRICT[sid]; ra = RNA[(dr, t)]
    for p in range(LO, HI + 1):
        bl = {l: A[sid][l][p] for l in A[sid] if 24 <= l <= 36 and A[sid][l][p] > 0}
        reg = ("N_CDS" if p <= 27620 else "N_stop_codon" if p <= NSTOP1 else "3UTR")
        rows.append(dict(position=p, region=reg, drug=dr, time_h=t, psite_formal_set=fa[p], psite_strict_core=sa[p],
                         density_fold_of_library_mean=fa[p] / MEAN_DEN[sid], CPM=fa[p] / CPM_DEN[sid],
                         psite_by_length=";".join(f"{l}:{int(v)}" for l, v in sorted(bl.items())), frame_rel_N=(p - NSTART) % 3,
                         mapping_quality="unique (STAR NH=1; MAPQ 255)", RNAseq_raw_depth=ra[p], RNAseq_depth_fold_of_library_mean=ra[p] / RNA_MEAN[(dr, t)]))
pos = pd.DataFrame(rows); pos.to_csv(T + r"/3UTR_position_level_Psites.tsv", sep="\t", index=False, float_format="%.5g")


hostlen = pd.read_csv(ROOT + r"/shared_riboseq_metrics/length_metrics_long.tsv", sep="\t", comment="#")
rowsq = []; LD = {}
LENS = list(range(24, 37))
for (dr, t), sid in LIBS.items():
    for rn in CAND_REGS:
        lo, hi = REG[rn]
        lens = {l: rng(A[sid][l], lo, hi).sum() if l in A[sid] else 0 for l in range(16, 41)}
        n_all = sum(lens.values()); LD[(sid, rn)] = np.array([lens[l] for l in LENS], float)
        arr = S_FORMAL[sid]; seg = rng(arr, lo, hi); n = seg.sum()
        ref = NSTART if rn not in OWN_FRAME else OWN_FRAME[rn]
        idx = np.arange(max(lo, 1), min(hi, L) + 1); fr = (idx - ref) % 3
        fc = np.array([seg[fr == k].sum() for k in range(3)])
        fcs = np.sort(seg)[::-1]
        seg3 = seg.copy(); seg3[np.argsort(seg)[::-1][:3]] = 0; seg1 = seg.copy(); seg1[np.argmax(seg)] = 0
        fc3 = np.array([seg3[fr == k].sum() for k in range(3)]); fc1 = np.array([seg1[fr == k].sum() for k in range(3)])
        med = np.median(np.repeat(list(lens), [lens[l] for l in lens])) if n_all else np.nan
        rowsq.append(dict(region=rn, library=LAB[sid], psites_formal_set=n, depth_class=depth_class(n),
                          frame0=fc[0] / n if n else np.nan, frame1=fc[1] / n if n else np.nan, frame2=fc[2] / n if n else np.nan,
                          frame_max=fc.max() / n if n else np.nan, frame_max_idx=int(fc.argmax()) if n else -1,
                          frame_max_after_top1_removed=fc1.max() / fc1.sum() if fc1.sum() else np.nan, frame_max_after_top3_removed=fc3.max() / fc3.sum() if fc3.sum() else np.nan,
                          top1_share=fcs[0] / n if n else np.nan, top3_share=fcs[:3].sum() / n if n else np.nan,
                          all_length_reads=n_all, median_len=med, frac_28_31=sum(lens[l] for l in range(28, 32)) / n_all if n_all else np.nan,
                          top_len=max(lens, key=lens.get) if n_all else np.nan, top_len_share=max(lens.values()) / n_all if n_all else np.nan))
qc = pd.DataFrame(rowsq)
ref_ld = {sid: LD[(sid, "N_CDS_body")] for sid in S_FORMAL}
qc["JSD_len_vs_N_CDS_body"] = [jsd(LD[(sid, r.region)], ref_ld[sid]) if LD[(sid, r.region)].sum() and ref_ld[sid].sum() else np.nan
                              for sid, r in zip([ {v: k for k, v in LAB.items()}[x] for x in qc.library], qc.itertuples())]
def conv(g):
    ok = (g.depth_class == "Measurable") & (g.frame_max >= 0.5) & (g.frame_max_after_top3_removed >= 0.5)
    return ok
qc["periodicity_candidate_flag(single_lib)"] = conv(qc)
def final_class(r):
    if r.depth_class == "Not assessable": return "Not assessable"
    if r.depth_class == "Low depth": return "Low depth"
    return "Measurable; periodic-like (single library)" if r["periodicity_candidate_flag(single_lib)"] else "Measurable; no convincing periodicity"
qc["evidence_class"] = qc.apply(final_class, axis=1)

cons = []
for rn in CAND_REGS:
    for dr in ("CHX", "HAR"):
        g = qc[(qc.region == rn) & (qc.library.str.startswith(dr)) & qc["periodicity_candidate_flag(single_lib)"]]
        cons.append((rn, dr, g.frame_max_idx.value_counts().max() if len(g) else 0))
cdf = {(a, b): c for a, b, c in cons}
qc["cross_library_consistent_(>=2_same_drug_same_frame)"] = [cdf[(r.region, r.library.split()[0])] >= 2 for r in qc.itertuples()]

hr = hostlen[(hostlen.domain == "host") if "domain" in hostlen.columns else slice(None)]
qc.to_csv(T + r"/3UTR_RPF_length_frame_QC.tsv", sep="\t", index=False, float_format="%.4g")
pd.DataFrame(LD.keys(), columns=["sample", "region"]).assign(**{f"len{l}": [LD[k][i] for k in LD] for i, l in enumerate(LENS)}).assign(
    library=lambda d: d["sample"].map(LAB)).drop(columns="sample").to_csv(T + r"/3UTR_RPF_length_distribution_counts.tsv", sep="\t", index=False)
print(qc[qc.region.isin(["N_CDS_body", "N_stop_pm30", "UTR3_all"] + UTR_ZONES)][["region", "library", "psites_formal_set", "depth_class", "frame0", "frame1", "frame2", "frame_max_after_top3_removed", "top_len_share", "evidence_class"]].round(2).to_string())


pk = []
for (dr, t), sid in LIBS.items():
    a = S_FORMAL[sid]
    for p in range(LO, HI + 1):
        c = a[p]
        if c < 10: continue
        w = np.concatenate([a[max(p - 50, 1):p], a[p + 1:min(p + 51, L) + 1]]); med = max(np.median(w), 1)
        if c >= 5 * med:
            nb = np.mean([a[p - 6], a[p - 3], a[p + 3], a[p + 6]])
            lens = {l: A[sid][l][p] for l in A[sid] if l <= 40 and A[sid][l][p] > 0}; tl = max(lens, key=lens.get); n = sum(lens.values())
            pk.append(dict(drug=dr, time_h=t, position=p, region=("N_CDS" if p <= 27622 else "N_stop_codon" if p <= NSTOP1 else "3UTR"), psites=c, fold_over_local_median=c / med,
                           density_fold_of_library_mean=c / MEAN_DEN[sid], top_len=tl, top_len_share=lens[tl] / n, median_len=float(np.median(np.repeat(list(lens), [lens[l] for l in lens]))),
                           codon_neighbours_over_peak=nb / c, frame_rel_N=(p - NSTART) % 3, RNA_fold=RNA[(dr, t)][p] / RNA_MEAN[(dr, t)]))
pk = pd.DataFrame(pk)
pk["single_length_isolated(D_like)"] = (pk.top_len_share >= 0.6) & (pk.codon_neighbours_over_peak <= 0.2)

def near(df, p): return df[(df.position >= p - 1) & (df.position <= p + 1)]
pk["n_same_drug_libs_present(+-1nt)"] = [near(pk[pk.drug == r.drug], r.position).time_h.nunique() for r in pk.itertuples()]
pk["opposite_drug_same_time_present(+-1nt)"] = [len(near(pk[(pk.drug != r.drug) & (pk.time_h == r.time_h)], r.position)) > 0 for r in pk.itertuples()]

def norm_at(dr, t, p): return S_FORMAL[LIBS[(dr, t)]][p] / MEAN_DEN[LIBS[(dr, t)]]
pk["HAR_over_CHX_norm_at_site"] = [(norm_at("HAR", r.time_h, r.position) + 1e-9) / (norm_at("CHX", r.time_h, r.position) + 1e-9) for r in pk.itertuples()]
pk.to_csv(T + r"/3UTR_local_peaks.tsv", sep="\t", index=False, float_format="%.4g")
print("peaks total", len(pk), pk.groupby(["region", "drug"]).size().to_dict())
print("3UTR peaks", pk[pk.region == "3UTR"].groupby("drug").position.nunique().to_dict(), "recurrent(>=2 same-drug libs)", int(((pk.region == "3UTR") & (pk["n_same_drug_libs_present(+-1nt)"] >= 2)).sum()),
      "D_like frac", round(pk[pk.region == "3UTR"]["single_length_isolated(D_like)"].mean(), 2))


def prof(a, lo, hi): return np.array([a[p] for p in range(lo, hi + 1)])
rows = []
for (dr, t), sid in LIBS.items():
    a = S_FORMAL[sid]; nb = rng(a, REG["N_CDS_body"][0], REG["N_CDS_body"][1]).mean()
    for nm, w in (("pm30", 30), ("pm60", 60), ("pm120", 120)):
        up = rng(a, NSTOP0 - w, NSTOP0 - 1); dn = rng(a, NSTOP1 + 1, NSTOP1 + w)
        def fr(seg, lo):
            idx = np.arange(lo, lo + len(seg)); f = (idx - NSTART) % 3; return [seg[f == k].sum() / seg.sum() if seg.sum() else np.nan for k in range(3)]
        fu, fd = fr(up, NSTOP0 - w), fr(dn, NSTOP1 + 1)
        rows.append(dict(library=LAB[sid], window=nm, up_psites=up.sum(), down_psites=dn.sum(), up_density_over_Nbody=up.mean() / nb, down_density_over_Nbody=dn.mean() / nb,
                         up_frame0=fu[0], up_frame1=fu[1], up_frame2=fu[2], down_frame0=fd[0], down_frame1=fd[1], down_frame2=fd[2],
                         down_depth_class=depth_class(dn.sum()), down_frame0_after_top3_removed=(lambda s: (s[(np.arange(len(s)) % 3 == 0)].sum() / s.sum()) if s.sum() else np.nan)(np.where(np.isin(np.arange(len(dn)), np.argsort(dn)[::-1][:3]), 0, dn))))

    pf = a[NSTOP0 - 60:NSTOP0 + 180 + 1]
    for b0 in range(-60, 120, 10):
        seg = a[NSTOP0 + b0:NSTOP0 + b0 + 10]
        rows.append(dict(library=LAB[sid], window=f"bin{b0:+d}..{b0 + 9:+d}", up_psites=np.nan, down_psites=seg.sum(), up_density_over_Nbody=np.nan, down_density_over_Nbody=seg.mean() / nb))
nst = pd.DataFrame(rows); nst.to_csv(T + r"/3UTR_N_stop_audit.tsv", sep="\t", index=False, float_format="%.4g")

pr = []
for (dr, t), sid in LIBS.items():
    a = S_FORMAL[sid]; w = a[NSTOP0 - 30:NSTOP0 + 31]; k = int(w.argmax()) - 30
    pr.append(dict(library=LAB[sid], argmax_rel_stop_first_nt=k, max_psites=w.max(), psites_at_minus3=a[NSTOP0 - 3], psites_at_minus6=a[NSTOP0 - 6], psites_at_0=a[NSTOP0], psites_at_plus3=a[NSTOP0 + 3],
                   total_pm30_up=a[NSTOP0 - 30:NSTOP0].sum(), total_pm30_down=a[NSTOP1 + 1:NSTOP1 + 31].sum()))
pd.DataFrame(pr).to_csv(T + r"/3UTR_N_stop_peak_positions.tsv", sep="\t", index=False, float_format="%.4g")
print(nst[nst.window.isin(["pm30", "pm60"])].round(3).to_string()); print(pd.DataFrame(pr).to_string())


STARTS = {"ATG": "AUG", "CTG": "CUG", "GTG": "GUG", "TTG": "UUG", "ATA": "AUA", "ATT": "AUU", "ATC": "AUC", "ACG": "ACG"}
STOPS = {"TAA", "TAG", "TGA"}
def enum_cands(lo, hi, minaa=4):
    out = []
    for s in range(lo, hi + 1):
        c = SEQ[s - 1:s + 2]
        if c not in STARTS: continue
        stop = None
        for k in range(s + 3, L - 1, 3):
            if SEQ[k - 1:k + 2] in STOPS: stop = k + 2; break
        closed = stop is not None; e = stop if closed else L
        n_aa = ((e - s + 1) // 3 - (1 if closed else 0))
        if n_aa >= minaa: out.append(dict(start=s, start_codon=STARTS[c], stop_end=e, closed=closed, len_aa=n_aa, frame_genome=(s - 1) % 3, frame_rel_N=(s - NSTART) % 3))
    return out
def strat(n): return "8-14" if n < 15 else "15-24" if n < 25 else "25-39" if n < 40 else ">=40"
cands = []
for reg, (lo, hi) in (("3UTR", UTR3), ("N_terminal60", (27566, 27625))):
    for d in enum_cands(lo, hi): d["start_region"] = reg; cands.append(d)
cd = pd.DataFrame(cands).sort_values(["start", "stop_end"]).reset_index(drop=True)
cd.insert(0, "candidate_id", [f"U3-{i + 1:02d}" for i in range(len(cd))])
cd["len_stratum_aa"] = [strat(n) if n >= 8 else "4-7 (below Fig5 cutoff)" for n in cd.len_aa]
cd["dist_from_N_stop_nt"] = cd.start - NSTOP1
cd["dist_from_genome_3p_end_nt"] = L - cd.start
cd["codon_seq"] = [SEQ[s - 1:s + 2] for s in cd.start]
ov = []
for r in cd.itertuples():
    o = cd[(cd.candidate_id != r.candidate_id) & (cd.start <= r.stop_end) & (cd.stop_end >= r.start)]
    ov.append(";".join(o.candidate_id.tolist()))
cd["overlaps_candidates"] = ov
cd["shares_stop_with"] = [";".join(cd[(cd.stop_end == r.stop_end) & (cd.candidate_id != r.candidate_id)].candidate_id) for r in cd.itertuples()]
cd["in_Fig5_previous_terminal_catalogue"] = False
prev = pd.read_csv(ROOT + r"/fig5_terminal_ORF/PEDV_terminal_ORF_candidates.tsv", sep="\t")
print(prev.columns.tolist()[:30])
cd.to_csv(T + r"/3UTR_candidate_ORFs.tsv", sep="\t", index=False)
print("candidates", len(cd), cd.groupby(["start_region", "len_stratum_aa"]).size().to_dict(), "unclosed", int((~cd.closed).sum()))

pickle.dump(dict(cd=cd), open(W + r"/cands.pkl", "wb"))
