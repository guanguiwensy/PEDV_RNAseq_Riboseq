# -*- coding: utf-8 -*-
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
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    exec(open(__import__('os').path.join(__import__('os').path.dirname(__REPO_FILE__), 'audit3utr.py'), encoding="utf-8").read())

# ---------- 已知翻译 ORF 的 frame0 中位数（逐库，判据 θ 参考）----------
KNOWN = {"S": (20655, 24713, 20610), "E": (25406, 25555, 25361), "M": (25653, 26243, 25608), "N": (26345, 27562, 26300)}
def frame0_of(a, lo, hi, ref):
    idx = np.arange(lo, hi + 1); seg = a[lo:hi + 1]; return seg[(idx - ref) % 3 == 0].sum() / seg.sum()
KF0 = {sid: float(np.median([frame0_of(S_FORMAL[sid], *v) for v in KNOWN.values()])) for sid in S_FORMAL}
print("known-ORF frame0 medians", {LAB[k]: round(v, 3) for k, v in KF0.items()})

# ---------- 候选逐库评价 ----------
def evaluate(cd_df, lo_note=""):
    hr, cr = [], []
    for c in cd_df.itertuples():
        s, e = c.start, c.stop_end
        for t in TIMES:
            sid = LIBS[("HAR", t)]; a = S_FORMAL[sid]; cs = S_FORMAL[LIBS[("CHX", t)]]
            win = a[s - 1:s + 2].sum(); local = a[max(s - 50, 1):min(s + 50, L) + 1].sum()
            flank = np.concatenate([a[max(s - 50, 1):s - 3], a[s + 4:min(s + 50, L) + 1]]); bg = flank.mean() if len(flank) else np.nan
            swe = (win / 3) / bg if bg and bg > 0 else (np.inf if win > 0 else np.nan)
            seg = a[max(s - 6, 1):s + 7]; am = int(np.argmax(seg)) + max(s - 6, 1) - s if seg.max() > 0 else None
            cwin = cs[s - 1:s + 2].sum(); cflank = np.concatenate([cs[max(s - 50, 1):s - 3], cs[s + 4:min(s + 50, L) + 1]]); cbg = cflank.mean() if len(cflank) else np.nan
            har_n = win / (3 * MEAN_DEN[sid]); chx_n = cwin / (3 * MEAN_DEN[LIBS[("CHX", t)]])
            flank_ratio = (bg / MEAN_DEN[sid]) / (cbg / MEAN_DEN[LIBS[("CHX", t)]]) if cbg and cbg > 0 and bg else np.nan
            shc = (har_n / chx_n) / flank_ratio if chx_n > 0 and flank_ratio and not np.isnan(flank_ratio) else np.nan
            evaluable = local >= 50; attrib = bool(evaluable and swe >= 5 and win >= 10 and am is not None and abs(am) <= 1)
            lens = {l: int(rng(A[sid][l], s - 1, s + 1).sum()) for l in HOST_SET[sid] if l in A[sid]}
            hr.append(dict(candidate_id=c.candidate_id, start=s, start_codon=c.start_codon, library=LAB[sid], time_h=t,
                           HAR_psites_start_pm1=win, HAR_local_pm50_total=local, local_background_per_nt=bg, SWE=swe, SHC_start_window_HAR_CHX_contrast=shc,
                           HAR_peak_rel_start_within_pm6=am, CHX_psites_start_pm1=cwin, evaluable=evaluable, start_attributable_HAR_response=attrib,
                           read_length_support=";".join(f"{l}:{n}" for l, n in lens.items() if n), n_lengths_supporting=sum(1 for n in lens.values() if n > 0)))
            sidc = LIBS[("CHX", t)]; ac = S_FORMAL[sidc]
            body = ac[s + 3:e - 2]; n = body.sum()
            idx = np.arange(s + 3, e - 2); inf = body[(idx - s) % 3 == 0].sum()
            b1 = body.copy()
            if len(b1): b1[np.argmax(b1)] = 0
            idx1 = idx; inf1 = b1[(idx1 - s) % 3 == 0].sum()
            b3 = body.copy()
            if len(b3): b3[np.argsort(b3)[::-1][:3]] = 0
            inf3 = b3[(idx - s) % 3 == 0].sum()
            fr = [body[(idx - s) % 3 == k].sum() for k in range(3)] if len(body) else [0, 0, 0]
            th = 0.8 * KF0[sidc]
            if n < 30: cl = "Not assessable"
            elif n < 100: cl = "Low depth"
            else:
                ok = (inf / n >= th) and (inf / n >= max(fr[1], fr[2]) / n) and (b1.sum() > 0 and inf1 / b1.sum() >= th)
                cl = "Measurable; candidate-associated periodicity" if ok else "Measurable; no convincing periodicity"
            cr.append(dict(candidate_id=c.candidate_id, start=s, library=LAB[sidc], time_h=t, body_nt=len(body), body_psites=n, body_psites_per_codon=n / max((e - s + 1) // 3, 1),
                           cand_frame_frac=inf / n if n else np.nan, other_frame1_frac=fr[1] / n if n else np.nan, other_frame2_frac=fr[2] / n if n else np.nan,
                           cand_frame_frac_after_top1_removed=inf1 / b1.sum() if b1.sum() else np.nan, cand_frame_frac_after_top3_removed=inf3 / b3.sum() if b3.sum() else np.nan, theta_0p8x_known_ORF_frame0=th, evidence_class=cl))
    return pd.DataFrame(hr), pd.DataFrame(cr)
HR, CR = evaluate(cd)
# 起始 → HAR 响应区聚类（相距 ≤6 nt 的 HAR 局部峰合并，区间 27,566-27,959）
harpk = pk[(pk.drug == "HAR") & (pk.position >= 27566)].position.sort_values().unique()
regions = []
for p in harpk:
    if regions and p - regions[-1][1] <= 6: regions[-1][1] = p
    else: regions.append([p, p])
def region_id(s):
    for i, (a, b) in enumerate(regions):
        if a - 6 <= s <= b + 6: return f"R{i + 1:02d}({a}-{b})"
    return ""
HR["HAR_response_region"] = HR.start.map(region_id)
HR.to_csv(T + r"/3UTR_HAR_candidate_responses.tsv", sep="\t", index=False, float_format="%.4g")
CR.to_csv(T + r"/3UTR_candidate_CHX_body_evidence.tsv", sep="\t", index=False, float_format="%.4g")
sumc = HR.groupby("candidate_id").agg(n_HAR_evaluable_libs=("evaluable", "sum"), n_start_attributable_libs=("start_attributable_HAR_response", "sum"), best_SWE=("SWE", "max"),
                                      HAR_response_region=("HAR_response_region", "first")).reset_index()
sumc = sumc.merge(CR.groupby("candidate_id").agg(n_CHX_assessable_libs=("evidence_class", lambda x: int((x != "Not assessable").sum())),
                                                  n_CHX_candidate_periodicity_libs=("evidence_class", lambda x: int((x == "Measurable; candidate-associated periodicity").sum())),
                                                  max_body_psites=("body_psites", "max")).reset_index(), on="candidate_id")
cd2 = cd.merge(sumc, on="candidate_id")
cd2["rank_key"] = list(zip(-cd2.n_start_attributable_libs, -cd2.n_CHX_candidate_periodicity_libs, -cd2.best_SWE.replace(np.inf, 1e9).fillna(0)))
cd2["rank_by_(attributable_HAR_timepoints,CHX_periodic_timepoints,best_SWE)"] = cd2.rank_key.rank(method="min").astype(int)
cd2 = cd2.drop(columns="rank_key").sort_values("rank_by_(attributable_HAR_timepoints,CHX_periodic_timepoints,best_SWE)")
# 与先前 3′ 末端目录交叉
prev = pd.read_csv(ROOT + r"/fig5_terminal_ORF/PEDV_terminal_ORF_candidates.tsv", sep="\t", comment="#")
print(prev.columns.tolist()[:25])
cd2.to_csv(T + r"/3UTR_candidate_ORFs.tsv", sep="\t", index=False, float_format="%.4g")
print("HAR response regions", len(regions), regions)
print(HR.groupby("library")[["evaluable", "start_attributable_HAR_response"]].sum())
print("candidates with >=1 attributable lib:", int((cd2.n_start_attributable_libs > 0).sum()), "of", len(cd2), "; with >=2:", int((cd2.n_start_attributable_libs >= 2).sum()))
print("CHX candidate periodicity libs>=1:", int((cd2.n_CHX_candidate_periodicity_libs > 0).sum()), ">=2:", int((cd2.n_CHX_candidate_periodicity_libs >= 2).sum()))
print(cd2.head(12)[["candidate_id", "start", "start_codon", "len_aa", "start_region", "n_HAR_evaluable_libs", "n_start_attributable_libs", "best_SWE", "n_CHX_candidate_periodicity_libs", "max_body_psites", "HAR_response_region"]].to_string())
pickle.dump(dict(cd2=cd2, HR=HR, CR=CR, regions=regions, KF0=KF0), open(W + r"/audit_b.pkl", "wb"))
