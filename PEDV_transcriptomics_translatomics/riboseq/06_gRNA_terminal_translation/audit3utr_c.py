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
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    exec(open(__import__('os').path.join(__import__('os').path.dirname(__REPO_FILE__), 'audit3utr_b.py'), encoding="utf-8").read())


cd2["in_frame_with_N_CDS(same_stop)"] = [(r.start_region == "N_terminal60" and r.frame_rel_N == 0) for r in cd2.itertuples()]
cd2.to_csv(T + r"/3UTR_candidate_ORFs.tsv", sep="\t", index=False, float_format="%.4g")


def region_row(name, lo, hi, sid, dr, t):
    a = S_FORMAL[sid]; seg = a[lo:hi + 1]; ra = RNA[(dr, t)][lo:hi + 1]
    return dict(region=name, library=LAB[sid], drug=dr, time_h=t, length_nt=hi - lo + 1, psites=seg.sum(), density_fold_of_library_mean=seg.mean() / MEAN_DEN[sid], CPM_per_nt=seg.mean() / CPM_DEN[sid],
                RNA_depth_fold_of_library_mean=ra.mean() / RNA_MEAN[(dr, t)], depth_class=depth_class(seg.sum()))
TREGS = {k: REG[k] for k in ["N_CDS_body", "N_stop_pm30", "UTR3_P(27626-27685)", "UTR3_Z1(27686-27776)", "UTR3_Z2(27777-27867)", "UTR3_Z3(27868-27959)", "UTR3_all", "M_body", "S_body"]}
TREGS["UTR5_all(1-292)"] = (1, 292)
tc = pd.DataFrame([region_row(n, lo, hi, sid, dr, t) for n, (lo, hi) in TREGS.items() for (dr, t), sid in LIBS.items()])
tc["RPF_over_RNA_descriptive"] = tc.density_fold_of_library_mean / tc.RNA_depth_fold_of_library_mean

pv = tc.pivot_table(index=["region", "time_h"], columns="drug", values="density_fold_of_library_mean").reset_index()
pv["HAR_over_CHX_norm_density"] = pv.HAR / pv.CHX
tc = tc.merge(pv[["region", "time_h", "HAR_over_CHX_norm_density"]], on=["region", "time_h"])
tc.insert(0, "level", "region")

cand_tc = HR.merge(CR[["candidate_id", "time_h", "body_psites", "cand_frame_frac", "evidence_class"]], on=["candidate_id", "time_h"])
cand_tc = cand_tc.rename(columns={"evidence_class": "CHX_body_evidence_class"}); cand_tc.insert(0, "level", "candidate")

man = pd.read_csv(ROOT + r"/GIb_rnaseq_reanalysis/sample_manifest_GIb.tsv", sep="\t")
man = man[man.sample_id.isin(SR.WHITELIST_LIBRARIES)]; SR.assert_exactly(man.sample_id, SR.RNA_VERO_ALIASES, "audit3utr_c RNA time course")
rr = []
for r in man.itertuples():
    try: d = rna_depth(r.sample_id)
    except Exception: continue
    m = d[1:L + 1].mean()
    for nm in ["N_CDS_body", "N_stop_pm30", "UTR3_P(27626-27685)", "UTR3_Z1(27686-27776)", "UTR3_Z2(27777-27867)", "UTR3_Z3(27868-27959)", "UTR3_all"]:
        lo, hi = REG[nm]; rr.append(dict(level="RNA_only", region=nm, treatment=r.treatment, time_h=r.time_h, RNA_depth_fold_of_library_mean=d[lo:hi + 1].mean() / m, rna_library_index=f"{r.treatment}{r.time_h}h-{len(rr)}"))
rdf = pd.DataFrame(rr)
rdf["rna_library_index"] = rdf.groupby(["region", "treatment", "time_h"]).cumcount() + 1
pd.concat([tc, cand_tc, rdf], ignore_index=True).to_csv(T + r"/3UTR_timecourse.tsv", sep="\t", index=False, float_format="%.4g")
print(tc[tc.region.isin(["N_CDS_body", "UTR3_all", "UTR5_all(1-292)"])].pivot_table(index=["region", "drug"], columns="time_h", values="density_fold_of_library_mean").round(3).to_string())
print(tc[tc.region.isin(["N_CDS_body", "UTR3_all", "UTR5_all(1-292)"])].pivot_table(index=["region", "drug"], columns="time_h", values="RPF_over_RNA_descriptive").round(3).to_string())


c5 = pd.DataFrame(enum_cands(1, 292, minaa=8)); c5.insert(0, "candidate_id", [f"U5-{i + 1:02d}" for i in range(len(c5))])
c5["start_region"] = "5UTR"
HR5, CR5 = evaluate(c5)
def peaks_in(lo, hi):
    out = []
    for (dr, t), sid in LIBS.items():
        a = S_FORMAL[sid]
        for p in range(lo, hi + 1):
            c = a[p]
            if c < 10: continue
            w = np.concatenate([a[max(p - 50, 1):p], a[p + 1:min(p + 51, L) + 1]]); med = max(np.median(w), 1)
            if c >= 5 * med:
                lens = {l: A[sid][l][p] for l in A[sid] if l <= 40 and A[sid][l][p] > 0}; tl = max(lens, key=lens.get); n = sum(lens.values())
                nb = np.mean([a[max(p - 6, 1)], a[max(p - 3, 1)], a[min(p + 3, L)], a[min(p + 6, L)]])
                out.append(dict(drug=dr, time_h=t, position=p, top_len_share=lens[tl] / n, nbr=nb / c))
    return pd.DataFrame(out)
pk5 = peaks_in(1, 292); pk3 = peaks_in(*UTR3)
def regions_of(df):
    ps = np.sort(df[df.drug == "HAR"].position.unique()); rg = []
    for p in ps:
        if rg and p - rg[-1][1] <= 6: rg[-1][1] = p
        else: rg.append([p, p])
    return rg
def cmp_row(name, lo, hi, cdf, hr, pkd):
    rg = regions_of(pkd)
    attrib_starts = hr[hr.start_attributable_HAR_response].start.unique()
    attrib_regions = {i for i, (a, b) in enumerate(rg) for s in attrib_starts if a - 6 <= s <= b + 6}
    return dict(region=name, length_nt=hi - lo + 1, seq_defined_ORFs_ge8aa_closed=int(((cdf.len_aa >= 8) & cdf.closed).sum()), of_which_AUG=int(((cdf.len_aa >= 8) & cdf.closed & (cdf.start_codon == "AUG")).sum()),
                evaluable_starts_ge1_HAR_lib=int(hr[hr.evaluable].start.nunique()), starts_with_start_attributable_HAR_response_ge1_lib=int(len(attrib_starts)),
                starts_attributable_ge2_libs=int((hr[hr.start_attributable_HAR_response].groupby("start").size() >= 2).sum()),
                distinct_HAR_peak_regions=len(rg), HAR_regions_containing_attributable_start=len(attrib_regions),
                local_peaks_unique_positions_CHX=int(pkd[pkd.drug == "CHX"].position.nunique()), local_peaks_unique_positions_HAR=int(pkd[pkd.drug == "HAR"].position.nunique()),
                frac_peaks_single_length_isolated=float(((pkd.top_len_share >= 0.6) & (pkd.nbr <= 0.2)).mean()) if len(pkd) else np.nan)
c3 = cd2[cd2.start_region == "3UTR"]; HR3 = HR[HR.candidate_id.isin(c3.candidate_id)]
cmpdf = pd.DataFrame([cmp_row("5UTR(1-292)", 1, 292, c5, HR5, pk5), cmp_row("3UTR(27626-27959)", UTR3[0], UTR3[1], c3[c3.len_aa >= 8], HR3, pk3)])
cmpdf = cmpdf.set_index("region").T.reset_index().rename(columns={"index": "metric"})
cmpdf["directly_comparable"] = ["yes"] * 3 + ["yes"] * 6 + ["yes", "yes"]
extra = pd.DataFrame([
    dict(metric="transcript context", **{"5UTR(1-292)": "leader present on gRNA and every sgRNA", "3UTR(27626-27959)": "3′ co-terminal; present on gRNA and every sgRNA"}, directly_comparable="descriptive only"),
    dict(metric="adjacent strong canonical signal", **{"5UTR(1-292)": "ORF1a AUG (293) immediately downstream", "3UTR(27626-27959)": "N stop codon (27623-27625) and N CDS immediately upstream"}, directly_comparable="NOT equivalent events"),
    dict(metric="reading-frame periodicity of regional P-sites", **{"5UTR(1-292)": "no annotated CDS frame (uORFs overlap ORF1a-frame only near 271-292)", "3UTR(27626-27959)": "no annotated CDS frame (N-frame continuation tested separately)"}, directly_comparable="not comparable (use candidate's own frame)"),
    dict(metric="first 1-9 nt of reference", **{"5UTR(1-292)": "filled from other reference, not measured; starts there flagged", "3UTR(27626-27959)": "no polyA in reference; genome end effects possible"}, directly_comparable="descriptive only")])
for c in ("5UTR(1-292)", "3UTR(27626-27959)"): cmpdf[c] = cmpdf[c].astype(object)
pd.concat([cmpdf, extra], ignore_index=True).to_csv(T + r"/5UTR_vs_3UTR_comparison.tsv", sep="\t", index=False)
print(pd.concat([cmpdf, extra], ignore_index=True).to_string())
HR5.to_csv(W + r"/HR5.tsv", sep="\t", index=False)
pickle.dump(dict(pk=pk, pk5=pk5, tc=tc, c5=c5, HR5=HR5, cmp=cmpdf), open(W + r"/audit_c.pkl", "wb"))

pp = pk[pk.region == "3UTR"].copy()
pp["zone"] = pd.cut(pp.position, [27625, 27685, 27776, 27867, 27959], labels=["P", "Z1", "Z2", "Z3"])
print(pp.groupby(["zone", "drug"], observed=True).position.nunique().unstack())
print("recurrent peak positions (>=2 same-drug libs):", pp[pp["n_same_drug_libs_present(+-1nt)"] >= 2].groupby("zone", observed=True).position.nunique().to_dict())
print("both drugs co-located same time:", int(pp["opposite_drug_same_time_present(+-1nt)"].sum()), "of", len(pp))
top = pk[pk.region.isin(["N_CDS", "3UTR"])].sort_values("psites", ascending=False).drop_duplicates("position").head(12)
print(top[["position", "region", "drug", "time_h", "psites", "top_len", "top_len_share", "codon_neighbours_over_peak", "frame_rel_N", "single_length_isolated(D_like)", "HAR_over_CHX_norm_at_site"]].round(2).to_string())
