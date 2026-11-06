#!/usr/bin/env python3
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
import os
from collections import defaultdict
import numpy as np
import pandas as pd

D = PEDV_P+"Fig2_Supp_junction_dynamics/data/"
GL = 27959
TIMES = [3, 5, 8, 12, 16]
CAN = ["S", "ORF3", "E", "M", "N"]
NCC = ["TRSL_dep_noncanonical", "TRSL_indep_distant", "TRSL_indep_local"]
DET_MIN = 3
BIN = 20
M = 1e6
rng = np.random.default_rng(20261010)

ps = pd.read_csv(D + "junction_event_sample_matrix_long.csv.gz")
cls = pd.read_csv(D + "per_sample_class_counts.csv")

den_path = D + "viral_fragments_mm2.tsv"
if os.path.exists(den_path):
    dn = pd.read_csv(den_path, sep="\t")
    cls = cls.merge(dn[["sample_id", "viral_frag_primary_q20"]].rename(columns={"viral_frag_primary_q20": "viral_fragments"}), on="sample_id")
    DEN_NOTE = "minimap2 primary alignments on PEDV_GIb, MAPQ>=20, distinct read names (fragments)"
else:
    st = pd.read_csv(PEDV_P+"GIb_rnaseq_reanalysis/figdata/fig1_counts.tsv", sep="\t",
                     names=["sample_id", "ver", "total_mapped", "pedv_mapped", "r1f", "r1r"])
    st = st[st.ver == "raw"][["sample_id", "pedv_mapped"]].rename(columns={"pedv_mapped": "viral_fragments"})
    cls = cls.merge(st, on="sample_id")
    DEN_NOTE = "STAR unique (MAPQ=255) PEDV fragments, non-deduplicated [fallback: server minimap2 count unavailable]"
cls["nc_total"] = cls[NCC].sum(axis=1)
cls = cls.sort_values(["MOI", "hpi", "treatment"]).reset_index(drop=True)


def pooled(df, cols):
    return df.groupby(["MOI", "hpi"])[cols + ["viral_fragments"]].sum().reset_index()


# ================= Panel A =================
A = cls[["sample_id", "MOI", "treatment", "hpi", "all_classified", "canonical_total", "nc_total"] + NCC + ["viral_fragments"]].copy()
A["f_canonical"] = A.canonical_total / A.all_classified
for c in NCC:
    A["f_" + c] = A[c] / A.all_classified
    A["nc_frac_" + c] = A[c] / A.nc_total
A.to_csv(D + "panelA_per_sample_class_fractions.csv", index=False)
pA = pooled(cls, ["all_classified", "canonical_total", "nc_total"] + NCC)
pA["f_canonical"] = pA.canonical_total / pA.all_classified
for c in NCC:
    pA["f_" + c] = pA[c] / pA.all_classified
    pA["nc_frac_" + c] = pA[c] / pA.nc_total
pA.to_csv(D + "panelA_pooled_CHX_HAR_by_MOI_hpi_DESCRIPTIVE.csv", index=False)
rows = []
for (moi, t), g in A.groupby(["MOI", "hpi"]):
    c = g[g.treatment == "CHX"].iloc[0]; h = g[g.treatment == "HAR"].iloc[0]
    for k, lab in [("f_canonical", "canonical")] + [("f_" + n, n) for n in NCC]:
        p1, p2 = c[k], h[k]; n1, n2 = c.all_classified, h.all_classified
        se = np.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
        rows.append(dict(MOI=moi, hpi=t, cls=lab, CHX=p1, HAR=p2, diff_pp=100 * (p1 - p2), binomial_SE_pp=100 * se,
                         diff_over_SE=(p1 - p2) / se if se > 0 else np.nan))
cons = pd.DataFrame(rows)
cons.to_csv(D + "panelA_CHX_vs_HAR_consistency.csv", index=False)

# ================= Panel B =================
Bd = cls[["sample_id", "MOI", "treatment", "hpi", "canonical_total", "viral_fragments"] + [f"canonical_{g}" for g in CAN]].copy()
for g in CAN:
    Bd["frac_" + g] = Bd[f"canonical_{g}"] / Bd.canonical_total
    Bd["per_M_viral_frag_" + g] = Bd[f"canonical_{g}"] / (Bd.viral_fragments / M)
assert np.allclose(Bd[[f"frac_{g}" for g in CAN]].sum(axis=1), 1)
Bd.to_csv(D + "panelB_canonical_sgRNA_per_sample.csv", index=False)
pB = pooled(cls, ["canonical_total"] + [f"canonical_{g}" for g in CAN])
for g in CAN:
    pB["frac_" + g] = pB[f"canonical_{g}"] / pB.canonical_total
    pB["per_M_viral_frag_" + g] = pB[f"canonical_{g}"] / (pB.viral_fragments / M)
pB.to_csv(D + "panelB_canonical_sgRNA_pooled_DESCRIPTIVE.csv", index=False)
zero_B = int((Bd[[f"canonical_{g}" for g in CAN]] == 0).sum().sum())

# ================= Panel C =================
Cd = cls[["sample_id", "MOI", "treatment", "hpi", "all_classified", "viral_fragments"] + NCC].copy()
for n in NCC:
    Cd["frac_all_" + n] = Cd[n] / Cd.all_classified
    Cd["per_M_viral_frag_" + n] = Cd[n] / (Cd.viral_fragments / M)
Cd.to_csv(D + "panelC_noncanonical_per_sample.csv", index=False)
pC = pooled(cls, ["all_classified"] + NCC)
for n in NCC:
    pC["frac_all_" + n] = pC[n] / pC.all_classified
    pC["per_M_viral_frag_" + n] = pC[n] / (pC.viral_fragments / M)
pC.to_csv(D + "panelC_noncanonical_pooled_DESCRIPTIVE.csv", index=False)
zero_C = int((Cd[NCC] == 0).sum().sum())


nc = ps[ps.junction_class != "canonical"].copy()
ev = nc.groupby(["donor", "acc", "junction_class"]).reads.sum().reset_index()


def chain(vals, w=5):
    vals = np.sort(np.unique(vals)); lab = {}; c = 0
    for i, v in enumerate(vals):
        if i and v - vals[i - 1] > w:
            c += 1
        lab[v] = c
    return lab


dep = ev[ev.junction_class == "TRSL_dep_noncanonical"].copy()
dep["cluster_key"] = "dep_" + dep.acc.map(chain(dep.acc.values)).astype(str)
ind = ev[ev.junction_class != "TRSL_dep_noncanonical"].copy().reset_index(drop=True)
par = np.arange(len(ind))


def find(i):
    while par[i] != i:
        par[i] = par[par[i]]
        i = par[i]
    return i


grid = defaultdict(list)
for i, r in enumerate(ind.itertuples()):
    grid[(r.junction_class, r.donor // 6, r.acc // 6)].append(i)
for i, r in enumerate(ind.itertuples()):
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for j in grid.get((r.junction_class, r.donor // 6 + dx, r.acc // 6 + dy), []):
                if j > i and abs(ind.donor[j] - r.donor) <= 5 and abs(ind.acc[j] - r.acc) <= 5:
                    a, b = find(i), find(j)
                    if a != b:
                        par[b] = a
ind["cluster_key"] = ["ind_" + str(find(i)) for i in range(len(ind))]
ev2 = pd.concat([dep, ind], ignore_index=True)
top = ev2.sort_values("reads", ascending=False).drop_duplicates("cluster_key")[["cluster_key", "donor", "acc", "junction_class"]]
top = top.rename(columns={"donor": "donor_rep", "acc": "acc_rep"})
agg = ev2.groupby("cluster_key").agg(reads_total=("reads", "sum"), n_exact_events=("reads", "size"),
                                    acc_min=("acc", "min"), acc_max=("acc", "max"),
                                    donor_min=("donor", "min"), donor_max=("donor", "max")).reset_index()
CL = top.merge(agg, on="cluster_key")
CL["cluster_id"] = CL.donor_rep.astype(str) + "|" + CL.acc_rep.astype(str) + "|+|" + CL.junction_class
CL["span_nt"] = CL.acc_rep - CL.donor_rep - 1
mapk = ev2[["donor", "acc", "junction_class", "cluster_key"]].merge(CL[["cluster_key", "cluster_id"]], on="cluster_key")
nc = nc.merge(mapk, on=["donor", "acc", "junction_class"], how="left")
assert nc.cluster_id.notna().all()
cs = nc.groupby(["cluster_id", "sample_id"]).reads.sum().unstack(fill_value=0)
cs = cs.reindex(columns=cls.sample_id.tolist(), fill_value=0)
nlib = (cs >= DET_MIN).sum(axis=1).rename("n_detected_libs").reset_index()
CL = CL.merge(nlib, on="cluster_id")
CL.sort_values("reads_total", ascending=False).to_csv(D + "noncanonical_clusters_all.csv", index=False)
cs.to_csv(D + "noncanonical_cluster_x_sample_counts.csv.gz")
isdep = CL.junction_class == "TRSL_dep_noncanonical"
print("clusters dep", int(isdep.sum()), "indep", int((~isdep).sum()), "max acceptor width", int((CL.acc_max - CL.acc_min).max()),
      "max donor width (indep)", int((CL.donor_max - CL.donor_min)[~isdep].max()))

cand = CL[(CL.reads_total >= 100) & (CL.n_detected_libs >= 5)].sort_values("reads_total", ascending=False)
sel = cand.head(30).copy()
sel.to_csv(D + "panelD_selected_clusters.csv", index=False)
cnt = cs.loc[sel.cluster_id]
cnt.to_csv(D + "panelD_selected_cluster_x_sample_counts.csv")
vf = cls.set_index("sample_id").viral_fragments.reindex(cnt.columns)
(cnt.div(vf / M, axis=1)).to_csv(D + "panelD_selected_cluster_x_sample_per_M_viral_frag.csv")
depth = cls.set_index("sample_id").all_classified.reindex(cnt.columns)
share = cs.sum(axis=1) / cls.all_classified.sum()
expct = pd.DataFrame(np.outer(share.loc[sel.cluster_id].values, depth.values), index=cnt.index, columns=cnt.columns)
expct.to_csv(D + "panelD_selected_expected_reads_given_depth.csv")
informative = ~((cnt == 0) & (expct < 3))
informative.to_csv(D + "panelD_selected_informative_mask.csv")
meta = cls[["sample_id", "MOI", "treatment", "hpi"]].set_index("sample_id")
dst = []
for cid in sel.cluster_id:
    r = cnt.loc[cid]
    mt = meta.assign(det=(r >= DET_MIN).reindex(meta.index).values)
    d = dict(cluster_id=cid, reads_total=int(r.sum()), n_libs_detected=int(mt.det.sum()),
             n_MOI_x_hpi_detected=int(mt[mt.det].groupby(["MOI", "hpi"]).ngroups),
             detected_libs_MOI0p1=int(mt[mt.MOI == 0.1].det.sum()), detected_libs_MOI1=int(mt[mt.MOI == 1.0].det.sum()))
    for t in TIMES:
        d[f"detected_libs_{t}hpi"] = int(mt[mt.hpi == t].det.sum())
    d["n_uninformative_cells"] = int((~informative.loc[cid]).sum())
    dst.append(d)
pd.DataFrame(dst).to_csv(D + "panelD_detection_summary.csv", index=False)
ov = {}
for moi in (0.1, 1.0):
    sid = cls[cls.MOI == moi].sample_id
    ov[moi] = cs[sid].sum(axis=1).sort_values(ascending=False).head(30).index
ov_n = len(set(ov[0.1]) & set(ov[1.0]))
pd.DataFrame({"top30_MOI0.1": list(ov[0.1]), "top30_MOI1": list(ov[1.0]),
              "in_both": [c in set(ov[1.0]) for c in ov[0.1]]}).to_csv(D + "panelD_top30_overlap_between_MOI.csv", index=False)


dmin = int(cls.all_classified.min())
sel_set = set(sel.cluster_id)
ncm = nc.groupby(["sample_id", "cluster_id"]).reads.sum().reset_index()
res = []
for s in cls.sample_id:
    g = ncm[ncm.sample_id == s]
    N = int(cls.set_index("sample_id").all_classified[s])
    p = min(1.0, dmin / N)
    insel = g.cluster_id.isin(sel_set).values
    det_all, det_sel = [], []
    for _ in range(100):
        k = rng.binomial(g.reads.values, p)
        det_all.append(int((k >= 1).sum()))
        det_sel.append(int(((k >= DET_MIN) & insel).sum()))
    res.append(dict(sample_id=s, classified_junction_reads=N, downsample_p=p, target_depth=dmin,
                    nc_clusters_detected_full=int((g.reads >= 1).sum()),
                    nc_clusters_detected_downsampled_mean=float(np.mean(det_all)),
                    nc_clusters_detected_downsampled_sd=float(np.std(det_all)),
                    panelD_rows_detected_full=int(((g.reads >= DET_MIN).values & insel).sum()),
                    panelD_rows_detected_downsampled_mean=float(np.mean(det_sel))))
pd.DataFrame(res).merge(meta.reset_index(), on="sample_id").to_csv(D + "panelD_downsampling_sensitivity.csv", index=False)


nb = int(np.ceil(GL / BIN))
rowsE = []
ps2 = ps.copy()
ps2["is_canonical"] = ps2.junction_class == "canonical"
for s, g in ps2.groupby("sample_id"):
    for kind, col in (("5prime_donor", "donor"), ("3prime_acceptor", "acc")):
        for cset, gg in (("noncanonical", g[~g.is_canonical]), ("canonical", g[g.is_canonical])):
            h = np.bincount(((gg[col].values - 1) // BIN).astype(int), weights=gg.reads.values, minlength=nb)
            for b in np.nonzero(h)[0]:
                rowsE.append((s, kind, cset, b * BIN + 1, (b + 1) * BIN, h[b]))
E = pd.DataFrame(rowsE, columns=["sample_id", "breakpoint", "set", "bin_start", "bin_end", "reads"]).merge(
    cls[["sample_id", "MOI", "treatment", "hpi", "viral_fragments"]], on="sample_id")
E["reads_per_M_viral_frag"] = E.reads / (E.viral_fragments / M)
E.to_csv(D + "panelE_breakpoint_histogram_per_sample_bin20.csv.gz", index=False)
pe = E.groupby(["MOI", "hpi", "breakpoint", "set", "bin_start", "bin_end"]).reads.sum().reset_index()
den = cls.groupby(["MOI", "hpi"]).viral_fragments.sum().rename("den").reset_index()
pe = pe.merge(den, on=["MOI", "hpi"])
pe["reads_per_M_viral_frag"] = pe.reads / (pe.den / M)
pe.to_csv(D + "panelE_breakpoint_histogram_pooled_DESCRIPTIVE.csv.gz", index=False)
chk_all = E.groupby("breakpoint").reads.sum()
with open(D + "s3_QC.txt", "w", encoding="utf-8") as f:
    f.write("# Step3 QC\n")
    f.write(f"denominator: {DEN_NOTE}\n")
    f.write(f"zero counts in panelB sgRNA cells: {zero_B}; zero counts in panelC class cells: {zero_C}\n")
    f.write(f"panelE: total junction reads {int(ps.reads.sum())}; 5' histogram sum {int(chk_all['5prime_donor'])}; "
            f"3' histogram sum {int(chk_all['3prime_acceptor'])}\n")
    f.write(f"downsample target depth (min classified junction reads): {dmin}\n")
    f.write(f"top30 overlap between MOI: {ov_n}/30\n")
    f.write(f"D candidates {len(cand)}; selected {len(sel)}\n\nCHX vs HAR consistency (|diff|/SE > 3 flagged):\n")
    f.write(cons[cons.diff_over_SE.abs() > 3].round(3).to_string() + "\n")
print(open(D + "s3_QC.txt", encoding="utf-8").read())
