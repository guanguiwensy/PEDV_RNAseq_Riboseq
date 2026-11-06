#!/usr/bin/env python3
"""分析一: Junction 四类统计 + 非经典 sgRNA 跨毒株同源/复现。
输入: _in/junc (gx_pe_junc.py 产物; 参考=LJX_merged_denovo / LLCPK_denovo), 冻结 GI-b: Fig2 事件矩阵/簇表 + GIb_atlas_v1.0 (phase5/_in/_work/sgRNA_atlas_GIb.tsv)。
复用: 簇规则(s03_highfreq.py: dep-nc 按 acceptor 链式<=5nt; 非依赖 donor&acceptor<=5nt 单连接), best_alignment_similarity (build_sgRNA_atlas_GIb.py 逻辑),
      junction_homology (module5_9 engine), 原项目对应 ±3 nt。分类规则 cls() 原样来自 p5_pe_junc.py(未改)。
输出: junction_composition_cross_strain.tsv, junction_recurrence_per_sample.tsv, junction_homology_GIb_LJX_LLCPK.tsv, noncanonical_sgRNA_similarity.tsv"""
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
import os, sys, json, gzip, collections, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gx_common import *
sys.path.insert(0, os.path.join(TASKS, "module5_9_tis_orf_infra", "engine"))
import pedv_orf_engine as E
IN = os.path.join(G, "_in"); F2 = os.path.join(TASKS, "Fig2_Supp_junction_dynamics", "data")
SAMPLES = [("SRR27751928", "LJX_P1", "LJX", "LJXm"), ("SRR27751927", "LJX_P2", "LJX", "LJXm"), ("SRR27751926", "LJX_P3", "LJX", "LJXm"), ("SRR31828424", "LLCPK_P10_1", "LLC-PK", "LLCPK")]
CLS = ["canonical", "TRSL_dep_noncanonical", "TRSL_indep_distant", "TRSL_indep_local"]
def norm(c): return "canonical" if c.startswith("TRSL_dep_canonical") else c
reads = []; viral = {}
for srr, lab, ds, ref in SAMPLES:
    d = pd.read_csv(os.path.join(IN, "junc", f"{srr}.{ref}.reads.tsv.gz"), sep="\t", dtype={"read": str})
    d["sample"] = lab; d["dataset"] = ds; d["cls"] = d["class"].map(norm)
    d["gene"] = d["class"].where(d["class"].str.startswith("TRSL_dep_canonical"), None).str.replace("TRSL_dep_canonical_", "", regex=False)
    reads.append(d); viral[lab] = json.load(open(os.path.join(IN, "junc", f"{srr}.{ref}.summary.json")))["primary_reads_in_prefilter"]
R = pd.concat(reads, ignore_index=True); n_nonhom = R[R.cls == "nonhomologous_end"].groupby("sample").size().to_dict()
R = R[R.cls != "nonhomologous_end"].copy()
# 旧参考(LJXc/LLCc)结果, 并列比较
OLD = {"LJX_P1": "SRR27751928.LJXc", "LJX_P2": "SRR27751927.LJXc", "LJX_P3": "SRR27751926.LJXc", "LLCPK_P10_1": "SRR31828424.LLCc"}
old = {}
for lab, key in OLD.items():
    d = pd.read_csv(os.path.join(ESV, "phase5", "_in5", "junc", key + ".reads.tsv.gz"), sep="\t"); d["cls"] = d["class"].map(norm); d = d[d.cls != "nonhomologous_end"]
    old[lab] = d
# 冻结 GI-b
M = pd.read_csv(os.path.join(F2, "junction_event_sample_matrix_long.csv.gz")); M = M.rename(columns={"acc": "acceptor", "junction_class": "cls"})
PA = pd.read_csv(os.path.join(F2, "panelA_per_sample_class_fractions.csv")); GIB_VIRAL = int(PA.viral_fragments.sum())
OC = pd.read_csv(os.path.join(F2, "noncanonical_clusters_all.csv")); GIB_TOT = int(M.reads.sum())
ATL = pd.read_csv(os.path.join(ESV, "phase5", "_in", "_work", "sgRNA_atlas_GIb.tsv"), sep="\t")
# ---------------- 1. 四类组成
rows = []
def comp(name, dataset, ref, df, vir, pooled_note=""):
    N = len(df)
    ev = df.groupby("cls").apply(lambda g: g[["donor", "acceptor"]].drop_duplicates().shape[0], include_groups=False).to_dict() if "donor" in df else {}
    Nd = sum(ev.get(c, 0) for c in CLS)
    for c in CLS:
        n = int((df.cls == c).sum()); e = ev.get(c, 0)
        rows.append(dict(dataset=dataset, sample=name, reference=ref, junction_class=c, junction_reads=n, distinct_junctions=e, pct_of_junction_reads=round(100 * n / N, 4) if N else np.nan,
                         pct_of_distinct_junctions=round(100 * e / Nd, 4) if Nd else np.nan, per_million_viral_reads=round(1e6 * n / vir, 2) if vir else np.nan,
                         total_junction_reads_sample=N, viral_reads_basis=vir, note=pooled_note))
for srr, lab, ds, ref in SAMPLES:
    comp(lab, ds, ref, R[R["sample"] == lab], viral[lab])
lj = R[R.dataset == "LJX"]; comp("LJX_pooled(P1+P2+P3)", "LJX", "LJXm", lj, sum(viral[l] for l in ("LJX_P1", "LJX_P2", "LJX_P3")), "三只合并统计; 不是三个独立重复")
gb = M.rename(columns={"acceptor": "acceptor"}).copy(); gb["donor"] = gb.donor
Ng = int(gb.reads.sum())
for c in CLS:
    g = gb[gb.cls == c]; n = int(g.reads.sum()); e = g[["donor", "acceptor"]].drop_duplicates().shape[0]
    rows.append(dict(dataset="GI-b(Vero,冻结Fig2)", sample="Vero_20libs_pooled", reference="PEDV_GIb v1.1", junction_class=c, junction_reads=n, distinct_junctions=e, pct_of_junction_reads=round(100 * n / Ng, 4),
                     pct_of_distinct_junctions=np.nan, per_million_viral_reads=round(1e6 * n / GIB_VIRAL, 2), total_junction_reads_sample=Ng, viral_reads_basis=GIB_VIRAL, note="原项目冻结矩阵(不重算); distinct=按 (donor,acceptor) 去重"))
tmp = pd.DataFrame(rows)
tot_d = tmp.groupby("sample").distinct_junctions.transform("sum"); m_ = tmp.dataset.str.startswith("GI-b")
tmp.loc[m_, "pct_of_distinct_junctions"] = (100 * tmp.distinct_junctions / tot_d)[m_].round(4)
for c in CLS: pass
for i, r in tmp.iterrows():
    if r["sample"] in old:
        o = old[r["sample"]]; tmp.loc[i, "old_ref_junction_reads"] = int((o.cls == r.junction_class).sum()); tmp.loc[i, "old_ref_pct"] = round(100 * (o.cls == r.junction_class).mean(), 4)
        tmp.loc[i, "old_ref_label"] = "LJXc" if r["sample"].startswith("LJX") else "LLCc"
write_tsv(tmp, "junction_composition_cross_strain.tsv")
# GI-b 每库(Vero 20 libs)四类比例, 供 Fig1
pa = PA[["sample_id", "MOI", "treatment", "hpi", "f_canonical", "f_TRSL_dep_noncanonical", "f_TRSL_indep_distant", "f_TRSL_indep_local", "all_classified"]]
write_tsv(pa, "junction_composition_GIb_per_library.tsv")
# ---------------- 2. 簇(GI-b 坐标), 同源归并 (s03 规则)
ev = R.groupby(["sample", "dataset", "cls", "gene", "donor", "acceptor", "donor_GIb", "acceptor_GIb"], dropna=False).size().rename("reads").reset_index(); ev["gene"] = ev.gene.fillna("")
Ns = R.groupby("sample").size().to_dict(); ev["N"] = ev["sample"].map(Ns)
def chain_1d(vals, w=5):
    vals = sorted(set(vals)); groups = []; cur = [vals[0]]
    for v in vals[1:]:
        if v - cur[-1] <= w: cur.append(v)
        else: groups.append(cur); cur = [v]
    groups.append(cur); return groups
cid = {}; n = 0
for k in ev[ev.cls == "canonical"][["gene", "donor_GIb", "acceptor_GIb"]].drop_duplicates().itertuples(index=False):
    n += 1; cid[("canonical", k.gene, k.donor_GIb, k.acceptor_GIb)] = (f"C{n:03d}", k.donor_GIb, k.donor_GIb, k.acceptor_GIb, k.acceptor_GIb)
dep = ev[ev.cls == "TRSL_dep_noncanonical"]; n = 0; dmap = {}
for grp in chain_1d(dep.acceptor_GIb.tolist()):
    n += 1
    for a in grp: dmap[a] = (f"D{n:03d}", 55, 72, min(grp), max(grp))
imap = {}; n = 0
for c in ("TRSL_indep_distant", "TRSL_indep_local"):
    sub = ev[ev.cls == c][["donor_GIb", "acceptor_GIb"]].drop_duplicates().values.tolist(); par = list(range(len(sub)))
    def find(i):
        while par[i] != i: par[i] = par[par[i]]; i = par[i]
        return i
    order = sorted(range(len(sub)), key=lambda i: sub[i][1])
    for ii, i in enumerate(order):
        for j in order[ii + 1:]:
            if sub[j][1] - sub[i][1] > 5: break
            if abs(sub[j][0] - sub[i][0]) <= 5: par[find(i)] = find(j)
    rt = collections.defaultdict(list)
    for i in range(len(sub)): rt[find(i)].append(i)
    for r_, mem in rt.items():
        n += 1; pre = "T" if c == "TRSL_indep_distant" else "L"
        ds_ = [sub[i][0] for i in mem]; as_ = [sub[i][1] for i in mem]
        for i in mem: imap[(c, sub[i][0], sub[i][1])] = (f"{pre}{n:04d}", min(ds_), max(ds_), min(as_), max(as_))
def assign(r):
    if r.cls == "canonical": return cid[("canonical", r.gene, r.donor_GIb, r.acceptor_GIb)]
    if r.cls == "TRSL_dep_noncanonical": return dmap[r.acceptor_GIb]
    return imap[(r.cls, r.donor_GIb, r.acceptor_GIb)]
a_ = ev.apply(assign, axis=1, result_type="expand"); a_.columns = ["cluster_id", "d_min", "d_max", "a_min", "a_max"]; ev = pd.concat([ev, a_], axis=1)
# 液位: 映射模糊性 (liftover class)
LF = pd.read_csv(os.path.join(G, "tables", "GII_to_GIb_coordinate_liftover.tsv"), sep="\t")
LFd = {ds: LF[LF.dataset == lab].set_index("gII_pos").mapping_class.to_dict() for ds, lab in (("LJX", "LJX_merged"), ("LLC-PK", "LLC-PK"))}
def amb(ds, pos): return LFd[ds].get(int(pos), "NA")
ev["donor_map_class"] = [amb(d, p) for d, p in zip(ev.dataset, ev.donor)]; ev["acc_map_class"] = [amb(d, p) for d, p in zip(ev.dataset, ev.acceptor)]
ev["mapping_ambiguous"] = ev.donor_map_class.str.startswith("ambiguous") | ev.acc_map_class.str.startswith("ambiguous") | ev.donor_map_class.str.startswith("unmappable") | ev.acc_map_class.str.startswith("unmappable")
ev.to_csv(os.path.join(IN, "gx_events_with_cluster.tsv"), sep="\t", index=False)
# 簇表
labs = [s[1] for s in SAMPLES]
def nsamp(df, cols, thr): return (df[cols] >= thr).sum(axis=1)
cl_rows = []
for c_, g in ev.groupby("cluster_id"):
    rep = g.sort_values("reads", ascending=False).iloc[0]
    r = dict(cluster_id=c_, junction_class=rep.cls, canonical_gene=rep.gene if rep.cls == "canonical" else "", rep_donor_GIb=int(rep.donor_GIb), rep_acceptor_GIb=int(rep.acceptor_GIb),
             donor_GIb_min=int(rep.d_min), donor_GIb_max=int(rep.d_max), acceptor_GIb_min=int(rep.a_min), acceptor_GIb_max=int(rep.a_max))
    for s in labs: r[f"reads|{s}"] = int(g[g["sample"] == s].reads.sum())
    r["any_mapping_ambiguous"] = bool(g.mapping_ambiguous.any())
    cl_rows.append(r)
CL = pd.DataFrame(cl_rows)
for lab_, ds in (("LJX", ["LJX_P1", "LJX_P2", "LJX_P3"]),):
    CL["LJX_pooled_reads"] = CL[[f"reads|{s}" for s in ds]].sum(axis=1)
    CL["LJX_n_piglets_ge1"] = nsamp(CL, [f"reads|{s}" for s in ds], 1); CL["LJX_n_piglets_ge3"] = nsamp(CL, [f"reads|{s}" for s in ds], 3)
CL["LLCPK_reads"] = CL["reads|LLCPK_P10_1"]
for s in labs: CL[f"per10k|{s}"] = (1e4 * CL[f"reads|{s}"] / Ns[s]).round(3)
CL["LJX_per10k_pooled"] = (1e4 * CL.LJX_pooled_reads / sum(Ns[s] for s in labs[:3])).round(3)
CL["LJX_recurrent(>=3reads in >=2 piglets)"] = CL.LJX_n_piglets_ge3 >= 2
CL["LLCPK_detected(>=3reads; single sample)"] = CL.LLCPK_reads >= 3
# ---------------- 3. 与冻结 GI-b 对应 (±3 nt, 同 s03 orig_match)
PAD = 3
def gib_match(r):
    if r.junction_class == "canonical":
        e = M[(M.cls == "canonical") & (M.canonical_gene == r.canonical_gene) & (M.acceptor == r.rep_acceptor_GIb) & (M.donor == r.rep_donor_GIb)]
        t = M[(M.cls == "canonical") & (M.canonical_gene == r.canonical_gene)]
        if len(e): return pd.Series(dict(GIb_Fig2_cluster=f"{r.canonical_gene}:{r.rep_donor_GIb}>{r.rep_acceptor_GIb}", GIb_Fig2_reads=int(e.reads.sum()), GIb_Fig2_libs=int(e.sample_id.nunique()), GIb_match_level="exact(donor,acceptor)"))
        return pd.Series(dict(GIb_Fig2_cluster="", GIb_Fig2_reads=0, GIb_Fig2_libs=0, GIb_match_level="none(canonical GI-b 事件无 exact)"))
    s = OC[OC.junction_class == r.junction_class]
    m = (s.acc_max >= r.rep_acceptor_GIb - PAD) & (s.acc_min <= r.rep_acceptor_GIb + PAD)
    if r.junction_class != "TRSL_dep_noncanonical": m &= (s.donor_max >= r.rep_donor_GIb - PAD) & (s.donor_min <= r.rep_donor_GIb + PAD)
    s = s[m]
    if not len(s): return pd.Series(dict(GIb_Fig2_cluster="", GIb_Fig2_reads=0, GIb_Fig2_libs=0, GIb_match_level="none"))
    s = s.assign(dist=(s.acc_rep - r.rep_acceptor_GIb).abs()).sort_values(["dist", "reads_total"], ascending=[True, False]); t = s.iloc[0]
    return pd.Series(dict(GIb_Fig2_cluster=t.cluster_id, GIb_Fig2_reads=int(t.reads_total), GIb_Fig2_libs=int(t.n_detected_libs), GIb_match_level=f"rep±{PAD}nt(n_candidate={len(s)})"))
CL = pd.concat([CL, CL.apply(gib_match, axis=1)], axis=1)
def atl(r):
    if r.junction_class != "TRSL_dep_noncanonical": return ""
    a = ATL[(ATL.type == "noncanonical") & (ATL.body_cluster_start.astype(int) - PAD <= r.rep_acceptor_GIb) & (ATL.body_cluster_end.astype(int) + PAD >= r.rep_acceptor_GIb)]
    return ";".join(a.sgRNA_id) if len(a) else ""
CL["GIb_atlas_v1.0_id"] = CL.apply(atl, axis=1)
write_tsv(CL.sort_values(["junction_class", "LJX_pooled_reads"], ascending=[True, False]), "junction_recurrence_per_sample.tsv")
# ---------------- 4. 以冻结 GI-b atlas 为主轴的同源表
SEQ_G = GIB
def flank(seq, pos, up, dn): return seq[max(0, pos - 1 - up):pos - 1 + dn]
def best_sim(a, b, max_shift=10):
    best = (0.0, 0, 0)
    for sh in range(-max_shift, max_shift + 1):
        m = t = 0
        for i in range(len(a)):
            j = i + sh
            if 0 <= j < len(b): t += 1; m += a[i] == b[j]
        if t >= 8 and m / t > best[0]: best = (m / t, t, sh)
    return best
def ident(a, b):
    n = min(len(a), len(b)); return round(sum(x == y for x, y in zip(a[:n], b[:n])) / n, 3) if n else np.nan
def trs_dist(s6, ref="CTAAAC"): return sum(x != y for x, y in zip(s6, ref)) if len(s6) == 6 else np.nan
reg = pd.read_csv(os.path.join(G, "tables", "GII_to_GIb_liftover_regions.tsv"), sep="\t")
indel_tab = {ds: reg[(reg.dataset == lab) & reg.type.isin(["ins_in_new", "del_in_new"])] for ds, lab in (("LJX", "LJX_merged"), ("LLC-PK", "LLC-PK"))}
def near_indel(ds, gpos, w=15):
    t = indel_tab[ds]; x = t[(t.GIb_pos_adjacent >= gpos - w) & (t.GIb_pos_adjacent <= gpos + w)]
    return ";".join(f"{r.type}:{r.length_nt}nt@GIb{r.GIb_pos_adjacent}" for r in x.itertuples())
def seq_evidence(ds, d, a, Gd, Ga):
    """GII 事件(d,a native) vs GI-b 同源事件(Gd,Ga) 的序列证据"""
    sq = SEQ[ds]
    acc_g = flank(GIB, Ga, 12, 12); acc_q = flank(sq, a, 12, 12)
    dn_g = GIB[max(0, Gd - 12):Gd]; dn_q = sq[max(0, d - 12):d]
    lead_q = sq[int(G2Q[ds][66]) - 26:int(G2Q[ds][66]) + 5] if G2Q[ds][66] else ""
    lead_g = flank(GIB, 66, 25, 5)
    sim_q = best_sim(lead_q, flank(sq, a, 25, 15))[0] if lead_q else np.nan; sim_g = best_sim(lead_g, flank(GIB, Ga, 25, 15))[0]
    return dict(acceptor_window24_identity=ident(acc_g, acc_q), leader_last12_identity=ident(dn_g, dn_q),
                TRS_B_6mer_GIb=GIB[Ga - 7:Ga - 1], TRS_B_6mer_GII=sq[a - 7:a - 1], TRS_B_dist_to_CTAAAC_GIb=trs_dist(GIB[Ga - 7:Ga - 1]), TRS_B_dist_to_CTAAAC_GII=trs_dist(sq[a - 7:a - 1]),
                TRS_similarity_GIb=round(sim_g, 3), TRS_similarity_GII=round(sim_q, 3) if sim_q == sim_q else np.nan,
                right_homology_GIb=E.junction_homology(GIB, Gd, Ga), right_homology_GII=E.junction_homology(sq, d, a),
                acceptor_flank_GIb_first12=GIB[Ga - 1:Ga + 11], acceptor_flank_GII_first12=sq[a - 1:a + 11], donor_flank_GIb_last12=dn_g, donor_flank_GII_last12=dn_q)
TRSL_core = {ds: SEQ[ds][int(G2Q[ds][67]) - 1:int(G2Q[ds][72])] if G2Q[ds][67] and G2Q[ds][72] else "NA" for ds in DATASETS}
TRSL_core["GI-b"] = GIB[66:72]
json.dump(TRSL_core, open(os.path.join(IN, "TRSL_core.json"), "w"))
# GI-b 总 reads 与 power
tot_by_sample = Ns
def status(reads_list, exp_list, ds):
    if ds == "LJX":
        ge3 = sum(x >= 3 for x in reads_list); anyr = sum(x >= 1 for x in reads_list)
        if ge3 >= 2: return "recurrent(>=3 reads in >=2 piglets)"
        if ge3 == 1: return "single_piglet_ge3"
        if anyr: return "subthreshold(1-2 reads)"
        return "not_detected_adequate_depth" if min(exp_list) >= 5 or np.mean(exp_list) >= 5 else "not_detected_low_depth(expected<5)"
    x = reads_list[0]
    if x >= 3: return "detected(>=3 reads; single sample)"
    if x >= 1: return "subthreshold(1-2 reads)"
    return "not_detected_adequate_depth" if exp_list[0] >= 5 else "not_detected_low_depth(expected<5)"
hrows = []; srows = []
Nlj = [Ns[s] for s in labs[:3]]; Nll = Ns["LLCPK_P10_1"]
def best_event(ds, cids):
    g = ev[(ev.dataset == ds) & (ev.cluster_id.isin(cids))]
    if not len(g): return None
    return g.groupby(["donor", "acceptor", "donor_GIb", "acceptor_GIb"]).reads.sum().sort_values(ascending=False).reset_index().iloc[0]
atl_nc = ATL[ATL.type == "noncanonical"]
# 对每个 atlas nc: 找到 Fig2 簇(给出冻结丰度)
for a in atl_nc.itertuples(index=False):
    rep = int(a.representative_breakpoint); lo = int(a.body_cluster_start) - PAD; hi = int(a.body_cluster_end) + PAD
    s = OC[(OC.junction_class == "TRSL_dep_noncanonical") & (OC.acc_max >= lo) & (OC.acc_min <= hi)]
    fr = (s.reads_total.sum() if len(s) else 0); frl = (s.n_detected_libs.max() if len(s) else 0)
    frac = fr / GIB_TOT
    cids = CL[(CL.junction_class == "TRSL_dep_noncanonical") & (CL.acceptor_GIb_max >= lo) & (CL.acceptor_GIb_min <= hi)].cluster_id.tolist()
    g = CL[CL.cluster_id.isin(cids)]
    lj_reads = [int(g[f"reads|{s_}"].sum()) for s_ in labs[:3]]; ll_reads = int(g["reads|LLCPK_P10_1"].sum())
    exp_lj = [frac * n_ for n_ in Nlj]; exp_ll = frac * Nll
    st_lj = status(lj_reads, exp_lj, "LJX"); st_ll = status([ll_reads], [exp_ll], "LLC-PK")
    row = dict(GIb_sgRNA_id=a.sgRNA_id, junction_class="TRSL_dep_noncanonical", GIb_region=a.genomic_region, GIb_rep_acceptor=rep, GIb_acceptor_range=f"{a.body_cluster_start}-{a.body_cluster_end}", GIb_atlas_confidence=a.confidence_level,
               GIb_atlas_positive_samples=int(a.positive_samples), GIb_Fig2_reads=int(fr), GIb_Fig2_libs_of20=int(frl), GIb_fraction_of_all_junction_reads=round(frac, 7),
               LJX_P1_reads=lj_reads[0], LJX_P2_reads=lj_reads[1], LJX_P3_reads=lj_reads[2], LJX_pooled_reads=sum(lj_reads), LLCPK_reads=ll_reads,
               LJX_expected_reads_P1_P2_P3=";".join(f"{x:.1f}" for x in exp_lj), LLCPK_expected_reads=round(exp_ll, 1), LJX_status=st_lj, LLCPK_status=st_ll,
               GII_cluster_ids=";".join(cids), LJX_per10k_pooled=round(1e4 * sum(lj_reads) / sum(Nlj), 3), LLCPK_per10k=round(1e4 * ll_reads / Nll, 3))
    for ds, tag in (("LJX", "LJX"), ("LLC-PK", "LLCPK")):
        be = best_event(ds, cids)
        if be is None:
            for k in ("native_donor", "native_acceptor", "GIb_donor_lifted", "GIb_acceptor_lifted", "breakpoint_shift_vs_GIb_rep", "shift_explanation", "mapping_ambiguity"): row[f"{tag}_{k}"] = ""
            continue
        shift = int(be.acceptor_GIb) - rep
        near = near_indel(ds, int(be.acceptor_GIb)); rh = E.junction_homology(SEQ[ds], int(be.donor), int(be.acceptor))
        expl = "none(0 nt)" if shift == 0 else (f"within_microhomology_ambiguity(right_homology={rh}nt)" if abs(shift) <= max(rh, 0) else (f"near_indel({near})" if near else f"within_cluster_range(±{PAD}nt, chain-merged)" if lo + PAD <= be.acceptor_GIb <= hi - PAD else "unexplained_shift"))
        amb_ = ev[(ev.dataset == ds) & (ev.donor == be.donor) & (ev.acceptor == be.acceptor)].mapping_ambiguous.any()
        row.update({f"{tag}_native_donor": int(be.donor), f"{tag}_native_acceptor": int(be.acceptor), f"{tag}_GIb_donor_lifted": int(be.donor_GIb), f"{tag}_GIb_acceptor_lifted": int(be.acceptor_GIb),
                    f"{tag}_breakpoint_shift_vs_GIb_rep": shift, f"{tag}_shift_explanation": expl, f"{tag}_mapping_ambiguity": "ambiguous_near_indel_or_unmappable" if amb_ else "unique"})
        # 序列证据 (GI-b 端以 atlas 代表 acceptor + 该事件 donor_GIb 近似)
        Gd = int(be.donor_GIb) if be.donor_GIb else 66
        ev_ = seq_evidence(ds, int(be.donor), int(be.acceptor), Gd, rep)
        srows.append(dict(GIb_sgRNA_id=a.sgRNA_id, dataset=ds, GIb_donor=Gd, GIb_acceptor=rep, native_donor=int(be.donor), native_acceptor=int(be.acceptor), event_reads=int(be.reads),
                          breakpoint_conserved="exact" if shift == 0 else f"shift{shift:+d}nt:{expl}", **ev_,
                          reads_per10k_GII=round(1e4 * be.reads / (Ns[labs[2]] if False else (sum(Nlj) if ds == "LJX" else Nll)), 3) if False else round(1e4 * (sum(lj_reads) if ds == "LJX" else ll_reads) / (sum(Nlj) if ds == "LJX" else Nll), 3),
                          GIb_per10k=round(1e4 * frac, 3), LJX_piglets_ge1=sum(x >= 1 for x in lj_reads) if ds == "LJX" else "", LJX_piglets_ge3=sum(x >= 3 for x in lj_reads) if ds == "LJX" else ""))
    # 分类
    lj_ok = st_lj.startswith("recurrent"); ll_ok = st_ll.startswith("detected"); lj_any = lj_ok or st_lj.startswith("single") ; ll_any = ll_ok
    lj_nd = st_lj.startswith("not_detected_adequate"); ll_nd = st_ll.startswith("not_detected_adequate")
    lj_na = st_lj.startswith("not_detected_low"); ll_na = st_ll.startswith("not_detected_low")
    if lj_ok and ll_ok: v = "shared_GIb_LJX_LLCPK"
    elif (lj_any and ll_any) : v = "shared_GIb_LJX(single piglet)_LLCPK"
    elif lj_ok and ll_nd: v = "GIb_LJX_only(LLCPK 足够深度未检出)"
    elif lj_ok and ll_na: v = "GIb_LJX(LLCPK 深度不足)"
    elif ll_ok and lj_nd: v = "GIb_LLCPK_only(LJX 足够深度未检出)"
    elif ll_ok and lj_na: v = "GIb_LLCPK(LJX 深度不足)"
    elif ll_ok or lj_any: v = "GIb+one_GII_partial"
    elif lj_nd and ll_nd: v = "GIb_only(两 GII 数据集足够深度未检出;潜在毒株特异)"
    elif lj_na and ll_na: v = "not_assessable(两 GII 数据集深度不足)"
    else: v = "GIb_only/partial_depth(见各自状态)"
    row["cross_strain_category"] = v
    hrows.append(row)
H = pd.DataFrame(hrows)
# canonical 5 个
for gene in ("S", "ORF3", "E", "M", "N"):
    g = CL[(CL.junction_class == "canonical") & (CL.canonical_gene == gene)]
    mm = M[(M.cls == "canonical") & (M.canonical_gene == gene)]
    row = dict(GIb_sgRNA_id=f"canonical_{gene}", junction_class="canonical", GIb_region=gene, GIb_Fig2_reads=int(mm.reads.sum()), GIb_fraction_of_all_junction_reads=round(mm.reads.sum() / GIB_TOT, 5),
               LJX_P1_reads=int(g["reads|LJX_P1"].sum()), LJX_P2_reads=int(g["reads|LJX_P2"].sum()), LJX_P3_reads=int(g["reads|LJX_P3"].sum()), LJX_pooled_reads=int(g.LJX_pooled_reads.sum()), LLCPK_reads=int(g.LLCPK_reads.sum()),
               LJX_status="recurrent(>=3 reads in >=2 piglets)" if (g[["reads|LJX_P1", "reads|LJX_P2", "reads|LJX_P3"]].sum() >= 3).sum() >= 2 else "see counts", LLCPK_status="detected" if g.LLCPK_reads.sum() >= 3 else "not_detected",
               GII_cluster_ids=";".join(g.cluster_id), cross_strain_category="canonical_all_three")
    hrows.append(row)
H = pd.DataFrame(hrows)
# GII 有而 atlas 无对应的 dep-nc / distant / local (GII-only candidates)
matched = set(";".join(H.GII_cluster_ids.fillna("")).split(";"))
extra = CL[(CL.junction_class != "canonical") & (~CL.cluster_id.isin(matched))].copy()
extra = extra[(extra.LJX_pooled_reads + extra.LLCPK_reads) >= 3]
ex = extra.copy()
ex_out = pd.DataFrame(dict(GIb_sgRNA_id="(无 atlas 对应)", junction_class=ex.junction_class, GIb_region="", GIb_rep_acceptor=ex.rep_acceptor_GIb, GIb_acceptor_range=ex.acceptor_GIb_min.astype(str) + "-" + ex.acceptor_GIb_max.astype(str),
                           GIb_Fig2_reads=ex.GIb_Fig2_reads, GIb_Fig2_libs_of20=ex.GIb_Fig2_libs, LJX_P1_reads=ex["reads|LJX_P1"], LJX_P2_reads=ex["reads|LJX_P2"], LJX_P3_reads=ex["reads|LJX_P3"],
                           LJX_pooled_reads=ex.LJX_pooled_reads, LLCPK_reads=ex.LLCPK_reads, GII_cluster_ids=ex.cluster_id, GIb_match_level=ex.GIb_match_level))
ex_out["cross_strain_category"] = np.where(ex_out.GIb_Fig2_reads > 0, "GII_detected;GIb_Fig2_cluster_present(非 atlas 层面)", "GII_only(冻结 GI-b 中无对应簇;潜在毒株特异)")
H2 = pd.concat([H, ex_out], ignore_index=True)
write_tsv(H2, "junction_homology_GIb_LJX_LLCPK.tsv")
S = pd.DataFrame(srows); write_tsv(S, "noncanonical_sgRNA_similarity.tsv")
print(H.cross_strain_category.value_counts()); print(H.LJX_status.value_counts()); print(H.LLCPK_status.value_counts())
print("n_atlas", len(atl_nc), "extra GII-only >=3 reads", len(ex_out), ex_out.cross_strain_category.value_counts().to_dict())
