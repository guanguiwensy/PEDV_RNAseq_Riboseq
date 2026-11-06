#!/usr/bin/env python3
"""Step 5: breakpoint hotspot dynamics(panel E 的数值依据) + Fig.2 hotspots a-f 的逐样本/逐 MOI 检出。描述性。"""
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
import numpy as np
import pandas as pd

D = PEDV_P+"Fig2_Supp_junction_dynamics/data/"
G2 = PEDV_P+"GIb_rnaseq_reanalysis/figures/"
TIMES = [3, 5, 8, 12, 16]
ps = pd.read_csv(D + "junction_event_sample_matrix_long.csv.gz")
cls = pd.read_csv(D + "per_sample_class_counts.csv")
import os
if os.path.exists(D + "viral_fragments_mm2.tsv"):
    dn = pd.read_csv(D + "viral_fragments_mm2.tsv", sep="\t")
    vf = dn.set_index("sample_id").viral_frag_primary_q20
else:
    st = pd.read_csv(PEDV_P+"GIb_rnaseq_reanalysis/figdata/fig1_counts.tsv", sep="\t",
                     names=["sample_id", "ver", "t", "pedv", "a", "b"])
    vf = st[st.ver == "raw"].set_index("sample_id").pedv
cls["vf"] = cls.sample_id.map(vf)
cls["nc_total"] = cls.TRSL_dep_noncanonical + cls.TRSL_indep_distant + cls.TRSL_indep_local
nc = ps[ps.junction_class != "canonical"]

# ---- Fig.2 hotspots a-f ----
hs = pd.read_csv(G2 + "Fig2A_hotspots_source_data.csv")
rows = []
for h in hs.itertuples():
    g = nc[(nc.donor == h.donor) & (nc.acc == h.acceptor)].set_index("sample_id").reads
    for r in cls.itertuples():
        n = int(g.get(r.sample_id, 0))
        rows.append(dict(label=h.label, donor=h.donor, acceptor=h.acceptor, sample_id=r.sample_id, MOI=r.MOI, treatment=r.treatment,
                         hpi=r.hpi, reads=n, reads_per_M_viral_frag=n / (r.vf / 1e6), share_of_nc_reads=n / r.nc_total))
H = pd.DataFrame(rows)
H.to_csv(D + "hotspots_a_f_per_sample.csv", index=False)
summ = H.groupby(["label", "donor", "acceptor", "MOI"]).agg(libs_detected_ge3=("reads", lambda x: int((x >= 3).sum())),
                                                          libs_detected_ge1=("reads", lambda x: int((x >= 1).sum())),
                                                          reads=("reads", "sum")).reset_index()
summ.to_csv(D + "hotspots_a_f_detection_by_MOI.csv", index=False)

# ---- 100-nt bin 层面的稳定热点: 汇总全部非 canonical, 取 5' 与 3' 各前 8 个互不相邻(>=300 nt)峰 ----
RB = 100
out = []
for bp, col in (("5prime_donor", "donor"), ("3prime_acceptor", "acc")):
    b = (nc[col] - 1) // RB
    t = nc.assign(b=b).groupby("b").reads.sum().sort_values(ascending=False)
    picked = []
    for bi, v in t.items():
        if all(abs(bi - p) >= 3 for p in picked):
            picked.append(bi)
        if len(picked) == 8:
            break
    for bi in picked:
        sub = nc[(b == bi)]
        for r in cls.itertuples():
            n = int(sub[sub.sample_id == r.sample_id].reads.sum())
            out.append(dict(breakpoint=bp, bin_start=int(bi * RB + 1), bin_end=int((bi + 1) * RB), sample_id=r.sample_id, MOI=r.MOI,
                            treatment=r.treatment, hpi=r.hpi, reads=n, reads_per_M_viral_frag=n / (r.vf / 1e6),
                            share_of_nc_reads=n / r.nc_total))
O = pd.DataFrame(out)
O.to_csv(D + "panelE_hotspot_bins_per_sample.csv", index=False)
pv = O.groupby(["breakpoint", "bin_start", "MOI", "hpi"]).agg(reads=("reads", "sum"), vf=("reads", "size")).reset_index()
# 读数占比随 hpi 的变化(描述): 每个热点 bin 占同 MOI×hpi 非 canonical reads 的比例
tot = cls.groupby(["MOI", "hpi"]).nc_total.sum().rename("nc").reset_index()
pv = pv.merge(tot, on=["MOI", "hpi"])
pv["share_of_nc_reads"] = pv.reads / pv.nc
w = pv.pivot_table(index=["breakpoint", "bin_start", "MOI"], columns="hpi", values="share_of_nc_reads")
w.to_csv(D + "panelE_hotspot_share_of_nc_reads_by_MOI_hpi.csv")
pd.set_option("display.width", 220)
print("== Fig2 hotspots a-f: detection by MOI (libs with >=3 reads)")
print(summ.pivot_table(index=["label", "donor", "acceptor"], columns="MOI", values=["libs_detected_ge3", "reads"]).to_string())
print("\n== top-bin share of nc reads by MOI x hpi (selected)")
print((w * 100).round(2).to_string())
