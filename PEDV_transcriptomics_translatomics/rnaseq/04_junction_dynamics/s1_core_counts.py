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
import sample_registry as SR
import re, os, sys, json
import numpy as np
import pandas as pd

BASE = PEDV_P+"GIb_rnaseq_reanalysis/"
OUT = PEDV_P+"Fig2_Supp_junction_dynamics/data/"
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', '03_classification_abundance'))
from fig1_moi_data import load_metadata

LEADER = (55, 72)
DIST = 1000
CANON = {"S": 20597, "ORF3": 24718, "E": 25342, "M": 25597, "N": 26285}
TOL = 5
RX = re.compile(r"(\d+)([MIDNSH=X])")


def anchors2(cig, gap):
    ops = [(int(l), o) for l, o in RX.findall(cig)]
    idx = [i for i, (l, o) in enumerate(ops) if o == "N" and l == gap]
    if not idx:
        return (0, 0)
    i = idx[0]
    la = sum(l for l, o in ops[:i] if o in "M=X")
    ba = 0
    for l, o in ops[i + 1:]:
        if o in "M=X":
            ba += l
        elif o in "NSH":
            break
    return la, ba


def classify(dn, ac):
    dep = LEADER[0] <= dn <= LEADER[1]
    can = ""
    for g, s in CANON.items():
        if abs(ac - s) <= TOL:
            can = g
    if dep:
        return ("canonical", can) if can else ("TRSL_dep_noncanonical", "")
    return ("TRSL_indep_" + ("local" if (ac - dn - 1) < DIST else "distant"), "")


def main():
    meta = load_metadata()
    ids = meta.sample_id.tolist()
    df = pd.read_csv(BASE + "junction_atlas/junction_events_raw.tsv.gz", sep="\t")
    d = df[df["sample"].isin(ids)].copy()
    aud = {"rows_20_samples": len(d)}
    d0 = d.copy()
    d = d[(~d.is_duplicate) & (~d.is_secondary) & (~d.is_supplementary) & (d.NH == 1)]
    aud["rows_after_flag_filters"] = len(d)
    ab = [anchors2(c, g) for c, g in zip(d.cigar, d.gap_len)]
    d["la"] = [a[0] for a in ab]; d["ba"] = [a[1] for a in ab]
    d["pass_qc"] = (d.mapq >= 20) & (d.la >= 15) & (d.ba >= 15)
    aud["rows_fail_qc"] = int((~d.pass_qc).sum())
    q = d[d.pass_qc].copy()
    q["donor"] = q.leader_end; q["acc"] = q.body_start
    aud["rows_pass_qc"] = len(q)
    # reads(fragments) with >1 distinct event
    u = q.drop_duplicates(["sample", "read_name", "donor", "acc"]).copy()
    ev_per_read = u.groupby(["sample", "read_name"]).size()
    aud["fragments_with_gt1_event"] = int((ev_per_read > 1).sum())
    aud["fragments_total_after_qc"] = int(len(ev_per_read))
    aud["junction_fragments_after_mate_dedup"] = len(u)

    f = pd.DataFrame({
        "rows_raw": d0.groupby("sample").size(),
        "rows_fail_qc": (~d.pass_qc).groupby(d["sample"]).sum(),
        "rows_pass_qc": q.groupby("sample").size(),
        "junction_fragments": u.groupby("sample").size(),
    }).reset_index().rename(columns={"sample": "sample_id"})
    f["qc_fail_frac"] = f.rows_fail_qc / f.rows_raw
    f.merge(meta[["sample_id", "MOI", "treatment", "hpi"]], on="sample_id").sort_values(["MOI", "hpi", "treatment"]).to_csv(
        OUT + "filter_audit_per_sample.csv", index=False)


    ps = u.groupby(["sample", "donor", "acc"]).size().rename("reads").reset_index()
    cl = ps[["donor", "acc"]].drop_duplicates().copy()
    k = [classify(a, b) for a, b in zip(cl.donor, cl.acc)]
    cl["junction_class"] = [x[0] for x in k]; cl["canonical_gene"] = [x[1] for x in k]
    cl["strand"] = "+"
    cl["event_id"] = (cl.donor.astype(str) + "|" + cl.acc.astype(str) + "|+|" + cl.junction_class +
                      np.where(cl.canonical_gene != "", "_" + cl.canonical_gene, ""))
    cl["span_nt"] = cl.acc - cl.donor - 1
    ps = ps.merge(cl, on=["donor", "acc"]).rename(columns={"sample": "sample_id"})
    ps = ps.merge(meta[["sample_id", "MOI", "treatment", "hpi"]], on="sample_id")
    ps.to_csv(OUT + "junction_event_sample_matrix_long.csv.gz", index=False)
    tot = ps.groupby("event_id").reads.sum().rename("reads_total")
    ns = ps.groupby("event_id").sample_id.nunique().rename("n_samples")
    et = cl.merge(tot, on="event_id").merge(ns, on="event_id").sort_values("reads_total", ascending=False)
    et.to_csv(OUT + "event_table.csv", index=False)

    piv = ps.assign(cls=np.where(ps.junction_class == "canonical", "canonical_" + ps.canonical_gene, ps.junction_class))\
        .pivot_table(index="sample_id", columns="cls", values="reads", aggfunc="sum", fill_value=0)
    piv = meta[["sample_id", "MOI", "treatment", "hpi"]].merge(piv.reset_index(), on="sample_id")
    for c in ["canonical_S", "canonical_ORF3", "canonical_E", "canonical_M", "canonical_N"]:
        assert c in piv
    piv["canonical_total"] = piv[[c for c in piv if c.startswith("canonical_")]].sum(axis=1)
    piv["all_classified"] = piv.canonical_total + piv.TRSL_dep_noncanonical + piv.TRSL_indep_distant + piv.TRSL_indep_local
    piv.sort_values(["MOI", "hpi", "treatment"]).to_csv(OUT + "per_sample_class_counts.csv", index=False)


    tot_all = int(ps.reads.sum())
    pooled = piv[["canonical_total", "TRSL_dep_noncanonical", "TRSL_indep_distant", "TRSL_indep_local"]].sum()
    lines = ["# Step1 对账: 本脚本 vs Figure 2 冻结结果", json.dumps(aud, indent=1)]
    lines.append(f"total junction reads (classified): {tot_all}   [Fig2: 2,674,258]")
    for k_, v in pooled.items():
        lines.append(f"  {k_:26s} {int(v):9d}  {100*v/tot_all:6.3f}%")
    lines.append("canonical by gene (%): " + ", ".join(f"{g} {100*piv['canonical_'+g].sum()/piv.canonical_total.sum():.2f}" for g in CANON))
    ev_main = pd.read_csv(BASE + "figures/Fig2_work/events_main.tsv", sep="\t")
    lines.append(f"Fig2 events_main.tsv: reads {int(ev_main.reads.sum())}, events {len(ev_main)}; mine events {len(et)}")
    m = et.merge(ev_main[["donor", "acceptor", "reads", "cls"]], left_on=["donor", "acc"], right_on=["donor", "acceptor"], how="outer",
                 indicator=True)
    lines.append("event-level merge: " + str(m._merge.value_counts().to_dict()) +
                 "; reads mismatch events: " + str(int((m.reads_total != m.reads).sum())))
    old = pd.read_csv(BASE + "figures/Fig2_per_sample_event_counts.csv.gz")
    a = old.groupby("sample").reads.sum()
    b = ps.groupby("sample_id").reads.sum()
    lines.append("per-sample total vs Fig2_per_sample_event_counts: max abs diff = " + str(int((a - b.reindex(a.index)).abs().max())))
    open(OUT + "s1_reconcile.txt", "w", encoding="utf-8").write("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
