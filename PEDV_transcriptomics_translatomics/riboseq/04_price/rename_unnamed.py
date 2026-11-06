# -*- coding: utf-8 -*-
"""Give every 'unnamed (not a Fig. 6 +1 family)' ORF a systematic name: <parent>+<frame>#<k>.
parent = canonical ORF the ORF overlaps (as in coordinate_class), frame = +1/+2 relative to the parent's reading frame, k = running index (ordered by start, then stop)
within (parent, frame) over the union of all ORFs in the four annotated_29nt tables. Non-overlapping ORFs: 5'UTR#k, 3'UTR#k, IG#k. Fig. 6 family names (S-x1 ...) and canonical/same-frame labels are kept.
Original annotated_29nt files are not modified; renamed copies go to annotated_29nt_named/."""
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
import os, re
import pandas as pd
R = PEDV_P+"orf_calling_GIb/_price_final"
SC = [SR.lib("RIBO_VERO_CHX_1MOI_5h"), SR.lib("RIBO_VERO_CHX_1MOI_8h"), SR.lib("RIBO_VERO_CHX_1MOI_16h"), "merged_CHX"]
T = {s: pd.read_csv(f"{R}/annotated_29nt/{s}.plus.orfs.annotated.tsv", sep="\t") for s in SC}
M = pd.read_csv(f"{R}/annotated_29nt_resultlevel_merge_of_3_samples.tsv", sep="\t")


def group_of(cls, rel):
    m = re.search(r"overlapping ([A-Za-z0-9]+)", str(cls))
    if m and rel in (1, 2, "1", "2", 1.0, 2.0, "+1", "+2"):
        return m.group(1), f"+{int(float(str(rel).replace('+', '')))}"
    if str(cls).startswith("5'"): return "5'UTR", ""
    if str(cls).startswith("3'"): return "3'UTR", ""
    if str(cls).startswith("intergenic"): return "IG", ""
    return None, None


def prep(d):
    d = d.copy()
    d["start"] = d["our_start_1based"]; d["stop"] = d["our_stop_end_1based"]
    d["rel"] = d["frame_relative_to_parent"].astype(str).str.replace("+", "", regex=False).str.replace(".0", "", regex=False)
    g = [group_of(c, r) for c, r in zip(d.coordinate_class, d.rel)]
    d["grp"] = [a for a, b in g]; d["frm"] = [b for a, b in g]
    return d


P = {s: prep(t) for s, t in T.items()}
U = pd.concat([p[p.ORF_name.str.startswith("unnamed")][["grp", "frm", "start", "stop"]] for p in P.values()]).drop_duplicates()
U = U.sort_values(["grp", "frm", "start", "stop"])
U["k"] = U.groupby(["grp", "frm"]).cumcount() + 1
name = {}
for _, r in U.iterrows():
    name[(r.start, r.stop)] = f"{r.grp}{r.frm}#{r.k}" if r.frm else f"{r.grp}#{r.k}"
assert len(name) == len(U)
os.makedirs(f"{R}/annotated_29nt_named", exist_ok=True)
for s, p in P.items():
    d = T[s].copy()
    new = [name[(a, b)] if str(n).startswith("unnamed") else n for n, a, b in zip(d.ORF_name, p.start, p.stop)]
    d["ORF_name_old"] = d.ORF_name; d["ORF_name"] = new
    d.to_csv(f"{R}/annotated_29nt_named/{s}.plus.orfs.annotated.tsv", sep="\t", index=False)
    print(s, "renamed", int(d.ORF_name_old.str.startswith("unnamed").sum()))
# merged-of-3 table keyed by stop: names via the start of the most significant library
M2 = M.copy()
ns = []
for _, r in M2.iterrows():
    if str(r.ORF_name).startswith("unnamed"):
        cand = []
        for t in ("CHX5", "CHX8", "CHX16"):
            if pd.notna(r[f"start_{t}"]) and (int(r[f"start_{t}"]), int(r.stop)) in name: cand.append((r[f"p_{t}"], name[(int(r[f"start_{t}"]), int(r.stop))]))
        ns.append(sorted(cand)[0][1] if cand else r.ORF_name)
    else: ns.append(r.ORF_name)
M2["ORF_name_old"] = M2.ORF_name; M2["ORF_name"] = ns
M2.to_csv(f"{R}/annotated_29nt_resultlevel_merge_of_3_samples_named.tsv", sep="\t", index=False)
print("remaining unnamed:", int(M2.ORF_name.str.startswith("unnamed").sum()), "of", len(M2))
print(M2[M2["n_libs_p<0.05"] >= 2][["stop", "ORF_name", "ORF_name_old", "n_libs_p<0.05"]].to_string())
mp = U.assign(ORF_name=[name[(a, b)] for a, b in zip(U.start, U.stop)]); mp.to_csv(f"{R}/_29nt_work/unnamed_to_systematic_name_map.tsv", sep="\t", index=False)
