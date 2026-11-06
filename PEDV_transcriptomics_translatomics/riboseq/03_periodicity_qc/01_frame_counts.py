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
import os, sys, json, gzip, math
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.dirname(HERE); TASK = os.path.dirname(OUT)
RGP = f"{TASK}/riboseq_GIb_mapping_psite"; SRV = f"{RGP}/_work/server/host_psite"
TAB = f"{OUT}/tables"; os.makedirs(TAB, exist_ok=True)
SAMPLES = [(SR.lib("RIBO_VERO_CHX_1MOI_5h"), "CHX", 5), (SR.lib("RIBO_VERO_CHX_1MOI_8h"), "CHX", 8), (SR.lib("RIBO_VERO_CHX_1MOI_16h"), "CHX", 16),
           (SR.lib("RIBO_VERO_HAR_1MOI_5h"), "HAR", 5), (SR.lib("RIBO_VERO_HAR_1MOI_8h"), "HAR", 8), (SR.lib("RIBO_VERO_HAR_1MOI_16h"), "HAR", 16)]
LS = list(range(20, 41))
rng = np.random.default_rng(20261012)
rd = lambda p: pd.read_csv(p, sep="\t", comment="#", low_memory=False)

def wilson(k, n, z=1.96):
    if n == 0: return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n; c = p + z * z / (2 * n); h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - h) / d, (c + h) / d)

off = rd(f"{RGP}/psite_offset_GIb.tsv")
rows = []; amb_in = {}
for sid, drug, t in SAMPLES:
    g = rd(f"{SRV}/{sid}/host_gene_cds_psite.tsv")
    mg = rd(f"{SRV}/{sid}/host_psite_metagene.tsv")
    qc = dict(l.rstrip().split("\t") for l in open(f"{SRV}/{sid}/host_psite_qc.tsv"))
    g["cds_mod3"] = g.cds_len_nt % 3

    tx_ids = g.transcript_id.unique(); tx_index = {t_: i for i, t_ in enumerate(tx_ids)}
    g["ti"] = g.transcript_id.map(tx_index)
    st = mg[(mg.anchor == "start")]; en = mg[(mg.anchor == "stop")]
    for L in LS:
        gl = g[g.L == L]
        N = np.array([gl.f0.sum(), gl.f1.sum(), gl.f2.sum()], dtype=np.int64)

        s = st[(st.L == L) & (st.rel >= 0) & (st.rel <= 44)]
        e = en[(en.L == L) & (en.rel >= -45) & (en.rel <= -1)]
        S = np.array([s[s.rel % 3 == f].n.sum() for f in range(3)], dtype=np.int64)
        E = np.array([e[e.rel % 3 == f].n.sum() for f in range(3)], dtype=np.int64)
        e30 = en[(en.L == L) & (en.rel >= -30) & (en.rel <= -1)]
        E30 = np.array([e30[e30.rel % 3 == f].n.sum() for f in range(3)], dtype=np.int64)
        body = N - S - E; body30 = N - S - E30

        bad = gl[gl.cds_mod3 != 0]; badfrac = (bad.f0.sum() + bad.f1.sum() + bad.f2.sum()) / max(N.sum(), 1)

        if N.sum() > 0 and len(gl) > 1:
            a = gl.f0.values.astype(float); tot = (gl.f0 + gl.f1 + gl.f2).values.astype(float)
            idx = rng.integers(0, len(gl), size=(1000, len(gl)))
            f0b = a[idx].sum(1) / np.maximum(tot[idx].sum(1), 1)
            lo_b, hi_b = np.percentile(f0b, [2.5, 97.5])
        else: lo_b = hi_b = np.nan
        for nm, X in (("full_CDS", N), ("CDS_body_45_45", body), ("CDS_body_45_30", body30)):
            n = int(X.sum()); f = X / n if n > 0 else np.array([np.nan] * 3)
            lo, hi = wilson(int(X[0]), n) if n > 0 else (np.nan, np.nan)
            o = off[(off.sample_id == sid) & (off.read_length == L)]
            rows.append(dict(sample=sid, treatment=drug, time_hpi=t, read_length=L, region=nm,
                             total_host_CDS_psites=n, N0=int(X[0]), N1=int(X[1]), N2=int(X[2]),
                             F0=f[0], F1=f[1], F2=f[2], F0_minus_one_third=f[0] - 1 / 3, F0_over_expected=f[0] * 3,
                             F0_CI95_wilson_low=lo, F0_CI95_wilson_high=hi,
                             F0_CI95_transcript_bootstrap_low=(lo_b if nm == "full_CDS" else np.nan),
                             F0_CI95_transcript_bootstrap_high=(hi_b if nm == "full_CDS" else np.nan),
                             top_to_second_frame_ratio=(np.sort(X)[::-1][0] / max(np.sort(X)[::-1][1], 1) if n > 0 else np.nan),
                             psite_offset_applied=(o.offset.iloc[0] if len(o) else np.nan), offset_confidence=(o.offset_confidence.iloc[0] if len(o) else ""),
                             length_set_category=(o.category_GIb.iloc[0] if len(o) else ""),
                             frac_psites_in_transcripts_CDSlen_not_mod3=badfrac, n_transcripts=int(len(gl))))
D = pd.DataFrame(rows)

for nm in D.region.unique():
    m = D.region == nm
    D.loc[m, "length_share_of_host_CDS_psites"] = D[m].groupby("sample").total_host_CDS_psites.transform(lambda x: x / x.sum())
D.to_pickle(f"{TAB}/_frame_counts.pkl")

# ---------------- QC ----------------
qc_lines = []
ok = ((D.N0 + D.N1 + D.N2) == D.total_host_CDS_psites).all(); qc_lines.append(f"N0+N1+N2 == total (all rows): {ok}")
v = D[D.total_host_CDS_psites > 0]; ok2 = np.allclose(v.F0 + v.F1 + v.F2, 1.0); qc_lines.append(f"F0+F1+F2 == 1 (rows with data): {ok2}")
neg = (D[["N0", "N1", "N2"]] < 0).any(axis=1).sum(); qc_lines.append(f"negative counts after window subtraction: {int(neg)}")
for ln in qc_lines: print(ln)
print("max frac_psites_in_transcripts_CDSlen_not_mod3:", D.frac_psites_in_transcripts_CDSlen_not_mod3.max().round(4))
pd.set_option("display.width", 250)
b = D[(D.region == "CDS_body_45_45") & D.read_length.between(25, 35)]
print(b.pivot(index="sample", columns="read_length", values="F0").round(2))
f = D[(D.region == "full_CDS") & D.read_length.between(25, 35)]
print(f.pivot(index="sample", columns="read_length", values="F0").round(2))
