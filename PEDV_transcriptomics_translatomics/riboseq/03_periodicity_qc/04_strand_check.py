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
import os, gzip, re
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.dirname(HERE); TASK = os.path.dirname(OUT)
GTF = os.path.join(os.path.dirname(TASK), "refdata_vero", "Chlorocebus_sabaeus.ChlSab1.1.116.gtf.gz")
SRV = f"{TASK}/riboseq_GIb_mapping_psite/_work/server/host_psite"
rx = re.compile(r'(\S+) "([^"]*)"'); strand = {}
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line.startswith("#"): continue
        f = line.split("\t")
        if f[2] == "CDS": strand[dict(rx.findall(f[8]))["transcript_id"]] = f[6]
rows = []
for sid in [SR.lib("RIBO_VERO_CHX_1MOI_5h"), SR.lib("RIBO_VERO_CHX_1MOI_8h"), SR.lib("RIBO_VERO_CHX_1MOI_16h"), SR.lib("RIBO_VERO_HAR_1MOI_5h"), SR.lib("RIBO_VERO_HAR_1MOI_8h"), SR.lib("RIBO_VERO_HAR_1MOI_16h")]:
    g = pd.read_csv(f"{SRV}/{sid}/host_gene_cds_psite.tsv", sep="\t"); g["strand"] = g.transcript_id.map(strand)
    for L in (29, 30):
        a = g[g.L == L]
        for st, b in a.groupby("strand"):
            rows.append(dict(sample=sid, L=L, strand=st, n=int(b.f0.sum() + b.f1.sum() + b.f2.sum()), F0=round(b.f0.sum() / (b.f0.sum() + b.f1.sum() + b.f2.sum()), 3)))
R = pd.DataFrame(R if False else rows)
R.to_csv(f"{OUT}/tables/strand_consistency_F0.csv", index=False)
print(R.pivot_table(index=["sample", "L"], columns="strand", values="F0").to_string())
