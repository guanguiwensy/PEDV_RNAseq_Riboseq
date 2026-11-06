#!/usr/bin/env python3
"""Supplementary Fig.1D: RNA-seq sample correlation (host gene expression)."""
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
import sample_registry as SR
import os
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

D = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
F = f"{D}/figures"
SF = f"{D}/figdata/suppfig1"
man = pd.read_csv(f"{D}/sample_manifest_GIb.tsv", sep="\t")
man = man[man.sample_id.isin(SR.LIB_RNA_VERO)][["sample_id", "treatment", "time_h"]]
man["__o"] = man.treatment.map({"CHX": 0, "HAR": 1})
man = man.sort_values(["__o", "time_h", "sample_id"]).reset_index(drop=True)
order = man.sample_id.tolist()
label = {r.sample_id: f"{r.treatment} {int(r.time_h)}h ({r.sample_id[-4:]})" for _, r in man.iterrows()}

# ---- featureCounts matrix ----
fc = pd.read_csv(f"{SF}/gene_counts.txt", sep="\t", comment="#")
bamcols = [c for c in fc.columns if c.endswith(".bam")]
samp = {c: "R260" + c.split("/R260")[1][:5] for c in bamcols}
fc = fc.rename(columns=samp)
cnt = fc.set_index("Geneid")[order]
length = fc.set_index("Geneid")["Length"]
is_viral = fc.set_index("Geneid")["Chr"].astype(str).str.contains("PEDV_GIb")
host_genes = is_viral[~is_viral].index
viral_genes = is_viral[is_viral].index


def tpm(counts, lengths):
    rpk = counts.div(lengths / 1000.0, axis=0)
    return rpk.div(rpk.sum(axis=0), axis=1) * 1e6


# host-only TPM (normalised within host genes so viral load does not shift host TPM)
tpm_host = tpm(cnt.loc[host_genes], length.loc[host_genes])
keep = (tpm_host > 1).any(axis=1)
tpm_host_f = tpm_host[keep]
log_host = np.log2(tpm_host_f + 1)
# all-feature TPM (host + viral) for viral-load dominance check
tpm_all = tpm(cnt, length)
keep_all = (tpm_all > 1).any(axis=1)
log_all = np.log2(tpm_all[keep_all] + 1)

pear_host = log_host.corr(method="pearson")
spear_host = log_host.corr(method="spearman")
pear_all = log_all.corr(method="pearson")

pear_host.to_csv(f"{SF}/Supplementary_Fig1D_RNAseq_Pearson_correlation.tsv", sep="\t")
spear_host.to_csv(f"{SF}/Supplementary_Fig1D_RNAseq_Spearman_correlation.tsv", sep="\t")
pear_all.to_csv(f"{SF}/Supplementary_Fig1D_RNAseq_Pearson_allfeatures.tsv", sep="\t")

# ---- clustering (reported, not used for main order) ----
d_ = 1 - pear_host.values
np.fill_diagonal(d_, 0)
Z = linkage(squareform(d_, checks=False), method="average")
clust_order = [pear_host.index[i] for i in leaves_list(Z)]

# ---- heatmap: fixed temporal order grouped by CHX/HAR ----
plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "font.size": 8})
CM = LinearSegmentedColormap.from_list("corr", ["#F7F7F7", "#FDDBC7", "#D6604D", "#B2182B"])
M = pear_host.loc[order, order].values
vmin = np.floor(M[~np.eye(len(order), dtype=bool)].min() * 100) / 100
mm = 1 / 25.4
fig, ax = plt.subplots(figsize=(130 * mm, 118 * mm))
im = ax.imshow(M, cmap=CM, vmin=vmin, vmax=1, aspect="equal")
labs = [label[s] for s in order]
ax.set_xticks(range(len(order)))
ax.set_yticks(range(len(order)))
ax.set_xticklabels(labs, rotation=90, fontsize=7)
ax.set_yticklabels(labs, fontsize=7)
ax.tick_params(length=2)
# 相关系数不逐格标注(20x20下字号会<7pt)，颜色已编码r，完整数值见源数据TSV
# CHX/HAR分隔
nchx = (man.treatment == "CHX").sum()
for pos in [nchx - 0.5]:
    ax.axhline(pos, color="#333333", lw=0.9)
    ax.axvline(pos, color="#333333", lw=0.9)
cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
cb.set_label("Pearson r\n(host log2(TPM+1))", fontsize=8)
cb.ax.tick_params(labelsize=7)
fig.tight_layout()
fig.savefig(f"{F}/Supplementary_Figure1D_RNAseq_correlation.png", dpi=300, facecolor="white")
fig.savefig(f"{F}/Supplementary_Figure1D_RNAseq_correlation.pdf", facecolor="white")
plt.close(fig)

# ---- summaries for report ----
off = M[~np.eye(len(order), dtype=bool)]
print("host genes kept:", keep.sum(), "of", len(host_genes))
print("viral genes:", len(viral_genes))
print("Pearson off-diagonal range: %.3f - %.3f (median %.3f)" % (off.min(), off.max(), np.median(off)))
# mean r to others per sample (outlier check)
mr = (pear_host.sum(axis=1) - 1) / (len(order) - 1)
print("mean r to others, lowest 3:")
print(mr.sort_values().head(3).round(3).to_string())
# host vs all-feature
pa = pear_all.loc[order, order].values
offa = pa[~np.eye(len(order), dtype=bool)]
print("all-feature Pearson off-diag range: %.3f - %.3f (median %.3f)" % (offa.min(), offa.max(), np.median(offa)))
print("mean |host - allfeature| off-diag: %.4f" % np.mean(np.abs(off - offa)))
print("clustering order:", [s[-4:] for s in clust_order])
# adjacent-time similarity within drug
for drug in ["CHX", "HAR"]:
    ss = man[man.treatment == drug].sample_id.tolist()
    sub = pear_host.loc[ss, ss].values
    print(f"{drug} within-drug off-diag mean r: %.3f" % sub[~np.eye(len(ss), dtype=bool)].mean())
# CHX vs HAR same-time
chx = man[man.treatment == "CHX"].set_index("time_h").sample_id
har = man[man.treatment == "HAR"].set_index("time_h").sample_id
for t in [3, 5, 8, 12, 16]:
    rs = [pear_host.loc[a, b] for a in chx.loc[[t]] for b in har.loc[[t]]]
    print(f"CHX-HAR {t}h mean r: %.3f" % np.mean(rs))
