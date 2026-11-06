#!/usr/bin/env python3
"""Supplementary Fig.1E: RNA-seq PCA (host gene expression, DESeq2-style VST)."""
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
import os, warnings
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
warnings.filterwarnings("ignore")
from pydeseq2.dds import DeseqDataSet

D = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
F = f"{D}/figures"
SF = f"{D}/figdata/suppfig1"
man = pd.read_csv(f"{D}/sample_manifest_GIb.tsv", sep="\t")
man = man[man.sample_id.isin(SR.LIB_RNA_VERO)][["sample_id", "treatment", "time_h"]].copy()
man["__o"] = man.treatment.map({"CHX": 0, "HAR": 1})
man = man.sort_values(["__o", "time_h", "sample_id"]).reset_index(drop=True)
man["rep_id"] = man.groupby(["treatment", "time_h"]).cumcount() + 1
order = man.sample_id.tolist()


def load(fn):
    fc = pd.read_csv(f"{SF}/{fn}", sep="\t", comment="#")
    bam = [c for c in fc.columns if c.endswith(".bam")]
    fc = fc.rename(columns={c: "R260" + c.split("/R260")[1][:5] for c in bam})
    return fc


host = load("gene_counts.txt")
host = host[~host.Chr.astype(str).str.contains("PEDV_GIb")]
hc = host.set_index("Geneid")[order].astype(int)
hl = host.set_index("Geneid")["Length"]
vir = load("viral_orf_counts.txt")
vc = vir.set_index("Geneid")[order].astype(int)
vl = vir.set_index("Geneid")["Length"]


def tpm(c, l):
    rpk = c.div(l / 1000.0, axis=0)
    return rpk.div(rpk.sum(axis=0), axis=1) * 1e6


# expressed host genes: TPM>1 in >=1 sample (same rule as Fig.1D)
th = tpm(hc, hl)
expr = (th > 1).any(axis=1)
hc_f = hc[expr]
n_host = int(expr.sum())


def vst_pca(counts_df):
    meta = man.set_index("sample_id").loc[counts_df.columns, ["treatment"]]
    dds = DeseqDataSet(counts=counts_df.T, metadata=meta, design="~1", quiet=True)
    dds.fit_size_factors()
    dds.vst(use_design=False)
    vst = pd.DataFrame(dds.layers["vst_counts"], index=counts_df.columns, columns=counts_df.index)
    X = vst.values
    Xc = X - X.mean(0)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    var = (S ** 2) / (S ** 2).sum() * 100
    pcs = U * S
    return pcs, var, vst


pcs, var, _ = vst_pca(hc_f)
coord = man.copy()
for i in range(4):
    coord[f"PC{i+1}"] = pcs[:, i]
coord = coord[["sample_id", "treatment", "time_h", "rep_id", "PC1", "PC2", "PC3", "PC4"]]
coord.columns = ["sample_id", "condition", "time_hpi", "rep_id", "PC1", "PC2", "PC3", "PC4"]
coord.to_csv(f"{SF}/Supplementary_Fig1E_RNAseq_PCA_coordinates.tsv", sep="\t", index=False)
vardf = pd.DataFrame({"PC": [f"PC{i+1}" for i in range(len(var))], "variance_pct": var})
vardf.head(10).to_csv(f"{SF}/Supplementary_Fig1E_RNAseq_PCA_variance.tsv", sep="\t", index=False)

# secondary: host + viral
allc = pd.concat([hc, vc])
alll = pd.concat([hl, vl])
ta = tpm(allc, alll)
expr_a = (ta > 1).any(axis=1)
pcs_a, var_a, _ = vst_pca(allc[expr_a])
# align sign for comparison
pd.DataFrame({"PC": [f"PC{i+1}" for i in range(len(var_a))], "variance_pct": var_a}).head(10).to_csv(
    f"{SF}/Supplementary_Fig1E_RNAseq_PCA_variance_hostviral.tsv", sep="\t", index=False)
ca = man.copy()
for i in range(4):
    ca[f"PC{i+1}"] = pcs_a[:, i]

# ---- plot ----
plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "font.size": 8})
TIMES = [3, 5, 8, 12, 16]
cmap = plt.get_cmap("viridis")
TCOL = {t: cmap(i / (len(TIMES) - 1)) for i, t in enumerate(TIMES)}
SHAPE = {"CHX": "o", "HAR": "s"}


def scatter(ax, C, v1, v2, pcx=1, pcy=2):
    cm = man.copy()
    cm["PCx"] = C[f"PC{pcx}"].values if isinstance(C, pd.DataFrame) else C[:, pcx-1]
    # draw thin connectors between the two libraries of the same condition x time (reading aid only)
    merged = man.copy()
    merged["x"] = C[f"PC{pcx}"].values
    merged["y"] = C[f"PC{pcy}"].values
    for (_, _), g in merged.groupby(["treatment", "time_h"]):
        if len(g) == 2:
            ax.plot(g.x, g.y, "-", color="#CCCCCC", lw=0.6, zorder=1)
    for _, r in merged.iterrows():
        ax.scatter(r.x, r.y, s=42, marker=SHAPE[r.treatment], facecolor=TCOL[r.time_h],
                   edgecolor="#333333", linewidths=0.5, zorder=3)
    ax.set_xlabel(f"PC{pcx} ({v1[pcx-1]:.1f}% variance)", fontsize=9)
    ax.set_ylabel(f"PC{pcy} ({v1[pcy-1]:.1f}% variance)", fontsize=9)
    ax.axhline(0, color="#EEEEEE", lw=0.6, zorder=0)
    ax.axvline(0, color="#EEEEEE", lw=0.6, zorder=0)
    for s_ in ["top", "right"]:
        ax.spines[s_].set_visible(False)
    ax.tick_params(labelsize=8)


mm = 1 / 25.4
fig, ax = plt.subplots(figsize=(100 * mm, 82 * mm))
scatter(ax, coord, var, var)
th_handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=TCOL[t], markeredgecolor="#333",
                     markersize=7, label=f"{t} h") for t in TIMES]
sh_handles = [Line2D([0], [0], marker=SHAPE[d], color="none", markerfacecolor="#BBBBBB", markeredgecolor="#333",
                     markersize=7, label=d) for d in ["CHX", "HAR"]]
leg1 = ax.legend(handles=th_handles, title="hpi", loc="center left", bbox_to_anchor=(1.01, 0.72),
                 frameon=False, fontsize=7.5, title_fontsize=8, handletextpad=0.2)
ax.add_artist(leg1)
ax.legend(handles=sh_handles, title="Condition", loc="center left", bbox_to_anchor=(1.01, 0.25),
          frameon=False, fontsize=7.5, title_fontsize=8, handletextpad=0.2)
fig.tight_layout()
fig.savefig(f"{F}/Supplementary_Figure1E_RNAseq_PCA.png", dpi=300, facecolor="white", bbox_inches="tight")
fig.savefig(f"{F}/Supplementary_Figure1E_RNAseq_PCA.pdf", facecolor="white", bbox_inches="tight")
plt.close(fig)

# ---- descriptive checks ----
from scipy.stats import pearsonr, spearmanr
print("host expressed genes:", n_host)
print("variance PC1-4:", [round(v, 2) for v in var[:4]])
print("hostviral variance PC1-4:", [round(v, 2) for v in var_a[:4]])
for pc in ["PC1", "PC2", "PC3"]:
    r, _ = spearmanr(coord[pc], coord.time_hpi)
    print(f"{pc} vs time_hpi Spearman rho = {r:.3f}")
# CHX vs HAR separation along each PC (mean diff in sd units)
for pc in ["PC1", "PC2", "PC3"]:
    a = coord[coord.condition == "CHX"][pc]
    b = coord[coord.condition == "HAR"][pc]
    pooled = coord[pc].std()
    print(f"{pc} CHX-HAR mean diff = {(a.mean()-b.mean()):.2f} ({(a.mean()-b.mean())/pooled:.2f} sd)")
# host-only vs host+viral: procrustes-ish correlation of PC1/PC2
def align(a, b):
    s = np.sign((a * b).sum()) or 1
    return s * b
for pc in [0, 1]:
    rr = pearsonr(pcs[:, pc], align(pcs[:, pc], pcs_a[:, pc]))[0]
    print(f"host vs host+viral PC{pc+1} |r| = {abs(rr):.4f}")
# pair (same cond x time) distance vs overall spread
import itertools
pair_d = []
for (_, _), g in coord.groupby(["condition", "time_hpi"]):
    if len(g) == 2:
        p = g[["PC1", "PC2"]].values
        pair_d.append(np.linalg.norm(p[0] - p[1]))
alld = [np.linalg.norm(coord[["PC1","PC2"]].values[i]-coord[["PC1","PC2"]].values[j])
        for i, j in itertools.combinations(range(20), 2)]
print(f"same cond×time pair dist (PC1-2): mean {np.mean(pair_d):.1f} vs all-pairs mean {np.mean(alld):.1f}")
