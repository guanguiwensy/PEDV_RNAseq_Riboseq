#!/usr/bin/env python3
"""Figure 1 B/C/D/F 按 MOI 重算的数据层(复用 GIb_rnaseq_reanalysis 的既有口径,不改任何既有文件)。

复用来源:
  - plot_Fig1DEF.py      : 唯一比对(MAPQ=255)fragment 计数、窗口/区域定义(refdata_GIb_v1.1 坐标)、未去重(raw)为主
  - plot_RNA00_read_composition.py : Panel B 的 PEDV/Host/Other 占比(分母=全部 QC 后 read pairs)
  - plot_figures_GIb.py  : Panel C 的 depth/总depth*1e6、101 nt 滑动平均
输入只读: GIb_rnaseq_reanalysis/{figdata,sample_manifest_GIb.tsv,figures}, <PEDV_P>/sample_info.txt
输出: Fig1_MOI_regrouped/data/*.csv
"""
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
import numpy as np
import pandas as pd

ROOT = PEDV_P+r""
SRC = f"{ROOT}/GIb_rnaseq_reanalysis"
OUT = f"{ROOT}/Fig1_MOI_regrouped"
DATA = f"{OUT}/data"
TIMES = [3, 5, 8, 12, 16]
GENOME_LEN = 27959

# 前缀 -> MOI。仅用于生成 metadata 表,图/源数据里不出现这些前缀。
PREFIX_MOI = {"Cv777": 0.1, "LZ": 1.0}


def load_metadata():
    si = pd.read_csv(f"{ROOT}/sample_info.txt", sep="\t")
    si = si[(si.assay == "RNA-seq") & si.production_id.between(SR.lib("RNA_VERO_CHX_0.1MOI_3h"), SR.lib("RNA_VERO_HAR_1MOI_16h"))].copy()
    assert len(si) == 20
    si["MOI"] = si.name_prefix.map(PREFIX_MOI)
    assert si.MOI.notna().all()
    meta = pd.DataFrame({
        "sample_id": si.production_id, "original_name": si.sample_name, "library_type": "RNA-seq",
        "virus_genotype": "GI-b", "MOI": si.MOI, "treatment": si.name_tag, "hpi": si.time_h.astype(int),
        "replicate": "NA", "has_biological_replicate": "no",
        "design_note": "full factorial 2 MOI x 2 treatment x 5 hpi, n=1 per cell"})
    # 与任务书给出的样本对照逐一核对
    moi01 = {SR.lib("RNA_VERO_CHX_0.1MOI_3h"), SR.lib("RNA_VERO_HAR_0.1MOI_3h"), SR.lib("RNA_VERO_CHX_0.1MOI_5h"), SR.lib("RNA_VERO_HAR_0.1MOI_5h"), SR.lib("RNA_VERO_CHX_0.1MOI_8h"), SR.lib("RNA_VERO_HAR_0.1MOI_8h"),
             SR.lib("RNA_VERO_CHX_0.1MOI_12h"), SR.lib("RNA_VERO_HAR_0.1MOI_12h"), SR.lib("RNA_VERO_CHX_0.1MOI_16h"), SR.lib("RNA_VERO_HAR_0.1MOI_16h")}
    chk = set(meta[meta.MOI == 0.1].sample_id)
    assert chk == moi01, chk ^ moi01
    assert (meta.groupby(["MOI", "treatment", "hpi"]).size() == 1).all()
    return meta.sort_values(["MOI", "hpi", "treatment"]).reset_index(drop=True)


REG = {"Leader": (10, 72, []), "ORF1a": (293, 12591, [(3298, 3499)]), "ORF1b": (12592, 20609, []),
       "S": (20610, 24757, []), "ORF3/E/M": (24758, 26288, [(24901, 25100)]), "N/3′UTR": (26300, 27959, [])}
CWIN = {"ORF1a": (1101, 12591), "ORF1b": (12592, 20609), "S": (20610, 24757), "N": (26300, 27625)}


def region_mean(dep, a, b, ex):
    m = np.zeros(len(dep) + 1, bool)
    m[a:b + 1] = True
    for x, y in ex:
        m[x:y + 1] = False
    return dep[m[1:]].mean()


def smooth(arr, win=101):
    return np.convolve(arr, np.ones(win) / win, mode="same")


def per_sample_table():
    meta = load_metadata()
    cnt = pd.read_csv(f"{SRC}/figdata/fig1_counts.tsv", sep="\t",
                      names=["sample_id", "ver", "total_mapped", "pedv_mapped", "r1_fwd", "r1_rev"])
    comp = pd.read_csv(f"{SRC}/figdata/read_composition_counts.tsv", sep="\t",
                       names=["sample_id", "total_pairs", "mapped_pairs", "pedv_only", "host_only", "cross", "unmapped"])
    rows = []
    for _, m in meta.iterrows():
        s = m.sample_id
        r = {"sample_id": s, "MOI": m.MOI, "treatment": m.treatment, "hpi": m.hpi}
        c = comp[comp.sample_id == s].iloc[0]
        r.update(B_total_QC_read_pairs=int(c.total_pairs), B_pedv_pairs=int(c.pedv_only), B_host_pairs=int(c.host_only),
                 B_other_pairs=int(c.total_pairs - c.pedv_only - c.host_only))
        r["B_pedv_frac"] = r.B_pedv_pairs / c.total_pairs if False else r["B_pedv_pairs"] / c.total_pairs
        r["B_host_frac"] = r["B_host_pairs"] / c.total_pairs
        r["B_other_frac"] = r["B_other_pairs"] / c.total_pairs
        for ver in ["raw", "dedup"]:
            k = cnt[(cnt.sample_id == s) & (cnt.ver == ver)].iloc[0]
            tag = "" if ver == "raw" else "_dedup"
            r[f"D_total_mapped_fragments{tag}"] = int(k.total_mapped)
            r[f"D_pedv_mapped_fragments{tag}"] = int(k.pedv_mapped)
            r[f"D_pedv_pct_of_mapped{tag}"] = 100 * k.pedv_mapped / k.total_mapped
            dep = pd.read_csv(f"{SRC}/figdata/fig1_depth/{s}.{ver}.depth.tsv", sep="\t", header=None)[2].values.astype(float)
            assert len(dep) == GENOME_LEN
            for g, (a, b, ex) in REG.items():
                v = region_mean(dep, a, b, ex) / (k.total_mapped / 1e6)
                r[f"F_{g}_cov_per_M_mapped{tag}"] = v
            if ver == "raw":
                r["D_r1_fwd_frac"] = k.r1_fwd / (k.r1_fwd + k.r1_rev)
            rpm = smooth(dep / dep.sum() * 1e6)
            for g, (a, b) in CWIN.items():
                r[f"C_{g}_median_RPM{tag}"] = float(np.median(rpm[a - 1:b]))
            r[f"C_peak_RPM{tag}"] = float(rpm.max())
            r[f"C_peak_pos{tag}"] = int(rpm.argmax() + 1)
        rows.append(r)
    return meta, pd.DataFrame(rows)


if __name__ == "__main__":
    os.makedirs(DATA, exist_ok=True)
    meta, S = per_sample_table()
    meta.to_csv(f"{DATA}/PEDV_sample_metadata_MOI_corrected_RNAseq20.tsv", sep="\t", index=False)
    S.to_csv(f"{DATA}/Fig1_per_sample_values_MOI.csv", index=False)
    print(S[["sample_id", "MOI", "treatment", "hpi", "B_pedv_frac", "D_pedv_pct_of_mapped"]].round(4).to_string())
