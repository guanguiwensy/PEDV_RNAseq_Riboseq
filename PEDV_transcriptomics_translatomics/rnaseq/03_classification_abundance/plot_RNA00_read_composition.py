#!/usr/bin/env python3
"""RNA00 / RNA00S: RNA-seq reads组成(PEDV/宿主/其他)按hpi。输入 read_composition_counts.tsv(服务器统计)+sample_manifest_GIb.tsv"""
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
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os
D=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
F=f"{D}/figures"
cnt=pd.read_csv(f"{D}/figdata/read_composition_counts.tsv",sep="\t",
    names=["sample_id","total_pairs","mapped_pairs","pedv_only","host_only","cross","unmapped"])
man=pd.read_csv(f"{D}/sample_manifest_GIb.tsv",sep="\t")[["sample_id","time_h"]]
man=man[man.sample_id.isin(SR.LIB_RNA_VERO)]
df=cnt.merge(man,on="sample_id")
assert len(df)==20
df["pedv_pairs"]=df.pedv_only; df["host_pairs"]=df.host_only
df["other_pairs"]=df.total_pairs-df.pedv_pairs-df.host_pairs   # 未比对 + 跨类(若有)
df["pedv_frac"]=df.pedv_pairs/df.total_pairs; df["host_frac"]=df.host_pairs/df.total_pairs; df["other_frac"]=df.other_pairs/df.total_pairs
assert np.allclose(df.pedv_frac+df.host_frac+df.other_frac,1)
assert (df.other_pairs==df.unmapped+df.cross).all()
df=df.sort_values(["time_h","sample_id"]).reset_index(drop=True)
cols=["sample_id","time_h","total_pairs","pedv_pairs","host_pairs","other_pairs","pedv_frac","host_frac","other_frac"]
df[cols].to_csv(f"{F}/RNA00S_read_composition_by_sample.csv",index=False)
# 时间点汇总
rows=[]
for t,g in df.groupby("time_h"):
    r={"time_h":t,"n":len(g)}
    for k in ["pedv","host","other"]:
        r[f"{k}_mean"]=g[f"{k}_frac"].mean(); r[f"{k}_min"]=g[f"{k}_frac"].min(); r[f"{k}_max"]=g[f"{k}_frac"].max()
    rows.append(r)
agg=pd.DataFrame(rows); agg["sum_of_means"]=agg.pedv_mean+agg.host_mean+agg.other_mean
agg.to_csv(f"{F}/RNA00_read_composition_by_hpi.csv",index=False)
C={"pedv":"#0072B2","host":"#009E73","other":"#BDBDBD"}
LAB={"pedv":"PEDV","host":"Host (Vero)","other":"Other"}
plt.rcParams.update({"font.family":"DejaVu Sans","pdf.fonttype":42,"font.size":8})
def style(ax):
    for s in ["top","right"]: ax.spines[s].set_visible(False)
    ax.set_ylim(0,1); ax.set_yticks([0,.25,.5,.75,1]); ax.tick_params(labelsize=8)
    ax.set_ylabel("Relative fraction of reads",fontsize=10)
def legend(ax):
    ax.legend(handles=[Patch(color=C[k],label=LAB[k]) for k in ["pedv","host","other"]],loc="lower center",
        bbox_to_anchor=(0.5,1.0),ncol=3,frameon=False,fontsize=8,handlelength=1,columnspacing=1.2)
# 主图
fig,ax=plt.subplots(figsize=(88/25.4,70/25.4))
x=np.arange(len(agg)); b=np.zeros(len(agg))
for k in ["pedv","host","other"]:
    ax.bar(x,agg[f"{k}_mean"],bottom=b,color=C[k],width=0.72,edgecolor="white",linewidth=0.4); b+=agg[f"{k}_mean"].values
ax.set_xticks(x); ax.set_xticklabels(agg.time_h.astype(int)); ax.set_xlabel("Hours post infection (hpi)",fontsize=10)
style(ax); legend(ax); fig.tight_layout()
for e,kw in [("png",{"dpi":300}),("pdf",{})]: fig.savefig(f"{F}/RNA00_read_composition_by_hpi.{e}",facecolor="white",**kw)
plt.close(fig)
# 补图: 20柱
fig,ax=plt.subplots(figsize=(120/25.4,70/25.4))
xs=[];pos=0;centers=[]
for t,g in df.groupby("time_h",sort=True):
    p=[pos+i for i in range(len(g))]; xs+=p; centers.append(np.mean(p)); pos+=len(g)+1
xs=np.array(xs); b=np.zeros(len(df))
for k in ["pedv","host","other"]:
    ax.bar(xs,df[f"{k}_frac"],bottom=b,color=C[k],width=0.8,edgecolor="white",linewidth=0.3); b+=df[f"{k}_frac"].values
ax.set_xticks(centers); ax.set_xticklabels(sorted(df.time_h.unique().astype(int))); ax.set_xlim(-0.8,xs.max()+0.8)
ax.set_xlabel("Hours post infection (hpi)",fontsize=10); style(ax); legend(ax); fig.tight_layout()
for e,kw in [("png",{"dpi":300}),("pdf",{})]: fig.savefig(f"{F}/RNA00S_read_composition_by_sample.{e}",facecolor="white",**kw)
print(agg.round(4).to_string()); print(df[cols].round(4).to_string())
print("cross total",df.cross.sum())
