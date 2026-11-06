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
import os, sys, numpy as np, pandas as pd
HERE=os.path.dirname(os.path.abspath(__file__)); G=os.path.dirname(HERE)
ANN=os.path.join(os.path.dirname(G),"strain_assembly_phylo_v1","tables","annot")
DS={"LJX_merged":"LJXmerged","LLC-PK":"LLCPK"}
PAD=5
rows=[];regs=[];maps={}
for ds,pre in DS.items():
    lf=pd.read_csv(os.path.join(ANN,f"{pre}_lift_to_GIb_v1.1.tsv"),sep="\t",dtype=str)
    ix=pd.read_csv(os.path.join(ANN,f"{pre}_indels_vs_GIb_v1.1.tsv"),sep="\t")
    L=len(lf); q=lf.new_pos.astype(int).values
    g=np.array([int(x) if x not in ("-","") else 0 for x in lf.GIb_pos]); st=lf.status.values
    mm=(st=="mismatch").astype(int); cs=np.concatenate([[0],np.cumsum(mm)])
    amb=np.zeros(L+2,dtype=bool); indel_nt=np.zeros(L+2,dtype=bool)
    for r in ix.itertuples(index=False):
        p=int(r.new_pos_1based); n=int(r.length)
        lo=max(1,p-PAD); hi=min(L,p+(n if r.type=="ins_in_new" else 0)+PAD); amb[lo:hi+1]=True
        regs.append(dict(dataset=ds,type=r.type,gII_pos_start=p,gII_pos_end=p+n-1 if r.type=="ins_in_new" else p,GIb_pos_adjacent=int(r[2]),length_nt=n,
                         note="GII 相对 GI-b 插入: 插入碱基无 GI-b 对应(unmappable)" if r.type=="ins_in_new" else "GII 相对 GI-b 缺失 %d nt: GI-b 这些位点在 GII 无对应"%n))

    nz=np.where(g>0)[0]
    if nz[0]>0: regs.append(dict(dataset=ds,type="5prime_unaligned_GII",gII_pos_start=1,gII_pos_end=int(nz[0]),GIb_pos_adjacent=int(g[nz[0]]),length_nt=int(nz[0]),note="GII 5′端无 GI-b 对应"))
    if g[nz[0]]>1: regs.append(dict(dataset=ds,type="5prime_GIb_not_covered",gII_pos_start=int(nz[0]+1),gII_pos_end=int(nz[0]+1),GIb_pos_adjacent=1,length_nt=int(g[nz[0]]-1),note="GI-b 1..%d 在该 GII 组装中缺失(组装5′端不完整); 涉及 leader 5′端坐标=unmappable"%(g[nz[0]]-1)))
    if g[nz[-1]]<27959: regs.append(dict(dataset=ds,type="3prime_GIb_not_covered",gII_pos_start=int(nz[-1]+1),gII_pos_end=int(nz[-1]+1),GIb_pos_adjacent=int(g[nz[-1]]),length_nt=int(27959-g[nz[-1]]),note="GII 3′端比 GI-b 短"))

    for i in range(L):
        lo=max(0,i-10); hi=min(L,i+11); loc=int(cs[hi]-cs[lo])
        if st[i]=="match": cls="unique_match"
        elif st[i]=="mismatch": cls="SNP_unique_column"
        else: cls="unmappable_insertion_or_unaligned_end"
        if cls!="unmappable_insertion_or_unaligned_end" and amb[i+1]: cls="ambiguous_near_indel(±%d nt)"%PAD
        rows.append((ds,i+1,lf.new_base[i],int(g[i]) if g[i] else "",lf.GIb_base[i] if g[i] else "",st[i],cls,loc))
    maps[ds]=g

    w=50
    for s in range(0,L-w,25):
        n=int(cs[s+w]-cs[s])
        if n>=10: regs.append(dict(dataset=ds,type="divergent_window(>=10 SNP/50nt)",gII_pos_start=s+1,gII_pos_end=s+w,GIb_pos_adjacent=int(g[s]) if g[s] else 0,length_nt=w,note=f"{n} 个 SNP/{w} nt: 局部同源较低, 该窗口内 junction 断点映射标 caution"))
out=pd.DataFrame(rows,columns=["dataset","gII_pos","gII_base","GIb_pos_v1.1","GIb_base","raw_status","mapping_class","SNP_within_±10nt"])
out.to_csv(os.path.join(G,"tables","GII_to_GIb_coordinate_liftover.tsv"),sep="\t",index=False)
rg=pd.DataFrame(regs);

rg.to_csv(os.path.join(G,"tables","GII_to_GIb_liftover_regions.tsv"),sep="\t",index=False)
print(out.groupby(["dataset","mapping_class"]).size()); print(rg.groupby(["dataset","type"]).size())
