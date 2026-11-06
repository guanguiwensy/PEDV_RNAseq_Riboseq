# -*- coding: utf-8 -*-
"""
Whole-genome +1-frame ORF-capacity screen (>40 aa). ORFs startable at AUG or the project's
near-cognate set, extending in the +1 frame (relative to each canonical parent) to first
in-frame stop, length >40 aa (i.e. >=41). Describe with FROZEN HAR D (unchanged) + CHX body.
No visualization, no binary calling, no composite score. Includes ORF1a/1b.
Outputs: plus1_ORF_start_variants_gt40aa.tsv, plus1_ORF_families_gt40aa.tsv
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
import gzip, csv, os
import numpy as np
from collections import defaultdict
RGP=PEDV_P+r"riboseq_GIb_mapping_psite"
LONG=RGP+r"/psite_counts_PEDV_GIb.tsv.gz"
REF=PEDV_P+r"refbuild_GIb/v1.1/virus/PEDV_GIb.fa"
OUT=PEDV_P+r"global_perturbation_screen/plus1_orf_capacity/tables"; os.makedirs(OUT,exist_ok=True)
VLEN=27959; T=10; EPS=1.0; WINS=[15,30,45,60]; PRIMARY=30
START_SET={"ATG","CTG","GTG","TTG","ATA","ATT","ATC","ACG"}; STOP={"TAA","TAG","TGA"}
SAMP={SR.lib("RIBO_VERO_CHX_1MOI_5h"):("CHX",5),SR.lib("RIBO_VERO_HAR_1MOI_5h"):("HAR",5),SR.lib("RIBO_VERO_CHX_1MOI_8h"):("CHX",8),
      SR.lib("RIBO_VERO_HAR_1MOI_8h"):("HAR",8),SR.lib("RIBO_VERO_CHX_1MOI_16h"):("CHX",16),SR.lib("RIBO_VERO_HAR_1MOI_16h"):("HAR",16)}
CHX={5:SR.lib("RIBO_VERO_CHX_1MOI_5h"),8:SR.lib("RIBO_VERO_CHX_1MOI_8h"),16:SR.lib("RIBO_VERO_CHX_1MOI_16h")}; HAR={5:SR.lib("RIBO_VERO_HAR_1MOI_5h"),8:SR.lib("RIBO_VERO_HAR_1MOI_8h"),16:SR.lib("RIBO_VERO_HAR_1MOI_16h")}; TPS=[5,8,16]
CHX_IDS=list(CHX.values()); HAR_IDS=list(HAR.values())
ORFS=[("ORF1a",293,12622),("ORF1b",12592,20613),("S",20610,24758),("ORF3",24758,25033),
      ("E",25361,25600),("M",25608,26288),("N",26300,27625)]
PRF_SITE=12595  # ORF1a/1b -1 PRF slippery region (ORF1b CGA start 12592)

seq="".join(l.strip() for l in open(REF) if not l.startswith(">")).upper().replace("U","T")
def cod(p): return "NA" if (p<1 or p+2>VLEN) else seq[p-1:p+2]
def rna(c): return c.replace("T","U")

# P-site frozen strict_core +strand
counts={s:np.zeros(VLEN+2) for s in SAMP}; pos_len=defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
with gzip.open(LONG,"rt",encoding="utf-8") as fh:
    idx=None
    for line in fh:
        if line.startswith("#"): continue
        if line.startswith("sample_id"): h=line.rstrip().split("\t"); idx={k:i for i,k in enumerate(h)}; continue
        p=line.rstrip("\n").split("\t")
        if p[idx["length_category"]]!="strict_core" or p[idx["strand"]]!="+": continue
        sid=p[idx["sample_id"]]
        if sid in SAMP: counts[sid][int(p[idx["position"]])]+=int(p[idx["psite_count"]]); pos_len[sid][int(p[idx["position"]])][int(p[idx["read_length"]])]+=int(p[idx["psite_count"]])
depth={s:int(counts[s].sum()) for s in SAMP}
cum={s:{g:np.cumsum(np.where(((np.arange(VLEN+2)-1)%3==g),counts[s],0.0)) for g in (0,1,2)} for s in SAMP}
anyps=sum(counts[s] for s in SAMP)
def gf(p): return (p-1)%3
def fs(s,g,a,b):
    a=max(a,1); b=min(b,VLEN); return 0.0 if b<a else float(cum[s][g][b]-cum[s][g][a-1])
def Dtp(ccx,chr_,ocx,ohr): return float(np.log2((chr_+1)/(ccx+1))-np.log2((ohr+1)/(ocx+1)))

# per-parent: candidate frame c, other, canon; D at every +1 position (window30) for local null
parent_frames={}; parent_nullD={}
for name,s,e in ORFS:
    gfp=gf(s); c=(gfp+1)%3; o=(gfp+2)%3
    parent_frames[name]=(gfp,c,o,s,e)
    pos=[p for p in range(s,e+1) if gf(p)==c]; nd=[]
    for p in pos:
        we=min(p+PRIMARY,e)
        dv=[];
        for tp in TPS:
            ccx=fs(CHX[tp],c,p,we); chr_=fs(HAR[tp],c,p,we); ocx=fs(CHX[tp],o,p,we); ohr=fs(HAR[tp],o,p,we)
            if min(ccx,chr_)>=T: dv.append(Dtp(ccx,chr_,ocx,ohr))
        nd.append((p, float(np.median(dv)) if dv else np.nan))
    parent_nullD[name]=nd

def local_null(name,p,Dval):
    nd=parent_nullD[name]; arr=np.array([(q,d) for q,d in nd],dtype=float)
    q=arr[:,0]; d=arr[:,1]; mask=(np.abs(q-p)>60)&(~np.isnan(d)); bg=d[mask]
    if bg.size<10 or np.isnan(Dval): return ("","","","")
    med=np.median(bg); mad=np.median(np.abs(bg-med))
    return (round(float(med),4),round(float(mad),4),round(float((bg<Dval).mean()*100),2),round(float((bg>=Dval).mean()),4))

def first_stop(p,c):
    q=p
    while q+2<=VLEN:
        if (q-1)%3==c and cod(q) in STOP: return q
        q+=3 if (q-1)%3==c else 1
    return None

def dual_rl(positions):
    dual=[SR.lib("RIBO_VERO_CHX_1MOI_5h"),SR.lib("RIBO_VERO_CHX_1MOI_8h"),SR.lib("RIBO_VERO_HAR_1MOI_8h"),SR.lib("RIBO_VERO_HAR_1MOI_16h")]
    return sum(1 for s in dual if sum(pos_len[s][p].get(29,0) for p in positions)>0 and sum(pos_len[s][p].get(30,0) for p in positions)>0)

SV=[]; fam_map=defaultdict(list)
for name,(gfp,c,o,ps,pe) in parent_frames.items():
    canon=gfp
    for p in range(ps,pe+1):
        if gf(p)!=c or cod(p) not in START_SET: continue
        st=first_stop(p,c)
        if st is None: continue
        aa=(st-p)//3
        if not (aa>40): continue
        stopc=cod(st); orf_nt=st-p+3
        # windows D
        wvals={}
        for W in WINS:
            we=min(p+W,st-1); dvw=[]; gt=[]
            for tp in TPS:
                ccx=fs(CHX[tp],c,p,we); chr_=fs(HAR[tp],c,p,we); ocx=fs(CHX[tp],o,p,we); ohr=fs(HAR[tp],o,p,we)
                dvw.append((tp,Dtp(ccx,chr_,ocx,ohr),min(ccx,chr_)>=T,ccx,chr_,ocx,ohr))
            wvals[W]=dvw
        pw=wvals[PRIMARY]; Ds={tp:d for tp,d,g,a,b,cc,dd in pw}; gd={tp:g for tp,d,g,a,b,cc,dd in pw}
        gated=[Ds[tp] for tp in TPS if gd[tp]]
        Dmed=float(np.median(gated)) if gated else np.nan; Dmean=float(np.mean(gated)) if gated else np.nan
        sign="".join(("+" if Ds[tp]>0 else("-" if Ds[tp]<0 else"0")) if gd[tp] else "." for tp in TPS)
        ln_med,ln_mad,ln_pct,ln_tail=local_null(name,p,Dmed)
        # strongest HAR peak in [p, min(p+60,st)] candidate frame
        peakpos=[q for q in range(p,min(p+60,st)) if gf(q)==c]
        hpool={q:sum(counts[HAR[tp]][q] for tp in TPS) for q in peakpos}
        pk=max(hpool,key=lambda q:hpool[q]) if hpool else p
        p29=sum(pos_len[s][pk].get(29,0) for s in SAMP); p30=sum(pos_len[s][pk].get(30,0) for s in SAMP)
        # CHX body metrics over [p, st-1]
        bcx=sum(fs(s,c,p,st-1) for s in CHX_IDS); bo=sum(fs(s,o,p,st-1) for s in CHX_IDS); bcan=sum(fs(s,canon,p,st-1) for s in CHX_IDS)
        codons=list(range(p,st-1,3)); cocc=[sum(counts[s][q] for s in CHX_IDS)>0 for q in codons]
        occ=float(np.mean(cocc)) if cocc else 0.0; longest=run=0
        for b in cocc:
            run=run+1 if b else 0; longest=max(longest,run)
        bodycov=float(np.mean([anyps[q]>0 for q in range(p,st)])) if st>p else 0.0
        cpool={q:float(sum(counts[s][q] for s in CHX_IDS)) for q in codons}; ctot=sum(cpool.values()); cord=sorted(cpool,key=lambda q:-cpool[q])
        ct1=(cpool[cord[0]]/ctot) if (cord and ctot>0) else 0.0; ct3=(sum(cpool[q] for q in cord[:3])/ctot) if(cord and ctot>0) else 0.0
        def fp_ex(ex): cb=bcx-sum(counts[s][q] for s in CHX_IDS for q in ex); return round(float(np.log2((cb+1)/(bo+1))),4)
        inpos=[q for q in range(p,min(p+30,st-1)+1) if gf(q)==c]
        leak=round(sum(fs(s,canon,max(1,p-150),min(VLEN,p+150)) for s in CHX_IDS)/max(1,301)/max(1,sum(depth[s] for s in CHX_IDS))*1e6,3)
        dqa=("ORF1a/1b_low_periodicity(5/6 libs);assessable_note:not_assessable_quality" if name in("ORF1a","ORF1b") else "assessable_quality")
        rec=dict(family_id=f"P1_{name}_stop{st}",start_position=p,start_codon=rna(cod(p)),stop_position=st,stop_codon=rna(stopc),
            ORF_length_nt=orf_nt,ORF_length_aa=aa,parent_canonical_ORF=name,overlap_type="plus1_internal",relative_frame=f"{name}:+1",
            genomic_region=f"{p}-{st+2}",data_quality_annotation=dqa,
            overlaps_ORF1a=("yes" if not(st< 293 or p>12622) else "no"),overlaps_ORF1b=("yes" if not(st<12592 or p>20613) else "no"),
            near_known_PRF_region=("yes" if abs(p-PRF_SITE)<=100 else "no"),distance_to_known_PRF_site=p-PRF_SITE,
            D5=round(Ds[5],4),D8=round(Ds[8],4),D16=round(Ds[16],4),D_median=("" if np.isnan(Dmed) else round(Dmed,4)),
            D_mean=("" if np.isnan(Dmean) else round(Dmean,4)),sign_pattern=sign,n_gated=len(gated),n_positive=sum(1 for tp in TPS if gd[tp] and Ds[tp]>0),
            local_null_median=ln_med,local_null_MAD=ln_mad,local_percentile=ln_pct,empirical_tail_probability=ln_tail,
            strongest_HAR_peak_position=pk,peak_offset_from_start=pk-p,top1_fraction=round((hpool.get(pk,0)/max(1,sum(hpool.values()))),4) if hpool else 0,
            primary_window=PRIMARY)
        for W in WINS:
            for tp,d,g,a,b,cc,dd in wvals[W]:
                rec[f"w{W}_D{tp}"]=round(d,4)
        for tp in TPS:
            _,d,g,a,b,ccx,ohr=0,0,0,0,0,0,0
            ccx=fs(CHX[tp],c,p,min(p+PRIMARY,st-1)); chr_=fs(HAR[tp],c,p,min(p+PRIMARY,st-1)); ocx=fs(CHX[tp],o,p,min(p+PRIMARY,st-1)); ohr=fs(HAR[tp],o,p,min(p+PRIMARY,st-1))
            rec[f"candidate_CHX{tp}"]=int(ccx); rec[f"candidate_HAR{tp}"]=int(chr_); rec[f"other_CHX{tp}"]=round(ocx,1); rec[f"other_HAR{tp}"]=round(ohr,1)
            rec[f"canonical_CHX{tp}"]=int(fs(CHX[tp],canon,p,min(p+PRIMARY,st-1))); rec[f"canonical_HAR{tp}"]=int(fs(HAR[tp],canon,p,min(p+PRIMARY,st-1)))
        rec.update(peak_candidate_count=int(hpool.get(pk,0)),peak_other_count=int(sum(counts[HAR[tp]][pk] for tp in TPS) and sum(fs(HAR[tp],o,pk,pk) for tp in TPS)) if False else int(sum(fs(HAR[tp],o,pk,pk) for tp in TPS)),
            peak_29nt_fraction=round(p29/max(1,p29+p30),3),peak_30nt_fraction=round(p30/max(1,p29+p30),3),
            CHX_candidate_body_count=int(bcx),CHX_other_body_count=int(bo),CHX_canonical_body_count=int(bcan),
            CHX_frame_preference=round(float(np.log2((bcx+1)/(bo+1))),4),occupied_codons=int(np.sum(cocc)),evaluable_codons=len(codons),
            occupancy_fraction=round(occ,4),longest_inframe_run=longest,body_coverage=round(bodycov,4),
            CHX_top1_fraction=round(ct1,4),CHX_top3_fraction=round(ct3,4),CHX_de_top1_frame_preference=(fp_ex(cord[:1]) if cord else ""),
            CHX_de_top3_frame_preference=(fp_ex(cord[:3]) if cord else ""),readlength_robustness=f"{dual_rl(inpos)}/4",leakage_level=leak,
            detectability_status=("below_count_floor" if not gated else "evaluable"))
        SV.append(rec); fam_map[rec["family_id"]].append(rec)
print("start variants (>40aa):",len(SV),"families:",len(fam_map))

# family table
def nf(x):
    try: return np.nan if x in("","nan") else float(x)
    except: return np.nan
FAM=[]
for fid,vs in fam_map.items():
    Ds=[(v["start_position"],nf(v["D_median"])) for v in vs if nf(v["D_median"])==nf(v["D_median"])]
    lps=[nf(v["local_percentile"]) for v in vs if nf(v["local_percentile"])==nf(v["local_percentile"])]
    best=max(Ds,key=lambda x:x[1]) if Ds else (None,np.nan)
    brow=[v for v in vs if v["start_position"]==best[0]] if best[0] else []
    brow=brow[0] if brow else vs[0]
    peaks=set(v["strongest_HAR_peak_position"] for v in vs)
    FAM.append(dict(family_id=fid,parent_ORF=vs[0]["parent_canonical_ORF"],overlap_type=vs[0]["overlap_type"],
        shared_stop=vs[0]["stop_position"],n_start_variants=len(vs),all_candidate_starts="|".join(str(v["start_position"]) for v in sorted(vs,key=lambda x:x["start_position"])),
        all_start_codons="|".join(sorted(set(v["start_codon"] for v in vs))),
        max_ORF_length_aa=max(v["ORF_length_aa"] for v in vs),min_ORF_length_aa=min(v["ORF_length_aa"] for v in vs),
        best_D_start=best[0] if best[0] else "",best_D_start_codon=brow["start_codon"],max_D=("" if np.isnan(best[1]) else round(best[1],4)),
        family_median_D=("" if not Ds else round(float(np.median([d for _,d in Ds])),4)),best_local_percentile=(max(lps) if lps else ""),
        n_positive_starts=sum(1 for _,d in Ds if d>0),n_cross_time_supported_starts=sum(1 for v in vs if int(v["n_gated"])>=2 and int(v["n_positive"])==int(v["n_gated"]) and int(v["n_positive"])>=2),
        shared_peak_position=("|".join(str(x) for x in sorted(peaks)) if len(peaks)>1 else list(peaks)[0]),
        shared_HAR_peak_yes_no=("yes" if len(peaks)==1 and len(vs)>1 else "no"),spread_of_supported_starts=max(v["start_position"] for v in vs)-min(v["start_position"] for v in vs),
        HAR_response_architecture=("focal" if nf(brow["top1_fraction"])>0.5 else("intermediate" if nf(brow["top1_fraction"])>0.25 else "distributed")),
        CHX_frame_preference=brow["CHX_frame_preference"],occupancy_fraction=brow["occupancy_fraction"],longest_inframe_run=brow["longest_inframe_run"],
        body_coverage=brow["body_coverage"],leakage_level=brow["leakage_level"],detectability_status=brow["detectability_status"],
        data_quality_annotation=brow["data_quality_annotation"],notes=""))
# ranks
def rank_by(key,rows,rev=True):
    vals=sorted(rows,key=lambda r:(nf(r[key]) if nf(r[key])==nf(r[key]) else -1e9),reverse=rev)
    for i,r in enumerate(vals): r["rank_by_"+key]=i+1
rank_by("best_local_percentile",FAM); rank_by("max_D",FAM)
for r in FAM: r["rank_by_local_response"]=r.pop("rank_by_best_local_percentile"); r["rank_by_D"]=r.pop("rank_by_max_D")

def wr(fn,rows,cols,note):
    with open(OUT+"/"+fn,"w",encoding="utf-8",newline="") as f:
        f.write(note+"\n"); w=csv.DictWriter(f,fieldnames=cols,delimiter="\t",restval=""); w.writeheader()
        w.writerows(rows)
    print(" ",fn,len(rows),"rows")
svcols=list(dict.fromkeys(k for r in SV for k in r))
NOTE="# whole-genome +1-frame ORF-capacity screen (>40aa; start=AUG/near-cognate{CUG,GUG,UUG,AUA,AUU,AUC,ACG}); FROZEN HAR D & strict_core; no score/calling"
# sort family by local percentile desc, max_D desc, cross-time desc, aa desc
FAMs=sorted(FAM,key=lambda r:(-(nf(r["best_local_percentile"]) if nf(r["best_local_percentile"])==nf(r["best_local_percentile"]) else -1),
    -(nf(r["max_D"]) if nf(r["max_D"])==nf(r["max_D"]) else -1e9),-int(r["n_cross_time_supported_starts"]),-int(r["max_ORF_length_aa"])))
famcols=list(dict.fromkeys(k for r in FAM for k in r))
wr("plus1_ORF_start_variants_gt40aa.tsv",sorted(SV,key=lambda r:r["start_position"]),svcols,NOTE)
wr("plus1_ORF_families_gt40aa.tsv",FAMs,famcols,NOTE)

# ===== report numbers =====
from collections import Counter
print("\n== SUMMARY ==")
print("families:",len(FAM)," start variants:",len(SV))
print("families by parent:",dict(Counter(r["parent_ORF"] for r in FAM)))
print("start variants by codon:",dict(Counter(r["start_codon"] for r in SV)))
aug=sum(1 for r in SV if r["start_codon"]=="AUG"); print(f"AUG starts: {aug}; near-cognate starts: {len(SV)-aug}")
print("families AUG-capable vs near-cognate-only:",
      sum(1 for r in FAM if "AUG" in r["all_start_codons"]),"AUG-capable;",sum(1 for r in FAM if "AUG" not in r["all_start_codons"]),"near-cognate-only")
fe=[r for r in FAM if r["best_local_percentile"]!="" ]
for thr in (90,95,99): print(f"families best_local_percentile>={thr}: {sum(1 for r in fe if nf(r['best_local_percentile'])>=thr)}")
print("families with >=1 start, >=2 timepoints D>0 (n_positive_starts>0 & cross-time):",sum(1 for r in FAM if int(r["n_cross_time_supported_starts"])>=1))
print("families with longCHX (longest_run>=15 & occ>=0.6):",sum(1 for r in FAM if nf(r["longest_inframe_run"])>=15 and nf(r["occupancy_fraction"])>=0.6))
