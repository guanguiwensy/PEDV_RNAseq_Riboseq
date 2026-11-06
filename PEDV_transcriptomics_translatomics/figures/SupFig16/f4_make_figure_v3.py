"""Supplementary Figure L REVISED v3 layout: page 1 = A (full width) / B + C3 (side by side) / D (+D2); page 2 = C1 over C2 (full width)."""
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
import csv,json,collections,textwrap,numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle,Polygon
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
B=PEDV_K+r""; O=B+"Supplementary_Figure_L/"; W=O+"work/"
plt.rcParams.update({"font.family":["Arial","DejaVu Sans"],"font.size":7,"axes.linewidth":0.6,"svg.fonttype":"none","pdf.fonttype":42,
    "figure.facecolor":"white","axes.facecolor":"white","xtick.major.width":0.5,"ytick.major.width":0.5,"xtick.major.size":2.5,"ytick.major.size":2.5})
RCOL={"ORF1a":"#0072B2","ORF1b":"#56B4E9","S":"#009E73","ORF3":"#CC79A7","E":"#D55E00","M":"#E69F00","N":"#117733"}
RED="#c1121f"; AUGC="#5e3c99"; NCC="#a6761d"
CDS=[("ORF1a",293,12622),("ORF1b",12592,20613),("S",20610,24758),("ORF3",24758,25033),("E",25361,25600),("M",25608,26288),("N",26300,27625)]
VLEN=27959; FW,FH=8.27,11.69; LM=0.45; WID=7.4
def rd(fn): return list(csv.DictReader(open(O+fn,encoding='utf-8-sig')))
A=rd("SourceData_PanelA_families_gt40.csv"); dens=rd("SourceData_PanelA_density_tracks.csv")
Bd=rd("SourceData_PanelB_region_density.csv"); Bs=rd("SourceData_PanelB_starts_per_family.csv"); Bst=json.load(open(W+"panelB_stats.json"))
feat=json.load(open(W+"featured.json")); NAME=feat['names']; C12=feat['C12_list']
nullr=list(csv.DictReader(open(W+"synonymous_null_by_family.csv",encoding='utf-8-sig')))
pg=rd("SourceData_PanelD_per_genome_completeness.csv"); Dm=rd("SourceData_PanelD_lineage_family_matrix.csv")
typ=list(csv.DictReader(open(W+"typing_S_based.tsv",encoding='utf-8'),delimiter='\t'))
lin=rd("PEDV_ORF_conservation_by_lineage.csv"); L={(r['family_id'],r['lineage_group']):r for r in lin}
import importlib.util
_sp=importlib.util.spec_from_file_location("f3",__import__('os').path.join(__import__('os').path.dirname(__REPO_FILE__),'f3_panelC.py')); f3=importlib.util.module_from_spec(_sp); _sp.loader.exec_module(f3)
def AX(fig,x,ytop,w,h,**k): return fig.add_axes([x/FW,1-(ytop+h)/FH,w/FW,h/FH],**k)
def TX(fig,x,y,s,**k): return fig.text(x/FW,1-y/FH,s,**k)
nc=[a for a in A if a['canonical_frame_family']=='no']
PP=PdfPages(O+"Supplementary_Figure_L_REVISED_v3.pdf")
# ============================================================ PAGE 1
fig=plt.figure(figsize=(FW,FH))
TX(fig,LM,0.12,"Supplementary Figure L  Sequence-defined candidate ORF families of PEDV and their retention across lineages",fontsize=9,fontweight="bold",va="top")
TX(fig,LM,0.42,"A",fontsize=10,fontweight="bold",va="top")
left=0.95; wA=WID-left+0.0      # axes x start (leaves room for frame labels)
axc=AX(fig,LM+left,0.62,wA,0.28); axc.set_xlim(0,VLEN); axc.set_ylim(0,1); axc.axis("off")
for n,s,e in CDS:
    axc.add_patch(Rectangle((s,0.1),e-s,0.8,facecolor=RCOL[n],edgecolor="none"))
    if e-s>1500: axc.text((s+e)/2,0.5,n,ha="center",va="center",color="white",fontsize=7,fontweight="bold")
    else: axc.text((s+e)/2,1.05,n,ha="center",va="bottom",fontsize=7,color=RCOL[n])
TX(fig,LM+left+0.1,0.97,"ORF1a/1b overlap 12,592–12,622; slippery heptamer 12,586–12,592 (known −1 PRF, not a new ORF)",fontsize=7,va="top",color="#444")
hl=[Rectangle((0,0),1,1,facecolor=AUGC,label="family contains AUG start"),Rectangle((0,0),1,1,facecolor=NCC,label="near-cognate starts only"),
    Rectangle((0,0),1,1,facecolor="white",edgecolor=RED,lw=1.1,label="Figure 6 featured (original, provisional)")]
fig.legend(handles=hl,loc="upper right",bbox_to_anchor=((LM+WID)/FW,1-0.30/FH),ncol=3,frameon=False,fontsize=7,handlelength=1.2,columnspacing=1.0)
lanes={}
for f in (0,1,2):
    fs=sorted([a for a in nc if int(a['genomic_frame_static'])==f],key=lambda a:int(a['start_min']))
    rows=[];pos={}
    for a in fs:
        s,e=int(a['start_min']),int(a['stop_position'])+2
        for i,end in enumerate(rows):
            if s>end+120: rows[i]=e; pos[a['family_id']]=i; break
        else: rows.append(e); pos[a['family_id']]=len(rows)-1
    lanes[f]=(len(rows),pos)
unit=0.19; y=1.2
FLAB={0:"frame 0"+chr(10)+"canonical: ORF1b"+chr(10)+"+1: S, M",1:"frame 1"+chr(10)+"canonical: ORF1a,"+chr(10)+"ORF3, E, N; +1: ORF1b",2:"frame 2"+chr(10)+"canonical: S, M"+chr(10)+"+1: ORF1a, ORF3, E, N"}
for f in (0,1,2):
    n,pos=lanes[f]; h=(n+2)*unit
    ax=AX(fig,LM+left,y,wA,h); ax.set_xlim(0,VLEN); ax.set_ylim(-0.5,n+1.5); ax.axis("off")
    ax.plot([0,VLEN],[-0.5,-0.5],color="#bbb",lw=0.5)
    cnt=0
    for a in sorted(nc,key=lambda a:int(a['start_min'])):
        if int(a['genomic_frame_static'])!=f: continue
        r=pos[a['family_id']]; s=int(a['start_min']); e=int(a['stop_position'])+2
        col=AUGC if a['contains_AUG']=='yes' else NCC
        isf=a['fig6_original_featured_provisional']=='yes'; isc=(a['featured']=='yes' and not isf)
        ax.add_patch(Rectangle((s,r-0.28),e-s,0.56,facecolor=col,edgecolor=RED if isf else "none",linewidth=1.1 if isf else 0))
        ax.plot([e,e],[r-0.4,r+0.4],color="#222",lw=0.6)
        if isf or isc:
            yl=n+0.2+0.8*(cnt%2); cnt+=1
            ax.plot([(s+e)/2]*2,[r+0.3,yl],color=RED if isf else "#555",lw=0.5)
            ax.text((s+e)/2,yl,a['figure6_name'],ha="center",va="bottom",fontsize=7,color=RED if isf else "#444",fontweight="bold" if isf else "normal")
    TX(fig,LM+left-0.08,y+h/2,FLAB[f],fontsize=7,ha="right",va="center")
    y+=h+0.04
xs=[int(d['window_center']) for d in dens]
for key,lab in (("families_per_kb","Families per kb"),("start_sites_per_kb","Start sites per kb")):
    ax=AX(fig,LM+left,y+0.05,wA,0.42)
    v=[float(d[key]) for d in dens]
    ax.fill_between(xs,v,step="mid",color="#7a7a7a",lw=0); ax.set_xlim(0,VLEN); ax.set_ylim(0,max(v)*1.1)
    ax.set_ylabel(lab,fontsize=7); ax.spines[["top","right"]].set_visible(False); ax.tick_params(labelsize=7)
    if key=="families_per_kb": ax.set_xticklabels([])
    else: ax.set_xlabel("Genome position, refdata_GIb_v1.1 (nt)",fontsize=7)
    y+=0.42+0.22
yB=y+0.3
# ---- B (left) and C3 (right)
TX(fig,LM,yB,"B",fontsize=10,fontweight="bold",va="top")
regs=[b['region'] for b in Bd if b['length_nt']]
bw=0.95; bh=1.1; by=yB+0.35
def barax(x0,key,title,fmt,ylab):
    ax=AX(fig,x0,by,bw,bh); v=[float(b[key]) for b in Bd if b['length_nt']]
    ax.bar(range(len(regs)),v,color=[RCOL[r] for r in regs],width=0.75)
    ax.set_xticks(range(len(regs))); ax.set_xticklabels(regs,rotation=90,fontsize=7)
    ax.set_title(title,fontsize=7,loc="left"); ax.spines[["top","right"]].set_visible(False); ax.tick_params(labelsize=7)
    ax.set_ylim(0,max(v)*1.18)
    for i,t in enumerate(v): ax.text(i,t,fmt(t),ha="center",va="bottom",fontsize=7)
    return ax
barax(LM+0.40,"families_per_kb","B1 Families per kb",lambda t:f"{t:.1f}",None)
barax(LM+0.40+bw+0.33,"start_sites_per_kb","B2 Start sites per kb",lambda t:f"{t:.0f}",None)
a3=AX(fig,LM+0.40+2*(bw+0.33),by,bw+0.05,bh)
order=regs+['UTR/gap']; rng=np.random.default_rng(3)
for i,r in enumerate(order):
    v=[int(b['n_start_variants_gt40']) for b in Bs if b['region']==r]
    if not v: continue
    a3.scatter(i+rng.uniform(-0.18,0.18,len(v)),v,s=7,color=RCOL.get(r,"#888"),edgecolor="black" if r=='N' else "none",linewidths=0.4,zorder=3)
    a3.plot([i-0.3,i+0.3],[np.median(v)]*2,color="black",lw=0.9)
a3.set_xticks(range(len(order))); a3.set_xticklabels(order,rotation=90,fontsize=7)
a3.set_title("B3 Starts per family",fontsize=7,loc="left"); a3.spines[["top","right"]].set_visible(False); a3.tick_params(labelsize=7)
TX(fig,LM,by+bh+0.62,textwrap.fill(f"N vs other CDS: families/kb {Bst['N_fam_per_kb']:.2f} vs {Bst['rest_fam_per_kb']:.2f} (rate ratio {Bst['fam_rate_ratio_N_over_rest']:.2f}, exact binomial P={Bst['fam_rate_binom_p']:.2f}); starts per family median {Bst['starts_per_family_median_N']:.0f} vs {Bst['starts_per_family_median_rest']:.0f} (Mann–Whitney P={Bst['mannwhitney_p']:.2f}): numerically higher, not significant. n={len(nc)} families; bar = median.",100),fontsize=7,va="top")
# C3
cx0=LM+4.7
TX(fig,cx0-0.5,yB,"C3",fontsize=10,fontweight="bold",va="top")
ax3=AX(fig,cx0,by,2.35,bh)
fidn={r['family_id']:r for r in nullr}
for r in nullr:
    f=r['family_id']; isf=f in C12
    ax3.scatter(float(r['expected_fraction_synonymous_null']),float(r['observed_fraction']),s=16 if isf else 6,color="white" if isf else "#8a8a8a",edgecolor="black" if isf else "none",linewidths=1.0,zorder=3 if isf else 2)
    if isf and float(r['observed_fraction'])<0.8: ax3.annotate(NAME[f],(float(r['expected_fraction_synonymous_null']),float(r['observed_fraction'])),xytext=(4,3),textcoords="offset points",fontsize=7,color=RED)
ax3.plot([0,1],[0,1],color="#333",lw=0.6,ls="--"); ax3.set_xlim(0,1.03); ax3.set_ylim(0,1.03)
ax3.set_xlabel("Expected under synonymous-substitution null",fontsize=7); ax3.set_ylabel("Observed retention",fontsize=7)
ax3.set_title(f"Observed vs null (n={len(nullr)}; rings = C1/C2 families)",fontsize=7,loc="left")
ax3.spines[["top","right"]].set_visible(False); ax3.tick_params(labelsize=7)
# ---- D
yD=by+bh+1.05
TX(fig,LM,yD,"D",fontsize=10,fontweight="bold",va="top")
fm={a['family_id']:a for a in A}
order=[a['family_id'] for a in A if a['canonical_frame_family']=='no']
rowsD=['GI-a','GI-b','GII-a','GII-b','GII-c','unresolved','ALL','N1_historical_G1b_label']
tot={}
for g in rowsD:
    tot[g]=len(typ) if g=='ALL' else (sum(1 for t in typ if t['N1_genotype_historical']=='GIb') if g.startswith('N1_') else sum(1 for t in typ if t['S_type']==g))
cm=LinearSegmentedColormap.from_list("pink_white_skyblue",["#4fb3e8","#a6d8f0","#ffffff","#f8bbd0","#f06292"])
hx=LM+1.55; hw=WID-1.55; hy=yD+0.62; rh=0.14
axs=AX(fig,hx,yD+0.38,hw,0.1); axs.set_xlim(0,len(order)); axs.set_ylim(0,1); axs.axis("off")
axh=AX(fig,hx,hy,hw,rh*len(rowsD))
nf=0
for j,fid in enumerate(order):
    p=fm[fid]['parent']; axs.add_patch(Rectangle((j,0),1,1,facecolor=RCOL.get(p,"#999"),edgecolor="none"))
    if fm[fid]['fig6_original_featured_provisional']=='yes':
        axs.add_patch(Polygon([[j+0.5,0],[j-0.5,-1.6],[j+1.5,-1.6]],closed=True,facecolor=RED,edgecolor="none"))
        lv=2.2+2.6*(nf%2); nf+=1
        axs.plot([j+0.5,j+0.5],[1,lv],color=RED,lw=0.4)
        axs.text(j+0.5,lv,fm[fid]['figure6_name'],fontsize=7,ha="center",va="bottom",color=RED)
for i,g in enumerate(rowsD[::-1]):
    for j,fid in enumerate(order):
        r=L[(fid,g)]; n=int(r['n_assessable'])
        if n<10 or r['fraction_retained']=='': axh.add_patch(Rectangle((j,i),1,1,facecolor="white",edgecolor="#bbb",linewidth=0.3,hatch="////"))
        else: axh.add_patch(Rectangle((j,i),1,1,facecolor=cm(float(r['fraction_retained'])),edgecolor="none"))
axh.set_xlim(0,len(order)); axh.set_ylim(0,len(rowsD)); axh.set_xticks([])
labs={'ALL':'all PEDV','N1_historical_G1b_label':"historical 'G1b' label"}
axh.set_yticks([i+0.5 for i in range(len(rowsD))]); axh.set_yticklabels([f"{labs.get(g,g)} (n={tot[g]})" for g in rowsD[::-1]],fontsize=7)
for lab in axh.get_yticklabels():
    if lab.get_text().startswith('GI-b'): lab.set_fontweight('bold')
gi=rowsD[::-1].index('GI-b'); axh.add_patch(Rectangle((0,gi),len(order),1,facecolor="none",edgecolor="black",linewidth=1.2))
for s in axh.spines.values(): s.set_visible(False)
axh.tick_params(length=0)
cb=AX(fig,hx,hy+rh*len(rowsD)+0.12,1.6,0.1)
cb.imshow(np.linspace(0,1,100)[None,:],aspect="auto",cmap=cm,extent=[0,1,0,1]); cb.set_yticks([]); cb.set_xticks([0,0.5,1]); cb.tick_params(labelsize=7)
TX(fig,hx+1.75,hy+rh*len(rowsD)+0.2,"Fraction of assessable genomes retaining the family (hatched: <10 assessable)",fontsize=7,va="center")
# D2
d2y=hy+rh*len(rowsD)+0.58
ax4=AX(fig,hx,d2y,hw,0.75)
grp=['GI-a','GI-b','GII-a','GII-b','GII-c','unresolved']; rng=np.random.default_rng(5)
for i,g in enumerate(grp):
    v=[float(r['fraction']) for r in pg if r['S_type']==g]
    if v:
        ax4.scatter(i+rng.uniform(-0.3,0.3,len(v)),v,s=3,color="black" if g=='GI-b' else "#9a9a9a",edgecolor="none",alpha=0.7)
        ax4.plot([i-0.35,i+0.35],[np.median(v)]*2,color="black",lw=1.0)
ax4.set_xticks(range(len(grp))); ax4.set_xticklabels([f"{g} (n={tot[g]})" for g in grp],fontsize=7)
ax4.set_ylabel("Fraction of assessable"+chr(10)+"families retained",fontsize=7); ax4.set_ylim(0.5,1.0)
ax4.spines[["top","right"]].set_visible(False); ax4.tick_params(labelsize=7)
ax4.set_title("D2  Per-genome completeness (bar = median; not an estimate of translational capacity)",fontsize=7,loc="left")
TX(fig,LM,FH-0.30,textwrap.fill("Lineage = nearest S-gene panel member (Guo 2019 anchors; GI-b n=30 on 5 anchors, GI-a n=23 on 3); S-gene lineage need not reflect the history of other loci (recombination). Retention = sequence retention of eligible start/stop architecture, not evidence of translation, expression or protein-level selection.",170),fontsize=7,va="bottom")
PP.savefig(fig)
for ext,kw in (("png",dict(dpi=300)),("jpg",dict(dpi=300,pil_kwargs={"quality":95})),("svg",{}),("pdf",{})):
    fig.savefig(O+f"Supplementary_Figure_L_REVISED_v3_p1.{ext}",**kw)
fig.savefig(O+"work/v3_p1.png",dpi=130); plt.close(fig)
# ============================================================ PAGE 2
fig=plt.figure(figsize=(FW,FH))
fams,c1,c2=f3.load(); n=len(fams); LG=f3.LG
TX(fig,LM,0.12,"Supplementary Figure L (continued)  Sequence retention of the families currently highlighted in Figure 6 (provisional list)",fontsize=9,fontweight="bold",va="top")
TX(fig,LM,0.42,"C",fontsize=10,fontweight="bold",va="top")
TX(fig,LM+0.5,0.45,"C1  Start-codon state of reference start variants (>40 aa) by S-gene lineage",fontsize=7,va="top")
x0=LM+1.0; wC=WID-1.0; y0=1.15; rowh=0.25
ax=AX(fig,x0,y0,wC,rowh*n)
colw=1.0; pad=0.17; ncol=len(LG)
for i,f in enumerate(fams):
    yy=n-1-i
    for j,g in enumerate(LG):
        r=c1[(f,g)]; left_=j*(colw+pad)
        for c,col in zip(f3.CLS,f3.CCOL):
            v=float(r[c+'_fraction']) if r[c+'_fraction']!='' else 0.0
            if v>0: ax.add_patch(Rectangle((left_,yy-0.34),v*colw,0.68,facecolor=col,edgecolor="#555" if col=="#ffffff" else "none",linewidth=0.4))
            left_+=v*colw
        ax.add_patch(Rectangle((j*(colw+pad),yy-0.34),colw,0.68,facecolor="none",edgecolor="#777",linewidth=0.4))
ax.set_xlim(-0.02,ncol*(colw+pad)-pad+0.02); ax.set_ylim(-0.6,n-0.4)
ax.set_yticks(range(n)); ax.set_yticklabels(fams[::-1],fontsize=7)
nn={g:c2[(fams[0],g)]['n_genomes'] for g in LG}
ax.set_xticks([j*(colw+pad)+colw/2 for j in range(ncol)]); ax.set_xticklabels([("All PEDV" if g=='ALL' else g)+f" n={nn[g]}" for g in LG],fontsize=7)
ax.xaxis.tick_top(); ax.tick_params(length=0)
for s in ('left','right','bottom','top'): ax.spines[s].set_visible(False)
ax.legend(handles=[Rectangle((0,0),1,1,facecolor=c,edgecolor="#555",linewidth=0.4,label=l) for c,l in zip(f3.CCOL,f3.CLAB)],loc="upper left",bbox_to_anchor=(0.02,-0.02),ncol=4,frameon=False,fontsize=7,handlelength=1.1,columnspacing=1.0)
TX(fig,LM,y0+rowh*n+0.75,"Fractions are over reference start variant × genome pairs; start variants of one family are not independent; an eligible start codon is not an initiation event.",fontsize=7,va="top")
# C2
yC2=y0+rowh*n+1.15
TX(fig,LM+0.5,yC2,"C2  Family-level retention (stop anchor + ≥1 eligible reference start; Wilson 95% CI; not-assessable excluded; n<10 omitted)",fontsize=7,va="top")
rh2=0.62
ax2=AX(fig,x0,yC2+0.25,wC,rh2*n)
for i,f in enumerate(fams):
    yy=n-1-i
    for j,g in enumerate(LG):
        r=c2[(f,g)]
        if r['retention_fraction']=='' or int(r['n_assessable'])<10: continue
        y2=yy+(j-2.5)*0.13; v=float(r['retention_fraction'])
        ax2.plot([float(r['wilson95_lo']),float(r['wilson95_hi'])],[y2,y2],color=f3.LCOL[g],lw=1.0)
        ax2.scatter(v,y2,s=26 if g in('GI-b','ALL') else 16,color=f3.LCOL[g],edgecolor="black" if g=='GI-b' else "none",linewidths=0.6,zorder=3,marker='D' if g=='ALL' else 'o')
ax2.set_ylim(-0.6,n-0.4); ax2.set_xlim(-0.03,1.03); ax2.set_yticks(range(n)); ax2.set_yticklabels(fams[::-1],fontsize=7)
for i in range(n): ax2.axhline(i-0.5,color="#e4e4e4",lw=0.5)
ax2.set_xlabel("Family ORF architecture retained / assessable genomes",fontsize=7)
ax2.spines[["top","right"]].set_visible(False); ax2.tick_params(labelsize=7)
ax2.legend(handles=[Line2D([],[],marker='D' if g=='ALL' else 'o',ls='none',color=f3.LCOL[g],markeredgecolor='black' if g=='GI-b' else 'none',ms=5,label=('All PEDV' if g=='ALL' else g)) for g in LG],
           loc="upper left",bbox_to_anchor=(0.0,-0.07),ncol=6,frameon=False,fontsize=7,handletextpad=0.3,columnspacing=1.2)
TX(fig,LM,FH-0.30,textwrap.fill("Families: provisional Figure 6 list (checkpoint-3 evidence architecture; Figure 6 not frozen), chosen by Figure 6 evidence class and not by conservation; order = genomic position of the most upstream start. No Ribo-seq quantity was used in any calculation; retention is not evidence of translation or protein function.",150),fontsize=7,va="bottom")
PP.savefig(fig); PP.close()
for ext,kw in (("png",dict(dpi=300)),("jpg",dict(dpi=300,pil_kwargs={"quality":95})),("svg",{}),("pdf",{})):
    fig.savefig(O+f"Supplementary_Figure_L_REVISED_v3_p2.{ext}",**kw)
fig.savefig(O+"work/v3_p2.png",dpi=130)
print("ok")
