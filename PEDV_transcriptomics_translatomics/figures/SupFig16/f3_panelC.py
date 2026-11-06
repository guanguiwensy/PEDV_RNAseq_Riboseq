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
import csv,json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
O=PEDV_K+r"Supplementary_Figure_L/"; W=O+"work/"
CLS=['exact_codon_retention','AUG_retention','near_cognate_retention','start_class_switching','loss_of_eligible_start','not_assessable']
CLAB=['exact non-AUG codon retained','AUG retained','near-cognate retained (different codon)','start-class switch','loss of eligible start','not assessable']
CCOL=['#4c4c4c','#5e3c99','#a6761d','#8fa8c8','#d9d9d9','#ffffff']
LG=['GI-a','GI-b','GII-a','GII-b','GII-c','ALL']
LCOL={'GI-a':'#7f8fa8','GI-b':'#000000','GII-a':'#a0a0a0','GII-b':'#b59c80','GII-c':'#7fa88a','ALL':'#c1121f'}
def load():
    sel=list(csv.DictReader(open(O+"SupL_C1C2_family_selection.csv",encoding='utf-8-sig')))
    c1={(r['family'],r['lineage']):r for r in csv.DictReader(open(O+"SourceData_PanelC1_start_codon_by_lineage.csv",encoding='utf-8-sig'))}
    c2={(r['family'],r['lineage']):r for r in csv.DictReader(open(O+"SourceData_PanelC2_family_retention_by_lineage.csv",encoding='utf-8-sig'))}
    return [r['family'] for r in sel],c1,c2
def draw(fig,x0,y0,w,h,lab_c='C'):
    fams,c1,c2=load()
    n=len(fams)
    wc1=0.62*w; wc2=0.26*w; gap=0.12*w
    ax1=fig.add_axes([x0,y0,wc1,h]); ax2=fig.add_axes([x0+wc1+gap,y0,wc2,h],sharey=None)
    # C1: rows=families, columns=lineages (mini stacked bars)
    colw=1.0; pad=0.25
    for i,f in enumerate(fams):
        y=n-1-i
        for j,g in enumerate(LG):
            r=c1[(f,g)]; left=j*(colw+pad)
            for c,col in zip(CLS,CCOL):
                v=float(r[c+'_fraction']) if r[c+'_fraction']!='' else 0.0
                if v>0: ax1.add_patch(Rectangle((left,y-0.32),v*colw,0.64,facecolor=col,edgecolor="#555" if col=="#ffffff" else "none",linewidth=0.4))
                left+=v*colw
            ax1.add_patch(Rectangle((j*(colw+pad),y-0.32),colw,0.64,facecolor="none",edgecolor="#777",linewidth=0.4))
    ax1.set_xlim(-0.1,len(LG)*(colw+pad)); ax1.set_ylim(-0.6,n-0.4)
    ax1.set_yticks(range(n)); ax1.set_yticklabels(fams[::-1],fontsize=7)
    ax1.set_xticks([j*(colw+pad)+colw/2 for j in range(len(LG))])
    nn={g:c2[(fams[0],g)]['n_genomes'] for g in LG}
    ax1.set_xticklabels([("All PEDV" if g=='ALL' else g)+chr(10)+"n="+nn[g] for g in LG],fontsize=7)
    ax1.xaxis.tick_top(); ax1.tick_params(length=0)
    for s in ('left','right','bottom','top'): ax1.spines[s].set_visible(False)
    ax1.set_title("C1  Start-codon state of reference start variants (>40 aa) by lineage",fontsize=7,loc="left",pad=22)
    # C2
    for i,f in enumerate(fams):
        y=n-1-i
        for j,g in enumerate(LG):
            r=c2[(f,g)]
            if r['retention_fraction']=='' or int(r['n_assessable'])<10: continue
            yy=y+(j-2.5)*0.11; v=float(r['retention_fraction'])
            ax2.plot([float(r['wilson95_lo']),float(r['wilson95_hi'])],[yy,yy],color=LCOL[g],lw=0.8)
            ax2.scatter(v,yy,s=14 if g in('GI-b','ALL') else 8,color=LCOL[g],edgecolor="black" if g=='GI-b' else "none",linewidths=0.5,zorder=3,marker='D' if g=='ALL' else 'o')
    ax2.set_ylim(-0.6,n-0.4); ax2.set_xlim(-0.03,1.03); ax2.set_yticks(range(n)); ax2.set_yticklabels([])
    ax2.set_xlabel("Family ORF architecture retained /"+chr(10)+"assessable genomes",fontsize=7)
    ax2.set_title("C2  Family-level retention",fontsize=7,loc="left",pad=22)
    ax2.spines[["top","right"]].set_visible(False); ax2.tick_params(labelsize=7)
    for i in range(n): ax2.axhline(i+0.5,color="#eee",lw=0.4)
    ax2.legend(handles=[Line2D([],[],marker='D' if g=='ALL' else 'o',ls='none',color=LCOL[g],markeredgecolor='black' if g=='GI-b' else 'none',ms=4,label=('All PEDV' if g=='ALL' else g)) for g in LG],
               loc="upper left",bbox_to_anchor=(-0.02,-0.20),ncol=3,frameon=False,fontsize=7,handletextpad=0.2,columnspacing=0.8)
    ax1.legend(handles=[Rectangle((0,0),1,1,facecolor=c,edgecolor="#555",linewidth=0.4,label=l) for c,l in zip(CCOL,CLAB)],loc="upper left",bbox_to_anchor=(0,-0.04),ncol=2,frameon=False,fontsize=7,handlelength=1.1)
    return ax1,ax2
if __name__=="__main__":
    plt.rcParams.update({"font.family":["Arial","DejaVu Sans"],"font.size":7,"axes.linewidth":0.6,"svg.fonttype":"none","pdf.fonttype":42,"figure.facecolor":"white","axes.facecolor":"white"})
    fig=plt.figure(figsize=(7.4,6.0))
    draw(fig,0.13,0.27,0.84,0.55)
    fig.text(0.01,0.97,"C",fontsize=10,fontweight="bold",va="top")
    fig.text(0.01,0.035,"Configured Figure 6 families."+chr(10)+"Start-codon state is not initiation; retention is not translation evidence. Wilson 95% intervals; n<10 omitted.",fontsize=7,va="bottom")
    for ext in ("pdf","svg","png"): fig.savefig(O+"Supplementary_Figure_16_PanelC1C2."+ext,dpi=300)
    fig.savefig(O+"work/c12.png",dpi=140)
