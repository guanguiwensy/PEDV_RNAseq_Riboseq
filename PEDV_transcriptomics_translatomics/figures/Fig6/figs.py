"""Stage 3 figures R6-A (small quantitative companion to the existing N1 three-frame heat maps), R6-B, R6-C, R6-D. Placeholder letters only in file names/titles.
Vector rules: data trimmed to axis limits beforehand, clip_on=False on all artists (no clipping masks relied upon), svg.fonttype='none', pdf.fonttype=42, font >= 7 pt."""
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
import os
import matplotlib; matplotlib.use('Agg'); import matplotlib.ticker
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
from r6lib import *
plt.rcParams.update({'svg.fonttype':'none','pdf.fonttype':42,'font.family':'DejaVu Sans','font.size':7.5,'axes.titlesize':8,'axes.labelsize':8,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'legend.fontsize':7.5,
  'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':0.7,'xtick.major.width':0.7,'ytick.major.width':0.7})
RCOL={"ORF1a":"#0072B2","ORF1b":"#56B4E9","S":"#009E73","ORF3":"#CC79A7","E":"#D55E00","M":"#E69F00","N":"#117733"}   # Supplementary Figure L palette
GR="#8a8a8a"; FR={'0':'#bdbdbd','+1':'#5e3c99','+2':'#a6761d'}
SD=OUT+'source_data/'; FD=OUT+'figures/'
os.makedirs(FD,exist_ok=True)
def finish(fig,name):
    for ax in fig.axes:
        for a in ax.get_children():
            try: a.set_clip_on(False)
            except Exception: pass
    for ext,kw in (('pdf',{}),('svg',{}),('png',{'dpi':300})): fig.savefig(FD+name+'.'+ext,**kw)
    plt.close(fig)
def lab(ax,t): ax.text(-0.16,1.06,t,transform=ax.transAxes,fontsize=10,fontweight='bold',va='bottom')
# ---------------- R6-D
T=pd.read_csv(SD+'G5_per_family_disruption_vs_null.csv'); t=T[(T.scope=='cluster_reps_579')&(T.criterion=='primary')]
fig,axs=plt.subplots(1,2,figsize=(7.6,3.3),gridspec_kw={'width_ratios':[1.6,1]})
ax=axs[0]; regs=['ORF1a','ORF1b','S','ORF3','M','N']; rng=np.random.default_rng(1)
for i,r in enumerate(regs):
    for off,role,col,ms in ((-0.2,'control_8to40',GR,3),(0.2,'target_gt40',RCOL[r],7)):
        v=t[(t.region==r)&(t.role==role)].excess_rate.values
        if len(v)==0: continue
        ax.scatter(i+off+rng.uniform(-0.14,0.14,len(v)),v,s=ms,color=col,lw=0,alpha=0.8)
        ax.plot([i+off-0.17,i+off+0.17],[np.median(v)]*2,color='k',lw=1.2)
ax.axhline(0,color='k',lw=0.6,ls='--'); ax.set_xticks(range(6)); ax.set_xticklabels(regs); ax.set_ylim(-0.3,1.0)
ax.set_ylabel('observed - expected disruption rate'); ax.set_title('Per family (579 cluster representatives)')
ax.legend(handles=[Line2D([],[],marker='o',color=GR,ls='',ms=3,label='control: 8-40 aa families (n=293)'),Line2D([],[],marker='o',color='#444',ls='',ms=5,label='target: >40 aa families (n=71), region colour')],frameon=False,loc='upper center',bbox_to_anchor=(0.5,-0.1),ncol=1,fontsize=7)
ax=axs[1]; rowsD=[]; labs=[]; y=0
for sc,scl in (('cluster_reps_579','579 clusters'),('full_1083','1,083 genomes')):
    for role,rl in (('target_gt40','targets'),('control_8to40','controls')):
        d=T[(T.scope==sc)&(T.criterion=='primary')&(T.role==role)]; n=len(d)
        sig=int(((d.q_lower<.05)&(d.obs_rate<d.null_rate)).sum()); low=int((d.obs_rate<d.null_rate).sum())-sig; hi=int((d.obs_rate>=d.null_rate).sum())
        ax.barh(y,sig/n,color='#2a7f62',lw=0); ax.barh(y,low/n,left=sig/n,color='#a8d5c2',lw=0); ax.barh(y,hi/n,left=(sig+low)/n,color='#d9d9d9',lw=0)
        ax.text(0.02,y,f'{sig+low}/{n} lower; {sig} BH',va='center',fontsize=7,color='white' if sig/n>0.3 else 'k')
        labs.append(f'{rl}, {scl}'); rowsD.append(dict(scope=sc,role=role,n=n,lower_than_null=sig+low,lower_BH_significant=sig,not_lower=hi)); y+=1
ax.set_yticks(range(4)); ax.set_yticklabels(labs,fontsize=7); ax.invert_yaxis(); ax.set_xlim(0,1); ax.set_xlabel('fraction of families'); ax.set_title('Fraction below null')
ax.legend(handles=[plt.Rectangle((0,0),1,1,color='#2a7f62',label='lower, BH q<0.05'),plt.Rectangle((0,0),1,1,color='#a8d5c2',label='lower, n.s.'),plt.Rectangle((0,0),1,1,color='#d9d9d9',label='not lower')],frameon=False,loc='upper center',bbox_to_anchor=(0.2,-0.22),ncol=3,fontsize=7,columnspacing=0.8,handlelength=1)
for a,tt in zip(axs,'AB'): lab(a,tt)
axs[1].texts[-1].set_x(-0.9); axs[0].texts[-1].set_x(-0.25)
fig.suptitle('R6-D  Candidate-ORF disruption vs random synonymous mutation (10,000 draws; BH)',fontsize=8,y=0.995)
fig.tight_layout(rect=(0,0,1,0.955)); finish(fig,'R6-D_candidate_vs_control_disruption')
t.to_csv(SD+'R6-D_panelA_per_family.csv',index=False,encoding='utf-8-sig'); pd.DataFrame(rowsD).to_csv(SD+'R6-D_panelB_counts.csv',index=False)
print('ok')
