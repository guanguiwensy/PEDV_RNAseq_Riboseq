"""R6-A, R6-B, R6-C (with chi-square p-values; C without ORF3/E) and R6-E (translation efficiency vs long-ORF capacity, n=5). Run stats_abc.py first.
Vector rules as in figs.py: clip_on=False, fonttype none/42, fonts >= 7 pt."""
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
import os, ast
import matplotlib; matplotlib.use('Agg'); import matplotlib.ticker
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
from r6lib import *
plt.rcParams.update({'svg.fonttype':'none','pdf.fonttype':42,'font.family':'DejaVu Sans','font.size':7.5,'axes.titlesize':8,'axes.labelsize':8,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'legend.fontsize':7.5,
  'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':0.7,'xtick.major.width':0.7,'ytick.major.width':0.7})
RCOL={"ORF1a":"#0072B2","ORF1b":"#56B4E9","S":"#009E73","ORF3":"#CC79A7","E":"#D55E00","M":"#E69F00","N":"#117733"}
FR={'0':'#bdbdbd','+1':'#5e3c99','+2':'#a6761d'}
SD=OUT+'source_data/'; FD=OUT+'figures/'; os.makedirs(FD,exist_ok=True)
def finish(fig,name):
    for ax in fig.axes:
        for a in ax.get_children():
            try: a.set_clip_on(False)
            except Exception: pass
    for ext,kw in (('pdf',{}),('svg',{}),('png',{'dpi':300})): fig.savefig(FD+name+'.'+ext,**kw)
    plt.close(fig)
def lab(ax,t,x=-0.16): ax.text(x,1.06,t,transform=ax.transAxes,fontsize=10,fontweight='bold',va='bottom')
def fp(p):
    if p<1e-300: return 'p<1e-300'
    if p>=0.001: return f'p={p:.3f}'.rstrip('0').rstrip('.') if p<0.1 else f'p={p:.2f}'
    e=int(np.floor(np.log10(p))); m=p/10**e; return f'p={m:.1f}e{e}'
def bracket(ax,x1,x2,y,text,h):
    ax.plot([x1,x1,x2,x2],[y,y+h,y+h,y],color='k',lw=0.7); ax.text((x1+x2)/2,y+h*1.15,text,ha='center',va='bottom',fontsize=7)
# ---------------- R6-A
g=pd.read_csv(SD+'G1_codon_density_by_frame_region.csv',dtype={'rel_frame':str}); g=g[g.region=='ALL']
A=pd.read_csv(SD+'R6-A_chi_square.csv')
def pA(cls,cmp_): return float(A[(A.codon_class==cls)&(A.comparison==cmp_)].p_chi2.iloc[0])
fig,axs=plt.subplots(1,3,figsize=(7.6,3.1))
cfg=(('stop (UAA/UAG/UGA)','stop','Stop codons',58,[(1,2,'frame +1 vs frame +2')]),('AUG','AUG','AUG',38,[(1,2,'frame +1 vs frame +2'),(0,1,'frame 0 vs frame +1'),(0,2,'frame 0 vs frame +2')]),
     ('AUG+7 near-cognate','AUG+7 near-cognate','AUG + 7 near-cognate',175,[(1,2,'frame +1 vs frame +2'),(0,1,'frame 0 vs frame +1'),(0,2,'frame 0 vs frame +2')]))
for ax,(c,cc,t,yl,brs) in zip(axs,cfg):
    d=g[g.codon_class==c].set_index('rel_frame')
    for i,fr in enumerate(('0','+1','+2')):
        r=d.loc[fr]; ax.bar(i,r.observed_per_kb,0.62,color=FR[fr],lw=0)
        ax.errorbar(i,r.observed_per_kb,yerr=[[r.observed_per_kb-r.obs_poisson95_lo_per_kb],[r.obs_poisson95_hi_per_kb-r.observed_per_kb]],color='k',lw=0.7,capsize=2)
    top=max(d.obs_poisson95_hi_per_kb.max(),1); step=yl*0.06; y0=top*1.1
    for j,(x1,x2,cmp_) in enumerate(sorted(brs,key=lambda b:b[1]-b[0])):
        bracket(ax,x1,x2,y0+j*step*1.7,'χ² '+fp(pA(cc,cmp_)),step*0.5)
    ax.set_xticks(range(3)); ax.set_xticklabels(['0','+1','+2']); ax.set_ylim(0,yl); ax.set_title(t); ax.set_xlabel('reading frame (relative to ORF start)')
axs[0].set_ylabel('codons per kb of CDS')
for a_,t_ in zip(axs,'ABC'): lab(a_,t_)
fig.text(0.5,0.0,'Bars: observed per kb (Poisson 95% CI). χ²: 2x2 test, hit vs non-hit triplets (8,975-8,980 triplets per frame),\nno continuity correction; overlapping frames are not independent, so p is descriptive.',ha='center',fontsize=7)
fig.suptitle('R6-A  Start/stop codon density by frame, chi-square tests',fontsize=8,y=0.995)
fig.tight_layout(rect=(0,0.11,1,0.95)); finish(fig,'R6-A_codon_density_by_frame')
# ---------------- R6-B
fam=pd.read_csv(SD+'G2_family_lengths_plus1_plus2.csv',dtype={'rel_frame':str}); hist=pd.read_csv(SD+'G2_length_hist_obs_vs_null.csv'); Bc=pd.read_csv(SD+'R6-B_chi_square.csv')
iso=isoforms(); iso=iso[iso.family_id.isin(fam.family_id)].merge(fam[['family_id','rel_frame']],on='family_id')
edges=np.arange(8,164,5)
def binned(v): return np.histogram(np.minimum(v,edges[-1]-1),bins=edges)[0]
fig,axs=plt.subplots(2,2,figsize=(7.6,5.2),sharex=True)
for j,rf in enumerate(('+1','+2')):
    v=fam[fam.rel_frame==rf].max_aa.values; ob=binned(v)
    hn=hist['null_mean_p%d'%(j+1)].values; aa=hist.aa.values
    nb=np.array([hn[(aa>=lo)&(aa<hi)].sum() for lo,hi in zip(edges[:-1],edges[1:])])
    ax=axs[0,j]; ax.bar(edges[:-1],ob,width=5,align='edge',color=FR[rf],lw=0,alpha=0.85)
    ax.step(np.r_[edges[:-1],edges[-1]],np.r_[nb,nb[-1]],where='post',color='#d62728',lw=1.2)
    ymax=float(max(ob.max(),nb.max()))*1.15; ax.axvline(40.5,color='k',lw=0.7,ls='--'); ax.set_ylim(0,ymax)
    ax.set_title(f'{rf} frame: stop-anchored families (n={len(v)})'); ax.text(43,ymax*0.96,f'>40 aa: {(v>40).sum()}  (max {v.max()} aa)',fontsize=7.5,va='top')
    vi=iso[iso.rel_frame==rf].peptide_length_aa.values; obi=binned(vi)
    ax2=axs[1,j]; ax2.bar(edges[:-1],obi,width=5,align='edge',color=FR[rf],lw=0,alpha=0.85); ax2.axvline(40.5,color='k',lw=0.7,ls='--'); ax2.set_ylim(0,520)
    ax2.set_title(f'{rf} frame: start isoforms (n={len(vi)})'); ax2.text(43,500,f'>40 aa: {(vi>40).sum()}  (max {vi.max()} aa)',fontsize=7.5,va='top')
for ax in axs[:,0]: ax.set_ylabel('count per 5-aa bin')
for ax in axs[1]: ax.set_xlabel('peptide length (aa, stop excluded; last bin = 158+ clipped)'); ax.set_xlim(8,163)
axs[0,1].legend(handles=[Line2D([],[],color='#d62728',lw=1.2,label='composition null (mean)'),Line2D([],[],color='k',lw=0.7,ls='--',label='40 aa')],frameon=False,loc='upper right',bbox_to_anchor=(1,0.78))
for a_,t_ in zip(axs.ravel(),'ABCD'): lab(a_,t_)
def row(unit,kind): r=Bc[(Bc.unit==unit)&(Bc.test.str.contains(kind))].iloc[0]; return r
rf2,rf3=row('family','2x2'),row('family','2x3'); ri2,ri3=row('start isoform','2x2'),row('start isoform','2x3')
fig.text(0.5,0.005,f'Chi-square, +1 vs +2 frame — families: 2x3 (8-20/21-40/>40 aa) χ²={rf3.chi2:.1f}, df=2, {fp(rf3.p)}; ≤40 vs >40 aa: χ²={rf2.chi2:.1f}, {fp(rf2.p)} (Fisher {fp(rf2.fisher_p)}).\nStart isoforms: 2x3 χ²={ri3.chi2:.1f}, {fp(ri3.p)}; ≤40 vs >40 aa: χ²={ri2.chi2:.1f}, {fp(ri2.p)}.',ha='center',fontsize=7.5)
fig.suptitle('R6-B  ORF length distribution, +1 vs +2 (top: families + null; bottom: isoforms), chi-square tests',fontsize=8,y=0.995)
fig.tight_layout(rect=(0,0.09,1,0.95)); finish(fig,'R6-B_length_histogram_plus1_plus2')
# ---------------- R6-C (ORF3 and E excluded)
Cc=pd.read_csv(SD+'R6-C_chi_square.csv'); c71=Cc.iloc[0]; c68=Cc.iloc[1]
cnt71=ast.literal_eval(c71.counts); cnt68=ast.literal_eval(c68.counts); O5=['ORF1a','ORF1b','S','M','N']
def garwood(k): return (0 if k==0 else stats.chi2.ppf(.025,2*k)/2, stats.chi2.ppf(.975,2*k+2)/2)
fig,ax=plt.subplots(figsize=(4.8,3.4))
for i,r in enumerate(O5):
    L=CDS[r][1]-CDS[r][0]+1; k=cnt71[r]; lo,hi=garwood(k)
    ax.bar(i,1000*k/L,0.65,color=RCOL[r],lw=0); ax.plot([i,i],[1000*lo/L,1000*hi/L],color='k',lw=0.7); ax.plot([i-0.15,i+0.15],[1000*hi/L]*2,color='k',lw=0.7)
    ax.text(i,1000*hi/L+0.25,f'n={k}',ha='center',fontsize=7.5)
    ax.plot(i,1000*cnt68[r]/L,marker='D',mfc='none',mec='k',ms=4,ls='',zorder=6)
ax.set_xticks(range(5)); ax.set_xticklabels([f'{r}\n{CDS[r][1]-CDS[r][0]+1}' for r in O5]); ax.set_ylim(0,16); ax.set_ylabel('>40-aa families per kb')
ax.set_title('Long alt-frame ORFs per kb (ORF3, E excluded)',pad=12)
ax.legend(handles=[Line2D([],[],marker='D',mfc='none',mec='k',ls='',ms=4,label='+1 frame only (68 families)'),Line2D([],[],color='k',lw=0.7,label='exact Poisson 95% CI (bars: +1 and +2, 70 families)')],frameon=False,loc='upper left',fontsize=7)
fig.text(0.5,0.0,f'χ² goodness of fit, rate proportional to length: χ²={c71.chi2_GOF_vs_length_proportional:.2f}, df=4, {fp(c71.p_chi2)} (Monte Carlo {fp(c71.p_monte_carlo_100000)}).\nM+N vs the other three ORFs: exact binomial {fp(c71.p_exact_binomial_M_plus_N_vs_length_share)}. Expected count for M is <5.',ha='center',fontsize=7)
lab(ax,'C',-0.14); fig.suptitle('R6-C  Long alternative-frame ORF capacity per canonical ORF',fontsize=8,y=0.995); fig.tight_layout(rect=(0,0.12,1,0.96)); finish(fig,'R6-C_capacity_per_ORF')
print('ok')
