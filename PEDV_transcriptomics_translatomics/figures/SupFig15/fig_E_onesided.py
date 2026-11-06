"""R6-E redone: translation efficiency (Ribo P-site density / RNA relative depth) vs long alt-frame ORF capacity, ORF3 included (n=6: ORF1a, ORF1b, S, ORF3, M, N; E excluded),
ONE-SIDED tests of the directional hypothesis 'higher translation efficiency -> more long alternative-frame ORFs per kb' (positive correlation). Exact permutation (720 permutations); Pearson one-sided on log10(TE)."""
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
import itertools, os
import matplotlib; matplotlib.use('Agg'); import matplotlib.ticker
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
from r6lib import *
plt.rcParams.update({'svg.fonttype':'none','pdf.fonttype':42,'font.family':'DejaVu Sans','font.size':7.5,'axes.titlesize':8,'axes.labelsize':8,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'legend.fontsize':7,
  'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':0.7,'xtick.major.width':0.7,'ytick.major.width':0.7})
RCOL={"ORF1a":"#0072B2","ORF1b":"#56B4E9","S":"#009E73","ORF3":"#CC79A7","E":"#D55E00","M":"#E69F00","N":"#117733"}
SD=OUT+'source_data/'; FD=OUT+'figures/'
P=pd.read_csv(SD+'G4_activity_vs_capacity_table.csv').set_index('region'); ORFS=['ORF1a','ORF1b','S','ORF3','M','N']
TE={'TE_CHX':P.CHX_density_ppm_per_nt/P.RNA_rel_depth_batch3_5_8_16h_median,'TE_HAR':P.HAR_P_per_nt/P.RNA_rel_depth_batch3_5_8_16h_median}
f=families(); ia=f[(f.primary_category=='internal alternative-frame ORF')&f.relative_frame_to_parent.isin([1,2])]
def fmt(p): return f'{p:.3f}' if p<0.1 else f'{p:.2f}'
rows=[]; rate={}
for T in (20,30,40):
    sub=ia[ia.max_peptide_aa>T]; k=np.array([(sub.fixed_parent_assignment==r).sum() for r in ORFS]); L=np.array([CDS[r][1]-CDS[r][0]+1 for r in ORFS]); rate[T]=1000*k/L
    for nm,te in TE.items():
        x=te.loc[ORFS].values; y=rate[T]; rho=stats.spearmanr(x,y).statistic
        perms=np.array([stats.spearmanr(x,p).statistic for p in itertools.permutations(y)])
        p1=float(np.mean(perms>=rho-1e-12)); p2=float(np.mean(np.abs(perms)>=abs(rho)-1e-12))
        pr=stats.pearsonr(np.log10(x),y); pp1=pr.pvalue/2 if pr.statistic>0 else 1-pr.pvalue/2
        rows.append(dict(threshold=f'>{T} aa',n_families=int(k.sum()),predictor=nm,n_ORFs=6,spearman_rho=round(rho,3),p_one_sided_exact=round(p1,4),p_two_sided_exact=round(p2,4),pearson_r_log10x=round(pr.statistic,3),p_pearson_one_sided=round(pp1,4),min_attainable_one_sided_p=round(1/720,4)))
D=pd.DataFrame(rows); D.to_csv(SD+'R6-E_TE_correlation_ORF3_included_one_sided.csv',index=False,encoding='utf-8-sig')
pd.DataFrame({'region':ORFS,'length_nt':[CDS[r][1]-CDS[r][0]+1 for r in ORFS],'TE_CHX':TE['TE_CHX'].loc[ORFS].values,'TE_HAR':TE['TE_HAR'].loc[ORFS].values,'fam_per_kb_gt20':rate[20],'fam_per_kb_gt30':rate[30],'fam_per_kb_gt40':rate[40]}).to_csv(SD+'R6-E_panel_data_ORF3_included.csv',index=False,encoding='utf-8-sig')
print(D.to_string())
fig,axs=plt.subplots(1,3,figsize=(9.8,3.6),gridspec_kw={'width_ratios':[1,1,1]})
for ax,nm,xl in zip(axs[:2],('TE_CHX','TE_HAR'),('TE (CHX Ribo density / RNA depth)','TE (HAR Ribo density / RNA depth)')):
    r=D[(D.threshold=='>40 aa')&(D.predictor==nm)].iloc[0]
    for o in ORFS:
        ax.scatter(TE[nm][o],rate[40][ORFS.index(o)],s=36,color=RCOL[o],edgecolor='k',lw=0.4,zorder=3)
        ax.annotate(o,(TE[nm][o],rate[40][ORFS.index(o)]),xytext=(4,-10) if o=='ORF3' else (4,3),textcoords='offset points',fontsize=7.5)
    ax.set_xscale('log'); ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v,p:'%g'%v)); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel(xl); ax.set_ylim(1.5,6.4)
    ax.text(0.04,0.97,f'Spearman ρ={r.spearman_rho:.2f}, one-sided p={fmt(r.p_one_sided_exact)} (exact)\nPearson r(log10 x)={r.pearson_r_log10x:.2f}, one-sided p={fmt(r.p_pearson_one_sided)}\nn=6 ORFs; min attainable one-sided p=0.0014',transform=ax.transAxes,va='top',fontsize=7)
axs[0].set_ylabel('>40-aa +1/+2 families per kb')
ax=axs[2]
for j,(nm,col) in enumerate((('TE_CHX','#1f77b4'),('TE_HAR','#d62728'))):
    for tt,T in enumerate((20,30,40)):
        r=D[(D.threshold==f'>{T} aa')&(D.predictor==nm)].iloc[0]; xx=tt+(j-0.5)*0.3; ax.bar(xx,r.spearman_rho,0.26,color=col,lw=0)
ax.axhline(0,color='k',lw=0.6); ax.set_xticks(range(3))
ax.set_xticklabels([f'>{T} aa\np='+'/'.join(fmt(D[(D.threshold==f">{T} aa")&(D.predictor==n)].iloc[0].p_one_sided_exact) for n in ('TE_CHX','TE_HAR')) for T in (20,30,40)],fontsize=7)
ax.set_ylim(-1,1.7); ax.set_ylabel('Spearman ρ with families per kb'); ax.set_title('Length threshold sensitivity (n=6)')
ax.legend(handles=[plt.Rectangle((0,0),1,1,color='#1f77b4',label='TE (CHX)'),plt.Rectangle((0,0),1,1,color='#d62728',label='TE (HAR)')],frameon=False,loc='upper right',ncol=2,fontsize=7)
for a_,t_ in zip(axs,'ABC'): a_.text(-0.16,1.06,t_,transform=a_.transAxes,fontsize=10,fontweight='bold',va='bottom')
fig.text(0.5,0.0,'One-sided hypothesis: higher translation efficiency, more long alt-frame ORFs per kb. ORF3 (276 nt, 1 family) included; E excluded.\nRibo-seq: strict_core 29/30 nt P-sites; RNA: batch3 5/8/16 h relative depth (non-strand-specific); N/M coverage includes nested sgRNA sharing.\nC: one-sided exact permutation p (CHX/HAR); no correction across predictors or thresholds.',ha='center',fontsize=7)
fig.suptitle('R6-E  Translation efficiency vs long alternative-frame ORF capacity (ORF1a, ORF1b, S, ORF3, M, N), one-sided',fontsize=8,y=0.995)
fig.tight_layout(rect=(0,0.16,1,0.95))
for a_ in fig.axes:
    for c in a_.get_children():
        try: c.set_clip_on(False)
        except Exception: pass
for ext,kw in (('pdf',{}),('svg',{}),('png',{'dpi':300})): fig.savefig(FD+'R6-E_translation_efficiency_vs_capacity.'+ext,**kw)
