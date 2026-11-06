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
import sys; sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from s4_panels import *
OUT=B+'figures/'; DAT=OUT+'Fig3_data/'
# ---- S3-1: Panel B, five metrics; vertical CI only (no horizontal error bars: metrics are deterministic per gene)
fig,axs=plt.subplots(1,5,figsize=(180/25.4,50/25.4))
cols=[('matched_of_12','Matched bases (of 12)'),('longest_contiguous','Longest contiguous run (nt)'),('mfe_a_kcal','MFE (a) (kcal/mol)'),('mfe_b_kcal','MFE (b) (kcal/mol)'),('maxid_matches','Best-shift matched bases')]
for ax,(col,lab) in zip(axs,cols):
    for g in GENES:
        r=A_T[A_T.gene==g].iloc[0]; p=PSR.loc[g]
        ax.errorbar(r[col],p.mean_reads_per_sample,yerr=[[p.mean_reads_per_sample-p.ci_lo],[p.ci_hi-p.mean_reads_per_sample]],fmt='o',ms=4,mfc=GENE_COL[g],mec='black',mew=0.5,ecolor='black',elinewidth=0.6,capsize=1.2)
        ax.annotate(g,(r[col],p.mean_reads_per_sample),xytext=(3,3),textcoords='offset points',fontsize=FS)
    ax.set_yscale('log'); ax.set_xlabel(lab,fontsize=7.5,labelpad=1); ax.grid(color='#EEEEEE',lw=0.4)
    r_=B3[B3.column==col].iloc[0]; ax.set_title('\u03c1 %+.2f, P %.2f'%(r_.rho,r_.P_exact_two_sided),fontsize=7.5)
axs[0].set_ylabel('Reads / sample (mean, 95% CI)',fontsize=7.5)
fig.subplots_adjust(left=0.07,right=0.99,bottom=0.26,top=0.86,wspace=0.35); savefig(fig,OUT+'FigS3_1_B_metrics'); plt.close(fig)
# ---- S3-5: unstratified vs span-stratified backgrounds (effect sizes)
MH_old=MHM[MHM.subset=='filtered (main)'].set_index('group'); MH_new=MHS3[MHS3.subset=='filtered (main)'].set_index('group')
MF_old=MFM[MFM.subset=='filtered (main)']; MF_new=MFS3[MFS3.subset=='filtered (main)']
fig,axs=plt.subplots(1,3,figsize=(180/25.4,52/25.4))
def pair(ax,rows,xl,ttl):
    for k,(lab,key,o,n) in enumerate(rows):
        ax.errorbar(o[0],k-0.12,xerr=[[o[0]-o[1]],[o[2]-o[0]]],fmt='o',ms=4,color=COL[key],mfc='white',capsize=1.5,lw=0.9)
        ax.errorbar(n[0],k+0.12,xerr=[[n[0]-n[1]],[n[2]-n[0]]],fmt='o',ms=4,color=COL[key],capsize=1.5,lw=0.9)
    ax.axvline(0,color='#555555',ls='--',lw=0.7); ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows]); ax.invert_yaxis(); ax.set_xlabel(xl,fontsize=7.5); ax.set_title(ttl,fontsize=7.5,loc='left')
gl=[('TRSL_dep_noncanonical','Dep.','dep'),('TRSL_indep_distant','Dist.','distant'),('TRSL_indep_local','Loc.','local')]
pair(axs[0],[(l,k,(MH_old.loc[g].delta_nt,MH_old.loc[g].delta_ci_lo,MH_old.loc[g].delta_ci_hi),(MH_new.loc[g].delta_nt,MH_new.loc[g].delta_ci_lo,MH_new.loc[g].delta_ci_hi)) for g,l,k in gl],'\u0394 microhomology (nt)','microhomology')
for ax,mode,ttl in [(axs[1],'b','MFE (b)'),(axs[2],'a','MFE (a)')]:
    o=MF_old[MF_old['mode']==mode].set_index('group'); n=MF_new[MF_new['mode']==mode].set_index('group')
    pair(ax,[(l,k,(o.loc[g].delta_kcal,o.loc[g].delta_ci_lo,o.loc[g].delta_ci_hi),(n.loc[g].delta_kcal,n.loc[g].delta_ci_lo,n.loc[g].delta_ci_hi)) for g,l,k in gl],'\u0394 MFE (kcal/mol)',ttl)
fig.text(0.5,0.01,'open = class-only pair shuffle; filled = shuffle within class and span quartile (95% CI)',ha='center',fontsize=FS)
fig.subplots_adjust(left=0.07,right=0.98,bottom=0.25,top=0.88,wspace=0.35); savefig(fig,OUT+'FigS3_5_background_strategies'); plt.close(fig)
# ---- S3-6: ECDF of microhomology and donor/acceptor positions
fig,axs=plt.subplots(1,3,figsize=(180/25.4,52/25.4))
for g,key in [('TRSL_dep_noncanonical','dep'),('TRSL_indep_distant','distant'),('TRSL_indep_local','local')]:
    o=np.sort(MH_E1[(MH_E1.group==g)&(MH_E1.kind=='observed')].mh_len.values); s=np.sort(pd.read_csv(RES+'mh_E1_spanshuf_%s.tsv'%g,sep='\t').mh_len.values)
    axs[0].step(o,np.arange(1,len(o)+1)/len(o),color=COL[key],lw=1.2); axs[0].step(s,np.arange(1,len(s)+1)/len(s),color=COL[key],lw=1.0,ls='--')
axs[0].set_xlim(0,8); axs[0].set_xlabel('Microhomology (nt)',fontsize=7.5); axs[0].set_ylabel('ECDF',fontsize=7.5); axs[0].set_title('solid observed; dashed shuffled',fontsize=7.5,loc='left')
for g,key in [('TRSL_dep_noncanonical','dep'),('TRSL_indep_distant','distant'),('TRSL_indep_local','local')]:
    d=M[(M.cls==g)&M.pass_filter]
    axs[1].hist(d.donor/1000,bins=np.arange(0,28.1,1),histtype='step',color=COL[key],lw=1.0); axs[2].hist(d.acceptor/1000,bins=np.arange(0,28.1,1),histtype='step',color=COL[key],lw=1.0)
axs[1].set_xlabel('Donor position (kb)',fontsize=7.5); axs[1].set_ylabel('Filtered junctions',fontsize=7.5); axs[2].set_xlabel('Acceptor position (kb)',fontsize=7.5); axs[1].set_yscale('log'); axs[2].set_yscale('log')
fig.subplots_adjust(left=0.07,right=0.99,bottom=0.24,top=0.88,wspace=0.35); savefig(fig,OUT+'FigS3_6_ECDF_and_positions'); plt.close(fig)
# ---- tables: read-weighted vs event-weighted (descriptive), unfiltered / Low-confidence observed means
rows=[]
for g in CLS3:
    for lab,sel in [('filtered (main)',M.pass_filter),('all events of class',M.pass_filter|True),('not passing filter',~M.pass_filter)]:
        d=M[(M.cls==g)&sel] if lab!='all events of class' else M[M.cls==g]
        if len(d)==0: continue
        rows.append(dict(group=g,set=lab,n=len(d),event_weighted_mean_mh=d.mh_len.mean(),read_weighted_mean_mh=np.average(d.mh_len,weights=d.reads),event_weighted_maxid_matches=(d.maxid*12).mean(),read_weighted_maxid_matches=np.average(d.maxid*12,weights=d.reads),
                         mean_mfe_b=d.mfe_b.mean(),mean_mfe_a=d.mfe_a.mean()))
pd.DataFrame(rows).to_csv(DAT+'supp_read_vs_event_weighted_and_unfiltered_descriptive.tsv',sep='\t',index=False)
print(pd.DataFrame(rows).round(3).to_string())
