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
from s3_plot_lib import *
from s3_lib import load_all, CLS3
M=load_all()
A_T=pd.read_csv(RES+'panelA_table.tsv',sep='\t')
B_T=pd.read_csv(RES+'panelB_spearman.tsv',sep='\t')
CF=pd.read_csv(RES+'panelC_frequencies.tsv',sep='\t')
D1=pd.read_csv(RES+'sim_curves_D1.tsv',sep='\t')
D2=pd.read_csv(RES+'sim_D2_points.tsv',sep='\t')
SIM=pd.read_csv(RES+'sim_effect_matches.tsv',sep='\t')
MHM=pd.read_csv(RES+'mh_main.tsv',sep='\t')
MH_E1=pd.read_csv(RES+'mh_E1_points.tsv',sep='\t')
MHC=pd.read_csv(RES+'mh_canonical_position_null.tsv',sep='\t')
MHS=pd.read_csv(RES+'mh_strata_samples.tsv',sep='\t')
MFM=pd.read_csv(RES+'mfe_main.tsv',sep='\t')
MF1=pd.read_csv(RES+'mfe_F1_points.tsv',sep='\t')
MF2=pd.read_csv(RES+'mfe_F2_spearman.tsv',sep='\t')
SEP_MFE=pd.read_csv(W+'separability_mfe.tsv',sep='\t')
GREY='#444444'
def legend_strip(c,y=1.0,x0=0.5,maxw=1e9):
    items=[('canonical','Canonical'),('dep','TRS-L-dep noncanonical'),('distant','TRS-L-indep distant'),('local','TRS-L-indep local'),('bg','Pair-/composition-shuffled background')]
    x=x0
    for k,lab in items:
        wlab=3.6+len(lab)*1.32+4.5
        if x+wlab>maxw+4.5 and x>x0: x=x0; y+=4.0
        c.fig.add_artist(Rectangle(((c.ox+x)/c.W,1-(c.oy+y+3.2)/c.H),2.6/c.W,2.6/c.H,fc=COL[k],ec='none',transform=c.fig.transFigure))
        c.text(x+3.6,y+1.4,lab,va='center',ha='left')
        x+=3.6+len(lab)*1.32+4.5
def note(c,x,y,s,**kw):
    kw.setdefault('color',GREY); c.text(x,y,s,va='top',ha='left',**kw)
# ------------------------------------------------------------------ A
def draw_A(c):
    c.letter('A'); c.title('Canonical TRS-L / TRS-B architecture',x=7)
    ax=c.ax(0,8,112,42); ax.set_xlim(0,112); ax.set_ylim(42,0); ax.axis('off')
    cw=3.25; x0=12.5
    ref=A_T.ref_12mer[0]
    ax.text(x0+3*cw,1.6,'core CTAAAC (67–72)',ha='left',va='center',fontsize=FS,color=GREY)
    ax.plot([x0+3*cw,x0+9*cw],[3.1,3.1],color=GREY,lw=0.8)
    yref=6.6
    ax.text(0.5,yref,'TRS-L',ha='left',va='center',fontsize=FS,fontweight='bold')
    for i,ch in enumerate(ref):
        core=3<=i<=8
        ax.add_patch(Rectangle((x0+i*cw,yref-2.1),cw,4.2,fc='#D9D9D9' if core else 'white',ec='#BDBDBD',lw=0.4))
        ax.text(x0+(i+.5)*cw,yref,ch,ha='center',va='center',fontsize=8,fontweight='bold')
    heads=['Matched\n(of 12)','Longest\nrun (nt)','Core\nmism.','MFE (a)\nkcal/mol','Reads\n(% can.)']
    xc=[x0+12*cw+2.5+i*11.4 for i in range(5)]
    for x,h in zip(xc,heads): ax.text(x+5.7,3.6,h,ha='center',va='center',fontsize=FS,linespacing=1.05)
    ax.plot([0,112],[9.6,9.6],color='black',lw=0.6)
    yb=[]
    for r,(i,row) in enumerate(A_T.iterrows()):
        y=13.4+r*5.7; g=row.gene; und=g=='ORF3'
        ax.text(0.5,y,g,ha='left',va='center',fontsize=8,fontweight='bold')
        for k,ch in enumerate(row.body_12mer):
            mt=ch==ref[k]
            fc=GENE_COL[g] if mt else 'white'
            if und: fc='#EDEDED' if mt else 'white'
            ax.add_patch(Rectangle((x0+k*cw,y-2.4),cw,4.8,fc=fc,ec='#BDBDBD',lw=0.4,hatch='////' if und else None))
            tc='white' if (mt and g in('E','M','N')) else ('#777777' if und else 'black')
            ax.text(x0+(k+.5)*cw,y,ch,ha='center',va='center',fontsize=8,color=tc)
        cm='n.d.' if und else str(int(row.core_mismatch_vs_CTAAAC))
        vals=[str(int(row.matched_of_12)),str(int(row.longest_contiguous)),cm,'%.1f'%row.mfe_a_kcal,'%.1f'%(100*row.canonical_share)]
        for x,v in zip(xc,vals): ax.text(x+5.7,y,v,ha='center',va='center',fontsize=FS+0.5)
    note(c,0.5,50.3,'Body windows aligned at the TRS-B core start (S GTAAAC differs from CTAAAC at its first base).')
    note(c,0.5,53.3,'ORF3: TRS-B core not resolved — window placed at breakpoint 24718 (hatched; "n.d.").')
    return A_T
# ------------------------------------------------------------------ B
def draw_B(c):
    c.letter('B'); c.title('Pairing strength vs canonical abundance',x=7)
    ax=c.ax(13,10,50,21)
    for _,r in A_T.iterrows():
        ax.scatter(r.mfe_a_kcal,r.canonical_reads,s=34,color=GENE_COL[r.gene],edgecolor='black',lw=0.5,zorder=3)
        ax.annotate(r.gene,(r.mfe_a_kcal,r.canonical_reads),xytext=(4,3),textcoords='offset points',fontsize=FS)
    ax.set_yscale('log'); ax.set_ylim(3e3,6e6); ax.set_xlim(-21,-10)
    ax.set_xlabel('MFE (a) (kcal/mol)',fontsize=7.5,labelpad=1); ax.set_ylabel('Canonical junction reads',fontsize=7.5)
    ra=B_T[B_T.column=='mfe_a_kcal'].iloc[0]
    ax.text(0.03,0.04,'Spearman ρ = %+.2f\nexact P = %.2f (n = 5)'%(ra.rho,ra.P_exact_two_sided),transform=ax.transAxes,fontsize=FS,va='bottom',linespacing=1.05)
    ax.grid(color='#EEEEEE',lw=0.4)
    tx=c.ax(0,39,68,18); tx.axis('off'); tx.set_xlim(0,68); tx.set_ylim(18,0)
    tx.text(0.5,1.4,'Metric (canonical, n = 5)',fontsize=FS,fontweight='bold',va='center'); tx.text(46,1.4,'ρ',fontsize=FS,fontweight='bold',ha='center',va='center'); tx.text(59,1.4,'exact P',fontsize=FS,fontweight='bold',ha='center',va='center')
    tx.plot([0,68],[3.0,3.0],color='black',lw=0.5)
    names={'maxid_matches':'Best-shift identity','longest_contig_at_best_shift':'Longest contiguous match','mfe_a_kcal':'MFE (a) — plotted','mfe_b_kcal':'MFE (b), direct pairing'}
    for i,(col,nm) in enumerate(names.items()):
        r=B_T[B_T.column==col].iloc[0]; y=5.3+i*2.9
        tx.text(0.5,y,nm,fontsize=FS,va='center'); tx.text(46,y,'%+.2f'%r.rho,fontsize=FS,ha='center',va='center'); tx.text(59,y,'%.2f'%r.P_exact_two_sided,fontsize=FS,ha='center',va='center')
    tx.text(0.5,17.0,'Smallest attainable two-sided P at n = 5: 0.017',fontsize=FS,va='center',color=GREY)
    return A_T[['gene','mfe_a_kcal','canonical_reads']]
# ------------------------------------------------------------------ C
def draw_C(c):
    c.letter('C'); c.title('Acceptor sequence context',x=7)
    labs=[('canonical','Canonical (n = 5)'),('TRSL_dep_noncanonical','TRS-L-dep noncanonical (n = 906)'),('TRSL_indep_distant','TRS-L-indep distant (n = 1,688)'),
          ('TRSL_indep_local','TRS-L-indep local (n = 1,479)'),('background_region_matched_positions','Region-matched random positions (n = {:,})')]
    nn=CF.groupby('group').n.first()
    x=44
    for nt in 'ACGU':
        c.fig.add_artist(Line2D([(c.ox+x)/c.W,(c.ox+x+4)/c.W],[1-(c.oy+5.0)/c.H]*2,color=NT_COL[nt],ls=NT_LS[nt],lw=1.2,transform=c.fig.transFigure))
        c.text(x+5,5.0,nt,va='center'); x+=8.5
    axes=[]
    for i,(g,lab) in enumerate(labs):
        y=9.0+i*10.0
        ax=c.ax(13,y+3.4,62,6.0)
        sub=CF[CF.group==g]
        for nt in 'ACGU':
            s=sub[sub.nucleotide==nt].sort_values('position')
            ax.plot(s.position,s.frequency,color=NT_COL[nt],ls=NT_LS[nt],lw=1.0)
        ax.axvline(0,color='#555555',lw=0.7,ls='--'); ax.set_xlim(-25,15); ax.set_ylim(0,1); ax.set_yticks([0,.5,1]); ax.set_yticklabels(['0','.5','1'])
        if i<4: ax.set_xticklabels([])
        c.text(13,y+1.6,lab.format(int(nn[g])) if '{' in lab else lab.replace('{:,}',''),va='center',ha='left',fontsize=FS)
        axes.append(ax)
    axes[-1].set_xlabel('Position relative to acceptor (nt)',fontsize=7.5,labelpad=1)
    c.text(2.4,33,'Nucleotide frequency',rotation=90,ha='center',va='center',fontsize=7.5)
    return CF
# ------------------------------------------------------------------ D
def draw_D(c):
    c.letter('D'); c.title('TRS-L similarity around the acceptor',x=7)
    c.text(11,9.6,'D1  position-dependent',fontsize=7.5,fontweight='bold',va='center')
    c.text(69,9.6,'D2  best shift (−20 to +20)',fontsize=7.5,fontweight='bold',va='center')
    ax=c.ax(11,17,48,37)
    sub=D1
    for g,key in [('canonical','canonical'),('TRSL_dep_noncanonical','dep'),('TRSL_indep_distant','distant'),('TRSL_indep_local','local')]:
        s=sub[(sub.group==g)&(sub.kind=='observed')].sort_values('offset')
        ax.plot(s.offset,s.mean_matches,color=COL[key],lw=1.3 if key!='canonical' else 1.0); ax.fill_between(s.offset,s.ci_lo,s.ci_hi,color=COL[key],alpha=0.12 if key=='canonical' else 0.25,lw=0)
    bgs=[sub[(sub.group==g)&(sub.kind=='composition-shuffled background')].sort_values('offset') for g in CLS3]
    yb=np.mean([b.mean_matches.values for b in bgs],axis=0); ax.plot(bgs[0].offset,yb,color=COL['bg'],lw=1.2,ls='--')
    ax.axvline(0,color='#555555',lw=0.7,ls='--'); ax.set_xlim(-20,20); ax.set_ylim(0,8); ax.set_yticks(range(0,9,2))
    ax.set_xlabel('Window offset relative to acceptor (nt)',fontsize=7.5,labelpad=1); ax.set_ylabel('Mean matched bases (of 12)',fontsize=7.5)
    ax.text(0.03,0.03,'bands: bootstrap 95% CI',transform=ax.transAxes,fontsize=FS,color=GREY)
    a2=c.ax(68,17,30,37)
    order=['canonical','TRSL_dep_noncanonical','TRSL_indep_distant','TRSL_indep_local']
    rng=np.random.default_rng(1)
    for i,g in enumerate(order):
        d=D2[D2.group==g]; key=GKEY[g]
        if g=='canonical':
            for _,r in d.iterrows(): a2.scatter(i+rng.uniform(-.12,.12),r.best_shift_matches,s=16,color=GENE_COL[r.canonical_gene],edgecolor='black',lw=0.4,zorder=5)
            a2.plot([i-.3,i+.3],[d.best_shift_matches.mean()]*2,color='black',lw=0.9)
        else:
            violin(a2,d.best_shift_matches,i,COL[key],width=0.8)
    violin(a2,D2[D2.group!='canonical'].background_one_shuffle,4,COL['bg'],width=0.8)
    a2.set_xticks(range(5)); a2.set_xticklabels(['Can.','Dep.','Dist.','Loc.','Bkg.']); a2.set_ylim(0,12.5); a2.set_yticks([0,3,6,9,12]); a2.set_xlim(-.6,4.6)
    a2.axhline(SIM[(SIM.group=='TRSL_dep_noncanonical')&SIM.background.str.startswith('comp')].bg_mean_matches.iloc[0],color=COL['bg'],lw=0.7,ls='--')
    a2.set_ylabel('Matched bases (of 12)',fontsize=7.5,labelpad=1)
    e=SIM[SIM.background.str.startswith('comp')].set_index('group')
    r=e.loc['canonical']
    a2.text(0,1.03,'Δ %+.1f'%(r.delta_matches),transform=a2.get_xaxis_transform(),ha='center',va='bottom',fontsize=FS,color=COL['canonical'],fontweight='bold',linespacing=1.0)
    rg=[e.loc[g].delta_matches for g in CLS3]
    a2.text(2.9,1.03,'Δ %+.2f to %+.2f\nnegligible effect size'%(min(rg),max(rg)),transform=a2.get_xaxis_transform(),ha='center',va='bottom',fontsize=FS,color='#333333',linespacing=1.0,style='italic')
    return pd.concat([D1.assign(panel='D1'),D2.assign(panel='D2')],ignore_index=True)
# ------------------------------------------------------------------ E
def draw_E(c):
    c.letter('E'); c.title('Donor–acceptor microhomology',x=7)
    c.text(12,8.6,'E1  by class',fontsize=7.5,fontweight='bold',va='center'); c.text(65,8.6,'E2  by recurrence',fontsize=7.5,fontweight='bold',va='center')
    ax=c.ax(12,21,48,26)
    order=[('canonical','Can.'),('TRSL_dep_noncanonical','Dep.'),('TRSL_indep_distant','Dist.'),('TRSL_indep_local','Loc.')]
    mm=MHM[MHM.subset=='filtered (main)'].set_index('group')
    rng=np.random.default_rng(3)
    for i,(g,lab) in enumerate(order):
        x=i
        if g=='canonical':
            d=M[M.canonical_primary]
            for _,r in d.iterrows(): ax.scatter(x-0.18+rng.uniform(-.08,.08),r.mh_len,s=16,color=GENE_COL[r.canonical_gene],edgecolor='black',lw=0.4,zorder=5)
            ref=MHC[(MHC.group=='canonical_primary5')&MHC.background.str.startswith('uniform')].iloc[0]
            ax.plot([x+0.05,x+0.4],[ref.null_mean]*2,color=COL['bg'],lw=1.4); ax.plot([x+0.22]*2,[ref.null_p2_5,ref.null_p97_5],color=COL['bg'],lw=1.0)
            delta=d.mh_len.mean()-ref.null_mean; fold=d.mh_len.mean()/ref.null_mean; p=ref.P
            nt_='%s\n%+.1f nt\n(%.1f×)'%(star(p),delta,fold)
        else:
            key=GKEY[g]
            dobs=MH_E1[(MH_E1.group==g)&(MH_E1.kind=='observed')].mh_len; dbg=MH_E1[(MH_E1.group==g)&(MH_E1.kind!='observed')].mh_len
            violin(ax,dobs,x-0.19,COL[key],width=0.36); violin(ax,dbg,x+0.19,COL['bg'],width=0.36)
            r=mm.loc[g]; ph=r.P_mean_holm_across_3_classes
            nt_='%s\n%+.2f nt\n(%.1f×)'%(star(ph),r.delta_nt,r.fold_mean)
            if g=='TRSL_indep_distant':
                r2=MHM[MHM.subset=='excl. donor<=100'].iloc[0]
                ax.scatter([x-0.19],[r2.obs_mean_mh],marker='D',s=14,facecolor='white',edgecolor='black',lw=0.7,zorder=6)
                ax.scatter([x+0.19],[r2.null_mean_mh],marker='D',s=14,facecolor='white',edgecolor='black',lw=0.7,zorder=6)
        ax.text(x,1.02,nt_,transform=ax.get_xaxis_transform(),ha='center',va='bottom',fontsize=FS,linespacing=1.0)
    ax.set_xticks(range(4)); ax.set_xticklabels([l for _,l in order]); ax.set_xlim(-.6,3.6); ax.set_ylim(0,14)
    ax.set_ylabel('Microhomology (nt)',fontsize=7.5,labelpad=1)
    note(c,0.5,53.6,'Colored: observed; gray: pair-shuffled.')
    note(c,0.5,56.6,'Open diamond: distant excl. donor ≤ 100.')
    a2=c.ax(67,21,21,26)
    bins=['1','2','3-4','>=5']; xs=np.arange(4)
    for g,key in [('TRSL_dep_noncanonical','dep'),('TRSL_indep_distant','distant'),('TRSL_indep_local','local')]:
        s=MHS[MHS.group==g].set_index('samples_bin').loc[bins]
        off={'dep':-0.09,'distant':0,'local':0.09}[key]
        a2.errorbar(xs+off,s.delta_nt,yerr=[s.delta_nt-s.delta_ci_lo,s.delta_ci_hi-s.delta_nt],color=COL[key],marker='o',ms=3.3,lw=1.0,capsize=1.5,elinewidth=0.8)
    a2.axhline(0,color='#555555',lw=0.7,ls='--'); a2.set_xticks(xs); a2.set_xticklabels(['1','2','3–4','≥5']); a2.set_xlim(-.4,3.4)
    a2.set_xlabel('Positive samples',fontsize=7.5,labelpad=1); a2.set_ylabel('Δ vs stratum shuffle (nt)',fontsize=7.5,labelpad=1)
    a2.set_ylim(-0.2,2.2)
    a2.text(0.5,1.02,'all events\n(unfiltered)',transform=a2.transAxes,ha='center',va='bottom',fontsize=FS,color=GREY,linespacing=1.0)
    return MH_E1
# ------------------------------------------------------------------ F
def draw_F(c):
    c.letter('F'); c.title('RNA–RNA interaction strength',x=7)
    c.text(13,8.6,'F1  MFE (b), direct pairing',fontsize=7.5,fontweight='bold',va='center'); c.text(63,8.6,'F2  vs reads',fontsize=7.5,fontweight='bold',va='center')
    ax=c.ax(13,21,44,26)
    order=[('canonical','Can.'),('TRSL_dep_noncanonical','Dep.'),('TRSL_indep_distant','Dist.'),('TRSL_indep_local','Loc.')]
    mb=MFM[(MFM['mode']=='b')&(MFM.subset=='filtered (main)')].set_index('group')
    rng=np.random.default_rng(4); nclip=0
    for i,(g,lab) in enumerate(order):
        x=i
        if g=='canonical':
            d=M[M.canonical_primary]
            for _,r in d.iterrows(): ax.scatter(x-0.18+rng.uniform(-.08,.08),r.mfe_b,s=16,color=GENE_COL[r.canonical_gene],edgecolor='black',lw=0.4,zorder=5)
            ref=SEP_MFE[(SEP_MFE.group=='canonical_primary5')&(SEP_MFE['mode']=='b')&SEP_MFE.null_type.str.startswith('acceptor')].iloc[0]
            ax.plot([x+0.05,x+0.4],[ref.null_mean]*2,color=COL['bg'],lw=1.4); ax.plot([x+0.22]*2,[ref.null_p2_5,ref.null_p97_5],color=COL['bg'],lw=1.0)
            nt_='ns\n%+.1f'%(d.mfe_b.mean()-ref.null_mean)
        else:
            key=GKEY[g]
            dobs=MF1[(MF1.group==g)&(MF1['mode']=='b')&(MF1.kind=='observed')].mfe.values; dbg=MF1[(MF1.group==g)&(MF1['mode']=='b')&(MF1.kind!='observed')].mfe.values
            nclip+=int((dobs<-30).sum()+(dbg<-30).sum())
            violin(ax,np.clip(dobs,-30,None),x-0.19,COL[key],width=0.36); violin(ax,np.clip(dbg,-30,None),x+0.19,COL['bg'],width=0.36)
            r=mb.loc[g]; nt_='%s\n%+.2f'%(star(r.P_holm_across_3_classes),r.delta_kcal)
            if g=='TRSL_indep_distant':
                r2=MFM[(MFM['mode']=='b')&(MFM.subset=='excl. donor<=100')].iloc[0]
                ax.scatter([x-0.19],[r2.obs_mean_mfe],marker='D',s=14,facecolor='white',edgecolor='black',lw=0.7,zorder=6)
                ax.scatter([x+0.19],[r2.null_mean],marker='D',s=14,facecolor='white',edgecolor='black',lw=0.7,zorder=6)
        ax.text(x,1.02,nt_,transform=ax.get_xaxis_transform(),ha='center',va='bottom',fontsize=FS,linespacing=1.0)
    ax.set_xticks(range(4)); ax.set_xticklabels([l for _,l in order]); ax.set_xlim(-.6,3.6); ax.set_ylim(-30,0)
    ax.set_ylabel('MFE (b) (kcal/mol)',fontsize=7.5,labelpad=1)
    note(c,0.5,53.6,'More negative = stronger pairing.')
    note(c,0.5,56.6,'Open diamond: distant excl. donor ≤ 100; %d value(s) < −30 shown at −30.'%nclip)
    a2=c.ax(66,21,22,26)
    f2=MF2[MF2['mode']=='b']; lines=[]
    for g,key,lab in [('TRSL_dep_noncanonical','dep','Dep.'),('TRSL_indep_distant','distant','Dist.'),('TRSL_indep_local','local','Loc.')]:
        d=M[(M.cls==g)&M.pass_filter]
        a2.scatter(d.mfe_b,d.reads,s=2.5,color=COL[key],alpha=0.35,lw=0,rasterized=True)
        fl=d[d.flag_leaderdonor_at_canon_site_in_distant]
        if len(fl): a2.scatter(fl.mfe_b,fl.reads,s=12,facecolor='none',edgecolor='black',lw=0.6,marker='s',zorder=5)
        r=f2[(f2.group==g)&(f2.subset=='filtered (main)')].iloc[0]; lines.append((lab,COL[key],r.rho))
    cp=M[M.canonical_primary]
    a2.scatter(cp.mfe_b,cp.reads,s=18,marker='D',color=COL['canonical'],edgecolor='black',lw=0.4,zorder=6)
    rc=f2[f2.group=='canonical_primary5'].iloc[0]; lines.insert(0,('Can.',COL['canonical'],rc.rho))
    a2.set_yscale('log'); a2.set_ylim(0.7,4e6); a2.set_xlim(-30,0)
    a2.set_xlabel('MFE (b) (kcal/mol)',fontsize=7.5,labelpad=1); a2.set_ylabel('Junction reads',fontsize=7.5,labelpad=1)
    for k,(lab,col,rho) in enumerate(lines):
        c.text(67,11.6+k*2.75,'%s ρ %+.2f'%(lab,rho),color=col,fontsize=FS,va='center',ha='left',fontweight='bold')
    c.text(67,19.0,'□ = leader-donor events at canonical-site acceptors',fontsize=FS,va='center',ha='left',color=GREY) if False else None
    return pd.DataFrame(dict(x=[0]))
# ------------------------------------------------------------------ G
def rbox(ax,x,y,w,h,fc,ec='#555555',ls='-',lw=0.7,r=0.6):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0,rounding_size=%g'%r,fc=fc,ec=ec,ls=ls,lw=lw))
def draw_G(c):
    c.letter('G'); c.title('Working model supported by A–F',x=7)
    ax=c.ax(0,7,180,37); ax.set_xlim(0,180); ax.set_ylim(37,0); ax.axis('off')
    rbox(ax,0.5,0.4,179,10.4,'#E8F0F8')
    ax.text(2,2.9,'1  Canonical: TRS-L/TRS-B complementarity',fontsize=7.5,fontweight='bold',va='center')
    ax.text(2,7.7,'matched bases 8.8 vs 6.6 in background (Δ +2.2); core mismatch 0–1 (ORF3 unresolved)',fontsize=FS,va='center',color='#333333')
    ax.add_patch(Rectangle((104,6.4),4,2.6,fc='#B8860B',ec='none')); ax.plot([108,168],[7.7]*2,color='#777777',lw=0.7)
    for cx,wd,fc in [(122,6,'#E1EBF5'),(136,5,'#EAE3F3'),(150,8,'#E0F0F0')]: ax.add_patch(Rectangle((cx,6.4),wd,2.6,fc=fc,ec='#9A9A9A',lw=0.4))
    ax.annotate('',xy=(124,6.2),xytext=(106,6.2),arrowprops=dict(arrowstyle='-|>',color=COL['canonical'],lw=1.0,connectionstyle='arc3,rad=-0.6'))
    ax.text(128,2.7,'TRS-L pairs with TRS-B copy',fontsize=FS,ha='left',va='center',color=COL['canonical'])
    rbox(ax,0.5,11.8,179,13.6,'#F3F3F3')
    ax.text(2,14.2,'2  Noncanonical junctions carry elevated microhomology; origin undetermined',fontsize=7.5,fontweight='bold',va='center')
    ax.text(2,19.3,'Δ +0.8 to +1.1 nt vs pair-shuffled; rises with recurrence, not depth',fontsize=FS,va='center',color='#333333')
    ax.text(2,22.9,'RT/PCR (library-preparation) chimeras not excluded',fontsize=FS,va='center',color='#333333')
    rbox(ax,100,16.6,38,7.6,'white',ec=COL['dep'])
    ax.text(119,18.9,'Biological',fontsize=FS,ha='center',va='center',fontweight='bold',color=COL['dep']); ax.text(119,22.0,'template switching',fontsize=FS,ha='center',va='center')
    rbox(ax,146,16.6,32,7.6,'white',ec='#555555')
    ax.text(162,18.9,'Library prep.',fontsize=FS,ha='center',va='center',fontweight='bold',color='#555555'); ax.text(162,22.0,'RT/PCR chimera',fontsize=FS,ha='center',va='center')
    ax.text(142,20.4,'or ?',fontsize=8,ha='center',va='center',fontweight='bold')
    rbox(ax,0.5,26.6,179,10.0,'white',ec='#555555',ls=(0,(4,2)))
    ax.text(2,29.0,'3  Candidate (not tested here): RNA structure / long-range proximity  ?',fontsize=7.5,fontweight='bold',va='center')
    ax.text(2,33.9,'MFE (b) differs from background only in the local class (Δ −0.7 kcal/mol)',fontsize=FS,va='center',color='#333333')
    ax.add_patch(Rectangle((112,33.6),4,2.2,fc='#B8860B',ec='none')); ax.plot([116,168],[34.7]*2,color='#777777',lw=0.7)
    ax.annotate('',xy=(150,33.4),xytext=(120,33.4),arrowprops=dict(arrowstyle='-',color='#555555',lw=0.9,ls=(0,(3,2)),connectionstyle='arc3,rad=-0.5'))
    ax.text(135,30.6,'?',fontsize=9,ha='center',va='center',fontweight='bold')
    return pd.DataFrame(dict(note=['schematic; no data']))
