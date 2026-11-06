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
from s3_lib import GENES
from s3_panels import *
from matplotlib.patches import FancyArrow
import matplotlib.colors as mcolors
B3=pd.read_csv(RES+'panelB_spearman_v3.tsv',sep='\t')
PSR=pd.read_csv(RES+'panelB_reads_per_sample.tsv',sep='\t').set_index('gene')
MHS3=pd.read_csv(RES+'mh_main_spanstrat.tsv',sep='\t')
MFS3=pd.read_csv(RES+'mfe_main_spanstrat.tsv',sep='\t')
MFP3=pd.read_csv(RES+'mfe_F1_points_spanstrat.tsv',sep='\t')
NTC={'A':'#4D4D4D','C':'#56B4E9','G':'#D55E00','U':'#CCCCCC'}
def fmt_reads(x):
    return '%.2fM'%(x/1e6) if x>=1e6 else ('%dk'%round(x/1e3) if x>=1e5 else '%.1fk'%(x/1e3))
# ---------------------------------------------------------------- A
def draw_A(c):
    c.letter('A'); c.title('Canonical TRS-L / TRS-B sequence architecture',x=7)
    ax=c.ax(0,8,112,42); ax.set_xlim(0,112); ax.set_ylim(42,0); ax.axis('off')
    cw=3.25; x0=12.5; ref=A_T.ref_12mer[0]
    ax.text(x0+3*cw,1.6,'core CTAAAC (67–72)',ha='left',va='center',fontsize=FS,color=GREY)
    ax.plot([x0+3*cw,x0+9*cw],[3.1,3.1],color=GREY,lw=0.8)
    yref=6.6
    ax.text(0.5,yref,'TRS-L',ha='left',va='center',fontsize=FS,fontweight='bold')
    for i,ch in enumerate(ref):
        core=3<=i<=8
        ax.add_patch(Rectangle((x0+i*cw,yref-2.1),cw,4.2,fc='#D9D9D9' if core else 'white',ec='#BDBDBD',lw=0.4))
        ax.text(x0+(i+.5)*cw,yref,ch,ha='center',va='center',fontsize=8,fontweight='bold')
    heads=['Matched\n(of 12)','Longest\nrun (nt)','Core\nmismatch','Junction\nreads']
    xc=[x0+12*cw+3+i*14.2 for i in range(4)]
    for x,h in zip(xc,heads): ax.text(x+5.5,3.6,h,ha='center',va='center',fontsize=FS,linespacing=1.05)
    ax.plot([0,112],[9.6,9.6],color='black',lw=0.6)
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
        cm='to be\nconfirmed' if und else str(int(row.core_mismatch_vs_CTAAAC))
        vals=[str(int(row.matched_of_12)),str(int(row.longest_contiguous)),cm,fmt_reads(row.canonical_reads)]
        for x,v in zip(xc,vals): ax.text(x+5.5,y,v,ha='center',va='center',fontsize=FS+(0 if und else 0.5),linespacing=0.95)
    note(c,0.5,50.3,'Body windows aligned at the TRS-B core start (S GTAAAC differs from CTAAAC at its first base).')
    note(c,0.5,53.3,'ORF3: TRS-B core to be confirmed — window placed at breakpoint 24718 (hatched).')
    return A_T[['gene','alignment_anchor','anchor_type','body_12mer','ref_12mer','matched_of_12','longest_contiguous','core_hexamer','core_mismatch_vs_CTAAAC','canonical_reads','canonical_share']]
# ---------------------------------------------------------------- B
def draw_B(c):
    c.letter('B'); c.title('Pairing metrics vs canonical abundance',x=7)
    tx=c.ax(0,8,68,26); tx.axis('off'); tx.set_xlim(0,68); tx.set_ylim(26,0)
    tx.text(0.5,1.4,'Metric (canonical, n = 5)',fontsize=FS,fontweight='bold',va='center'); tx.text(47,1.4,'ρ',fontsize=FS,fontweight='bold',ha='center',va='center'); tx.text(59,1.4,'exact P',fontsize=FS,fontweight='bold',ha='center',va='center')
    tx.plot([0,68],[3.0,3.0],color='black',lw=0.5)
    labs={'matched_of_12':'Matched bases (of 12)†','longest_contiguous':'Longest contiguous run†','mfe_a_kcal':'MFE (a), nascent-strand','mfe_b_kcal':'MFE (b), direct pairing'}
    for i,(col,nm) in enumerate(labs.items()):
        r=B3[B3.column==col].iloc[0]; y=5.3+i*2.9
        tx.text(0.5,y,nm,fontsize=FS,va='center'); tx.text(47,y,'%+.2f'%r.rho,fontsize=FS,ha='center',va='center'); tx.text(59,y,'%.2f'%r.P_exact_two_sided,fontsize=FS,ha='center',va='center')
    tx.text(0.5,17.6,'† tied values (average ranks). Exact P over all 120',fontsize=FS,va='center',color=GREY)
    tx.text(0.5,20.5,'permutations; smallest attainable two-sided P = 0.017.',fontsize=FS,va='center',color=GREY)
    tx.text(0.5,23.4,'None supported; descriptive only (n = 5).',fontsize=FS,va='center',color=GREY,fontweight='bold')
    ax=c.ax(13,38,50,13)
    for g in GENES:
        r=A_T[A_T.gene==g].iloc[0]; p=PSR.loc[g]
        ax.errorbar(r.mfe_a_kcal,p.mean_reads_per_sample,yerr=[[p.mean_reads_per_sample-p.ci_lo],[p.ci_hi-p.mean_reads_per_sample]],fmt='o',ms=4.5,mfc=GENE_COL[g],mec='black',mew=0.5,ecolor='black',elinewidth=0.7,capsize=1.5,zorder=3)
        ax.annotate(g,(r.mfe_a_kcal,p.mean_reads_per_sample),xytext=(4,3),textcoords='offset points',fontsize=FS)
    ax.set_yscale('log'); ax.set_ylim(2e2,3e5); ax.set_xlim(-21,-10)
    ax.set_xlabel('MFE (a) (kcal/mol)',fontsize=7.5,labelpad=1); ax.set_ylabel('Reads / sample',fontsize=7.5,labelpad=1)
    ax.grid(color='#EEEEEE',lw=0.4)
    return B3
# ---------------------------------------------------------------- C
def draw_C(c):
    c.letter('C'); c.title('Acceptor sequence context',x=7)
    labs=[('canonical','Canonical (n = 5)'),('TRSL_dep_noncanonical','TRS-L-dep noncanonical (n = 906)'),('TRSL_indep_distant','TRS-L-indep distant (n = 1,688)'),
          ('TRSL_indep_local','TRS-L-indep local (n = 1,479)'),('background_region_matched_positions','Region-matched random positions (n = {:,})')]
    nn=CF.groupby('group').n.first()
    x=42
    for nt in 'ACGU':
        c.fig.add_artist(Rectangle(((c.ox+x)/c.W,1-(c.oy+5.9)/c.H),2.4/c.W,2.4/c.H,fc=NTC[nt],ec='#777777',lw=0.3,transform=c.fig.transFigure))
        c.text(x+3.4,4.7,nt,va='center'); x+=8
    axes=[]
    for i,(g,lab) in enumerate(labs):
        y=9.0+i*10.0
        ax=c.ax(13,y+3.4,62,6.0); sub=CF[CF.group==g]
        bottom=np.zeros(41)
        for nt in 'ACGU':
            s=sub[sub.nucleotide==nt].sort_values('position')
            ax.bar(s.position,s.frequency,bottom=bottom,width=1.0,color=NTC[nt],ec='none',lw=0)
            bottom+=s.frequency.values
        ax.axvline(-0.5,color='black',lw=0.8,ls='--'); ax.set_xlim(-25.5,15.5); ax.set_ylim(0,1); ax.set_yticks([0,.5,1]); ax.set_yticklabels(['0','.5','1'])
        if i<4: ax.set_xticklabels([])
        c.text(13,y+1.6,lab.format(int(nn[g])) if '{' in lab else lab,va='center',ha='left',fontsize=FS)
        axes.append(ax)
    axes[-1].set_xlabel('Position relative to acceptor (nt); dashed = junction',fontsize=7.5,labelpad=1)
    c.text(2.4,33,'Frequency (stacked)',rotation=90,ha='center',va='center',fontsize=7.5)
    return CF
# ---------------------------------------------------------------- E
def draw_E(c):
    c.letter('E'); c.title('Donor–acceptor microhomology',x=7)
    c.text(12,8.6,'E1  by class',fontsize=7.5,fontweight='bold',va='center'); c.text(65,8.6,'E2  by recurrence',fontsize=7.5,fontweight='bold',va='center')
    ax=c.ax(12,21,48,26)
    order=[('canonical','Can.'),('TRSL_dep_noncanonical','Dep.'),('TRSL_indep_distant','Dist.'),('TRSL_indep_local','Loc.')]
    mm=MHS3[MHS3.subset=='filtered (main)'].set_index('group'); rng=np.random.default_rng(3)
    for i,(g,lab) in enumerate(order):
        x=i
        if g=='canonical':
            d=M[M.canonical_primary]
            for _,r in d.iterrows(): ax.scatter(x-0.18+rng.uniform(-.08,.08),r.mh_len,s=16,color=GENE_COL[r.canonical_gene],edgecolor='black',lw=0.4,zorder=5)
            ref=MHC[(MHC.group=='canonical_primary5')&MHC.background.str.startswith('uniform')].iloc[0]
            ax.plot([x+0.05,x+0.4],[ref.null_mean]*2,color=COL['bg'],lw=1.4); ax.plot([x+0.22]*2,[ref.null_p2_5,ref.null_p97_5],color=COL['bg'],lw=1.0)
            nt_='%s\n%+.1f nt\n(%.1f×)'%(star(ref.P),d.mh_len.mean()-ref.null_mean,d.mh_len.mean()/ref.null_mean)
        else:
            key=GKEY[g]; dobs=MH_E1[(MH_E1.group==g)&(MH_E1.kind=='observed')].mh_len
            dbg=pd.read_csv(RES+'mh_E1_spanshuf_%s.tsv'%g,sep='\t').mh_len
            violin(ax,dobs,x-0.19,COL[key],width=0.36); violin(ax,dbg,x+0.19,COL['bg'],width=0.36)
            r=mm.loc[g]; nt_='%s\n%+.2f nt\n(%.1f×)'%(star(r.P_mean_holm_across_3_classes),r.delta_nt,r.fold_mean)
            if g=='TRSL_indep_distant':
                r2=MHS3[MHS3.subset=='excl. donor<=100'].iloc[0]
                ax.scatter([x-0.19],[r2.obs_mean_mh],marker='D',s=14,facecolor='white',edgecolor='black',lw=0.7,zorder=6)
                ax.scatter([x+0.19],[r2.null_mean_mh],marker='D',s=14,facecolor='white',edgecolor='black',lw=0.7,zorder=6)
        ax.text(x,1.02,nt_,transform=ax.get_xaxis_transform(),ha='center',va='bottom',fontsize=FS,linespacing=1.0)
    ax.set_xticks(range(4)); ax.set_xticklabels([l for _,l in order]); ax.set_xlim(-.6,3.6); ax.set_ylim(0,14)
    ax.set_ylabel('Microhomology (nt)',fontsize=7.5,labelpad=1)
    note(c,0.5,53.6,'Gray: shuffled within class and span quartile.')
    note(c,0.5,56.6,'Open diamond: distant excl. donor ≤ 100.')
    a2=c.ax(67,21,21,26)
    bins=['1','2','3-4','>=5']; xs=np.arange(4)
    for g,key in [('TRSL_dep_noncanonical','dep'),('TRSL_indep_distant','distant'),('TRSL_indep_local','local')]:
        s=MHS[MHS.group==g].set_index('samples_bin').loc[bins]; off={'dep':-0.09,'distant':0,'local':0.09}[key]
        a2.errorbar(xs+off,s.delta_nt,yerr=[s.delta_nt-s.delta_ci_lo,s.delta_ci_hi-s.delta_nt],color=COL[key],marker='o',ms=3.3,lw=1.0,capsize=1.5,elinewidth=0.8)
    a2.axhline(0,color='#555555',lw=0.7,ls='--'); a2.set_xticks(xs); a2.set_xticklabels(['1','2','3–4','≥5']); a2.set_xlim(-.4,3.4)
    a2.set_xlabel('Positive samples',fontsize=7.5,labelpad=1); a2.set_ylabel('Δ vs stratum shuffle (nt)',fontsize=7.5,labelpad=1); a2.set_ylim(-0.2,2.2)
    a2.text(0.5,1.02,'all events\n(unfiltered)',transform=a2.transAxes,ha='center',va='bottom',fontsize=FS,color=GREY,linespacing=1.0)
    return pd.concat([MH_E1.assign(panel='E1_observed_and_unstratified_draw'),MHS.assign(panel='E2')],ignore_index=True)
# ---------------------------------------------------------------- F
def draw_F(c):
    c.letter('F'); c.title('Direct plus-strand pairing (MFE b)',x=7)
    c.text(13,8.6,'F1  vs shuffled background',fontsize=7.5,fontweight='bold',va='center'); c.text(63,8.6,'F2  vs reads',fontsize=7.5,fontweight='bold',va='center')
    ax=c.ax(13,21,44,26)
    order=[('canonical','Can.'),('TRSL_dep_noncanonical','Dep.'),('TRSL_indep_distant','Dist.'),('TRSL_indep_local','Loc.')]
    mb=MFS3[(MFS3['mode']=='b')&(MFS3.subset=='filtered (main)')].set_index('group'); rng=np.random.default_rng(4); nclip=0
    for i,(g,lab) in enumerate(order):
        x=i
        if g=='canonical':
            d=M[M.canonical_primary]
            for _,r in d.iterrows(): ax.scatter(x-0.18+rng.uniform(-.08,.08),r.mfe_b,s=16,color=GENE_COL[r.canonical_gene],edgecolor='black',lw=0.4,zorder=5)
            ref=SEP_MFE[(SEP_MFE.group=='canonical_primary5')&(SEP_MFE['mode']=='b')&SEP_MFE.null_type.str.startswith('acceptor')].iloc[0]
            ax.plot([x+0.05,x+0.4],[ref.null_mean]*2,color=COL['bg'],lw=1.4); ax.plot([x+0.22]*2,[ref.null_p2_5,ref.null_p97_5],color=COL['bg'],lw=1.0)
            nt_='ns\n%+.1f'%(d.mfe_b.mean()-ref.null_mean)
        else:
            key=GKEY[g]; dobs=MF1[(MF1.group==g)&(MF1['mode']=='b')&(MF1.kind=='observed')].mfe.values
            dbg=MFP3[(MFP3.group==g)&(MFP3['mode']=='b')].mfe.values
            nclip+=int((dobs<-30).sum()+(dbg<-30).sum())
            violin(ax,np.clip(dobs,-30,None),x-0.19,COL[key],width=0.36); violin(ax,np.clip(dbg,-30,None),x+0.19,COL['bg'],width=0.36)
            r=mb.loc[g]; nt_='%s\n%+.2f'%(star(r.P_holm_across_3_classes),r.delta_kcal)
            if g=='TRSL_indep_distant':
                r2=MFS3[(MFS3['mode']=='b')&(MFS3.subset=='excl. donor<=100')].iloc[0]
                ax.scatter([x-0.19],[r2.obs_mean_mfe],marker='D',s=14,facecolor='white',edgecolor='black',lw=0.7,zorder=6)
                ax.scatter([x+0.19],[r2.null_mean],marker='D',s=14,facecolor='white',edgecolor='black',lw=0.7,zorder=6)
        ax.text(x,1.02,nt_,transform=ax.get_xaxis_transform(),ha='center',va='bottom',fontsize=FS,linespacing=1.0)
    ax.set_xticks(range(4)); ax.set_xticklabels([l for _,l in order]); ax.set_xlim(-.6,3.6); ax.set_ylim(-30,0)
    ax.set_ylabel('MFE (b) (kcal/mol)',fontsize=7.5,labelpad=1)
    note(c,0.5,53.6,'Δ = observed − shuffled (kcal/mol).')
    note(c,0.5,56.6,'Local: donor/acceptor are neighbours (confound).')
    a2=c.ax(66,21,22,26)
    f2=MF2[MF2['mode']=='b']; lines=[]
    for g,key,lab in [('TRSL_dep_noncanonical','dep','Dep.'),('TRSL_indep_distant','distant','Dist.'),('TRSL_indep_local','local','Loc.')]:
        d=M[(M.cls==g)&M.pass_filter]
        a2.scatter(d.mfe_b,d.reads,s=2.5,color=COL[key],alpha=0.35,lw=0,rasterized=True)
        fl=d[d.flag_leaderdonor_at_canon_site_in_distant]
        if len(fl): a2.scatter(fl.mfe_b,fl.reads,s=12,facecolor='none',edgecolor='black',lw=0.6,marker='s',zorder=5)
        r=f2[(f2.group==g)&(f2.subset=='filtered (main)')].iloc[0]; lines.append((lab,COL[key],r.rho))
    cp=M[M.canonical_primary]; a2.scatter(cp.mfe_b,cp.reads,s=18,marker='D',color=COL['canonical'],edgecolor='black',lw=0.4,zorder=6)
    rc=f2[f2.group=='canonical_primary5'].iloc[0]; lines.insert(0,('Can.',COL['canonical'],rc.rho))
    a2.set_yscale('log'); a2.set_ylim(0.7,4e6); a2.set_xlim(-30,0)
    a2.set_xlabel('MFE (b) (kcal/mol)',fontsize=7.5,labelpad=1); a2.set_ylabel('Junction reads',fontsize=7.5,labelpad=1)
    for k,(lab,col,rho) in enumerate(lines): c.text(67,11.6+k*2.75,'%s ρ %+.2f'%(lab,rho),color=col,fontsize=FS,va='center',ha='left',fontweight='bold')
    return pd.DataFrame(dict(group=[l[0] for l in lines],rho_within_class=[l[2] for l in lines]))
# ---------------------------------------------------------------- G
from s4_g import draw_G
