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
from s4_panels import *
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties
from matplotlib.transforms import Affine2D
import matplotlib.colors as mcolors
# ---- final palette: three colours (blue canonical, orange noncanonical in three tints, gray background)
C5={'canonical':'#2171B5','dep':'#B35806','distant':'#E69F00','local':'#F6D5A0','bg':'#8C8C8C'}
ORD5=[('canonical','Can.'),('TRSL_dep_noncanonical','Dep.'),('TRSL_indep_distant','Dist.'),('TRSL_indep_local','Loc.'),('bg','Bkg.')]
GK={'canonical':'canonical','TRSL_dep_noncanonical':'dep','TRSL_indep_distant':'distant','TRSL_indep_local':'local','bg':'bg'}
SIMM=pd.read_csv(RES+'sim_merged_effect.tsv',sep='\t').iloc[0]; SIMP=pd.read_csv(RES+'sim_merged_points.tsv',sep='\t')
NN={'canonical':5,'TRSL_dep_noncanonical':906,'TRSL_indep_distant':1688,'TRSL_indep_local':1479}
def title2(c,s,x=7):
    c.text(x,1.0,s,fontsize=9,fontweight='bold',ha='left',va='top',linespacing=1.0)
def legend5(c,y=1.0,x0=1.0,maxw=1e9):
    items=[('canonical','Canonical'),('dep','TRS-L-dependent'),('distant','TRS-L-independent distant'),('local','TRS-L-independent local'),('bg','Shuffled / random background')]
    x=x0
    for k,lab in items:
        if x+3.6+len(lab)*1.3>maxw and x>x0: x=x0; y+=4.2
        c.fig.add_artist(Rectangle(((c.ox+x)/c.W,1-(c.oy+y+3.2)/c.H),2.6/c.W,2.6/c.H,fc=C5[k],ec='#777777' if k=='local' else 'none',lw=0.4,transform=c.fig.transFigure))
        c.text(x+3.6,y+1.4,lab,va='center',ha='left'); x+=3.6+len(lab)*1.3+4.0
# ---------------------------------------------------------------- A
def draw_A(c):
    c.letter('A'); title2(c,'Canonical TRS-L/TRS-B\nsequence architecture')
    ax=c.ax(0,11,90,46); ax.set_xlim(0,90); ax.set_ylim(46,0); ax.axis('off')
    cw=4.3; x0=14; ref=A_T.ref_12mer[0]; ys=[5.0]+[12.5+k*6.6 for k in range(5)]
    def cell(x,y,ch,fc,tc,bold=False,base=None):
        if base: ax.add_patch(Rectangle((x,y-2.7),cw,5.4,fc=base,ec='white',lw=0.4))
        ax.add_patch(Rectangle((x+(0.45 if base else 0),y-2.7+(0.45 if base else 0)),cw-(0.9 if base else 0),5.4-(0.9 if base else 0),fc=fc,ec='#D0D0D0' if not base else 'none',lw=0.4))
        ax.text(x+cw/2,y,ch,ha='center',va='center',fontsize=8.5,color=tc,fontweight='bold' if bold else 'normal')
    ax.text(0.5,ys[0],'TRS-L',fontsize=8,fontweight='bold',va='center')
    for i,ch in enumerate(ref): cell(x0+i*cw,ys[0],ch,'white','black',True,'#BDBDBD' if 3<=i<=8 else None)
    for r,(_,row) in enumerate(A_T.iterrows()):
        y=ys[r+1]; g=row.gene; und=g=='ORF3'
        ax.text(0.5,y,g,fontsize=8,fontweight='bold',va='center')
        for k,ch in enumerate(row.body_12mer):
            mt=ch==ref[k]; core=3<=k<=8 and not und
            cell(x0+k*cw,y,ch,C5['canonical'] if mt else 'white','white' if mt else '#444444',mt,'#BDBDBD' if core else None)
        if und: ax.text(x0+12*cw+1.5,y,'core undetermined',fontsize=7,va='center',ha='left',color='#555555',style='italic')
    return A_T[['gene','alignment_anchor','anchor_type','body_12mer','ref_12mer','matched_of_12','longest_contiguous','core_hexamer','core_mismatch_vs_CTAAAC','canonical_reads','canonical_share']]
# ---------------------------------------------------------------- B (logos + best-shift)
def logo_freqs(g):
    sub=CF[CF.group==g]; n=int(sub.n.iloc[0]); F=np.zeros((41,4))
    for j,nt in enumerate('ACGU'): F[:,j]=sub[sub.nucleotide==nt].sort_values('position').frequency.values
    return F,n
def merged_freqs():
    Fs=[];ns=[]
    for g in ['TRSL_dep_noncanonical','TRSL_indep_distant','TRSL_indep_local']:
        F,n=logo_freqs(g); Fs.append(F*n); ns.append(n)
    return sum(Fs)/sum(ns),sum(ns)
def info_content(F,n):
    with np.errstate(divide='ignore',invalid='ignore'):
        H=-np.nansum(np.where(F>0,F*np.log2(F),0),axis=1)
    en=3/(2*np.log(2)*n)
    return np.clip(2-(H+en),0,None)
def draw_logo(ax,F,n,color):
    ic=info_content(F,n); fp=FontProperties(family='DejaVu Sans',weight='bold')
    for i in range(41):
        pos=-25+i; base=0.0
        for j in np.argsort(F[i]):
            h=F[i,j]*ic[i]
            if h<1e-4: continue
            tp=TextPath((0,0),'ACGU'[j],size=1,prop=fp); bb=tp.get_extents()
            tr=Affine2D().translate(-bb.x0,-bb.y0).scale(0.92/bb.width,h/bb.height).translate(pos-0.46,base)
            ax.add_patch(PathPatch(tr.transform_path(tp),fc=color,ec='none',lw=0)); base+=h
    return ic
def draw_B(c):
    c.letter('B'); title2(c,'Canonical TRS similarity is lost\nin noncanonical junctions')
    Fm,nm=merged_freqs(); rows=[('canonical',logo_freqs('canonical'),'Canonical (n = 5)','canonical'),('nc',(Fm,nm),'Noncanonical (n = {:,})'.format(nm),'dep'),
        ('bg',logo_freqs('background_region_matched_positions'),'Random positions (n = {:,})','bg')]
    ics=[]
    for i,(k,(F,n),lab,col) in enumerate(rows):
        y=11+i*13.2; ax=c.ax(12,y+2.8,38,8.2)
        ic=draw_logo(ax,F,n,C5[col]); ics.append(pd.DataFrame(dict(group=k,n=n,position=range(-25,16),information_bits=ic)))
        ax.set_xlim(-25.5,15.5); ax.set_ylim(0,2); ax.set_yticks([0,1,2]); ax.axvline(-0.5,color='black',lw=0.7,ls='--')
        if i<2: ax.set_xticklabels([])
        c.text(12,y+1.2,(lab.format(n) if '{' in lab else lab),va='center',ha='left',fontsize=FS)
    ax.set_xticks([-20,-10,0,10]); ax.set_xlabel('Position relative to acceptor (nt)',fontsize=7.5,labelpad=1)
    c.text(2.2,33,'Information content (bits)',rotation=90,ha='center',va='center',fontsize=7.5)
    a2=c.ax(67,16,21,32)
    can=D2[D2.group=='canonical']; rng=np.random.default_rng(1)
    a2.scatter(rng.uniform(-.1,.1,len(can)),can.best_shift_matches,s=16,color=C5['canonical'],edgecolor='black',lw=0.4,zorder=5); a2.plot([-.3,.3],[can.best_shift_matches.mean()]*2,color='black',lw=0.9)
    violin(a2,SIMP.best_shift_matches,1,C5['dep'],width=0.8); violin(a2,SIMP.background_one_shuffle,2,C5['bg'],width=0.8)
    a2.set_xticks([0,1,2]); a2.set_xticklabels(['Can.\n5','Non.\n{:,}'.format(len(SIMP)),'Bkg.\n{:,}'.format(len(SIMP))]); a2.set_xlim(-.6,2.6); a2.set_ylim(0,12.5); a2.set_yticks([0,3,6,9,12])
    a2.set_ylabel('Best-shift TRS-L similarity\n(matched bases of 12)',fontsize=7.5,labelpad=1)
    e=pd.read_csv(RES+'sim_effect_matches.tsv',sep='\t'); ec=e[(e.group=='canonical')&e.background.str.startswith('comp')].iloc[0]
    a2.text(-0.15,1.02,'Δ %+.1f'%ec.delta_matches,transform=a2.get_xaxis_transform(),ha='center',va='bottom',fontsize=FS,color=C5['canonical'],fontweight='bold')
    a2.text(1.25,1.02,'Δ %+.2f'%SIMM.delta_matches,transform=a2.get_xaxis_transform(),ha='center',va='bottom',fontsize=FS,color=C5['dep'],fontweight='bold')
    return pd.concat(ics,ignore_index=True)
# ---------------------------------------------------------------- C
def _five(ax,get,ylim,ylabel,annot,dobs_fun=None):
    pass
def draw_C(c):
    c.letter('C'); title2(c,'Donor–acceptor microhomology\nat junction sites')
    ax=c.ax(13,24,74,30); mm=MHS3[MHS3.subset=='filtered (main)'].set_index('group'); rng=np.random.default_rng(3); bgs=[]
    for i,(g,lab) in enumerate(ORD5):
        if g=='canonical':
            d=M[M.canonical_primary]; ax.scatter(i+rng.uniform(-.12,.12,len(d)),d.mh_len,s=16,color=C5['canonical'],edgecolor='black',lw=0.4,zorder=5)
            ref=MHC[(MHC.group=='canonical_primary5')&MHC.background.str.startswith('uniform')].iloc[0]; ax.plot([i+0.25,i+0.45],[ref.null_mean]*2,color=C5['bg'],lw=1.6)
            nt_='%s\n%+.1f nt\n(%.1f×)'%(star(ref.P),d.mh_len.mean()-ref.null_mean,d.mh_len.mean()/ref.null_mean)
        elif g=='bg':
            bg=np.concatenate([pd.read_csv(RES+'mh_E1_spanshuf_%s.tsv'%q,sep='\t').mh_len.values for q in CLS3]); violin(ax,bg,i,C5['bg'],width=0.8); nt_=''
        else:
            dobs=MH_E1[(MH_E1.group==g)&(MH_E1.kind=='observed')].mh_len; violin(ax,dobs,i,C5[GK[g]],width=0.8)
            r=mm.loc[g]; nt_='%s\n%+.2f nt\n(%.1f×)'%(star(r.P_mean_holm_across_3_classes),r.delta_nt,r.fold_mean)
            if g=='TRSL_indep_distant':
                r2=MHS3[MHS3.subset=='excl. donor<=100'].iloc[0]; ax.scatter([i],[r2.obs_mean_mh],marker='D',s=14,facecolor='white',edgecolor='black',lw=0.7,zorder=6)
        ax.text(i,1.02,nt_,transform=ax.get_xaxis_transform(),ha='center',va='bottom',fontsize=FS,linespacing=1.0)
    ax.set_xticks(range(5)); ax.set_xticklabels(['Can.\n5','Dep.\n906','Dist.\n1,688','Loc.\n1,479','Bkg.']); ax.set_xlim(-.6,4.6); ax.set_ylim(0,14)
    ax.set_ylabel('Microhomology (nt)',fontsize=7.5,labelpad=1)
    return MH_E1
# ---------------------------------------------------------------- D
def draw_D(c):
    c.letter('D'); title2(c,'Predicted RNA pairing does not explain\njunction abundance')
    ax=c.ax(13,24,40,30); mb=MFS3[(MFS3['mode']=='b')&(MFS3.subset=='filtered (main)')].set_index('group'); rng=np.random.default_rng(4); nclip=0
    for i,(g,lab) in enumerate(ORD5):
        if g=='canonical':
            d=M[M.canonical_primary]; ax.scatter(i+rng.uniform(-.12,.12,len(d)),d.mfe_b,s=16,color=C5['canonical'],edgecolor='black',lw=0.4,zorder=5)
            ref=SEP_MFE[(SEP_MFE.group=='canonical_primary5')&(SEP_MFE['mode']=='b')&SEP_MFE.null_type.str.startswith('acceptor')].iloc[0]; ax.plot([i+0.25,i+0.45],[ref.null_mean]*2,color=C5['bg'],lw=1.6)
            nt_='ns\n%+.1f'%(d.mfe_b.mean()-ref.null_mean)
        elif g=='bg':
            bg=np.concatenate([MFP3[(MFP3.group==q)&(MFP3['mode']=='b')].mfe.values for q in CLS3]); nclip+=int((bg<-30).sum()); violin(ax,np.clip(bg,-30,None),i,C5['bg'],width=0.8); nt_=''
        else:
            dobs=MF1[(MF1.group==g)&(MF1['mode']=='b')&(MF1.kind=='observed')].mfe.values; nclip+=int((dobs<-30).sum()); violin(ax,np.clip(dobs,-30,None),i,C5[GK[g]],width=0.8)
            r=mb.loc[g]; nt_='%s\n%+.2f'%(star(r.P_holm_across_3_classes),r.delta_kcal)
        ax.text(i,1.02,nt_,transform=ax.get_xaxis_transform(),ha='center',va='bottom',fontsize=FS,linespacing=1.0)
    ax.set_xticks(range(5)); ax.set_xticklabels(['Can.','Dep.','Dist.','Loc.','Bkg.']); ax.set_xlim(-.6,4.6); ax.set_ylim(-30,0)
    ax.set_ylabel('MFE (b) (kcal/mol)',fontsize=7.5,labelpad=1)
    a2=c.ax(63,24,25,30); f2=MF2[MF2['mode']=='b']; lines=[]
    for g,key,lab in [('TRSL_dep_noncanonical','dep','Dep.'),('TRSL_indep_distant','distant','Dist.'),('TRSL_indep_local','local','Loc.')]:
        d=M[(M.cls==g)&M.pass_filter]; a2.scatter(d.mfe_b,d.reads,s=2.5,color=C5[key],alpha=0.4 if key!='local' else 0.7,lw=0,rasterized=True)
        r=f2[(f2.group==g)&(f2.subset=='filtered (main)')].iloc[0]; lines.append((lab,C5[key],r.rho,len(d)))
    cp=M[M.canonical_primary]; a2.scatter(cp.mfe_b,cp.reads,s=18,marker='D',color=C5['canonical'],edgecolor='black',lw=0.4,zorder=6)
    rc=f2[f2.group=='canonical_primary5'].iloc[0]; lines.insert(0,('Can.',C5['canonical'],rc.rho,5))
    a2.set_yscale('log'); a2.set_ylim(0.7,4e6); a2.set_xlim(-30,0); a2.set_xlabel('MFE (b) (kcal/mol)',fontsize=7.5,labelpad=1); a2.set_ylabel('Junction reads',fontsize=7.5,labelpad=1)
    for k,(lab,col,rho,n) in enumerate(lines):
        c.text(63,12.6+k*2.9,'%s n=%s ρ %+.2f'%(lab,format(n,','),rho),color='#222222' if lab=='Loc.' else col,fontsize=FS,va='center',ha='left',fontweight='bold')
    return pd.DataFrame(dict(group=[l[0] for l in lines],n=[l[3] for l in lines],rho_within_class=[l[2] for l in lines]))
# ---------------------------------------------------------------- E
def junction(ax,x,y,color,marks,mcol):
    # observed junction structure: leader (left) and body (right) with skipped sequence (dashed) and joined product (below)
    ax.add_patch(Rectangle((x,y),8,2.4,fc='#B8B8B8',ec='none')); ax.plot([x+8,x+30],[y+1.2]*2,color='#777777',lw=0.7,ls=':'); ax.add_patch(Rectangle((x+30,y),20,2.4,fc='#B8B8B8',ec='none'))
    ax.add_patch(Rectangle((x,y+6.2),8,2.4,fc='#B8B8B8',ec='none')); ax.add_patch(Rectangle((x+8,y+6.2),20,2.4,fc='#B8B8B8',ec='none'))
    ax.plot([x+8,x+8],[y+3.4,y+5.2],color=color,lw=0.8); ax.plot([x+30,x+30],[y+3.4,y+5.2],color=color,lw=0.8)
    for k in range(marks):
        ax.add_patch(Rectangle((x+5.5-1.1*k,y),1.0,2.4,fc=mcol,ec='none')); ax.add_patch(Rectangle((x+31.0+1.1*k,y),1.0,2.4,fc=mcol,ec='none'))
def draw_E(c):
    c.letter('E'); c.title('A continuum of sequence constraint in PEDV template switching',x=7)
    ax=c.ax(0,8,180,28); ax.set_xlim(0,180); ax.set_ylim(28,0); ax.axis('off')
    zones=[(0.5,'Canonical: TRS-guided',C5['canonical']),(60.5,'Noncanonical',C5['dep']),(120.5,'TRS-poor / unresolved','#555555')]
    for x,ttl,col in zones:
        ax.add_patch(Rectangle((x,0.5),59,17.5,fc='none',ec=col,lw=0.9,ls='-' if col!='#555555' else (0,(4,2))))
        ax.text(x+29.5,3.0,ttl,fontsize=8,fontweight='bold',ha='center',va='center',color=col)
    junction(ax,5,7.0,C5['canonical'],6,C5['canonical']); junction(ax,65,7.0,C5['dep'],2,C5['dep']); junction(ax,125,7.0,'#555555',0,'#8C8C8C')
    grad=np.linspace(0,1,256)[None,:]; cmap=mcolors.LinearSegmentedColormap.from_list('g',[C5['canonical'],'#D9D9D9'])
    ax.imshow(grad,extent=[1,172,26.8,23.6],cmap=cmap,aspect='auto',zorder=1)
    ax.add_patch(Polygon([[172,22.4],[178.5,25.2],[172,28.0]],closed=True,fc='#D9D9D9',ec='none'))
    ax.text(1.5,20.6,'Strong',fontsize=7.5,fontweight='bold',color=C5['canonical'],va='center',ha='left'); ax.text(171,20.6,'Weak',fontsize=7.5,fontweight='bold',color='#555555',va='center',ha='right')
    ax.text(90,20.6,'Strength of canonical sequence constraint',fontsize=7.5,fontweight='bold',va='center',ha='center')
    ax.set_xlim(0,180); ax.set_ylim(28,0)
    return pd.DataFrame(dict(note=['schematic of observed junction structure; no data']))
