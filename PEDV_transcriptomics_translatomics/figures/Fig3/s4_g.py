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
from s3_panels import *
import matplotlib.colors as mcolors
def draw_G(c):
    c.letter('G'); c.title('Conservative working model (what A–F support)',x=7)
    ax=c.ax(0,7,180,37); ax.set_xlim(0,180); ax.set_ylim(37,0); ax.axis('off')
    zones=[(0.5,'#E8F0F8','Canonical: TRS-guided',COL['canonical']),(60.5,'#F3F3F3','Noncanonical junctions','#7A5A08'),(120.5,'#FFFFFF','TRS-poor / unresolved','#555555')]
    for x,fc,ttl,col in zones:
        rbox(ax,x,0.4,59,26.0,fc,ec=col if col!='#7A5A08' else COL['dep'],ls='-' if fc!='#FFFFFF' else (0,(4,2)))
        ax.text(x+29.5,2.9,ttl,fontsize=7.5,fontweight='bold',ha='center',va='center',color=col)
    ax.add_patch(Rectangle((4,11.3),4,2.4,fc='#B8860B',ec='none')); ax.plot([8,56],[12.5]*2,color='#777777',lw=0.7)
    for cx,wd,fc in [(30,6,'#E1EBF5'),(40,5,'#EAE3F3'),(48,7,'#E0F0F0')]: ax.add_patch(Rectangle((cx,11.3),wd,2.4,fc=fc,ec='#9A9A9A',lw=0.4))
    ax.annotate('',xy=(32,11.0),xytext=(6,11.0),arrowprops=dict(arrowstyle='-|>',color=COL['canonical'],lw=1.0,connectionstyle='arc3,rad=-0.25'))
    ax.text(30.5,19.2,'TRS-L / TRS-B complementarity:\nmatched 8.8 vs 6.6 bases (Δ +2.2),\ncore mismatch 0–1',fontsize=FS,ha='center',va='center',color='#333333',linespacing=1.05)
    ax.text(30.5,24.4,'ORF3 TRS-B to be confirmed',fontsize=FS,ha='center',va='center',color=GREY)
    ax.text(90,7.0,'TRS-L similarity lost (Δ ≤ +0.19)',fontsize=FS,ha='center',va='center',color='#333333')
    ax.text(90,10.4,'microhomology elevated (+0.8 to +1.1 nt)',fontsize=FS,ha='center',va='center',color='#333333')
    ax.text(90,13.8,'origin undetermined:',fontsize=FS,ha='center',va='center',fontweight='bold')
    rbox(ax,62,16.0,26,9.4,'white',ec=COL['dep']); ax.text(75,18.3,'Biological',fontsize=FS,ha='center',va='center',fontweight='bold',color='#7A5A08'); ax.text(75,22.2,'template switching\nat short homology',fontsize=FS,ha='center',va='center',linespacing=0.95)
    rbox(ax,92,16.0,26,9.4,'white',ec='#555555'); ax.text(105,18.3,'Library prep.',fontsize=FS,ha='center',va='center',fontweight='bold',color='#555555'); ax.text(105,22.2,'RT/PCR chimera at\nshort homology',fontsize=FS,ha='center',va='center',linespacing=0.95)
    ax.text(90,20.7,'or',fontsize=7.5,ha='center',va='center',fontweight='bold')
    ax.text(150,7.6,'no further feature resolved',fontsize=FS,ha='center',va='center',color='#333333')
    ax.text(150,12.2,'direct pairing MFE (b): background,\nexcept local (proximity-confounded)',fontsize=FS,ha='center',va='center',color='#333333',linespacing=1.05)
    ax.text(150,18.4,'RNA structure / proximity ?',fontsize=FS,ha='center',va='center',color=GREY,style='italic'); ax.text(150,22.0,'(candidate, untested)',fontsize=FS,ha='center',va='center',color=GREY)
    grad=np.linspace(0,1,256)[None,:]
    cmap=mcolors.LinearSegmentedColormap.from_list('g',[COL['canonical'],'#BDBDBD'])
    ax.imshow(grad,extent=[2,172,36.0,33.0],cmap=cmap,aspect='auto',zorder=1)
    ax.add_patch(Polygon([[172,31.6],[178.5,34.5],[172,37.4]],closed=True,fc='#BDBDBD',ec='none'))
    ax.text(2,29.8,'Strength of canonical sequence constraint:',fontsize=FS,va='center',ha='left',fontweight='bold')
    ax.text(84,29.8,'Strong',fontsize=FS,va='center',ha='center',color=COL['canonical'],fontweight='bold')
    ax.text(172,29.8,'Weak',fontsize=FS,va='center',ha='right',color='#555555',fontweight='bold')
    ax.set_xlim(0,180); ax.set_ylim(37,0)
    return pd.DataFrame(dict(note=['schematic; statements restricted to Supported = yes in the decision table']))
