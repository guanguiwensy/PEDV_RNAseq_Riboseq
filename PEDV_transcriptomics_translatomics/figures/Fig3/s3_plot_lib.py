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
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, FancyArrowPatch, Polygon
from matplotlib.lines import Line2D
from common import W, B
RES=W+'res/'
plt.rcParams.update({'font.family':'Arial','pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':0.6,'xtick.major.width':0.6,'ytick.major.width':0.6,
    'xtick.major.size':2.5,'ytick.major.size':2.5,'font.size':7,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,
    'axes.spines.top':False,'axes.spines.right':False,'mathtext.default':'regular'})
FS=7
# ---------------- project palette (group -> hex -> meaning), consistent with Fig2
COL={'canonical':'#2171B5','dep':'#B8860B','distant':'#4DAF4A','local':'#984EA3','bg':'#8C8C8C'}
GENE_COL={'S':'#9ECAE1','ORF3':'#6BAED6','E':'#4292C6','M':'#2171B5','N':'#08306B'}
NT_COL={'A':'#000000','C':'#56B4E9','G':'#D55E00','U':'#8C8C8C'}
NT_LS={'A':'-','C':'--','G':'-','U':':'}
GROUPS=[('canonical','Can.','Canonical'),('TRSL_dep_noncanonical','Dep.','TRS-L-dep noncanonical'),
        ('TRSL_indep_distant','Dist.','TRS-L-indep distant'),('TRSL_indep_local','Loc.','TRS-L-indep local')]
GKEY={'canonical':'canonical','TRSL_dep_noncanonical':'dep','TRSL_indep_distant':'distant','TRSL_indep_local':'local'}
class Canvas:
    def __init__(self,fig,W_,H_,ox=0,oy=0):
        self.fig=fig; self.W=W_; self.H=H_; self.ox=ox; self.oy=oy
    def ax(self,x,y,w,h,**kw):
        return self.fig.add_axes([(self.ox+x)/self.W,1-(self.oy+y+h)/self.H,w/self.W,h/self.H],**kw)
    def text(self,x,y,s,**kw):
        kw.setdefault('fontsize',FS)
        return self.fig.text((self.ox+x)/self.W,1-(self.oy+y)/self.H,s,**kw)
    def letter(self,s):
        self.text(0.5,0.5,s,fontsize=13,fontweight='bold',ha='left',va='top')
    def title(self,s,x=6.5):
        self.text(x,1.2,s,fontsize=8,fontweight='bold',ha='left',va='top')
def violin(ax,data,pos,color,width=0.8,alpha=0.55,box=True,lw=0.6):
    data=np.asarray(data,float)
    if len(data)<2 or np.ptp(data)==0:
        return
    v=ax.violinplot([data],positions=[pos],widths=width,showextrema=False)
    for b in v['bodies']:
        b.set_facecolor(color); b.set_edgecolor(color); b.set_alpha(alpha); b.set_linewidth(lw)
    if box:
        q=np.percentile(data,[25,50,75]); iqr=q[2]-q[0]
        lo=data[data>=q[0]-1.5*iqr].min(); hi=data[data<=q[2]+1.5*iqr].max()
        ax.add_patch(Rectangle((pos-width*0.13,q[0]),width*0.26,q[2]-q[0],fc='white',ec='black',lw=0.5,zorder=3))
        ax.plot([pos-width*0.13,pos+width*0.13],[q[1],q[1]],color='black',lw=0.8,zorder=4)
        ax.plot([pos,pos],[lo,q[0]],color='black',lw=0.5,zorder=3); ax.plot([pos,pos],[q[2],hi],color='black',lw=0.5,zorder=3)
def star(p):
    return '***' if p<=0.001 else ('**' if p<=0.01 else ('*' if p<=0.05 else 'ns'))
def savefig(fig,base,dpi=300):
    for ext in ['png','pdf','svg']:
        fig.savefig(base+'.'+ext,dpi=dpi,facecolor='white')
