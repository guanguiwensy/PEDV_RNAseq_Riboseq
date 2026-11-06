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
import sys; sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from s4_panels import *
OUT=B+'figures/'
LAY={'A':(0,6,112,56),'B':(112,6,68,56),'C':(0,64,78,68),'D':(78,64,102,68),'E':(0,134,90,58),'F':(90,134,90,58),'G':(0,194,180,44)}
FN={'A':draw_A,'B':draw_B,'C':draw_C,'D':draw_D,'E':draw_E,'F':draw_F,'G':draw_G}
SRC={}
def combined():
    Wm,Hm=180,238
    fig=plt.figure(figsize=(Wm/25.4,Hm/25.4))
    cv0=Canvas(fig,Wm,Hm,0,0); legend_strip(cv0,y=0.8,x0=1)
    for k,(x,y,w,h) in LAY.items():
        cv=Canvas(fig,Wm,Hm,x,y); SRC[k]=FN[k](cv)
    savefig(fig,OUT+'Fig3_combined'); plt.close(fig)
def singles():
    for k,(x,y,w,h) in LAY.items():
        ex=0 if k in('A','G') else 8
        fig=plt.figure(figsize=(w/25.4,(h+ex+4)/25.4))
        cv0=Canvas(fig,w,h+ex+4,0,0); legend_strip(cv0,y=0.8,x0=1,maxw=w) if ex else None
        cv=Canvas(fig,w,h+ex+4,0,ex)
        FN[k](cv); savefig(fig,OUT+'Fig3%s'%k); plt.close(fig)
if __name__=='__main__':
    combined()
    if len(sys.argv)>1 and sys.argv[1]=='all': singles()
    for k,v in SRC.items():
        v.to_csv(OUT+'Fig3%s_source_data.csv'%k,index=False)
    print('ok')
