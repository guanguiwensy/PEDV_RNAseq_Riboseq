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
import pandas as pd, numpy as np
B=PEDV_P+'GIb_rnaseq_reanalysis/'
W=B+'figures/Fig3_work/'
seq=''.join(l.strip() for l in open(B+'../refbuild_GIb/v1.1/virus/PEDV_GIb.fa') if not l.startswith('>')).upper()
GL=len(seq)
SEED=42
NPERM=10000
CLS_NAME={'canonical':'Canonical','TRSL_dep_noncanonical':'TRS-L-dep noncanonical',
          'TRSL_indep_distant':'TRS-L-indep distant','TRSL_indep_local':'TRS-L-indep local'}
def load_pool():
    m=pd.read_csv(W+'events_main_with_tier.tsv',sep='\t')
    hp=m.event_type.isin(['High-confidence','Probable'])
    m['tier']=m.event_type.replace({'canonical':'canonical(not_tiered)','internal_genomic_donor(classB)':'classB(no_tier)'})
    p=pd.Series(False,index=m.index)
    p[m.cls=='canonical']=True
    p[m.cls=='TRSL_dep_noncanonical']=hp[m.cls=='TRSL_dep_noncanonical']
    ind=m.cls.isin(['TRSL_indep_distant','TRSL_indep_local'])
    p[ind]=(m.n_samples>=2)[ind]
    m['pass_filter']=p
    m['filter_rule']=np.where(m.cls=='canonical','none(canonical)',np.where(m.cls=='TRSL_dep_noncanonical','High+Probable','>=2 positive samples (+MAPQ/anchor, Fig2)'))
    m['hotspot_27620_27682']=(m.donor==27620)&(m.acceptor==27682)
    m['canonical_primary']=False
    c=m[m.cls=='canonical']
    m.loc[c.loc[c.groupby('canonical_gene').reads.idxmax()].index,'canonical_primary']=True
    m['junction_id']=m.donor.astype(str)+'_'+m.acceptor.astype(str)
    return m
