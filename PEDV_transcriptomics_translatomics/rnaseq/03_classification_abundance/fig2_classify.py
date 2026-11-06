#!/usr/bin/env python3
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
# Fig2 junction classification (GI-b reference, 20 pooled RNA-seq samples)
import re, sys, json
import pandas as pd, numpy as np
BASE=PEDV_P+'GIb_rnaseq_reanalysis/'
OUT=BASE+'figures/Fig2_work/'
LEADER=(int(sys.argv[1]),int(sys.argv[2])) if len(sys.argv)>2 else (55,71)
DIST=int(sys.argv[3]) if len(sys.argv)>3 else 1000
TAG=sys.argv[4] if len(sys.argv)>4 else 'main'
CANON={'S':20597,'ORF3':24718,'E':25342,'M':25597,'N':26285}; TOL=5
CDS={'ORF1a':(293,12622),'ORF1b':(12592,20613),'S':(20610,24758),'ORF3':(24758,25033),
     'E':(25361,25600),'M':(25608,26288),'N':(26300,27625)}
seq=''.join(l.strip() for l in open(BASE+'../refbuild_GIb/v1.1/virus/PEDV_GIb.fa') if not l.startswith('>'))
GL=len(seq)
man=pd.read_csv(BASE+'sample_manifest_GIb.tsv',sep='\t')
s20=man[man.batch=='batch3'].sample_id.tolist()
df=pd.read_csv(BASE+'junction_atlas/junction_events_raw.tsv.gz',sep='\t')
d=df[df['sample'].isin(s20)].copy()
d=d[(~d.is_duplicate)&(~d.is_secondary)&(~d.is_supplementary)&(d.NH==1)]
n_rows=len(d)
R=re.compile(r'(\d+)([MIDNSH=X])')
def anchors(cig,gap):
    ops=[(int(l),o) for l,o in R.findall(cig)]
    idx=[i for i,(l,o) in enumerate(ops) if o=='N' and l==gap]
    if not idx: return (0,0)
    i=idx[0]
    la=0
    for l,o in reversed(ops[:i]):
        if o in 'M=X': la+=l
        elif o in 'NSH': break
    ba=0
    for l,o in ops[i+1:]:
        if o in 'M=X': ba+=l
        elif o in 'NSH': break
    return la,ba
# leader anchor in the extractor summed all M before N; use same convention (sum, not break) for consistency
def anchors2(cig,gap):
    ops=[(int(l),o) for l,o in R.findall(cig)]
    idx=[i for i,(l,o) in enumerate(ops) if o=='N' and l==gap]
    if not idx: return (0,0)
    i=idx[0]
    la=sum(l for l,o in ops[:i] if o in 'M=X')
    ba=0
    for l,o in ops[i+1:]:
        if o in 'M=X': ba+=l
        elif o in 'NSH': break
    return la,ba
ab=[anchors2(c,g) for c,g in zip(d.cigar,d.gap_len)]
d['la']=[a[0] for a in ab]; d['ba']=[a[1] for a in ab]
qc=(d.mapq>=20)&(d.la>=15)&(d.ba>=15)
n_qcfail=int((~qc).sum())
d=d[qc].copy()
d['donor']=d.leader_end; d['acc']=d.body_start
n_qc_rows=len(d)
u=d.drop_duplicates(['sample','read_name','donor','acc']).copy()   # mates spanning same junction counted once
n_frag=len(u)
ev=u.groupby(['donor','acc']).agg(reads=('read_name','size'),n_samples=('sample','nunique')).reset_index()
per_sample=u.groupby(['sample','donor','acc']).size().rename('reads').reset_index()

def cds_of(pos):
    hits=[(s,g) for g,(s,e) in CDS.items() if s<=pos<=e]
    if not hits: return None
    s,g=max(hits)   # overlap: later-starting (downstream) ORF
    return g
def phase(pos,g): return (pos-CDS[g][0])%3
def first_atg(a):
    i=seq.find('ATG',a-1); return i+1 if i>=0 else None
def region(pos):
    g=cds_of(pos)
    if g: return g
    if pos<293: return "5'UTR"
    if pos>27625: return "3'UTR"
    return 'intergenic'
rows=[]
for r in ev.itertuples():
    dn,ac=r.donor,r.acc
    dep=LEADER[0]<=dn<=LEADER[1]
    canon=None
    for g,s in CANON.items():
        if abs(ac-s)<=TOL: canon=g
    dist=ac-dn-1
    y=cds_of(ac); x=cds_of(dn)
    frame='NA'; orf_cons='NA'; aug_any='NA'
    if not (dep and canon):
        if y is None:
            frame='noncoding'
            p=first_atg(ac)
            orf_cons='noncoding_acceptor'
            aug_any = 'first_AUG_is_CDS_start' if any(p==v[0] for v in CDS.values()) else 'other'
        else:
            p=first_atg(ac)
            if p is not None and p<=CDS[y][1]:
                ph=(p-CDS[y][0])%3
                frame='in' if ph==0 else 'out'
                orf_cons=('full_length' if p==CDS[y][0] else 'N_trunc_internal_AUG') if ph==0 else 'frameshift'
            else:
                frame='no_AUG_in_acceptor_CDS'; orf_cons='no_AUG_in_acceptor_CDS'
                # scanning ribosome continues past CDS end: is there any AUG downstream in genome, and its frame vs covering CDS
                if p is None: aug_any='none_in_genome'
                else:
                    z=cds_of(p)
                    aug_any = 'AUG_outside_CDS' if z is None else ('AUG_in_frame_with_'+z if (p-CDS[z][0])%3==0 else 'AUG_offframe_in_'+z)
    if dep:
        cls='canonical' if canon else 'TRSL_dep_noncanonical'
        if canon: cat=f'canonical_{canon}'
        else: cat=cls+'_'+frame
    else:
        cls='TRSL_indep_'+('local' if dist<DIST else 'distant')
        cat=None
    # legacy donor-frame continuation (supplementary column)
    cont='NA'
    if not dep:
        cont='noncoding' if (x is None or y is None) else ('in' if (phase(dn,x)+1)%3==phase(ac,y) else 'out')
        cat=cls+'_'+cont
    rows.append(dict(donor=dn,acceptor=ac,reads=r.reads,n_samples=r.n_samples,gap_nt=dist,
        trsl_dependent=dep,canonical_gene=canon or '',cls=cls,frame=frame,category=cat,
        donor_region=region(dn),acceptor_region=region(ac),orf_consequence=orf_cons,downstream_AUG_note=aug_any,donor_frame_continuation=cont))
E=pd.DataFrame(rows)
E.to_csv(OUT+f'events_{TAG}.tsv',sep='\t',index=False)
per_sample.to_csv(OUT+'per_sample_events.tsv.gz',sep='\t',index=False)
tot=E.reads.sum()
cat=E.groupby('category').reads.sum().sort_values(ascending=False)
info=dict(leader=LEADER,dist=DIST,rows_20samples_reliable=n_rows,rows_fail_qc=n_qcfail,rows_qc=n_qc_rows,
          fragment_junctions=n_frag,total=int(tot),events=len(E))
print(json.dumps(info))
out=pd.DataFrame({'reads':cat,'pct':cat/tot*100}); print(out.round(3).to_string()); print('sum',out.pct.sum())
