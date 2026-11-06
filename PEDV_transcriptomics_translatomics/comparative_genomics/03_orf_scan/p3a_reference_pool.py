"""Reference-derived candidate pool on refdata_GIb_v1.1 (all 3 genomic frames, 8 start codons, >=8 aa, complete ORFs).
Fixed-assignment rule (declared before analysis): parent = annotated CDS containing the start codon; if >1 CDS contains it
(ORF1a/ORF1b 12592-12622, S/ORF3 24758) prefer the CDS in which the ORF is in relative frame +1, then +2, then 0, then genomic order.
Full overlap annotation is kept in `overlapping_annotations`."""
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
import sys,csv,collections,json
sys.path[0:0] = [REPO_ROOT + '/comparative_genomics/02_alignment_homology_lineage', REPO_ROOT + '/comparative_genomics/03_orf_scan', REPO_ROOT + '/comparative_genomics/04_retention_null_stats']
from scan_core import *
B=PEDV_K+r""; O=B+"Supplementary_Figure_L/"
ref="".join(l.strip() for l in open(PEDV_P+r"refbuild_GIb/v1.1/virus/PEDV_GIb.fa") if l[0]!='>').upper()
CDS=[("ORF1a",293,12622),("ORF1b",12592,20613),("S",20610,24758),("ORF3",24758,25033),("E",25361,25600),("M",25608,26288),("N",26300,27625)]
PRF=(12586,12622)
comp,trunc=scan_seq(ref)
def gframe(p): return (p-1)%3
CSTOP={b-2:(n,a,b) for n,a,b in CDS}   # first base of each annotated CDS stop codon
def classify(p,s):
    end=s+2; f=gframe(p)
    ov=[n for n,a,b in CDS if a<=end and b>=p]
    if s in CSTOP:   # stop anchor is the stop codon of an annotated CDS (same frame by construction)
        n,a,b=CSTOP[s]
        if p==a: cat='annotated canonical CDS start'
        elif p<a: cat='in-frame extension'
        else: cat='in-frame internal initiation'
        return n,0,cat,ov,('upstream of '+n if p<a else n)
    cont=[(n,a,b) for n,a,b in CDS if a<=p<=b]
    def rel(a): return (f-gframe(a))%3
    if cont:
        pref=sorted(cont,key=lambda x:({1:0,2:1,0:2}[rel(x[1])],x[1]))[0]
        n,a,b=pref; r=rel(a)
        if r==0:
            if p==a: cat='annotated canonical CDS start'
            elif s==b-2: cat='in-frame internal initiation'
            else: cat='in-frame within CDS (stop differs)'
        else: cat='internal alternative-frame ORF'
        parent=n; relf=r; reg=n
    else:
        down=[(n,a,b) for n,a,b in CDS if a>p]
        ovl=[(n,a,b) for n,a,b in down if a<=end]
        if ovl:
            n,a,b=ovl[0]; r=rel(a)
            if r==0 and s==b-2: cat='in-frame extension'
            else: cat='upstream overlapping ORF'
            parent=n; relf=r; reg='5UTR/gap upstream of '+n
        elif down:
            cat='upstream ORF'; parent='-'; relf=None; reg=('5UTR' if p<293 else 'inter-CDS gap')
            parent=down[0][0]; relf=None
        else:
            cat='intergenic-UTR ORF'; parent='-'; relf=None; reg='3UTR'
    return parent,relf,cat,ov,reg
iso=[]
for r in comp:
    if r['aa']<8: continue
    p,s=r['start'],r['stop']
    parent,relf,cat,ov,reg=classify(p,s)
    iso.append(dict(start=p,stop=s,genomic_frame=r['frame'],aa=r['aa'],start_codon=r['start_codon'].replace('T','U'),stop_codon=r['stop_codon'].replace('T','U'),
        parent=parent,relative_frame=('' if relf is None else {0:'0',1:'+1',2:'+2'}[relf]),category=cat,region=reg,overlapping=';'.join(ov),
        prf_overlap='yes' if (p<=PRF[1] and s+2>=PRF[0]) else 'no',
        aug='yes' if r['start_codon']=='ATG' else 'no'))
fam=collections.defaultdict(list)
for r in iso: fam[(r['genomic_frame'],r['stop'])].append(r)
fams=[]
for (f,s),v in sorted(fam.items(),key=lambda kv:min(x['start'] for x in kv[1])):
    v.sort(key=lambda x:x['start']); top=v[0]
    fams.append(dict(frame=f,stop=s,isos=v,n=len(v),maxaa=max(x['aa'] for x in v),minaa=min(x['aa'] for x in v),first=top))
# stable IDs (frame-agnostic, coordinate-based; NOT renamed against Figure 6 names)
for k,F in enumerate(fams,1):
    F['fid']=f"SDF_f{F['frame']}_stop{F['stop']}"
    F['tier']='gt40' if F['maxaa']>40 else ('ge20' if F['maxaa']>=20 else 'ge8')
    for x in F['isos']: x['family']=F['fid']
print('isoforms',len(iso),'families',len(fams),collections.Counter(F['tier'] for F in fams))
print('gt40 isoforms',sum(1 for r in iso if r['aa']>40),' families with max>40',sum(1 for F in fams if F['maxaa']>40))
print(collections.Counter(F['first']['category'] for F in fams if F['maxaa']>40))
print(collections.Counter(F['first']['parent'] for F in fams if F['maxaa']>40))
json.dump(dict(fams=[{k:v for k,v in F.items()} for F in fams],trunc=trunc),open(O+"work/refpool.json","w"))
def wr(fn,rows,cols):
    with open(O+fn,"w",encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols); w.writeheader(); w.writerows(rows)
irows=[]
for F in fams:
    for x in F['isos']:
        irows.append(dict(family_id=F['fid'],start_position=x['start'],start_codon=x['start_codon'],stop_position=x['stop'],stop_codon=x['stop_codon'],
          peptide_length_aa=x['aa'],orf_length_nt=x['aa']*3+3,genomic_frame_static=x['genomic_frame'],relative_frame_to_parent=x['relative_frame'],
          fixed_parent_assignment=x['parent'],category=x['category'],region=x['region'],overlapping_annotations=x['overlapping'],
          overlaps_PRF_or_ORF1ab_overlap_region=x['prf_overlap'],is_AUG='yes' if x['aug']=='yes' else 'no',
          length_tier_gt40=('yes' if x['aa']>40 else 'no'),length_eq40=('yes' if x['aa']==40 else 'no'),length_tier_ge20=('yes' if x['aa']>=20 else 'no')))
wr("PEDV_reference_ORF_isoforms.csv",irows,list(irows[0]))
frows=[]
for F in fams:
    cats=collections.Counter(x['category'] for x in F['isos'])
    frows.append(dict(family_id=F['fid'],genomic_frame_static=F['frame'],shared_stop_position=F['stop'],stop_codon=F['isos'][0]['stop_codon'],
      n_start_variants=F['n'],n_AUG_starts=sum(1 for x in F['isos'] if x['aug']=='yes'),n_near_cognate_starts=sum(1 for x in F['isos'] if x['aug']=='no'),
      start_positions='|'.join(str(x['start']) for x in F['isos']),start_codons='|'.join(sorted({x['start_codon'] for x in F['isos']})),
      max_peptide_aa=F['maxaa'],min_peptide_aa=F['minaa'],length_tier=F['tier'],fixed_parent_assignment=F['first']['parent'],
      relative_frame_to_parent=F['first']['relative_frame'],primary_category=F['first']['category'],categories_all=';'.join(f"{k}:{v}" for k,v in cats.items()),
      region=F['first']['region'],overlapping_annotations=F['first']['overlapping'],overlaps_PRF_or_ORF1ab_overlap_region=('yes' if any(x['prf_overlap']=='yes' for x in F['isos']) else 'no'),
      canonical_stop_family=('yes' if F['stop'] in CSTOP else 'no'),sharing_stop_with=(CSTOP[F['stop']][0] if F['stop'] in CSTOP else ''),
      starts_span_regions='|'.join(sorted({x['region'] for x in F['isos']})),
      has_eq40_isoform=('yes' if any(x['aa']==40 for x in F['isos']) else 'no')))
wr("PEDV_reference_ORF_families.csv",frows,list(frows[0]))
trows=[dict(start_position=t['start'],genomic_frame_static=t['frame'],start_codon=t['start_codon'].replace('T','U'),open_aa_to_genome_end=t['aa_open'],status='boundary-truncated/incomplete') for t in trunc]
wr("PEDV_reference_boundary_truncated_starts.csv",trows,list(trows[0]))
