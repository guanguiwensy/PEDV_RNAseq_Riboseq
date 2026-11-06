"""Homology-mapped family assignments: for every genome x reference family, evaluate stop-anchor retention and start classes
using per-genome alignment-derived coordinate maps (mapper.py). Not-assessable (end-incomplete / ambiguous base) is kept
separate from structural loss and excluded from denominators downstream."""
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
import sys,csv,json,collections,numpy as np,warnings
warnings.filterwarnings('ignore')
B=PEDV_K+r""; O=B+"Supplementary_Figure_L/"; W=O+"work/"
z=np.load(W+"refmaps.npz"); ids=list(z['ids']); M=z['M']
gen={};cur=None
for l in open(B+"N1_conservation_analysis/results/PEDV_genomes_all.fasta",encoding='utf-8'):
    if l[0]=='>': cur=l[1:].split()[0]; gen[cur]=[]
    else: gen[cur].append(l.strip().upper())
gen={k:''.join(v) for k,v in gen.items()}
pool=json.load(open(W+"refpool.json"))['fams']
STOPS={"TAA","TAG","TGA"}; STARTS={"ATG","CTG","GTG","TTG","ATA","ATT","ATC","ACG"}
ACGT=set("ACGT")
typing={r['accession']:r for r in csv.DictReader(open(W+"typing_S_based.tsv",encoding='utf-8'),delimiter='\t')}
out=[]; iso_out=[]
SA=collections.Counter()
ST=typing  # alias

for gi,a in enumerate(ids):
    g=gen[a]; m=M[gi].tolist(); Lg=len(g)
    # per-position codon annotation in genome
    isstop=bytearray(Lg); isstart=bytearray(Lg); amb=bytearray(Lg)
    for p in range(Lg-2):
        c=g[p:p+3]
        if c in STOPS: isstop[p]=1
        elif c in STARTS: isstart[p]=1
        elif not (c[0] in ACGT and c[1] in ACGT and c[2] in ACGT): amb[p]=1
    def trip(p):   # ref 0-based codon start -> (status,gpos)
        x,y,w=m[p],m[p+1],m[p+2]
        if x==-2 or y==-2 or w==-2: return 'END',None
        if x==-1 or y==-1 or w==-1: return 'DEL',None
        if y!=x+1 or w!=y+1: return 'DISCONT',None
        if amb[x]: return 'AMB',None
        return 'OK',x
    for F in pool:
        s0=F['stop']-1
        st,sx=trip(s0)
        base=dict(accession=a,family_id=F['fid'])
        n_cls=collections.Counter(); n40=collections.Counter(); hom40=0; hom40_na=0; hom40_tot=sum(1 for x in F['isos'] if x['aa']>40)
        if st in('END','AMB'):
            status8=status40='not_assessable'; reason='end_incomplete_or_unaligned' if st=='END' else 'ambiguous_base_in_stop_codon'
            maxaa=''; nany=''
            for x in F['isos']:
                n_cls['not_assessable']+=1
                if x['aa']>40: n40['not_assessable']+=1; SA[(F['fid'],x['start'],typing[a]['S_type'],'not_assessable')]+=1
        elif st in('DEL','DISCONT') or not isstop[sx]:
            status8=status40='lost'; reason={'DEL':'stop_codon_region_deleted','DISCONT':'indel_within_stop_codon','OK':'stop_codon_substituted_nonstop'}[st]
            maxaa=''; nany=''
            for x in F['isos']:
                n_cls['loss_of_eligible_start']+=1
                if x['aa']>40: n40['loss_of_eligible_start']+=1; SA[(F['fid'],x['start'],typing[a]['S_type'],'loss_of_eligible_start')]+=1
        else:
            f=sx%3
            # previous in-frame stop
            q=sx-3; ambig=False
            while q>=0 and not isstop[q]:
                if amb[q]: ambig=True
                q-=3
            lo=q+3   # first codon of open stretch
            starts=[p for p in range(lo,sx,3) if isstart[p]]
            aas=[(sx-p)//3 for p in starts]
            maxaa=max(aas) if aas else 0; nany=len(starts)
            # ambiguity upstream hiding possible stop? only matters if no start found
            def status(thr):
                if any(x>=thr for x in aas): return 'retained'
                return 'unresolved' if (ambig and lo<=0) else 'lost'
            status8=status('8') if False else ('retained' if any(x>=8 for x in aas) else ('unresolved' if ambig else 'lost'))
            status40='retained' if any(x>40 for x in aas) else ('unresolved' if ambig else 'lost')
            reason='stop_retained' if status8=='retained' else ('no_eligible_start_in_open_stretch' if status8=='lost' else 'ambiguity_in_open_stretch')
            sset=set(starts)
            for x in F['isos']:
                rp=x['start']-1
                sst,sp=trip(rp)
                if sst in('END','AMB'): cl='not_assessable'
                elif sst in('DEL','DISCONT'): cl='loss_of_eligible_start'
                else:
                    rc=x['start_codon'].replace('U','T'); gc=g[sp:sp+3]
                    inframe_open=(sp%3==f) and lo<=sp<sx
                    if not inframe_open: cl='loss_of_eligible_start'
                    elif gc in STARTS:
                        if gc==rc: cl='AUG_retention' if rc=='ATG' else 'exact_codon_retention'
                        elif gc!='ATG' and rc!='ATG': cl='near_cognate_retention'
                        else: cl='start_class_switching'
                    else: cl='loss_of_eligible_start'
                n_cls[cl]+=1
                if x['aa']>40: n40[cl]+=1; SA[(F['fid'],x['start'],typing[a]['S_type'],cl)]+=1
                if x['aa']>40:
                    if cl in('exact_codon_retention','AUG_retention','near_cognate_retention','start_class_switching'): hom40+=1
                    elif cl=='not_assessable': hom40_na+=1
        hom_status=('NA' if st in('END','AMB') else ('lost' if (st in('DEL','DISCONT') or not isstop[sx]) else ('retained' if hom40>0 else ('NA' if hom40_na>0 and hom40_na==hom40_tot else 'lost'))))
        base.update(dict(hom_status_gt40=hom_status if hom40_tot>0 else '',n_hom_starts_retained_gt40=hom40,status_ge8=status8,status_gt40=status40,reason=reason,max_peptide_aa_in_genome=maxaa,n_eligible_starts_in_open_stretch=nany,
            n_ref_starts=F['n'],n_ref_starts_gt40=sum(1 for x in F['isos'] if x['aa']>40),**{k+'_gt40':n40.get(k,0) for k in ('exact_codon_retention','AUG_retention','near_cognate_retention','start_class_switching','loss_of_eligible_start','not_assessable')},**{k:n_cls.get(k,0) for k in ('exact_codon_retention','AUG_retention','near_cognate_retention','start_class_switching','loss_of_eligible_start','not_assessable')}))
        out.append(base)
cols=list(out[0])
with open(O+"PEDV_family_homology_mapping.csv","w",encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['S_type']+cols); w.writeheader()
    for r in out: r['S_type']=typing[r['accession']]['S_type']; w.writerow(r)
print(len(out))
c=collections.Counter(r['status_gt40'] for r in out); print(c)

# per start variant x lineage table (aa>40 isoforms)
iso={(F['fid'],x['start']):x for F in pool for x in F['isos'] if x['aa']>40}
lin=['GI-a','GI-b','GII-a','GII-b','GII-c','unresolved']
tot=collections.defaultdict(int)
for (f,s,l,c),n in SA.items(): tot[(f,s,l)]+=n
rows=[]
for (f,s),x in sorted(iso.items(),key=lambda kv:(kv[0][1],kv[0][0])):
    for l in lin:
        d=dict(family_id=f,start_position=s,reference_start_codon=x['start_codon'],reference_peptide_aa=x['aa'],lineage=l,n_genomes_in_lineage=sum(1 for a in typing if typing[a]['S_type']==l))
        n=tot.get((f,s,l),0)
        for c in ('exact_codon_retention','AUG_retention','near_cognate_retention','start_class_switching','loss_of_eligible_start','not_assessable'): d[c]=SA.get((f,s,l,c),0)
        d['n_assessable_denominator']=n-d['not_assessable']
        rows.append(d)
with open(O+"SupL_start_codon_conservation.tsv","w",encoding='utf-8',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0]),delimiter='	'); w.writeheader(); w.writerows(rows)
print('start table rows',len(rows))
