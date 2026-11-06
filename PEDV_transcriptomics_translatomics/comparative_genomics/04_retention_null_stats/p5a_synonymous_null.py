"""Synonymous-substitution null (generalises the idea of N1 scripts/phase2g_step4_null1_unbiased.py: the protein
sequence of the parent CDS is held fixed, only synonymous codon changes are introduced).
Differences from Null-1: (i) works for any parent CDS and any family lying inside it; (ii) the number of changes per
replicate equals each genome's observed number of synonymous codon changes in that parent CDS vs the reference;
(iii) alternatives are single-nt synonymous changes weighted by the empirical transition/transversion ratio.
Statistic: homologous-start family retention (stop codon intact AND >=1 reference start with aa>40 whose codon is still
eligible and with no in-frame stop between start and stop). Applied only to internal alternative-frame families inside a
parent CDS; all other families are descriptive only."""
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
from Bio.Data.CodonTable import unambiguous_dna_by_id
B=PEDV_K+r""; O=B+"Supplementary_Figure_L/"; W=O+"work/"
ref="".join(l.strip() for l in open(PEDV_P+r"refbuild_GIb/v1.1/virus/PEDV_GIb.fa") if l[0]!='>').upper()
CDS={"ORF1a":(293,12622),"ORF1b":(12592,20613),"S":(20610,24758),"ORF3":(24758,25033),"E":(25361,25600),"M":(25608,26288),"N":(26300,27625)}
tbl=unambiguous_dna_by_id[1].forward_table
STOP={"TAA","TAG","TGA"}; STARTS={"ATG","CTG","GTG","TTG","ATA","ATT","ATC","ACG"}
z=np.load(W+"refmaps.npz"); ids=list(z['ids']); M=z['M']
gen={};cur=None
for l in open(B+"N1_conservation_analysis/results/PEDV_genomes_all.fasta",encoding='utf-8'):
    if l[0]=='>': cur=l[1:].split()[0]; gen[cur]=[]
    else: gen[cur].append(l.strip().upper())
gen={k:''.join(v) for k,v in gen.items()}
pool=json.load(open(W+"refpool.json"))['fams']
codon_list=[a+b+c for a in "ACGT" for b in "ACGT" for c in "ACGT"]
code={s:i for i,s in enumerate(codon_list)}
aa_of=np.array([tbl.get(codon_list[i],'*') for i in range(64)])
isstop=np.array([codon_list[i] in STOP for i in range(64)]); isstart=np.array([codon_list[i] in STARTS for i in range(64)])
enc=np.full(256,255,np.uint8)
for i,c in enumerate("ACGT"): enc[ord(c)]=i
def transition(a,b): return {a,b} in ({'A','G'},{'C','T'})
syn_counts={}; ts=tv=0
for pn,(a,b) in CDS.items():
    n=(b-a+1)//3
    p0=np.arange(a-1,a-1+3*n,3)
    rc=np.array([code[ref[p:p+3]] for p in p0])
    for gi,acc in enumerate(ids):
        m=M[gi]; m0,m1,m2=m[p0],m[p0+1],m[p0+2]
        ok=(m0>=0)&(m1==m0+1)&(m2==m1+1)
        e=enc[np.frombuffer(gen[acc].encode(),np.uint8)]
        idx=np.where(ok,m0,0); Lm=len(e)-1
        x,y,w=e[idx],e[np.minimum(idx+1,Lm)],e[np.minimum(idx+2,Lm)]
        ok&=(x<4)&(y<4)&(w<4)
        gc=x.astype(int)*16+y.astype(int)*4+w.astype(int)
        gc=np.where(ok,gc,0)
        diff=ok&(gc!=rc)
        syn=diff&(aa_of[gc]==aa_of[rc])&(~isstop[gc])
        syn_counts[(acc,pn)]=int(syn.sum())
        if pn=='ORF1a' and gi%5==0:
            for j in np.nonzero(syn)[0][:100]:
                r=codon_list[rc[j]]; q=codon_list[gc[j]]
                d=[(u,v) for u,v in zip(r,q) if u!=v]
                if len(d)==1:
                    if transition(*d[0]): ts+=1
                    else: tv+=1
kappa=ts/max(1,tv); print('ts',ts,'tv',tv,'ts/tv',round(kappa,2))
alts={}
for s in code:
    if s in STOP or tbl.get(s) in ('M','W'): continue
    L=[]
    for k in range(3):
        for nb in "ACGT":
            if nb==s[k]: continue
            t=s[:k]+nb+s[k+1:]
            if t in STOP or tbl.get(t)!=tbl.get(s): continue
            L.append((code[t],kappa if transition(s[k],nb) else 1.0))
    if L: alts[code[s]]=L
rng=np.random.default_rng(20261012)
refcode=enc[np.frombuffer(ref.encode(),np.uint8)].astype(np.int64)
ASS=[F for F in pool if F['isos'][0]['category']=='internal alternative-frame ORF' and F['isos'][0]['parent'] in CDS and F['maxaa']>40]
print('families for null',len(ASS),dict(collections.Counter(F['isos'][0]['parent'] for F in ASS)))
def fam_ok(arr,F):
    s0=F['stop']-1
    if not isstop[arr[s0]*16+arr[s0+1]*4+arr[s0+2]]: return False
    st=[x['start']-1 for x in F['isos'] if x['aa']>40]
    smin=min(st)
    cod=arr[smin:s0+3].reshape(-1,3); cv=cod[:,0]*16+cod[:,1]*4+cod[:,2]
    stp=np.nonzero(isstop[cv])[0]; inner=stp[stp<len(cv)-1]
    last=inner.max() if len(inner) else -1
    for p in st:
        k=(p-smin)//3
        if k>last and isstart[cv[k]]: return True
    return False
R=40
nullp={}
for pn,(a,b) in CDS.items():
    fams=[F for F in ASS if F['isos'][0]['parent']==pn]
    if not fams: continue
    n=(b-a+1)//3-1
    p0=np.arange(a-1,a-1+3*n,3)
    cand=np.array([i for i in range(n) if code[ref[p0[i]:p0[i]+3]] in alts])
    ks=sorted({syn_counts[(acc,pn)] for acc in ids})
    # coarse binning of k to keep runtime bounded: simulate at bin centres of width 5, look up by nearest centre
    centres=sorted({5*round(k/5) for k in ks})
    print(pn,'unique k',len(ks),'sim bins',len(centres),'candidate codons',len(cand),flush=True)
    for kc in centres:
        kk=min(kc,len(cand)); cnt=collections.Counter()
        for r in range(R):
            arr=refcode.copy()
            sel=rng.choice(cand,kk,replace=False)
            for i in sel:
                pos=p0[i]; c=code[ref[pos:pos+3]]; Lc=alts[c]
                wts=np.array([x[1] for x in Lc]); j=rng.choice(len(Lc),p=wts/wts.sum())
                t=codon_list[Lc[j][0]]
                arr[pos]=enc[ord(t[0])]; arr[pos+1]=enc[ord(t[1])]; arr[pos+2]=enc[ord(t[2])]
            for F in fams:
                if fam_ok(arr,F): cnt[F['fid']]+=1
        nullp[(pn,kc)]={F['fid']:cnt[F['fid']]/R for F in fams}
rows=list(csv.DictReader(open(O+"PEDV_family_homology_mapping.csv",encoding='utf-8-sig')))
obs=collections.defaultdict(lambda:[0,0,0.0,0.0,[]])
fm={F['fid']:F for F in ASS}
for r in rows:
    F=fm.get(r['family_id'])
    if not F or r['hom_status_gt40'] not in('retained','lost'): continue
    pn=F['isos'][0]['parent']; k=syn_counts[(r['accession'],pn)]; kc=5*round(k/5)
    p=nullp[(pn,kc)][F['fid']]
    o=obs[F['fid']]; o[1]+=1; o[0]+= (r['hom_status_gt40']=='retained'); o[2]+=p; o[3]+=p*(1-p); o[4].append(p)
out=[]
rngp=np.random.default_rng(7)
for fid,(k,n,sp,var,pl) in obs.items():
    pa=np.array(pl); sim=(rngp.random((2000,len(pa)))<pa).sum(1)
    out.append(dict(family_id=fid,parent=fm[fid]['isos'][0]['parent'],n_assessable=n,n_retained_observed=k,observed_fraction=round(k/n,4),
        expected_fraction_synonymous_null=round(sp/n,4),expected_count=round(sp,1),null_sd_count=round(var**0.5,2),empirical_percentile_of_observed=round(float((sim<=k).mean()),4),parent_cds=fm[fid]['isos'][0]['parent'],z=round((k-sp)/max(var,1e-9)**0.5,2),null_replicates_per_bin=R))
with open(O+"work/synonymous_null_by_family.csv","w",encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
import statistics
print('families',len(out),'median obs',statistics.median(o['observed_fraction'] for o in out),'median exp',statistics.median(o['expected_fraction_synonymous_null'] for o in out))
print('obs>exp',sum(1 for o in out if o['observed_fraction']>o['expected_fraction_synonymous_null']),'obs<exp',sum(1 for o in out if o['observed_fraction']<o['expected_fraction_synonymous_null']))
json.dump(dict(kappa=kappa,ts=ts,tv=tv),open(O+"work/null_params.json","w"))
