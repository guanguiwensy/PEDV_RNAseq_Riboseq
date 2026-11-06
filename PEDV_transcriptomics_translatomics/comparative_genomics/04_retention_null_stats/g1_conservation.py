"""G1 conservation: fraction of assessable genomes in which a reference alt-frame stop triplet is still a stop at the homologous position;
expected value under the synonymous-substitution null (analytic, exact over the outcomes of the <=2 codons spanned; k = each genome's own synonymous codon changes in that CDS; same ts/tv kappa as Sup L).
Assessable = all 3 positions mapped contiguously to unambiguous A/C/G/T. Unassessable (gap/end/IUPAC/indel junction -> NOT counted as non-conserved; internal deletions are counted as 'lost' only when mapped as deletion)."""
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
from r6lib import *
z=np.load(OUT+'work/proj.npz',allow_pickle=True); P=z['P']; contig=z['contig']; ids=list(z['ids'])
cl=pd.read_csv(SL+'work/dedup_clusters.tsv',sep='\t'); rep=cl.set_index('accession').loc[ids,'is_cluster_representative'].eq('yes').values
kappa=json.load(open(SL+'work/null_params.json'))['kappa']
OVL=(12592,12622)
def ts(a,b): return {a,b} in ({'A','G'},{'C','T'})
alts={}
for s in code:
    if s in STOPS or tbl.get(s) in ('M','W'): continue
    L=[]
    for k in range(3):
        for nb in 'ACGT':
            if nb==s[k]: continue
            t=s[:k]+nb+s[k+1:]
            if t in STOPS or tbl.get(t)!=tbl.get(s): continue
            L.append((code[t],kappa if ts(s[k],nb) else 1.0))
    if L: tot=sum(w for _,w in L); alts[code[s]]=[(c,w/tot) for c,w in L]
body={}
for n,(a,b) in CDS.items():
    lo,hi=a,b-3
    if n=='ORF1a': hi=OVL[0]-1
    if n=='ORF1b': lo=OVL[1]+1
    body[n]=(lo,hi)
# per genome syn counts per CDS (codon-level, as Sup L p5a)
G=P.shape[0]; K={}; NC={}
for n,(a,b) in CDS.items():
    nco=(b-a+1)//3-1; p0=np.arange(a-1,a-1+3*nco,3)
    rc=np.array([code[ref[p:p+3]] for p in p0])
    ok=contig[:,p0]&contig[:,p0+1]&(P[:,p0]<4)&(P[:,p0+1]<4)&(P[:,p0+2]<4)
    gc=P[:,p0].astype(int)*16+P[:,p0+1]*4+P[:,p0+2]
    syn=ok&(gc!=rc[None,:])&(aa_of[np.where(ok,gc,0)]==aa_of[rc][None,:])&(~isstop[np.where(ok,gc,0)])
    K[n]=syn.sum(1); NC[n]=sum(1 for c in rc if c in alts)
rows=[]; sites=[]
for n in ORDER:
    a,b=CDS[n]; lo,hi=body[n]
    for p in range(lo,hi-1):
        f=(p-a)%3
        if f==0: continue
        cod=ref[p-1:p+2]
        if cod not in STOPS: continue
        if not set(ref[p-1:p+2])<=set('ACGT'): continue
        # genome state
        tri=P[:,p-1:p+2]; ok=(tri<4).all(1)&contig[:,p-1]&contig[:,p]
        gc=tri[:,0].astype(int)*16+tri[:,1]*4+tri[:,2]
        kept=ok&isstop[np.where(ok,gc,0)]
        # codon grid: j = (p-a)//3 ; triplet spans codons j and j+1
        j=(p-a)//3; c0=code[ref[a-1+3*j:a+3*j+2]]; c1=code[ref[a+3*(j+1)-1:a+3*(j+1)+2]]
        cand0=c0 in alts; cand1=c1 in alts
        # outcome coefficients
        def stopfrac(nt0,nt1,alt0,alt1):
            return None
        def tri_of(x0,x1):
            s=codon_list[x0]+codon_list[x1]; return s[f:f+3]
        c01=sum(w for c,w in alts.get(c1,[]) if tri_of(c0,c)in STOPS)
        c10=sum(w for c,w in alts.get(c0,[]) if tri_of(c,c1)in STOPS)
        c11=sum(w0*w1 for cc0,w0 in alts.get(c0,[]) for cc1,w1 in alts.get(c1,[]) if tri_of(cc0,cc1)in STOPS)
        q=np.minimum(K[n],NC[n])/NC[n]
        q0=q if cand0 else 0*q; q1=q if cand1 else 0*q
        Pexp=(1-q0)*(1-q1)+(1-q0)*q1*c01+q0*(1-q1)*c10+q0*q1*c11
        sites.append(dict(region=n,rel_frame=f'+{f}',pos=p,stop=cod,n_assess=int(ok.sum()),n_kept=int(kept.sum()),cons=kept.sum()/max(ok.sum(),1),
            n_assess_cl=int((ok&rep).sum()),n_kept_cl=int((kept&rep).sum()),exp=float(Pexp[ok].mean()) if ok.any() else np.nan,exp_cl=float(Pexp[ok&rep].mean()) if (ok&rep).any() else np.nan))
S=pd.DataFrame(sites); S.to_csv(OUT+'source_data/G1_stop_site_conservation_per_site.csv',index=False,encoding='utf-8-sig')
def summ(df,label):
    out=[]
    for (r,f),d in df.groupby(['region','rel_frame']):
        out.append(dict(scope=label,region=r,rel_frame=f,n_stop_sites=len(d),pooled_kept=int(d.n_kept.sum()),pooled_assessable=int(d.n_assess.sum()),pooled_conservation=d.n_kept.sum()/d.n_assess.sum(),
          mean_site_conservation=d.cons.mean(),median_site_conservation=d.cons.median(),frac_sites_ge99=(d.cons>=.99).mean(),
          pooled_expected_syn_null=(d.exp*d.n_assess).sum()/d.n_assess.sum(),n_assessable_genomes_range=f'{int(d.n_assess.min())}-{int(d.n_assess.max())}'))
    return out
rows=[]
rows+=summ(S,'all_1083_genomes')
S2=S.copy(); S2['n_assess']=S2.n_assess_cl; S2['n_kept']=S2.n_kept_cl; S2['cons']=S2.n_kept/S2.n_assess.clip(lower=1); S2['exp']=S2.exp_cl
rows+=summ(S2,'cluster_representatives_579')
A=S.assign(region='ALL'); rows+=summ(A,'all_1083_genomes'); A2=S2.assign(region='ALL'); rows+=summ(A2,'cluster_representatives_579')
R=pd.DataFrame(rows); R.to_csv(OUT+'source_data/G1_stop_conservation_by_frame_region.csv',index=False,encoding='utf-8-sig')
pd.set_option('display.width',250); print(R.round(4).to_string())
# paired nonparametric comparison of site-level conservation +1 vs +2 (sites are not independent of phylogeny; descriptive)
from scipy import stats
a=S[S.rel_frame=='+1'].cons; b=S[S.rel_frame=='+2'].cons
print('MWU site cons +2 vs +1',stats.mannwhitneyu(b,a,alternative='two-sided'))
