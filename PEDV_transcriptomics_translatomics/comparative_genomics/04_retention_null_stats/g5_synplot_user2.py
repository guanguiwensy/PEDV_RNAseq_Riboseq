"""SynPlot2-style test on the USER'S de-duplicated gene alignments (<private-workstation-path> read-only) using the USER'S tree-based pair lists (pairs.txt, NJ peripheral traversal).
Per gene: (1) codon frame of the alignment is detected (best protein match of alignment rows to the GI-b ORF over frame offsets 0/1/2); (2) for every listed pair, per codon observed synonymous substitutions
and expected = syn sites x (pair-wide S_d/S_s), binomial variance; summed over pairs -> obs, exp, var per alignment codon (SynPlot2 logic; recomputed because the stored ORF1ab frasynplot_plot.txt is out of frame);
(3) alignment codons are mapped to GI-b v1.1 host-ORF codons by protein alignment (best rows, union fill); (4) candidate-span z and same-gene same-width window permutation p; controls 8-40 aa; positive control PRF.
Validation: recomputed obs/exp are compared to the user's frasynplot_plot.txt where its frame is correct."""
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
from scipy import stats
from Bio.Seq import Seq
from Bio import Align
U=PEDV_R+"PEDV_ribo/"
GENES={'ORF1a':('ORF1a/PEDV_ORF1a_polyprotein_extracted_aln_remove_seq.fas','ORF1a/pairs.txt','ORF1a/frasynplot_plot.txt'),
 'ORF1b':('ORF1ab/PEDV_ORF1ab_polyprotein_tail_extracted_aln_remove_seq.fas','ORF1ab/pairs.txt','ORF1ab/frasynplot_plot.txt'),
 'S':('spike/PEDV_spike_protein_extracted_aln_remove_seq.fas','spike/pairs.txt','spike/frasynplot_plot.txt'),
 'ORF3':('ORF3/PEDV_ORF3_extracted_aln.fasta','ORF3/pairs.txt','ORF3/frasynplot_plot.txt'),
 'E':('envelope/PEDV_envelope_protein_extracted_aln.fasta','envelope/pairs.txt','envelope/frasynplot_plot.txt'),
 'M':('membrane/PEDV_membrane_protein_extracted_aln.fasta','membrane/pairs.txt','membrane/frasynplot_plot.txt'),
 'N':('nucleocapsid/PEDV_nucleocapsid_protein_extracted_aln_remove_seq.fas','nucleocapsid/pairs.txt','nucleocapsid/frasynplot_plot.txt')}
def readfa(p):
    ids,seqs=[],[]
    for l in open(p,encoding='utf-8',errors='replace'):
        l=l.rstrip()
        if l.startswith('>'): ids.append(l[1:]); seqs.append([])
        elif l: seqs[-1].append(l)
    return ids,[''.join(s).upper() for s in seqs]
ns=np.zeros(64); nd=np.zeros((64,64),np.int8); issyn=np.zeros((64,64),bool)
for c in range(64):
    s=codon_list[c]
    if s in STOPS: continue
    for k in range(3):
        for nb in 'ACGT':
            if nb==s[k]: continue
            t=s[:k]+nb+s[k+1:]
            if t not in STOPS and tbl.get(t)==tbl.get(s): ns[c]+=1/3
for a in range(64):
    for b in range(64):
        nd[a,b]=sum(x!=y for x,y in zip(codon_list[a],codon_list[b]))
        issyn[a,b]=(a!=b) and not isstop[a] and not isstop[b] and aa_of[a]==aa_of[b]
al=Align.PairwiseAligner(); al.mode='local'; al.match_score=2; al.mismatch_score=-1; al.open_gap_score=-4; al.extend_gap_score=-1
cmap={'A':0,'C':1,'G':2,'T':3}
res={}; qc=[]
for g,(fa,pr,pl) in GENES.items():
    ids_,seqs=readfa(U+fa); L=len(seqs[0]); a,b=CDS[g]; gb=str(Seq(ref[a-1:b-3]).translate())
    def prot(s,f):
        n=(len(s)-f)//3; return ''.join(str(Seq(s[f+3*i:f+3*i+3]).translate()) if set(s[f+3*i:f+3*i+3])<=set('ACGT') else 'X' for i in range(n))
    sc={f:max(al.score(prot(s,f).replace('*','X'),gb) for s in seqs[:12]) for f in (0,1,2)}
    f=max(sc,key=sc.get); ncod=(L-f)//3
    # codon code matrix (rows x ncod), -1 for gap/ambiguous
    C=np.full((len(seqs),ncod),-1,np.int16)
    for r,s in enumerate(seqs):
        for i in range(ncod):
            t=s[f+3*i:f+3*i+3]
            if set(t)<=set('ACGT'): C[r,i]=cmap[t[0]]*16+cmap[t[1]]*4+cmap[t[2]]
    # pairs
    idx={sid[:50]:i for i,sid in enumerate(ids_)}
    pairs=[]
    for l in open(U+pr,encoding='utf-8',errors='replace'):
        if l.startswith('#') or not l.strip(): continue
        p=l.split();
        if len(p)>=2 and p[0][:50] in idx and p[1][:50] in idx: pairs.append((idx[p[0][:50]],idx[p[1][:50]]))
    pairs=list(dict.fromkeys(pairs))
    obs=np.zeros(ncod); ex=np.zeros(ncod); var=np.zeros(ncod)
    for i,j in pairs:
        v=(C[i]>=0)&(C[j]>=0); ci=np.where(v,C[i],0); cj=np.where(v,C[j],0); v&=~isstop[ci]&~isstop[cj]
        sy=np.where(v&issyn[ci,cj],nd[ci,cj],0).astype(float); sp=np.where(v,(ns[ci]+ns[cj])/2,0.0); Ss=sp.sum()
        if Ss==0: continue
        ps=sy.sum()/Ss; obs+=sy; ex+=sp*ps; var+=sp*ps*(1-ps)
    # mapping alignment codon -> GI-b codon (union over best rows)
    scores=[]
    for r,s in enumerate(seqs):
        scores.append((al.score(prot(s,f).replace('*','X'),gb),r))
    scores.sort(reverse=True); m=np.full(ncod,-1); used=set()
    for sc_,r in scores[:12]:
        pr_=prot(seqs[r],f).replace('*','X'); aln=al.align(pr_,gb)[0]
        for (s0,s1),(t0,t1) in zip(aln.aligned[0],aln.aligned[1]):
            for k in range(s1-s0):
                i=s0+k; kk=t0+k+1
                if m[i]<0 and kk not in used: m[i]=kk; used.add(kk)
    nco=len(gb)
    O=np.full(nco+1,np.nan); E_=np.full(nco+1,np.nan); V=np.full(nco+1,np.nan)
    for i in range(ncod):
        if m[i]>0: O[m[i]]=obs[i]; E_[m[i]]=ex[i]; V[m[i]]=var[i]
    with np.errstate(divide='ignore',invalid='ignore'): zc=(obs-ex)/np.sqrt(np.where(var>0,var,np.nan))
    phi=max(1.0,float(np.nanmean(zc**2)))
    userplot=np.loadtxt(U+pl); cmp_=np.nan
    if f==0 and len(userplot)==ncod:
        m_=(userplot[:,2]>0)&(ex>0); cmp_=float(np.corrcoef(userplot[m_,1],obs[m_])[0,1]) if m_.sum()>10 else np.nan
    res[g]=dict(obs=O,ex=E_,var=V,phi=phi)
    qc.append(dict(gene=g,n_sequences=len(seqs),n_pairs=len(pairs),alignment_columns=L,codon_frame_offset=f,frame_scores=str(sc),aln_codons=ncod,GIb_codons=nco,GIb_codons_mapped=int((~np.isnan(O[1:])).sum()),
        obs_total=float(obs.sum()),exp_total=float(ex.sum()),user_plot_obs_total=float(userplot[:,1].sum()),corr_obs_with_user_plot=cmp_,phi=round(phi,2)))
Q=pd.DataFrame(qc); Q.to_csv(OUT+'source_data/G5_synplot_user_mapping_qc.csv',index=False); pd.set_option('display.width',250); print(Q.drop(columns=['frame_scores']).to_string())
def stat(g,lo,hi):
    r=res[g]; a,b=CDS[g]
    k0=max(1,(lo-a)//3+1); k1=min(len(r['obs'])-1,(hi-a)//3+1)
    idx=np.array([k for k in range(k0,k1+1) if not np.isnan(r['obs'][k])])
    if len(idx)<5: return None
    d=(r['obs'][idx]-r['ex'][idx]).sum(); v=r['var'][idx].sum()
    if v<=0: return None
    z=d/np.sqrt(r['phi']*v)
    ok=~np.isnan(r['obs'][1:]); dd=np.where(ok,r['obs'][1:]-r['ex'][1:],0.0); vv=np.where(ok,r['var'][1:],0.0)
    cd=np.r_[0,np.cumsum(dd)]; cv=np.r_[0,np.cumsum(vv)]; span=k1-k0+1
    dw=cd[span:]-cd[:-span]; vw=cv[span:]-cv[:-span]; good=vw>0
    zw=dw[good]/np.sqrt(r['phi']*vw[good])
    return dict(n_codons=len(idx),obs=float(r['obs'][idx].sum()),exp=float(r['ex'][idx].sum()),obs_over_exp=float(r['obs'][idx].sum()/max(r['ex'][idx].sum(),1e-9)),z=float(z),p_lower=float(stats.norm.cdf(z)),p_window_perm=float((1+(zw<=z).sum())/(1+len(zw))),host_phi=r['phi'])
f_=families(); iso=isoforms(); isog=iso.groupby('family_id')
ia=f_[(f_.primary_category=='internal alternative-frame ORF')&f_.fixed_parent_assignment.isin(CDS)].copy()
rows=[]
for _,F in ia.iterrows():
    g=F.fixed_parent_assignment; d=isog.get_group(F.family_id); thr=41 if F.max_peptide_aa>40 else 8
    st=int(d[d.peptide_length_aa>=thr].start_position.min()); stop=int(d.stop_position.iloc[0])+2
    s=stat(g,st,stop)
    if s is None: continue
    rows.append(dict(family_id=F.family_id,role='target_gt40' if F.max_peptide_aa>40 else 'control_8to40',region=g,rel_frame=f'+{int(F.relative_frame_to_parent)}',max_aa=F.max_peptide_aa,span_start=st,span_end=stop,**s))
T=pd.DataFrame(rows)
def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); r=p[o]*n/(np.arange(n)+1); r=np.minimum.accumulate(r[::-1])[::-1]; q=np.empty(n); q[o]=np.minimum(r,1); return q
for role in T.role.unique():
    m=T.role==role; T.loc[m,'q_lower']=bh(T.loc[m,'p_lower']); T.loc[m,'q_window_perm']=bh(T.loc[m,'p_window_perm'])
T.to_csv(OUT+'source_data/G5_synplot_user_per_family.csv',index=False,encoding='utf-8-sig')
pc=[]
for g,lo,hi in (('ORF1b',12592,12700),('ORF1a',12540,12622)):
    s=stat(g,lo,hi); pc.append(dict(element='PRF',host=g,lo=lo,hi=hi,**(s or {})))
pd.DataFrame(pc).to_csv(OUT+'source_data/G5_synplot_user_positive_control_PRF.csv',index=False)
S=[]
for role,d in T.groupby('role'):
    S.append(dict(role=role,n=len(d),n_obs_lt_exp=int((d.obs_over_exp<1).sum()),n_p_lt_005=int((d.p_lower<.05).sum()),n_BH_lt_005=int((d.q_lower<.05).sum()),n_winperm_p_lt_005=int((d.p_window_perm<.05).sum()),n_winperm_BH_lt_005=int((d.q_window_perm<.05).sum()),median_z=d.z.median(),median_obs_over_exp=d.obs_over_exp.median()))
S=pd.DataFrame(S); S.to_csv(OUT+'source_data/G5_synplot_user_summary.csv',index=False)
t=T[T.role=='target_gt40']; c=T[T.role=='control_8to40']
print(S.round(3).to_string()); print('MWU z targets vs controls',stats.mannwhitneyu(t.z,c.z).pvalue); print(pd.DataFrame(pc).round(4).to_string())
print(t.groupby('region').agg(n=('z','size'),lt=('obs_over_exp',lambda x:(x<1).sum()),med_z=('z','median'),sig=('p_lower',lambda x:(x<.05).sum()),perm=('p_window_perm',lambda x:(x<.05).sum())))
print(t.sort_values('p_window_perm').head(10)[['family_id','region','n_codons','obs','exp','obs_over_exp','z','p_lower','q_lower','p_window_perm','q_window_perm']].round(4).to_string())

# ---------------- pooled test (power): sum obs / sum exp over all candidate spans vs the same span widths placed at random in the same gene (10,000 draws, seed 20261022)
rng=np.random.default_rng(20261022)
def pooled(role):
    d=T[T.role==role]; spans=[]
    for _,x in d.iterrows():
        g=x.region; a,b=CDS[g]; k0=max(1,(int(x.span_start)-a)//3+1); k1=min(len(res[g]['obs'])-1,(int(x.span_end)-a)//3+1); spans.append((g,k0,k1-k0+1))
    cum={}
    for g in res:
        o=np.nan_to_num(res[g]['obs'][1:]); e=np.nan_to_num(res[g]['ex'][1:]); ok=(~np.isnan(res[g]['obs'][1:])).astype(float)
        cum[g]=(np.r_[0,np.cumsum(o)],np.r_[0,np.cumsum(e)],np.r_[0,np.cumsum(ok)])
    def ratio(starts):
        O=E=0.0
        for (g,k0,w),s in zip(spans,starts):
            c=cum[g]; O+=c[0][s-1+w]-c[0][s-1]; E+=c[1][s-1+w]-c[1][s-1]
        return O/E
    obs_r=ratio([k0 for _,k0,_ in spans])
    nul=np.empty(10000)
    for r in range(10000):
        st=[int(rng.integers(1,len(res[g]['obs'])-w+1)) for g,_,w in spans]; nul[r]=ratio(st)
    return dict(role=role,n_spans=len(spans),pooled_obs_over_exp=float(obs_r),null_mean=float(nul.mean()),null_p2_5=float(np.percentile(nul,2.5)),null_p97_5=float(np.percentile(nul,97.5)),p_lower=float((1+(nul<=obs_r).sum())/10001),p_upper=float((1+(nul>=obs_r).sum())/10001),seed=20261022,n_draws=10000)
PL=pd.DataFrame([pooled('target_gt40'),pooled('control_8to40')]); PL.to_csv(OUT+'source_data/G5_synplot_user_pooled.csv',index=False); print(PL.round(4).to_string())
