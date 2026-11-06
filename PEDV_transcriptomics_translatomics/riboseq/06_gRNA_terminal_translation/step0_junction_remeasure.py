"""Step 0: re-measure canonical leader-body junctions on GI-b (refdata_GIb_v1.1).
Input: GIb_rnaseq_reanalysis/junction_atlas/junction_events_raw.tsv.gz (STAR split reads on GI-b), PEDV_GIb.fa
Convention: leader_end = last leader base (1-based), body_start = first body base (1-based), as reported by STAR (left-flush)."""
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
import sample_registry as SR
import pandas as pd, numpy as np, gzip, sys, collections
B=PEDV_P+''
OUT=B+'fig5_terminal_ORF/'
ev=pd.read_csv(B+'GIb_rnaseq_reanalysis/junction_atlas/junction_events_raw.tsv.gz',sep='\t',
    usecols=['sample','leader_start','leader_end','body_start','gap_index_in_read','mapq','NH','is_duplicate','is_secondary','is_supplementary'])
seq=''.join(l.strip() for l in open(B+'refbuild_GIb/v1.1/virus/PEDV_GIb.fa') if not l.startswith('>')).upper()
print('events',len(ev),'genome',len(seq))
SITES={'S':20597,'ORF3':24718,'E':25342,'M':25597,'N':26285}
ev=ev[ev.leader_end.between(40,80)]
ev=ev[~ev.is_secondary&~ev.is_supplementary]
def gene_of(b):
    for g,p in SITES.items():
        if abs(b-p)<=12: return g
    return None
ev['gene']=ev.body_start.map(gene_of)
c=ev[ev.gene.notna()].copy()
print('canonical-window events',len(c))
man=pd.read_csv(B+'GIb_rnaseq_reanalysis/sample_manifest_GIb.tsv',sep='\t')
man=man[man.sample_id.isin(SR.WHITELIST_LIBRARIES)].set_index('sample_id'); SR.assert_exactly(man.index, SR.RNA_VERO_ALIASES, 'step0 junction remeasure')
def ob(i): return seq[i-1]  # 1-based
def amb_window(d,b):
    """equivalent junctions (d+k,b+k): ref[d+1..d+k]==ref[b..b+k-1]; (d-k,b-k): ref[d-k+1..d]==ref[b-k..b-1]"""
    r=0
    while ob(d+r+1)==ob(b+r): r+=1
    l=0
    while ob(d-l)==ob(b-l-1): l+=1
    return l,r
rows=[];persamp=[]
for g,p in SITES.items():
    cg=c[c.gene==g]
    pair=cg.groupby(['leader_end','body_start']).size().sort_values(ascending=False)
    (md,mb)=pair.index[0]
    dmode=cg.leader_end.value_counts(); bmode=cg.body_start.value_counts()
    l,r=amb_window(md,mb)
    # per sample modes
    ps=cg.groupby("sample")[["leader_end","body_start"]].apply(lambda x: pd.Series({'n':len(x),
        'donor_mode':x.leader_end.mode().iat[0],'body_mode':x.body_start.mode().iat[0],
        'pair_mode_frac':x.groupby(['leader_end','body_start']).size().max()/len(x)}))
    ps['gene']=g; persamp.append(ps.reset_index())
    top=(pair.iloc[0])
    rows.append(dict(sgRNA=g,donor_leader_end_mode=md,body_acceptor_first_base_mode=mb,
        reads_total=len(cg),reads_modal_pair=int(top),modal_pair_frac=round(top/len(cg),4),
        donor_mode_alone=int(dmode.index[0]),donor_mode_frac=round(dmode.iloc[0]/len(cg),4),
        body_mode_alone=int(bmode.index[0]),body_mode_frac=round(bmode.iloc[0]/len(cg),4),
        n_samples_with_reads=ps.shape[0],n_samples_donor_mode_equal=int((ps.donor_mode==md).sum()),
        n_samples_body_mode_equal=int((ps.body_mode==mb).sum()),
        min_reads_per_sample=int(ps.n.min()),median_reads_per_sample=float(ps.n.median()),
        ambiguity_shift_left_nt=l,ambiguity_shift_right_nt=r,ambiguity_window_width_nt=l+r+1,
        left_most_equivalent=f'{md-l}|{mb-l}',right_most_equivalent=f'{md+r}|{mb+r}',
        other_pairs_top3=';'.join(f'{a}|{b}:{n}' for (a,b),n in pair.iloc[1:4].items())))
res=pd.DataFrame(rows); ps=pd.concat(persamp)
ps['time_h']=ps['sample'].map(man.time_h); ps['drug']=ps['sample'].map(man.treatment)
ps.to_csv(OUT+'_work/junction_per_sample.tsv',sep='\t',index=False)
res.to_csv(OUT+'_work/junction_summary_raw.tsv',sep='\t',index=False)
pd.set_option('display.width',250,'display.max_columns',40)
print(res.T)
print(ps.groupby('gene').agg(nsamp=('sample','nunique'),donor_vals=('donor_mode',lambda x:sorted(set(x))),body_vals=('body_mode',lambda x:sorted(set(x))),minfrac=('pair_mode_frac','min')))
# sequences around
for g,p in SITES.items():
    md=int(res[res.sgRNA==g].donor_leader_end_mode.iat[0]);mb=int(res[res.sgRNA==g].body_acceptor_first_base_mode.iat[0])
    print(g,'leader',seq[md-20:md],'|',seq[mb-1:mb+19],'  upstream-of-body',seq[mb-21:mb-1])

# ---------------- final frozen table ----------------
def mm(a,b): return sum(x!=y for x,y in zip(a,b))
def best_core(lo,hi):
    best=min(((mm(seq[i-1:i+5],'CTAAAC'),i) for i in range(lo,hi)),key=lambda t:(t[0],t[1]))
    return best
old={'S':('62|62','20597'),'ORF3':('61|61','24718'),'E':('55|55','25340'),'M':('66|66','25597'),'N':('66|66','26285')}
gstart={'S':20610,'ORF3':24758,'E':25361,'M':25608,'N':26300}
fin=[]
for _,r in res.iterrows():
    g=r.sgRNA; d=int(r.donor_leader_end_mode); b=int(r.body_acceptor_first_base_mode)
    cg=c[c.gene==g]; wr=int(r.ambiguity_shift_right_nt)
    sh=cg.groupby('sample').apply(lambda x:((x.leader_end==d)&(x.body_start==b)).mean(),include_groups=False)
    chim=seq[:d]+'|'+seq[b-1:b+24]
    lo=b-6; hi=b+wr+30
    bc=best_core(b+wr-5,b+wr+40)  # 6-mers starting at/after right-most junction
    note=''
    if g=='ORF3': note='junction坐标实测(97.4%读段,23/23样本一致);TRS-B核心(CTAAAC<=1错配)在junction下游仍未找到,最佳6-mer见best_downstream_core_*;TRS-B序列身份待确认,不影响转录本嵌合序列'
    if g=='M': note='三峰microheterogeneity(同一gap=25531,donor 63/66/67对body 25594/25597/25598);众数66|25597占约42%,三峰嵌合序列互不相同(差3nt);下游以众数建模,其余两峰作敏感性'
    if g=='E': note='旧TRS.tsv写55|25340(实为次峰,约2.5%);实测众数57|25342(约97%)'
    fin.append(dict(sgRNA=g,junction_status='实测(GIb v1.1, STAR split reads)' if g!='ORF3' else '坐标实测;TRS-B身份待确认',
        leader_end_donor_left_aligned=d,body_start_acceptor_left_aligned=b,gap_nt=b-d-1,
        reads_total=int(r.reads_total),reads_modal_pair=int(r.reads_modal_pair),modal_pair_frac=r.modal_pair_frac,
        n_samples=23,n_samples_modal_pair_is_sample_mode=int(r.n_samples_body_mode_equal),
        per_sample_modal_frac_min=round(sh.min(),4),per_sample_modal_frac_median=round(sh.median(),4),
        min_reads_per_sample=int(r.min_reads_per_sample),median_reads_per_sample=r.median_reads_per_sample,
        ambiguity_left_nt=int(r.ambiguity_shift_left_nt),ambiguity_right_nt=wr,ambiguity_window_width=int(r.ambiguity_window_width_nt),
        rightmost_equivalent=r.right_most_equivalent,
        homology_leader=f'{d+1}-{d+wr}' if wr else 'none',homology_body=f'{b}-{b+wr-1}' if wr else 'none',
        chimera_leader20__body25=seq[:d][-20:]+'__'+seq[b-1:b+24],
        other_top_pairs=r.other_pairs_top3,old_TRS_tsv_donor_body=f'{old[g][0]} / {old[g][1]}',
        canonical_main_ORF_start=gstart[g],
        best_downstream_core_mismatch=bc[0],best_downstream_core_pos=bc[1],best_downstream_core_seq=seq[bc[1]-1:bc[1]+5],
        note=note))
f=pd.DataFrame(fin)
hdr='# GIb_junction_remeasured v1.0 FROZEN | refdata_GIb_v1.1 | source: GIb_rnaseq_reanalysis/junction_atlas/junction_events_raw.tsv.gz (STAR split reads, 23 RNA-seq libs, no dedup) | convention: 1-based, donor=last leader base, acceptor=first body base, STAR left-flush (left-aligned) | script: step0_junction_remeasure.py\n'
with open(OUT+'GIb_junction_remeasured.tsv','w',encoding='utf-8') as fh:
    fh.write(hdr); f.to_csv(fh,sep='\t',index=False)
# M variants
m=c[c.gene=='M']; print(m.groupby(['leader_end','body_start']).size().sort_values(ascending=False).head(4))
msh=m[m.body_start.isin([25594,25597,25598])].groupby(['sample','body_start']).size().unstack(fill_value=0); msh=msh.div(msh.sum(1),axis=0); print(msh.describe().loc[['min','50%','max']])
print(f.T)
