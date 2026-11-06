"""G2: stop-anchored family max-peptide-length distribution, +1 vs +2 (relative frame, Sup L rule A) with composition-null reference."""
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
f=families(); iso=isoforms()
d=f[f.relative_frame_to_parent.isin([1,2])&f.fixed_parent_assignment.isin(CDS)].copy()
d['rel_frame']=d.relative_frame_to_parent.map({1:'+1',2:'+2'})
d.rename(columns={'fixed_parent_assignment':'region','shared_stop_position':'stop','max_peptide_aa':'max_aa','min_peptide_aa':'min_aa','n_start_variants':'n_start_isoforms'})[['family_id','region','rel_frame','stop','max_aa','min_aa','n_start_isoforms','n_AUG_starts','primary_category']].to_csv(OUT+'source_data/G2_family_lengths_plus1_plus2.csv',index=False,encoding='utf-8-sig')
a=d[d.rel_frame=='+1'].max_peptide_aa; b=d[d.rel_frame=='+2'].max_peptide_aa
res=dict(n_plus1=len(a),n_plus2=len(b),median_plus1=float(a.median()),median_plus2=float(b.median()),mean_plus1=float(a.mean()),mean_plus2=float(b.mean()),
  mwu_p=float(stats.mannwhitneyu(a,b).pvalue),ks_p=float(stats.ks_2samp(a,b).pvalue),frac_gt20_plus1=float((a>20).mean()),frac_gt20_plus2=float((b>20).mean()),frac_gt40_plus1=float((a>40).mean()),frac_gt40_plus2=float((b>40).mean()),
  max_plus1=int(a.max()),max_plus2=int(b.max()))
z=np.load(OUT+'work/null_replicates.npz'); o=np.load(OUT+'work/null_observed.npz')
for rfi,rf in enumerate(('+1','+2')):
    nf=z['nfam'][:,:,rfi,:].sum(1)   # reps x tiers (ge8, ge20, gt40)
    for ti,t in ((1,'ge20'),(2,'gt40')):
        fr=nf[:,ti]/np.maximum(nf[:,0],1); ob=o['nfam'][0,:,rfi,ti].sum()/o['nfam'][0,:,rfi,0].sum()
        res[f'null_frac_{t}_of_ge8_{rf}']=[float(fr.mean()),float(np.percentile(fr,2.5)),float(np.percentile(fr,97.5)),float(ob)]
# null histogram (pooled over regions, per replicate mean) vs observed for the figure later
h=z['hist'].sum(1)/int(z['nrep'])  # (region,rf,aa)->mean per replicate; sum regions below
hn=np.stack([ (z['hist'][:,rf,:].sum(0))/int(z['nrep']) for rf in range(2)]); ho=np.stack([o['hist'][:,rf,:].sum(0) for rf in range(2)])
pd.DataFrame({'aa':np.arange(hn.shape[1]),'null_mean_p1':hn[0],'obs_p1':ho[0],'null_mean_p2':hn[1],'obs_p2':ho[1]}).query('aa>=8').to_csv(OUT+'source_data/G2_length_hist_obs_vs_null.csv',index=False)
json.dump(res,open(OUT+'source_data/G2_summary.json','w'),indent=1); print(json.dumps(res,indent=1))
