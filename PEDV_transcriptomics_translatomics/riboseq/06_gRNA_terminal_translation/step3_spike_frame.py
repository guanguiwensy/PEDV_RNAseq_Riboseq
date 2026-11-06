"""Step 3: region-level spike / frame-consistency / read-length characterisation (5' leader, 5'UTR, 3'UTR; all libraries)."""
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
import pandas as pd, numpy as np
B = PEDV_P+''
OUT = B + 'fig5_terminal_ORF/'
RGP = B + 'riboseq_GIb_mapping_psite/'
SM = [(SR.lib("RIBO_VERO_CHX_1MOI_5h"), 5, 'CHX'), (SR.lib("RIBO_VERO_HAR_1MOI_5h"), 5, 'HAR'), (SR.lib("RIBO_VERO_CHX_1MOI_8h"), 8, 'CHX'), (SR.lib("RIBO_VERO_HAR_1MOI_8h"), 8, 'HAR'), (SR.lib("RIBO_VERO_CHX_1MOI_16h"), 16, 'CHX'), (SR.lib("RIBO_VERO_HAR_1MOI_16h"), 16, 'HAR')]
lo = pd.read_csv(RGP + 'psite_counts_PEDV_GIb.tsv.gz', sep='\t', comment='#', low_memory=False)
lo = lo[lo.strand == '+'].copy()
lo['read_length'] = pd.to_numeric(lo.read_length, errors='coerce')
lo = lo.dropna(subset=['read_length']); lo['read_length'] = lo.read_length.astype(int)
rls = pd.read_csv(RGP + 'readlength_sets_GIb.tsv', sep='\t', comment='#')
rls['read_length'] = pd.to_numeric(rls.read_length, errors='coerce'); rls = rls.dropna(subset=['read_length']); rls['read_length'] = rls.read_length.astype(int)
core = rls[rls.in_strict_core.astype(str) == 'True'][['sample_id', 'read_length']]
sc = lo.merge(core, on=['sample_id', 'read_length'])
sc['end5'] = sc.position - sc.offset
REG = {'leader_1_72': (1, 72), 'UTR5_1_292': (1, 292), 'UTR3_27626_27959': (27626, 27959), 'N_CDS_ref_26345_27580': (26345, 27580)}
rows = []
for sid, t, d in SM:
    x = sc[sc.sample_id == sid]
    tot_lib = x.psite_count.sum()
    f29_lib = x[x.read_length == 29].psite_count.sum() / tot_lib
    for reg, (a, b) in REG.items():
        g = x[(x.position >= a) & (x.position <= b)]
        pc = g.groupby('position').psite_count.sum().sort_values(ascending=False)
        n = pc.sum()
        fr = g.groupby(g.position % 3).psite_count.sum()
        fr = fr.reindex([0, 1, 2], fill_value=0)
        # frame share after removing top1 / top3 positions
        def frame_after(k):
            gg = g[~g.position.isin(pc.index[:k])]
            f = gg.groupby(gg.position % 3).psite_count.sum().reindex([0, 1, 2], fill_value=0)
            return f.max() / f.sum() if f.sum() else np.nan
        skew_ok = sid != SR.lib("RIBO_VERO_HAR_1MOI_5h")
        top1 = pc.index[0]
        f29_top1 = g[(g.position == top1) & (g.read_length == 29)].psite_count.sum() / g[g.position == top1].psite_count.sum()
        f29_reg = g[g.read_length == 29].psite_count.sum() / n
        rows.append(dict(sample=sid, time_h=t, drug=d, region=reg, n_psites=int(n), ppm=n / tot_lib * 1e6,
                         top1_pos=int(top1), top1_frac=pc.iloc[0] / n, top3_frac=pc.iloc[:3].sum() / n,
                         max_frame_frac=fr.max() / n, max_frame_frac_deTop1=frame_after(1), max_frame_frac_deTop3=frame_after(3),
                         resid_deTop1=int(n - pc.iloc[0]), resid_deTop3=int(n - pc.iloc[:3].sum()),
                         f29_region=f29_reg if skew_ok else np.nan, f29_top1pos=f29_top1 if skew_ok else np.nan, f29_library=f29_lib if skew_ok else np.nan,
                         length_skew_testable=skew_ok))
df = pd.DataFrame(rows)
df.to_csv(OUT + '_work/region_spike_frame_characterization.tsv', sep='\t', index=False)
pd.set_option('display.width', 250, 'display.max_columns', 30)
print(df.round(3).to_string())
# 5'-end anchoring of the position 21/22 pile
pile = sc[(sc.position <= 40)]
rows = []
for sid, t, d in SM:
    g = pile[pile.sample_id == sid]
    e = g.groupby('end5').psite_count.sum().sort_values(ascending=False)
    top = g[g.position.isin([21, 22])]
    ln = top.groupby('read_length').psite_count.sum()
    rows.append(dict(sample=sid, time_h=t, drug=d, strict_core_psites_pos21_22=int(top.psite_count.sum()),
                     top_5prime_end=int(e.index[0]), top_5prime_end_frac_of_pos1_40=e.iloc[0] / e.sum(),
                     read_lengths_at_pos21_22=';'.join(f'{k}:{int(v)}' for k, v in ln.items())))
all_len = lo[(lo.position.between(21, 22)) & (lo.read_length.between(25, 33))].merge(
    pd.DataFrame(SM, columns=['sample_id', 't', 'd']), on='sample_id')
alll = all_len.assign(end5=all_len.position - all_len.offset).groupby(['sample_id', 'read_length', 'end5']).psite_count.sum().reset_index()
alll.to_csv(OUT + '_work/pos21_22_pile_all_lengths_end5.tsv', sep='\t', index=False)
pd.DataFrame(rows).to_csv(OUT + '_work/pos21_22_pile_summary.tsv', sep='\t', index=False)
print(pd.DataFrame(rows).round(3).to_string())
