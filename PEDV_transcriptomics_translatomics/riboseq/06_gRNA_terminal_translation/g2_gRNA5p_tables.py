"""gRNA 5' analysis, part 2: peak-sharing, 21/22 and hotspot characterisation, final evidence table, descriptive categories."""
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
import json, numpy as np, pandas as pd
B = PEDV_P+''
OUT = B + 'fig5_terminal_ORF/gRNA_5prime/'
RGP = B + 'riboseq_GIb_mapping_psite/'
ev = pd.read_pickle(OUT + '_ev.pkl')
cat = pd.read_pickle(OUT + '_cat.pkl')
meta = json.load(open(OUT + '_meta.json'))
HOT = {int(k): v for k, v in meta['HOT'].items()}; MST = meta['MST']
TPS = [5, 8, 16]
SAMP = {5: (SR.lib("RIBO_VERO_CHX_1MOI_5h"), SR.lib("RIBO_VERO_HAR_1MOI_5h")), 8: (SR.lib("RIBO_VERO_CHX_1MOI_8h"), SR.lib("RIBO_VERO_HAR_1MOI_8h")), 16: (SR.lib("RIBO_VERO_CHX_1MOI_16h"), SR.lib("RIBO_VERO_HAR_1MOI_16h"))}
SM = [(SR.lib("RIBO_VERO_CHX_1MOI_5h"), 5, 'CHX'), (SR.lib("RIBO_VERO_HAR_1MOI_5h"), 5, 'HAR'), (SR.lib("RIBO_VERO_CHX_1MOI_8h"), 8, 'CHX'), (SR.lib("RIBO_VERO_HAR_1MOI_8h"), 8, 'HAR'), (SR.lib("RIBO_VERO_CHX_1MOI_16h"), 16, 'CHX'), (SR.lib("RIBO_VERO_HAR_1MOI_16h"), 16, 'HAR')]

# ---- sharing of strongest peaks across candidates ----
for t in TPS:
    pk = ev[f'HAR_peak_pos_{t}h']
    ev[f'n_candidates_with_same_strongest_peak_{t}h'] = pk.map(pk.value_counts())
    ev[f'candidates_sharing_peak_{t}h'] = [','.join(ev.candidate_id[pk == p].str[-3:]) if not np.isnan(p) else '' for p in pk]

# ---- hotspot tables ----
lo = pd.read_csv(RGP + 'psite_counts_PEDV_GIb.tsv.gz', sep='\t', comment='#', low_memory=False)
lo = lo[lo.strand == '+'].copy()
lo['read_length'] = pd.to_numeric(lo.read_length, errors='coerce'); lo = lo.dropna(subset=['read_length']); lo['read_length'] = lo.read_length.astype(int)
rls = pd.read_csv(RGP + 'readlength_sets_GIb.tsv', sep='\t', comment='#')
rls['read_length'] = pd.to_numeric(rls.read_length, errors='coerce'); rls = rls.dropna(subset=['read_length']); rls['read_length'] = rls.read_length.astype(int)
core = rls[rls.in_strict_core.astype(str) == 'True'][['sample_id', 'read_length']]
sc = lo.merge(core, on=['sample_id', 'read_length'])
name = {s: f'{d}{t}' for s, t, d in SM}
reg = sc[sc.position <= 443]
tot = reg.groupby('sample_id').psite_count.sum()
HP = [21, 22, 23, 31, 102, 103, 104, 112, 113, 271, 289, 290, 292, 293, 294]
rows = []
for p in HP:
    o = dict(position=p, group=HOT.get(p, 'shoulder/none'), rel_frame_to_ORF1a=(p - MST) % 3, offset_from_ORF1a_AUG=p - MST,
             n_candidate_windows_containing=int(((ev.start_position - 15 <= p) & (p <= ev.start_position + 45)).sum()),
             n_candidate_windows_in_candidate_frame=int((((ev.start_position - 15 <= p) & (p <= ev.start_position + 45)) & ((p - ev.start_position) % 3 == 0)).sum()))
    for s, t, d in SM:
        g = reg[(reg.sample_id == s) & (reg.position == p)]
        o[f'{d}{t}_n'] = int(g.psite_count.sum()); o[f'{d}{t}_share_pct'] = 100 * g.psite_count.sum() / tot[s]
        o[f'{d}{t}_frac_len28_29_30'] = (f"{g[g.read_length == 28].psite_count.sum() / g.psite_count.sum():.2f}/" if g.psite_count.sum() else 'NA/') + \
                                       (f"{g[g.read_length == 29].psite_count.sum() / g.psite_count.sum():.2f}/{g[g.read_length == 30].psite_count.sum() / g.psite_count.sum():.2f}" if g.psite_count.sum() else 'NA/NA')
    o['n_libraries_share_ge_1pct'] = sum(o[f'{d}{t}_share_pct'] >= 1 for s, t, d in SM)
    o['recurs_in_both_drugs'] = any(o[f'CHX{t}_share_pct'] >= 1 for t in TPS) and any(o[f'HAR{t}_share_pct'] >= 1 for t in TPS)
    rows.append(o)
hot = pd.DataFrame(rows)
hot.to_csv(OUT + 'hotspot_recurrence_table.tsv', sep='\t', index=False)
print(hot[['position', 'group', 'rel_frame_to_ORF1a', 'n_candidate_windows_containing', 'n_candidate_windows_in_candidate_frame', 'n_libraries_share_ge_1pct', 'recurs_in_both_drugs']].to_string())

# ---- A: 21/22 region, read-end distribution ----
rowsA = []
for s, t, d in SM:
    m = pd.read_csv(RGP + f'_work/server/stats/{s}/pedv_5p_matrix.tsv.gz', sep='\t'); m = m[m.strand == '+']
    x = m[m.pos5.between(1, 40)]
    p10 = m[m.pos5 == 10]; p9 = m[m.pos5 == 9]
    allL = m.groupby('L').n.sum()
    orf = m[m.pos5.between(279, 283)]            # reads whose P-site lies at 291-295 (offset ~12): ORF1a start region
    end3 = (p10.assign(e3=p10.pos5 + p10.L - 1)).groupby('e3').n.sum().sort_values(ascending=False)
    ppm_pile = sc[(sc.sample_id == s) & sc.position.between(21, 23)].psite_count.sum() / tot[s] * 1e6
    ppm_orf = sc[(sc.sample_id == s) & sc.position.between(289, 296)].psite_count.sum() / tot[s] * 1e6
    rowsA.append(dict(library=f'{d} {t} h', reads_with_5prime_end_1_40=int(x.n.sum()), frac_5prime_end_at_10=p10.n.sum() / x.n.sum(), frac_5prime_end_at_9=p9.n.sum() / x.n.sum(),
                      frac_5prime_end_9_10=(p10.n.sum() + p9.n.sum()) / x.n.sum(),
                      softclip5_zero_frac_at_pos10=p10[p10.sc5 == 0].n.sum() / max(1, p10.n.sum()),
                      L28_frac_at_pos10=p10[p10.L == 28].n.sum() / max(1, p10.n.sum()), L29_30_frac_at_pos10=p10[p10.L.isin([29, 30])].n.sum() / max(1, p10.n.sum()),
                      top3prime_ends_at_pos10=';'.join(f'{int(k)}:{int(v)}' for k, v in end3.head(3).items()),
                      L28_frac_genome_wide=allL.get(28, 0) / allL.sum(), L28_frac_ORF1a_start_region=orf[orf.L == 28].n.sum() / max(1, orf.n.sum()),
                      L29_30_frac_ORF1a_start_region=orf[orf.L.isin([29, 30])].n.sum() / max(1, orf.n.sum()),
                      ppm_P21_23=ppm_pile, ppm_P289_296=ppm_orf))
A = pd.DataFrame(rowsA)
for t in TPS:
    h = A[A.library == f'HAR {t} h'].iloc[0]; c = A[A.library == f'CHX {t} h'].iloc[0]
    A.loc[A.library == f'HAR {t} h', 'HAR_over_CHX_ppm_P21_23'] = h.ppm_P21_23 / c.ppm_P21_23
    A.loc[A.library == f'HAR {t} h', 'HAR_over_CHX_ppm_P289_296'] = h.ppm_P289_296 / c.ppm_P289_296
A.to_csv(OUT + 'hotspot_21_22_characterization.tsv', sep='\t', index=False)
pd.set_option('display.width', 300, 'display.max_columns', 30)
print(A.round(3).T.to_string())

# ---- descriptive categories (manual, documented; no score) ----
CATS = {
 'gRNA_5p_001': ('C', 'strongest candidate-frame HAR peak is the 21-23 5-prime-end pile-up (shoulder); no independent start-proximal peak'),
 'gRNA_5p_002': ('C', 'strongest HAR peak is the 21-23 pile-up shoulder; CHX body reads moderate (8 occupied codons at 8 h) but candidate frame not dominant'),
 'gRNA_5p_003': ('A', 'HAR peak sits on the start codon itself (offset 0) at 8 h and 16 h (79/164 P-sites), outside every hotspot (3% hotspot-derived); CHX body at 8 h has 11/23 occupied codons; the 16 h CHX body is 98% hotspot 112 and is not counted as support. Caveat: CHX also has a peak on this start (8 h) and HAR/CHX window ratio is only 0.6 (8 h) and 1.3 (16 h), so HAR enrichment over CHX is not shown; 5 h has too few reads (6)'),
 'gRNA_5p_004': ('C', 'HAR peak = hotspot 103 at all three time points (70-82% of candidate-frame window reads); after masking the strongest peak moves to the start of another candidate (AUU@46)'),
 'gRNA_5p_005': ('C', 'HAR peak = hotspot 112/113 at 5 h and 8 h; 16 h peak is AUG@95 (a different, downstream-nested start)'),
 'gRNA_5p_006': ('C', 'HAR peak = hotspot 103 at all three time points (89-92% hotspot-derived)'),
 'gRNA_5p_007': ('C', 'HAR peak = hotspot 103 at all three time points (94-96% hotspot-derived)'),
 'gRNA_5p_008': ('A', 'AUG; at 16 h the strongest HAR peak (170) sits on its own start codon (offset 0, HAR/CHX 3.3) after hotspots are masked; at 5 h and 8 h the strongest peak is hotspot 112/113; CHX body at 8 h has 9/12 occupied codons, run of 6, 48% hotspot-derived. Single time point for the start-centred peak; not a reproducible-initiation claim'),
 'gRNA_5p_009': ('C', 'nested downstream of AUG@95 in the same ORF: HAR peaks are hotspot 112/113 (5/8 h) or the AUG@95 start (16 h, 6 nt upstream); CHX body is shared with AUG@95'),
 'gRNA_5p_010': ('C', 'nested in the AUG@95 ORF; peaks are hotspot 112/113 or the upstream AUG@95 start'),
 'gRNA_5p_011': ('C', 'nested in the AUG@95 ORF; peaks are hotspot 112/113 or the upstream AUG@95 start'),
 'gRNA_5p_012': ('C', 'HAR strongest peak = hotspot 103 shoulder (102) at 8 h and 16 h; CHX body reads very low (7-13)'),
 'gRNA_5p_013': ('D', 'HAR strongest peak is the 102 shoulder of hotspot 103; CHX body reads 0-9; too little signal'),
 'gRNA_5p_014': ('B', 'CHX body spread over 10/39 codons at 16 h (81 reads, span 99 nt, no hotspot contribution); HAR peaks small (3-54 P-sites) and shared with neighbouring candidates (162/168)'),
 'gRNA_5p_015': ('B', 'CHX body 8-10 occupied codons at 8 h and 16 h, no hotspot contribution; HAR peak small and not at the start'),
 'gRNA_5p_016': ('B', 'same as gRNA_5p_015 (shared body and shared 160 peak)'),
 'gRNA_5p_017': ('B', 'CHX body 9-10 occupied codons, no hotspot contribution; HAR peak small and shared with gRNA_5p_014/018 (162)'),
 'gRNA_5p_018': ('B', 'CHX body 9 occupied codons at 8 h and 16 h, no hotspot contribution; HAR peak (162, 54 P-sites at 16 h) is shared with gRNA_5p_014/017'),
 'gRNA_5p_019': ('D', '11-50 body reads, 2 occupied codons at most; HAR peaks 3-24 P-sites'),
 'gRNA_5p_020': ('B', 'CHX body 6-8 occupied codons over 51-75 nt, no hotspot contribution; HAR peak small'),
 'gRNA_5p_021': ('D', 'low counts (19-29 body reads; HAR peaks 7-54) and no start-aligned structure'),
 'gRNA_5p_022': ('D', 'low counts (16-28 body reads; HAR peaks 4-52) with unstable peak position'),
 'gRNA_5p_023': ('C', 'peaks and body are hotspot 271 (own start P-site, 98-100% of body reads) and, from 8 h, the ORF1a start region (292); cannot separate from the hotspot'),
 'gRNA_5p_024': ('C', 'canonical_start_confounded: strongest HAR peak is ORF1a AUG (290/293) at every time point (99-100% hotspot-derived)'),
}
ev['descriptive_category'] = ev.candidate_id.map(lambda c: CATS[c][0]); ev['category_rationale'] = ev.candidate_id.map(lambda c: CATS[c][1])

def note(r):
    bits = []
    for t in TPS:
        bits.append(f"{t}h:{r[f'peak_attribution_class_{t}h']}@{'NA' if np.isnan(r[f'HAR_peak_pos_{t}h']) else int(r[f'HAR_peak_pos_{t}h'])}(n={int(r[f'HAR_peak_height_{t}h'])})")
    return '; '.join(bits)
famm = cat.set_index('candidate_id').nested_starts_in_same_ORF
ev['ORF_family_members'] = ev.candidate_id.map(famm)
ev['category_rationale'] = ev.category_rationale + ' | same ORF (shared stop, shared body): ' + ev.ORF_family_members.str.replace('gRNA_5p_', '', regex=False)
ev['attribution_note'] = ev.apply(note, axis=1)
ev.to_pickle(OUT + '_ev2.pkl')

# ---- final evidence table (requested names first, exploratory last) ----
out = pd.DataFrame({'candidate_id': ev.candidate_id, 'start': ev.start_position, 'codon': ev.start_codon, 'stop': ev.stop_position, 'stop_codon': ev.stop_codon,
                    'length_aa': ev.length_aa, 'class': ev.ORF_class, 'frame_to_ORF1a': ev.relative_frame_to_ORF1a, 'distance_to_ORF1a_start': ev.distance_to_ORF1a_start,
                    'start_type': ev.start_type, 'n_starts_in_same_ORF': ev.n_starts_in_same_ORF})
for t in TPS:
    out[f'HAR{t}_candidate_reads'] = ev[f'HAR{t}_candidate_reads']
for t in TPS:
    for k in ['HAR_total_reads', ]: pass
for t in TPS:
    out[f'HAR{t}_total_reads'] = ev[f'HAR{t}_total_reads']
for t in TPS:
    for k in [f'HAR_peak_pos_{t}h', f'HAR_peak_offset_{t}h', f'HAR_peak_height_{t}h', f'HAR_local_background_{t}h', f'HAR_peak_focality_{t}h']: out[k] = ev[k]
    out[f'HAR_frame_fraction_{t}h'] = ev[f'HAR{t}_frame_fraction']
    out[f'HAR_frac_ORF1a_frame_{t}h'] = ev[f'HAR{t}_frac_ORF1a']; out[f'HAR_frac_third_frame_{t}h'] = ev[f'HAR{t}_frac_third']
    out[f'CHXwindow_frame_fraction_{t}h'] = ev[f'CHXwin{t}_frac_candidate']; out[f'CHXwindow_frac_ORF1a_frame_{t}h'] = ev[f'CHXwin{t}_frac_ORF1a']; out[f'CHXwindow_frac_third_frame_{t}h'] = ev[f'CHXwin{t}_frac_third']
    out[f'HAR_CHX_ratio_{t}h'] = ev[f'HAR_CHX_ratio_{t}h']
    for tag in ('_m15p30', '_m15p60'): out[f'HAR{t}_candidate_reads{tag}'] = ev[f'HAR{t}_candidate_reads{tag}']; out[f'HAR{t}_frame_fraction{tag}'] = ev[f'HAR{t}_frame_fraction{tag}']
for t in TPS:
    for k in [f'CHX_body_reads_{t}h', f'CHX_body_frac_candidate_{t}h', f'CHX_body_frac_ORF1a_{t}h', f'CHX_body_frac_third_{t}h', f'occupied_codons_{t}h', f'evaluable_codons_{t}h', f'occupancy_fraction_{t}h',
              f'longest_inframe_run_{t}h', f'first_occupied_codon_{t}h', f'last_occupied_codon_{t}h', f'occupied_span_nt_{t}h', f'top1_fraction_{t}h', f'top3_fraction_{t}h', f'de_top1_occupancy_{t}h', f'de_top3_occupancy_{t}h',
              f'de_top1_reads_{t}h', f'CHX_body_frac_from_hotspots_{t}h', f'occupied_codons_excl_hotspot_codons_{t}h', f'body_top1_pos_{t}h']:
        out[k] = ev[k]
for t in TPS:
    for k in [f'peak_attribution_class_{t}h', f'distance_to_candidate_start_{t}h', f'distance_to_ORF1a_start_{t}h', f'distance_to_nearest_other_candidate_start_{t}h', f'n_candidates_with_same_strongest_peak_{t}h',
              f'HAR_cand_frac_from_hotspots_{t}h', f'HAR_masked_peak_pos_{t}h', f'HAR_masked_peak_height_{t}h', f'HAR_masked_peak_attribution_{t}h', f'HAR_anyframe_peak_pos_{t}h', f'HAR_anyframe_peak_class_{t}h']:
        out[k] = ev[k]
for k in ['nearest_major_hotspot', 'distance_to_hotspot', 'hotspots_in_primary_window', 'canonical_start_confounded', 'five_prime_end_confounded', 'recurrent_hotspot_confounded', 'attribution_note',
          'HAR_peak_positions_5_8_16', 'HAR_peak_position_stable_within3nt_n_tp', 'ORF_family_members', 'descriptive_category', 'category_rationale']:
    out[k] = ev[k]
for c in ev.columns:
    if c.startswith('exploratory_'): out[c] = ev[c]
hdr = ("# PEDV_gRNA_5prime_candidate_evidence v1.0 | genomic RNA only (continuous reference; no sgRNA reference dependence) | ORF1a AUG verified from annotation = 293-295 | P-sites: RGP_v1.3 strict_core, + strand | "
       "primary window start-15..+45, sensitivity -15/+30 and -15/+60 | continuous metrics only: no composite score, no translated/not-translated call | descriptive_category A-D = manual summary, not a score | "
       "hotspot definition: hotspot_definition.json | exploratory_prev_* = previous round fixed thresholds, retained as descriptive variables only | 5/8/16 h are distinct contexts, not replicates\n")
with open(OUT + 'PEDV_gRNA_5prime_candidate_evidence.tsv', 'w', encoding='utf-8') as f:
    f.write(hdr); out.to_csv(f, sep='\t', index=False, float_format='%.4g')
print(out.shape, out.descriptive_category.value_counts().to_dict())
