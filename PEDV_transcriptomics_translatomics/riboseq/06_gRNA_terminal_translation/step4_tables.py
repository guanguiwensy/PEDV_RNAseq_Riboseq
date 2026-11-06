"""Step 4: write the two deliverable tables (schema of TASK §10; per-timepoint columns kept separate, no composite score)."""
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
import pandas as pd, numpy as np
B = PEDV_P+''
OUT = B + 'fig5_terminal_ORF/'
d = pd.read_csv(OUT + '_work/candidates_with_evidence.tsv', sep='\t', keep_default_na=False, na_values=[''])
TPS = [5, 8, 16]
ren = {}
for t in TPS:
    ren[f'HAR_start_cand_{t}h'] = f'HAR{t}_start'
    ren[f'CHX_start_cand_{t}h'] = f'CHX{t}_start'
    ren[f'HAR_CHX_{t}h'] = f'HAR_CHX{t}'
d = d.rename(columns=ren)
d['sgRNA_junction_support'] = d.sgRNA_junction_support
d['transcript_abundance_proxy'] = d.transcript_abundance_proxy_junction_ppm
d['RNA_local_coverage_5h_8h_16h'] = d.apply(lambda r: ';'.join(f"{r[f'RNA_local_cov_{t}h']:.1f}" for t in TPS), axis=1)
d['distance_to_leader_body_junction'] = np.nan
junc = {'gRNA': np.nan}
jn = pd.read_csv(OUT + 'GIb_junction_remeasured.tsv', sep='\t', comment='#').set_index('sgRNA')
for k in jn.index:
    junc[k] = jn.loc[k, 'leader_end_donor_left_aligned']
d['distance_to_leader_body_junction'] = d.apply(lambda r: r.transcript_start - junc[r.transcript_id] if not np.isnan(junc[r.transcript_id]) else np.nan, axis=1)
d['specific_5prime_UTR_region'] = np.where((d.region == '5prime') & ~d.leader_derived, True, False)
first = ['transcript_id', 'transcript_class', 'canonical_main_ORF', 'leader_body_junction', 'junction_status', 'candidate_id', 'ORF_class', 'region',
         'genomic_start', 'genomic_stop', 'transcript_start', 'transcript_stop', 'start_codon', 'start_type', 'length_nt', 'length_aa', 'length_bin_aa',
         'relative_frame', 'distance_to_main_start', 'distance_to_main_stop', 'overlaps_main_ORF', 'overlap_length',
         'is_5prime_proximal_on_this_transcript', 'n_upstream_AUG_on_transcript', 'n_upstream_anyStart_on_transcript',
         'leader_derived', 'crosses_leader_body_junction', 'junction_proximal', 'distance_to_leader_body_junction', 'specific_5prime_UTR_region',
         'start_in_filled_region_1_9', 'needs_confirmation', 'confirmation_reason', 'M_variant_sensitive',
         'locus_key', 'n_transcripts_identical_ORF', 'transcripts_identical_ORF', 'start_site_n_transcripts']
har = ['HAR5_start', 'HAR8_start', 'HAR16_start', 'CHX5_start', 'CHX8_start', 'CHX16_start', 'HAR_CHX5', 'HAR_CHX8', 'HAR_CHX16',
       'strongest_HAR_tp', 'strongest_HAR_peak_position', 'peak_offset_from_start', 'n_HAR_positive_timepoints']
for t in TPS:
    har += [f'HAR_start_s30_cand_{t}h', f'HAR_start_s60_cand_{t}h', f'HAR_start_total_{t}h', f'HAR_start_main_{t}h', f'HAR_start_other_{t}h',
            f'HAR_cand_frame_frac_{t}h', f'HAR_peak_pos_{t}h', f'HAR_peak_off_{t}h', f'HAR_peak_n_{t}h', f'HAR_peak_is_main_start_{t}h',
            f'HAR_win_top1_frac_{t}h', f'HAR_win_top1_pos_{t}h', f'HAR_win_top3_frac_{t}h', f'HAR_win_resid_deTop1_{t}h', f'HAR_win_resid_deTop3_{t}h',
            f'HAR_win_f29_{t}h', f'HAR_win_f29_lib_{t}h', f'HAR_win_includes_pos21_22_{t}h', f'HAR_start_signal_{t}h']
chx = ['n_CHX_supported_timepoints', 'n_CHX_supported_timepoints_deTop1']
for t in TPS:
    chx += [f'CHX_body_cand_{t}h', f'CHX_body_main_{t}h', f'CHX_body_other_{t}h', f'CHX_body_cand_share_{t}h', f'occupied_codons_{t}h', f'evaluable_codons_{t}h',
            f'n_body_codons_{t}h', f'occupancy_fraction_{t}h', f'longest_inframe_run_{t}h', f'body_coverage_{t}h', f'first_occupied_codon_{t}h',
            f'last_occupied_codon_{t}h', f'occupied_span_nt_{t}h', f'top1_fraction_{t}h', f'top3_fraction_{t}h', f'de_top1_occupancy_{t}h',
            f'de_top3_occupancy_{t}h', f'body_top1_pos_{t}h', f'CHX_body_f29_{t}h', f'CHX_body_f29_lib_{t}h', f'CHX_body_cand_deTop1_{t}h',
            f'CHX_body_cand_share_deTop1_{t}h', f'CHX_body_signal_{t}h', f'CHX_body_signal_deTop1_{t}h']
rna = ['RNA_local_coverage_5h_8h_16h', 'sgRNA_junction_support', 'transcript_abundance_proxy']
cols = first + har + chx + rna + ['ORF_nt_seq']
d = d[cols].copy()
hdr = ("# PEDV_terminal_ORF_candidates v1.0 | refdata_GIb_v1.1 | junctions: GIb_junction_remeasured.tsv (frozen; ORF3 TRS-B identity pending) | P-sites: RGP_v1.3 strict_core, + strand | "
       "HAR/CHX window = start-15..+45 nt (transcript coords; _s30/_s60 = -15/+30, -15/+60) | no composite score, no binary translated call | "
       "descriptor flags *_signal_*: HAR = cand-frame reads>=10, peak in [-3,+30], peak n>=5, cand-frame fraction>=0.5, HAR/CHX ppm ratio>=2, peak not at canonical main start; "
       "CHX body = cand-frame reads>=10 and cand-frame share>=0.5 (deTop1 = after removing the top codon) | 5h/8h/16h are distinct contexts, not replicates; only 5h is a fully paired CHX/HAR preparation | "
       "HAR/CHX ratios at 8h/16h carry preparation confounding | 3prime rows repeat the same 33 loci on each transcript (identical P-site data); is_5prime_proximal differs\n")
with open(OUT + 'PEDV_terminal_ORF_candidates.tsv', 'w', encoding='utf-8') as f:
    f.write(hdr)
    d.to_csv(f, sep='\t', index=False, float_format='%.4g')
# ---------------- summary by transcript ----------------
order = ['gRNA', 'S', 'ORF3', 'E', 'M', 'N']
rows = []
for t in order:
    x = d[d.transcript_id == t]
    five = x[x.region == '5prime']; three = x[x.region == '3prime']
    rows.append(dict(transcript_id=t, canonical_main_ORF=x.canonical_main_ORF.iat[0], junction_status=x.junction_status.iat[0],
                     n_uORF=(five.ORF_class == 'uORF').sum(), n_uoORF=(five.ORF_class == 'uoORF').sum(),
                     n_upstream_extension=(five.ORF_class == 'upstream_inframe_extension').sum(), n_3prime_terminal_ORF=len(three),
                     n_AUG_candidates=(x.start_type == 'AUG').sum(), n_near_cognate_candidates=(x.start_type == 'near-cognate').sum(),
                     n_AUG_5prime=(five.start_type == 'AUG').sum(), n_AUG_3prime=(three.start_type == 'AUG').sum(),
                     median_length_aa=x.length_aa.median(), max_length_aa=x.length_aa.max(),
                     median_length_aa_5prime_excl_extension=five[five.ORF_class != 'upstream_inframe_extension'].length_aa.median(),
                     n_5prime_proximal=int(x.is_5prime_proximal_on_this_transcript.sum() if t != 'N' else five.is_5prime_proximal_on_this_transcript.sum()),
                     n_HAR_supported_5prime=(five.n_HAR_positive_timepoints >= 1).sum(), n_CHX_body_supported_5prime=(five.n_CHX_supported_timepoints >= 1).sum(),
                     n_both_5prime=((five.n_HAR_positive_timepoints >= 1) & (five.n_CHX_supported_timepoints >= 1)).sum(),
                     n_HAR_supported_3prime=(three.n_HAR_positive_timepoints >= 1).sum(), n_CHX_body_supported_3prime=(three.n_CHX_supported_timepoints >= 1).sum(),
                     n_both_3prime=((three.n_HAR_positive_timepoints >= 1) & (three.n_CHX_supported_timepoints >= 1)).sum(),
                     n_HAR_supported=(x.n_HAR_positive_timepoints >= 1).sum(), n_CHX_body_supported=(x.n_CHX_supported_timepoints >= 1).sum(),
                     n_start_in_filled_1_9=int(x.start_in_filled_region_1_9.sum())))
s = pd.DataFrame(rows)
with open(OUT + 'PEDV_terminal_ORF_summary_by_transcript.tsv', 'w', encoding='utf-8') as f:
    f.write('# PEDV_terminal_ORF_summary_by_transcript v1.0 | counts are candidate rows per transcript; 3prime columns repeat the same 33 loci on every transcript | n_*_supported use descriptor flags (>=1 timepoint), not binary translation calls\n')
    s.to_csv(f, sep='\t', index=False)
pd.set_option('display.width', 250, 'display.max_columns', 40)
print(s.T)
print(len(d), d.shape)
