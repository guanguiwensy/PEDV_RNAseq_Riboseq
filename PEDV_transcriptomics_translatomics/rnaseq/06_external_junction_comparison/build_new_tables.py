"""Tables for the two candidate-ORF figures under the reclassified scheme.
Counted classes: truncated (frame 'in'), frameshift (frame 'out'), out-of-frame ORF (no AUG in acceptor CDS / noncoding,
  but the sgRNA carries an ATG...stop ORF not in canonical frame). Excluded from statistics: translates_canonical.
ORF identity = start ATG (1-based). Truncated / frameshift start = first downstream AUG of the acceptor (as before).
Out-of-frame ORF start = first ATG at or after the acceptor in the sgRNA that is followed by an in-frame stop inside it.
"""
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
import sample_registry as SR
from pathlib import Path
import pandas as pd

BASE = Path(PEDV_P+'')
TAB = BASE / 'Noncanonical_reclass' / 'tables'
seq = ''.join(l.strip() for l in open(BASE / 'refbuild_GIb/virus/PEDV_GIb.fa') if not l.startswith('>')).upper()
STOP = {'TAA', 'TAG', 'TGA'}
CDS = {'ORF1a': (293, 12622), 'ORF1b': (12592, 20613), 'S': (20610, 24758), 'ORF3': (24758, 25033),
       'E': (25361, 25600), 'M': (25608, 26288), 'N': (26300, 27625)}
CLS_KEY = {'truncated (counted)': 'truncated', 'frameshift (counted)': 'frameshift',
           'out-of-frame ORF (counted)': 'out-of-frame ORF', 'translates_canonical (excluded)': 'excluded'}


def gene_of(p):
    hit = [k for k, (s, t) in CDS.items() if s <= p <= t]
    return max(hit, key=lambda k: CDS[k][0]) if hit else 'intergenic/UTR'


def first_atg(a):
    i = seq.find('ATG', a - 1)
    return i + 1 if i >= 0 else None


def orf_start_oof(D, A):
    """first body ATG (genome position >= A) with an in-frame stop inside the molecule."""
    mol = seq[:D] + seq[A - 1:]
    for i in range(D, len(mol) - 2):          # body part of the molecule
        if mol[i:i + 3] != 'ATG':
            continue
        j = i + 3
        while j + 3 <= len(mol):
            if mol[j:j + 3] in STOP:
                return i - D + A              # genome position of the ATG
            j += 3
    return None


rec = pd.read_csv(TAB / 'TRSL_dep_noncanonical_reclassified.csv')
met = pd.read_csv(BASE / 'Noncanonical_TOP_sgRNA_TRS' / 'tables' / 'T_all_TRSL_dep_noncanonical_metrics.csv')
ev = rec.merge(met[['donor', 'acceptor', 'trs_identity_of_12', 'microhomology_nt', 'mfe_b_kcal', 'mfe_a_kcal']],
               on=['donor', 'acceptor'], how='left')
ev['cls'] = ev.new_class.map(CLS_KEY)
ev['orf_start'] = None
for idx, r in ev.iterrows():
    if r.cls in ('truncated', 'frameshift'):
        ev.at[idx, 'orf_start'] = first_atg(int(r.acceptor))
    elif r.cls == 'out-of-frame ORF':
        ev.at[idx, 'orf_start'] = orf_start_oof(int(r.donor), int(r.acceptor))
ev['start_codon'] = ev.orf_start.map(lambda p: seq[int(p) - 1:int(p) + 2] if pd.notna(p) else None)
ev['start_gene'] = ev.orf_start.map(lambda p: gene_of(int(p)) if pd.notna(p) else None)
ev.to_csv(TAB / 'events_with_new_ORF_assignment.csv', index=False)

# per-library support
ps = pd.read_csv(BASE / 'GIb_rnaseq_reanalysis/figures/Fig2_work/per_sample_events.tsv.gz', sep='\t').rename(columns={'acc': 'acceptor'})
meta = pd.read_csv(BASE.parent / 'final_figure_reproducibility/00_sample_metadata/PEDV_sample_metadata_MOI_corrected.tsv', sep='\t')
meta = meta[meta.sample_id.between(SR.lib("RNA_VERO_CHX_0.1MOI_3h"), SR.lib("RNA_VERO_HAR_1MOI_16h"))].copy()
meta['MOI'] = meta.MOI.astype(float)
ps = ps.merge(meta[['sample_id', 'MOI', 'treatment', 'hpi']], left_on='sample', right_on='sample_id')
cnt = ev[ev.cls != 'excluded'].dropna(subset=['orf_start'])
ps = ps.merge(cnt[['donor', 'acceptor', 'cls', 'orf_start', 'start_codon', 'start_gene']], on=['donor', 'acceptor'])
ps['orf_start'] = ps.orf_start.astype(int)
ps.to_csv(TAB / 'per_library_support_counted.csv', index=False)

# ORF ranking (one row per start ATG within class)
rank = (ps.groupby(['cls', 'orf_start', 'start_codon', 'start_gene'])
        .agg(n_supporting_junctions=('donor', lambda s: len(set(zip(s, ps.loc[s.index, 'acceptor'])))),
             n_donors=('donor', 'nunique'), reads=('reads', 'sum'), n_libraries_detected=('sample', 'nunique'))
        .reset_index())
rank = rank.sort_values(['cls', 'n_supporting_junctions', 'reads'], ascending=[True, False, False])
rank['rank_within_class'] = rank.groupby('cls').cumcount() + 1
rank.to_csv(TAB / 'T1_new_ORF_ranking.csv', index=False)

# time course per ORF x MOI x treatment x hpi
ps['condition'] = 'MOI' + ps.MOI.map('{:g}'.format) + '_' + ps.treatment + '_' + ps.hpi.astype(str) + 'h'
tc = (ps.groupby(['cls', 'orf_start', 'start_codon', 'start_gene', 'MOI', 'treatment', 'hpi'])
      .agg(reads=('reads', 'sum'), n_supporting_junctions=('donor', 'size')).reset_index())
tc.to_csv(TAB / 'T3_new_time_course.csv', index=False)

# summary by class
summ = ev.groupby('new_class').agg(events=('reads', 'size'), reads=('reads', 'sum')).reset_index()
distinct = rank.groupby('cls').size().rename('distinct_ORFs').reset_index()
distinct.to_csv(TAB / 'T4_new_distinct_ORFs.csv', index=False)
summ.to_csv(TAB / 'T4_new_class_summary.csv', index=False)
print(summ.to_string(index=False))
print(distinct.to_string(index=False))
print(rank.groupby('cls').head(3).to_string(index=False))
