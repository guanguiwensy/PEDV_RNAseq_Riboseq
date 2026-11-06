"""Step 1: transcript model + terminal ORF enumeration (sequence only, no translation evidence).
Transcripts: gRNA, S, ORF3, E, M, N. Junctions from frozen GIb_junction_remeasured.tsv
(donor = last leader base, acceptor = first body base)."""
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
import pandas as pd, numpy as np, hashlib
B = PEDV_P+''
OUT = B + 'fig5_terminal_ORF/'
seq = ''.join(l.strip() for l in open(B + 'refbuild_GIb/v1.1/virus/PEDV_GIb.fa') if not l.startswith('>')).upper()
G = len(seq)
jn = pd.read_csv(OUT + 'GIb_junction_remeasured.tsv', sep='\t', comment='#').set_index('sgRNA')
MAIN = {'gRNA': ('ORF1a', 293, 12622), 'S': ('S', 20610, 24758), 'ORF3': ('ORF3', 24758, 25033),
        'E': ('E', 25361, 25600), 'M': ('M', 25608, 26288), 'N': ('N', 26300, 27625)}
for g, (n, s, e) in MAIN.items():
    assert seq[s - 1:s + 2] == 'ATG' and seq[e - 3:e] in ('TAA', 'TAG', 'TGA'), (g, seq[s - 1:s + 2], seq[e - 3:e])
START = {'ATG', 'CTG', 'GTG', 'TTG', 'ATA', 'ATT', 'ATC', 'ACG'}
STOP = {'TAA', 'TAG', 'TGA'}
FILL = 9    # task convention: genome 1-9 treated as filled / unmeasured
TAIL = 60   # window (nt) at the end of N CDS in which overlapping 3' ORFs are still enumerated


def tx(name, variant=None):
    if name == 'gRNA':
        return seq, list(range(1, G + 1)), 0
    d, a = int(jn.loc[name, 'leader_end_donor_left_aligned']), int(jn.loc[name, 'body_start_acceptor_left_aligned'])
    if variant:
        d, a = variant
    return seq[:d] + seq[a - 1:], list(range(1, d + 1)) + list(range(a, G + 1)), d


def enumerate_orfs(name, variant=None):
    s, m, donor = tx(name, variant)
    L = len(s)
    _, ms, me = MAIN[name]
    mst, mend = m.index(ms) + 1, m.index(me) + 1
    nst = m.index(MAIN['N'][1]) + 1
    nend = m.index(MAIN['N'][2]) + 1
    out = []
    for i in range(L - 2):
        cod = s[i:i + 3]
        if cod not in START:
            continue
        t0 = i + 1
        if t0 < mst:
            side = '5'
        elif t0 >= nend - TAIL + 1:
            side = '3'
        else:
            continue
        j = i
        stop = None
        while j + 3 <= L:
            if s[j:j + 3] in STOP:
                stop = j
                break
            j += 3
        if stop is None:
            out.append(dict(unterminated=True, side=side, t0=t0))
            continue
        t1 = stop + 3
        naa = (stop - i) // 3
        frame = (t0 - mst) % 3
        if side == '5':
            if frame == 0 and t1 == mend:
                cls = 'upstream_inframe_extension'
            else:
                cls = 'uORF' if t1 < mst else 'uoORF'
            ov = max(0, min(t1, mend) - mst + 1)
        else:
            if t0 <= nend and (t0 - nst) % 3 == 0:
                continue          # in-frame inside N: N truncation, not a terminal ORF
            cls = '3prime_terminal_ORF'
            ov = max(0, nend - t0 + 1)
        out.append(dict(unterminated=False, side=side, ORF_class=cls, t0=t0, t1=t1, naa=naa, codon=cod,
                        frame=frame, ov=ov, nt=s[i:t1], g0=m[t0 - 1], g1=m[t1 - 1], mst=mst, mend=mend,
                        donor=donor, nend=nend))
    return out, s, m, mst, mend, nend


# alternative M junction variants (sensitivity)
alt = {}
for v in [(63, 25594), (67, 25598)]:
    o = enumerate_orfs('M', v)[0]
    alt.setdefault('M', []).append({(r['t0'], r['nt']) for r in o if not r['unterminated']})

rows, unterm = [], []
for name in MAIN:
    o, s, m, mst, mend, nend = enumerate_orfs(name)
    aug_pos = [i + 1 for i in range(mst - 1) if s[i:i + 3] == 'ATG']
    any_pos = [i + 1 for i in range(mst - 1) if s[i:i + 3] in START]
    jstat = '坐标实测;TRS-B身份待确认' if name == 'ORF3' else ('n/a(gRNA无junction)' if name == 'gRNA' else '实测')
    jtxt = 'NA' if name == 'gRNA' else f"{int(jn.loc[name, 'leader_end_donor_left_aligned'])}|{int(jn.loc[name, 'body_start_acceptor_left_aligned'])}"
    for r in o:
        if r['unterminated']:
            unterm.append((name, r['side']))
            continue
        t0 = r['t0']
        if r['side'] == '5':
            n_aug = sum(1 for p in aug_pos if p < t0)
            n_any = sum(1 for p in any_pos if p < t0)
            prox = n_aug == 0
            dms, dme = r['mst'] - t0, np.nan
        else:
            n_aug = n_any = np.nan
            prox = name == 'N'
            dms, dme = np.nan, t0 - nend - 1
        crosses = name != 'gRNA' and t0 <= r['donor'] < r['t1']
        lead = (t0 <= r['donor']) if name != 'gRNA' else (t0 <= 72)
        msens = ''
        if name == 'M':
            msens = any((t0, r['nt']) not in ak for ak in alt['M'])
        rows.append(dict(
            transcript_id=name, transcript_class='genomic' if name == 'gRNA' else 'sgRNA',
            canonical_main_ORF=MAIN[name][0], leader_body_junction=jtxt, junction_status=jstat,
            ORF_class=r['ORF_class'], region='5prime' if r['side'] == '5' else '3prime',
            genomic_start=r['g0'], genomic_stop=r['g1'], transcript_start=t0, transcript_stop=r['t1'],
            start_codon=r['codon'].replace('T', 'U'), start_type='AUG' if r['codon'] == 'ATG' else 'near-cognate',
            length_nt=r['t1'] - t0 + 1, length_aa=r['naa'], relative_frame=r['frame'],
            distance_to_main_start=dms, distance_to_main_stop=dme,
            overlaps_main_ORF=r['ov'] > 0, overlap_length=r['ov'],
            is_5prime_proximal_on_this_transcript=prox, n_upstream_AUG_on_transcript=n_aug,
            n_upstream_anyStart_on_transcript=n_any, leader_derived=bool(lead),
            crosses_leader_body_junction=bool(crosses),
            junction_proximal=bool(abs(t0 - r['donor']) <= 15) if name != 'gRNA' else False,
            start_in_filled_region_1_9=bool(r['g0'] <= FILL), M_variant_sensitive=msens, ORF_nt_seq=r['nt']))
df = pd.DataFrame(rows)
df = df[df.length_aa >= 8].copy()
df['length_bin_aa'] = pd.cut(df.length_aa, [7, 19, 39, 79, 10 ** 6], labels=['8-19', '20-39', '40-79', '>=80'])
df['locus_key'] = [hashlib.md5(f"{a}|{b}|{c}".encode()).hexdigest()[:8]
                   for a, b, c in zip(df.genomic_start, df.ORF_nt_seq, df.region)]
df['n_transcripts_identical_ORF'] = df.groupby('locus_key').transcript_id.transform('nunique')
order = list(MAIN)
df['transcripts_identical_ORF'] = df.groupby('locus_key').transcript_id.transform(
    lambda x: ','.join(sorted(set(x), key=order.index)))
df['universal_on_all_6'] = df.n_transcripts_identical_ORF == 6
df['needs_confirmation'] = df.start_in_filled_region_1_9
df['confirmation_reason'] = np.where(df.start_in_filled_region_1_9,
                                     'start in reference position 1-9 (filled; true gRNA 5-prime end unknown)',
                                     np.where(df.junction_status.str.contains('待确认'), 'ORF3 TRS-B identity pending', ''))
df['transcript_id'] = pd.Categorical(df.transcript_id, order)
df = df.sort_values(['region', 'transcript_id', 'transcript_start']).reset_index(drop=True)
df['transcript_id'] = df.transcript_id.astype(str)
df['candidate_id'] = df.transcript_id + '_' + df.region.str[0] + 'p_' + (df.groupby(['transcript_id', 'region']).cumcount() + 1).map('{:03d}'.format)
df['start_site_key'] = df.region + ':' + df.genomic_start.astype(str)
df['start_site_n_transcripts'] = df.groupby('start_site_key').transcript_id.transform('nunique')
df['start_site_transcripts'] = df.groupby('start_site_key').transcript_id.transform(lambda x: ','.join(sorted(set(x), key=order.index)))
# sub-threshold AUG uORFs (<8 aa), for transparency only
sub = []
for name in MAIN:
    s_, m_, d_ = tx(name)
    mst_ = m_.index(MAIN[name][1]) + 1
    for i in range(mst_ - 1):
        if s_[i:i + 3] == 'ATG':
            j = i
            while j + 3 <= len(s_) and s_[j:j + 3] not in STOP:
                j += 3
            sub.append(dict(transcript_id=name, genomic_start=m_[i], length_aa=(j - i) // 3, stop_found=j + 3 <= len(s_)))
pd.DataFrame(sub).to_csv(OUT + '_work/upstream_AUG_all_lengths.tsv', sep='	', index=False)
print(pd.DataFrame(sub))
u5 = df[df.region == '5prime']
print('5prime distinct start sites', u5.start_site_key.nunique(), ' shared by >=2 transcripts', (u5.drop_duplicates('start_site_key').start_site_n_transcripts >= 2).sum())
df.to_csv(OUT + '_work/catalog_sequence_only.tsv', sep='\t', index=False)

pd.set_option('display.width', 220, 'display.max_columns', 30)
print('unterminated excluded:', pd.Series(unterm).value_counts().to_dict() if unterm else 0)
print('rows', len(df), ' unique loci', df.locus_key.nunique())
print(df.pivot_table(index='transcript_id', columns='ORF_class', values='candidate_id', aggfunc='count', fill_value=0, observed=False).reindex(order))
print(df.pivot_table(index='transcript_id', columns='start_type', values='candidate_id', aggfunc='count', fill_value=0, observed=False).reindex(order))
print(df.groupby('ORF_class').length_aa.describe())
print(pd.crosstab(df.ORF_class, df.length_bin_aa))
u = df.drop_duplicates('locus_key')
print('UNIQUE LOCI'); print(pd.crosstab(u.ORF_class, u.start_type)); print(pd.crosstab(u.ORF_class, u.length_bin_aa))
print('filled1-9', df.start_in_filled_region_1_9.sum(), ' 5p-proximal', df.groupby('transcript_id').is_5prime_proximal_on_this_transcript.sum().to_dict())
print('M variant sens', (df.M_variant_sensitive == True).sum(), ' universal rows', df.universal_on_all_6.sum(), ' cross junction', df.crosses_leader_body_junction.sum())
print(df.groupby(['region', 'transcript_id'], observed=True).size())
