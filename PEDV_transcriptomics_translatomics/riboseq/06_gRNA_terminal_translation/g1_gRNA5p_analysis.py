"""gRNA 5' terminal translation: candidate catalogue + continuous evidence + peak attribution (no sgRNA, no binary calls).
Reference: PEDV_GIb refdata_GIb_v1.1 (genomic RNA, continuous). Reads: RGP_v1.3 strict_core P-sites, + strand."""
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
import os, re, json, numpy as np, pandas as pd
B = PEDV_P+''
OUT = B + 'fig5_terminal_ORF/gRNA_5prime/'
os.makedirs(OUT, exist_ok=True)
RGP = B + 'riboseq_GIb_mapping_psite/'
seq = ''.join(l.strip() for l in open(B + 'refbuild_GIb/v1.1/virus/PEDV_GIb.fa') if not l.startswith('>')).upper()
G = len(seq)
# ---- verify ORF1a start from annotation (not hard-coded) ----
gtf = [l.split('\t') for l in open(B + 'refbuild_GIb/v1.1/virus/PEDV_GIb.annotation.gtf') if not l.startswith('#')]
sc = [g for g in gtf if g[2] == 'start_codon' and 'ORF1a"' in g[8]][0]
cds = [g for g in gtf if g[2] == 'CDS' and 'ORF1a"' in g[8]][0]
MST = int(sc[3]); MEND = int(cds[4])
assert int(sc[4]) == MST + 2 and seq[MST - 1:MST + 2] == 'ATG', (sc[3], seq[MST - 1:MST + 2])
ANN = dict(ORF1a_start_codon=f'{MST}-{MST + 2}', ORF1a_CDS=f'{MST}-{MEND}', ORF1a_start_codon_seq=seq[MST - 1:MST + 2],
           analysis_window=f'1-{MST + 150}')
WIN_END = MST + 150
print(ANN)
START = {'ATG', 'CTG', 'GTG', 'TTG', 'ATA', 'ATT', 'ATC', 'ACG'}
STOP = {'TAA', 'TAG', 'TGA'}
CT = {}
b = 'TCAG'
aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
for i, a in enumerate(b):
    for j, bb in enumerate(b):
        for k, c in enumerate(b):
            CT[a + bb + c] = aas[16 * i + 4 * j + k]


def tr(s): return ''.join(CT[s[i:i + 3]] for i in range(0, len(s) - 2, 3))


# ---------------- catalogue ----------------
rows = []
for i in range(MST - 1):
    cod = seq[i:i + 3]
    if cod not in START: continue
    j = i
    while j + 3 <= G and seq[j:j + 3] not in STOP: j += 3
    st = i + 1; stop_first = j + 1; stop_end = j + 3
    naa = (j - i) // 3
    frame = (st - MST) % 3
    if frame == 0 and stop_end == MEND: cls = 'upstream_inframe_extension'
    else: cls = 'uORF' if stop_end < MST else 'uoORF'
    rows.append(dict(start_position=st, start_codon=cod.replace('T', 'U'), start_type='AUG' if cod == 'ATG' else 'near-cognate',
                     stop_codon_start=stop_first, stop_position=stop_end, stop_codon=seq[j:j + 3].replace('T', 'U'), length_aa=naa,
                     relative_frame_to_ORF1a=frame, ORF_class=cls, distance_to_ORF1a_AUG=MST - st,
                     overlap_with_ORF1a=max(0, min(stop_end, MEND) - MST + 1), nt_seq=seq[i:stop_end], aa_seq=tr(seq[i:j])))
allorf = pd.DataFrame(rows)
special = allorf[(allorf.start_position == 22)].copy()
cat = allorf[allorf.length_aa >= 8].copy().sort_values('start_position').reset_index(drop=True)
cat['candidate_id'] = ['gRNA_5p_%03d' % (i + 1) for i in range(len(cat))]
cat['length_nt'] = cat.stop_position - cat.start_position + 1
cat['ORF_family_stop'] = cat.stop_position
fam = cat.groupby('stop_position').candidate_id.transform(lambda x: ','.join(x))
cat['nested_starts_in_same_ORF'] = fam
cat['n_starts_in_same_ORF'] = cat.groupby('stop_position').candidate_id.transform('count')
cat['sequence_aa'] = np.where(cat.length_aa > 60, cat.aa_seq.str[:60] + '...(+%d aa to ORF1a stop)' % 0, cat.aa_seq)
cat.loc[cat.length_aa > 60, 'sequence_aa'] = cat.aa_seq.str[:60] + '...[' + cat.length_aa.astype(str) + ' aa total]'
cat['sequence_nt'] = np.where(cat.length_nt > 200, cat.nt_seq.str[:90] + '...', cat.nt_seq)
print('catalogue', len(cat), cat.ORF_class.value_counts().to_dict(), 'ORF families (distinct stops):', cat.stop_position.nunique())
cat_cols = ['candidate_id', 'start_position', 'start_codon', 'start_type', 'stop_position', 'stop_codon', 'length_aa', 'relative_frame_to_ORF1a', 'ORF_class',
            'distance_to_ORF1a_AUG', 'overlap_with_ORF1a', 'n_starts_in_same_ORF', 'nested_starts_in_same_ORF', 'sequence_aa', 'sequence_nt']
cat[cat_cols].to_csv(OUT + 'gRNA_5prime_candidate_catalogue.tsv', sep='\t', index=False)
special[['start_position', 'start_codon', 'stop_position', 'stop_codon', 'length_aa', 'relative_frame_to_ORF1a', 'ORF_class', 'distance_to_ORF1a_AUG', 'aa_seq', 'nt_seq']].assign(
    note='special unresolved 5-prime-end feature; below the 8 aa threshold, NOT part of the formal catalogue').to_csv(OUT + 'special_feature_AUG22_3aa.tsv', sep='\t', index=False)
print(special[['start_position', 'start_codon', 'length_aa', 'aa_seq', 'nt_seq']].to_string())

# ---------------- P-site data ----------------
SAMP = {5: (SR.lib("RIBO_VERO_CHX_1MOI_5h"), SR.lib("RIBO_VERO_HAR_1MOI_5h")), 8: (SR.lib("RIBO_VERO_CHX_1MOI_8h"), SR.lib("RIBO_VERO_HAR_1MOI_8h")), 16: (SR.lib("RIBO_VERO_CHX_1MOI_16h"), SR.lib("RIBO_VERO_HAR_1MOI_16h"))}
TPS = [5, 8, 16]
lo = pd.read_csv(RGP + 'psite_counts_PEDV_GIb.tsv.gz', sep='\t', comment='#', low_memory=False)
lo = lo[lo.strand == '+'].copy()
lo['read_length'] = pd.to_numeric(lo.read_length, errors='coerce'); lo = lo.dropna(subset=['read_length']); lo['read_length'] = lo.read_length.astype(int)
rls = pd.read_csv(RGP + 'readlength_sets_GIb.tsv', sep='\t', comment='#')
rls['read_length'] = pd.to_numeric(rls.read_length, errors='coerce'); rls = rls.dropna(subset=['read_length']); rls['read_length'] = rls.read_length.astype(int)
core = rls[rls.in_strict_core.astype(str) == 'True'][['sample_id', 'read_length']]
sc_ = lo.merge(core, on=['sample_id', 'read_length'])
P = {}
for s, g in sc_.groupby('sample_id'):
    a = np.zeros(G + 2); np.add.at(a, g.position.values, g.psite_count.values); P[s] = a
TOT = {s: P[s].sum() for s in P}
man = pd.read_csv(B + 'GIb_rnaseq_reanalysis/sample_manifest_GIb.tsv', sep='\t').set_index('sample_id')
man = man[man.index.isin(SR.WHITELIST_LIBRARIES)]; SR.assert_exactly(man.index, SR.RNA_VERO_ALIASES, 'gRNA RNA depth libraries')   # final sample table only
RNA = {}
for t in TPS:
    arr = np.zeros(G + 2); ids = man[man.time_h == t].index
    for i in ids:
        d = pd.read_csv(B + f'GIb_rnaseq_reanalysis/coverage/{i}.depth.tsv', sep='\t', header=None); arr[d[1].values] += d[2].values
    RNA[t] = arr / len(ids)
RNA_ALL = sum(RNA[t] for t in TPS) / 3

# ---------------- hotspots (data-driven rule + the user-named 103) ----------------
reg = np.arange(1, WIN_END + 1)
share = {k: P[k][reg] / P[k][reg].sum() for k in P}
rule_pos = sorted(set(int(p) for k in share for p in reg[share[k] >= 0.02]))
HOT = {}
for p in rule_pos + [103]:
    if p in (21, 22): g = '5prime_end_pileup'
    elif 285 <= p <= 296: g = 'ORF1a_start_region'
    elif p in (112, 113): g = 'hotspot_112_113'
    elif p == 103: g = 'hotspot_103'
    else: g = f'hotspot_{p}'
    HOT[p] = g
for p in range(285, 297): HOT[p] = 'ORF1a_start_region'
HOTCORE = dict(HOT)
for p0 in list(HOT):                       # +-1 nt shoulders (offset assignment is only good to ~1 nt; pile reads with 5' end 9-11 land on P-sites 21-23)
    if HOT[p0] == 'ORF1a_start_region': continue
    for q in (p0 - 1, p0 + 1):
        HOT.setdefault(q, HOT[p0])
GROUPS = sorted(set(HOT.values()), key=lambda g: min(p for p in HOT if HOT[p] == g))
print(f'hotspot positions (>=2pct of 1-{WIN_END} P-sites in >=1 library, plus 103):', rule_pos, GROUPS)
json.dump({'rule': f'>=2 percent of strict_core P-sites in 1-{WIN_END} in at least one of 6 libraries; 103 added by request; ORF1a start region fixed as 285-296 (-8..+3); other groups carry +-1 nt shoulders',
           'core_positions': {int(k): v for k, v in HOTCORE.items()}, 'positions_with_shoulders': {int(k): v for k, v in HOT.items()}}, open(OUT + 'hotspot_definition.json', 'w'), ensure_ascii=False, indent=1)
CLSMAP = {'5prime_end_pileup': '5prime_end_pileup', 'ORF1a_start_region': 'ORF1a_canonical_start'}


def hot_class(p):
    g = HOT.get(int(p))
    if g is None: return None
    return CLSMAP.get(g, 'recurrent_shared_hotspot')


starts = cat.start_position.values; stops = cat.stop_position.values
famstop = dict(zip(cat.candidate_id, cat.stop_position))


def attribute(cid, st, peak, famrows):
    """returns class, ..."""
    if peak is None or np.isnan(peak): return 'no_reads'
    peak = int(peak)
    hc = hot_class(peak)
    off = peak - st
    in_own = -3 <= off <= 30
    if hc: return hc + ('|also_in_own_window' if in_own else '')
    if in_own:
        closer = [s2 for s2 in famrows if s2 > st and peak - s2 >= -3]
        return 'nested_family_start' if closer else 'own_start'
    for s2 in starts:
        if s2 != st and -3 <= peak - s2 <= 30 and (peak - s2) % 3 == 0: return 'other_candidate_start'
    return 'unattributed'


def fr_fracs(vec, ts, st):
    rel = (ts - st) % 3
    ar = (MST - st) % 3
    x = vec[ts]; tot = x.sum()
    f = [x[rel == k].sum() / tot if tot else np.nan for k in range(3)]
    fc = f[0]
    if ar == 0: return fc, fc, np.nan, f, ar      # candidate frame == ORF1a frame: ORF1a value = candidate value, no distinct third frame
    third = [k for k in range(3) if k not in (0, ar)][0]
    return fc, f[ar], f[third], f, ar


rows = []; contrib = []
for _, r in cat.iterrows():
    st = int(r.start_position); cid = r.candidate_id
    famrows = [int(x) for x in cat[cat.stop_position == r.stop_position].start_position]
    o = dict(candidate_id=cid, start_position=st, start_codon=r.start_codon, stop_position=int(r.stop_position), stop_codon=r.stop_codon,
             length_aa=int(r.length_aa), ORF_class=r.ORF_class, relative_frame_to_ORF1a=int(r.relative_frame_to_ORF1a),
             distance_to_ORF1a_start=int(r.distance_to_ORF1a_AUG), start_type=r.start_type, n_starts_in_same_ORF=int(r.n_starts_in_same_ORF))
    # nearest hotspot to candidate start
    hp = np.array(sorted(HOTCORE)); dd = hp - st
    nh = hp[np.argmin(np.abs(dd))]
    o['nearest_major_hotspot'] = HOTCORE[nh] + f'@{nh}'; o['distance_to_hotspot'] = int(nh - st)
    others = [s for s in starts if s != st]
    o['distance_to_nearest_other_candidate_start_nt'] = int(min(abs(s - st) for s in others))
    ts = np.arange(max(1, st - 15), min(G, st + 45) + 1)
    o['hotspots_in_primary_window'] = ';'.join(f'{p}' for p in sorted(HOT) if ts[0] <= p <= ts[-1])
    conf = {'canonical': [], 'fiveend': [], 'hot': []}
    for t in TPS:
        chx, har = SAMP[t]
        # ---------------- HAR windows ----------------
        for tag, (a, b2) in {'': (-15, 45), '_m15p30': (-15, 30), '_m15p60': (-15, 60)}.items():
            w = np.arange(max(1, st + a), min(G, st + b2) + 1)
            rel = (w - st) % 3
            x = P[har][w]
            o[f'HAR{t}_candidate_reads{tag}'] = int(x[rel == 0].sum()); o[f'HAR{t}_total_reads{tag}'] = int(x.sum())
            o[f'HAR{t}_frame_fraction{tag}'] = x[rel == 0].sum() / x.sum() if x.sum() else np.nan
        w = ts; rel = (w - st) % 3; x = P[har][w]; cf = rel == 0
        hc = int(x[cf].sum()); cc = int(P[chx][w][cf].sum())
        o[f'CHXwin{t}_candidate_reads'] = cc
        o[f'HAR_CHX_ratio_{t}h'] = ((hc + 0.5) / TOT[har]) / ((cc + 0.5) / TOT[chx])
        for lab, vec in (('HAR', P[har]), ('CHXwin', P[chx])):
            fc, fa, f3, fl, ar = fr_fracs(vec, w, st)
            o[f'{lab}{t}_frac_candidate'] = fc; o[f'{lab}{t}_frac_ORF1a'] = fa; o[f'{lab}{t}_frac_third'] = f3
        # cand-frame peak
        vv = x * cf
        if vv.max() > 0:
            idx = np.where(vv == vv.max())[0]; k = idx[np.argmin(np.abs(w[idx] - st))]
            pk = int(w[k]); ph = int(vv.max())
            bg = (vv.sum() - ph) / max(1, cf.sum() - 1)
            o[f'HAR_peak_pos_{t}h'] = pk; o[f'HAR_peak_offset_{t}h'] = pk - st; o[f'HAR_peak_height_{t}h'] = ph
            o[f'HAR_local_background_{t}h'] = bg; o[f'HAR_peak_over_background_{t}h'] = ph / bg if bg > 0 else np.inf
            o[f'HAR_peak_focality_{t}h'] = ph / vv.sum()
        else:
            pk = np.nan; o[f'HAR_peak_pos_{t}h'] = np.nan; o[f'HAR_peak_offset_{t}h'] = np.nan; o[f'HAR_peak_height_{t}h'] = 0
            o[f'HAR_local_background_{t}h'] = 0.0; o[f'HAR_peak_over_background_{t}h'] = np.nan; o[f'HAR_peak_focality_{t}h'] = np.nan
        att = attribute(cid, st, pk, famrows)
        o[f'peak_attribution_class_{t}h'] = att
        o[f'distance_to_candidate_start_{t}h'] = (pk - st) if not np.isnan(pk) else np.nan
        o[f'distance_to_ORF1a_start_{t}h'] = (pk - MST) if not np.isnan(pk) else np.nan
        o[f'distance_to_nearest_other_candidate_start_{t}h'] = int(min(abs(pk - s) for s in others)) if not np.isnan(pk) else np.nan
        # strongest any-frame peak
        if x.max() > 0:
            ka = int(np.argmax(x)); pa = int(w[ka]); o[f'HAR_anyframe_peak_pos_{t}h'] = pa; o[f'HAR_anyframe_peak_height_{t}h'] = int(x.max())
            o[f'HAR_anyframe_peak_in_candidate_frame_{t}h'] = bool(cf[ka]); o[f'HAR_anyframe_peak_class_{t}h'] = hot_class(pa) or 'non_hotspot'
        else:
            o[f'HAR_anyframe_peak_pos_{t}h'] = np.nan; o[f'HAR_anyframe_peak_height_{t}h'] = 0; o[f'HAR_anyframe_peak_in_candidate_frame_{t}h'] = False; o[f'HAR_anyframe_peak_class_{t}h'] = 'no_reads'
        # hotspot contribution & masked residual (candidate frame)
        cfpos = w[cf]; cfv = x[cf]
        ismask = np.array([int(p) in HOT for p in cfpos])
        o[f'HAR_cand_frac_from_hotspots_{t}h'] = cfv[ismask].sum() / cfv.sum() if cfv.sum() else np.nan
        mv = cfv * (~ismask)
        o[f'HAR_cand_reads_hotspot_masked_{t}h'] = int(mv.sum())
        if mv.max() > 0:
            kk = int(np.argmax(mv)); mp = int(cfpos[kk]); o[f'HAR_masked_peak_pos_{t}h'] = mp; o[f'HAR_masked_peak_height_{t}h'] = int(mv.max())
            o[f'HAR_masked_peak_offset_{t}h'] = mp - st; o[f'HAR_masked_peak_attribution_{t}h'] = attribute(cid, st, mp, famrows)
            o[f'HAR_masked_peak_focality_{t}h'] = mv.max() / mv.sum()
        else:
            o[f'HAR_masked_peak_pos_{t}h'] = np.nan; o[f'HAR_masked_peak_height_{t}h'] = 0; o[f'HAR_masked_peak_offset_{t}h'] = np.nan
            o[f'HAR_masked_peak_attribution_{t}h'] = 'no_reads'; o[f'HAR_masked_peak_focality_{t}h'] = np.nan
        for g in GROUPS:
            gs = np.array([HOT.get(int(p)) == g for p in cfpos]); contrib.append(dict(candidate_id=cid, tp=t, assay='HAR_window_cand_frame', group=g, reads=int(cfv[gs].sum()), total=int(cfv.sum())))
        conf['canonical'].append(att.startswith('ORF1a_canonical_start')); conf['fiveend'].append(att.startswith('5prime_end_pileup')); conf['hot'].append(att.startswith('recurrent_shared_hotspot'))
        # ---------------- CHX body ----------------
        end = (MST - 1) if r.ORF_class == 'upstream_inframe_extension' else int(r.stop_position) - 3
        cod = np.arange(st, end - 1, 3); cod = cod[cod + 2 <= end]; nco = len(cod)
        bts = np.arange(st, end + 1); vC = P[chx]
        xb = vC[bts]; rel = (bts - st) % 3
        o[f'CHX_body_reads_{t}h'] = int(xb[rel == 0].sum()); o[f'CHX_body_total_reads_{t}h'] = int(xb.sum())
        fc, fa, f3, fl, ar = fr_fracs(vC, bts, st)
        o[f'CHX_body_frac_candidate_{t}h'] = fc; o[f'CHX_body_frac_ORF1a_{t}h'] = fa; o[f'CHX_body_frac_third_{t}h'] = f3
        cv = vC[cod] if nco else np.array([]); ev = (RNA[t][cod] >= 5) if nco else np.array([], bool); occ = cv > 0
        n_occ = int(occ.sum()); n_ev = int(ev.sum())
        o[f'occupied_codons_{t}h'] = n_occ; o[f'evaluable_codons_{t}h'] = n_ev; o[f'body_codons_{t}h'] = nco
        o[f'occupancy_fraction_{t}h'] = n_occ / n_ev if n_ev else np.nan
        run = br = 0
        for z in occ:
            run = run + 1 if z else 0; br = max(br, run)
        o[f'longest_inframe_run_{t}h'] = br
        if n_occ:
            ix = np.where(occ)[0]; so = np.sort(cv)[::-1]
            o[f'first_occupied_codon_{t}h'] = int(ix[0] + 1); o[f'last_occupied_codon_{t}h'] = int(ix[-1] + 1); o[f'occupied_span_nt_{t}h'] = int((ix[-1] - ix[0] + 1) * 3)
            o[f'top1_fraction_{t}h'] = so[0] / cv.sum(); o[f'top3_fraction_{t}h'] = so[:3].sum() / cv.sum()
            o[f'de_top1_occupancy_{t}h'] = (n_occ - 1) / max(1, n_ev - 1); o[f'de_top3_occupancy_{t}h'] = max(0, n_occ - 3) / max(1, n_ev - 3)
            o[f'de_top1_reads_{t}h'] = int(cv.sum() - so[0])
            o[f'body_top1_pos_{t}h'] = int(cod[int(np.argmax(cv))]); o[f'body_top1_attribution_{t}h'] = hot_class(cod[int(np.argmax(cv))]) or 'non_hotspot'
        else:
            for kx in ['first_occupied_codon', 'last_occupied_codon', 'occupied_span_nt', 'top1_fraction', 'top3_fraction', 'de_top1_occupancy', 'de_top3_occupancy', 'de_top1_reads', 'body_top1_pos', 'body_top1_attribution']: o[f'{kx}_{t}h'] = np.nan
        hm = np.array([int(p) in HOT for p in cod]) if nco else np.array([], bool)
        o[f'CHX_body_frac_from_hotspots_{t}h'] = cv[hm].sum() / cv.sum() if nco and cv.sum() else np.nan
        o[f'occupied_codons_excl_hotspot_codons_{t}h'] = int((occ & ~hm).sum()) if nco else 0
        o[f'CHX_body_reads_hotspot_masked_{t}h'] = int(cv[~hm].sum()) if nco else 0
        for g in GROUPS:
            gs = np.array([HOT.get(int(p)) == g for p in cod]) if nco else np.array([], bool)
            contrib.append(dict(candidate_id=cid, tp=t, assay='CHX_body_cand_frame', group=g, reads=int(cv[gs].sum()) if nco else 0, total=int(cv.sum()) if nco else 0))
        o[f'RNA_local_cov_{t}h'] = float(RNA[t][ts].mean())
    o['canonical_start_confounded'] = ';'.join(f'{t}h' for t, c in zip(TPS, conf['canonical']) if c) or 'no'
    o['five_prime_end_confounded'] = ';'.join(f'{t}h' for t, c in zip(TPS, conf['fiveend']) if c) or 'no'
    o['recurrent_hotspot_confounded'] = ';'.join(f'{t}h' for t, c in zip(TPS, conf['hot']) if c) or 'no'
    # time dimension descriptors
    pks = [o[f'HAR_peak_pos_{t}h'] for t in TPS]
    pk_ok = [p for p in pks if not np.isnan(p)]
    o['HAR_peak_positions_5_8_16'] = '/'.join('NA' if np.isnan(p) else str(int(p)) for p in pks)
    o['HAR_peak_position_stable_within3nt_n_tp'] = max([sum(abs(q - p) <= 3 for q in pk_ok) for p in pk_ok]) if pk_ok else 0
    rows.append(o)
ev = pd.DataFrame(rows)
pd.DataFrame(contrib).to_csv(OUT + 'hotspot_contribution_long.tsv', sep='\t', index=False)

# exploratory legacy descriptors (previous round thresholds; NOT used for any classification)
leg = pd.read_csv(B + 'fig5_terminal_ORF/PEDV_terminal_ORF_candidates.tsv', sep='\t', comment='#', keep_default_na=False, na_values=[''])
leg = leg[(leg.transcript_id == 'gRNA') & (leg.region == '5prime')].set_index('genomic_start')
for t in TPS:
    for nm in ('HAR_start_signal', 'CHX_body_signal', 'CHX_body_signal_deTop1'):
        ev[f'exploratory_prev_{nm}_{t}h'] = [leg.loc[s, f'{nm}_{t}h'] if s in leg.index else np.nan for s in ev.start_position]
ev['exploratory_note'] = 'exploratory_prev_* = last round fixed thresholds (reads>=10, peak -3..+30, peak>=5, frac>=0.5, ratio>=2); descriptive only'
ev.to_pickle(OUT + '_ev.pkl')
cat.to_pickle(OUT + '_cat.pkl')
json.dump(dict(MST=MST, MEND=MEND, HOT={int(k): v for k, v in HOT.items()}, GROUPS=GROUPS, ANN=ANN), open(OUT + '_meta.json', 'w'))
print(ev.shape)
