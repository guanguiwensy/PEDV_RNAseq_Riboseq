"""Step 2: HAR start-window / CHX body / spike / RNA-context evidence for every catalogued candidate.
Inputs: _work/catalog_sequence_only.tsv, RGP_v1.3 P-site long table (strict_core = 29/30 nt per library set), RNA-seq depth.
No composite score, no binary translated call. 'signal' flags are descriptors with fixed thresholds written in THRESH below."""
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
import pandas as pd, numpy as np, glob, os
B = PEDV_P+''
OUT = B + 'fig5_terminal_ORF/'
RGP = B + 'riboseq_GIb_mapping_psite/'
seq = ''.join(l.strip() for l in open(B + 'refbuild_GIb/v1.1/virus/PEDV_GIb.fa') if not l.startswith('>')).upper()
G = len(seq)
SAMP = {5: (SR.lib("RIBO_VERO_CHX_1MOI_5h"), SR.lib("RIBO_VERO_HAR_1MOI_5h")), 8: (SR.lib("RIBO_VERO_CHX_1MOI_8h"), SR.lib("RIBO_VERO_HAR_1MOI_8h")), 16: (SR.lib("RIBO_VERO_CHX_1MOI_16h"), SR.lib("RIBO_VERO_HAR_1MOI_16h"))}  # (CHX, HAR)
TPS = [5, 8, 16]
THRESH = dict(har_min_reads=10, har_peak_lo=-3, har_peak_hi=30, har_peak_min=5, har_frame_frac=0.5, har_chx_ratio=2.0,
              chx_min_reads=10, chx_frame_share=0.5, rna_eval_depth=5)
BASE29 = {SR.lib("RIBO_VERO_CHX_1MOI_8h"): 0.67, SR.lib("RIBO_VERO_HAR_1MOI_8h"): 0.41, SR.lib("RIBO_VERO_HAR_1MOI_16h"): 0.57, SR.lib("RIBO_VERO_CHX_1MOI_5h"): 0.42}  # fractions quoted in task (29 nt share)

# ---------- load P-sites (strict_core, + strand) ----------
lo = pd.read_csv(RGP + 'psite_counts_PEDV_GIb.tsv.gz', sep='\t', comment='#', low_memory=False)
lo = lo[lo.strand == '+']
lo['read_length'] = pd.to_numeric(lo.read_length, errors='coerce')
lo = lo.dropna(subset=['read_length']); lo['read_length'] = lo.read_length.astype(int)
rls = pd.read_csv(RGP + 'readlength_sets_GIb.tsv', sep='\t', comment='#')
rls['read_length'] = pd.to_numeric(rls.read_length, errors='coerce')
rls = rls.dropna(subset=['read_length']); rls['read_length'] = rls.read_length.astype(int)
core = rls[rls.in_strict_core.astype(str) == 'True'][['sample_id', 'read_length']]
lo = lo.merge(core, on=['sample_id', 'read_length'])
P, P29, P30 = {}, {}, {}
for s, g in lo.groupby('sample_id'):
    a = np.zeros(G + 2); a29 = np.zeros(G + 2); a30 = np.zeros(G + 2)
    for rl, gg in g.groupby('read_length'):
        np.add.at(a, gg.position.values, gg.psite_count.values)
        if rl == 29: np.add.at(a29, gg.position.values, gg.psite_count.values)
        if rl == 30: np.add.at(a30, gg.position.values, gg.psite_count.values)
    P[s], P29[s], P30[s] = a, a29, a30
TOT = {s: P[s].sum() for s in P}
# cross-check with wide table
w = pd.read_csv(RGP + 'psite_counts_PEDV_GIb_wide.tsv.gz', sep='\t', comment='#')
w = w[w.strand == '+']
for s in P:
    assert abs(w[w.sample_id == s].strict_core.sum() - TOT[s]) < 1, s
LIB29 = {s: P29[s].sum() / (P29[s].sum() + P30[s].sum()) for s in P}
print('lib 29-nt fraction among strict_core:', {s: round(v, 3) for s, v in LIB29.items()})
CAN_SKEW = {SR.lib("RIBO_VERO_HAR_1MOI_5h"): False}  # single-length libs: no skew test

# ---------- RNA depth (pooled over libs of the time point; non-strand-specific, gRNA+all sgRNA summed) ----------
man = pd.read_csv(B + 'GIb_rnaseq_reanalysis/sample_manifest_GIb.tsv', sep='\t').set_index('sample_id')
man = man[man.index.isin(SR.WHITELIST_LIBRARIES)]; SR.assert_exactly(man.index, SR.RNA_VERO_ALIASES, 'gRNA RNA depth libraries')   # final sample table only
RNA = {}
for t in TPS:
    ids = man[man.time_h == t].index
    arr = np.zeros(G + 2)
    for i in ids:
        d = pd.read_csv(B + f'GIb_rnaseq_reanalysis/coverage/{i}.depth.tsv', sep='\t', header=None)
        arr[d[1].values] += d[2].values
    RNA[t] = arr / len(ids)
jn = pd.read_csv(OUT + 'GIb_junction_remeasured.tsv', sep='\t', comment='#').set_index('sgRNA')
tot_j = jn.reads_modal_pair.sum()


def tx(name):
    if name == 'gRNA':
        return seq, list(range(1, G + 1))
    d, a = int(jn.loc[name, 'leader_end_donor_left_aligned']), int(jn.loc[name, 'body_start_acceptor_left_aligned'])
    return seq[:d] + seq[a - 1:], list(range(1, d + 1)) + list(range(a, G + 1))


MAIN = {'gRNA': 293, 'S': 20610, 'ORF3': 24758, 'E': 25361, 'M': 25608, 'N': 26300}
cat = pd.read_csv(OUT + '_work/catalog_sequence_only.tsv', sep='\t', keep_default_na=False, na_values=[''])
cache = {}
for name in MAIN:
    s_, m_ = tx(name)
    mm = np.array([0] + m_)  # 1-based transcript -> genomic
    cache[name] = dict(m=mm, L=len(m_), mst=m_.index(MAIN[name]) + 1,
                       P={k: np.concatenate([[0], P[k][mm[1:]]]) for k in P},
                       P29={k: np.concatenate([[0], P29[k][mm[1:]]]) for k in P},
                       P30={k: np.concatenate([[0], P30[k][mm[1:]]]) for k in P},
                       R={t: np.concatenate([[0], RNA[t][mm[1:]]]) for t in TPS})


def frames(ts, t0, mst):
    c = (ts - t0) % 3 == 0
    mn = (ts - mst) % 3 == 0
    return c, mn


def spike_stats(vec29, vec30, v, ts, k):
    """signal-structure descriptors over positions ts (array of transcript positions) for sample k, all frames"""
    x = v[ts]
    tot = x.sum()
    if tot == 0:
        return dict(n=0, top1=np.nan, top1_pos=np.nan, top3=np.nan, resid1=np.nan, resid3=np.nan, f29=np.nan)
    o = np.argsort(-x)
    t1 = x[o[0]] / tot
    t3 = x[o[:3]].sum() / tot
    f29 = vec29[ts].sum() / max(1.0, vec29[ts].sum() + vec30[ts].sum())
    return dict(n=int(tot), top1=t1, top1_pos=int(ts[o[0]]), top3=t3, resid1=int(tot - x[o[0]]), resid3=int(tot - x[o[:3]].sum()), f29=f29)


rows = []
for _, r in cat.iterrows():
    c = cache[r.transcript_id]
    L, mst = c['L'], c['mst']
    t0, t1 = int(r.transcript_start), int(r.transcript_stop)
    o = dict(candidate_id=r.candidate_id)
    best = None
    for t in TPS:
        chx, har = SAMP[t]
        # ---- start window (HAR focus) ----
        for tag, (a, b) in {'': (-15, 45), 's30': (-15, 30), 's60': (-15, 60)}.items():
            ts = np.arange(max(1, t0 + a), min(L, t0 + b) + 1)
            cf, mf = frames(ts, t0, mst)
            for lab, k in (('HAR', har), ('CHX', chx)):
                v = c['P'][k]
                x = v[ts]
                o[f'{lab}_start{tag and "_" + tag}_cand_{t}h'] = int(x[cf].sum())
                if tag == '':
                    o[f'{lab}_start_total_{t}h'] = int(x.sum())
                    o[f'{lab}_start_main_{t}h'] = int(x[mf & ~cf].sum())
                    o[f'{lab}_start_other_{t}h'] = int(x[~cf & ~mf].sum())
        ts = np.arange(max(1, t0 - 15), min(L, t0 + 45) + 1)
        cf, mf = frames(ts, t0, mst)
        hc = o[f'HAR_start_cand_{t}h']; cc = o[f'CHX_start_cand_{t}h']
        o[f'HAR_CHX_{t}h'] = ((hc + 0.5) / TOT[har]) / ((cc + 0.5) / TOT[chx])
        o[f'HAR_cand_ppm_{t}h'] = hc / TOT[har] * 1e6
        o[f'CHX_cand_ppm_{t}h'] = cc / TOT[chx] * 1e6
        tot_h = o[f'HAR_start_total_{t}h']
        o[f'HAR_cand_frame_frac_{t}h'] = hc / tot_h if tot_h else np.nan
        # HAR peak in cand frame within window
        vv = c['P'][har][ts] * cf
        if vv.max() > 0:
            idx = np.where(vv == vv.max())[0]
            cand_off = ts[idx] - t0
            pk = ts[idx[np.argmin(np.abs(cand_off))]]
            o[f'HAR_peak_pos_{t}h'] = int(c['m'][pk]); o[f'HAR_peak_off_{t}h'] = int(pk - t0); o[f'HAR_peak_n_{t}h'] = int(vv.max())
        else:
            o[f'HAR_peak_pos_{t}h'] = np.nan; o[f'HAR_peak_off_{t}h'] = np.nan; o[f'HAR_peak_n_{t}h'] = 0
        # spike descriptors in the HAR start window (all frames)
        sp = spike_stats(c['P29'][har], c['P30'][har], c['P'][har], ts, har)
        o[f'HAR_win_top1_frac_{t}h'] = sp['top1']; o[f'HAR_win_top1_pos_{t}h'] = int(c['m'][int(sp['top1_pos'])]) if sp['n'] else np.nan
        o[f'HAR_win_top3_frac_{t}h'] = sp['top3']; o[f'HAR_win_resid_deTop1_{t}h'] = sp['resid1']; o[f'HAR_win_resid_deTop3_{t}h'] = sp['resid3']
        o[f'HAR_win_f29_{t}h'] = sp['f29'] if har not in CAN_SKEW else np.nan
        o[f'HAR_win_f29_lib_{t}h'] = LIB29[har] if har not in CAN_SKEW else np.nan
        o[f'HAR_win_includes_pos21_22_{t}h'] = bool(any(g in (21, 22) for g in c['m'][ts]))
        o[f'HAR_peak_is_main_start_{t}h'] = bool((not np.isnan(o[f'HAR_peak_pos_{t}h'])) and o[f'HAR_peak_pos_{t}h'] == c['m'][mst])
        pos_ok = (not np.isnan(o[f'HAR_peak_off_{t}h'])) and THRESH['har_peak_lo'] <= o[f'HAR_peak_off_{t}h'] <= THRESH['har_peak_hi']
        o[f'HAR_start_signal_{t}h'] = bool(hc >= THRESH['har_min_reads'] and pos_ok and o[f'HAR_peak_n_{t}h'] >= THRESH['har_peak_min']
                                           and o[f'HAR_cand_frame_frac_{t}h'] >= THRESH['har_frame_frac'] and o[f'HAR_CHX_{t}h'] >= THRESH['har_chx_ratio']
                                           and not o[f'HAR_peak_is_main_start_{t}h'])
        if best is None or o[f'HAR_cand_ppm_{t}h'] > best[0]:
            best = (o[f'HAR_cand_ppm_{t}h'], t)
        # ---- CHX body ----
        end = (mst - 1) if r.ORF_class == 'upstream_inframe_extension' else (t1 - 3)
        end = min(end, L)
        cod = np.arange(t0, end - 1, 3) if end >= t0 + 2 else np.array([], int)
        cod = cod[cod + 2 <= end]
        nco = len(cod)
        vC = c['P'][chx]
        bts = np.arange(t0, end + 1)
        bc, bm = frames(bts, t0, mst)
        xb = vC[bts]
        cand_reads = int(xb[bc].sum()); main_reads = int(xb[bm & ~bc].sum()); oth_reads = int(xb[~bc & ~bm].sum())
        o[f'CHX_body_cand_{t}h'] = cand_reads; o[f'CHX_body_main_{t}h'] = main_reads; o[f'CHX_body_other_{t}h'] = oth_reads
        tot_b = cand_reads + main_reads + oth_reads
        o[f'CHX_body_cand_share_{t}h'] = cand_reads / tot_b if tot_b else np.nan
        cv = vC[cod] if nco else np.array([])
        ev = (c['R'][t][cod] >= THRESH['rna_eval_depth']) if nco else np.array([], bool)
        occ = (cv > 0)
        n_occ = int(occ.sum()); n_ev = int(ev.sum())
        o[f'occupied_codons_{t}h'] = n_occ; o[f'evaluable_codons_{t}h'] = n_ev; o[f'n_body_codons_{t}h'] = nco
        o[f'occupancy_fraction_{t}h'] = n_occ / n_ev if n_ev else np.nan
        run = best_run = 0
        for z in occ:
            run = run + 1 if z else 0
            best_run = max(best_run, run)
        o[f'longest_inframe_run_{t}h'] = best_run
        o[f'body_coverage_{t}h'] = float((xb > 0).mean()) if len(xb) else np.nan
        if n_occ:
            idxs = np.where(occ)[0]
            o[f'first_occupied_codon_{t}h'] = int(idxs[0] + 1); o[f'last_occupied_codon_{t}h'] = int(idxs[-1] + 1)
            o[f'occupied_span_nt_{t}h'] = int((idxs[-1] - idxs[0] + 1) * 3)
            so = np.sort(cv)[::-1]
            o[f'top1_fraction_{t}h'] = so[0] / cv.sum(); o[f'top3_fraction_{t}h'] = so[:3].sum() / cv.sum()
            o[f'de_top1_occupancy_{t}h'] = (n_occ - 1) / max(1, n_ev - 1); o[f'de_top3_occupancy_{t}h'] = max(0, n_occ - 3) / max(1, n_ev - 3)
            o[f'body_top1_pos_{t}h'] = int(c['m'][cod[int(np.argmax(cv))]])
        else:
            for kx in ['first_occupied_codon', 'last_occupied_codon', 'occupied_span_nt', 'top1_fraction', 'top3_fraction', 'de_top1_occupancy', 'de_top3_occupancy', 'body_top1_pos']:
                o[f'{kx}_{t}h'] = np.nan
        sp = spike_stats(c['P29'][chx], c['P30'][chx], vC, bts, chx)
        o[f'CHX_body_win_top1_frac_{t}h'] = sp['top1']; o[f'CHX_body_f29_{t}h'] = sp['f29'] if chx not in CAN_SKEW else np.nan
        o[f'CHX_body_f29_lib_{t}h'] = LIB29[chx] if chx not in CAN_SKEW else np.nan
        mx = cv.max() if nco else 0
        o[f'CHX_body_cand_deTop1_{t}h'] = int(cand_reads - mx)
        o[f'CHX_body_cand_share_deTop1_{t}h'] = (cand_reads - mx) / (tot_b - mx) if tot_b - mx > 0 else np.nan
        o[f'CHX_body_signal_deTop1_{t}h'] = bool((cand_reads - mx) >= THRESH['chx_min_reads'] and tot_b - mx > 0 and (cand_reads - mx) / (tot_b - mx) >= THRESH['chx_frame_share'])
        o[f'CHX_body_signal_{t}h'] = bool(cand_reads >= THRESH['chx_min_reads'] and o[f'CHX_body_cand_share_{t}h'] >= THRESH['chx_frame_share'])
        o[f'RNA_local_cov_{t}h'] = float(c['R'][t][ts].mean())
    o['strongest_HAR_tp'] = best[1]
    o['strongest_HAR_peak_position'] = o[f'HAR_peak_pos_{best[1]}h']; o['peak_offset_from_start'] = o[f'HAR_peak_off_{best[1]}h']
    o['n_HAR_positive_timepoints'] = sum(o[f'HAR_start_signal_{t}h'] for t in TPS)
    o['n_CHX_supported_timepoints_deTop1'] = sum(o[f'CHX_body_signal_deTop1_{t}h'] for t in TPS)
    o['n_CHX_supported_timepoints'] = sum(o[f'CHX_body_signal_{t}h'] for t in TPS)
    o['sgRNA_junction_support'] = int(jn.loc[r.transcript_id, 'reads_modal_pair']) if r.transcript_id != 'gRNA' else np.nan
    o['transcript_abundance_proxy_junction_ppm'] = jn.loc[r.transcript_id, 'reads_modal_pair'] / tot_j * 1e6 if r.transcript_id != 'gRNA' else np.nan
    rows.append(o)
ev = pd.DataFrame(rows)
res = cat.merge(ev, on='candidate_id')
res['THRESH'] = ''
res.to_csv(OUT + '_work/candidates_with_evidence.tsv', sep='\t', index=False)
print(len(res))
pd.set_option('display.width', 250, 'display.max_columns', 40)
for reg in ['5prime', '3prime']:
    d = res[res.region == reg]
    u = d.drop_duplicates('locus_key')
    print(reg, 'rows', len(d), 'loci', len(u))
    for nm, dd in (('rows', d), ('unique loci', u)):
        har = dd.n_HAR_positive_timepoints >= 1; chx = dd.n_CHX_supported_timepoints >= 1
        print(f'  {nm}: HAR>=1tp {har.sum()}, CHX>=1tp {chx.sum()}, both {(har & chx).sum()}, HAR>=2tp {(dd.n_HAR_positive_timepoints >= 2).sum()}, CHX>=2tp {(dd.n_CHX_supported_timepoints >= 2).sum()}')
for reg in ['5prime','3prime']:
    u = res[res.region==reg].drop_duplicates('locus_key'); print(reg,'CHX deTop1>=1tp loci',(u.n_CHX_supported_timepoints_deTop1>=1).sum(),' >=2tp',(u.n_CHX_supported_timepoints_deTop1>=2).sum())
print(res.groupby(['region', 'ORF_class']).agg(n=('candidate_id', 'count'), har=('n_HAR_positive_timepoints', lambda x: (x >= 1).sum()), chx=('n_CHX_supported_timepoints', lambda x: (x >= 1).sum())))
