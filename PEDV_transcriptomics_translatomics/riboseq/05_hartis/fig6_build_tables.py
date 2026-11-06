# -*- coding: utf-8 -*-
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
import os, numpy as np, pandas as pd
H = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.abspath(os.path.join(H, '..')); TB = os.path.join(OUT, 'tables'); os.makedirs(TB, exist_ok=True)
B = os.path.abspath(os.path.join(H, '..', '..', '..'))
R = pd.read_csv(os.path.join(B, 'HARTIS_v3', 'G3_pedv_initiation', 'tables', 'HARTIS_v3_PEDV_start_level.tsv.gz'), sep='\t', dtype={'frame': str, 'read_length': str})
R = R[R['sid'].isin(SR.WHITELIST_LIBRARIES)].copy(); SR.assert_whitelisted(R['sid'].unique(), 'HARTIS start-level table')   # final sample table only; HARTIS fits are per library (independent)
C = pd.read_csv(os.path.join(B, 'HARTIS_v3', 'G3_pedv_initiation', 'tables', 'candidates_plus1_plus2.tsv'), sep='\t', dtype={'frame': str})
NM = pd.read_csv(os.path.join(B, 'global_perturbation_screen', 'plus1_orf_capacity', 'panelC_family_name_mapping.tsv'), sep='\t')
PM = pd.read_csv(os.path.join(B, 'orf_calling_GIb', '_price_final', 'annotated_29nt_resultlevel_merge_of_3_samples_named.tsv'), sep='\t')
PL = {'CHX5': SR.lib("RIBO_VERO_CHX_1MOI_5h"), 'CHX8': SR.lib("RIBO_VERO_CHX_1MOI_8h"), 'CHX16': SR.lib("RIBO_VERO_CHX_1MOI_16h")}
CANON = [('ORF1a', 293, 12622), ('ORF1b', 12592, 20613), ('S', 20610, 24758), ('ORF3', 24758, 25033), ('E', 25361, 25600), ('M', 25608, 26288), ('N', 26300, 27625)]
checks = []


fam_name = dict(zip(NM.original_family_id, NM.displayed_name))
checks.append(('69 个 +1 家族显示名唯一', NM.displayed_name.is_unique))
checks.append(('家族 id ↔ 共同终止一一对应', NM.shared_stop.is_unique and NM.original_family_id.is_unique))
c1 = C[C.frame == '+1']
checks.append(('HARTIS +1 家族集合 = 名称映射集合', set(c1.family_id) == set(NM.original_family_id)))
checks.append(('每个家族只有一个终止位点', int(C.groupby('family_id').stop_position.nunique().max()) == 1))
checks.append(('每个家族只有一个父 ORF 与一个框', int(C.groupby('family_id')[['parent', 'frame']].nunique().max().max()) == 1))
C['family_name'] = C.family_id.map(fam_name).fillna(C.family_id)
R['family_name'] = R.family_id.map(fam_name).fillna(R.family_id)


def price_lib(sid):
    d = pd.read_csv(os.path.join(B, 'orf_calling_GIb', '_price_final', 'annotated_29nt_named', f'{sid}.plus.orfs.annotated.tsv'), sep='\t')
    p = d['p value'].values.astype(float); n = len(p); o = np.argsort(p); q = np.empty(n); prev = 1.0
    for i in range(n - 1, -1, -1):
        prev = min(prev, p[o[i]] * n / (i + 1)); q[o[i]] = prev
    d['q_BH_within_file'] = q
    return d[['our_start_1based', 'our_stop_end_1based', 'p value', 'q_BH_within_file', 'Codon']].rename(columns={'p value': 'p'})
PLIB = {k: price_lib(v) for k, v in PL.items()}; PMERGED = price_lib('merged_CHX')
for k, sid in PL.items():
    rep = {SR.lib("RIBO_VERO_CHX_1MOI_5h"): 4, SR.lib("RIBO_VERO_CHX_1MOI_8h"): 8, SR.lib("RIBO_VERO_CHX_1MOI_16h"): 7}[sid]
    checks.append((f'PRICE {k} 文件内 BH q<0.1 个数与 PRICE_29nt_REPORT 一致（{rep}）', int((PLIB[k].q_BH_within_file < 0.1).sum()) == rep))
checks.append(('PRICE merged_CHX 文件内 BH q<0.1 个数 = 6', int((PMERGED.q_BH_within_file < 0.1).sum()) == 6))
fam_stop_end = {fid: s + 2 for fid, s in zip(NM.original_family_id, NM.shared_stop)}
fam_stop_end.update({f: s + 2 for f, s in C[C.frame == '+2'].groupby('family_id').stop_position.first().items()})
PRICE_rows = []
for fid, se in fam_stop_end.items():
    r = dict(family_id=fid, PRICE_stop_end=se)
    m = PM[PM.stop == se]
    r['PRICE_called'] = len(m) > 0
    if len(m) > 1: r['PRICE_note'] = 'multiple_PRICE_rows_same_stop'
    if len(m):
        m = m.iloc[0]; r['PRICE_ORF_name'] = m.ORF_name; r['PRICE_coordinate_class'] = m.coordinate_class; r['PRICE_frame_relative_to_parent'] = m.frame_relative_to_parent
        ps = {}
        for k in PL:
            st, p = m[f'start_{k}'], m[f'p_{k}']
            if pd.notna(st):
                x = PLIB[k][(PLIB[k].our_start_1based == int(st)) & (PLIB[k].our_stop_end_1based == se)]
                q = float(x.q_BH_within_file.iloc[0]) if len(x) else np.nan
                ps[k] = (int(st), float(p), q)
                r[f'PRICE_start_{k}'] = int(st); r[f'PRICE_p_{k}'] = float(p); r[f'PRICE_q_{k}'] = q
        if ps:
            kbest = min(ps, key=lambda k: ps[k][1])
            r['PRICE_best_lib'] = kbest; r['PRICE_best_start'] = ps[kbest][0]; r['PRICE_best_p'] = ps[kbest][1]; r['PRICE_best_q_BH_within_file'] = ps[kbest][2]
            r['PRICE_best_p_reported_zero'] = ps[kbest][1] == 0
            r['PRICE_n_libs_called'] = len(ps); r['PRICE_n_libs_p05'] = sum(v[1] < 0.05 for v in ps.values())
            r['PRICE_starts_differ_between_libs'] = len({v[0] for v in ps.values()}) > 1
        mg = PMERGED[PMERGED.our_stop_end_1based == se]
        if len(mg):
            mg = mg.sort_values('p').iloc[0]; r['PRICE_mergedCHX_start'] = int(mg.our_start_1based); r['PRICE_mergedCHX_p'] = float(mg.p); r['PRICE_mergedCHX_q'] = float(mg.q_BH_within_file)
    PRICE_rows.append(r)
PR = pd.DataFrame(PRICE_rows)


key = ['frame', 'family_id', 'family_name', 'parent', 'start_position', 'start_codon', 'stop_position', 'ORF_length_aa', 'candidate_frame_equals_canonical']
S = C.copy()
def summarize(d, drug, window, prefix):
    x = d[(d.drug == drug) & (d.window == window)]
    out = {}
    ins = x[x.status == 'in_support']; cap = x[x.status == 'capped_reference']
    libs_ins = ins.groupby('label').p_emp.min()
    out[f'{prefix}_n_libs_valid_in_support'] = len(libs_ins)
    out[f'{prefix}_n_libs_p05_in_support'] = int((libs_ins <= 0.05).sum())
    if len(ins):
        b = ins.loc[ins.p_emp.idxmin()]
        out[f'{prefix}_min_nominal_p'] = b.p_emp; out[f'{prefix}_Z_at_min_p'] = b.Z_cand; out[f'{prefix}_min_p_source'] = f'{b.label}|{b.method}'
        out[f'{prefix}_K_at_min_p'] = b.K
    libs_cap = cap.groupby('label').p_cap.min()
    out[f'{prefix}_n_libs_capped'] = len(libs_cap); out[f'{prefix}_n_libs_pcap05'] = int((libs_cap <= 0.05).sum())
    if len(cap):
        b = cap.loc[cap.p_cap.idxmin()]
        out[f'{prefix}_capped_min_p'] = b.p_cap; out[f'{prefix}_capped_Zcap_at_min'] = b.Z_cap; out[f'{prefix}_capped_source'] = f'{b.label}|{b.method}'
    out[f'{prefix}_status_counts'] = ';'.join(f'{a}:{b}' for a, b in x.status.value_counts().items())
    return out
rows = []
for k, d in R.groupby(['frame', 'family_id', 'start_position']):
    r = {}
    for drug in ('HAR', 'CHX'):
        r.update(summarize(d, drug, 'W1', f'{drug}_W1')); r.update(summarize(d, drug, 'W5', f'{drug}_W5'))
    r['frame'], r['family_id'], r['start_position'] = k
    rows.append(r)
ST = S.merge(pd.DataFrame(rows), on=['frame', 'family_id', 'start_position'], how='left')
ST = ST.merge(PR[['family_id', 'PRICE_called', 'PRICE_best_p', 'PRICE_best_start', 'PRICE_best_lib']], on='family_id', how='left')
ST['PRICE_best_start_is_this_start'] = ST.PRICE_best_start == ST.start_position
ST['PRICE_note'] = 'PRICE 为 ORF（终止位点）级证据，此处仅作家族映射参考，不是本起点的 P'
ST['p_value_type'] = 'raw nominal (HARTIS: minimum nominal P across libraries×lengths, uncalibrated, not a joint P; PRICE: family-level, see family table)'
ST.to_csv(os.path.join(OUT, 'HARTIS_PRICE_start_summary.tsv'), sep='\t', index=False, float_format='%.5g')


F = []
for fid, d in ST.groupby('family_id'):
    r = dict(family_id=fid, family_name=d.family_name.iloc[0], frame=d.frame.iloc[0], parent_region=d.parent.iloc[0], stop_position=int(d.stop_position.iloc[0]),
             n_start_variants=len(d), ORF_length_aa_range=f'{d.ORF_length_aa.min()}-{d.ORF_length_aa.max()}',
             candidate_frame_equals_canonical=';'.join(sorted(set(d.candidate_frame_equals_canonical.dropna().astype(str)))) or '')
    for drug in ('HAR', 'CHX'):
        pre = f'{drug}_W1'
        e = d[d[f'{pre}_min_nominal_p'].notna()]
        x = R[(R.family_id == fid) & (R.drug == drug) & (R.window == 'W1')]
        libs_valid = x[x.status == 'in_support'].label.nunique(); libs_sup = x[(x.status == 'in_support') & (x.p_emp <= 0.05)].label.nunique()
        r[f'{drug}_W1_n_starts_with_valid_p'] = len(e)
        r[f'{drug}_W1_n_libs_valid'] = libs_valid; r[f'{drug}_W1_n_libs_support_p05'] = libs_sup
        if len(e):
            b = e.loc[e[f'{pre}_min_nominal_p'].idxmin()]
            r[f'{drug}_W1_best_start'] = int(b.start_position); r[f'{drug}_W1_best_start_codon'] = b.start_codon
            r[f'{drug}_W1_best_p'] = b[f'{pre}_min_nominal_p']; r[f'{drug}_W1_Z_at_best'] = b[f'{pre}_Z_at_min_p']; r[f'{drug}_W1_best_source'] = b[f'{pre}_min_p_source']
        ec = d[d[f'{pre}_capped_min_p'].notna()]
        r[f'{drug}_W1_n_starts_capped'] = len(ec)
        if len(ec):
            b = ec.loc[ec[f'{pre}_capped_min_p'].idxmin()]
            r[f'{drug}_W1_capped_best_start'] = int(b.start_position); r[f'{drug}_W1_capped_best_p'] = b[f'{pre}_capped_min_p']; r[f'{drug}_W1_capped_Zcap'] = b[f'{pre}_capped_Zcap_at_min']
            r[f'{drug}_W1_capped_n_libs_pcap05'] = int(R[(R.family_id == fid) & (R.drug == drug) & (R.window == 'W1') & (R.status == 'capped_reference') & (R.p_cap <= 0.05)].label.nunique())
        e5 = d[d[f'{drug}_W5_min_nominal_p'].notna()]
        if len(e5):
            b = e5.loc[e5[f'{drug}_W5_min_nominal_p'].idxmin()]; r[f'{drug}_W5_best_p'] = b[f'{drug}_W5_min_nominal_p']; r[f'{drug}_W5_best_start'] = int(b.start_position)
        e5c = d[d[f'{drug}_W5_capped_min_p'].notna()]
        if len(e5c):
            b = e5c.loc[e5c[f'{drug}_W5_capped_min_p'].idxmin()]; r[f'{drug}_W5_capped_best_p'] = b[f'{drug}_W5_capped_min_p']; r[f'{drug}_W5_capped_best_start'] = int(b.start_position)
    F.append(r)
F = pd.DataFrame(F).merge(PR, on='family_id', how='left')
F['n_valid_libraries_total'] = F.HAR_W1_n_libs_valid + F.CHX_W1_n_libs_valid
F['PRICE_frame_mismatch'] = F.PRICE_called.fillna(False) & (F.PRICE_frame_relative_to_parent.fillna(-1) != F.frame.str[1].astype(int))
def har_status(r):
    if pd.notna(r.get('HAR_W1_best_p')): return 'support_p05' if r.HAR_W1_best_p <= 0.05 else 'not_significant'
    if pd.notna(r.get('HAR_W1_capped_best_p')): return 'capped_only(exploratory)'
    return 'not_assessable'
def price_status(r):
    if not bool(r.PRICE_called) or pd.isna(r.get('PRICE_best_p')): return 'not_called(not_assessable)'
    if r.PRICE_frame_mismatch: return 'frame_mismatch(canonical_ORF)'
    return 'support_p05' if r.PRICE_best_p < 0.05 else 'not_significant'
F['HAR_W1_status'] = F.apply(har_status, axis=1); F['PRICE_status'] = F.apply(price_status, axis=1)
def joint(r):
    h, p = r.HAR_W1_status, r.PRICE_status
    if h == 'support_p05' and p == 'support_p05': return 'PRICE_and_HAR_nominal'
    if h == 'support_p05' and p == 'not_significant': return 'HAR_only_nominal'
    if h == 'not_significant' and p == 'support_p05': return 'PRICE_only_nominal'
    if h == 'not_significant' and p == 'not_significant': return 'neither_nominal(both_assessed)'
    if h == 'support_p05': return 'HAR_nominal_PRICE_not_assessable'
    if p == 'support_p05': return 'PRICE_nominal_HAR_not_assessable'
    return 'not_assessable(at_least_one_line)'
F['joint_category_exploratory'] = F.apply(joint, axis=1)
F['exploratory_candidate'] = ((F.HAR_W1_status == 'support_p05') | (F.PRICE_status == 'support_p05') | (F.CHX_W1_best_p <= 0.05)
                              | (F.HAR_W1_capped_n_libs_pcap05.fillna(0) >= 2))
F['p_value_type'] = 'raw nominal; HARTIS = minimum nominal P over starts×libraries×lengths (best-start selection bias; not calibrated, not joint); PRICE best = min raw P over 3 CHX libraries; q = BH within each PRICE file'
F['interpretation_note'] = np.where(F.candidate_frame_equals_canonical != '', '候选框与经典 ORF 读框相同：' + F.candidate_frame_equals_canonical, '')
F = F.sort_values(['frame', 'parent_region', 'stop_position'])
F.to_csv(os.path.join(OUT, 'HARTIS_PRICE_family_summary.tsv'), sep='\t', index=False, float_format='%.5g')
checks.append(('PRICE 与家族读框不一致的家族', ';'.join(F[F.PRICE_frame_mismatch].family_name)))
checks.append(('PRICE 报告 p = 0 的家族', ';'.join(F[F.PRICE_best_p_reported_zero.fillna(False)].family_name)))
pd.DataFrame(checks, columns=['check', 'result']).to_csv(os.path.join(TB, 'data_checks.tsv'), sep='\t', index=False)
CT = pd.crosstab(F[F.frame == '+1'].HAR_W1_status, F[F.frame == '+1'].PRICE_status, margins=True)
CT.to_csv(os.path.join(TB, 'family_crosstab_plus1_HAR_W1_vs_PRICE.tsv'), sep='\t')
pd.set_option('display.width', 250)
for c in checks: print(c)
print(CT.to_string()); print(F[F.frame == '+1'].joint_category_exploratory.value_counts().to_string())
