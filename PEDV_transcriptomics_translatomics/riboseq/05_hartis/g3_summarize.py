# -*- coding: utf-8 -*-
"""HARTIS_v3 PEDV 起始占据结果汇总（只用 g3_hartis_pedv_initiation 的输出；不新增检验）。
- 正式结果只来自 status = in_support 的 EMP_Zs（单读长）与 J_POOL（W5）；capped_reference 单列，探索性。
- 文库层面的重复：同一文库的多个读长不是独立证据，按"文库"计数——某库任一合格读长（或 J_POOL）p_emp ≤ 0.05 即该库命中。
- 家族：同父 ORF、同框、同终止；家族统计量为家族内最小 p（只作描述；宿主伪家族零分布未建立），q 为 BH（描述性，multiplicity_not_certified）。
"""
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
import os, numpy as np, pandas as pd
H = os.path.dirname(os.path.abspath(__file__)); TB = os.path.join(H, '..', 'tables'); R0 = os.path.abspath(os.path.join(H, '..', '..', '..'))
R = pd.read_csv(os.path.join(TB, 'HARTIS_v3_PEDV_start_level.tsv.gz'), sep='\t', dtype={'frame': str, 'read_length': str})
R = R[R['sid'].isin(SR.WHITELIST_LIBRARIES)].copy(); SR.assert_whitelisted(R['sid'].unique(), 'HARTIS start-level table')   # final sample table only; HARTIS fits are per library (independent)
NM = pd.read_csv(os.path.join(R0, 'global_perturbation_screen', 'plus1_orf_capacity', 'panelC_family_name_mapping.tsv'), sep='\t')
nm = dict(zip(NM.original_family_id, NM.displayed_name))
R['family_name'] = R.family_id.map(nm).fillna(R.family_id)
R['canonical_note'] = np.select([R.candidate_frame_equals_canonical.fillna('') == 'S', R.candidate_frame_equals_canonical.fillna('') == 'ORF3',
                                 R.candidate_frame_equals_canonical.fillna('') == 'M'],
                                ['= 经典 S 读框（S 起始/上游）', '= 经典 ORF3 读框（ORF3 起始/上游延伸）', '= 经典 M 读框（M 上游 N 端延伸）'], '')
R['lib'] = R.label
ok = R.status == 'in_support'


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n); prev = 1.0
    for i in range(n - 1, -1, -1):
        prev = min(prev, p[o[i]] * n / (i + 1)); q[o[i]] = prev
    return q


# 1. 起点 × 框 × 窗口：按文库的命中
key = ['frame', 'family_id', 'family_name', 'parent', 'start_position', 'start_codon', 'canonical_note', 'window']
S = []
for k, d in R.groupby(key, dropna=False):
    r = dict(zip(key, k))
    for drug in ('CHX', 'HAR'):
        x = d[(d.drug == drug) & ok.loc[d.index]]
        libs = x.groupby('lib').agg(minp=('p_emp', 'min'), maxZ=('Z_cand', 'max'))
        libs['minp_cp'] = x[x.transfer_flag == 'checks_passed'].groupby('lib').p_emp.min().reindex(libs.index)   # 只在迁移检查通过的行内取最小 p
        r[f'{drug}_libs_assessed'] = len(libs); r[f'{drug}_libs_p05'] = int((libs.minp <= .05).sum())
        r[f'{drug}_libs_p05_checks_passed'] = int((libs.minp_cp <= .05).sum())
        r[f'{drug}_min_p'] = libs.minp.min() if len(libs) else np.nan; r[f'{drug}_max_Z'] = libs.maxZ.max() if len(libs) else np.nan
        r[f'{drug}_hit_libs'] = ';'.join(libs.index[libs.minp <= .05])
        c = d[(d.drug == drug) & (d.status == 'capped_reference')]
        r[f'{drug}_capped_n'] = c.lib.nunique(); r[f'{drug}_capped_maxZcap'] = c.Z_cap.max() if len(c) else np.nan
        r[f'{drug}_capped_libs_pcap05'] = c[c.p_cap <= .05].lib.nunique()
    r['status_summary'] = ';'.join(f'{a}:{b}' for a, b in d.status.value_counts().items())
    S.append(r)
S = pd.DataFrame(S); S.to_csv(os.path.join(TB, 'start_level_library_replication.tsv'), sep='\t', index=False, float_format='%.4g')

# 2. 家族 × 方法 × 窗口：家族内最小 p、BH（描述）
F = []
for (fr, w, meth, lib), d in R[ok].groupby(['frame', 'window', 'method', 'lib']):
    g = d.groupby(['family_id', 'family_name']).agg(n_starts=('p_emp', 'size'), min_p=('p_emp', 'min'), max_Z=('Z_cand', 'max'),
                                                    drug=('drug', 'first'), transfer=('transfer_flag', 'first')).reset_index()
    g['min_p_bonf_within_family'] = (g.min_p * g.n_starts).clip(upper=1)
    g['q_BH_descriptive'] = bh(g.min_p_bonf_within_family.values)
    g['frame'] = fr; g['window'] = w; g['method'] = meth; g['lib'] = lib; g['n_families_in_set'] = len(g)
    F.append(g)
F = pd.concat(F); F['flag'] = 'multiplicity_not_certified;tail_uncertified_below_0.01'
F.to_csv(os.path.join(TB, 'family_level_in_support.tsv'), sep='\t', index=False, float_format='%.4g')

# 3. 状态总览（每个起点在每个方法下的状态）
ST = R.groupby(['frame', 'window', 'parent', 'status']).size().unstack(fill_value=0)
ST.to_csv(os.path.join(TB, 'status_overview_by_parent.tsv'), sep='\t')
pd.set_option('display.width', 260); pd.set_option('display.max_rows', 200)
print(ST.to_string())
rep = S[(S.CHX_libs_p05 + S.HAR_libs_p05) >= 1]
print('\n至少 1 个文库 p ≤ 0.05 的起点 × 窗口：', len(rep))
cols = ['frame', 'window', 'family_name', 'start_position', 'start_codon', 'canonical_note', 'CHX_libs_assessed', 'CHX_libs_p05', 'CHX_hit_libs', 'HAR_libs_assessed', 'HAR_libs_p05', 'HAR_hit_libs', 'CHX_min_p', 'HAR_min_p']
print(S[(S.CHX_libs_p05 >= 2) | (S.HAR_libs_p05 >= 2) | ((S.CHX_libs_p05 >= 1) & (S.HAR_libs_p05 >= 1))][cols].sort_values(['frame', 'window', 'start_position']).to_string(index=False))
print(F[F.q_BH_descriptive <= 0.1][['frame', 'window', 'method', 'lib', 'family_name', 'n_starts', 'min_p', 'q_BH_descriptive', 'transfer']].to_string(index=False))
