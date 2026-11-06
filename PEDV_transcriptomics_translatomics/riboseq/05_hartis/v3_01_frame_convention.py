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
import os, json, gzip, re, collections, numpy as np, pandas as pd
H = os.path.dirname(os.path.abspath(__file__)); V3 = os.path.join(H, '..'); R = os.path.join(V3, '..')
RG = os.path.join(R, 'riboseq_GIb_mapping_psite'); ST = os.path.join(RG, '_work', 'server', 'stats')
T = os.path.join(V3, 'tables'); os.makedirs(T, exist_ok=True)
OFF = {s: {int(k): v for k, v in d.items()} for s, d in json.load(open(os.path.join(RG, '_work', 'offset_table_GIb.json'))).items()}
SAMPLES = [SR.lib("RIBO_VERO_CHX_1MOI_5h"), SR.lib("RIBO_VERO_HAR_1MOI_5h"), SR.lib("RIBO_VERO_CHX_1MOI_8h"), SR.lib("RIBO_VERO_HAR_1MOI_8h"), SR.lib("RIBO_VERO_CHX_1MOI_16h"), SR.lib("RIBO_VERO_HAR_1MOI_16h")]
ORFS = [('ORF1a', 293, 12622), ('ORF1b', 12592, 20613), ('S', 20610, 24758), ('ORF3', 24758, 25033),
        ('E', 25361, 25600), ('M', 25608, 26288), ('N', 26300, 27625)]
GL = 27959; TRIM = 45


def rd(p):
    return pd.read_csv(p, sep='\t', comment='#')



cover = np.zeros(GL + 2, int)
for _, a, b in ORFS: cover[a:b + 1] += 1
body = {}
for n, a, b in ORFS:
    p = np.arange(a + TRIM, b - TRIM + 1); body[n] = (a, set(p[cover[p] == 1].tolist()))
PV = rd(os.path.join(RG, 'periodicity_virus.tsv'))
rowsA = []
for s in SAMPLES:
    m = pd.read_csv(os.path.join(ST, s, 'pedv_5p_matrix.tsv.gz'), sep='\t'); m = m[m.strand == '+']
    m = m.groupby(['L', 'pos5'], as_index=False).n.sum()
    for n, (a, bset) in body.items():
        mb = m[m.pos5.isin(bset)]
        for L, d in mb.groupby('L'):
            o = OFF[s].get(int(min(L, 150)))
            if o is None: continue
            fr = (d.pos5.values + o - a) % 3
            c = np.array([d.n.values[fr == k].sum() for k in range(3)]); t = c.sum()
            ref = PV[(PV.sample_id == s) & (PV.domain == f'PEDV_{n}') & (PV.read_length.astype(str) == str(L))]
            if len(ref) == 0: continue
            r = ref.iloc[0]
            rowsA.append(dict(sample_id=s, ORF=n, L=int(L), n_recomputed=int(t), n_table=int(r['n_reads(5′端在CDS主体)']),
                              F0=round(c[0] / t, 3), F1=round(c[1] / t, 3), F2=round(c[2] / t, 3),
                              F0_tab=r.psite_frame0, F1_tab=r.psite_frame1, F2_tab=r.psite_frame2))
A = pd.DataFrame(rowsA)
A['match'] = (A.n_recomputed == A.n_table) & (np.abs(A.F0 - A.F0_tab) < 6e-4) & (np.abs(A.F1 - A.F1_tab) < 6e-4) & (np.abs(A.F2 - A.F2_tab) < 6e-4)
A.to_csv(os.path.join(T, 'frame_check_virus_vs_periodicity_virus.tsv'), sep='\t', index=False)
print('A virus rows', len(A), 'all match', A.match.all(), 'mismatch', (~A.match).sum())


PH = rd(os.path.join(RG, 'periodicity_host.tsv'))
rowsB = []
for s in SAMPLES:
    f5 = rd(os.path.join(ST, s, 'frame5.tsv')); f5 = f5[f5.domain == 'host']
    g = pd.read_csv(os.path.join(RG, '_work', 'server', 'host_psite', s, 'host_gene_cds_psite.tsv'), sep='\t')
    for L in range(25, 35):
        o = OFF[s].get(L)
        d = f5[f5.L == L]
        if o is None or d.empty: continue
        c5 = d[['f0', 'f1', 'f2']].values.sum(0)
        c = np.array([c5[(k - o) % 3] for k in range(3)])
        t = c.sum()
        ref = PH[(PH.sample_id == s) & (PH.read_length.astype(str) == str(L))]
        gg = g[g.L == L][['f0', 'f1', 'f2']].values.sum(0); tg = gg.sum()
        r = ref.iloc[0] if len(ref) else None
        rowsB.append(dict(sample_id=s, L=L, offset=o, n_B1=int(t), n_table=int(r['n_reads(5′端在CDS主体)']) if r is not None else -1,
                          F0_B1=round(c[0] / t, 3), F1_B1=round(c[1] / t, 3), F2_B1=round(c[2] / t, 3),
                          F0_tab=r.psite_frame0 if r is not None else np.nan, F1_tab=r.psite_frame1 if r is not None else np.nan,
                          F2_tab=r.psite_frame2 if r is not None else np.nan,
                          n_B2_wholeCDS=int(tg), F0_B2=round(gg[0] / tg, 3), F1_B2=round(gg[1] / tg, 3), F2_B2=round(gg[2] / tg, 3),
                          F1_over_F1F2_tab=round(r.psite_frame1 / (r.psite_frame1 + r.psite_frame2), 3) if r is not None else np.nan,
                          F1_over_F1F2_B2=round(gg[1] / (gg[1] + gg[2]), 3)))
B = pd.DataFrame(rowsB)
B['B1_match'] = (B.n_B1 == B.n_table) & (np.abs(B.F0_B1 - B.F0_tab) < 6e-4) & (np.abs(B.F1_B1 - B.F1_tab) < 6e-4)
B.to_csv(os.path.join(T, 'frame_check_host_vs_periodicity_host.tsv'), sep='\t', index=False)
print('B host rows', len(B), 'B1 all match', B.B1_match.all())
pd.set_option('display.width', 250)
print(B[B.L.isin([28, 29, 30, 31])][['sample_id', 'L', 'F0_tab', 'F0_B1', 'F0_B2', 'F1_over_F1F2_tab', 'F1_over_F1F2_B2', 'B1_match']].to_string())


gtf = os.path.join(R, '..', 'refdata_vero', 'Chlorocebus_sabaeus.ChlSab1.1.116.gtf.gz')
strand = {}
for line in gzip.open(gtf, 'rt'):
    if line.startswith('#'): continue
    f = line.split('\t')
    if f[2] == 'transcript' and 'Ensembl_canonical' in f[8]:
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1); strand[t] = f[6]
rowsS = []
for s in SAMPLES:
    g = pd.read_csv(os.path.join(RG, '_work', 'server', 'host_psite', s, 'host_gene_cds_psite.tsv'), sep='\t')
    g = g[g.L == 29]; g['strand'] = g.transcript_id.map(strand)
    for st, d in g.groupby('strand'):
        c = d[['f0', 'f1', 'f2']].values.sum(0)
        rowsS.append(dict(sample_id=s, L=29, strand=st, n_genes=d.transcript_id.nunique(), n=int(c.sum()),
                          F0=round(c[0] / c.sum(), 3), F1=round(c[1] / c.sum(), 3), F2=round(c[2] / c.sum(), 3),
                          F1_over_F1F2=round(c[1] / (c[1] + c[2]), 3)))
S = pd.DataFrame(rowsS); S.to_csv(os.path.join(T, 'frame_check_host_by_strand_29nt.tsv'), sep='\t', index=False); print(S.to_string())


PS = pd.read_csv(os.path.join(RG, 'psite_counts_PEDV_GIb.tsv.gz'), sep='\t', comment='#',
                 usecols=['sample_id', 'strand', 'position', 'read_length', 'psite_count'])
PS = PS[(PS.strand == '+') & (PS.read_length == 29)]
rowsC = []
for s in SAMPLES:
    d = PS[PS.sample_id == s]
    for lab, lo, hi in [('ORF1a_body_293+45..12585', 338, 12585), ('ORF1b_body_12623..20568', 12623, 20568)]:
        x = d[(d.position >= lo) & (d.position <= hi)]
        fr = (x.position.values - 293) % 3; c = np.array([x.psite_count.values[fr == k].sum() for k in range(3)])
        rowsC.append(dict(sample_id=s, region=lab, frame_reference='ORF1a (start 293)', n=int(c.sum()),
                          F0=round(c[0] / c.sum(), 3), F1=round(c[1] / c.sum(), 3), F2=round(c[2] / c.sum(), 3)))
C = pd.DataFrame(rowsC); C.to_csv(os.path.join(T, 'frame_check_ORF1b_in_ORF1a_frame.tsv'), sep='\t', index=False); print(C.to_string())
print('(12592-293)%3 =', (12592 - 293) % 3)
