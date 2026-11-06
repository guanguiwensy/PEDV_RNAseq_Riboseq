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
import os, gzip, re, collections, pandas as pd
H = os.path.dirname(os.path.abspath(__file__)); V3 = os.path.join(H, '..'); ROOT = os.path.join(V3, '..', '..')
GTF = os.path.join(ROOT, 'refdata_vero', 'Chlorocebus_sabaeus.ChlSab1.1.116.gtf.gz'); FA = os.path.join(ROOT, 'refdata_vero', 'Chlorocebus_sabaeus.ChlSab1.1.dna.toplevel.fa.gz')
TX = pd.read_csv(os.path.join(V3, '_server', SR.lib("RIBO_VERO_CHX_1MOI_8h"), 'transcripts.tsv'), sep='\t')
rx = re.compile(r'(\S+) "([^"]*)"')
ex = collections.defaultdict(list); info = {}
for line in gzip.open(GTF, 'rt'):
    if line.startswith('#'): continue
    f = line.split('\t')
    if f[2] != 'exon' or 'Ensembl_canonical' not in f[8]: continue
    d = dict(rx.findall(f[8]))
    if d.get('transcript_biotype') != 'protein_coding': continue
    ex[d['transcript_id']].append((int(f[3]), int(f[4]))); info[d['transcript_id']] = (f[0], f[6])
need = collections.defaultdict(list)
for t in TX.transcript_id:
    c, s = info[t]; need[c].append(t)
comp = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
seqs = {}


def flush(chrom, s):
    for t in need.get(chrom, []):
        st = info[t][1]; el = sorted(ex[t], key=lambda x: x[0])
        parts = [s[b - 1:e] for b, e in el]; q = ''.join(parts).upper()
        if st == '-': q = q.translate(comp)[::-1]
        seqs[t] = q


chrom = None; buf = []
for line in gzip.open(FA, 'rt'):
    if line[0] == '>':
        if chrom in need: flush(chrom, ''.join(buf))
        chrom = line[1:].split()[0]; buf = []
    elif chrom in need: buf.append(line.rstrip())
if chrom in need: flush(chrom, ''.join(buf))
TX['seq'] = TX.transcript_id.map(seqs)
TX['len_ok'] = TX.seq.str.len() == TX.tx_len
TX['start_codon'] = [s[c:c + 3] if isinstance(s, str) else None for s, c in zip(TX.seq, TX.cds_start_tx)]
print('transcripts', len(TX), 'seq found', TX.seq.notna().sum(), 'length match', TX.len_ok.sum(), 'ATG at CDS start', (TX.start_codon == 'ATG').sum())
TX[['tx_idx', 'transcript_id', 'gene_name', 'strand', 'tx_len', 'cds_start_tx', 'cds_end_tx', 'cds_len_nt', 'split_cds_gap_le10', 'start_codon', 'seq']].to_csv(
    os.path.join(V3, 'tables', 'g2_tx_seq.tsv.gz'), sep='\t', index=False)
