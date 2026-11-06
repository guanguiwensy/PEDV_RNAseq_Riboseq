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
# keep ONLY query_length==29 alignments from the existing md BAM (same logic as pilot p29_filter.py), any sample. Output tree $O/price_29nt/all4
import pysam, os, json, sys
from collections import Counter
O = os.environ['O']; P = f'{O}/price_29nt/all4'; S = sys.argv[1]
md = f'{O}/price_813/bam_md/{S}.md.bam' if S == SR.lib("RIBO_VERO_CHX_1MOI_16h") else f'{O}/bam_md/{S}.md.bam'
src = pysam.AlignmentFile(md, 'rb', threads=4)
hv = pysam.AlignmentFile(f'{P}/virus/{S}.PEDV.29nt.bam', 'wb', header=src.header)
hh = pysam.AlignmentFile(f'{P}/host/pricebam_29nt/{S}.all.29nt.bam', 'wb', header=src.header)
c = Counter(); lenv = Counter(); lenh = Counter()
for r in src.fetch(until_eof=True):
    if r.is_unmapped: continue
    c['n_mapped_records'] += 1
    isv = r.reference_name == 'PEDV_GIb'
    ql = r.query_length
    (lenv if isv else lenh)[ql] += 1
    if ql != 29: continue
    if isv:
        hv.write(r); c['virus_29nt'] += 1; st = '-' if r.is_reverse else '+'
        c['virus_29nt_' + ('minus' if st == '-' else 'plus')] += 1
        if r.reference_length == 29: c['virus_29nt_refspan29'] += 1
    else:
        hh.write(r); c['host_29nt'] += 1
        if r.reference_length == 29: c['host_29nt_refspan29'] += 1
hv.close(); hh.close(); src.close()
pysam.index(f'{P}/virus/{S}.PEDV.29nt.bam'); pysam.index(f'{P}/host/pricebam_29nt/{S}.all.29nt.bam')
c['virus_all_lengths'] = sum(lenv.values()); c['host_all_lengths'] = sum(lenh.values())
c['virus_len16_50'] = sum(v for k, v in lenv.items() if k is not None and k <= 50); c['host_len16_50'] = sum(v for k, v in lenh.items() if k is not None and k <= 50)
d = dict(c); d['virus_len_hist'] = dict(sorted(lenv.items())); d['host_len_hist'] = dict(sorted(lenh.items()))
json.dump(d, open(f'{P}/logs/filter_counts_29nt_{S}.json', 'w'))
print(S, dict(c), flush=True)
