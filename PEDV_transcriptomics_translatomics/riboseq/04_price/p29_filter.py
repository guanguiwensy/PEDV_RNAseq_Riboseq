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
# PRICE 29-nt pilot: keep ONLY alignments with query_length == 29 from the existing md BAM of RIBO_VERO_CHX_1MOI_16h (CHX16h).
# Same logic as p813_filter.py except the length rule (==29 instead of <=50). New output tree: $O/price_29nt
import pysam, os, json
from collections import Counter
O = os.environ['O']; P = f'{O}/price_29nt'; S = SR.lib("RIBO_VERO_CHX_1MOI_16h")
md = f'{O}/price_813/bam_md/{S}.md.bam'
src = pysam.AlignmentFile(md, 'rb', threads=4)
hv = pysam.AlignmentFile(f'{P}/virus/{S}.PEDV.29nt.bam', 'wb', header=src.header)
hh = pysam.AlignmentFile(f'{P}/host/pricebam_29nt/{S}.all.29nt.bam', 'wb', header=src.header)
n_in = 0; nv = nh = 0
lenv = Counter(); lenh = Counter()
for r in src.fetch(until_eof=True):
    if r.is_unmapped: continue
    n_in += 1
    isv = r.reference_name == 'PEDV_GIb'
    ql = r.query_length
    (lenv if isv else lenh)[ql] += 1   # full length histogram (all mapped records) for reference
    if ql != 29: continue
    if isv: hv.write(r); nv += 1
    else: hh.write(r); nh += 1
hv.close(); hh.close(); src.close()
pysam.index(f'{P}/virus/{S}.PEDV.29nt.bam'); pysam.index(f'{P}/host/pricebam_29nt/{S}.all.29nt.bam')
json.dump({'n_mapped_records': n_in, 'virus_29nt': nv, 'host_29nt': nh,
           'virus_len_hist': dict(sorted(lenv.items())), 'host_len_hist': dict(sorted(lenh.items()))}, open(f'{P}/logs/filter_counts_29nt_{S}.json', 'w'))
print(S, n_in, 'virus29', nv, 'host29', nh, flush=True)
