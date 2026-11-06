"""Reclassify TRS-L-dependent noncanonical junctions (1,173) by the ORFs their sgRNA can encode.
Rules (agreed with the user):
  translates_canonical (excluded): sgRNA contains an ATG-initiated, in-frame-stop ORF that starts inside a canonical
      CDS in the canonical reading frame. These stay noncanonical sgRNAs but are NOT counted in plots/statistics.
  truncated (counted): first downstream AUG in frame with the acceptor's CDS (Fig. 2 frame == 'in').
  frameshift (counted): first downstream AUG out of frame with the acceptor's CDS (Fig. 2 frame == 'out').
  out-of-frame ORF (counted): acceptor with no AUG in its CDS (or noncoding acceptor) whose sgRNA contains an
      ATG...stop ORF, none of which is in canonical frame.
  no ORF: no ATG...stop ORF in the sgRNA (none observed here).
ORF search: any ATG with an in-frame stop inside the molecule (leader donor + body from acceptor to genome end). No minimum
length was imposed (open question for the user).
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
from pathlib import Path
import pandas as pd

BASE = Path(PEDV_P+'')
OUT = BASE / 'Noncanonical_reclass' / 'tables'
OUT.mkdir(parents=True, exist_ok=True)
seq = ''.join(l.strip() for l in open(BASE / 'refbuild_GIb/virus/PEDV_GIb.fa') if not l.startswith('>')).upper()
STOP = {'TAA', 'TAG', 'TGA'}
CDS = {'ORF1a': (293, 12622), 'ORF1b': (12592, 20613), 'S': (20610, 24758), 'ORF3': (24758, 25033),
       'E': (25361, 25600), 'M': (25608, 26288), 'N': (26300, 27625)}


def orfs(mol):
    out = []
    for i in range(len(mol) - 2):
        if mol[i:i + 3] != 'ATG':
            continue
        j, n = i + 3, 1
        while j + 3 <= len(mol):
            if mol[j:j + 3] in STOP:
                out.append((i, n))
                break
            j += 3
            n += 1
    return out


ev = pd.read_csv(BASE / 'GIb_rnaseq_reanalysis/figures/Fig2_work/events_main.tsv', sep='\t')
dep = ev[ev.cls == 'TRSL_dep_noncanonical'].copy().reset_index(drop=True)
assert len(dep) == 1173 and dep.reads.sum() == 32646
new = []
for r in dep.itertuples():
    if r.frame in ('in', 'out'):
        new.append('truncated (counted)' if r.frame == 'in' else 'frameshift (counted)')
        continue
    D, A = int(r.donor), int(r.acceptor)
    mol = seq[:D] + seq[A - 1:]
    tags = set()
    for i, _n in orfs(mol):
        g = (i + 1) if i < D else (i - D + A)
        hit = [k for k, (s, t) in CDS.items() if s <= g <= t]
        if hit:
            k = max(hit, key=lambda k: CDS[k][0])
            tags.add('canon_inframe' if (g - CDS[k][0]) % 3 == 0 else 'canon_outframe')
        else:
            tags.add('noncds')
    if 'canon_inframe' in tags:
        new.append('translates_canonical (excluded)')
    elif tags:
        new.append('out-of-frame ORF (counted)')
    else:
        new.append('no ORF')
dep['new_class'] = new
dep['counted_in_plots'] = ~dep.new_class.str.startswith('translates_canonical') & (dep.new_class != 'no ORF')
dep[['donor', 'acceptor', 'reads', 'n_samples', 'frame', 'orf_consequence', 'new_class', 'counted_in_plots']].to_csv(
    OUT / 'TRSL_dep_noncanonical_reclassified.csv', index=False)
summ = dep.groupby('new_class').agg(events=('reads', 'size'), reads=('reads', 'sum')).reset_index()
summ.to_csv(OUT / 'reclass_summary.csv', index=False)
print(summ.to_string(index=False))
print('counted events', int(dep.counted_in_plots.sum()), 'reads', int(dep[dep.counted_in_plots].reads.sum()))
