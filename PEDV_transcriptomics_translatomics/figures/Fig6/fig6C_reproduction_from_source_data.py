"""Fig 6C -- POST-HOC REPRODUCTION SCRIPT (written after submission; NOT the original plotting program).

The original bar-plot script of the submitted Fig 6C was not found in the project archive. This script redraws
the panel from the stored source table `R6-E_panel_data_ORF3_included.csv` (column `fam_per_kb_gt40`:
alternative-frame ORF families > 40 aa per kb of parent CDS; +1 and +2 families, ORF3 included, 71 families).
The table itself is produced by figures/SupFig15/fig_E_onesided.py (statistics are NOT recomputed here).

Layout (colours, order, 1-decimal value labels, rotated region labels) was read off the submitted figure.
The final check compares the plotted values (1 decimal) with the six numbers printed on the submitted figure.
Usage: python fig6C_reproduction_from_source_data.py [OUTDIR]
"""
import os, sys, csv
_d = os.path.dirname(os.path.abspath(__file__))
while _d != os.path.dirname(_d) and not os.path.isdir(os.path.join(_d, "config")):
    _d = os.path.dirname(_d)
sys.path.insert(0, os.path.join(_d, "config"))
from pedv_paths import PEDV_K
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SRC = PEDV_K + "Result6_genomics_figures/source_data/R6-E_panel_data_ORF3_included.csv"
OUT = sys.argv[1] if len(sys.argv) > 1 else "out_fig6CD"
os.makedirs(OUT, exist_ok=True)
RCOL = {"ORF1a": "#0072B2", "ORF1b": "#56B4E9", "S": "#009E73", "ORF3": "#CC79A7", "M": "#E69F00", "N": "#117733"}
SUBMITTED = {"ORF1a": 2.6, "ORF1b": 2.2, "S": 2.9, "ORF3": 3.6, "M": 4.4, "N": 3.8}   # printed on the submitted Fig 6C
ORDER = ["ORF1a", "ORF1b", "S", "ORF3", "M", "N"]

rows = {r["region"]: r for r in csv.DictReader(open(SRC, encoding="utf-8-sig"))}
vals = [float(rows[k]["fam_per_kb_gt40"]) for k in ORDER]

plt.rcParams.update({"font.family": ["Arial", "DejaVu Sans"], "pdf.fonttype": 42, "svg.fonttype": "none"})
fig, ax = plt.subplots(figsize=(2.6, 3.6))
ax.bar(range(6), vals, color=[RCOL[k] for k in ORDER], width=0.85)
for i, v in enumerate(vals):
    ax.text(i, v + 0.05, "%.1f" % v, ha="center", va="bottom", fontsize=11)
ax.set_xticks(range(6)); ax.set_xticklabels(ORDER, rotation=90, fontsize=11)
ax.set_yticks([0, 2, 4]); ax.set_ylim(0, 5)
ax.set_ylabel("Families per kb", fontsize=11)
ax.set_title(">40 aa", fontsize=11)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("pdf", "svg", "png"):
    fig.savefig(os.path.join(OUT, "Fig6C_reproduction." + ext), dpi=300)

with open(os.path.join(OUT, "Fig6C_reproduction_source_data.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f); w.writerow(["region", "length_nt", "families_per_kb_gt40", "label_on_submitted_figure"])
    for k, v in zip(ORDER, vals): w.writerow([k, rows[k]["length_nt"], repr(v), SUBMITTED[k]])
bad = [(k, v, SUBMITTED[k]) for k, v in zip(ORDER, vals) if round(v, 1) != SUBMITTED[k]]
print("Fig 6C values:", dict(zip(ORDER, [round(v, 3) for v in vals])))
print("CHECK vs numbers printed on submitted figure:", "OK (6/6)" if not bad else "MISMATCH %s" % bad)
sys.exit(1 if bad else 0)
