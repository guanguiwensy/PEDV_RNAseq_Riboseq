"""Fig 6D -- POST-HOC REPRODUCTION SCRIPT (written after submission; NOT the original plotting program).

The original script of the submitted Fig 6D (+1-frame family map) was not found. This script redraws the panel from
`Supplementary_Figure_L/SourceData_PanelA_families_gt40.csv` (relative_frame == +1: 68 families, the number printed on
Fig 6B). No statistic is computed; the table was produced by comparative_genomics/04_retention_null_stats.

What is DATA: one box per +1-frame family spanning start_min .. stop_position (nt of the PEDV_GIb genome), colour by
`contains_AUG` (purple: family contains an AUG start; gold: near-cognate starts only), black tick at the stop codon.
What is INFERRED from the submitted figure (cosmetic only, no effect on data): (i) the two lanes alternate in genomic order,
(ii) red outlines mark ORF1a-x3, N-x1 and N-x2, (iii) the canonical-ORF bar on top uses the Fig 6 region colours.
Usage: python fig6D_reproduction_from_source_data.py [OUTDIR]
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
from matplotlib.patches import Rectangle

SRC = PEDV_K + "Supplementary_Figure_L/SourceData_PanelA_families_gt40.csv"
OUT = sys.argv[1] if len(sys.argv) > 1 else "out_fig6CD"
os.makedirs(OUT, exist_ok=True)
RCOL = {"ORF1a": "#0072B2", "ORF1b": "#56B4E9", "S": "#009E73", "ORF3": "#CC79A7", "E": "#D55E00", "M": "#E69F00", "N": "#117733"}
CDS = [("ORF1a", 293, 12622), ("ORF1b", 12592, 20613), ("S", 20610, 24758), ("ORF3", 24758, 25033), ("E", 25361, 25600), ("M", 25608, 26288), ("N", 26300, 27625)]
HILITE = {"ORF1a-x3", "N-x1", "N-x2"}
PUR, GOLD = "#5e3c99", "#a6761d"

fams = [r for r in csv.DictReader(open(SRC, encoding="utf-8-sig")) if r["relative_frame"] in ("+1", "1", "1.0")]
fams.sort(key=lambda r: int(r["start_min"]))
n_aug = sum(r["contains_AUG"] == "yes" for r in fams)
assert len(fams) == 68, "expected 68 +1 families (Fig 6B), got %d" % len(fams)

plt.rcParams.update({"font.family": ["Arial", "DejaVu Sans"], "pdf.fonttype": 42, "svg.fonttype": "none"})
fig, ax = plt.subplots(figsize=(12.4, 2.6))
for nm, a, b in CDS:
    ax.add_patch(Rectangle((a, 3.0), b - a, 0.6, color=RCOL[nm], lw=0))
    if nm in ("ORF1a", "ORF1b", "S", "M", "N"):
        ax.text((a + b) / 2, 3.3, nm, ha="center", va="center", color="white", fontsize=11, fontweight="bold")
    elif nm in ("ORF3", "E"):
        ax.text(a + (b - a) / 2 + (-150 if nm == "ORF3" else 150), 3.75, nm, ha=("right" if nm == "ORF3" else "left"), va="bottom", color=RCOL[nm], fontsize=11)
ax.text(-300, 1.35, "+1 frame", ha="right", va="center", fontsize=12)
rows = []
for i, r in enumerate(fams):
    lane = 1 if i % 2 == 0 else 0          # alternating lanes (inferred layout)
    y = 1.6 - 0.9 * lane
    a, s = int(r["start_min"]), int(r["stop_position"])
    col = PUR if r["contains_AUG"] == "yes" else GOLD
    hl = r["figure6_name"] in HILITE
    ax.add_patch(Rectangle((a, y - 0.3), s - a, 0.6, facecolor=col, edgecolor="#d62728" if hl else "none", lw=1.5 if hl else 0, zorder=3))
    ax.plot([s, s], [y - 0.4, y + 0.4], color="k", lw=1.0, zorder=4)
    rows.append([r["figure6_name"], r["family_id"], a, s, r["contains_AUG"], lane, "yes" if hl else "no"])
ax.add_patch(Rectangle((11000, -0.55), 500, 0.3, color=PUR, lw=0)); ax.text(11700, -0.4, "family contains AUG start", va="center", fontsize=12)
ax.add_patch(Rectangle((17500, -0.55), 500, 0.3, color=GOLD, lw=0)); ax.text(18200, -0.4, "near-cognate starts only", va="center", fontsize=12)
ax.set_xlim(-300, 28100); ax.set_ylim(-0.8, 4.2); ax.axis("off")
fig.tight_layout()
for ext in ("pdf", "svg", "png"):
    fig.savefig(os.path.join(OUT, "Fig6D_reproduction." + ext), dpi=300)
with open(os.path.join(OUT, "Fig6D_reproduction_source_data.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f); w.writerow(["figure6_name", "family_id", "start_min", "stop_position", "contains_AUG", "lane(inferred)", "red_outline(inferred)"]); w.writerows(rows)
print("Fig 6D: %d +1-frame families; %d contain AUG, %d near-cognate only" % (len(fams), n_aug, len(fams) - n_aug))
