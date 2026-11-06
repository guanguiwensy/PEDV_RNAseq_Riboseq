#!/usr/bin/env python3
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
import os, csv, datetime
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.dirname(HERE); TASK = os.path.dirname(OUT)
SH = f"{TASK}/shared_riboseq_metrics"; WORK = f"{TASK}/riboseq_GIb_mapping_psite/_work/server"
QC = f"{TASK}/batch23_riboseq_qc_mapping_output/batch23_riboseq_qc_mapping_summary.tsv"
FIG = f"{OUT}/figures"; os.makedirs(FIG, exist_ok=True)
MM = 1 / 25.4; NOW = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
rd = lambda p: pd.read_csv(p, sep="\t", comment="#", low_memory=False)
SAMPLES = [(SR.lib("RIBO_VERO_CHX_1MOI_5h"), "CHX", 5), (SR.lib("RIBO_VERO_CHX_1MOI_8h"), "CHX", 8), (SR.lib("RIBO_VERO_CHX_1MOI_16h"), "CHX", 16),
           (SR.lib("RIBO_VERO_HAR_1MOI_5h"), "HAR", 5), (SR.lib("RIBO_VERO_HAR_1MOI_8h"), "HAR", 8), (SR.lib("RIBO_VERO_HAR_1MOI_16h"), "HAR", 16)]
IDS = [s for s, _, _ in SAMPLES]; DRUG = {s: d for s, d, _ in SAMPLES}; TIME = {s: t for s, _, t in SAMPLES}
DRUGC = {"CHX": "#000000", "HAR": "#D55E00"}; HOST = "#009E73"; PEDV = "#0072B2"
CAT = [("Adapter/length-discarded", "adapter", "#D9D9D9"), ("rRNA/ncRNA (bowtie2)", "rRNA", "#8C8C8C"), ("Host unique", "host", HOST),
       ("PEDV unique", "pedv", PEDV), (">1 locus", "multi", "#CC79A7"), ("Unmapped: too short", "short", "#F0E442"), ("Unmapped: other", "other", "#56B4E9")]
MARK = {5: "o", 8: "s", 16: "^"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.linewidth": 0.7, "svg.fonttype": "none", "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
                     "xtick.labelsize": 7, "ytick.labelsize": 7})


F5 = rd(f"{SH}/sample_funnel.tsv").set_index("sample_id")
rowsA = []
for s in IDS:
    if s != SR.lib("RIBO_VERO_CHX_1MOI_16h"):
        r = F5.loc[s]
        n = dict(adapter=r.n_adapter_length_discarded, rRNA=r.n_rRNA_ncRNA_discarded, host=r.n_host_unique, pedv=r.n_PEDV_unique,
                 multi=r.n_more_than_one_locus, short=r.n_unmapped_too_short, other=r.n_unmapped_other)
        raw = r.raw_r1_reads; myco = r.myco_pct_all_of_nonrRNA
    else:
        q = rd(QC).set_index("seq_id").loc[s]
        raw, trim, nonr = int(q.raw_r1_reads), int(q.trimmed_reads), int(q.non_rRNA_reads)
        fin = {}
        for l in open(f"{WORK}/align/{s}/{s}.Log.final.out", encoding="utf-8"):
            if "|" in l: k, v = l.split("|", 1); fin[k.strip()] = v.strip()
        inp = int(fin["Number of input reads"]); uniq = int(fin["Uniquely mapped reads number"])
        mm_ = int(fin["Number of reads mapped to too many loci"]); short = int(fin["Number of reads unmapped: too short"])
        other = int(fin["Number of reads unmapped: other"]) + int(fin["Number of reads unmapped: too many mismatches"])
        rl = pd.read_csv(f"{WORK}/stats/{s}/readlen_aln.tsv", sep="\t")
        host, pedv = int(rl.host_n.sum()), int(rl.pedv_plus_n.sum() + rl.pedv_minus_n.sum())
        assert inp == nonr, (inp, nonr); assert host + pedv == uniq, (host, pedv, uniq)
        assert uniq + mm_ + short + other == inp, (uniq, mm_, short, other, inp)
        n = dict(adapter=raw - trim, rRNA=trim - nonr, host=host, pedv=pedv, multi=mm_, short=short, other=other)

        bt = open(f"{WORK}/myco/{s}.bt2.log", encoding="utf-8").read()
        one = int([l for l in bt.splitlines() if "aligned exactly 1 time" in l][0].split()[0]); more = int([l for l in bt.splitlines() if "aligned >1 times" in l][0].split()[0])
        myco = 100 * (one + more) / nonr
    tot = sum(n.values()); assert int(tot) == int(raw), (s, tot, raw)
    for name, key, _ in CAT: rowsA.append(dict(sample=s, treatment=DRUG[s], time_hpi=TIME[s], category=name, reads=int(n[key]), raw_r1_reads=int(raw), fraction_of_raw=n[key] / raw))
    rowsA.append(dict(sample=s, treatment=DRUG[s], time_hpi=TIME[s], category="Mycoplasma (% of non-rRNA reads, bowtie2 all alignments)", reads=np.nan, raw_r1_reads=int(raw), fraction_of_raw=myco / 100))
A = pd.DataFrame(rowsA)


LS = rd(f"{SH}/length_summary.tsv")
rows = []
for s in IDS:
    for dom in ("host", "PEDV_plus"):
        if s != SR.lib("RIBO_VERO_CHX_1MOI_16h"):
            r = LS[(LS.sample_id == s) & (LS.domain == dom)].iloc[0]; med, f2931, n = int(r.median_length), float(r.frac_29_31), int(r.n_unique_reads)
        else:
            rl = pd.read_csv(f"{WORK}/stats/{s}/readlen_aln.tsv", sep="\t").sort_values("L")
            col = "host_n" if dom == "host" else "pedv_plus_n"; v = rl[col].values.astype(float); L = rl.L.values; tot = v.sum()
            cum = np.cumsum(v) / tot; med = int(L[np.searchsorted(cum, 0.5)]); f2931 = v[(L >= 29) & (L <= 31)].sum() / tot; n = int(tot)
        rows.append(dict(sample=s, treatment=DRUG[s], time_hpi=TIME[s], molecule=("Host" if dom == "host" else "PEDV (+)"), n_unique_reads=n, median_length_nt=med, fraction_29_31nt=f2931))
BC = pd.DataFrame(rows)

def wr(name, df, title):
    head = [f"# {title}", "# version: SUP12v2_1.0", f"# generated: {NOW} (local)",
            f"# libraries: CHX 5 h {SR.lib('RIBO_VERO_CHX_1MOI_5h')}, CHX 8 h {SR.lib('RIBO_VERO_CHX_1MOI_8h')}, CHX 16 h {SR.lib('RIBO_VERO_CHX_1MOI_16h')}, HAR 5 h {SR.lib('RIBO_VERO_HAR_1MOI_5h')}, HAR 8 h {SR.lib('RIBO_VERO_HAR_1MOI_8h')}, HAR 16 h {SR.lib('RIBO_VERO_HAR_1MOI_16h')} (one library per condition)",
            "# denominators: fraction_of_raw = reads / raw R1 reads; length statistics use all uniquely mapped reads of the molecule (PEDV (+) strand), length = aligned block length",
            "# script: scripts/make_supfig12_v2_ABC.py"]
    with open(f"{FIG}/{name}", "w", encoding="utf-8", newline="") as f:
        f.write("\n".join(head) + "\n"); df.to_csv(f, index=False, float_format="%.6g")
wr("panelA_read_fate_source_data.csv", A, "Panel A: fate of raw R1 reads and mycoplasma estimate")
wr("panelBC_length_summary_source_data.csv", BC, "Panels B/C: median aligned length and fraction of footprints at 29-31 nt, host vs PEDV (+)")


W, H = 180.0, 142.0
fig = plt.figure(figsize=(W * MM, H * MM))
def axmm(l, t, w, h): return fig.add_axes([l / W, 1 - (t + h) / H, w / W, h / H])
fig.text(1 / W, 1 - 1 / H, "A", fontsize=13, fontweight="bold", va="top")
axA = axmm(22, 8, 128, 44)
for i, s in enumerate(IDS):
    left = 0
    for name, key, col in CAT:
        v = A[(A["sample"] == s) & (A.category == name)].fraction_of_raw.iloc[0]
        axA.barh(i, v, left=left, color=col, height=0.7, linewidth=0); left += v
axA.set_ylim(len(IDS) - 0.4, -0.6); axA.set_xlim(0, 1); axA.set_yticks(range(len(IDS)))
axA.set_yticklabels([f"{DRUG[s]} {TIME[s]} h" for s in IDS])
for t_, s in zip(axA.get_yticklabels(), IDS): t_.set_color(DRUGC[DRUG[s]])
axA.set_xlabel("Fraction of raw reads", fontsize=8); axA.tick_params(axis="y", labelsize=8)
mp = A[A.category.str.startswith("Mycoplasma")].set_index("sample").fraction_of_raw * 100
for i, s in enumerate(IDS): axA.text(1.02, i, f"{mp[s]:.3f}" if mp[s] < 1 else f"{mp[s]:.2f}", va="center", fontsize=7.5, transform=axA.get_yaxis_transform())
axA.text(1.02, -0.62, "Mycoplasma\n(% of non-rRNA)", fontsize=7, va="bottom", ha="left", transform=axA.get_yaxis_transform(), linespacing=1.1)
hl = [Rectangle((0, 0), 1, 1, fc=c) for _, _, c in CAT]
fig.legend(hl, [n for n, _, _ in CAT], loc="upper left", bbox_to_anchor=(22 / W, 1 - 64 / H), ncol=4, frameon=False, fontsize=7, handlelength=1.2, columnspacing=1.2)

fig.text(1 / W, 1 - 78 / H, "B", fontsize=13, fontweight="bold", va="top")
fig.text(92 / W, 1 - 78 / H, "C", fontsize=13, fontweight="bold", va="top")
def paired(ax, drug, col, ylim, yt, fmt, gap):
    labs = []
    for s in [s for s in IDS if DRUG[s] == drug]:
        h = BC[(BC["sample"] == s) & (BC.molecule == "Host")].iloc[0][col]; p = BC[(BC["sample"] == s) & (BC.molecule == "PEDV (+)")].iloc[0][col]
        dx = {5: -0.05, 8: 0.0, 16: 0.05}[TIME[s]]
        ax.plot([0 + dx, 1 + dx], [h, p], color="#8C8C8C", lw=0.9, zorder=1)
        ax.scatter([0 + dx], [h], marker=MARK[TIME[s]], s=26, color=HOST, zorder=3, lw=0)
        ax.scatter([1 + dx], [p], marker=MARK[TIME[s]], s=26, color=PEDV, zorder=3, lw=0)
        labs.append([p, fmt(p - h)])
    labs.sort(key=lambda z: z[0]); ys = [z[0] for z in labs]
    for i in range(1, len(ys)):
        if ys[i] - ys[i - 1] < gap: ys[i] = ys[i - 1] + gap
    for y_, z in zip(ys, labs): ax.text(1.18, y_, z[1], fontsize=7, va="center")
    ax.set_xlim(-0.3, 1.7); ax.set_ylim(*ylim); ax.set_yticks(yt)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Host", "PEDV"])
    for t_, c in zip(ax.get_xticklabels(), (HOST, PEDV)): t_.set_color(c)
    ax.set_title(f"{drug}\nΔ = PEDV − Host", fontsize=7.5, color=DRUGC[drug], pad=3, linespacing=1.1)
for k, drug in enumerate(("CHX", "HAR")):
    ax = axmm(14 + k * 36, 88, 28, 44)
    paired(ax, drug, "median_length_nt", (27, 33), [27, 28, 29, 30, 31, 32, 33], lambda d: f"{int(round(d)):+d}" if round(d) != 0 else "0", 0.35)
    if k == 0: ax.set_ylabel("Median length (nt)", fontsize=8)
    else: ax.tick_params(labelleft=False)
for k, drug in enumerate(("CHX", "HAR")):
    ax = axmm(104 + k * 30, 88, 23, 44)
    paired(ax, drug, "fraction_29_31nt", (0.3, 0.7), [0.3, 0.4, 0.5, 0.6, 0.7], lambda d: f"{d:+.2f}", 0.035)
    if k == 0: ax.set_ylabel("Fraction in 29–31 nt", fontsize=8)
    else: ax.tick_params(labelleft=False)
fig.legend([Line2D([], [], ls="", marker=MARK[t], color="#444444", ms=5.5) for t in (5, 8, 16)], ["5 h", "8 h", "16 h"], loc="center left",
           bbox_to_anchor=(160 / W, 1 - 110 / H), frameon=False, fontsize=7.5, handletextpad=0.3, labelspacing=0.8)
for ext, kw in (("pdf", {}), ("svg", {}), ("png", dict(dpi=300))):
    try: fig.savefig(f"{FIG}/SupFig12_v2_ABC.{ext}", facecolor="white", **kw)
    except PermissionError: print("LOCKED", ext)
plt.close(fig)
pd.set_option("display.width", 220)
print(A[A.category.isin(["Host unique", "PEDV unique", ">1 locus"]) | A.category.str.startswith("Myco")].pivot(index="sample", columns="category", values="fraction_of_raw").round(4))
print(BC.pivot(index="sample", columns="molecule", values=["median_length_nt", "fraction_29_31nt"]).round(3))
