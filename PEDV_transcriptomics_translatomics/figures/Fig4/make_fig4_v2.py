#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 4 v2：样本 = CHX 5/8 h (RIBO_VERO_CHX_1MOI_5h, RIBO_VERO_CHX_1MOI_8h)、CHX 16 h = RIBO_VERO_CHX_1MOI_16h（第二批）、HAR 5/8/16 h (RIBO_VERO_HAR_1MOI_5h/62/64)；
读长集合 = 每个文库中通过宿主判据的读长（25–35 nt；offset high、≥5,000 CDS-body P-site、F0≥0.5、最高帧/次高帧≥1.5）。
数据全部来自已同步到本地的 RGP_v1.3 / 服务器一遍统计，不上服务器、不重算 offset（RIBO_VERO_CHX_1MOI_16h 的 offset 由已冻结的 Module 1 方法 B 判定函数逐行复刻）。"""
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
import os, sys, csv, gzip, math, json, collections, datetime
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.dirname(HERE); TASK = os.path.dirname(OUT)
RGP = f"{TASK}/riboseq_GIb_mapping_psite"; WORK = f"{RGP}/_work/server"; FRM = f"{TASK}/read_length_frame_periodicity"
FIG = f"{OUT}/figures"; SD = f"{OUT}/source_data"; os.makedirs(FIG, exist_ok=True); os.makedirs(SD, exist_ok=True)
MM = 1 / 25.4; NOW = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.linewidth": 0.7, "svg.fonttype": "none", "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
                     "xtick.labelsize": 7, "ytick.labelsize": 7})
rdt = lambda p: pd.read_csv(p, sep="\t", comment="#", low_memory=False)
SAMPLES = [(SR.lib("RIBO_VERO_CHX_1MOI_5h"), "CHX", 5), (SR.lib("RIBO_VERO_CHX_1MOI_8h"), "CHX", 8), (SR.lib("RIBO_VERO_CHX_1MOI_16h"), "CHX", 16),
           (SR.lib("RIBO_VERO_HAR_1MOI_5h"), "HAR", 5), (SR.lib("RIBO_VERO_HAR_1MOI_8h"), "HAR", 8), (SR.lib("RIBO_VERO_HAR_1MOI_16h"), "HAR", 16)]
IDS = [s for s, _, _ in SAMPLES]; DRUG = {s: d for s, d, _ in SAMPLES}; TIME = {s: t for s, _, t in SAMPLES}
PAL = {"CHX": "#000000", "HAR": "#D55E00", "RNA": "#7F7F7F", "host": "#009E73", "PEDV": "#0072B2"}
CCOL = {"CDS": "#08519C", "5UTR": "#9ECAE1", "3UTR": "#E69F00", "other_exon": "#D8C9A3"}
CNAME = {"CDS": "CDS", "5UTR": "5′UTR", "3UTR": "3′UTR", "other_exon": "Other exon / ncRNA"}
GENES = [("5′UTR", 1, 292, "#EDEDED", 1), ("ORF1a", 293, 12622, "#F3E3EC", 1), ("ORF1b", 12592, 20613, "#E1EBF5", 2),
         ("S", 20610, 24758, "#E6F1E3", 1), ("ORF3", 24758, 25033, "#FBF0D9", 1), ("E", 25361, 25600, "#EAE3F3", 1),
         ("M", 25608, 26288, "#F8E3E3", 1), ("N", 26300, 27625, "#E0F0F0", 1), ("3′UTR", 27626, 27959, "#EDEDED", 1)]
GLEN = 27959
J_LO, J_HI = 7, 19

# ============================================================ 1. 读长集合（宿主判据）
main = pd.read_csv(f"{FRM}/read_length_frame_counts_and_fractions.csv", comment="#")
SR.assert_table_whitelisted(main, "sample", "read_length_frame_counts"); SR.assert_exactly(set(main["sample"]), SR.RIBO_VERO_ALIASES, "read_length_frame_counts")
B = main[main.region == "CDS_body_45_45"].copy()   # all six whitelisted libraries come from the same (whitelist-only) table
B = B[B.read_length.between(25, 35)]
SETL = {s: sorted(B[(B["sample"] == s) & B.host_criteria_pass].read_length.astype(int)) for s in IDS}
OFFS = {}
for s in IDS:
    o = B[B["sample"] == s].set_index("read_length").psite_offset_applied
    OFFS[s] = {int(L): int(v) for L, v in o.items() if not pd.isna(v)}
print("length sets:", SETL)
pd.DataFrame([dict(sample=s, treatment=DRUG[s], time_hpi=TIME[s], lengths=";".join(map(str, SETL[s])), offsets=";".join(f"{L}:{OFFS[s][L]}" for L in SETL[s]))
              for s in IDS]).to_csv(f"{SD}/length_sets_host_criteria.csv", index=False)

# ============================================================ 2. RIBO_VERO_CHX_1MOI_16h 的 offset（Module 1 方法 B，只读 stats）
def peak_stats(x, L):
    x = x[1:L + 1]; n = int(x.sum())
    if n == 0: return None, 0, 0, 0.0, 0.0
    w = x[J_LO - 1:J_HI]; k = int(np.argmax(w)); srt = np.sort(w); pk = int(x[J_LO + k - 1])
    return J_LO + k - 1, n, pk, pk / n, (srt[-1] / srt[-2] if srt[-2] > 0 else np.inf)
def offsets_from_stats(sid):
    d = f"{WORK}/stats/{sid}"
    hs = collections.defaultdict(lambda: np.zeros(600, np.int64))
    for r in csv.DictReader(open(f"{d}/host_startpos.tsv"), delimiter="\t"): hs[int(r["L"])][int(r["j"])] += int(r["n"])
    rl = {int(r["L"]): r for r in csv.DictReader(open(f"{d}/readlen_aln.tsv"), delimiter="\t")}
    off = {}
    for L in rl:
        o = peak_stats(hs[L], L)[0]; off[L] = o
    return off, rl
OFF813_ALL, RL813 = offsets_from_stats(SR.lib("RIBO_VERO_CHX_1MOI_16h"))
for L in SETL[SR.lib("RIBO_VERO_CHX_1MOI_16h")]: assert OFF813_ALL[L] == OFFS[SR.lib("RIBO_VERO_CHX_1MOI_16h")][L], (L, OFF813_ALL[L], OFFS[SR.lib("RIBO_VERO_CHX_1MOI_16h")][L])

# ============================================================ 3. Panel A：读长分布
dist = rdt(f"{RGP}/riboseq_readlength_distribution.tsv"); dist["L"] = pd.to_numeric(dist.aligned_length, errors="coerce")
def length_series(sid, kind):
    x = np.arange(18, 51); y = np.zeros(len(x))
    if sid == SR.lib("RIBO_VERO_CHX_1MOI_16h"):
        tot = sum(int(r["host_n" if kind == "host" else "pedv_plus_n"]) for r in RL813.values())
        for L, r in RL813.items():
            if 18 <= L <= 50: y[L - 18] = int(r["host_n" if kind == "host" else "pedv_plus_n"]) / tot
    else:
        d = dist[dist.sample_id == sid]; col = "host_reads" if kind == "host" else "PEDV_plus_reads"
        tot = d[col].sum()
        for _, r in d.dropna(subset=["L"]).iterrows():
            if 18 <= int(r.L) <= 50: y[int(r.L) - 18] = r[col] / tot
    return x, y
rowsA = []
def draw_A(axs):
    for ax, drug in zip(axs, ("CHX", "HAR")):
        for key, nm in (("host", "Host"), ("PEDV", "PEDV (+)")):
            ys = []
            for s in [s for s in IDS if DRUG[s] == drug]:
                x, y = length_series(s, "host" if key == "host" else "pedv"); ys.append(y)
                ax.plot(x, y, color=PAL[key], lw=0.6, alpha=0.55)
                for L, v in zip(x, y): rowsA.append(dict(sample=s, treatment=drug, time_hpi=TIME[s], molecule=nm, aligned_length=int(L), fraction=v))
            ax.plot(x, np.mean(ys, 0), color=PAL[key], lw=1.6, label=nm)
        ax.set_xlim(18, 50); ax.set_xticks([20, 25, 30, 35, 40, 45, 50]); ax.set_title(drug, fontsize=8, color=PAL[drug], fontweight="bold")
        ax.set_xlabel("Aligned footprint length (nt)", fontsize=8)
    axs[0].set_ylabel("Fraction of reads", fontsize=8)
    ym = max(a.get_ylim()[1] for a in axs)
    for a in axs: a.set_ylim(0, ym)
    axs[1].legend(frameon=False, fontsize=7, loc="upper right")

# ============================================================ 4. Panel B：三读框（宿主 CDS-body，所选读长）
rna = rdt(f"{TASK}/fig4_riboseq_qc/figures/Fig4C1_frame_fractions_source_data.csv") if False else pd.read_csv(f"{TASK}/fig4_riboseq_qc/figures/Fig4C1_frame_fractions_source_data.csv", comment="#")
rna = rna[(rna.category == "RNA-seq") & rna.sample_id.isin(SR.WHITELIST_LIBRARIES)]   # RNA-seq rows of the final sample table only
SR.assert_whitelisted(rna.sample_id, "Fig4B RNA-seq frame table")
Bset = B.set_index(["sample", "read_length"])
rowsB = []
FR = {}
for s in IDS:
    N = np.zeros(3)
    for L in SETL[s]: N += Bset.loc[(s, L), ["N0", "N1", "N2"]].values.astype(float)
    FR[s] = N / N.sum()
    rowsB.append(dict(sample=s, treatment=DRUG[s], time_hpi=TIME[s], lengths=";".join(map(str, SETL[s])), n_CDS_body_psites=int(N.sum()), frame0=FR[s][0], frame1=FR[s][1], frame2=FR[s][2]))
RNAF = rna[["frame0", "frame1", "frame2"]].values.astype(float)
def draw_B(ax):
    cats = [("CHX", np.array([FR[s] for s in IDS if DRUG[s] == "CHX"])), ("HAR", np.array([FR[s] for s in IDS if DRUG[s] == "HAR"])), ("RNA", RNAF)]
    w = 0.25
    for j, (c, a) in enumerate(cats):
        for i in range(3):
            x = i + (j - 1) * w
            ax.bar(x, a[:, i].mean(), width=w * 0.9, color=PAL[c], alpha=0.25 if c != "CHX" else 0.18, edgecolor=PAL[c], lw=0.7)
            ax.scatter(np.full(len(a), x) + np.linspace(-0.06, 0.06, len(a)), a[:, i], s=7, color=PAL[c], zorder=3, lw=0)
    ax.axhline(1 / 3, color="#8C8C8C", lw=0.7, ls="--"); ax.set_ylim(0, 1)
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["Frame 0", "Frame 1", "Frame 2"]); ax.set_ylabel("Fraction of host CDS P-sites", fontsize=8)
    ax.legend([Line2D([], [], marker="s", ls="", color=PAL[c]) for c in ("CHX", "HAR", "RNA")],
              [f"CHX Ribo (n={len(cats[0][1])})", f"HAR Ribo (n={len(cats[1][1])})", f"RNA-seq (n={len(RNAF)})"], frameon=False, fontsize=7, loc="upper right")

# ============================================================ 5. Panel C：起始/终止 metagene（所选读长，每 10^6 个 CDS P-site）
rowsC = []; MGS = {}
for s in IDS:
    mg = rdt(f"{WORK}/host_psite/{s}/host_psite_metagene.tsv"); g = rdt(f"{WORK}/host_psite/{s}/host_gene_cds_psite.tsv")
    Ls = SETL[s]; den = g[g.L.isin(Ls)][["f0", "f1", "f2"]].values.sum()
    for an in ("start", "stop"):
        m = mg[(mg.anchor == an) & mg.L.isin(Ls)].groupby("rel").n.sum().reindex(range(-50, 51), fill_value=0)
        MGS[(s, an)] = m.values / den * 1e6
        for r_, v in m.items(): rowsC.append(dict(sample=s, treatment=DRUG[s], time_hpi=TIME[s], anchor=an, rel_nt=r_, psite_count=int(v), per_million_set_CDS_psites=v / den * 1e6, denominator_CDS_psites=int(den)))
# 校验：同一函数复现旧 strict_core 版本（060）
old = rdt(f"{RGP}/host_metagene_GIb.tsv"); o60 = old[(old.sample_id == SR.lib("RIBO_VERO_CHX_1MOI_5h")) & (old.length_group == "strict_core") & (old.anchor == "start")].set_index("rel_nt").per_million_strict_core_CDS_psites
mg = rdt(f"{WORK}/host_psite/{SR.lib('RIBO_VERO_CHX_1MOI_5h')}/host_psite_metagene.tsv"); g = rdt(f"{WORK}/host_psite/{SR.lib('RIBO_VERO_CHX_1MOI_5h')}/host_gene_cds_psite.tsv")
chk = mg[(mg.anchor == "start") & mg.L.isin([29, 30])].groupby("rel").n.sum() / g[g.L.isin([29, 30])][["f0", "f1", "f2"]].values.sum() * 1e6
common = [i for i in range(-50, 51) if i in o60.index and i in chk.index]
print("metagene reproduction check (strict_core 060): max abs diff =", float(np.max(np.abs(np.array([o60[i] for i in common]) - np.array([chk[i] for i in common])))))
def draw_C(axs):
    xs = np.arange(-50, 51)
    for c, an in enumerate(("start", "stop")):
        top = 0
        for r_, drug in enumerate(("CHX", "HAR")):
            ax = axs[r_][c]; ys = []
            for s in [s for s in IDS if DRUG[s] == drug]:
                y = MGS[(s, an)]; ys.append(y); ax.plot(xs, y, color=PAL[drug], lw=0.6, alpha=0.5)
            ax.plot(xs, np.mean(ys, 0), color=PAL[drug], lw=1.4)
            ax.axvline(0, color="#8C8C8C", lw=0.7, ls="--"); ax.set_xlim(-50, 50); ax.set_xticks([-50, -25, 0, 25, 50])
            top = max(top, ax.get_ylim()[1])
            ax.text(0.04, 0.95, drug, transform=ax.transAxes, fontsize=7.5, color=PAL[drug], fontweight="bold", va="top")
            if c == 0: ax.set_ylabel("P-sites / 10⁶", fontsize=8)
            if r_ == 1: ax.set_xlabel("From start codon (nt)" if an == "start" else "From stop codon (nt)", fontsize=8)
            else: ax.set_xticklabels([])
        for r_ in range(2): axs[r_][c].set_ylim(0, top)

# ============================================================ 6. Panel D：宿主基因组注释（全部唯一比对 reads 的 P-site，外显子内重新归一）
Uall = {r["sample_id"]: r for r in csv.DictReader(open(f"{TASK}/riboseq_qc_figs_R1_R3/_work/host_annotation/host_annotation_unstranded_all.tsv"), delimiter="\t")}
u813 = list(csv.DictReader(open(f"{TASK}/riboseq_qc_figs_R1_R3/_work/host_annotation/{SR.lib('RIBO_VERO_CHX_1MOI_16h')}.host_annotation_unstranded.tsv"), delimiter="\t"))[0]
Uall[SR.lib("RIBO_VERO_CHX_1MOI_16h")] = u813
Uall = {k: v for k, v in Uall.items() if k in SR.WHITELIST_LIBRARIES}   # final sample table only
rman = {r["sample_id"]: r for r in csv.DictReader(open(f"{TASK}/fig4_riboseq_qc/_work/rna_manifest.tsv"), delimiter="\t") if r["sample_id"] in SR.WHITELIST_LIBRARIES}   # final sample table only
EXC = ["CDS", "5UTR", "3UTR", "other_exon"]
def exfrac(s):
    v = np.array([int(Uall[s][c]) for c in EXC], float); return v / v.sum()
rowsD = []
def draw_D(ax):
    x = 0; xs = []; labs = []
    for grp in ("Ribo", "RNA"):
        if grp == "Ribo":
            units = [(s, f"{TIME[s]}h", DRUG[s], exfrac(s)) for s in IDS]
        else:
            units = []
            for drug in ("CHX", "HAR"):
                for t in (5, 8, 16):
                    ss = [s for s in rman if rman[s]["treatment"] == drug and int(rman[s]["time_h"]) == t and s in Uall]
                    units.append((",".join(ss), f"{t}h", drug, np.mean([exfrac(s) for s in ss], 0)))
        for s, lab, drug, f in units:
            b = 0
            for c, v in zip(EXC, f):
                ax.bar(x, v, bottom=b, color=CCOL[c], width=0.8, edgecolor="white", lw=0.3); b += v
            rowsD.append(dict(assay=grp, sample=s, treatment=drug, time_hpi=lab, **{c: v for c, v in zip(EXC, f)}))
            xs.append(x); labs.append(lab); x += 1
        x += 0.8
    ax.set_xticks(xs); ax.set_xticklabels(labs, fontsize=7); ax.set_ylim(0, 1); ax.set_xlim(-0.7, x - 0.9)
    ax.set_ylabel("Fraction of P-sites", fontsize=8)
    ax.text(2.5, 1.03, "Ribo-seq", ha="center", va="bottom", fontsize=8); ax.text(9.9, 1.03, "RNA-seq", ha="center", va="bottom", fontsize=8, color=PAL["RNA"])
    for x0, x1, drug in ((-0.4, 2.4, "CHX"), (2.6, 5.4, "HAR"), (6.2, 9.0, "CHX"), (9.2, 12.0, "HAR")):
        ax.plot([x0, x1], [-0.17, -0.17], color=PAL[drug], lw=0.8, transform=ax.get_xaxis_transform(), clip_on=False)
        ax.text((x0 + x1) / 2, -0.19, drug, ha="center", va="top", fontsize=7.5, color=PAL[drug], transform=ax.get_xaxis_transform())
    ax.legend([Rectangle((0, 0), 1, 1, fc=CCOL[c]) for c in EXC], [CNAME[c] for c in EXC], frameon=False, fontsize=7, ncol=2, loc="upper center",
              bbox_to_anchor=(0.5, -0.36), handlelength=1.0, columnspacing=1.0)

# ============================================================ 7. Panel E：PEDV 占据（所选读长；另出全部读长版）
def psite_vec(sid, Ls=None, include_na=False):
    """+ 链逐碱基 P-site 计数（1..27959）。Ls=None → 全部有 offset 的读长。"""
    c = np.zeros(GLEN + 1)
    if sid != SR.lib("RIBO_VERO_CHX_1MOI_16h"):
        d = pd.read_csv(f"{RGP}/psite_counts_PEDV_GIb.tsv.gz", sep="\t", comment="#", low_memory=False) if not hasattr(psite_vec, "_d") else psite_vec._d
        psite_vec._d = d
        d = d[(d.sample_id == sid) & (d.strand == "+")]
        if Ls is not None: d = d[d.read_length.isin(Ls)]
        np.add.at(c, d.position.values, d.psite_count.values)
        if include_na:
            na = pd.read_csv(f"{RGP}/pedv_5p_offsetNA_reads_GIb.tsv.gz", sep="\t", comment="#"); na = na[(na.sample_id == sid) & (na.strand == "+")]
            np.add.at(c, na.pos5.values, na.n.values)
        return c[1:]
    off = OFF813_ALL; NA = 0
    with gzip.open(f"{WORK}/stats/{sid}/pedv_5p_matrix.tsv.gz", "rt") as f:
        next(f)
        for l in f:
            st, L, p5, sc, n = l.rstrip("\n").split("\t"); L = int(L); p5 = int(p5); n = int(n)
            if st != "+": continue
            o = off.get(min(L, 150))
            if o is None:
                if include_na and Ls is None: c[p5] += n
                continue
            if Ls is not None and L not in Ls: continue
            c[p5 + o] += n
    with gzip.open(f"{WORK}/stats/{sid}/pedv_spliced5p.tsv.gz", "rt") as f:
        next(f)
        for l in f:
            st, L, sc, coords, n = l.rstrip("\n").split("\t"); L = int(L); n = int(n)
            if st != "+": continue
            o = off.get(min(L, 150)); cs = coords.split(",")
            if o is None or o >= len(cs):
                if include_na and Ls is None: c[int(cs[0])] += n
                continue
            if Ls is not None and L not in Ls: continue
            c[int(cs[o])] += n
    return c[1:]
VE = {s: psite_vec(s, SETL[s]) for s in IDS}
VA = {s: psite_vec(s, None, True) for s in IDS}
# 校验：旧 6 库口径（全部读长 + offset=NA）的最大值与既有 Figure 4E 一致（060: 67,029）
print("max per-base (all lengths) 060 =", int(VA[SR.lib("RIBO_VERO_CHX_1MOI_5h")].max()), "(Fig4E: 67,029)")
def draw_E(fig, left, top, w, h, V, gap=2.0, label_unit=1e4):
    th = (h - 5 * gap - 15) / 6 ; ax_list = []
    for i, (s, drug, t) in enumerate(SAMPLES):
        y0 = top + i * (th + gap)
        ax = fig.add_axes([left / FW, 1 - (y0 + th) / FH, w / FW, th / FH]); ax_list.append(ax)
        v = V[s]; xs = np.arange(1, GLEN + 1); keep = v >= 1
        ax.vlines(xs[keep], 0.5, v[keep], colors=PAL[drug], linewidth=0.25)
        ax.set_yscale("log"); ax.set_ylim(1, 1e6); ax.set_xlim(1, GLEN + 1); ax.set_yticks([1, 1e2, 1e4, 1e6]); ax.minorticks_off()
        ax.set_yticklabels(["10⁰", "10²", "10⁴", "10⁶"], fontsize=7); ax.tick_params(axis="x", bottom=False, labelbottom=False)
        ax.text(0.5, 0.97, f"{drug} {t}h", transform=ax.transAxes, ha="center", va="top", fontsize=8.5)
        ax.text(1.005, 0.5, f"{v.max()/label_unit:.1f}", transform=ax.transAxes, fontsize=7.5, va="center", ha="left")
        if i == 0: ax.text(1.005, 1.02, "max (10⁴)", transform=ax.transAxes, fontsize=7.5, va="bottom", ha="right")
        if i == 3: ax.set_ylabel("P-sites per base (log₁₀)", fontsize=8, labelpad=2, y=1.0)
    gy = top + 6 * (th + gap) - gap + 2.0; gh = 4.0
    gax = fig.add_axes([left / FW, 1 - (gy + 2 * gh + 1) / FH, w / FW, (2 * gh) / FH]); gax.set_xlim(1, GLEN + 1); gax.set_ylim(0, 2 * gh)
    for sp in ("left", "right", "top"): gax.spines[sp].set_visible(False)
    gax.set_yticks([]); rowy = {1: gh, 2: 0}
    for g, s_, e_, col, r_ in GENES:
        gax.add_patch(Rectangle((s_, rowy[r_] + 0.4), e_ - s_ + 1, gh - 0.8, fc=col, ec="#777777", lw=0.4))
        if g in ("ORF1a", "ORF1b", "S", "N", "M"): gax.text((s_ + e_) / 2, rowy[r_] + gh / 2, g, ha="center", va="center", fontsize=7.5)
    gax.text(1, rowy[2] + gh / 2, "5′UTR", ha="left", va="center", fontsize=7); gax.text(25150, rowy[2] + gh / 2, "ORF3", ha="right", va="center", fontsize=7)
    gax.text(25480, rowy[2] + gh / 2, "E", ha="center", va="center", fontsize=7); gax.text(GLEN, rowy[2] + gh / 2, "3′UTR", ha="right", va="center", fontsize=7)
    gax.text(12300, rowy[2] + gh / 2, "31-nt overlap", ha="right", va="center", fontsize=7, color="#444444")
    gax.plot([12607, 12607], [rowy[2] + 0.4, rowy[1] + gh - 0.4], color="#444444", lw=0.6)
    gax.set_xticks([1, 5000, 10000, 15000, 20000, 25000]); gax.set_xticklabels(["1", "5,000", "10,000", "15,000", "20,000", "25,000"])
    gax.set_xlabel("Position on PEDV genome (nt, 1-based; 27,959 nt)", fontsize=8, labelpad=2)
    return ax_list

# ============================================================ 8. 组图（A–E，与现有 Fig 4 排版一致）
FW, FH = 180.0, 220.0
def combined(name, V):
    global FW, FH
    FW, FH = 180.0, 246.0
    fig = plt.figure(figsize=(FW * MM, FH * MM))
    def ax_mm(l, t, w, h): return fig.add_axes([l / FW, 1 - (t + h) / FH, w / FW, h / FH])
    for L_, T_, s_ in ((1, 1, "A"), (112, 1, "B"), (1, 62, "C"), (102, 62, "D"), (1, 134, "E")):
        fig.text(L_ / FW, 1 - T_ / FH, s_, fontsize=13, fontweight="bold", ha="left", va="top")
    draw_A([ax_mm(14, 7, 40, 38), ax_mm(62, 7, 40, 38)])
    draw_B(ax_mm(128, 7, 48, 38))
    draw_C([[ax_mm(22, 68, 30, 22), ax_mm(66, 68, 30, 22)], [ax_mm(22, 94, 30, 22), ax_mm(66, 94, 30, 22)]])
    draw_D(ax_mm(118, 70, 58, 40))
    draw_E(fig, 20, 142, 146, 96, V, gap=3.0)
    for ext, kw in (("pdf", {}), ("svg", {}), ("png", dict(dpi=300))):
        try: fig.savefig(f"{FIG}/{name}.{ext}", facecolor="white", **kw)
        except PermissionError: print("LOCKED", name, ext)
    plt.close(fig)
combined("Fig4_v2_combined", VE)
combined("Fig4_v2_combined_E_all_lengths", VA)

# ---- 单面板导出 ----
def single(name, w, h, drawfun, **kw):
    global FW, FH
    FW, FH = w, h
    fig = plt.figure(figsize=(w * MM, h * MM)); drawfun(fig)
    for ext in ("pdf", "svg", "png"):
        try: fig.savefig(f"{FIG}/{name}.{ext}", facecolor="white", dpi=300 if ext == "png" else None, bbox_inches="tight")
        except PermissionError: print("LOCKED", name, ext)
    plt.close(fig)
single("Fig4v2_A_length_distribution", 100, 45, lambda f: draw_A([f.add_axes([0.12, 0.2, 0.38, 0.65]), f.add_axes([0.58, 0.2, 0.38, 0.65])]))
single("Fig4v2_B_frame_fractions", 70, 50, lambda f: draw_B(f.add_axes([0.16, 0.14, 0.8, 0.78])))
single("Fig4v2_C_metagene", 80, 62, lambda f: draw_C([[f.add_axes([0.12, 0.56, 0.38, 0.38]), f.add_axes([0.6, 0.56, 0.38, 0.38])],
                                                       [f.add_axes([0.12, 0.12, 0.38, 0.38]), f.add_axes([0.6, 0.12, 0.38, 0.38])]]))
single("Fig4v2_D_genomic_annotation", 90, 60, lambda f: draw_D(f.add_axes([0.1, 0.34, 0.88, 0.58])))
single("Fig4v2_E_PEDV_occupancy", 170, 118, lambda f: draw_E(f, 14, 3, 140, 105, VE))
single("Fig4v2_E_PEDV_occupancy_all_lengths", 170, 118, lambda f: draw_E(f, 14, 3, 140, 105, VA))

# ---- 源数据 ----
def wr(name, df, title, extra=()):
    head = [f"# {title}", "# version: FIG4v2_1.0", f"# generated: {NOW} (local)",
            f"# libraries: CHX 5 h {SR.lib('RIBO_VERO_CHX_1MOI_5h')}, CHX 8 h {SR.lib('RIBO_VERO_CHX_1MOI_8h')}, CHX 16 h {SR.lib('RIBO_VERO_CHX_1MOI_16h')}, HAR 5 h {SR.lib('RIBO_VERO_HAR_1MOI_5h')}, HAR 8 h {SR.lib('RIBO_VERO_HAR_1MOI_8h')}, HAR 16 h {SR.lib('RIBO_VERO_HAR_1MOI_16h')} (one library per condition)",
            "# length sets: per library, 25-35 nt lengths meeting host criteria (offset confidence high, >=5000 CDS-body P-sites, F0>=0.5, top/second frame >=1.5) -- see length_sets_host_criteria.csv",
            "# host reference ChlSab1.1; virus reference PEDV_GIb (27,959 nt); coordinates 1-based", "# script: scripts/make_fig4_v2.py"] + [f"# {e}" for e in extra]
    with open(f"{SD}/{name}", "w", encoding="utf-8", newline="") as f:
        f.write("\n".join(head) + "\n"); df.to_csv(f, index=False, float_format="%.6g")
wr("panelA_length_distribution.csv", pd.DataFrame(rowsA), "Panel A: aligned footprint length distribution (host, PEDV (+)); denominator = all uniquely mapped reads of that molecule")
wr("panelB_frame_fractions.csv", pd.DataFrame(rowsB), "Panel B: host CDS-body (45/45) P-site frame fractions over the selected lengths; RNA-seq = pseudo P-site (Fig4C1 source data, unchanged)")
wr("panelC_metagene.csv", pd.DataFrame(rowsC), "Panel C: host start/stop codon P-site metagene over the selected lengths, per 10^6 host CDS P-sites of the same lengths")
wr("panelD_genomic_annotation.csv", pd.DataFrame(rowsD), "Panel D: genomic annotation of host P-sites (all uniquely mapped host reads, all lengths; exonic categories renormalised); RNA-seq = mean of the two libraries per condition")
vc = []
for s in IDS:
    for name, V in (("selected_lengths", VE), ("all_lengths_plus_offsetNA", VA)):
        v = V[s]; vc.append(dict(sample=s, treatment=DRUG[s], time_hpi=TIME[s], variant=name, total_psites=int(v.sum()), max_per_base=int(v.max()), argmax_position=int(v.argmax()) + 1))
wr("panelE_summary.csv", pd.DataFrame(vc), "Panel E: summary of PEDV per-base P-site tracks (+ strand)")
rowsE = []
for s in IDS:
    for p in range(GLEN): rowsE.append((s, DRUG[s], TIME[s], p + 1, int(VE[s][p]), int(VA[s][p])))
wr("panelE_per_base_counts.csv", pd.DataFrame(rowsE, columns=["sample", "treatment", "time_hpi", "position", "psites_selected_lengths", "psites_all_lengths_plus_offsetNA"]),
   "Panel E: per-base (+ strand) PEDV P-site counts, selected lengths and all lengths")
json.dump({"length_sets": SETL, "B_frames": {s: FR[s].tolist() for s in IDS}}, open(f"{OUT}/_numbers.json", "w"), indent=1)
print(pd.DataFrame(vc).to_string()); print(pd.DataFrame(rowsB).round(3).to_string())
