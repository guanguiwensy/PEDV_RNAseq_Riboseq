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
import os, json, datetime
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.dirname(HERE); TASK = os.path.dirname(OUT)
TAB = f"{OUT}/tables"; FIG = f"{OUT}/figures"; os.makedirs(FIG, exist_ok=True)
RGP = f"{TASK}/riboseq_GIb_mapping_psite"; SH = f"{TASK}/shared_riboseq_metrics"; SRV = f"{RGP}/_work/server/host_psite"
MM = 1 / 25.4; NOW = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
rd = lambda p: pd.read_csv(p, sep="\t", comment="#", low_memory=False)
SAMPLES = [(SR.lib("RIBO_VERO_CHX_1MOI_5h"), "CHX", 5), (SR.lib("RIBO_VERO_CHX_1MOI_8h"), "CHX", 8), (SR.lib("RIBO_VERO_CHX_1MOI_16h"), "CHX", 16),
           (SR.lib("RIBO_VERO_HAR_1MOI_5h"), "HAR", 5), (SR.lib("RIBO_VERO_HAR_1MOI_8h"), "HAR", 8), (SR.lib("RIBO_VERO_HAR_1MOI_16h"), "HAR", 16)]
DRUGC = {"CHX": "#000000", "HAR": "#D55E00"}; FC = ["#7B3294", "#56B4E9", "#999999"]; SC = "#F0E442"
LS = list(range(25, 36))
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.linewidth": 0.7, "svg.fonttype": "none", "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
                     "xtick.labelsize": 7, "ytick.labelsize": 7})

D = pd.read_pickle(f"{TAB}/_frame_counts.pkl")
off = rd(f"{RGP}/psite_offset_GIb.tsv")
q = rd(f"{RGP}/length_qc_GIb.tsv"); q["L"] = pd.to_numeric(q.read_length, errors="coerce")
q = q.dropna(subset=["L"]).astype({"L": int})
sets = rd(f"{RGP}/readlength_sets_GIb.tsv")
dist = rd(f"{SH}/length_metrics_long.tsv")
aln = rd(f"{RGP}/riboseq_GIb_alignment_summary.tsv").set_index("sample_id")


MIN_N = 5000; MIN_F0 = 0.50; MIN_RATIO = 1.5
body = D[D.region == "CDS_body_45_45"].copy()
body["offset_high"] = body.offset_confidence == "high"
body["host_criteria_pass"] = body.offset_high & (body.total_host_CDS_psites >= MIN_N) & (body.F0 >= MIN_F0) & (body.top_to_second_frame_ratio >= MIN_RATIO)
qq = q[["sample_id", "L", "QC_pass_C1_C2_C3", "C2_periodicity", "C3_frame0", "PEDV_periodicity_ratio(C2指标)", "PEDV_correct_frame_enrichment"]]
qq.columns = ["sample", "read_length", "PEDV_QC_pass_C1_C2_C3_frozen", "PEDV_C2_periodicity_frozen", "PEDV_C3_frame0_frozen", "PEDV_periodicity_ratio", "PEDV_frame0_enrichment"]
body = body.merge(qq, on=["sample", "read_length"], how="left")
body["in_strict_core"] = body.length_set_category == "strict_core"
body["extended_candidate"] = body.host_criteria_pass & body.PEDV_QC_pass_C1_C2_C3_frozen.fillna(False).astype(bool) | body.in_strict_core
flags = body[["sample", "read_length", "host_criteria_pass", "PEDV_QC_pass_C1_C2_C3_frozen", "PEDV_C2_periodicity_frozen", "PEDV_C3_frame0_frozen",
              "in_strict_core", "extended_candidate"]]


amb = pd.read_csv(f"{TAB}/canonical_cds_frame_ambiguity.csv"); amb_ids = set(amb[amb.ambiguous_nt > 0].transcript_id)
sens = []
for sid, drug, t in SAMPLES:
    g = rd(f"{SRV}/{sid}/host_gene_cds_psite.tsv")
    for L in range(20, 41):
        a = g[g.L == L]; b = a[~a.transcript_id.isin(amb_ids)]
        n_all = (a.f0 + a.f1 + a.f2).sum(); n_b = (b.f0 + b.f1 + b.f2).sum()
        sens.append(dict(sample=sid, read_length=L, F0_full_CDS=a.f0.sum() / n_all if n_all else np.nan,
                         F0_full_CDS_excl_ambiguous_transcripts=b.f0.sum() / n_b if n_b else np.nan, psites_removed=int(n_all - n_b)))
sens = pd.DataFrame(sens); sens["abs_change"] = (sens.F0_full_CDS - sens.F0_full_CDS_excl_ambiguous_transcripts).abs()


out = D.merge(flags, on=["sample", "read_length"], how="left")
hd = ["# Reading-frame distribution of host CDS P-sites by sample and aligned read length",
      "# version: FRM_v1.0", f"# generated: {NOW} (local)",
      "# input: riboseq_GIb_mapping_psite/_work/server/host_psite/<sample>/{host_gene_cds_psite,host_psite_metagene}.tsv (RGP_v1.3 server pass host_psite_pass.py, GI-b BAM, STAR NH=1); psite_offset_GIb.tsv; length_qc_GIb.tsv",
      "# P-site = alignment 5' end + calibrated offset of the same library and length (no re-estimation); host ChlSab1.1 Ensembl canonical protein_coding (19,165 transcripts); frame = (P-site transcript coordinate - CDS start) mod 3; exon junctions handled in spliced transcript coordinates; strand handled",
      "# regions: full_CDS = all GTF CDS positions (stop codon excluded); CDS_body_45_45 = full CDS minus first 45 and last 45 nt (project definition, = periodicity_host); CDS_body_45_30 = minus first 45 and last 30 nt",
      "# CDS_body is obtained by subtracting the start (rel 0..44) and stop (rel -45..-1 / -30..-1) windows of the pooled metagene; stop-window frame assumes CDS length = 0 mod 3 (fraction of P-sites in transcripts violating this is given)",
      "# F0 CI: Wilson 95% (binomial; assumes independent P-sites, too narrow) and transcript-level bootstrap 95% (full_CDS only; 1000 resamples)",
      "# host_criteria_pass (body): offset_confidence==high & >=5000 body P-sites & F0>=0.50 & top/second frame ratio>=1.5; PEDV QC flags are the frozen C1-C3 values; no reads are removed by any flag"]
with open(f"{OUT}/read_length_frame_counts_and_fractions.csv", "w", encoding="utf-8", newline="") as f:
    f.write("\n".join(hd) + "\n"); out.to_csv(f, index=False, float_format="%.6g")


body_idx = body.set_index(["sample", "read_length"])
full_idx = D[D.region == "full_CDS"].set_index(["sample", "read_length"])
rows = []; setdefs = {}
for sid, drug, t in SAMPLES:
    sc = sorted(off[(off.sample_id == sid) & (off.category_GIb == "strict_core")].read_length.astype(int))
    ext = sorted(body[(body["sample"] == sid) & body.extended_candidate & body.read_length.between(25, 35)].read_length)
    hostonly = sorted(body[(body["sample"] == sid) & body.host_criteria_pass].read_length)
    exp_ = sorted(sets[(sets.sample_id == sid) & (sets.in_expanded.astype(str).str.lower() == "true")].read_length.astype(int)) if "in_expanded" in sets.columns else []
    S = {"strict_core (current main set)": sc, "fixed_29_30 (uniform)": [29, 30], "extended_candidate (strict_core + host & PEDV criteria)": ext,
         "host_criteria_only (reference)": hostonly, "frozen_expanded (RGP_v1.3)": exp_,
         "all_25_35 (reference, not proposed)": LS}
    dl = dist[(dist.sample_id == sid)]
    th = int(dl[dl.domain == "host"].reads.sum() + 0); tp = int(dl[dl.domain == "PEDV_plus"].reads.sum())

    th = int(aln.loc[sid, "host_unique_reads"]); tp = int(aln.loc[sid, "PEDV_plus_reads"])
    for nm, Ls in S.items():
        hr = int(dl[(dl.domain == "host") & dl.read_length.isin(Ls)].reads.sum()); pr = int(dl[(dl.domain == "PEDV_plus") & dl.read_length.isin(Ls)].reads.sum())
        nb = sum(int(body_idx.loc[(sid, L), "total_host_CDS_psites"]) for L in Ls if (sid, L) in body_idx.index)
        n0 = sum(int(body_idx.loc[(sid, L), "N0"]) for L in Ls if (sid, L) in body_idx.index)
        nf = sum(int(full_idx.loc[(sid, L), "total_host_CDS_psites"]) for L in Ls if (sid, L) in full_idx.index)
        nf0 = sum(int(full_idx.loc[(sid, L), "N0"]) for L in Ls if (sid, L) in full_idx.index)
        rows.append(dict(sample=sid, treatment=drug, time_hpi=t, length_set=nm, lengths=";".join(map(str, Ls)),
                         host_reads_retained=hr, host_unique_reads_total=th, host_retention_fraction=hr / th,
                         PEDV_plus_reads_retained=pr, PEDV_plus_unique_reads_total=tp, PEDV_retention_fraction=pr / tp,
                         host_CDS_body_psites=nb, host_CDS_body_F0=(n0 / nb if nb else np.nan),
                         host_full_CDS_psites=nf, host_full_CDS_F0=(nf0 / nf if nf else np.nan)))
R = pd.DataFrame(rows)
hd2 = ["# Read-length filter retention audit (no data in Figures 4-6 changed)", "# version: FRM_v1.0", f"# generated: {NOW} (local)",
       "# retention = uniquely mapped reads (host; PEDV + strand) at the listed aligned lengths / all uniquely mapped reads of that molecule (all lengths incl. >=150 nt bucket)",
       "# host_*_F0 = frame-0 fraction of host CDS P-sites in the listed lengths (CDS_body = project definition 45/45; full_CDS = all CDS)",
       "# extended_candidate = strict_core of that library plus lengths 25-35 nt passing host criteria (offset high, >=5000 body P-sites, F0>=0.50, top/second>=1.5) AND the frozen PEDV C1-C3 QC",
       "# these are audit sets only; they are NOT used for any figure or table of the main analysis"]
with open(f"{OUT}/read_length_filter_retention_audit.csv", "w", encoding="utf-8", newline="") as f:
    f.write("\n".join(hd2) + "\n"); R.to_csv(f, index=False, float_format="%.6g")


bb = body[body.read_length.between(25, 35)]
share = bb.pivot(index="sample", columns="read_length", values="total_host_CDS_psites")
share = share.div(share.sum(axis=1), axis=0)
ref = share.mean(axis=0)
f0 = bb.pivot(index="sample", columns="read_length", values="F0")
comp = pd.DataFrame({"observed_F0_25_35": (share * f0).sum(axis=1), "standardised_F0_common_composition": (ref * f0).sum(axis=1)})
comp.to_csv(f"{TAB}/composition_standardised_F0.csv")


def n_lbl(n): return f"{n/1e6:.1f}M" if n >= 1e6 else (f"{n/1e3:.0f}k" if n >= 1e3 else str(int(n)))
def panel_grid(name, mode):
    W, H = 185.0, 200.0
    fig = plt.figure(figsize=(W * MM, H * MM))
    PW, PH, LM, GX, TOP = 78.0, 36.0, 15.0, 11.0, 7.0
    for sid, drug, t in SAMPLES:
        r = [5, 8, 16].index(t); c = 0 if drug == "CHX" else 1
        l = LM + c * (PW + GX); tp = TOP + r * (PH + 24)
        ax = fig.add_axes([l / W, 1 - (tp + PH) / H, PW / W, PH / H])
        b = body[(body["sample"] == sid) & body.read_length.isin(LS)].sort_values("read_length")
        x = np.arange(len(LS))
        if mode == "grouped":
            for k, (col, fc) in enumerate(zip(["F0", "F1", "F2"], FC)):
                ax.bar(x + (k - 1) * 0.27, b[col].values, width=0.25, color=fc, linewidth=0, label=f"Frame {k}")
        else:
            ax.bar(x, b.F0.values, width=0.62, color=FC[0], linewidth=0)
            lo, hi = b.F0_CI95_wilson_low.values, b.F0_CI95_wilson_high.values
            ax.errorbar(x, b.F0.values, yerr=[b.F0.values - lo, hi - b.F0.values], fmt="none", ecolor="black", elinewidth=0.6, capsize=1.2, capthick=0.6)
        ax.axhline(1 / 3, color="#555555", lw=0.8, ls="--", zorder=0)
        ax.set_xlim(-0.7, len(LS) - 0.3); ax.set_ylim(0, 1.0); ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0]); ax.set_yticklabels(["0", "0.25", "0.5", "0.75", "1"])
        ax.set_xticks(x); ax.set_xticklabels([str(L) for L in LS])
        ax.set_title(f"{drug} {t} h", color=DRUGC[drug], fontweight="bold", fontsize=8, pad=2)
        if c == 0: ax.set_ylabel("Fraction of host CDS-body P-sites" if mode == "grouped" else "Frame-0 fraction (F0)", fontsize=7.5)
        else: ax.tick_params(labelleft=False)
        ax.set_xlabel("Aligned footprint length (nt)", fontsize=7.5, labelpad=1)

        for i, rr in enumerate(b.itertuples()):
            ax.text(i, -0.30, n_lbl(rr.total_host_CDS_psites), fontsize=7, ha="center", va="top", rotation=0, transform=ax.transData)
            if rr.in_strict_core:
                ax.plot([i], [0.965], ls="", marker="s", ms=4.5, mfc=SC, mec="black", mew=0.5, zorder=5)
            if rr.host_criteria_pass:
                ax.plot([i], [0.885], ls="", marker="o", ms=3.2, mfc="black", mec="black", zorder=5)
        ax.text(-0.75, -0.30, "n", fontsize=7, ha="right", va="top", style="italic")
    if mode == "grouped":
        hl = [Rectangle((0, 0), 1, 1, fc=c_) for c_ in FC]
        labs = ["Frame 0", "Frame 1", "Frame 2"]
    else:
        hl = [Rectangle((0, 0), 1, 1, fc=FC[0]), Line2D([], [], color="black", lw=0.6, marker="|", ms=5)]
        labs = ["Frame-0 fraction", "95% CI (Wilson)"]
    hl += [Line2D([], [], color="#555555", lw=0.8, ls="--"), Line2D([], [], ls="", marker="s", ms=4.5, mfc=SC, mec="black", mew=0.5),
           Line2D([], [], ls="", marker="o", ms=3.2, mfc="black", mec="black")]
    labs += ["1/3 (no frame preference)", "strict_core length of this library", "host criteria met: offset high, n ≥ 5,000, F0 ≥ 0.5, top/second ≥ 1.5"]
    fig.legend(hl, labs, loc="lower left", bbox_to_anchor=(LM / W, 0.0), ncol=2, frameon=False, fontsize=7, handlelength=1.4, columnspacing=1.2)
    for ext, kw in (("pdf", {}), ("svg", {}), ("png", dict(dpi=300))):
        try: fig.savefig(f"{OUT}/{name}.{ext}", facecolor="white", **kw)
        except PermissionError: print("LOCKED", name, ext)
    plt.close(fig)
panel_grid("read_length_frame_periodicity_grouped_bars", "grouped")
panel_grid("read_length_frame0_fraction_bars", "f0")


pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print("sensitivity (exclude %d ambiguous transcripts): max abs change in F0 (25-35 nt) = %.4f; max psites removed = %d" %
      (len(amb_ids), sens[sens.read_length.between(25, 35)].abs_change.max(), sens.psites_removed.max()))
sens.to_csv(f"{TAB}/ambiguity_sensitivity.csv", index=False, float_format="%.6g")
print(comp.round(3))
print(R[R.length_set.str.startswith(("strict", "fixed", "extended", "all"))].pivot(index="sample", columns="length_set", values=["host_retention_fraction", "PEDV_retention_fraction"]).round(3).to_string())
print(R[R.length_set.str.startswith(("strict", "extended", "fixed"))].pivot(index="sample", columns="length_set", values="host_CDS_body_F0").round(3).to_string())
print(body[body.extended_candidate & ~body.in_strict_core][["sample", "read_length", "F0", "total_host_CDS_psites", "offset_confidence"]].to_string())

d = body[body.read_length.between(27, 33)].groupby(["treatment", "read_length"]).F0.agg(["mean", "min", "max"]).round(2)
print(d.unstack(0).to_string())

fb = D[D.read_length.between(25, 35)].pivot_table(index=["sample", "read_length"], columns="region", values="F0")
print("max |F0 body - F0 full| =", (fb.CDS_body_45_45 - fb.full_CDS).abs().max().round(3), " max |45_30 - 45_45| =", (fb.CDS_body_45_30 - fb.CDS_body_45_45).abs().max().round(3))

ph = rd(f"{RGP}/periodicity_host.tsv"); ph["L"] = pd.to_numeric(ph.read_length, errors="coerce")
ph = ph[ph.L.between(25, 35)].drop(columns=["read_length"]).rename(columns={"sample_id": "sample", "L": "read_length"})
cmp_ = body.merge(ph[["sample", "read_length", "psite_frame0", "n_reads(5′端在CDS主体)"]], on=["sample", "read_length"])
cmp_["diff"] = cmp_.F0 - cmp_.psite_frame0
cmp_[["sample", "read_length", "F0", "psite_frame0", "diff", "total_host_CDS_psites", "n_reads(5′端在CDS主体)"]].to_csv(f"{TAB}/body_F0_vs_periodicity_host.csv", index=False, float_format="%.6g")
print("body F0 (true P-site, subtraction) vs periodicity_host (5'-end in body): max |diff| =", cmp_["diff"].abs().max().round(3), " mean diff =", cmp_["diff"].mean().round(4))
