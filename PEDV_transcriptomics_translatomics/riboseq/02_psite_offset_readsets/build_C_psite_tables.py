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
import os, sys, json, gzip, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

OFF = {sid: {int(k): v for k, v in d.items()} for sid, d in json.load(open(f"{OUT}/_work/offset_table_GIb.json")).items()}
SETROWS = [r for r in csv.DictReader((l for l in open(f"{OUT}/readlength_sets_GIb.tsv", encoding="utf-8") if not l.startswith("#")), delimiter="\t")]
SETS = collections.defaultdict(dict)
for r in SETROWS:
    if r["read_length"] == "≥150": continue
    SETS[r["sample_id"]][int(r["read_length"])] = r
def in_set(sid, L, k):
    r = SETS[sid].get(L)
    if k == "all_lengths": return OFF[sid].get(L) is not None
    if k == "LPF_ge35": return L >= LONG_MIN
    return bool(r) and r[f"in_{k}"] == "True"
SETNAMES = ["strict_core", "fixed_29_30", "expanded", "short_RPF", "LPF_ge35", "all_lengths"]

cover = np.zeros(VLEN + 2, np.int16)
for _, a, b in ORFS: cover[a:b + 1] += 1
def orf_frame_str(p):
    s = [f"{n}:{(p - a) % 3}" for n, a, b in ORFS if a <= p <= b]
    return ";".join(s) if s else "-"
OFS = {p: orf_frame_str(p) for p in range(1, VLEN + 1)}

pf = lambda f, off: np.array([f[(k - off) % 3] for k in range(3)])
f3 = lambda x: f"{x:.3f}"
long_rows_written = 0; NA_ROWS = []
wide = {}; checks = []; per_v = []; per_h = []; na_report = []
prov_extra = ["offset table: psite_offset_GIb.tsv (RGP_v1.0；与 Module 1 M1_v1.0 逐样本×逐读长 0 差异)",
              "P-site: + 链 = 比对 5′端 + offset；− 链 = 比对 5′端(reference_end) − offset；5′端 22 nt 内有剪接跳跃者沿比对坐标取第 offset+1 个碱基；"
              "5′ soft-clip 不计入（B 口径，比对 5′端）",
              "scripts: scripts/build_C_psite_tables.py；R1 scripts/gib_bam_stats.py（pedv_5p_matrix / pedv_spliced5p）"] + bam_prov()

LONG = f"{OUT}/psite_counts_PEDV_GIb.tsv.gz"
fo = gzip.open(LONG, "wt", encoding="utf-8", newline="")
fo.write(provenance("psite_counts_PEDV_GIb：sample × strand × position × read_length（碱基分辨率长表）", prov_extra +
         ["columns: genome_frame=(position-1)%3；orf_frame=该位点所在 canonical ORF 及相对 ORF 起点的读框（重叠区列多个；- 表示 ORF 外）；"
          "length_category=该样本该读长的冻结类别（分析口径，非过滤）"]))
fo.write("sample_id\ttime_h\tdrug\tstrand\tposition\tread_length\toffset\tlength_category\tpsite_count\tgenome_frame\torf_frame\n")
for sid in IDS:
    m = META[sid]; off = OFF[sid]
    C = collections.Counter(); CD = collections.Counter(); nNA = collections.Counter(); ntot = 0
    with gzip.open(f"{WORK}/stats/{sid}/pedv_5p_matrix.tsv.gz", "rt") as f:
        next(f)
        for l in f:
            st, L, p5, sc, n = l.rstrip("\n").split("\t"); L = int(L); p5 = int(p5); n = int(n); ntot += n
            o = off.get(min(L, 150))
            if o is None:
                nNA[L] += n; NA_ROWS.append({**mcols(sid), "strand": st, "read_length": L, "pos5": p5, "spliced_5p_coords": "-", "n": n}); continue
            p = p5 + o if st == "+" else p5 - o
            C[(st, p, L)] += n; CD[(st, L, p5, p)] = 1
    with gzip.open(f"{WORK}/stats/{sid}/pedv_spliced5p.tsv.gz", "rt") as f:
        next(f)
        for l in f:
            st, L, sc, coords, n = l.rstrip("\n").split("\t"); L = int(L); n = int(n); ntot += n
            o = off.get(min(L, 150)); cs = coords.split(",")
            if o is None or o >= len(cs):
                nNA[L] += n; NA_ROWS.append({**mcols(sid), "strand": st, "read_length": L, "pos5": cs[0], "spliced_5p_coords": coords, "n": n}); continue
            p = int(cs[o]); C[(st, p, L)] += n; CD[(st, L, coords, p)] = 1
    A = {r["domain"]: r for r in rd(f"{WORK}/stats/{sid}/aln_metrics.tsv")}
    assert ntot == int(A["pedv_all"]["reads"]), (sid, ntot)
    na_report.append({**mcols(sid), "PEDV_reads_total": ntot, "PEDV_reads_with_psite": ntot - sum(nNA.values()),
                      "PEDV_reads_offset_NA": sum(nNA.values()), "offset_NA_lengths": ";".join(f"{L}:{n}" for L, n in sorted(nNA.items())) or "-"})

    cat = lambda L: SETS[sid][L]["category"] if L in SETS[sid] else ("long_protected_fragment" if L >= LONG_MIN else "not_included")
    for (st, p, L), n in sorted(C.items(), key=lambda x: (x[0][0], x[0][1], x[0][2])):
        gf = (p - 1) % 3 if st == "+" else "NA(−链)"
        fo.write(f"{sid}\t{m['time_h']}\t{m['drug']}\t{st}\t{p}\t{L}\t{off[min(L,150)]}\t{cat(L)}\t{n}\t{gf}\t{OFS.get(p,'-') if st=='+' else '-'}\n")
        long_rows_written += 1

    Wd = collections.defaultdict(lambda: collections.Counter())
    for (st, p, L), n in C.items():
        for k in SETNAMES:
            if in_set(sid, L, k): Wd[(st, p)][k] += n
    for (st, L, key, p), one in CD.items():
        if in_set(sid, L, "strict_core"): Wd[(st, p)]["strict_core_dedup(敏感性)"] += 1
    wide[sid] = Wd
    for k in SETNAMES:
        a = sum(n for (st, p, L), n in C.items() if in_set(sid, L, k)); b = sum(v[k] for v in Wd.values())
        checks.append({**mcols(sid), "set": k, "long_table_sum": a, "wide_table_sum": b, "match": a == b})

    F5 = collections.defaultdict(lambda: np.zeros(3, np.int64))
    for r in rd(f"{WORK}/stats/{sid}/frame5.tsv"): F5[(r["domain"], int(r["L"]))] += [int(r["f0"]), int(r["f1"]), int(r["f2"])]
    doms = ["PEDV_all"] + [f"PEDV_{n}" for n, a, b in ORFS]
    for dom in doms + ["host"]:
        Ls = sorted({L for (d, L) in F5 if d == dom})
        for L in Ls:
            o = off.get(L)
            if o is None: continue
            p = pf(F5[(dom, L)], o); t = p.sum()
            row = {**mcols(sid), "domain": dom.replace("PEDV_all", "PEDV_all_7ORF"), "read_length": L, "offset": o,
                   "length_category": cat(L), "n_reads(5′端在CDS主体)": int(t),
                   "psite_frame0": f3(p[0] / t), "psite_frame1": f3(p[1] / t), "psite_frame2": f3(p[2] / t), "phase0_pct": f"{100*p[0]/t:.2f}"}
            (per_h if dom == "host" else per_v).append(row)
        for k in ("strict_core", "expanded", "fixed_29_30", "all_lengths"):
            p = np.zeros(3, np.int64)
            for L in Ls:
                if in_set(sid, L, k) and off.get(L) is not None: p += pf(F5[(dom, L)], off[L])
            t = p.sum()
            if t == 0: continue
            row = {**mcols(sid), "domain": dom.replace("PEDV_all", "PEDV_all_7ORF"), "read_length": f"pooled_{k}", "offset": "per-length",
                   "length_category": k, "n_reads(5′端在CDS主体)": int(t),
                   "psite_frame0": f3(p[0] / t), "psite_frame1": f3(p[1] / t), "psite_frame2": f3(p[2] / t), "phase0_pct": f"{100*p[0]/t:.2f}"}
            (per_h if dom == "host" else per_v).append(row)
fo.close()


WCOLS = SETNAMES + ["strict_core_dedup(敏感性)"]
with gzip.open(f"{OUT}/psite_counts_PEDV_GIb_wide.tsv.gz", "wt", encoding="utf-8", newline="") as f:
    f.write(provenance("psite_counts_PEDV_GIb_wide：sample × strand × position，各读长集合合并（每个读长用自己的 offset）", prov_extra +
            ["all_lengths = 该样本所有有 offset 的读长（含 LPF）；strict_core_dedup = 去重敏感性（键=链,读长,5′端坐标或剪接坐标串），主结果不去重"]))
    f.write("sample_id\ttime_h\tdrug\tstrand\tposition\t" + "\t".join(WCOLS) + "\tgenome_frame\torf_frame\n")
    for sid in IDS:
        m = META[sid]
        for (st, p) in sorted(wide[sid]):
            v = wide[sid][(st, p)]
            f.write(f"{sid}\t{m['time_h']}\t{m['drug']}\t{st}\t{p}\t" + "\t".join(str(v[k]) for k in WCOLS) +
                    f"\t{(p-1)%3 if st=='+' else 'NA(−链)'}\t{OFS.get(p,'-') if st=='+' else '-'}\n")
H = lambda t: provenance(t, prov_extra)
write(f"{OUT}/psite_table_crosscheck.tsv", checks, H("长表 vs 宽表互校（各集合逐样本求和）"))
write(f"{OUT}/psite_offset_NA_report.tsv", na_report, H("offset=NA 的 PEDV reads 计数：这些读长在宿主起始密码子上无覆盖，P-site 无法定；"
      "按本任务红线不丢弃，逐条 5′端见 pedv_5p_offsetNA_reads_GIb.tsv.gz"))
write(f"{OUT}/pedv_5p_offsetNA_reads_GIb.tsv.gz", NA_ROWS, H("offset=NA 读长的 PEDV reads（按 5′端坐标保留，不进 P-site 表，不丢弃）"))
write(f"{OUT}/periodicity_virus.tsv", per_v, H("periodicity_virus：PEDV(+) 5′端位于各 canonical ORF 唯一主体区（两端各去 45 nt、去重叠区）的 reads，"
      "按 offset 换算的 P-site 读框（相对该 ORF 自身读框）；只呈现，不作翻译与否结论"))
write(f"{OUT}/periodicity_host.tsv", per_h, H("periodicity_host：宿主 canonical protein_coding CDS 主体（两端各去 45 nt）P-site 读框；phase0_pct=frame0 百分比"))
print("long rows", long_rows_written, "checks all match:", all(c["match"] for c in checks))
for r in na_report: print(r)
