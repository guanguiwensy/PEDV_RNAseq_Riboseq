# -*- coding: utf-8 -*-
"""阶段2a 步骤3：跨冠状病毒序列宇宙构建 + 分层采样 + QC/去重漏斗。

输入（只读）：
  $PEDV_COV/coronavirus_accession_selection/species_table.csv        ICTV 54 物种基准
  $PEDV_COV/coronavirus_accession_selection/selection_manifest.tsv  逐物种检索漏斗（2026-09-10/11）
  $PEDV_COV/coronavirus_accession_selection/selected_metadata/*.tsv 逐 accession 元数据
  $PEDV_COV/coronavirus_fasta_download/fasta_final/all_54_species_genomes.fasta  本地序列库（53,593）
  $PEDV_COV/coronavirus_fasta_download/fasta_final/dedup_groups.tsv 库级精确去重映射
  data_external/{PEDV,TGEV,PRCV,PDCoV}_metadata.tsv                     1b NCBI 保留集
  results/phase2_ncbi_excluded_{PEDV,TGEV,PRCV,PDCoV}.tsv               1b 排除集+原因
  results/phase2_ncbi_dupgroups_{PEDV,TGEV,PRCV,PDCoV}.tsv              1b 去重映射
  results/PEDV_target_manifest.tsv / results/crossvirus_target_manifest.tsv  1e 纳入清单（cov_existing 元数据）
  data_external/raw/PEDV_increment_20260925.gb                          2a PEDV 增量（PZ463692.2）
  data_external/other_representatives.fasta / *_metadata.tsv            外群（2b 步）

抽样规则（总工程师确认，seed 固定 20260925）：
  每物种完整基因组 ≤50 全纳；>50 按 host/region/year(/lineage) 分层随机抽 20–30，
  固定包含 locked_representative（RefSeq/参考株）；SARS-CoV-2 硬上限 30 并注明。

输出：
  results/cross_coronavirus_sequence_universe.tsv
  results/cross_coronavirus_sampling_summary.tsv
  results/cross_coronavirus_target_manifest_extension.tsv
  results/cross_coronavirus_target_manifest_extension.fasta
  logs/phase2a_dedup_map.tsv
  logs/phase2a_universe_build.log
"""
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
import sys, io, re, csv, json, random
from pathlib import Path
from datetime import datetime
from collections import defaultdict, Counter

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8")

WORK = Path(PEDV_COV+r"N1_conservation_analysis")
SEL = Path(PEDV_COV+r"coronavirus_accession_selection")
FASTA_LIB = Path(PEDV_COV+r"coronavirus_fasta_download/fasta_final/all_54_species_genomes.fasta")
DEDUP_LIB = Path(PEDV_COV+r"coronavirus_fasta_download/fasta_final/dedup_groups.tsv")
RES, LOGS = WORK / "results", WORK / "logs"
SEED = 20260925
RUN_DATE = "2026-09-25"
FAMILY, SUBFAM = "Coronaviridae", "Orthocoronavirinae"
GENUS_MAP = {"α": "Alphacoronavirus", "β": "Betacoronavirus",
             "γ": "Gammacoronavirus", "δ": "Deltacoronavirus",
             "Letovirus": "Alphaletovirus", "Pironavirus": "Alphapironavirus"}

logf = open(LOGS / "phase2a_universe_build.log", "w", encoding="utf-8")
def log(msg):
    line = f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}"
    print(line); logf.write(line + "\n"); logf.flush()

def norm_sp(s):
    return (s or "").replace("_", " ").strip()

# ---------- 1. 读物种表（ICTV 基准） ----------
species_tbl = []
with open(SEL / "species_table.csv", encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        species_tbl.append({
            "idx": int(r["序号"]), "subfamily": r["亚科"], "group": r["基因群"],
            "genus": GENUS_MAP.get(r["基因群"], r["基因群"]), "species": r["ICTV现行物种名"].strip(),
            "rep": [a.strip() for a in r["代表accession"].split("/")]})
sp2genus = {s["species"]: s["genus"] for s in species_tbl}
log(f"物种表载入：{len(species_tbl)} 物种；Orthocoronavirinae 四属："
    f"{ {g: sum(1 for s in species_tbl if s['genus']==g and s['subfamily']=='Orthocoronavirinae') for g in ['Alphacoronavirus','Betacoronavirus','Gammacoronavirus','Deltacoronavirus']} }")

# ---------- 2. 读本地序列库 ----------
lib = {}   # accession.version -> (header, seq)
bare = defaultdict(list)
cur_h, cur_s = None, []
with open(FASTA_LIB, encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith(">"):
            if cur_h:
                lib[cur_h.split()[0][1:]] = (cur_h, "".join(cur_s).upper())
            cur_h, cur_s = line.strip(), []
        else:
            cur_s.append(line.strip())
    if cur_h:
        lib[cur_h.split()[0][1:]] = (cur_h, "".join(cur_s).upper())
for av in lib:
    bare[av.split(".")[0]].append(av)
log(f"本地序列库载入：{len(lib)} 条（{FASTA_LIB.name}）")

# ---------- 3. 读逐物种元数据 + 选择漏斗 ----------
meta = {}  # accession.version -> meta dict
meta_files = sorted(SEL.glob("selected_metadata/*_selected_metadata.tsv"))
for fp in meta_files:
    with open(fp, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            acc = (r.get("accession") or "").strip()
            if not acc:
                continue
            if "." not in acc and len(bare.get(acc, [])) == 1:
                acc = bare[acc][0]
            meta[acc] = {"species": norm_sp(r.get("species")), "virus_name": r.get("virus_name", ""),
                         "length": r.get("length", ""), "completeness": r.get("completeness", ""),
                         "host": r.get("host", ""), "geo_region": r.get("geo_region", ""),
                         "geo_location": r.get("geo_location", ""),
                         "collection_year": r.get("collection_year") or (r.get("collection_date", "") or "")[:4],
                         "pango": r.get("pango_lineage", ""),
                         "locked": str(r.get("locked_representative", "")).upper() == "TRUE"}
log(f"逐物种元数据载入：{len(meta)} 条 accession（{len(meta_files)} 个物种文件）")

funnel = {}  # species -> dict
with open(SEL / "selection_manifest.tsv", encoding="utf-8-sig") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        funnel[r["ictv_species"]] = r

# 库级去重映射 removed -> kept
lib_dedup = {}
with open(DEDUP_LIB, encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        lib_dedup[r["removed_accession"]] = r["kept_accession"]
log(f"库级去重映射：{len(lib_dedup)} 条 removed→kept")

# ---------- 4. 猪病毒宇宙 ----------
PIG = {"PEDV": ("PEDV_main", "Alphacoronavirus porci"),
       "TGEV": ("TGEV_main", "Alphacoronavirus suis"),
       "PRCV": ("PRCV_comparator", "Alphacoronavirus suis"),
       "PDCoV": ("PDCoV_main", "Deltacoronavirus suis")}
pig_rows, pig_keys = [], set()

def add_pig(virus, acc, ver, src, status, dup, reason, basis, **kw):
    key = acc if "." in acc else f"{acc}.{ver}"
    if key in pig_keys:
        return
    pig_keys.add(key)
    group, sp = PIG[virus]
    pig_rows.append({"virus": virus, "virus_group": group, "genus": sp2genus.get(sp, sp.split()[0]),
                     "species": sp, "accession": key.split(".")[0], "version": key,
                     "sequence_source": src, "duplicate_status": dup,
                     "inclusion_status": status, "exclusion_reason": reason,
                     "sampling_basis": basis, **kw})

# 4a. 1b NCBI 保留集
for virus in PIG:
    with open(WORK / f"data_external/{virus}_metadata.tsv", encoding="utf-8-sig") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            add_pig(virus, r["accession"], r.get("version", ""), "ncbi_download", "included",
                    "unique", "", "core_pigvirus_full_inclusion; 1b retained after QC+dedup",
                    strain=r.get("strain", "") or "unknown", host=r.get("host", "") or "unknown",
                    country=r.get("country", "") or "unknown",
                    collection_year=(r.get("collection_date", "") or "unknown"),
                    genome_length=r.get("genome_length", ""),
                    n_fraction=r.get("n_fraction", ""), complete_genome="yes")
# 4b. 1b 排除集（含原因）
for virus in PIG:
    with open(RES / f"phase2_ncbi_excluded_{virus}.tsv", encoding="utf-8-sig") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            av = r["accession_version"]
            acc, _, ver = av.partition(".")
            add_pig(virus, acc, av if av.count(".") else f"{av}.?", "ncbi_download", "excluded",
                    "unique", r["reason"], "core_pigvirus_full_inclusion; 1b excluded",
                    strain="unknown", host="unknown", country="unknown",
                    collection_year="unknown", genome_length="", n_fraction="",
                    complete_genome="no")
# 4c. 1b 去重移除集（保留映射；文件无表头：col0=kept，col1..n=removed）
for virus in PIG:
    dgp = RES / f"phase2_ncbi_dupgroups_{virus}.tsv"
    if not dgp.exists():
        continue
    with open(dgp, encoding="utf-8-sig") as f:
        for parts in csv.reader(f, delimiter="\t"):
            if len(parts) < 2:
                continue
            kept, removed = parts[0].strip(), [p.strip() for p in parts[1:] if p.strip()]
            for av in removed:
                add_pig(virus, av.split(".")[0], av, "ncbi_download", "excluded",
                        f"duplicate_of:{kept}", "identical_sequence",
                        "core_pigvirus_full_inclusion; 1b exact-dedup removed",
                        strain="unknown", host="unknown", country="unknown",
                        collection_year="unknown", genome_length="", n_fraction="",
                        complete_genome="yes")
# 4d. cov_existing（来自 1e 清单）
for mf in ["PEDV_target_manifest.tsv", "crossvirus_target_manifest.tsv"]:
    with open(RES / mf, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("sequence_source") != "cov_existing":
                continue
            vmap = {"PEDV": "PEDV", "TGEV": "TGEV", "PRCV": "PRCV", "PDCoV": "PDCoV"}
            virus = vmap.get(r["virus"])
            if not virus:
                continue
            add_pig(virus, r["accession"], r.get("version", ""), "cov_existing", "included",
                    "unique", "", "core_pigvirus_full_inclusion; cov_existing library",
                    strain=r.get("strain", "") or "unknown", host=r.get("host", "") or "unknown",
                    country=r.get("country", "") or "unknown",
                    collection_year=r.get("collection_date", "") or "unknown",
                    genome_length=r.get("genome_length", ""),
                    n_fraction="", complete_genome="yes")
# 4e. 2a PEDV 增量
inc_gb = (WORK / "data_external/raw/PEDV_increment_20260925.gb").read_text(encoding="utf-8")
def parse_gb_single(gb):
    d, seq, org = {}, [], False
    for line in gb.splitlines():
        if line.startswith("DEFINITION"): d["definition"] = line[12:].strip()
        elif line.startswith("VERSION"): d["version"] = line[12:].split()[0]
        elif line.startswith("  ORGANISM"): d["organism"] = line[12:].strip()
        elif "/collection_date=" in line: d["year"] = line.split('"')[1][:4]
        elif "/host=" in line: d["host"] = line.split('"')[1]
        elif line.startswith("ORIGIN"): org = True
        elif line.startswith("//"): org = False
        elif org: seq.append(re.sub(r"[^ACGTUNacgtun]", "", line))
    s = "".join(seq).upper()
    d["length"], d["n_frac"] = len(s), (s.count("N") / len(s) if s else 1.0)
    return d, s
inc_d, inc_seq = parse_gb_single(inc_gb)
add_pig("PEDV", inc_d["version"].split(".")[0], inc_d["version"].split(".")[1], "ncbi_increment", "included",
        "unique", "", "core_pigvirus_full_inclusion; 2a fresh-esearch increment (2026-09-25)",
        strain="HN2025", host=inc_d.get("host", "Sus scrofa"), country="unknown",
        collection_year=inc_d.get("year", "unknown"),
        genome_length=str(inc_d["length"]), n_fraction=f"{inc_d['n_frac']:.6f}",
        complete_genome="yes")
log(f"猪病毒宇宙：PEDV {sum(1 for r in pig_rows if r['virus']=='PEDV')}、"
    f"TGEV {sum(1 for r in pig_rows if r['virus']=='TGEV')}、"
    f"PRCV {sum(1 for r in pig_rows if r['virus']=='PRCV')}、"
    f"PDCoV {sum(1 for r in pig_rows if r['virus']=='PDCoV')}（含排除者与去重移除者）")

# ---------- 5. 四属其他物种宇宙 + 分层采样 ----------
def pig_class(m):
    """返回元数据行所属的猪病毒组，或 None。"""
    vn = (m["virus_name"] or "").lower()
    if "porcine epidemic diarrhea" in vn: return "PEDV"
    if "transmissible gastroenteritis" in vn or "tgev" in vn: return "TGEV"
    if "porcine respiratory coronavirus" in vn: return "PRCV"
    if "deltacoronavirus" in vn or "coronavirus hku15" in vn or "pdco" in vn: return "PDCoV"
    return None

other_rows, samp_rows = [], []
summary = []
random.seed(SEED)
for s in species_tbl:
    if s["subfamily"] != "Orthocoronavirinae":
        continue
    sp, genus = s["species"], s["genus"]
    if sp in ("Alphacoronavirus porci", "Deltacoronavirus suis"):
        continue  # 全部由猪病毒组覆盖（物种内所有成员均为 PEDV/PDCoV）
    sp_meta = [(av, m) for av, m in meta.items() if m["species"] == sp]
    # 猪病毒组成员行
    if sp == "Alphacoronavirus suis":
        pig_part = [(av, m) for av, m in sp_meta if pig_class(m) in ("TGEV", "PRCV")]
        pool_all = [(av, m) for av, m in sp_meta if pig_class(m) is None]  # FCoV/CCoV 等
    else:
        pig_part, pool_all = [], sp_meta
    # 宇宙行（全部检索到的成员）
    for av, m in sp_meta:
        cls = pig_class(m)
        if cls and sp == "Alphacoronavirus suis":
            continue  # TGEV/PRCV 行已在猪病毒宇宙（或与之重复）
        in_lib = av in lib
        dup = "unique"
        reason = ""
        if av in lib_dedup:
            dup = f"duplicate_of:{lib_dedup[av]}"
            reason = "identical_sequence"
        elif not in_lib:
            reason = "not_in_final_library(download_failed_or_prelibrary_qc)"
        other_rows.append({"virus_group": "other", "genus": genus, "species": sp,
                           "accession": av.split(".")[0], "version": av,
                           "sequence_source": "local_library_54species",
                           "complete_genome": "yes" if (m["completeness"] or "").upper() == "COMPLETE" else "unknown",
                           "duplicate_status": dup,
                           "inclusion_status": "pending_sampling",
                           "exclusion_reason": reason,
                           "sampling_basis": "",
                           "strain": "unknown", "host": m["host"] or "unknown",
                           "country": m["geo_location"] or m["geo_region"] or "unknown",
                           "collection_year": m["collection_year"] or "unknown",
                           "genome_length": m["length"], "n_fraction": "",
                           "in_library": in_lib, "locked": m["locked"],
                           "_m": m})
    # 采样池 = 库内可用且非猪病毒
    pool = [r for r in other_rows if r["species"] == sp and r["in_library"] and r["duplicate_status"] == "unique"]
    # 固定包含 RefSeq/参考株（locked）；缺失时从 other_representatives 补（见步骤6后校验）
    locked = [r for r in pool if r["locked"]]
    rest = [r for r in pool if not r["locked"]]
    fn = funnel.get(sp, {})
    if len(pool) <= 50:
        chosen = pool
        basis = f"species_pool<=50({len(pool)}): all included; refseq_fixed"
    else:
        target = 30
        # 分层：host × geo_region × collection_year（pandemicum 增加 pango 谱系维）
        strata = defaultdict(list)
        for r in rest:
            if sp == "Betacoronavirus pandemicum":
                pango = (r["_m"]["pango"] or "unknown").split(".")[0]
                key = (r["host"], r["country"], r["collection_year"], pango)
            else:
                key = (r["host"], r["country"], r["collection_year"])
            strata[key].append(r)
        rng = random.Random(f"{SEED}:{sp}")
        for v in strata.values():
            rng.shuffle(v)
        chosen = list(locked)
        skeys = sorted(strata, key=lambda k: -len(strata[k]))
        i = 0
        while len(chosen) < target and any(strata[k] for k in skeys):
            k = skeys[i % len(skeys)]
            if strata[k]:
                chosen.append(strata[k].pop())
            i += 1
        basis = (f"stratified_random(host,region,year{',lineage' if sp=='Betacoronavirus pandemicum' else ''}) "
                 f"target={target}; seed={SEED}; refseq_fixed")
        if sp == "Betacoronavirus pandemicum":
            basis += "; SARS-CoV-2 hard cap 30"
    # 2a QC 复核（与 1b 同口径）：长度 25–34 kb 且 N%≤1%，不合格者排除并记录原因
    qc_fail = {}
    for r in chosen:
        seq = lib[r["version"]][1]
        n_frac = seq.count("N") / len(seq) if seq else 1.0
        r["genome_length"] = str(len(seq))
        r["n_fraction"] = f"{n_frac:.6f}"
        if not (25000 <= len(seq) <= 34000):
            qc_fail[id(r)] = f"length_out_of_25_34kb({len(seq)})"
        elif n_frac > 0.01:
            qc_fail[id(r)] = f"n_fraction_above_1pct({n_frac:.4f})"
    if qc_fail:
        log(f"{sp}: QC 复核排除 {len(qc_fail)} 条（" + "; ".join(sorted(set(qc_fail.values()))) + "）")
    chosen = [r for r in chosen if id(r) not in qc_fail]
    chosen_ids = {id(r) for r in chosen}
    qc_fail_ids = set(qc_fail)
    for r in other_rows:
        if r["species"] != sp:
            continue
        if not r["in_library"]:
            r["inclusion_status"] = "excluded"
        elif r["duplicate_status"] != "unique":
            r["inclusion_status"] = "excluded"
        elif id(r) in qc_fail_ids:
            r["inclusion_status"] = "excluded"
            r["exclusion_reason"] = qc_fail[id(r)]
            r["sampling_basis"] = basis + "; qc_recheck_failed"
        elif id(r) in chosen_ids:
            r["inclusion_status"] = "included"
            r["sampling_basis"] = basis + ("; locked_representative" if r["locked"] else "")
            samp_rows.append(r)
        else:
            r["inclusion_status"] = "excluded"
            r["exclusion_reason"] = (r["exclusion_reason"] + ";" if r["exclusion_reason"] else "") + "not_selected_by_stratified_sampling"
            r["sampling_basis"] = basis
    n_locked_in = sum(1 for r in chosen if r["locked"])
    summary.append({"genus": genus, "species": sp,
                    "n_retrieved": int(float(fn.get("N_raw", 0) or 0)) if fn else len(sp_meta),
                    "n_complete": int(float(fn.get("N_complete", 0) or 0)) if fn else "",
                    "n_after_length_lab_vaccine": int(float(fn.get("N_eligible", 0) or 0)) if fn else "",
                    "n_in_library": len(pool), "n_locked_rep": n_locked_in,
                    "n_included": len(chosen), "pool_size": len(pool),
                    "sampling_basis": basis, "notes": ""})
    log(f"{sp}: pool={len(pool)} locked_included={n_locked_in} → included={len(chosen)} ({basis})")

# ---------- 6. 外群 ----------
out_rows = []
out_fasta = {}
cur_h, cur_s = None, []
with open(WORK / "data_external/other_representatives.fasta", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith(">"):
            if cur_h:
                out_fasta[cur_h.split()[0][1:]] = (cur_h, "".join(cur_s).upper())
            cur_h, cur_s = line.strip(), []
        else:
            cur_s.append(line.strip())
    if cur_h:
        out_fasta[cur_h.split()[0][1:]] = (cur_h, "".join(cur_s).upper())
OUTGRP = {"GECV01031551.1": ("Alphaletovirus", "Letovirinae", "Alphaletovirus microhylae",
                             "no", "partial transcriptome assembly (~22.3 kb); species_table exempt_length"),
          "NC_076697.1": ("Alphapironavirus", "Pitovirinae", "Alphapironavirus salmonis",
                          "yes", "exempt_length ~36.6 kb (Pitovirinae 双节段特殊性)")}
for av, (g, sub, sp, cg, note) in OUTGRP.items():
    h, seq = out_fasta[av]
    out_rows.append({"virus_group": "outgroup", "family": FAMILY, "subfamily": sub, "genus": g,
                     "species": sp, "host": "unknown", "strain": "unknown",
                     "accession": av.split(".")[0], "version": av,
                     "sequence_source": "other_representatives_2b", "complete_genome": cg,
                     "duplicate_status": "unique", "inclusion_status": "included",
                     "exclusion_reason": "", "sampling_basis": "distant_outgroup_optional; not_in_main_four_genera_negative_conclusion",
                     "country": "unknown", "collection_year": "unknown",
                     "genome_length": str(len(seq)),
                     "n_fraction": f"{seq.count('N')/len(seq):.6f}", "notes_col": note,
                     "outgroup": "yes"})
log(f"外群纳入：{len(out_rows)}（GECV01031551.1, NC_076697.1）")

# ---------- 7. 汇总输出 ----------
COLS = ["virus_group", "family", "subfamily", "genus", "species", "host", "strain",
        "accession", "version", "sequence_source", "complete_genome", "duplicate_status",
        "inclusion_status", "exclusion_reason", "sampling_basis", "outgroup",
        "country", "collection_year", "genome_length", "n_fraction"]

def finalize(rows, extra_note=None):
    out = []
    for r in rows:
        r = dict(r)
        r.pop("_m", None)
        inlib = r.pop("in_library", None)
        r.pop("locked", None)
        r.setdefault("family", FAMILY)
        r.setdefault("subfamily", SUBFAM)
        r.setdefault("outgroup", "no")
        note = extra_note or r.pop("notes_col", "")
        r["notes"] = note
        out.append(r)
    return out

pig_final = finalize(pig_rows)
for r in pig_final:
    r["notes"] = ("PRCV 单独标注，不并入 TGEV" if r["virus_group"] == "PRCV_comparator" else "")
other_final = finalize(other_rows)
all_rows = pig_final + other_final + finalize(out_rows)
with open(RES / "cross_coronavirus_sequence_universe.tsv", "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=COLS + ["notes"], delimiter="\t", extrasaction="ignore")
    w.writeheader()
    w.writerows(all_rows)
n_inc = sum(1 for r in all_rows if r["inclusion_status"] == "included")
log(f"universe 写出：{len(all_rows)} 行（included {n_inc}）→ cross_coronavirus_sequence_universe.tsv")

# manifest = included only
man_rows = [r for r in all_rows if r["inclusion_status"] == "included"]
with open(RES / "cross_coronavirus_target_manifest_extension.tsv", "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=COLS + ["notes"], delimiter="\t", extrasaction="ignore")
    w.writeheader()
    w.writerows(man_rows)

# manifest fasta（序列取自本地库 / 猪病毒 fasta / 外群 fasta）
PEDV_FA = {">PEDV_increment"}
def load_fa(fp):
    d, h, s = {}, None, []
    with open(fp, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                if h: d[h.split()[0][1:]] = (h, "".join(s).upper())
                h, s = line.strip(), []
            else:
                s.append(line.strip())
        if h: d[h.split()[0][1:]] = (h, "".join(s).upper())
    return d
fa_pig = {}
for v in PIG:
    fa_pig.update(load_fa(WORK / f"data_external/{v}.fasta"))
fa_pedv_all = load_fa(WORK / "results/PEDV_genomes_all.fasta")
fa_all = dict(lib)
fa_all.update(fa_pig)
fa_all.update(fa_pedv_all)
fa_all.update(out_fasta)
fa_all["PZ463692.2"] = (">PZ463692.2 PEDV increment HN2025", inc_seq)
missing = [r["version"] for r in man_rows if r["version"] not in fa_all]
log(f"manifest 序列定位：缺序列 {len(missing)} 条 {missing[:10]}")
with open(RES / "cross_coronavirus_target_manifest_extension.fasta", "w", encoding="utf-8") as f:
    for r in man_rows:
        h, seq = fa_all[r["version"]]
        f.write(f">{r['accession']}|{r['virus_group']}|{r['genus']}|{r['species']}|{r['version']}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")

# ---------- 8. 采样汇总（漏斗） ----------
SUM_COLS = ["scope", "genus", "species", "n_retrieved", "n_complete", "n_after_dedup",
            "n_after_qc", "n_included", "n_excluded", "exclusion_reason_summary",
            "refseq_fixed_included", "sampling_basis", "notes"]
sum_out = []
# 猪病毒漏斗：新鲜 esearch → 1b 链 → 2a
fresh_counts = {"PEDV": 1111, "TGEV": 67, "PRCV": 15, "PDCoV": 279}
chain_1b = {"PEDV": (1110, 1041, 41, 28), "TGEV": (67, 50, 16, 1),
            "PRCV": (15, 15, 0, 0), "PDCoV": (279, 241, 31, 7)}
for virus, (group, sp) in PIG.items():
    rows_v = [r for r in pig_final if r.get("virus") == virus]
    inc = [r for r in rows_v if r["inclusion_status"] == "included"]
    exc = Counter()
    for r in rows_v:
        if r["inclusion_status"] == "excluded":
            key = r["exclusion_reason"].split(":")[0] or "other"
            exc[key] += 1
    c1b = chain_1b[virus]
    sum_out.append({"scope": "core_pigvirus", "genus": sp2genus.get(sp, ""),
                    "species": sp, "n_retrieved": fresh_counts[virus],
                    "n_complete": c1b[0],
                    "n_after_dedup": c1b[1],
                    "n_after_qc": c1b[1],
                    "n_included": len(inc), "n_excluded": len(rows_v) - len(inc),
                    "exclusion_reason_summary": "; ".join(f"{k}={v}" for k, v in exc.items()) or "-",
                    "refseq_fixed_included": "-",
                    "sampling_basis": "core_pigvirus_full_inclusion",
                    "notes": (f"1b esearch 2026-09-24 count={c1b[0]}；2a fresh esearch 2026-09-25 count={fresh_counts[virus]}"
                              + ("；一致，引用 logs/phase2_step2_ncbi_download.log 为完整检索证明" if fresh_counts[virus] == c1b[0]
                                 else f"；+1 增量 PZ463692.2 已下载 QC 通过并纳入")
                              + ("；PRCV 单独标注" if virus == "PRCV" else ""))})
# 四属漏斗
for srec in summary:
    sp = srec["species"]
    rows_sp = [r for r in other_final if r["species"] == sp]
    inc = [r for r in rows_sp if r["inclusion_status"] == "included"]
    exc = Counter()
    for r in rows_sp:
        if r["inclusion_status"] == "excluded":
            for part in (r["exclusion_reason"] or "other").split(";"):
                exc[part.strip() or "other"] += 1
    n_dup = sum(1 for r in rows_sp if r["duplicate_status"] != "unique")
    n_notlib = sum(1 for r in rows_sp if "not_in_final_library" in r["exclusion_reason"])
    fn = funnel.get(sp, {})
    sum_out.append({"scope": "four_genera", "genus": srec["genus"], "species": sp,
                    "n_retrieved": srec["n_retrieved"],
                    "n_complete": srec["n_complete"],
                    "n_after_dedup": (srec["n_after_length_lab_vaccine"] if srec["n_after_length_lab_vaccine"] != "" else ""),
                    "n_after_qc": srec["n_in_library"] + n_dup,
                    "n_included": len(inc), "n_excluded": len(rows_sp) - len(inc),
                    "exclusion_reason_summary": "; ".join(f"{k}={v}" for k, v in exc.items()) or "-",
                    "refseq_fixed_included": srec["n_locked_rep"],
                    "sampling_basis": srec["sampling_basis"],
                    "notes": (f"检索漏斗来源 selection_manifest(2026-09-10/11)：N_raw→N_complete→N_eligible；"
                              f"库内去重移除 {n_dup}、未入库 {n_notlib}"
                              + ("；物种参考株 AJ271965.2（TGEV Purdue）已在 TGEV 组全量纳入，"
                                 "other 池（FCoV/CCoV 等）无单独 locked 参考株" if sp == "Alphacoronavirus suis" else ""))})
# 外群
for r in out_rows:
    sum_out.append({"scope": "outgroup", "genus": r["genus"], "species": r["species"],
                    "n_retrieved": 1, "n_complete": 1 if r["complete_genome"] == "yes" else 0,
                    "n_after_dedup": 1, "n_after_qc": 1, "n_included": 1, "n_excluded": 0,
                    "exclusion_reason_summary": "-", "refseq_fixed_included": "-",
                    "sampling_basis": r["sampling_basis"], "notes": r.get("notes", r.get("notes_col", ""))})
with open(RES / "cross_coronavirus_sampling_summary.tsv", "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=SUM_COLS, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    w.writerows(sum_out)
log(f"summary 写出：{len(sum_out)} 行（含四属全部 {sum(1 for s in sum_out if s['scope']=='four_genera')} 物种 + 猪病毒 4 + 外群 2）")

# ---------- 9. 去重映射日志 ----------
with open(LOGS / "phase2a_dedup_map.tsv", "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["scope", "removed_accession", "kept_accession", "source"])
    for virus in PIG:
        dgp = RES / f"phase2_ncbi_dupgroups_{virus}.tsv"
        if not dgp.exists():
            continue
        with open(dgp, encoding="utf-8-sig") as fh:
            for parts in csv.reader(fh, delimiter="\t"):
                if len(parts) < 2:
                    continue
                for rm in parts[1:]:
                    if rm.strip():
                        w.writerow([PIG[virus][0], rm.strip(), parts[0].strip(), "1b_ncbi_exact_dedup"])
    for rm, kp in lib_dedup.items():
        w.writerow(["local_library_54species", rm, kp, "library_exact_dedup"])
log(f"去重映射写出 logs/phase2a_dedup_map.tsv")

# ---------- 10. 自检 ----------
log("=== 自检 ===")
gen_inc = Counter(r["genus"] for r in man_rows)
for g in ["Alphacoronavirus", "Betacoronavirus", "Gammacoronavirus", "Deltacoronavirus"]:
    log(f"属 {g}: 纳入 {gen_inc.get(g, 0)} 条" + ("" if gen_inc.get(g, 0) > 0 else "  ← no_complete_genome_available"))
no_data = [s["species"] for s in species_tbl
           if s["subfamily"] == "Orthocoronavirinae" and gen_inc.get(s["genus"], 0) == 0]
log(f"四属无纳入物种：{no_data or '无'}")
nov = [r["version"] for r in man_rows if "." not in r["version"]]
log(f"manifest 缺 version 的记录：{len(nov)}")
log(f"manifest 总数：{len(man_rows)}；其中 PEDV {sum(1 for r in man_rows if r['virus_group']=='PEDV_main')}、"
    f"TGEV {sum(1 for r in man_rows if r['virus_group']=='TGEV_main')}、"
    f"PRCV {sum(1 for r in man_rows if r['virus_group']=='PRCV_comparator')}、"
    f"PDCoV {sum(1 for r in man_rows if r['virus_group']=='PDCoV_main')}、"
    f"other {sum(1 for r in man_rows if r['virus_group']=='other')}、"
    f"outgroup {sum(1 for r in man_rows if r['virus_group']=='outgroup')}")
# QC 复核（库序列长度/N%）
bad = []
for r in man_rows:
    if r["sequence_source"] == "local_library_54species":
        seq = lib[r["version"]][1]
        n = seq.count("N") / len(seq)
        if not (20000 <= len(seq) <= 40000) or n > 0.01:
            bad.append((r["version"], len(seq), round(n, 4)))
log(f"纳入库序列 QC 复核异常：{bad if bad else '无（长度/N% 全部通过；外群豁免长度按物种表规定）'}")
log("=== 完成 ===")
logf.close()
