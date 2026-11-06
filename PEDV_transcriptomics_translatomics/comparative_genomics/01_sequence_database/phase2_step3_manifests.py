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
import sys, io, json, csv, re
from pathlib import Path
from datetime import datetime

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8")

COV = Path(PEDV_COV+r"")
WORK = COV / "N1_conservation_analysis"
EXIST = WORK / "results/phase2_existing_target_seqs.tsv"
EXT = WORK / "data_external"
FINAL_LIST = COV / "coronavirus_fasta_download/fasta_final/final_accession_list.txt"
REP54 = COV / "orf_analysis/tables/representative_54_genomes.tsv"
PEDV_OUT = WORK / "results/PEDV_target_manifest.tsv"
CROSS_OUT = WORK / "results/crossvirus_target_manifest.tsv"
LOG = WORK / "logs/phase2_step3_manifests.log"

logf = open(LOG, "w", encoding="utf-8")
def log(msg):
    line = f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}"
    print(line); logf.write(line + "\n"); logf.flush()

def base_acc(a):
    return a.split(".")[0]

GROUP = {"PEDV": ("PEDV", "PEDV_main"), "TGEV": ("TGEV", "TGEV_main"),
         "PRCV": ("PRCV", "PRCV_comparator"), "PDCoV": ("PDCoV", "PDCoV_main")}

def main():
    log("=== 阶段2 步骤3：建立目标清单 ===")

    # ---- cov_existing ----
    existing = {}   # base_acc -> row(dict)
    with open(EXIST, encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["qc_status"] != "pass":
                continue
            existing[base_acc(r["accession"])] = r
    log(f"cov_existing QC pass 目标序列：{len(existing)}")


    ncbi = {}       # virus -> {"by_version": {acc_ver: row}, "by_base": {base: [rows]}}
    for virus in ("PEDV", "TGEV", "PRCV", "PDCoV"):
        by_version, by_base = {}, {}
        p = EXT / f"{virus}_metadata.tsv"
        with open(p, encoding="utf-8") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                av = f'{r["accession"]}.{r["version"]}' if r["version"] != "NA" else r["accession"]
                r["_av"] = av
                by_version[av] = r
                by_base.setdefault(base_acc(r["accession"]), []).append(r)
        ncbi[virus] = {"by_version": by_version, "by_base": by_base}
        multi = sum(1 for v in by_base.values() if len(v) > 1)
        log(f"ncbi_download {virus}: {len(by_version)} 条（其中 {multi} 个 base accession 含多版本）")


    manifest_cols = ["virus", "species", "genotype", "genotype_status", "strain",
                     "accession", "version", "genome_length", "sequence_source",
                     "sequence_available", "n1_status", "virus_group",
                     "country", "collection_date", "host", "qc_status",
                     "recombination_flag", "ncbi_current_version", "note"]
    merged = []     # PEDV+TGEV+PRCV+PDCoV
    for virus in ("PEDV", "TGEV", "PRCV", "PDCoV"):
        vname, group = GROUP[virus]
        nv = ncbi[virus]["by_version"]
        nb = ncbi[virus]["by_base"]
        local_rows = [r for k, r in existing.items() if r["virus"] == vname]
        matched_ncbi = set()
        rows_v = []

        for e in local_rows:
            acc = e["accession"]
            n = nv.get(acc)
            ver = acc.split(".", 1)[1] if "." in acc else "NA"
            if n:
                matched_ncbi.add(acc)
                ncbi_ver, note = n["version"], "also_in_ncbi_download"
            else:
                cands = nb.get(base_acc(acc), [])
                if cands:
                    n = sorted(cands, key=lambda x: x["_av"])[-1]
                    ncbi_ver, note = n["version"], "version_mismatch_with_ncbi"
                else:
                    n = None
                    ncbi_ver, note = "NA", "not_in_current_ncbi_query"
            recom = n["recombination_flag"] if n else "unknown"
            gnote = (n["genotype_note"] if n else "") or ""
            rows_v.append((e["species"], e["strain"], acc, ver, e["genome_length"],
                           e["geo_location"], e["collection_date"], e["host"],
                           "cov_existing", ncbi_ver, note, recom, gnote,
                           e.get("n1_status", "pending_stage3")))

        for av, n in nv.items():
            if av in matched_ncbi:
                continue
            if base_acc(n["accession"]) in {base_acc(a) for a in matched_ncbi}:
                continue
            if any(base_acc(n["accession"]) == base_acc(r[2]) for r in rows_v):
                continue
            rows_v.append((n["organism"], n["strain"] or n["isolate"] or "Unknown",
                           n["accession"], n["version"], n["genome_length"],
                           n["country"], n["collection_date"], n["host"],
                           "ncbi_download", n["version"], "", n["recombination_flag"],
                           n["genotype_note"] or "", "pending_stage3"))
        log(f"{vname}: 合并后 {len(rows_v)} 条（本地 {len(local_rows)}，其中与 NCBI 重合 "
            f"{len(matched_ncbi)}；纯 NCBI 新增 {len(rows_v) - len(local_rows)}）")
        for species, strain, acc, ver, glen, country, cdate, host, source,\
                ncbi_ver, note, recom, gnote, n1s in rows_v:
            gnote = gnote or ""
            genotype, gstatus = "unknown", "unknown"
            if gnote.startswith("record_annotation:"):
                m = re.search(r"genotype\s*[:=]?\s*([A-Za-z0-9\-/]+)", gnote, re.I)
                token = m.group(1) if m else ""

                if token and token.upper() not in ("INDEL", "NA", "UNKNOWN", "NONE"):
                    genotype, gstatus = token, "verified"
            if acc == "PP999774.1":
                n1s = "reference_annotated"
            merged.append({
                "virus": vname, "species": species,
                "genotype": genotype, "genotype_status": gstatus, "strain": strain,
                "accession": acc, "version": ver, "genome_length": glen,
                "sequence_source": source, "sequence_available": "yes",
                "n1_status": n1s, "virus_group": group,
                "country": country, "collection_date": cdate, "host": host,
                "qc_status": "pass", "recombination_flag": recom,
                "ncbi_current_version": ncbi_ver, "note": note,
                "_gnote": gnote,
            })


    final_set = set(x.strip() for x in open(FINAL_LIST, encoding="utf-8") if x.strip())
    rep_rows = []
    with open(REP54, encoding="utf-8") as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            acc = r.get("accession") or list(r.values())[0]
            sp = r.get("species") or list(r.values())[1]
            if any(k in sp for k in ("porci", "Alphacoronavirus suis", "Deltacoronavirus suis")):
                continue
            if acc in final_set:
                rep_rows.append((acc, sp))
    log(f"other 代表基因组（rep-54 去除目标物种、且存在于本库）：{len(rep_rows)}")
    for acc, sp in rep_rows:
        merged.append({
            "virus": "other_coronavirus", "species": sp,
            "genotype": "unknown", "genotype_status": "unknown", "strain": "representative",
            "accession": acc, "version": acc.split(".", 1)[1] if "." in acc else "NA",
            "genome_length": "NA", "sequence_source": "cov_existing",
            "sequence_available": "yes", "n1_status": "pending_stage3",
            "virus_group": "other", "country": "NA", "collection_date": "NA",
            "host": "NA", "qc_status": "pass", "recombination_flag": "unknown",
            "ncbi_current_version": "NA", "note": "rep-54 species representative (in local final library)",
            "_gnote": "",
        })


    rep_ext = EXT / "other_representatives_metadata.tsv"
    if rep_ext.exists():
        n_ext = 0
        with open(rep_ext, encoding="utf-8") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                acc = f'{r["accession"]}.{r["version"]}' if r["version"] != "NA" else r["accession"]
                merged.append({
                    "virus": "other_coronavirus",
                    "species": r["species"] or r["organism"],
                    "genotype": "unknown", "genotype_status": "unknown",
                    "strain": r["strain"] or r["isolate"] or "representative",
                    "accession": acc, "version": r["version"],
                    "genome_length": r["genome_length"],
                    "sequence_source": "ncbi_download",
                    "sequence_available": "yes",
                    "n1_status": "pending_stage3",
                    "virus_group": "other",
                    "country": "NA", "collection_date": "NA",
                    "host": r["organism"],
                    "qc_status": "pass" if r["qc_status"] == "pass" else "fail",
                    "recombination_flag": "unknown",
                    "ncbi_current_version": r["version"],
                    "note": r.get("note", ""),
                    "_gnote": "",
                })
                n_ext += 1
        log(f"step2b 补充的其他物种代表：{n_ext} 条")


    n_other_local = 0
    for k, e in existing.items():
        if e["virus_group"] != "other":
            continue
        acc = e["accession"]
        merged.append({
            "virus": e["virus"], "species": e["species"],
            "genotype": "unknown", "genotype_status": "unknown",
            "strain": e["strain"], "accession": acc,
            "version": acc.split(".", 1)[1] if "." in acc else "NA",
            "genome_length": e["genome_length"], "sequence_source": "cov_existing",
            "sequence_available": "yes", "n1_status": "pending_stage3",
            "virus_group": "other", "country": e["geo_location"],
            "collection_date": e["collection_date"], "host": e["host"],
            "qc_status": "pass", "recombination_flag": "unknown",
            "ncbi_current_version": "NA", "note": "alphacoronavirus-1 member (TGEV 近缘对照)",
            "_gnote": "",
        })
        n_other_local += 1
    log(f"本地位 Alphacoronavirus-1 其他成员（近缘对照）：{n_other_local} 条")


    def write_tsv(path, rows):
        with open(path, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=manifest_cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for r in rows:
                w.writerow(r)
    pedv_rows = [r for r in merged if r["virus"] == "PEDV"]
    write_tsv(PEDV_OUT, pedv_rows)
    log(f"PEDV_target_manifest: {len(pedv_rows)} 行 → {PEDV_OUT}")

    cross_rows = [r for r in merged if r["virus"] != "PEDV"]
    write_tsv(CROSS_OUT, cross_rows)
    log(f"crossvirus_target_manifest: {len(cross_rows)} 行 → {CROSS_OUT}")


    stats = {}
    for r in merged:
        key = (r["virus"], r["sequence_source"])
        stats[key] = stats.get(key, 0) + 1

    other_local = sum(1 for k, r in existing.items() if r["virus_group"] == "other")
    stats[("other_alphacoronavirus1_members(cov_existing)", "cov_existing")] = other_local
    with open(WORK / "logs/phase2_step3_run_stats.json", "w", encoding="utf-8") as f:
        json.dump({f"{k[0]}|{k[1]}": v for k, v in sorted(stats.items())}, f,
                  ensure_ascii=False, indent=2)
    for k, v in sorted(stats.items()):
        log(f"  {k[0]} [{k[1]}] = {v}")
    log("=== 步骤3 完成 ===")
    logf.close()

if __name__ == "__main__":
    main()
