# -*- coding: utf-8 -*-
"""阶段2 步骤2：NCBI 定点补充下载 PEDV / TGEV / PRCV / PDCoV 完整基因组。

检索式沿用项目 species_table.csv 的病毒级口径；仅下载目标病毒，不重下整个冠状病毒库。
原始 GenBank 存 data_external/raw/，解析后输出 FASTA + metadata TSV，并做：
  筛选（complete、病毒名明确、排除 partial/合成/专利/载体、长度 25–34kb、N<=1%）
  → 完全相同序列精确去重 → QC 标记 → recombination_flag（仅基于记录注释，不删除）。

输出：
  data_external/raw/{virus}.gb
  data_external/{virus}.fasta / {virus}_metadata.tsv
  logs/phase2_step2_ncbi_download.log
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
import sys, io, os, re, time, json, urllib.parse, urllib.request
from pathlib import Path
from datetime import datetime

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8")

WORK = Path(PEDV_COV+r"N1_conservation_analysis")
EXT = WORK / "data_external"
RAW = EXT / "raw"
LOG = WORK / "logs/phase2_step2_ncbi_download.log"
for d in (EXT, RAW, WORK / "logs"):
    d.mkdir(parents=True, exist_ok=True)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ACCESS_DATE = "2026-09-24"
BATCH = 200          # efetch genbank 每批 id 数
SLEEP = 0.45         # 无 API key，限速 3 req/s 以内

# 病毒级检索式（与 species_table.csv 同一口径，病毒名替换为单一病毒）
QUERIES = {
    "PEDV": ('("Porcine epidemic diarrhea virus"[Organism] OR '
             '"Porcine epidemic diarrhea virus"[Title]) AND biomol_genomic[PROP] AND '
             '25000:34000[SLEN] AND ("complete genome"[Title] OR "complete sequence"[Title]) '
             'NOT partial[Title]'),
    "TGEV": ('("Transmissible gastroenteritis virus"[Organism] OR '
             '"Transmissible gastroenteritis virus"[Title]) AND biomol_genomic[PROP] AND '
             '25000:34000[SLEN] AND ("complete genome"[Title] OR "complete sequence"[Title]) '
             'NOT partial[Title]'),
    "PRCV": ('("Porcine respiratory coronavirus"[Organism] OR '
             '"Porcine respiratory coronavirus"[Title]) AND biomol_genomic[PROP] AND '
             '25000:34000[SLEN] AND ("complete genome"[Title] OR "complete sequence"[Title]) '
             'NOT partial[Title]'),
    "PDCoV": ('("Porcine deltacoronavirus"[Organism] OR "Porcine deltacoronavirus"[Title] OR '
              '"Porcine coronavirus HKU15"[Organism] OR "Porcine coronavirus HKU15"[Title] OR '
              '"Swine deltacoronavirus"[Title] OR "Coronavirus HKU15"[Title]) AND '
              'biomol_genomic[PROP] AND 25000:34000[SLEN] AND '
              '("complete genome"[Title] OR "complete sequence"[Title]) NOT partial[Title]'),
}

EXCLUDE_TITLE = re.compile(r"partial|synthetic|patent|vector|cloning|clone|construct|"
                           r"unknown( isolate| host| source)|unverified", re.I)

logf = open(LOG, "w", encoding="utf-8")
def log(msg):
    line = f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}"
    print(line); logf.write(line + "\n"); logf.flush()

def get(url, params, binary=False, retries=4):
    qs = urllib.parse.urlencode(params)
    full = f"{url}?{qs}"
    for i in range(retries):
        try:
            req = urllib.request.Request(full, headers={"User-Agent": "pedv-n1-conservation/1.0"})
            with urllib.request.urlopen(req, timeout=120) as r:
                data = r.read()
            return data if binary else data.decode("utf-8", errors="replace")
        except Exception as e:
            log(f"  请求失败({i+1}/{retries}) {url}: {e}")
            time.sleep(3 * (i + 1))
    raise RuntimeError(f"多次重试仍失败: {full}")

def esearch(query):
    txt = get(f"{EUTILS}/esearch.fcgi", {
        "db": "nucleotide", "term": query, "retmax": 100000,
        "retmode": "json", "tool": "pedv_n1", "email": "NA"})
    j = json.loads(txt)
    return j["esearchresult"]["idlist"], int(j["esearchresult"]["count"])

def efetch_gb(idlist):
    return get(f"{EUTILS}/efetch.fcgi", {
        "db": "nucleotide", "id": ",".join(idlist), "rettype": "gb",
        "retmode": "text", "tool": "pedv_n1", "email": "NA"})

# ---- 极简 GenBank 解析（避免依赖 BioPython 的 LOCUS 容错差异）----
def parse_gb_records(text):
    recs, cur = [], []
    for line in text.splitlines():
        if line.startswith("LOCUS") and cur:
            recs.append("\n".join(cur)); cur = []
        cur.append(line)
    if cur:
        recs.append("\n".join(cur))
    return recs

def parse_record(rec):
    def field(name):
        m = re.search(rf"^{name}\s+(.+)$", rec, re.M)
        return m.group(1).strip() if m else ""
    def clean(s):
        # GenBank 限定符值可跨行（如 isolation_source），折叠所有空白为单空格，保证 TSV 一行一条
        return " ".join((s or "").split())
    acc_v = field("VERSION") or field("ACCESSION").split()[0]
    definition = field("DEFINITION")
    organism = ""
    m = re.search(r"^  ORGANISM\s+(.+)$", rec, re.M)
    if m: organism = m.group(1).strip()
    quals = {}
    for q in re.finditer(r"/(\w+)=\"([^\"]*)\"", rec):
        quals.setdefault(q.group(1), []).append(clean(q.group(2)))
    # 序列
    m = re.search(r"^ORIGIN\s*$", rec, re.M)
    seq = ""
    if m:
        seq = re.sub(r"[^ACGTUNacgtun]", "", rec[m.end():])
    return {
        "accession_version": acc_v,
        "definition": clean(definition),
        "organism": clean(organism),
        "strain": clean((quals.get("strain") or [""])[0]),
        "isolate": clean((quals.get("isolate") or [""])[0]),
        "collection_date": clean((quals.get("collection_date") or [""])[0]),
        "country": clean((quals.get("country") or [""])[0]),
        "host": clean((quals.get("host") or [""])[0]),
        "isolation_source": clean((quals.get("isolation_source") or [""])[0]),
        "note": clean("; ".join(quals.get("note", [])[:3])),
        "title_flags_recombination": bool(re.search(r"recombinant", rec, re.I)
                                          and re.search(r"recombinant (virus|strain|isolate|PEDV|genome)",
                                                        rec, re.I)),
        "sequence": seq.upper().replace("U", "T"),
    }

def main():
    log("=== 阶段2 步骤2：NCBI 定点补充下载 ===")
    log(f"数据库：NCBI Nucleotide（E-utilities esearch/efetch）；访问日期 {ACCESS_DATE}；"
        f"E-utilities 不暴露数据库版本号，按访问日期留痕")
    log(f"筛选条件：complete genome/sequence；病毒名明确；排除 partial/合成/专利/载体/克隆/来源不明；"
        f"长度 25–34 kb；N 比例 ≤1%；完全相同序列精确去重；不删除重组株，仅标记 recombination_flag")
    log(f"download executed locally (access date {ACCESS_DATE})")

    run_stats = {}
    reuse_raw = os.environ.get("PHASE2_REUSE_RAW") == "1"
    for virus, query in QUERIES.items():
        log(f"--- {virus} ---")
        raw_path = RAW / f"{virus}.gb"
        all_recs, failed = [], []
        if reuse_raw and raw_path.exists() and raw_path.stat().st_size > 0:
            log(f"PHASE2_REUSE_RAW=1：跳过 esearch/efetch，直接从 {raw_path} 重新解析")
            all_recs = parse_gb_records(raw_path.read_text(encoding="utf-8"))
            count = len(all_recs)
            log(f"重新解析记录数：{len(all_recs)}")
        else:
            log(f"检索式：{query}")
            ids, count = esearch(query)
            log(f"esearch 命中 count={count}，取回 id={len(ids)}")
            time.sleep(SLEEP)
            with open(raw_path, "w", encoding="utf-8") as out:
                for i in range(0, len(ids), BATCH):
                    batch = ids[i:i + BATCH]
                    try:
                        txt = efetch_gb(batch)
                        out.write(txt if txt.endswith("\n") else txt + "\n")
                        all_recs.extend(parse_gb_records(txt))
                    except Exception as e:
                        failed.append((i, i + len(batch), str(e)))
                    time.sleep(SLEEP)
        log(f"efetch/解析完成：{len(all_recs)} 条记录（raw: {raw_path}）；失败批次：{len(failed)}")
        for f in failed:
            log(f"  失败批次 id[{f[0]}:{f[1]}]: {f[2]}")

        # 解析 + 筛选 + QC + 去重
        kept, excluded = [], []
        for rec in all_recs:
            r = parse_record(rec)
            title = r["definition"]
            reason = None
            if EXCLUDE_TITLE.search(title):
                reason = "标题含排除关键词(partial/合成/专利/载体/来源不明)"
            elif r["organism"] and not re.search(
                    {"PEDV": r"epidemic diarrhea",
                     "TGEV": r"gastroenteritis virus|TGEV",
                     "PRCV": r"respiratory coronavirus",
                     "PDCoV": r"deltacoronavirus|HKU15"}[virus],
                    r["organism"], re.I):
                reason = f"物种名不匹配: {r['organism']}"
            else:
                ln = len(r["sequence"])
                nfrac = (r["sequence"].count("N") / ln) if ln else 1.0
                if not (25000 <= ln <= 34000):
                    reason = f"长度 {ln} 超出 25–34 kb"
                elif nfrac > 0.01:
                    reason = f"N 比例 {nfrac:.3f} > 1%"
            if reason:
                excluded.append((r["accession_version"], reason))
            else:
                kept.append(r)
        # 完全相同序列精确去重（保留 accession 字典序最小者）
        by_seq = {}
        for r in kept:
            by_seq.setdefault(r["sequence"], []).append(r)
        dedup, dup_groups = [], []
        for seq, rs in by_seq.items():
            rs.sort(key=lambda x: x["accession_version"])
            dedup.append(rs[0])
            if len(rs) > 1:
                dup_groups.append([x["accession_version"] for x in rs])
        log(f"{virus}: 解析 {len(all_recs)} → 筛选排除 {len(excluded)} → "
            f"完全相同序列去重 -{len(kept) - len(dedup)} → 保留 {len(dedup)}")
        for acc, reason in excluded[:20]:
            log(f"  排除 {acc}: {reason}")
        if len(excluded) > 20:
            log(f"  …其余 {len(excluded) - 20} 条排除记录见 results/phase2_ncbi_excluded_{virus}.tsv")
        with open(WORK / f"results/phase2_ncbi_excluded_{virus}.tsv", "w", encoding="utf-8") as out:
            out.write("accession_version\treason\n")
            for acc, reason in excluded:
                out.write(f"{acc}\t{reason}\n")
        if dup_groups:
            with open(WORK / f"results/phase2_ncbi_dupgroups_{virus}.tsv", "w", encoding="utf-8") as out:
                for g in dup_groups:
                    out.write("\t".join(g) + "\n")

        # 输出 FASTA + metadata
        fasta_path = EXT / f"{virus}.fasta"
        meta_path = EXT / f"{virus}_metadata.tsv"
        cols = ["virus", "accession", "version", "definition", "organism", "strain", "isolate",
                "collection_date", "country", "host", "isolation_source", "genome_length",
                "n_fraction", "recombination_flag", "genotype_note", "source_url",
                "sequence_source", "download_date"]
        with open(fasta_path, "w", encoding="utf-8") as fa, \
             open(meta_path, "w", encoding="utf-8") as mt:
            mt.write("\t".join(cols) + "\n")
            for r in sorted(dedup, key=lambda x: x["accession_version"]):
                ln = len(r["sequence"])
                nfrac = r["sequence"].count("N") / ln if ln else 1.0
                acc = r["accession_version"].split(".")[0]
                fa.write(f">{r['accession_version']} {virus} {r['strain'] or r['isolate'] or ''}\n")
                for i in range(0, ln, 60):
                    fa.write(r["sequence"][i:i + 60] + "\n")
                # 基因型：仅当记录明确注释时记 verified，否则 unknown（推断留给后续阶段）
                g_note = ""
                m = re.search(r"genotype\s*[:=\s]?\s*([A-Za-z0-9\-/]+)", r["note"] + " " + r["definition"], re.I)
                if m and re.search(r"genotype", r["note"], re.I):
                    g_note = f"record_annotation:{m.group(0)[:80]}"
                mt.write("\t".join([
                    virus, acc, r["accession_version"].split(".", 1)[1] if "." in r["accession_version"] else "NA",
                    r["definition"], r["organism"], r["strain"], r["isolate"],
                    r["collection_date"], r["country"], r["host"], r["isolation_source"],
                    str(ln), f"{nfrac:.6f}",
                    "yes" if r["title_flags_recombination"] else "not_annotated",
                    g_note,
                    f"https://www.ncbi.nlm.nih.gov/nuccore/{acc}",
                    "ncbi_download", ACCESS_DATE,
                ]) + "\n")
        run_stats[virus] = {"esearch_count": count, "fetched": len(all_recs),
                            "excluded": len(excluded), "dedup_removed": len(kept) - len(dedup),
                            "final": len(dedup), "failed_batches": len(failed)}
        log(f"输出：{fasta_path} / {meta_path}")

    with open(WORK / "logs/phase2_step2_run_stats.json", "w", encoding="utf-8") as f:
        json.dump(run_stats, f, ensure_ascii=False, indent=2)
    log("=== 步骤2 完成 ===")
    log(json.dumps(run_stats, ensure_ascii=False))
    logf.close()

if __name__ == "__main__":
    main()
