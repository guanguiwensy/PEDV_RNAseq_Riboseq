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
import sys, io, json, time, urllib.parse, urllib.request
from pathlib import Path
from datetime import datetime

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8")

WORK = Path(PEDV_COV+r"N1_conservation_analysis")
LOG = WORK / "logs/phase2a_fresh_esearch.log"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
RUN_DATE = "2026-09-25"

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


COUNTS_1B = {"PEDV": 1110, "TGEV": 67, "PRCV": 15, "PDCoV": 279}

logf = open(LOG, "w", encoding="utf-8")
def log(msg):
    line = f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}"
    print(line); logf.write(line + "\n"); logf.flush()

def get(url, params, retries=4):
    full = f"{url}?{urllib.parse.urlencode(params)}"
    for i in range(retries):
        try:
            req = urllib.request.Request(full, headers={"User-Agent": "pedv-n1-extension/2a"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read().decode("utf-8", errors="replace")
        except Exception as e:
            log(f"  请求失败({i+1}/{retries}) {url}: {e}")
            time.sleep(3 * (i + 1))
    raise RuntimeError(f"多次重试仍失败: {full}")

def esearch(query):
    txt = get(f"{EUTILS}/esearch.fcgi",
              {"db": "nucleotide", "term": query, "retmax": 100000, "retmode": "json"})
    j = json.loads(txt)
    return j["esearchresult"]["idlist"], int(j["esearchresult"]["count"])

def accession_of(uid):
    txt = get(f"{EUTILS}/esummary.fcgi",
              {"db": "nucleotide", "id": uid, "retmode": "json"})
    j = json.loads(txt)
    d = j["result"][uid]
    return d.get("accessionversion") or f"{d.get('accession')}.{d.get('version')}"

log(f"=== 阶段2a：新鲜 esearch 计数核对（只计数不下载） ===")
log(f"数据库：NCBI Nucleotide（E-utilities esearch，v1，经 esummary 取 accession.version）；"
    f"检索日期 {RUN_DATE}；E-utilities 不暴露数据库版本号，按访问日期留痕")
result = {}
for virus, query in QUERIES.items():
    log(f"--- {virus} ---")
    log(f"检索式：{query}")
    ids, count = esearch(query)
    log(f"esearch 命中 count={count}（1b 记录 {COUNTS_1B[virus]}）")
    time.sleep(0.5)

    accs = []
    for i in range(0, len(ids), 100):
        batch_txt = get(f"{EUTILS}/esummary.fcgi",
                        {"db": "nucleotide", "id": ",".join(ids[i:i+100]), "retmode": "json"})
        j = json.loads(batch_txt)
        for uid in j["result"]["uids"]:
            d = j["result"][uid]
            accs.append(d.get("accessionversion") or f"{d.get('accession')}.{d.get('version')}")
        time.sleep(0.5)
    (WORK / f"logs/phase2a_esearch_ids_{virus}.txt").write_text(
        "\n".join(sorted(accs)) + "\n", encoding="utf-8")
    log(f"accession.version 取回 {len(accs)} 条 → logs/phase2a_esearch_ids_{virus}.txt")
    result[virus] = {"count": count, "n_ids": len(ids), "n_accs": len(accs)}

log("=== 核对结论 ===")
for virus in QUERIES:
    c1b, c2a = COUNTS_1B[virus], result[virus]["count"]
    if c2a == c1b:
        log(f"{virus}: 新命中数 {c2a} 与 1b 一致 → 引用 logs/phase2_step2_ncbi_download.log 作为完整检索证明")
    else:
        log(f"{virus}: 新命中数 {c2a} ≠ 1b 的 {c1b} → 需比对 accession 清单，仅下载增量")
log(json.dumps(result, ensure_ascii=False))
log("=== 完成 ===")
logf.close()
