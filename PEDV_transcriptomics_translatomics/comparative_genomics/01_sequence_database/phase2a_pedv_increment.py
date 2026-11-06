# -*- coding: utf-8 -*-
"""阶段2a 步骤2：PEDV 增量下载（PZ463692.2，2026-09-25 新鲜 esearch 新增）。

背景：1b（2026-09-24）PEDV esearch count=1110；2026-09-25 新鲜核对 count=1111，
新增 accession PZ463692.2（其余 1110 条 accession.version 完全一致，见
logs/phase2a_esearch_ids_PEDV.txt 比对）。按任务规定只下载增量。

下载后做与 1b 相同的 QC：长度 25–34 kb；N 比例 ≤1%；标题排除 partial/合成/专利/载体；
与现有 PEDV 1083 条（results/PEDV_target_manifest.tsv + data_external 库）精确序列去重比对。

输出：
  data_external/raw/PEDV_increment_20260925.gb
  logs/phase2a_pedv_increment.log
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
import sys, io, re, time, json, urllib.parse, urllib.request
from pathlib import Path
from datetime import datetime

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8")

WORK = Path(PEDV_COV+r"N1_conservation_analysis")
EXT = WORK / "data_external"
LOG = WORK / "logs/phase2a_pedv_increment.log"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NEW_ACC = "PZ463692.2"

EXCLUDE_TITLE = re.compile(r"partial|synthetic|patent|vector|cloning|clone|construct|"
                           r"unknown( isolate| host| source)|unverified", re.I)

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

log("=== 阶段2a：PEDV 增量下载（新增 accession 1 条） ===")
log(f"数据库：NCBI Nucleotide efetch；日期 2026-09-25；工具 E-utilities v1；增量对象 {NEW_ACC}")
gb = get(f"{EUTILS}/efetch.fcgi",
         {"db": "nucleotide", "id": NEW_ACC, "rettype": "gb", "retmode": "text"})
raw_path = EXT / "raw/PEDV_increment_20260925.gb"
raw_path.write_text(gb, encoding="utf-8")
log(f"原始 GenBank 已存 {raw_path}（{len(gb)} 字符）")

# 解析单条记录
acc = ver = definition = organism = strain = country = date = host = ""
seq_parts, in_origin = [], False
seq_len_gb = None
for line in gb.splitlines():
    if line.startswith("DEFINITION"):
        definition = line[12:].strip()
    elif line.startswith("ACCESSION"):
        acc = line[12:].split()[0]
    elif line.startswith("VERSION"):
        ver = line[12:].split()[0]
    elif line.startswith("SOURCE") or line.startswith("  ORGANISM"):
        if line.startswith("  ORGANISM"):
            organism = line[12:].strip()
    elif "/strain=" in line:
        strain = line.split('/strain="')[1].split('"')[0]
    elif "/country=" in line:
        country = line.split('/country="')[1].split('"')[0]
    elif "/collection_date=" in line:
        date = line.split('/collection_date="')[1].split('"')[0]
    elif "/host=" in line:
        host = line.split('/host="')[1].split('"')[0]
    elif line.startswith("ORIGIN"):
        in_origin = True
    elif line.startswith("//"):
        in_origin = False
    elif in_origin:
        seq_parts.append(re.sub(r"[^ACGTUNacgtun]", "", line))
seq = "".join(seq_parts).upper()
n_frac = seq.count("N") / len(seq) if seq else 1.0
log(f"解析：{ver} | len={len(seq)} | N%={n_frac:.4f} | definition={definition}")

# QC（与 1b 同口径）
fail = []
if EXCLUDE_TITLE.search(definition):
    fail.append("标题含排除关键词")
if not (25000 <= len(seq) <= 34000):
    fail.append(f"长度 {len(seq)} 超出 25–34 kb")
if n_frac > 0.01:
    fail.append(f"N 比例 {n_frac:.2%} > 1%")
if "Porcine epidemic diarrhea virus" not in organism:
    fail.append(f"物种名不匹配: {organism}")
if fail:
    log(f"QC 未通过：{'；'.join(fail)} → 标记 excluded，不纳入")
else:
    log("QC 通过（长度/N%/标题/物种名）")
    # 精确序列去重比对：data_external/PEDV.fasta + PEDV_genomes_all
    dup_of = None
    for fp, tag in [(EXT / "PEDV.fasta", "ncbi_1041"),
                    (WORK / "results/PEDV_genomes_all.fasta", "all_1083")]:
        cur, name = [], ""
        with open(fp, encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith(">"):
                    if cur and "".join(cur).upper() == seq:
                        dup_of = (tag, name)
                        break
                    cur, name = [], line[1:].split()[0]
                else:
                    cur.append(line.strip())
            else:
                if cur and "".join(cur).upper() == seq:
                    dup_of = (tag, name)
        if dup_of:
            break
    if dup_of:
        log(f"精确去重：与 {dup_of[0]} 中 {dup_of[1]} 完全相同 → 标记 duplicate_removed，保留映射")
    else:
        log("精确去重：无相同序列 → 可作为新纳入序列")

log(json.dumps({"accession": acc, "version": ver, "definition": definition,
                "organism": organism, "strain": strain, "country": country,
                "collection_date": date, "host": host, "length": len(seq),
                "n_fraction": round(n_frac, 6), "qc": fail or "pass",
                "duplicate_of": dup_of[1] if dup_of else None},
               ensure_ascii=False))
log("=== 完成 ===")
logf.close()
