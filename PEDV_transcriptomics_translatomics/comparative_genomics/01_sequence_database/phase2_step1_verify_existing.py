# -*- coding: utf-8 -*-
"""阶段2 步骤1：逐一核实已有 54 物种数据中 PEDV/TGEV/PRCV/PDCoV 的存在情况并留痕。

输入（只读）：
  $PEDV_COV/coronavirus_fasta_download/fasta_final/all_54_species_genomes.fasta
  $PEDV_COV/coronavirus_fasta_download/fasta_final/final_accession_list.txt
  $PEDV_COV/coronavirus_accession_selection/selected_metadata/{3个目标物种}.tsv
  $PEDV_COV/coronavirus_accession_selection/all_54_species_selected_accessions.tsv
  $PEDV_COV/orf_analysis/tables/representative_54_genomes.tsv

输出：
  results/phase2_existing_target_seqs.tsv   （分类+QC 后的本地位目标序列表）
  logs/phase2_step1_verify_existing.log     （逐项核查留痕）
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
import sys, io, os, re
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8")

COV = Path(PEDV_COV+r"")
WORK = COV / "N1_conservation_analysis"
FINAL_FASTA = COV / "coronavirus_fasta_download/fasta_final/all_54_species_genomes.fasta"
FINAL_LIST = COV / "coronavirus_fasta_download/fasta_final/final_accession_list.txt"
SEL = COV / "coronavirus_accession_selection/all_54_species_selected_accessions.tsv"
META_DIR = COV / "coronavirus_accession_selection/selected_metadata"
REP54 = COV / "orf_analysis/tables/representative_54_genomes.tsv"
OUT_TSV = WORK / "results/phase2_existing_target_seqs.tsv"
LOG = WORK / "logs/phase2_step1_verify_existing.log"

SPECIES_FILES = {
    "21": ("Alphacoronavirus porci", META_DIR / "Alphacoronavirus_porci_selected_metadata.tsv"),
    "28": ("Alphacoronavirus suis", META_DIR / "Alphacoronavirus_suis_selected_metadata.tsv"),
    "52": ("Deltacoronavirus suis", META_DIR / "Deltacoronavirus_suis_selected_metadata.tsv"),
}

logf = open(LOG, "w", encoding="utf-8")
def log(msg):
    line = f"[{__import__('datetime').datetime.now():%Y-%m-%d %H:%M:%S}] {msg}"
    print(line)
    logf.write(line + "\n"); logf.flush()

def classify(virus_name: str):
    """返回 (virus, virus_group)。病毒名规则基于 virus_name 字段，不看文件名。"""
    v = (virus_name or "").strip()
    vl = v.lower()
    if "epidemic diarrhea" in vl:
        return "PEDV", "PEDV_main"
    if "porcine respiratory coronavirus" in vl:
        return "PRCV", "PRCV_comparator"
    if "transmissible gastroenteritis" in vl or vl.startswith("tgev"):
        return "TGEV", "TGEV_main"
    # PDCoV 曾用名：Porcine coronavirus HKU15 / Coronavirus HKU15 / Swine deltacoronavirus
    if ("deltacoronavirus" in vl) or ("hku15" in vl) or ("pdcov" in vl):
        return "PDCoV", "PDCoV_main"
    # species suis 中的其他成员（CCoV/FCoV/FIPV/Alphacoronavirus 1/SeCoV）
    if any(k in vl for k in ("feline coronavirus", "feline infectious peritonitis",
                             "canine coronavirus", "alphacoronavirus 1",
                             "swine enteric coronavirus")):
        return v.split(" strain")[0].split(" isolate")[0], "other"
    return v.split(" strain")[0].split(" isolate")[0], "other"

def parse_strain(virus_name: str):
    v = (virus_name or "").strip()
    m = re.search(r"(?:strain|isolate)\s+(.+)$", v, re.I)
    if m:
        return m.group(1).strip()
    # 形如 "Porcine epidemic diarrhea virus L00721/GER/2014"：尾巴当作株系线索
    for pref in ("Porcine epidemic diarrhea virus", "Porcine deltacoronavirus",
                 "Deltacoronavirus PDCoV", "Swine deltacoronavirus",
                 "Porcine coronavirus HKU15", "TGEV", "Porcine respiratory coronavirus"):
        if v.startswith(pref) and len(v) > len(pref) + 1:
            return v[len(pref):].strip()
    return "Unknown"

def main():
    log("=== 阶段2 步骤1：已有54物种数据中目标病毒逐一核查 ===")
    log(f"核查对象文件：{FINAL_FASTA}")
    log(f"筛选关键字（病毒名层面）：epidemic diarrhea / transmissible gastroenteritis / "
        f"porcine respiratory coronavirus / deltacoronavirus / HKU15 / PDCoV / TGEV")

    # 1) FASTA 头核查（病毒全名关键词；头格式为 ">ACCESSION.VERSION 物种名"）
    headers = []
    with open(FINAL_FASTA, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                headers.append(line[1:].strip())
    log(f"FASTA 序列总数（去重后唯一基因组）：{len(headers)}")
    kw = re.compile(r"epidemic diarrhea|gastroenteritis|respiratory coronavirus|"
                    r"deltacoronavirus|porci|suis", re.I)
    # 头部只有物种拉丁名，需要结合 accession 清单与 metadata 判断
    hit_headers = [h for h in headers if kw.search(h)]
    log(f"头部含 porci/suis/deltacoronavirus 等物种名的序列数：{len(hit_headers)}")

    # 2) 最终 accession 清单
    final_set = set()
    with open(FINAL_LIST, encoding="utf-8") as fh:
        for line in fh:
            a = line.strip()
            if a:
                final_set.add(a)
    log(f"final_accession_list.txt 条目数：{len(final_set)}")

    # 3) 三个目标物种：已选 accession vs 最终保留
    rows = []
    for sid, (sp, meta_path) in SPECIES_FILES.items():
        sel_acc = []
        with open(SEL, encoding="utf-8") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f and f[0] == sid:
                    sel_acc.append(f[2])
        in_final = [a for a in sel_acc if a in final_set]
        not_in_final = [a for a in sel_acc if a not in final_set]
        log(f"物种#{sid} {sp}: 已选 {len(sel_acc)} 条，进入最终库 {len(in_final)} 条，"
            f"未进入 {len(not_in_final)} 条（下载失败/精确去重剔除，详见 dedup_groups.tsv）")
        # metadata
        with open(meta_path, encoding="utf-8") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            idx = {c: i for i, c in enumerate(hdr)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < len(hdr):
                    continue
                acc = f[idx["accession"]]
                rows.append({
                    "accession": acc,
                    "species": f[idx["species"]],
                    "virus_name": f[idx["virus_name"]],
                    "length": f[idx["length"]],
                    "host": f[idx["host"]],
                    "collection_date": f[idx["collection_date"]],
                    "geo_location": f[idx["geo_location"]],
                    "is_lab_host": f[idx["is_lab_host"]],
                    "is_vaccine_strain": f[idx["is_vaccine_strain"]],
                    "selected": f[idx["selected"]] if "selected" in idx else "True",
                    "in_final": acc in final_set,
                })
    log(f"三个目标物种 metadata 记录合计：{len(rows)}")

    # 4) 病毒名分类
    cls_count = {}
    for r in rows:
        v, g = classify(r["virus_name"])
        r["virus"], r["virus_group"] = v, g
        cls_count[(v, g)] = cls_count.get((v, g), 0) + 1
    for (v, g), c in sorted(cls_count.items(), key=lambda x: -x[1]):
        log(f"  病毒名分类: {v} [{g}] = {c} 条（metadata 层面，含未进入最终库者）")

    # 5) FASTA 序列级 QC（长度 25–34 kb、N 比例）
    seqinfo = {}
    cur_h, cur_len, cur_n = None, 0, 0
    want = {r["accession"] for r in rows}
    with open(FINAL_FASTA, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_h and cur_h in want:
                    seqinfo[cur_h] = (cur_len, cur_n)
                cur_h = line[1:].split()[0]
                cur_len, cur_n = 0, 0
            else:
                s = line.strip()
                cur_len += len(s)
                cur_n += s.upper().count("N")
        if cur_h and cur_h in want:
            seqinfo[cur_h] = (cur_len, cur_n)
    log(f"成功从 FASTA 取回序列信息的目标 accession 数：{len(seqinfo)} / {len(want)}")

    # 6) rep-54 代表基因组（判断 n1_status 可用的注释）
    rep_acc = {}
    if REP54.exists():
        with open(REP54, encoding="utf-8") as fh:
            rh = fh.readline().rstrip("\n").split("\t")
            ri = {c: i for i, c in enumerate(rh)}
            a_col = "accession" if "accession" in ri else rh[0]
            s_col = "species" if "species" in ri else (rh[1] if len(rh) > 1 else rh[0])
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= len(rh):
                    rep_acc[f[ri[a_col]]] = f[ri[s_col]]
        for a, sp in rep_acc.items():
            if "porci" in sp or sp == "Alphacoronavirus suis" or "Deltacoronavirus suis" in sp:
                log(f"  目标物种代表基因组：{a} ({sp})")
    log(f"rep-54 代表基因组总数：{len(rep_acc)}")

    # 7) 输出 TSV
    cols = ["virus", "virus_group", "species", "virus_name", "strain", "accession",
            "genome_length", "n_count", "n_fraction", "qc_status", "host",
            "collection_date", "geo_location", "is_lab_host", "is_vaccine_strain",
            "sequence_source", "sequence_available", "n1_status", "note"]
    kept, fail = 0, 0
    with open(OUT_TSV, "w", encoding="utf-8", newline="") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            if not r["in_final"]:
                continue
            if r["accession"] not in seqinfo:
                fail += 1
                continue
            ln, n = seqinfo[r["accession"]]
            nfrac = n / ln if ln else 1.0
            qc = "pass" if (25000 <= ln <= 34000 and nfrac <= 0.01) else "fail"
            if qc == "fail":
                fail += 1
            else:
                kept += 1
            n1_status = "pending_stage3"
            note = ""
            if r["accession"] == "PP999774.1":
                n1_status = "reference_annotated"
                note = "LZ202401 参考株，N1 坐标 26434..26772（阶段1已核验）"
            elif r["accession"] in rep_acc:
                n1_status = "orf_annotation_available"
                note = "rep-54 代表基因组，已有 ORF 注释表可用"
            out.write("\t".join([
                r["virus"], r["virus_group"], r["species"], r["virus_name"],
                parse_strain(r["virus_name"]), r["accession"], str(ln), str(n),
                f"{nfrac:.6f}", qc, r["host"], r["collection_date"],
                r["geo_location"], r["is_lab_host"], r["is_vaccine_strain"],
                "cov_existing", "yes", n1_status, note,
            ]) + "\n")
    log(f"QC pass：{kept} 条；QC fail 或缺失：{fail} 条")
    log(f"输出：{OUT_TSV}")
    log("=== 步骤1 完成 ===")
    logf.close()

if __name__ == "__main__":
    main()
