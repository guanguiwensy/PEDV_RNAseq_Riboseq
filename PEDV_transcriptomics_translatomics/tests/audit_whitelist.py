#!/usr/bin/env python3
"""Final sample-whitelist audit (automatic part).

Input : the per-script file-access log of a guarded run (tests/chain_results.jsonl: for every script the data files it read/wrote, recorded by the
        runtime path guard) and the data roots (PEDV_P / PEDV_K).
Check : (1) the runtime guard never fired (a non-whitelisted library in a file PATH stops the program);
        (2) every data file that was READ is scanned for non-whitelisted experimental library numbers / SRR runs INSIDE the file
            (text and .gz tables). Files with hits are reported with the hit counts. A hit does not by itself mean contamination: row-separable
            inputs (one row per library) are filtered by the script right after loading and the script asserts the remaining set; the audit
            therefore reports 'non_whitelist_content' and the script-level filter separately (see FIGURES below).
Output: tests/FINAL_SAMPLE_WHITELIST_AUDIT_FILES.tsv (one row per script x input file) and the figure-level table is written by
        tests/build_final_audit.py from this plus the manual figure map.
Usage : PEDV_P=... PEDV_K=... python tests/audit_whitelist.py tests/chain_results.jsonl tests/FINAL_SAMPLE_WHITELIST_AUDIT_FILES.tsv
Exit  : 1 if the guard fired in any script, 0 otherwise.
"""
import os, sys, json, gzip, csv
_d = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(_d), "config"))
import sample_registry as SR

ROOTS = {"P": os.environ.get("PEDV_P", ""), "K": os.environ.get("PEDV_K", ""), "R": os.environ.get("PEDV_R", "")}
MAXMB = 400


def resolve(rel):
    rel = rel.replace("\\", "/")
    parts = rel.strip("/").split("/", 1)
    parts[0] = parts[0].upper()   # the guard log stores lower-cased Windows paths
    if parts[0] in ROOTS and ROOTS[parts[0]]:
        return os.path.join(ROOTS[parts[0]], parts[1]) if len(parts) > 1 else ROOTS[parts[0]]
    return None


_CACHE = {}


def scan_file(p):
    if p in _CACHE:
        return _CACHE[p]
    hits = {}
    try:
        if not os.path.isfile(p) or os.path.getsize(p) > MAXMB * 1e6:
            res = None
        else:
            if p.endswith(".gz"):
                f = gzip.open(p, "rt", encoding="utf-8", errors="replace")
            elif p.endswith((".tsv", ".csv", ".txt", ".json", ".md", ".fa", ".gtf", ".bed")):
                f = open(p, encoding="utf-8", errors="replace")
            else:
                f = None   # binary (pickle, images): not scanned
            if f is None:
                res = None
            else:
                with f:
                    for line in f:
                        for b in SR.scan_text(line):
                            hits[b] = hits.get(b, 0) + 1
                res = hits
    except Exception as e:
        res = None
    _CACHE[p] = res
    return res


def main(log, out):
    rows = []
    guard_fired = 0
    for l in open(log, encoding="utf-8"):
        r = json.loads(l)
        for d in r["detail"]:
            if d["forbidden"]:
                guard_fired += 1
            for rel in d["reads"]:
                p = resolve(rel)
                if p is None:
                    continue
                h = scan_file(p)
                rows.append([r["test"], d["script"], rel, "binary_or_unscanned" if h is None else ("non_whitelist_content" if h else "clean"),
                             ";".join("%s:%d" % kv for kv in sorted((h or {}).items()))])
    with open(out, "w", encoding="utf-8", newline="") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n")
        w.writerow(["stage", "script", "input_file(root-relative)", "content_status", "non_whitelist_ids:count"])
        w.writerows(rows)
    n_hit = sum(1 for x in rows if x[3] == "non_whitelist_content")
    print("scripts with guard violations:", guard_fired, "| input files scanned:", len(rows), "| files with non-whitelist content:", n_hit)
    return 1 if guard_fired else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
