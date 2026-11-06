# tests/

* `results.tsv` – ledger of the first-round local tests T01-T22 (sandbox runs; superseded where different by the isolated re-run).
* `isolated_run_results.tsv`, `isolated_run_log.txt` – authoritative re-run from an independent copy of this repository, sandbox without any script files, with an audit hook for reads of the original tree.
* `repro_commands.sh` – the commands used for the tests. It expects a sandbox in which the needed inputs have been copied under `$PEDV_P`, `$PEDV_K`, `$PEDV_R`; comparison is by `cmp`/`diff` against the originals, which are never written to.
* **`FINAL_SAMPLE_WHITELIST_AUDIT.tsv`** (figure/panel level: aliases, SRR, inputs, recomputation, plot run, old-vs-new differences, unfinished steps), `FINAL_SAMPLE_WHITELIST_AUDIT_FILES.tsv` (content scan of every input file read), `chain_results.jsonl` (guarded rebuild chain log), `audit_whitelist.py`, `build_final_audit.py` — authoritative for the final 30-library sample set. The first-round ledger `results.tsv`, `isolated_run_results.tsv` and `TESTED_VS_RELEASED.txt` describe earlier packaging states and are kept as history.
* Result files themselves are not stored here (they would duplicate the data tables).
