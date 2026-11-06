# Release blockers (after the final sample-whitelist correction)

Missing original plotting scripts, Illustrator edits and the absence of pixel-identical figures are documented limitations, not blockers.

## A. Must be decided/fixed before the results can be called "all traceable to the 30 libraries"
| # | Item | State |
|---|---|---|
| 1 | **Fig 5** | analysis recomputed with whitelisted libraries; independent plot exists but the submitted y-scale and the red/blue signal rule are not reproduced -> authors must confirm the rule and decide whether Fig 5 panels are redrawn (the CHX 16 h dependent values and 124 terminal-candidate calls changed, see the diff report) |
| 2 | **Sup 13** | RNA-depth based columns changed (RNA partners replaced by whitelisted libraries) -> check the panel/legend against the new tables |
| 3 | **Sup 4, 8, 9** | not re-run; the scripts that process unapproved external datasets are withheld -> needs re-implementation (approved external runs only) and server-processed reads |
| 4 | **Sup 17, Sup 18** | not re-run (server outputs); code restricted to whitelisted libraries |

## B. Decisions to record (not blockers)
* Licence: MIT (`LICENSE`, generic copyright line).
* BioProject `PRJNA1554149` is registered but not public (README states it); the PK1 libraries are in `metadata/sample.txt` (confirm accession scope before release).
* Fig 6B/6F wording, Fig 3 methods text, post-hoc Fig 6C/6D scripts: see `UNRESOLVED_ISSUES.md`.

## C. Documented limitations
No class C reproduction (raw reads -> counts); no pixel comparison; HARTIS G1/G2 and PRICE are server steps (rows of non-whitelisted libraries are dropped at load, per-library fits are independent); several stored multi-library inputs are filtered explicitly; no biological replicates.

## D. Verified in this round
* Registry + always-on path guard; 0 guard violations in the full rebuild chain; all stages run (see `tests/chain_results.jsonl`, `docs/REPRODUCTION_TEST_REPORT.md`).
* Fig 1-4, Fig 6, Sup 2-7, 10-12, 15, 16: recomputed values identical to the submitted ones.
* Sup 15 background totals reproduced exactly from whitelisted P-sites.
