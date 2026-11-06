# Server-side steps that were NOT re-run (and whose scripts must be supplied or confirmed)

Nothing here was executed or modified in this round. Paths are relative to the project's server working tree; user/host-specific prefixes were removed.

| Step | Needed by | Script / file (relative) | In this repo? | Output consumed locally |
|---|---|---|---|---|
| Ribo-seq trimming + STAR (oh40) + QC per library | Fig 4-6, Sup 10-13 | `analysis_batch23/Riboseq_qc_mapping/scripts/run_one.sh` | server only | BAMs, one-pass stats (P-sites, host metagene) |
| RNA-seq trimming + STAR (oh149) + minimap2 junction mapping | Fig 1-3, Sup 2-9 | `analysis_sgRNA_junction/scripts/run_qc_align.sh`; `rnaseq/01_alignment/run_align_GIb.sh` (in repo) | partially | `viral_fragments_mm2.tsv`, junction tables |
| Server denominator of mapped viral fragments | Sup 3, 6, 7 | `s2_server_denominator.sh` | in repo (`rnaseq/04_junction_dynamics`) | used as stored input in test T04 |
| Reference build (merged Vero + PEDV_GIb index) | all | `refbuild_GIb/server/*` | server only | STAR indices `..._oh149`, `..._oh40` |
| Python/Java environment for PRICE | Fig 6E/F, Sup 17 | `env_orf.sh` (server home, `tmp/`) | server only | — |
| PRICE host training + virus run (`-filter 29:29`) | Fig 6E/F, Sup 17 | PRICE 1.0.5 / Gedi 1.0.6f, parameters in `riboseq/04_price/README` | wrapper scripts in repo | `*.orfs.tsv` annotated outputs |
| HARTIS G1 (host per-position extraction), G2 (background model), G3 (start-level run) | Fig 6E/F | `HARTIS_v3` G1–G3 steps | G3/G4 table builders in repo; G1/G2 server only | stored G3 start-level table (used in T10) |
| Public RNA-seq processing (PRJNA1069035, PRJNA1203018) | Sup 4, 8, 9 | strain assembly + minimap2 scripts, `rnaseq/06_external_junction_comparison/` | partially | per-dataset junction tables |
| PK1 Ribo-seq mapping + pig-host P-site offset | Sup 18 | `PK1_psite` pipeline | server only | PK1 P-site long table |

## Request to the data owner
Before release please provide (or confirm unnecessary): the six "server only" items above, a `sha256sum` list for the stored intermediate inputs used in the tests (they are listed in `tests/results.tsv` by name), and, after release of BioProject PRJNA1554149, the run-level SRA accessions.
