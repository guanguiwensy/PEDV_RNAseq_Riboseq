# Sensitive-content and implicit-dependency scan (final package)

Scan = regex search over every text file for: workstation/server paths and drive letters, user names, hostnames/IPs, ssh/scp/job-scheduler calls, e-mail addresses, credentials/keys/tokens (including typical key prefixes), personal names, internal wording; plus file-size/type checks and the runtime file-open audit of the isolated test run.

| Check | Result |
|---|---|
| BAM / SAM / FASTQ / CRAM / index files | none |
| Files > 1 MB | none (package < 4 MB) |
| Windows/Unix absolute paths of the workstation or servers in code | none. Original paths appear only as `<ORIG_*>` tokens in `docs/PATH_REWRITE_LOG.tsv` and `docs/path_rewrite_diffs/`. Remaining drive-letter-like hits are URLs (`https://...`), a regex in a shell script and escape sequences |
| User names, host names, IP addresses, e-mail addresses | none |
| Credentials, API keys, tokens | none (the word `token` occurs as a variable name in `phase2_step3_manifests.py`) |
| Personal names | removed (two docstrings now say "the PI"; logged) |
| Description of private infrastructure | the log message about a failed ssh login to a compute server was replaced by a neutral message (logged in `CURATION_PATCHES.tsv`); internal folder names removed from configuration text |
| Internal words | only `User-Agent: pedv-n1-*` strings (generic) and Chinese docstrings that mention task names; left unchanged to keep archived code identical to what produced the figures |
| Third-party data/code | none included. The one script of unclear provenance/licence (`plot_Fig5_recolored.py`) is **excluded** (`EXCLUDED_FROM_PUBLIC_PACKAGE.md`) |
| Unpublished data | no sequencing data; sample identifiers (R260…/XN260…) appear in metadata and code — they belong to BioProject PRJNA1554149, which is registered but not yet public (authors to confirm that publishing the IDs before release is acceptable) |

## Implicit dependency on the historical project directories
* `config/pedv_paths.py` has **no default pointing at any historical directory**: without environment variables or `config/paths.local.yaml` the roots are `<repo>/data_root/{P,K,R,cov}` (test T31: a script started without configuration fails on `<repo>/data_root/K/...` and opens nothing in the original tree).
* All isolated tests (`tests/isolated_run_results.tsv`) ran under a hook that records every file opened below the original project tree outside the sandbox: **0 reads** in all tests. The sandbox itself contained **no script files** (so sibling modules had to come from this repository).
* Remaining by-design coupling: scripts expect the *relative layout inside the data roots* (e.g. `GIb_rnaseq_reanalysis/figures/Fig3_work/`) and several read stored server-derived intermediates (listed in `IO_DEPENDENCY_MODULES.md` and `SERVER_SCRIPTS_REQUIRED.md`).
