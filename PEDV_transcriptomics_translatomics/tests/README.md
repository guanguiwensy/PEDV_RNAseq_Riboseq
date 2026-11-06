# Tests

- `repro_commands.sh` lists the analysis entry points used for reproduction. It requires populated `PEDV_P`, `PEDV_K`, and `PEDV_R` data roots.
- `audit_whitelist.py` scans a runtime file-access log and its referenced inputs for sample identifiers outside the registry.

Generated logs and historical packaging reports are intentionally excluded from the repository.
