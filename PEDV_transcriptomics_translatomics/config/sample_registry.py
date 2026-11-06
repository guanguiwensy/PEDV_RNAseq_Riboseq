"""Single sample registry of the PEDV analyses.

Only the 30 experimental libraries listed in ``metadata/sample.txt`` (BioProject PRJNA1554149) may enter RNA-seq / Ribo-seq results.
Real library numbers, SRR accessions and FASTQ names live ONLY in ``metadata/sample.txt``; analysis code uses the biological aliases, e.g.

    RNA_VERO_CHX_0.1MOI_16h   RNA_VERO_HAR_1MOI_8h   RIBO_VERO_CHX_1MOI_16h   RIBO_PK1_HAR_1MOI_16h

and resolves file names / table keys through this module (no renaming of FASTQ/BAM files is needed).

Guard: ``assert_whitelisted(ids, context)`` stops the program (SystemExit) if any experimental library that is not in the table
(or any unapproved external run) appears in an analysis input. It never falls back to old data.
"""
import csv
import os
import re

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
SAMPLE_TABLE = os.environ.get("PEDV_SAMPLE_TABLE") or os.path.join(_REPO, "metadata", "sample.txt")

# external public runs that the authors approved for the comparative RNA-seq analyses (Sup. Fig. 4/8/9): 3 piglet intestine runs (PRJNA1069035) + 1 LLC-PK run (PRJNA1203018)
APPROVED_EXTERNAL_RUNS = frozenset({"SRR27751926", "SRR27751927", "SRR27751928", "SRR31828424"})

# any token that looks like an experimental library number of this project
_LIB_PATTERN = re.compile(r"(?<![A-Za-z0-9])(R2[0-9]{7}|XN2[0-9]{7}|R_(?:CHX|HAR)_[0-9]+h|XN_(?:CHX|HAR)_[0-9]+h)(?![A-Za-z0-9])")
_SRR_PATTERN = re.compile(r"(?<![A-Za-z0-9])(SRR[0-9]{7,9})(?![A-Za-z0-9])")


class WhitelistError(SystemExit):
    pass


def _load():
    rows = list(csv.DictReader(open(SAMPLE_TABLE, encoding="utf-8"), delimiter="\t"))
    if len(rows) != 30:
        raise WhitelistError("sample table must contain exactly 30 libraries, found %d (%s)" % (len(rows), SAMPLE_TABLE))
    return rows


_ROWS = _load()
SAMPLES = {r["alias"]: r for r in _ROWS}                      # alias -> record
ALIASES = tuple(r["alias"] for r in _ROWS)
_LIB2ALIAS = {r["original_library_id"]: r["alias"] for r in _ROWS}
_SRR2ALIAS = {r["sra_run_accession"]: r["alias"] for r in _ROWS}
WHITELIST_LIBRARIES = frozenset(_LIB2ALIAS)
WHITELIST_SRR = frozenset(_SRR2ALIAS)


def lib(alias):
    """library number used as key in stored tables / file names (resolved from the registry)"""
    try:
        return SAMPLES[alias]["original_library_id"]
    except KeyError:
        raise WhitelistError("unknown sample alias %r (not in %s)" % (alias, SAMPLE_TABLE))


def srr(alias):
    return SAMPLES[alias]["sra_run_accession"]


def alias_of(library_id):
    try:
        return _LIB2ALIAS[library_id]
    except KeyError:
        raise WhitelistError("library %r is NOT in the final sample table -> refusing to continue" % library_id)


def find(assay=None, cell=None, treatment=None, moi=None, hpi=None):
    """aliases matching the given properties (assay 'RNA-seq'/'Ribo-seq', cell 'Vero'/'PK-1', treatment 'CHX'/'HAR', moi 0.1/1, hpi int)"""
    out = []
    for r in _ROWS:
        if assay and r["assay"] != assay: continue
        if cell and r["cell_line"] != cell: continue
        if treatment and r["treatment"] != treatment: continue
        if moi is not None and float(r["moi"]) != float(moi): continue
        if hpi is not None and int(r["hpi"]) != int(hpi): continue
        out.append(r["alias"])
    return out


def ribo_vero(treatment, hpi):
    """alias of the Vero Ribo-seq library (MOI 1)"""
    return "RIBO_VERO_%s_1MOI_%dh" % (treatment, hpi)


def rna_vero(treatment, moi, hpi):
    m = ("%g" % float(moi))
    return "RNA_VERO_%s_%sMOI_%dh" % (treatment, m, hpi)


# convenience: the six Vero Ribo-seq libraries in analysis order (CHX5, HAR5, CHX8, HAR8, CHX16, HAR16)
RIBO_VERO_ALIASES = tuple(ribo_vero(t, h) for h in (5, 8, 16) for t in ("CHX", "HAR"))
RNA_VERO_ALIASES = tuple(a for a in ALIASES if a.startswith("RNA_VERO_"))
# Vero RNA-seq partner of each Ribo-seq library (MOI 1, same drug and time)
RNA_PARTNER = {ribo_vero(t, h): rna_vero(t, 1, h) for h in (5, 8, 16) for t in ("CHX", "HAR")}
LIB_RIBO_VERO = tuple(lib(a) for a in RIBO_VERO_ALIASES)
LIB_RNA_VERO = tuple(lib(a) for a in RNA_VERO_ALIASES)


def scan_text(text):
    """experimental library numbers / SRR accessions found in a string that are not allowed"""
    bad = {m for m in _LIB_PATTERN.findall(text) if m not in WHITELIST_LIBRARIES}
    bad |= {m for m in _SRR_PATTERN.findall(text) if m not in WHITELIST_SRR and m not in APPROVED_EXTERNAL_RUNS}
    return sorted(bad)


def assert_whitelisted(ids, context=""):
    """stop if any identifier is an experimental library (or SRR run) outside the final table / approved external list"""
    bad = sorted({str(i) for i in ids if str(i) not in WHITELIST_LIBRARIES and str(i) not in WHITELIST_SRR and str(i) not in APPROVED_EXTERNAL_RUNS
                  and (_LIB_PATTERN.fullmatch(str(i)) or _SRR_PATTERN.fullmatch(str(i)))})
    if bad:
        raise WhitelistError("NON-WHITELISTED sample(s) in analysis input %s: %s. Regenerate the upstream table from the final sample set; "
                             "there is no fallback to old data." % (("(" + context + ")") if context else "", ", ".join(bad)))


def assert_table_whitelisted(df, column, context=""):
    """pandas helper: the values of a sample column must all be whitelisted (and non-empty set)"""
    vals = set(df[column].astype(str))
    assert_whitelisted(vals, context or column)
    return df


def assert_columns_whitelisted(columns, context=""):
    """wide tables whose column names are library numbers"""
    assert_whitelisted([c for c in columns], context)


def assert_exactly(ids, aliases, context=""):
    """the identifiers must be exactly the libraries of the given aliases"""
    want = {lib(a) for a in aliases}
    got = set(str(i) for i in ids)
    assert_whitelisted(got, context)
    if got != want:
        raise WhitelistError("%s: expected libraries %s, found %s" % (context, sorted(want), sorted(got)))


# ---------------------------------------------------------------------------------------------------------------------------------
# Runtime path guard (always on): any attempt to OPEN a file whose path contains an experimental library number (or SRR run) that is not
# in the final table stops the program. Content-level checks are done by assert_whitelisted()/assert_exactly() in the loaders and by
# tests/audit_whitelist.py. There is no fallback and no switch to disable the guard.
# ---------------------------------------------------------------------------------------------------------------------------------
import sys as _sys

_GUARD_INSTALLED = False


def _path_hook(event, args):
    if event != "open":
        return
    p = args[0]
    if not isinstance(p, (str, bytes, os.PathLike)):
        return
    try:
        s = os.fsdecode(p)
    except Exception:
        return
    bad = scan_text(s.replace("\\", "/"))
    if bad:
        raise WhitelistError("refusing to open %r: path contains non-whitelisted sample(s) %s (final sample table: %s)" % (s, ", ".join(bad), SAMPLE_TABLE))


def install_path_guard():
    global _GUARD_INSTALLED
    if not _GUARD_INSTALLED:
        _sys.addaudithook(_path_hook)
        _GUARD_INSTALLED = True


install_path_guard()


def _cli(argv):
    """shell helper:  python sample_registry.py libs ribo_vero | libs rna_vero | libs all | lib ALIAS | srr ALIAS | table"""
    if not argv or argv[0] == "table":
        for a in ALIASES:
            print(a, lib(a), srr(a))
    elif argv[0] == "libs":
        sel = {"ribo_vero": LIB_RIBO_VERO, "rna_vero": LIB_RNA_VERO, "all": tuple(lib(a) for a in ALIASES)}[argv[1]]
        print(" ".join(sel))
    elif argv[0] == "lib":
        print(lib(argv[1]))
    elif argv[0] == "srr":
        print(srr(argv[1]))
    else:
        raise SystemExit(_cli.__doc__)


if __name__ == "__main__":
    import sys
    _cli(sys.argv[1:])
