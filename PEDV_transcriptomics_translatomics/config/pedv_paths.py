"""Resolve the repository's configurable data roots."""

import os

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))).replace("\\", "/")


def _load_yaml():
    p = os.path.join(REPO_ROOT, "config", "paths.local.yaml")
    d = {}
    if os.path.exists(p):
        with open(p, encoding="utf-8") as f:
            for line in f:
                line = line.split("#", 1)[0].strip()
                if ":" in line:
                    k, v = line.split(":", 1)
                    d[k.strip()] = v.strip().strip("'\"")
    return d


_Y = _load_yaml()


def _get(name, default):
    v = os.environ.get(name) or _Y.get(name) or default
    return v.replace("\\", "/")


def _slash(v):
    return v if v.endswith("/") else v + "/"


PEDV_P = _slash(_get("PEDV_P", REPO_ROOT + "/data_root/P"))
PEDV_K = _slash(_get("PEDV_K", REPO_ROOT + "/data_root/K"))
PEDV_R = _slash(_get("PEDV_R", REPO_ROOT + "/data_root/R"))
PEDV_COV = _slash(_get("PEDV_COV", REPO_ROOT + "/data_root/cov"))
PEDV_SERVER_ROOT = _get("PEDV_SERVER_ROOT", REPO_ROOT + "/data_root/server").rstrip("/")


def _script_dirs():
    out = []
    for dp, dn, fn in os.walk(REPO_ROOT):
        dn[:] = [d for d in dn if d not in ("config", "tests", "docs", "__pycache__", ".git")]
        if any(f.endswith(".py") for f in fn):
            out.append(dp.replace("\\", "/"))
    return sorted(out)


SCRIPT_DIRS = _script_dirs()


def _origin():
    """Map repository scripts to their data-tree working locations."""
    d = {}
    p = os.path.join(REPO_ROOT, "config", "script_origin.tsv")
    roots = {"P": PEDV_P, "K": PEDV_K, "R": PEDV_R}
    if os.path.exists(p):
        with open(p, encoding="utf-8") as f:
            for line in f.read().splitlines()[1:]:
                a = line.split("	")
                if len(a) == 3:
                    d[a[0]] = roots[a[1]] + a[2]
    return d


ORIGIN = _origin()


def check_sup16_featured(feat):
    """Sup 16: the 9 featured families are a FIXED input (config/sup16_featured_families.json); abort if a working file disagrees."""
    import json
    cfg = json.load(open(os.path.join(REPO_ROOT, "config", "sup16_featured_families.json"), encoding="utf-8"))
    fixed = cfg["featured_family_ids"]
    got = list(feat["featured"])
    if got != fixed:
        raise SystemExit("work/featured.json['featured'] differs from the fixed Sup 16 list (config/sup16_featured_families.json): "
                         "run p6_aggregate.py before s1_C12_featured_and_data.py, never edit the list by hand. got=%s fixed=%s" % (got, fixed))
    return fixed
