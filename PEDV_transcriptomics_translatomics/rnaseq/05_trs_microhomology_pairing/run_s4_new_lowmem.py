"""Low-memory runner for the Fig. 3 span-stratified permutation analysis.

The runner replaces the large process pool in ``s4_new.py`` with a serial pool
or a pool of at most four workers. The statistical model, seed, and 10,000
permutation draws are unchanged.

Usage: python run_s4_new_lowmem.py [--workers 1]
Requires ViennaRNA Python bindings.
"""
import os, sys
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_k, "1")
_here = os.path.dirname(os.path.abspath(__file__))
workers = 1
if "--workers" in sys.argv:
    workers = int(sys.argv[sys.argv.index("--workers") + 1])
assert 1 <= workers <= 4, "refusing to start more than 4 processes"
import multiprocessing as _mp


class _SerialPool:
    def __init__(self, *a, **k): pass
    def __enter__(self): return self
    def __exit__(self, *e): return False
    def map(self, fn, it, chunksize=None): return [fn(x) for x in it]
    def imap(self, fn, it, chunksize=1): return (fn(x) for x in it)


_RealPool = _mp.Pool
_mp.Pool = (lambda *a, **k: _SerialPool()) if workers == 1 else (lambda *a, **k: _RealPool(workers))
sys.path.insert(0, _here)
import common  # unmodified module: SEED / NPERM definition
assert common.SEED == 42 and common.NPERM == 10000, (common.SEED, common.NPERM)
print("config: SEED=%d NPERM=%d workers=%d background=span-quartile-stratified pair-swap chain" % (common.SEED, common.NPERM, workers), flush=True)
sys.argv = [os.path.join(_here, "s4_new.py")]
import runpy
runpy.run_path(sys.argv[0], run_name="__main__")
