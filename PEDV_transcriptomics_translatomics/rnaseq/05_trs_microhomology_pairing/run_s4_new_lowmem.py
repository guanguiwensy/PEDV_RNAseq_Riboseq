"""Low-memory / low-process runner for the FULL Fig 3 span-stratified permutation analysis (s4_new.py).

Added during release curation (NOT part of the original analysis code).  Motivation: the original s4_new.py opens
`Pool(26)`; on a 16 GB Windows laptop this exhausted memory / handles (multiprocessing pickling error) and the run crashed.

What this runner changes: ONLY the process model.  `multiprocessing.Pool` is replaced by an in-process (serial) pool or,
with --workers N (N<=4), a small real pool.  Nothing else is touched: background model (span-quartile-stratified pair-swap
chain, 'sweeps=2, burn=30'), random seed (SEED=42 from common.py), number of retained permutation draws (NPERM=10,000)
and all statistics are executed by the unmodified s4_new.py / metrics_mfe.py.  The runner asserts SEED==42 and NPERM==10000
and prints them, so a changed configuration cannot go unnoticed.

Cost (measured on the dev laptop, ~0.17 ms per RNA.duplexfold call): ~1.1 M donor x acceptor pairs x 2 modes -> a few minutes
for the energy matrices, plus the swap chains; peak memory < 1 GB with the default --workers 1.
Usage:  python run_s4_new_lowmem.py [--workers 1]
Requires ViennaRNA python bindings (import RNA, 2.7.2 used for the submitted numbers).
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
