"""Sequence and coordinate helpers for external junction comparisons."""

from pathlib import Path

import numpy as np
from Bio import Align


CANON = {"S": 20597, "ORF3": 24718, "E": 25342, "M": 25597, "N": 26285}
LEADER_WIN = (55, 72)


def read_fa(path):
    """Return the uppercase sequence from a single-record FASTA file."""
    with Path(path).open(encoding="utf-8") as handle:
        return "".join(line.strip() for line in handle if not line.startswith(">")).upper()


def cigar_map(query_fasta, reference_fasta):
    """Build one-based query-to-reference and reference-to-query coordinate maps."""
    query = read_fa(query_fasta)
    reference = read_fa(reference_fasta)
    aligner = Align.PairwiseAligner(mode="fogsaa")
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -1
    alignment = aligner.align(query, reference)[0]

    q2r = np.zeros(len(query) + 2, dtype=np.int64)
    r2q = np.zeros(len(reference) + 2, dtype=np.int64)
    query_blocks, reference_blocks = alignment.aligned
    for (q0, q1), (r0, r1) in zip(query_blocks, reference_blocks):
        if q1 - q0 != r1 - r0:
            raise ValueError("unexpected unequal aligned blocks")
        size = q1 - q0
        qpos = np.arange(q0 + 1, q0 + size + 1)
        rpos = np.arange(r0 + 1, r0 + size + 1)
        q2r[qpos] = rpos
        r2q[rpos] = qpos
    return q2r, r2q
