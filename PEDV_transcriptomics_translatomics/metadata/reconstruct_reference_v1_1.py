#!/usr/bin/env python3
"""Rebuild the project reference `PEDV_GIb` (refdata_GIb_v1.1, 27,959 nt) from the public GenBank record KP162057.1.

Fact established during release preparation: v1.1 is colinear with KP162057.1 (same length, no indel) and differs at exactly 8 positions.
Positions are 1-based; the table lists KP162057.1 -> v1.1.

Usage:
    1. download KP162057.1 as FASTA from NCBI Nucleotide (efetch ... rettype=fasta) -> KP162057.1.fa
    2. python reconstruct_reference_v1_1.py KP162057.1.fa PEDV_GIb.fa
The script verifies the sequence of the input record and of the output against the SHA-256 values recorded when the reference was frozen.
Exit code 0 = reconstructed reference is identical to the one used for all analyses.
"""
import sys, hashlib

SNV = [(57, "C", "T"), (1700, "T", "C"), (2326, "T", "G"), (9599, "C", "T"), (20773, "A", "G"), (21574, "C", "T"), (24279, "G", "T"), (27959, "A", "C")]
KP_SEQ_SHA256 = "102ce18bee8bd6cb7fefbb0fd2f2b621cba8a520485ae8aa076b1bc830100491"   # KP162057.1, upper case, no whitespace
V11_SEQ_SHA256 = "fef2fcb2e00b0756a55499958f8eed8bc485e9c3a901409bd4271a9a7cf60baf"  # refdata_GIb_v1.1 sequence
V11_FILE_SHA256 = "48c36027c9b24ddcdbbec34537dea126ade378ad4094e5d72054f697c8412b92"  # PEDV_GIb.fa as used (header '>PEDV_GIb', 70 nt per line)


def read_fasta(p):
    seq = []
    for l in open(p):
        if not l.startswith(">"):
            seq.append(l.strip())
    return "".join(seq).upper()


def main(src, dst):
    s = read_fasta(src)
    h = hashlib.sha256(s.encode()).hexdigest()
    if len(s) != 27959 or h != KP_SEQ_SHA256:
        raise SystemExit("input is not the expected KP162057.1 sequence (len=%d sha=%s). Note: a later GenBank version would differ." % (len(s), h))
    t = list(s)
    for pos, old, new in SNV:
        assert t[pos - 1] == old, (pos, t[pos - 1], old)
        t[pos - 1] = new
    t = "".join(t)
    assert hashlib.sha256(t.encode()).hexdigest() == V11_SEQ_SHA256, "reconstructed sequence differs from refdata_GIb_v1.1"
    with open(dst, "w", newline="\n") as f:
        f.write(">PEDV_GIb\n")
        for i in range(0, len(t), 70):
            f.write(t[i:i + 70] + "\n")
    fh = hashlib.sha256(open(dst, "rb").read()).hexdigest()
    print("sequence identical to refdata_GIb_v1.1: OK; file sha256 %s (%s)" % (fh, "matches" if fh == V11_FILE_SHA256 else "differs only in file formatting"))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2])
