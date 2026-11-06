"""Junction sequence helpers used by the cross-strain analysis."""


def junction_homology(sequence, donor, acceptor, max_nt=100):
    """Count identical bases to the right of two one-based junction boundaries."""
    sequence = sequence.upper()
    donor_index = int(donor)
    acceptor_index = int(acceptor) - 1
    limit = min(max_nt, len(sequence) - donor_index, len(sequence) - acceptor_index)
    length = 0
    while length < limit and sequence[donor_index + length] == sequence[acceptor_index + length]:
        length += 1
    return length
