"""
sequence_utils.py

Summary
-------
Utilities for DNA sequence manipulation and random subsequence extraction,
designed for dataset creation and taxonomic labeling workflows.

Functions
---------
extract_subseqs
    Extracts (random) subsequences from a given DNA sequence, with optional
    reproducibility.
get_complement
    Returns the reverse complement of a DNA sequence, including ambiguous
    IUPAC bases.

Notes
-----
Designed for large-scale, robust, and reproducible DNA data workflows.
"""
import random
import math


def extract_subseqs(
    seq: str,
    n: int,
    min_len: int,
    max_len: int,
    rng: random.Random | None = None
) -> list[str]:
    """
    Extract n subsequences from the input DNA sequence, with length between
    min_len and max_len.

    The extraction strategy depends on the length of the input sequence:

    - If `len(seq) >= 2 * n * max_len`: Extract n non-overlapping random
      windows of max_len.
    - If `n * max_len <= len(seq) < 2 * n * max_len`: Extract n non-adjacent,
      nearly non-overlapping windows of max_len.
    - If `len(seq) < n * max_len`: Extract overlapping subsequences of random
      length between min_len and max_len, adding reverse complements as needed
      for diversity.

    Parameters
    ----------
    seq : str
        Input DNA sequence.
    n : int
        Number of subsequences to extract.
    min_len : int
        Minimum subsequence length.
    max_len : int
        Maximum subsequence length.
    rng : random.Random, optional
        Random number generator for reproducibility. Defaults to Python's
        random module.

    Returns
    -------
    list of str
        List of extracted subsequences.

    Notes
    -----
    The function tries to maximize sequence diversity and minimize overlap
    unless sequence length requires otherwise.
    When extracting overlapping subsequences, both the sequence and its reverse
    complement are added, if not already present.

    Examples
    --------
    >>> import random
    >>> seq = "ATGCGTACGTTAGCTAGCTAGCGTACGTTAGC" * 5
    >>> rng = random.Random(42)
    >>> subs = extract_subseqs(seq, n=10, min_len=10, max_len=20, rng=rng)
    >>> print(subs[0])
    'ATGCGTACGT...'
    """
    if n <= 0:
        raise ValueError("n must be positive")
    if min_len > max_len:
        raise ValueError("min_len must be <= max_len")
    if isinstance(rng, int):
        rng = random.Random(rng)
    elif rng is None:
        rng = random
    subseqs = []
    if len(seq) < min_len:
        return []
    # extract non overlapping subsequences of max_len randomly
    if len(seq) >= 2 * n * max_len:
        blacklist = []
        while len(subseqs) < n:
            idx = rng.randrange(0, len(seq) - max_len + 1)
            try_again = False
            for pos in blacklist:
                if idx >= pos and idx <= pos + max_len:
                    try_again = True
                    break
            if try_again:
                continue
            else:
                blacklist.append(idx)
                subseqs.append(seq[idx: idx + max_len])
    # extract non adjacent overlapping subsequences of max_len (window between)
    elif len(seq) < 2 * n * max_len and len(seq) >= n * max_len:
        rest = (len(seq) / max_len) - n
        rest_bases = int(math.floor(len(seq) / rest))
        window_bases = int((rest_bases / n) - 1)
        left_start = 0
        right_start = len(seq) - max_len
        if n % 2 == 0:
            operations = int(n / 2)
        else:
            operations = int((n - 1) / 2)
        for i in range(operations):
            subseqs.append(seq[left_start: left_start + max_len])
            left_start += max_len + window_bases
            subseqs.append(seq[right_start: right_start + max_len])
            right_start -= max_len + window_bases
        if n % 2 != 0:
            mid_seq = int(math.floor(len(seq) / 2))
            mid_max_len = int(math.floor(max_len / 2))
            subseqs.append(seq[mid_seq - mid_max_len: mid_seq + mid_max_len])
    # extract overlapping subsequences of length between min_len and max_len
    else:
        iter = 0
        iter_limit = n * 2
        while len(subseqs) < n:
            iter += 1
            if iter >= iter_limit:
                max_len_subseq = max(subseqs, key=len)
                max_len_subseqs = (i for i, x in enumerate(subseqs)
                                   if len(x) == max_len_subseq)
                subseqs.pop(max_len_subseqs.pop())
            subseq_len = rng.randrange(min_len, min(max_len, len(seq)))
            idx = rng.randrange(0, len(seq) - subseq_len + 1)
            subseq = seq[idx: (idx + subseq_len)]
            complement = get_complement(subseq)
            if subseq not in subseqs:
                subseqs.append(subseq)
            if complement not in subseqs:
                subseqs.append(complement)
        # compensate double append when n is odd
        if len(subseqs) > n:
            subseqs.pop()
    return subseqs


def get_complement(seq: str) -> str:
    """
    Return the reverse complement sequence of a DNA string using IUPAC
    notation.

    Parameters
    ----------
    seq : str
        Input DNA sequence (A, T, C, G, or IUPAC ambiguity codes).

    Returns
    -------
    str
        Reverse complement sequence.

    Notes
    -----
    Ambiguous bases (e.g., Y, R, N) are handled according to IUPAC standards.
    The returned sequence is always reversed relative to input.

    Example
    -------
    >>> get_complement('ATGCY')
    'RGCAT'
    """
    # IUPAC notation https://doi.org/10.1021/bi00822a023
    complementary_bases = {
        "A": "T", "T": "A", "C": "G", "G": "C",
        "Y": "R", "R": "Y", "W": "W", "S": "S",
        "K": "M", "M": "K", "D": "H", "H": "D",
        "V": "B", "B": "V", "N": "N"
    }
    return "".join(
        [complementary_bases.get(i.upper(), "N") for i in seq]
    )[::-1]
