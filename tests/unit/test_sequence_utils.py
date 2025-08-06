import random
import pytest
from dataset import sequence_utils


def test_extract_subseqs_basic():
    seq = "ACGT" * 100
    rng = random.Random(42)
    result = sequence_utils.extract_subseqs(
        seq, n=5, min_len=10, max_len=20, rng=rng)
    assert len(result) == 5
    for subseq in result:
        assert 10 <= len(subseq) <= 20


def test_extract_subseqs_reproducible():
    seq = "ACGT" * 100
    rng1 = random.Random(42)
    rng2 = random.Random(42)
    result1 = sequence_utils.extract_subseqs(
        seq, n=5, min_len=10, max_len=20, rng=rng1)
    result2 = sequence_utils.extract_subseqs(
        seq, n=5, min_len=10, max_len=20, rng=rng2)
    assert result1 == result2


def test_extract_subseqs_zero_n():
    seq = "ACGTACGTACGT"
    with pytest.raises(ValueError):
        sequence_utils.extract_subseqs(seq, n=0, min_len=3, max_len=5)


def test_extract_subseqs_minlen_gt_maxlen():
    seq = "ACGTACGTACGT"
    with pytest.raises(ValueError):
        sequence_utils.extract_subseqs(seq, n=1, min_len=10, max_len=5)


def test_extract_subseqs_sequence_shorter_than_minlen():
    seq = "ACG"
    result = sequence_utils.extract_subseqs(seq, n=1, min_len=5, max_len=10)
    assert result == []


def test_get_complement_all_iupac():
    inp = "ATGCYRWSKMBDHVN"
    expected = "NBDHVKMSWYRGCAT"
    assert sequence_utils.get_complement(inp) == expected


def test_get_complement_lowercase():
    inp = "atgcyrwskmbdhvn"
    expected = "NBDHVKMSWYRGCAT"
    assert sequence_utils.get_complement(inp) == expected


def test_get_complement_empty():
    assert sequence_utils.get_complement("") == ""


def test_get_complement_non_iupac():
    inp = "ATGCX"
    expected = "NGCAT"
    assert sequence_utils.get_complement(inp) == expected
