from Bio import pairwise2
import numpy as np
from graphalign import DNAAlphabet, get_align_func, get_score_mat, SequenceGraph
import pytest


def get_bio_align_func(gap_open, score_matrix, gap_extend=None):
    match_score = score_matrix[0, 0]
    mismatch_score = score_matrix[0, 1]

    def align_func(seq_a, seq_b):
        seq_a = "".join(DNAAlphabet.to_str[c] for c in seq_a)
        seq_b = "".join(DNAAlphabet.to_str[c] for c in seq_b)
        res = pairwise2.align.globalms(seq_a, seq_b, match_score, mismatch_score,
                                       gap_open, gap_extend)[0]
        return (res[0], res[1], res[2])

    return align_func


@pytest.fixture
def align_affine():
    return get_align_func(-3, get_score_mat(-1), -1, return_seq=True)


@pytest.fixture
def align_graph():
    return get_align_func(-3, get_score_mat(-1), -1, use_graphs=True, return_seq=True)


@pytest.fixture
def bio_align():
    return get_bio_align_func(-3, get_score_mat(-1), -1)


def test_align_snp_del(bio_align, align_affine):
    a = [1, 2, 0, 3, 2, 1, 0]
    b = [2, 0, 2, 2, 1, 0]
    assert bio_align(a, b) == align_affine(a, b)


def test_small(bio_align, align_affine):
    a = [DNAAlphabet.to_num[c] for c in 'TTTATGACCA']
    b = [DNAAlphabet.to_num[c] for c in 'ACA']
    assert bio_align(a, b) == align_affine(a, b)


@pytest.mark.skip
def test_random(bio_align, align_affine):
    np.random.seed(1000)
    for _ in range(100):
        a = list(np.random.randint(0, 4, 100))
        b = list(np.random.randint(0, 4, 90))
        assert bio_align(a, b) == align_affine(a, b)
