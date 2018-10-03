from fixtures import *
# import pytest


def test_align(align):
    a = [1, 2, 0, 3, 2, 1, 0]
    b = [1, 2, 0, 3, 2, 1, 0]
    score = align(a,  b)
    assert score == len(a)


def test_align_snp(align):
    a = [1, 2, 0, 3, 2, 1, 0]
    b = [1, 2, 0, 2, 2, 1, 0]
    score = align(a,  b)
    assert score == len(a)-2


def test_align_del(align):
    a = [1, 2, 0, 3, 2, 1, 0]
    b = [1, 2, 0, 2, 1, 0]
    score = align(a,  b)
    assert score == len(b)-1


def test_align_snp_del(align):
    a = [1, 2, 0, 3, 2, 1, 0]
    b = [   2, 0, 2, 2, 1, 0]
    score = align(a,  b)
    assert score == len(b)-3
