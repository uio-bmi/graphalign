from graphalign import get_align_func, get_score_mat


def test_align():
    a = [1, 2, 0, 3, 2, 1, 0]
    b = [1, 2, 0, 3, 2, 1, 0]
    align = get_align_func(-1, get_score_mat(-1))
    score = align(a,  b)
    assert score == len(a)
