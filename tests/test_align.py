from graph_align import align


def test_align():
    a = "ACGTACGT"
    b = "ACGTACGT"
    score = align(a, b)
    assert score == len(a)
