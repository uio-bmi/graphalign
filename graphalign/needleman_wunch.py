import numpy as np
from collections import namedtuple

Alphabet = namedtuple("Alphabet", ["to_str", "to_num"])

_dna_chars = ["A", "C", "G", "T"]
DNAAlphabet = Alphabet(_dna_chars, {c: i for i, c in enumerate(_dna_chars)})


def get_score_mat(mismatch_score, alphabet_size=4):
    scores = mismatch_score*np.ones((alphabet_size, alphabet_size))
    for i in range(alphabet_size):
        scores[i, i] = 1
    return scores


def get_align_func(gap_open, score_matrix, gap_extend=None):
    if gap_extend is None:
        gap_extend = gap_open

    def get_comb_scores(seq_a, seq_b):
        return np.array([[score_matrix[c_a][c_b] for c_b in seq_b] for c_a in seq_a])

    def align(seq_a, seq_b):
        comb_scores = get_comb_scores(seq_a, seq_b)
        matrix = np.zeros((len(seq_a)+1, len(seq_b)+1))
        open_matrix_a = np.zeros((len(seq_a)+1, len(seq_b)+1))
        open_matrix_b = np.zeros((len(seq_a)+1, len(seq_b)+1))
        matrix[1:, 0] = gap_open+gap_extend*np.arange(len(seq_a))
        open_matrix_a[1:, 0] = gap_open+gap_extend*np.arange(len(seq_a))
        open_matrix_a[0, :] = -100
        open_matrix_b[:, 0] = -100
        open_matrix_b[0, 1:] = gap_open+gap_extend*np.arange(len(seq_b))

        matrix[0, 1:] = gap_open+gap_extend*np.arange(len(seq_b))
        for i in range(1, len(seq_a)+1):
            for j in range(1, len(seq_b)+1):
                scores = [matrix[i-1, j]+gap_open,
                          matrix[i, j-1]+gap_open,
                          matrix[i-1, j-1]+comb_scores[i-1, j-1]]
                open_matrix_a[i, j] = max(matrix[i-1, j]+gap_open,
                                          open_matrix_a[i-1, j]+gap_extend)
                open_matrix_b[i, j] = max(matrix[i, j-1]+gap_open,
                                          open_matrix_b[i, j-1] + gap_extend)
                matrix[i, j] = max(max(scores), open_matrix_a[i, j], open_matrix_b[i, j])

        print(matrix)
        return matrix[-1, -1]

    return align
