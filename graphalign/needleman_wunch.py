import numpy as np
from collections import namedtuple
from itertools import chain

from .sequencegraph import get_prev_func, naive_graph


Alphabet = namedtuple("Alphabet", ["to_str", "to_num"])

_dna_chars = ["A", "C", "G", "T"]
DNAAlphabet = Alphabet(_dna_chars, {c: i for i, c in enumerate(_dna_chars)})


def get_score_mat(mismatch_score, alphabet_size=4):
    scores = mismatch_score*np.ones((alphabet_size, alphabet_size))
    for i in range(alphabet_size):
        scores[i, i] = 1
    return scores


def get_align_func(gap_open, score_matrix, gap_extend=None, use_graphs=True):
    if gap_extend is None:
        gap_extend = gap_open

    def get_comb_scores(seq_a, seq_b):
        return np.array([[score_matrix[c_a][c_b] for c_b in seq_b] for c_a in seq_a])

    def init_matrices(len_a, len_b):
        matrix = np.zeros((len_a+1, len_b+1))
        open_matrix_a = np.zeros((len_a+1, len_b+1))
        open_matrix_b = np.zeros((len_a+1, len_b+1))
        matrix[1:, 0] = gap_open+gap_extend*np.arange(len_a)
        open_matrix_a[1:, 0] = gap_open+gap_extend*np.arange(len_a)
        open_matrix_b[0, 1:] = gap_open+gap_extend*np.arange(len_b)
        matrix[0, 1:] = gap_open+gap_extend*np.arange(len_b)
        open_matrix_a[0, :] = -100
        open_matrix_b[:, 0] = -100
        return matrix, open_matrix_a, open_matrix_b

    def align(seq_a, seq_b):
        comb_scores = get_comb_scores(seq_a, seq_b)
        matrix, open_matrix_a, open_matrix_b = init_matrices(
            len(seq_a), len(seq_b))
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

    def graph_align(graph_a, graph_b):
        seq_a = graph_a.sequences
        seq_b = graph_b.sequences
        comb_scores = get_comb_scores(seq_a, seq_b)
        matrix, open_matrix_a, open_matrix_b = init_matrices(
            len(seq_a), len(seq_b))
        get_prev_a = get_prev_func(graph_a)
        get_prev_b = get_prev_func(graph_b)
        backtrack_matrices = np.zeros((3, len(seq_a)+1, len(seq_b)+1, 3), dtype="int")

        def get_best_linear(i, j):
            indels = list(chain(((prev_i, j) for prev_i in prev_is),
                                ((i, prev_j) for prev_j in prev_js)))
            indel_scores = [matrix[_i, _j]+gap_open for _i, _j in indels]
            matches = [(prev_i, prev_j) for prev_i in prev_is for prev_j in prev_js]
            match_scores = [matrix[_i, _j]+comb_scores[_i, _j] for _i, _j in matches]
            tmp = list(chain(zip(indel_scores, indels),
                             zip(match_scores, matches)))
            return max(tmp)

        def get_best_a(prev_is, j):
            new_open = ((matrix[prev_i, j] + gap_open, (prev_i, j, 0)) for prev_i in prev_is)
            old_open = ((open_matrix_a[prev_i, j] + gap_extend, (prev_i, j, 1)) for prev_i in prev_is)
            return max(chain(new_open, old_open))

        def get_best_b(i, prev_js):
            new_open = ((matrix[i, prev_j] + gap_open, (i, prev_j, 0)) for prev_j in prev_js)
            old_open = ((open_matrix_b[i, prev_j] + gap_extend, (i, prev_j, 2)) for prev_j in prev_js)
            return max(chain(new_open, old_open))

        for i in range(1, len(seq_a)+1):
            for j in range(1, len(seq_b)+1):
                prev_is = get_prev_a(i)
                prev_js = get_prev_b(j)
                score, ij = get_best_linear(i, j)
                ijk = (ij[0], ij[1], 0)
                a_score, a_ijk = get_best_a(prev_is, j)
                b_score, b_ijk = get_best_b(i, prev_js)
                open_matrix_a[i, j] = a_score
                backtrack_matrices[1, i, j] = a_ijk
                open_matrix_b[i, j] = b_score
                backtrack_matrices[2, i, j] = b_ijk
                m_score, m_ijk = max(((score, ijk),
                                     (a_score, a_ijk), (b_score, b_ijk)))
                backtrack_matrices[0, i, j] = m_ijk
                matrix[i, j] = m_score

#                 open_matrix_a[i, j] = max(max(matrix[prev_i, j] + gap_open,
#                                               open_matrix_a[prev_i, j]+gap_extend)
#                                           for prev_i in prev_is)
# 
#                 open_matrix_b[i, j] = max(max(matrix[i, prev_j]+gap_open,
#                                               open_matrix_b[i, prev_j] + gap_extend)
#                                           for prev_j in prev_js)
#                 matrix[i, j] = max(max_indel, max_match,
#                                    open_matrix_a[i, j], open_matrix_b[i, j])

        print(matrix)
        return matrix[-1, -1]

    if use_graphs:
        return lambda seq_a, seq_b: graph_align(naive_graph(seq_a),
                                                naive_graph(seq_b))
    return align
