import numpy as np


def get_score_mat(mismatch_score, alphabet_size=4):
    scores = mismatch_score*np.ones((alphabet_size, alphabet_size))
    for i in range(alphabet_size):
        scores[i, i] = 1
    return scores


def get_align_func(gap_open, score_matrix, gap_extend=None):
    def get_comb_scores(seq_a, seq_b):
        return np.array([[score_matrix[c_a][c_b] for c_b in seq_b] for c_a in seq_a])

    def align(seq_a, seq_b):
        comb_scores = get_comb_scores(seq_a, seq_b)
        matrix = np.zeros((len(seq_a)+1, len(seq_b)+1))
        matrix[:, 0] = -gap_open*np.arange(len(seq_a)+1)
        matrix[0, :] = -gap_open*np.arange(len(seq_b)+1)
        for i in range(1, len(seq_a)+1):
            for j in range(1, len(seq_b)+1):
                scores = [matrix[i-1, j]+gap_open, matrix[i, j-1]+gap_open, matrix[i-1, j-1]+comb_scores[i-1, j-1]]
                matrix[i, j] = max(scores)
        print(matrix)
        return matrix[-1, -1]

    return align
