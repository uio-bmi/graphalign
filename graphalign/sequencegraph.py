from collections import namedtuple

SequenceGraph = namedtuple("SequenceGraph",
                           ["sequences", "node_offsets", "adj_list"])


def get_prev_func(graph):
    return lambda i: [i-1]
