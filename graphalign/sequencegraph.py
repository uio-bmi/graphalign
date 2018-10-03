from collections import namedtuple, defaultdict
import numpy as np


SequenceGraph = namedtuple("SequenceGraph",
                           ["sequences", "node_offsets", "adj_list"])


def naive_graph(sequence):
    if isinstance(sequence, SequenceGraph):
        return sequence
    n_nodes = len(sequence)//3
    node_offsets = list(i*3 for i in range(n_nodes))
    adj_list = {i: [i+1] for i in range(n_nodes-1)}
    return SequenceGraph(sequence, node_offsets, adj_list)


def _get_reverse_adj_list(adj_list):
    reverse_adj_list = defaultdict(list)
    for from_node, to_nodes in adj_list.items():
        for to_node in to_nodes:
            reverse_adj_list[to_node].append(from_node)
    return reverse_adj_list


def get_prev_func(graph):
    node_ids = np.zeros(len(graph.sequences)+1)
    node_ids[graph.node_offsets] = 1
    node_ids = np.cumsum(node_ids)-1
    reverse_adj_list = _get_reverse_adj_list(graph.adj_list)

    def get_prev(i):
        if i == 0:
            return []
        if node_ids[i-1] != node_ids[i]:
            return [graph.node_offsets[prev_node+1]-1
                    for prev_node in reverse_adj_list[node_ids[i]]]
        return [i-1]

    return get_prev
