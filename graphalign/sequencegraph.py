from collections import namedtuple, defaultdict
import numpy as np


SequenceGraph = namedtuple("SequenceGraph",
                           ["sequences", "node_offsets", "adj_list"])

Alignment = namedtuple("Alignment", ["seq_a", "seq_b"])


def alignment_to_sequencegraph(alignment):
    seq_a, seq_b = alignment
    print(seq_a, seq_b)
    adj_list = defaultdict(set)
    node_offsets = []
    sequences = []
    iter_a, iter_b = (iter(seq_a), iter(seq_b))
    cur_a, cur_b = (next(iter_a), next(iter_b))
    print(list(zip(seq_a, seq_b)))
    pairs = iter(zip(seq_a, seq_b))
    cur_a, cur_b = next(pairs)
    i = 0
    prev_a, prev_b = (None, None)
    while True:
        cur_node = len(node_offsets)
        node_offsets.append(i)
        try:
            if cur_a == cur_b:
                adj_list[prev_a].add(cur_node)
                adj_list[prev_b].add(cur_node)
                prev_a = cur_node
                prev_b = cur_node
                while cur_a == cur_b:
                    print("MATCH",  cur_a, cur_b)
                    sequences.append(cur_a)
                    i += 1
                    cur_a, cur_b = next(pairs)
                continue
            elif cur_a == "-":
                adj_list[prev_a].add(cur_node)
                adj_list[prev_b].add(cur_node)
                prev_b = cur_node
                while cur_a == "-":
                    sequences.append(cur_b)
                    i += 1
                    cur_a, cur_b = next(pairs)
                continue
            elif cur_b == "-":
                adj_list[prev_a].add(cur_node)
                adj_list[prev_b].add(cur_node)
                prev_a = cur_node
                while cur_b == "-":
                    sequences.append(cur_a)
                    i += 1
                    cur_a, cur_b = next(pairs)
                continue
            else:
                adj_list[prev_a].add(cur_node)
                adj_list[prev_b].add(cur_node)
                adj_list[prev_a].add(cur_node+1)
                adj_list[prev_b].add(cur_node+1)
                prev_a = cur_node
                prev_b = cur_node+1
                sequences.append(cur_a)
                i += 1
                node_offsets.append(i)
                sequences.append(cur_b)
                i += 1
                cur_a, cur_b = next(pairs)
        except StopIteration:
            break
    del adj_list[None]
    return SequenceGraph(sequences, node_offsets, adj_list)


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
