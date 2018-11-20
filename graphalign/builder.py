from collections import defaultdict, namedtuple
from itertools import chain
import numpy as np
from .sequencegraph import SequenceGraph
from .datastructs import Interval, Position
from .topological_sort import topological_sort
import logging
BuilderStruct = namedtuple("BS", ["nodes", "adj_list"])


def adj_integrity(func):
    def new_func(self, *args, **kwargs):
        ret = func(self, *args, **kwargs)
        for from_node, to_nodes in self._adj_list.items():
            assert all(from_node in self._reverse_adj_list[to_node]
                       for to_node in to_nodes), (self._adj_list, self._reverse_adj_list)
        for from_node, to_nodes in self._reverse_adj_list.items():
            assert all(from_node in self._adj_list[to_node]
                       for to_node in to_nodes), (self._adj_list, self._reverse_adj_list)
        assert all(len(v) == len(set(v)) for v in self._adj_list.values())
        assert all(len(v) == len(set(v)) for v in self._reverse_adj_list.values())
        return ret
    return new_func


class GraphBuilder:

    @adj_integrity
    def __init__(self, nodes={}, adj_list={}):
        self._nodes = nodes
        self._adj_list = defaultdict(list, adj_list)
        self._reverse_adj_list = defaultdict(list)
        [self._reverse_adj_list[to_node].append(from_node)
         for from_node, to_nodes in self._adj_list.items()
         for to_node in to_nodes]
        print(self._reverse_adj_list)
        self._max_node = max(self._nodes) if self._nodes else -1

    def __eq__(self, other):
        return self._nodes == other._nodes and self._adj_list == other._adj_list

    def __repr__(self):
        pass

    def __add__(self, other):
        node_offset = self._max_node+1
        nodes = {node_offset+node_id: seq
                 for node_id, seq in other._nodes.items()}
        nodes.update(self._nodes)
        adj_list = {node_offset+node_id: [node_offset+n for n in nexts]
                    for node_id, nexts in other._adj_list}
        adj_list.update(self._adj_list)
        return self.__class__(nodes, adj_list)

    def size(self):
        return self._max_node + 1

    def add_node(self, data, node_id=None):
        if node_id is None:
            node_id = self._max_node+1
        self._nodes[node_id] = data
        self._max_node = max(node_id, self._max_node)
        return node_id

    @adj_integrity
    def add_edge(self, from_node, to_node):
        self._adj_list[from_node].append(to_node)
        self._reverse_adj_list[to_node].append(from_node)

    @adj_integrity
    def set_edge(self, from_node, to_node):
        logging.info("Setting edge %s to %s", from_node, to_node)
        for node in self._adj_list[from_node]:
            self._reverse_adj_list[node].remove(from_node)
        self._adj_list[from_node] = [to_node]
        self._reverse_adj_list[to_node] = [from_node]

    @adj_integrity
    def add_edges(self, from_nodes, to_nodes):
        logging.debug("Adding edges: %s -> %s", from_nodes, to_nodes)
        for node in from_nodes:
            self._adj_list[node].extend(n for n in to_nodes if n not in self._adj_list[node])
        for node in to_nodes:
            self._reverse_adj_list[node].extend(n for n in from_nodes if n not in self._reverse_adj_list[node])

    @adj_integrity
    def remove_node(self, node_id):
        del self._nodes[node_id]
        from_nodes = self._reverse_adj_list[node_id]
        to_nodes = self._adj_list[node_id]
        logging.debug("Removing node %s from lists %s, %s", node_id, from_nodes, to_nodes)
        for from_node in from_nodes:
            assert node_id in self._adj_list[from_node], (from_node, self._adj_list[from_node], node_id)
            self._adj_list[from_node].remove(node_id)

        for to_node in to_nodes:
            assert node_id in self._reverse_adj_list[to_node], (to_node, self._reverse_adj_list[to_node], node_id)
            self._reverse_adj_list[to_node].remove(node_id)
        del self._adj_list[node_id]
        del self._reverse_adj_list[node_id]

    def node_size(self, node_id):
        return len(self._nodes[node_id])

    def add_sequence_graph(self, sequence_graph):
        seqs = [sequence_graph.sequences[start:end] for start, end in
                zip(sequence_graph.node_offsets,
                    chain(sequence_graph.node_offsets, [len(sequence_graph)]))]
        prev_max = self._max_node
        [self.add_node(seq) for seq in seqs]
        [self.add_edge([from_node+prev_max], [to_node+prev_max for to_node in to_nodes])
         for from_node, to_nodes in sequence_graph._adj_list.items()]
        return self._max_node

    def to_sequence_graph(self):
        node_sequence = topological_sort(self)
        lookup = {node: i for i, node in enumerate(node_sequence)}
        sequences = [self._nodes[node_id] for node_id in node_sequence]
        node_offsets = chain([0], np.cumsum([len(seq) for seq in sequences])[:-1])
        adj_list = {lookup[from_node]: [lookup[to] for to in to_nodes]
                    for from_node, to_nodes in self._adj_list.items()}
        return SequenceGraph(list(chain(sequences)), node_offsets, adj_list)

    def to_struct(self):
        return BuilderStruct(self._nodes, {k: list(sorted(v)) for k, v in self._adj_list.items() if v})

    def to_topology(self):
        topology = defaultdict(list)
        for node_id, seq in self._nodes.items():
            if not self._adj_list[node_id]:
                continue
            topology[tuple(seq)].append(
                tuple(sorted(self._nodes[n] for n in self._adj_list[node_id])))
        for v in topology.values():
            v.sort()
        return topology

    def add_snp(self, pos, seq):
        node_id, offset = pos
        next_node = None
        if offset > 0:
            node_id = self._split_node(node_id, offset)
        if self.node_size(node_id) > 1:
            next_node = self._split_node(node_id, 1)
        self._copy_node(node_id, [seq])
        return next_node

    def add_snps(self, node_id, offsets, seqs):
        cur_offset = 0
        for offset, snp in zip(offsets, seqs):
            node_id = self.add_snp(Position(node_id, offset-cur_offset),
                                   snp)
            cur_offset = offset+1

    def add_deletion(self, interval):
        start_node = interval.nodes[0]
        end_node = interval.nodes[-1]
        if interval.end < self.node_size(end_node):
            self._split_node(end_node, interval.end)
        if interval.start > 0:
            start_node = self._split_node(start_node, interval.start)
        if end_node == interval.nodes[0]:
            end_node = start_node
        self.add_edges(self._reverse_adj_list[start_node],
                       self._adj_list[end_node])

    def _is_valid_position(self, position):
        return 0 <= position.offset < self.node_size(position.node)

    def add_insertion(self, position, sequence):
        """Add insertion before position"""
        assert self._is_valid_position(position), position
        new_node = self.add_node(sequence)
        node = position.node
        if position.offset == 0:
            self.add_edges(self._reverse_adj_list[node], [new_node])
            self.add_edge(new_node, node)  # TODO: Maybe include all nexts
            return new_node
        after_node = self._split_node(*position)
        self.add_edge(node, new_node)
        self.add_edge(new_node, after_node)

    def _split_node(self, node_id, split):
        """Split node before position"""
        # Update seq
        orig_seq = self._nodes[node_id]
        self._nodes[node_id] = orig_seq[:split]
        self._max_node += 1
        self._nodes[self._max_node] = orig_seq[split:]

        # Update edges
        self.add_edges([self._max_node], self._adj_list[node_id])
        self.set_edge(node_id, self._max_node)
        return self._max_node

    def _copy_node(self, node_id, seq):
        new_id = self.add_node(seq)
        self.add_edges(self._reverse_adj_list[node_id], [new_id])
        self.add_edges([new_id], self._adj_list[node_id])


class GraphMerger(GraphBuilder):
    def _get_offsets(self, interval):
        if len(interval.nodes) == 1:
            return [interval.end-interval.start]
        start_len = self.node_size(interval.nodes[0])-interval.start
        offsets = [start_len] + [self.node_size(n) for n in interval.nodes[1:-1]]+[interval.end]
        return offsets

    def _split_interval(self, interval, offsets):
        logging.debug("Split %s at %s", interval, offsets)
        nodes = iter(interval.nodes)
        cur_node = next(nodes)
        if interval.start > 0:
            cur_node = self._split_node(cur_node, interval.start)
        cur_node_size = self.node_size(cur_node)
        cur_offset = 0
        path = [cur_node]
        cur_node_size = self.node_size(cur_node)
        for i, offset in enumerate(offsets):
            cur_offset += offset
            while cur_offset > cur_node_size:
                cur_offset -= cur_node_size
                cur_node = next(nodes)
                path.append(cur_node)
                cur_node_size = self.node_size(cur_node)
            if 0 < cur_offset < cur_node_size:
                cur_node = self._split_node(cur_node, cur_offset)
                if i != len(offsets)-1:
                    path.append(cur_node)
                cur_offset = 0
                cur_node_size = self.node_size(cur_node)
        return path
        # TODO: return path with new nodes

    def _merge_nodes(self, node_a, node_b):
        logging.debug("Merging nodes %s and %s", node_a, node_b)
        seq_a = self._nodes[node_a]
        seq_b = self._nodes[node_b]
        assert self.node_size(node_a) == self.node_size(node_b)
        self.add_edges([node_a], self._adj_list[node_b])
        self.add_edges(self._reverse_adj_list[node_b], [node_a])
        self.remove_node(node_b)
        tmp = [(i, pair[1])
               for i, pair in enumerate(zip(seq_a, seq_b))
               if pair[0] != pair[1]]
        if not tmp:
            return
        offsets, seqs = zip(*tmp)
        self.add_snps(node_a, offsets, seqs)

    def _merge_paths(self, path_a, path_b):
        logging.info("Merging paths %s and %s", path_a, path_b)
        assert len(path_a) == len(path_b)
        [self._merge_nodes(*nodes) for nodes in zip(path_a, path_b)]

    def merge_intervals(self, interval_a, interval_b):
        offsets_a = self._get_offsets(interval_a)
        offsets_b = self._get_offsets(interval_b)
        new_path_a = self._split_interval(interval_a, offsets_b)
        new_path_b = self._split_interval(interval_b, offsets_a)
        logging.debug(new_path_a)
        logging.debug(new_path_b)
        self._merge_paths(new_path_a, new_path_b)


def merge_graphs(graph_a, graph_b, intervals_a, intervals_b):
    builder = GraphMerger()
    b_node_offset = builder.add_sequence_graph(graph_a)
    builder.add_sequence_graph(graph_b)
    intervals_b = [Interval(i.start, i.end, [node+b_node_offset for node in i.nodes])
                   for i in intervals_b]
    for interval_pair in reversed(zip(intervals_a, intervals_b)):
        builder.merge_intervals(*interval_pair)

    return builder.to_sequence_graph()
