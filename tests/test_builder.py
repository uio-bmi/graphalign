import pytest

from graphalign.builder import GraphBuilder, BuilderStruct, GraphMerger
from graphalign import *


@pytest.fixture
def short_ref():
    return [0, 1, 2, 3, 3, 2, 1, 0]


@pytest.fixture
def short_ref2():
    return [3, 3, 2, 3, 3, 2, 3, 3]


@pytest.fixture
def merge_input():
    return GraphMerger({0: short_ref(), 1: short_ref2()})


@pytest.fixture
def short_graph():
    return GraphBuilder({0: short_ref()})


def test_add_snp(short_graph, short_ref):
    short_graph.add_snp((0, 3), 2)
    nodes = {0: short_ref[:3], 1: [3], 2: short_ref[4:], 3: [2]}
    adj_list = {0: [1, 3], 1: [2], 3: [2]}
    assert short_graph.to_struct() == BuilderStruct(nodes, adj_list)


def test_add_del(short_graph, short_ref):
    short_graph.add_deletion(Interval(2, 5, [0]))
    nodes = {0: short_ref[:2], 2: short_ref[2:5], 1: short_ref[5:]}
    adj_list = {0: [1, 2], 2: [1]}
    assert short_graph.to_struct() == BuilderStruct(nodes, adj_list)


def test_add_insertion(short_graph, short_ref):
    short_graph.add_insertion(Position(0, 4), [1, 1, 1])
    nodes = {0: short_ref[:4], 2: short_ref[4:], 1: [1, 1, 1]}
    adj_list = {0: [1, 2], 1: [2]}
    assert short_graph.to_struct() == BuilderStruct(nodes, adj_list)


def test_simple_merge(merge_input):
    merge_input.merge_intervals(Interval(2, 6, [0]),
                                Interval(2, 6, [1]))
    nodes = {0: [0, 1],
             1: [3, 3],
             2: [2, 3, 3, 2],
             3: [1, 0],
             5: [3, 3]}
    adj_list = {0: [2], 1: [2], 2: [3, 5]}
    assert merge_input.to_struct() == BuilderStruct(nodes, adj_list)


def test_divided_merge(short_ref, short_ref2):
    adj_list = {i: [i+1] for i in range(4)}
    del adj_list[1]
    builder = GraphMerger({0: short_ref[:4], 1: short_ref[4:],
                           2: short_ref2[:3], 3: short_ref2[3:6],
                           4: short_ref2[6:8]}, adj_list)
    
    builder.merge_intervals(Interval(2, 2, [0, 1]),
                            Interval(2, 3, [2, 3]))
    nodes = {0: short_ref[0:2],
             2: short_ref2[0:2],
             5: short_ref[2:3],
             6: short_ref[3:4],
             1: short_ref[4:6],
             7: short_ref[6:],
             4: short_ref2[6:]}
    adj_list = {0: [5], 2: [5], 5: [6], 6: [1], 1: [4, 7]}
    assert builder.to_struct() == BuilderStruct(nodes, adj_list)
