"""Characterization tests for acore.network_analysis: what each function does today."""

import networkx as nx
import pytest

pytest.importorskip("community")
pytest.importorskip("snf")

import acore.network_analysis as na


def two_cliques() -> nx.Graph:
    """Nodes 0-3 and 4-7 each fully connected, joined by the single edge 3-4."""
    graph = nx.Graph()
    graph.add_edges_from(
        (a, b)
        for block in (range(4), range(4, 8))
        for a in block
        for b in block
        if a < b
    )
    graph.add_edge(3, 4)
    nx.set_edge_attributes(graph, 1.0, "weight")
    nx.set_edge_attributes(graph, 1.0, "width")
    return graph


def assert_cliques_split(partition: dict) -> None:
    assert set(partition) == set(range(8))
    assert len({partition[node] for node in range(4)}) == 1
    assert len({partition[node] for node in range(4, 8)}) == 1
    assert partition[0] != partition[4]


def test_louvain_partitions_split_two_cliques_joined_by_one_edge():
    # Given two cliques joined by a single edge
    graph = two_cliques()

    # When the Louvain partition is computed on the edge weights
    partition = na.get_louvain_partitions(graph, "weight")

    # Then each clique is one community, numbered from 0
    assert_cliques_split(partition)
    assert set(partition.values()) == {0, 1}


@pytest.mark.parametrize(
    "algorithm", ["louvain", "greedy_modularity", "affinity_propagation"]
)
def test_network_communities_split_two_cliques_with_each_algorithm(algorithm):
    # Given two cliques joined by a single edge
    graph = two_cliques()

    # When communities are found with one of the algorithms
    communities = na.get_network_communities(
        graph, {"values": "weight", "communities_algorithm": algorithm}
    )

    # Then each clique is one community
    assert_cliques_split(communities)


def test_network_communities_default_to_louvain_and_record_the_choice():
    # Given settings that name no algorithm
    graph = two_cliques()
    settings = {"values": "weight"}

    # When communities are found
    communities = na.get_network_communities(graph, settings)

    # Then Louvain ran, and the settings dictionary now says so
    assert_cliques_split(communities)
    assert settings["communities_algorithm"] == "louvain"


def test_most_central_edge_returns_the_node_with_the_highest_centrality():
    # Given a star whose centre is the most central node
    graph = nx.star_graph(4)
    nx.set_edge_attributes(graph, 1.0, "width")

    # When the most central edge is asked for
    result = na.most_central_edge(graph)

    # Then it is the centre node, not a centrality value
    assert result == 0
