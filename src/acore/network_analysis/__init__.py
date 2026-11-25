try:
    import community
    import networkx as nx
    import snf
except ImportError as e:
    raise ImportError(
        "Error importing network analysis dependencies (community/python-louvain, networkx, snf/snfpy). "
        "Install them with: pip install 'acore[all]'"
        f"\n\nError: {e}"
    ) from e

from typing import Any

import numpy as np
import pandas as pd
from sklearn import cluster
from sklearn.cluster import AffinityPropagation

from acore import utils


def get_network_communities(graph: nx.Graph, args: dict[str, Any]) -> dict[Any, int]:
    """
    Finds communities in a graph using different methods. For more information on the methods visit:

        - https://networkx.github.io/documentation/latest/reference/algorithms/generated/networkx.algorithms.community.modularity_max.greedy_modularity_communities.html
        - https://networkx.github.io/documentation/networkx-2.0/reference/algorithms/generated/networkx.algorithms.community.asyn_lpa.asyn_lpa_communities.html
        - https://networkx.github.io/documentation/latest/reference/algorithms/generated/networkx.algorithms.community.centrality.girvan_newman.html
        - https://networkx.github.io/documentation/latest/reference/generated/networkx.convert_matrix.to_pandas_adjacency.html

    :param nx.Graph graph: networkx graph whose nodes are grouped into communities.
    :param dict args: settings, read by key. ``values`` is the name of the edge
                      attribute that holds the edge weight; it is required by every
                      algorithm except ``girvan_newman`` and ``affinity_propagation``
                      (which reads the edge attribute ``width`` instead).
                      ``communities_algorithm`` picks the method and is one of
                      ``'louvain'``, ``'greedy_modularity'``,
                      ``'asyn_label_propagation'``, ``'girvan_newman'`` and
                      ``'affinity_propagation'``; it defaults to ``'louvain'`` and is
                      added to ``args`` when absent.
    :return: Dictionary mapping each node to the number of the community it belongs
             to (from 0 to the number of communities minus one).
    :rtype: dict
    """
    if "communities_algorithm" not in args:
        args["communities_algorithm"] = "louvain"

    if args["communities_algorithm"] == "louvain":
        communities = get_louvain_partitions(graph, args["values"])
    elif args["communities_algorithm"] == "greedy_modularity":
        gcommunities = nx.algorithms.community.greedy_modularity_communities(
            graph, weight=args["values"]
        )
        communities = utils.generator_to_dict(gcommunities)
    elif args["communities_algorithm"] == "asyn_label_propagation":
        gcommunities = nx.algorithms.community.label_propagation.asyn_lpa_communities(
            graph, args["values"]
        )
        communities = utils.generator_to_dict(gcommunities)
    elif args["communities_algorithm"] == "girvan_newman":
        gcommunities = nx.algorithms.community.girvan_newman(
            graph, most_valuable_edge=most_central_edge
        )
        communities = utils.generator_to_dict(gcommunities)
    elif args["communities_algorithm"] == "affinity_propagation":
        adjacency = nx.to_pandas_adjacency(graph, weight="width")
        nodes = list(adjacency.columns)
        communities = AffinityPropagation().fit(adjacency.values).labels_
        communities = {nodes[i]: communities[i] for i in range(len(communities))}

    return communities


def get_snf_clusters(
    data_tuples: list[tuple[np.ndarray, str]],
    num_clusters: int | None = None,
    metric: str = "euclidean",
    k: int = 5,
    mu: float = 0.5,
) -> tuple[np.ndarray, list[int], int, float]:
    """
    Cluster samples based on Similarity Network Fusion (SNF) (ref: https://www.ncbi.nlm.nih.gov/pubmed/24464287)

    :param list data_tuples: one ``(data, metric)`` tuple per dataset, where ``data``
                              is a samples-by-features array and ``metric`` the
                              distance metric used to compare its samples. Every
                              dataset must hold the same samples in the same order.
    :param int num_clusters: number of clusters to be identified, if None, the algorithm finds the best number based on SNF algorithm (recommended)
    :param str metric: distance metric used to calculate the sample similarity network
    :param int k: number of neighbors used to measure local affinity (KNN)
    :param float mu: normalization factor to scale similarity kernel when constructing affinity matrix
    :return: Tuple of 1) fused_aff: the fused sample-by-sample affinity matrix,
             2) fused_labels: cluster labels, numbered from 1, in sample order,
             3) num_clusters: number of clusters, 4) silhouette: average silhouette score
    :rtype: tuple
    """
    affinities = []
    for d, m in data_tuples:
        affinities += [snf.make_affinity(d, metric=m, K=k, mu=mu)]
    fused_aff = snf.snf(affinities, K=k)
    if num_clusters is None:
        num_clusters, _second = snf.get_n_clusters(fused_aff)
    fused_labels = cluster.spectral_clustering(fused_aff, n_clusters=num_clusters)
    fused_labels = [i + 1 for i in fused_labels]
    silhouette = snf.metrics.silhouette_score(fused_aff, fused_labels)

    return (fused_aff, fused_labels, num_clusters, silhouette)


def most_central_edge(G: nx.Graph) -> Any:
    """
    Compute the eigenvector centrality for the graph G, and finds the highest value.

    This is the ``most_valuable_edge`` callable that ``girvan_newman`` takes in
    :func:`get_network_communities`. It ranks nodes, using the edge attribute ``width``
    as weight, and returns the node with the highest centrality, not the centrality value.

    :param nx.Graph G: networkx graph whose edges have a ``width`` attribute.
    :return: The node with the highest eigenvector centrality.
    :rtype: Any
    """
    centrality = nx.eigenvector_centrality_numpy(G, weight="width")

    return max(centrality, key=centrality.get)


def get_louvain_partitions(G: nx.Graph, weight: str) -> dict[Any, int]:
    """
    Computes the partition of the graph nodes which maximises the modularity (or try..) using the Louvain heuristices. For more information visit https://python-louvain.readthedocs.io/en/latest/api.html.

    :param nx.Graph G: networkx graph which is decomposed.
    :param str weight: name of the edge attribute to use as weight.
    :return: The partition: a dictionary mapping each node to its community, with
             communities numbered from 0 to the number of communities minus one.
    :rtype: dict
    """
    partition = community.best_partition(G, weight=weight)

    return partition


def run_snf(
    df_dict: dict[str, pd.DataFrame],
    index: str,
    num_clusters: int | None = None,
    distance_metric: str = "euclidean",
    k_affinity: int = 5,
    mu_affinity: float = 0.5,
) -> tuple[pd.DataFrame, np.ndarray, pd.DataFrame, float]:
    """
    Runs Similarity Network Fusion: integration of multiple omics datasets to identify
    similar samples (clusters) (ref: https://www.ncbi.nlm.nih.gov/pubmed/24464287).
    We make use of the pyton version SNFpy (https://github.com/rmarkello/snfpy)

    :param dict df_dict: dictionary of datasets to be used (i.e {'rnaseq': rnaseq_data, 'proteomics': proteomics_data}).
                         Each dataset has one row per sample, the sample identifier
                         in the column ``index`` and the other columns as features.
    :param str index: name of the column that identifies a sample in every dataset; only the
                      samples present in all datasets are used.
    :param int num_clusters: number of clusters to be identified, if None, the algorithm finds the best number based on SNF algorithm (recommended)
    :param str distance_metric: distance metric used to calculate the sample similarity network
    :param int k_affinity: number of neighbors used to measure local affinity (KNN)
    :param float mu_affinity: normalization factor to scale similarity kernel when constructing affinity matrix
    :return: Tuple of 1) feature_df: SNF features and mutual information score (MIscore),
             2) fused_aff: the fused sample-by-sample similarity matrix,
             3) fused_labels: a dataframe with one row per sample and the cluster in
             the column ``cluster``, 4) silhouette: silhouette score
    :rtype: tuple
    """
    datasets = []
    dataset_labels = []
    common_samples = set()
    for dtype in df_dict:
        dataset_labels.append(dtype)
        df = df_dict[dtype]
        df = df.set_index(index)
        datasets.append(df)
        if len(common_samples) > 1:
            common_samples = common_samples.intersection(df.index)
        else:
            common_samples.update(df.index.tolist())

    data_tuples = [
        (d.loc[list(common_samples)].values, distance_metric) for d in datasets
    ]

    fused_aff, fused_labels, num_clusters, silhouette_score = get_snf_clusters(
        data_tuples, num_clusters, metric=distance_metric, k=k_affinity, mu=mu_affinity
    )

    fused_labels = pd.DataFrame(fused_labels, index=common_samples, columns=["cluster"])

    snf_features = snf.metrics.rank_feature_by_nmi(
        data_tuples, fused_aff, K=k_affinity, mu=mu_affinity, n_clusters=num_clusters
    )

    feature_df = pd.DataFrame(columns=["MIscore"])
    indexes = [df.columns for df in datasets]
    for i, dtype in enumerate(snf_features):
        df = pd.DataFrame(dtype, index=indexes[i], columns=["MIscore"]).sort_values(
            by="MIscore", ascending=False
        )
        df["dataset"] = dataset_labels[i]
        feature_df = feature_df.append(df)

    feature_df = feature_df.sort_values(by="MIscore", ascending=False)

    return feature_df, fused_aff, fused_labels, silhouette_score
