import community
import networkx as nx
import pandas as pd
import snf
from sklearn import cluster
from sklearn.cluster import AffinityPropagation

from acore import utils


def get_network_communities(graph, args):
    """
    Finds communities in a graph using different methods. For more information on the methods visit:

        - https://networkx.github.io/documentation/latest/reference/algorithms/generated/networkx.algorithms.community.modularity_max.greedy_modularity_communities.html
        - https://networkx.github.io/documentation/networkx-2.0/reference/algorithms/generated/networkx.algorithms.community.asyn_lpa.asyn_lpa_communities.html
        - https://networkx.github.io/documentation/latest/reference/algorithms/generated/networkx.algorithms.community.centrality.girvan_newman.html
        - https://networkx.github.io/documentation/latest/reference/generated/networkx.convert_matrix.to_pandas_adjacency.html

    :param graph graph: networkx graph
    :param dict args: config file arguments
    :return: Dictionary of nodes and which community they belong to (from 0 to number of communities).
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


def get_snf_clusters(data_tuples, num_clusters=None, metric="euclidean", k=5, mu=0.5):
    """
    Cluster samples based on Similarity Network Fusion (SNF) (ref: https://www.ncbi.nlm.nih.gov/pubmed/24464287)

    :param df_tuples: list of (dataset,metric) tuples
    :param index: how the datasets can be merged (common columns)
    :param num_clusters: number of clusters to be identified, if None, the algorithm finds the best number based on SNF algorithm (recommended)
    :param distance_metric: distance metric used to calculate the sample similarity network
    :param k: number of neighbors used to measure local affinity (KNN)
    :param mu: normalization factor to scale similarity kernel when constructing affinity matrix
    :return tuple: 1) fused_aff: affinity clustered samples, 2) fused_labels: cluster labels,
                    3) num_clusters: number of clusters, 4) silhouette: average silhouette score
    """
    affinities = []
    for d, m in data_tuples:
        affinities += [snf.make_affinity(d, metric=m, K=k, mu=mu)]
    fused_aff = snf.snf(affinities, K=k)
    if num_clusters is None:
        num_clusters, second = snf.get_n_clusters(fused_aff)
    fused_labels = cluster.spectral_clustering(fused_aff, n_clusters=num_clusters)
    fused_labels = [i + 1 for i in fused_labels]
    silhouette = snf.metrics.silhouette_score(fused_aff, fused_labels)

    return (fused_aff, fused_labels, num_clusters, silhouette)


def most_central_edge(G):
    """
    Compute the eigenvector centrality for the graph G, and finds the highest value.

    :param graph G: networkx graph
    :return: Highest eigenvector centrality value.
    :rtype: float
    """
    centrality = nx.eigenvector_centrality_numpy(G, weight="width")

    return max(centrality, key=centrality.get)


def get_louvain_partitions(G, weight):
    """
    Computes the partition of the graph nodes which maximises the modularity (or try..) using the Louvain heuristices. For more information visit https://python-louvain.readthedocs.io/en/latest/api.html.

    :param graph G: networkx graph which is decomposed.
    :param str weight: the key in graph to use as weight.
    :return: The partition, with communities numbered from 0 to number of communities.
    :rtype: dict
    """
    partition = community.best_partition(G, weight=weight)

    return partition


def run_snf(
    df_dict,
    index,
    num_clusters=None,
    distance_metric="euclidean",
    k_affinity=5,
    mu_affinity=0.5,
):
    """
    Runs Similarity Network Fusion: integration of multiple omics datasets to identify
    similar samples (clusters) (ref: https://www.ncbi.nlm.nih.gov/pubmed/24464287).
    We make use of the pyton version SNFpy (https://github.com/rmarkello/snfpy)

    :param df_dict: dictionary of datasets to be used (i.e {'rnaseq': rnaseq_data, 'proteomics': proteomics_data})
    :param index: how the datasets can be merged (common columns)
    :param num_clusters: number of clusters to be identified, if None, the algorithm finds the best number based on SNF algorithm (recommended)
    :param distance_metric: distance metric used to calculate the sample similarity network
    :param k_affinity: number of neighbors used to measure local affinity (KNN)
    :param mu_ffinity: normalization factor to scale similarity kernel when constructing affinity matrix
    :return tuple: 1) feature_df: SNF features and mutual information score (MIscore), 2) fused_aff: adjacent similarity matrix, 3)fused_labels: cluster labels,
                    4) silhouette: silhouette score
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
    i = 0
    for dtype in snf_features:
        df = pd.DataFrame(dtype, index=indexes[i], columns=["MIscore"]).sort_values(
            by="MIscore", ascending=False
        )
        df["dataset"] = dataset_labels[i]
        i += 1
        feature_df = feature_df.append(df)

    feature_df = feature_df.sort_values(by="MIscore", ascending=False)

    return feature_df, fused_aff, fused_labels, silhouette_score
