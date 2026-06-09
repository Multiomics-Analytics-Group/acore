# acore.network_analysis package

### get_network_communities(graph, args)

Finds communities in a graph using different methods. For more information on the methods visit:

> - [https://networkx.github.io/documentation/latest/reference/algorithms/generated/networkx.algorithms.community.modularity_max.greedy_modularity_communities.html](https://networkx.github.io/documentation/latest/reference/algorithms/generated/networkx.algorithms.community.modularity_max.greedy_modularity_communities.html)
> - [https://networkx.github.io/documentation/networkx-2.0/reference/algorithms/generated/networkx.algorithms.community.asyn_lpa.asyn_lpa_communities.html](https://networkx.github.io/documentation/networkx-2.0/reference/algorithms/generated/networkx.algorithms.community.asyn_lpa.asyn_lpa_communities.html)
> - [https://networkx.github.io/documentation/latest/reference/algorithms/generated/networkx.algorithms.community.centrality.girvan_newman.html](https://networkx.github.io/documentation/latest/reference/algorithms/generated/networkx.algorithms.community.centrality.girvan_newman.html)
> - [https://networkx.github.io/documentation/latest/reference/generated/networkx.convert_matrix.to_pandas_adjacency.html](https://networkx.github.io/documentation/latest/reference/generated/networkx.convert_matrix.to_pandas_adjacency.html)
* **Parameters:**
  * **graph** (*graph*) – networkx graph
  * **args** ([*dict*](https://docs.python.org/3/library/stdtypes.html#dict)) – config file arguments
* **Returns:**
  Dictionary of nodes and which community they belong to (from 0 to number of communities).

### get_snf_clusters(data_tuples, num_clusters=None, metric='euclidean', k=5, mu=0.5)

Cluster samples based on Similarity Network Fusion (SNF) (ref: [https://www.ncbi.nlm.nih.gov/pubmed/24464287](https://www.ncbi.nlm.nih.gov/pubmed/24464287))

* **Parameters:**
  * **df_tuples** – list of (dataset,metric) tuples
  * **index** – how the datasets can be merged (common columns)
  * **num_clusters** – number of clusters to be identified, if None, the algorithm finds the best number based on SNF algorithm (recommended)
  * **distance_metric** – distance metric used to calculate the sample similarity network
  * **k** – number of neighbors used to measure local affinity (KNN)
  * **mu** – normalization factor to scale similarity kernel when constructing affinity matrix
* **Return tuple:**
  1) fused_aff: affinity clustered samples, 2) fused_labels: cluster labels,
  3) num_clusters: number of clusters, 4) silhouette: average silhouette score

### most_central_edge(G)

Compute the eigenvector centrality for the graph G, and finds the highest value.

* **Parameters:**
  **G** (*graph*) – networkx graph
* **Returns:**
  Highest eigenvector centrality value.
* **Return type:**
  [float](https://docs.python.org/3/library/functions.html#float)

### get_louvain_partitions(G, weight)

Computes the partition of the graph nodes which maximises the modularity (or try..) using the Louvain heuristices. For more information visit [https://python-louvain.readthedocs.io/en/latest/api.html](https://python-louvain.readthedocs.io/en/latest/api.html).

* **Parameters:**
  * **G** (*graph*) – networkx graph which is decomposed.
  * **weight** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – the key in graph to use as weight.
* **Returns:**
  The partition, with communities numbered from 0 to number of communities.
* **Return type:**
  [dict](https://docs.python.org/3/library/stdtypes.html#dict)

### run_snf(df_dict, index, num_clusters=None, distance_metric='euclidean', k_affinity=5, mu_affinity=0.5)

Runs Similarity Network Fusion: integration of multiple omics datasets to identify
similar samples (clusters) (ref: [https://www.ncbi.nlm.nih.gov/pubmed/24464287](https://www.ncbi.nlm.nih.gov/pubmed/24464287)).
We make use of the pyton version SNFpy ([https://github.com/rmarkello/snfpy](https://github.com/rmarkello/snfpy))

* **Parameters:**
  * **df_dict** – dictionary of datasets to be used (i.e {‘rnaseq’: rnaseq_data, ‘proteomics’: proteomics_data})
  * **index** – how the datasets can be merged (common columns)
  * **num_clusters** – number of clusters to be identified, if None, the algorithm finds the best number based on SNF algorithm (recommended)
  * **distance_metric** – distance metric used to calculate the sample similarity network
  * **k_affinity** – number of neighbors used to measure local affinity (KNN)
  * **mu_ffinity** – normalization factor to scale similarity kernel when constructing affinity matrix
* **Return tuple:**
  1) feature_df: SNF features and mutual information score (MIscore), 2) fused_aff: adjacent similarity matrix, 3)fused_labels: cluster labels,
  4) silhouette: silhouette score
