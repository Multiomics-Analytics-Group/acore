# %% [markdown]
# # Decomposition
#
# Using a metabolomics dataset ([MTBLS13411](../../../example_data/MTBLS13411.md)) as an
# example, we show how to use the decomposition methods
# for the four groups present in the data. The groups are derived based on the sample
# names, but you could also merge them from available metadata.
#
# - Principal Component Analysis (PCA)
# - Uniform Manifold Approximation and Projection (UMAP)
# - t-Distributed Stochastic Neighbor Embedding (t-SNE)
#
# The focus of acore is not plotting, but getting the data in the right format for
# plotting. Nonetheless, we show some examples of how to plot the results. VueCore
# offers more examples.

# %%
# %pip install acore

# %%
from pprint import pprint
from typing import Optional

import matplotlib.pyplot as plt
import pandas as pd

import acore.exploratory_analysis as ea
from acore.types.exploratory_analysis import (
    AnnotationResult,
    ThreeComponentSchema,
    ThreeLoadingsSchema,
    TwoComponentSchema,
    TwoLoadingsSchema,
    TwoVariance,
    VarianceRatio,
)

# %% [markdown]
# Utility function for plotting


# %%
# relace wit vuecore plot
def make_plot(
    embeddings,
    x: str,
    y: str,
    annotation: Optional[dict[str, str]] = None,
    group: str = "group",
    **kwargs,
):
    """Utility function for static plot of dimensionality reductions."""
    fig, ax = plt.subplots()
    for i, (group, group_df) in enumerate(embeddings.groupby("group")):
        ax = group_df.rename(columns=map_names).plot.scatter(
            x=x,
            y=y,
            label=group,
            c=f"C{i}",
            ax=ax,
        )
    if annotation is not None:
        _ = ax.set(ylabel=annotation.y_title, xlabel=annotation.x_title)
    return fig, ax


# %% [markdown]
# ## Load metabolomics example data

# %%
data = (
    "https://raw.githubusercontent.com/Multiomics-Analytics-Group/acore/"
    "refs/heads/main/"
    "example_data/MTBLS13311/MTBLS13411_processed_data.csv"
)
data = pd.read_csv(data, index_col=0)
# specific to this data, we shorten some column names for better readability
data.columns = data.columns.str.split("(").str[-1].str.replace(")", "")
data

# %% [markdown]
# We add the group here based on the sample names. Alternatively you could merge it from
# the avilable metadata.

# %%
data["group"] = data.index.str.split("-").str[0]
data["group"].value_counts()

# %% [markdown]
# ## Principal Component Analysis (PCA)
# Show first two principal components of the data.

# %%
# map_names gives the column names for the plot axes (which default to "x" and "y")
map_names = {
    "value": "feature_communiality",
    "x": "PC1",
    "y": "PC2",
}
results_dfs, annotation = ea.decomposition.run_pca(
    data, drop_cols=[], annotation_cols=[], group="group", components=2, dropna=True
)
pcs, loadings, var_explained = results_dfs

# %% [markdown]
# See how much variance is explained by the first two components and validate that
# they adhere to the expected format:

# %%
pprint(var_explained)
TwoVariance(pd.Series(var_explained, index=["PC1", "PC2"]))

# %% [markdown]
# Show the annotation information for plotting and validate that
# they adhere to the expected format:

# %%
pprint(annotation)
annotation = AnnotationResult(**annotation)
annotation

# %% [markdown]
# Show what was computed and validate that
# they adhere to the expected format:
# - first two principal components of the samples
# - loadings for the features on the first two components
#
# We rename the columns for better readability.

# %%
TwoComponentSchema(pcs).rename(columns=map_names).head(10)

# %% [markdown]
# The feature communality of the loading is the absolute length of the projection.
# So the features listed first here contribute the most to the two first components,
# therefore driving the PCA separation.

# %%
TwoLoadingsSchema(loadings).rename(columns=map_names).head(10)

# %% [markdown]
# Show the PCA plot of the first two components. VueCore understands the
# annotation information (ToDo)

# %%
# replace wit vuecore plot
fig, ax = make_plot(pcs, annotation=annotation, **map_names)

# %% [markdown]
# ### More than two components
# With more than two components, the results contain an additional `z` column (and
# `3`, `4`, ... for further components). The annotation then also has a `z_title`.
# The `ThreeComponentSchema` and `ThreeLoadingsSchema` validate these results, and
# `VarianceRatio` checks that the explained variances are fractions between 0 and 1.

# %%
results_dfs_3, annotation_3 = ea.decomposition.run_pca(
    data, drop_cols=[], annotation_cols=[], group="group", components=3, dropna=True
)
pcs_3, loadings_3, var_explained_3 = results_dfs_3
VarianceRatio(pd.Series(var_explained_3, index=["PC1", "PC2", "PC3"]))
AnnotationResult(**annotation_3)

# %%
ThreeComponentSchema(pcs_3)

# %%
ThreeLoadingsSchema(loadings_3)

# %% [markdown]
# ## Uniform Manifold Approximation and Projection (UMAP)
# Visualize UMAP low-dimensional embedding of the data.
# This uses the `umap-learn` package, which is documented with examples at
# [umap-learn.readthedocs.io](https://umap-learn.readthedocs.io).

# %%
# map_names gives the column names for the plot axes (which default to "x" and "y")
map_names = {
    "x": "UMAP1",
    "y": "UMAP2",
}
result, annotation = ea.decomposition.run_umap(
    data,
    drop_cols=["sample", "subject"],
    group="group",
    n_neighbors=10,
    min_dist=0.3,
    metric="cosine",
    dropna=True,
)

# %%
annotation = AnnotationResult(**annotation)
annotation

# %%
fig, ax = make_plot(result["umap"], **map_names)
TwoComponentSchema(result["umap"]).rename(columns=map_names)

# %% [markdown]
# Make sure to check the parameters and tutorials annotations in the API docs at
# [umap-learn.readthedocs.io](https://umap-learn.readthedocs.io).


# %% [markdown]
# ## t-Distributed Stochastic Neighbor Embedding (t-SNE)
#
# Visualize the t-SNE low-dimensional embedding of the data, using
# [`sklearn.manifold.TSNE`](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html).
# The default `barnes_hut` method can differ slightly between CPU architectures, so we
# use the `exact` method and a `random_state` for reproducible results here.

# %%
map_names = {
    "x": "t-SNE1",
    "y": "t-SNE2",
}
result, annotation = ea.decomposition.run_tsne(
    data,
    drop_cols=["sample", "subject"],
    group="group",
    perplexity=10,
    method="exact",
    random_state=42,
    dropna=True,
)

# %%
annotation = AnnotationResult(**annotation)
annotation

# %%
fig, ax = make_plot(result["tsne"], **map_names)
TwoComponentSchema(result["tsne"]).rename(columns=map_names)

# %% [markdown]
# Done.
