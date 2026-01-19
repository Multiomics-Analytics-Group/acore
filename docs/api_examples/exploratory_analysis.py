# %% [markdown]
# # Exploratory Analysis

# %%
# %pip install acore

# %%
from typing import Optional

import matplotlib.pyplot as plt
import pandas as pd

import acore.exploratory_analysis as ea
from acore.types.exploratory_analysis import (
    AnnotationResult,
    TwoComponentSchema,
    TwoLoadingsSchema,
    TwoVariance,
)


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


data = (
    "https://raw.githubusercontent.com/Multiomics-Analytics-Group/acore/"
    "refs/heads/main/"
    "example_data/MTBLS13311/MTBLS13411_processed_data.csv"
)
data = pd.read_csv(data, index_col=0)
data

# %% [markdown]
# We add the group here based on the sample names. Alternatively you could merge it from
# the avilable metadata.

# %%
data["group"] = data.index.str.split("-").str[0]
data["group"].value_counts()

# %% [markdown]
# Show first two principal components of the data.

# %%
map_names = {
    "value": "feature_communiality",
    "x": "PC1",
    "y": "PC2",
}
results_dfs, annotation = ea.run_pca(
    data, drop_cols=[], annotation_cols=[], group="group", components=2, dropna=True
)
pcs, loadings, var_explained = results_dfs

# %%
TwoVariance(pd.Series(var_explained, index=["PC1", "PC2"]))

# %%
annotation = AnnotationResult(**annotation)  # .model_dump()
annotation

# %%
fig, ax = make_plot(pcs, annotation=annotation, **map_names)

# %% [markdown]
# Show what was computed:

# %%
TwoComponentSchema(pcs).rename(columns=map_names)

# %%
TwoLoadingsSchema(loadings).rename(columns=map_names)

# %% [markdown]
# Visualize UMAP low-dimensional embedding of the data.

# %%
map_names = {
    "x": "UMAP1",
    "y": "UMAP2",
}
result, annotation = ea.run_umap(
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
fig, ax = make_plot(result["umap"], annotation=annotation, **map_names)
TwoComponentSchema(result["umap"]).rename(columns=map_names)

# %% [markdown]
# Make sure to check the parameter annotations in the API docs.


# %% [markdown]
# ## Correlation analysis
#
# ### Coefficient of variation
# - as of now does an internal data transformation

# %%
res = ea.get_coefficient_variation(data=data, group="group")
res

# %%
map_names = {"x": "mean", "y": "coef_of_var", "group": "group"}
fig, ax = make_plot(res, **map_names)
