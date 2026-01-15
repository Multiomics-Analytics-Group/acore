# %% [markdown]
# # Exploratory Analysis

# %%
# %pip install acore

# %%
import matplotlib.pyplot as plt
import pandas as pd

import acore.exploratory_analysis as ea
from acore.types.exploratory_analysis import (
    AnnotationResult,
    TwoComponentSchema,
    TwoLoadingsSchema,
    TwoVariance,
)

data = pd.DataFrame(
    {
        "group": ["A", "A", "B", "B"],
        "protein1": [1.4, 2.2, 5.3, 4.2],
        "protein2": [5.6, 0.3, 2.1, 8.1],
        "protein3": [9.1, 10.01, 11.2, 12.9],
    }
)

data = (
    "https://raw.githubusercontent.com/Multiomics-Analytics-Group/acore/"
    "refs/heads/add_metabolomics_data/"
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
fig, ax = plt.subplots()
for i, (group, group_df) in enumerate(pcs.groupby("group")):
    ax = group_df.rename(columns=map_names).plot.scatter(
        x="PC1",
        y="PC2",
        label=group,
        c=f"C{i}",
        ax=ax,
    )
_ = ax.set(ylabel=annotation.y_title, xlabel=annotation.x_title)

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
TwoComponentSchema(result["umap"]).rename(columns=map_names)

# %%
AnnotationResult(**annotation)

# %% [markdown]
# Make sure to check the parameter annotations in the API docs.
