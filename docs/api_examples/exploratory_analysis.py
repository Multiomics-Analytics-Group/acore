# %% [markdown]
# # Exploratory Analysis

# %%
# %pip install acore

# %%
import pandas as pd

import acore.exploratory_analysis as ea

data = pd.DataFrame(
    {
        "group": ["A", "A", "B", "B"],
        "protein1": [1.4, 2.2, 5.3, 4.2],
        "protein2": [5.6, 0.3, 2.1, 8.1],
        "protein3": [9.1, 10.01, 11.2, 12.9],
    }
)

# %% [markdown]
# Show first two principal components of the data.

# %%
result_dfs, annotation = ea.run_pca(
    data, drop_cols=[], annotation_cols=[], group="group", components=2, dropna=True
)

# %% [markdown]
# Show what was computed:

# %%
result_dfs[0]

# %%
result_dfs[1]

# %%
result_dfs[2]

# %%
annotation

# %% [markdown]
# Visualize UMAP low-dimensional embedding of the data.

# %%
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
result["umap"]

# %%
annotation

# %% [markdown]
# Make sure to check the parameter annotations in the API docs.
