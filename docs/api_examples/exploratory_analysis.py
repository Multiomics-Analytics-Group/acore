# %% [markdown]
# # Exploratory Analysis

# %%
# %pip install acore

# %%
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import seaborn as sns

import acore.correlation_analysis as ca
import acore.exploratory_analysis as ea
from acore.types.exploratory_analysis import (
    AnnotationResult,
    TwoComponentSchema,
    TwoLoadingsSchema,
    TwoVariance,
)

# %% [markdown]
# Utility function for plotting


# %%
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
results_dfs, annotation = ea.run_pca(
    data, drop_cols=[], annotation_cols=[], group="group", components=2, dropna=True
)
pcs, loadings, var_explained = results_dfs

# %% [markdown]
# See how much variance is explained by the first two components and validate that
# they adhere to the expected format:

# %%
TwoVariance(pd.Series(var_explained, index=["PC1", "PC2"]))

# %% [markdown]
# Show the annotation information for plotting and validate that
# they adhere to the expected format:

# %%
annotation = AnnotationResult(**annotation)
annotation

# %% [markdown]
# Make the PCA plot:

# %%
fig, ax = make_plot(pcs, annotation=annotation, **map_names)

# %% [markdown]
# Show what was computed and validate that
# they adhere to the expected format:
# - first two principal components of the samples
# - loadings for the features on the first two components
#
# We rename the columns for better readability.

# %%
TwoComponentSchema(pcs).rename(columns=map_names)

# %% [markdown]
# The feature communality of the loading is the absolute length of the projection.
# So the features listed first here contribute the most to the two first components,
# therefore driving the PCA separation.

# %%
TwoLoadingsSchema(loadings).rename(columns=map_names)

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
# Make sure to check the parameters and tutorials annotations in the API docs at
# [umap-learn.readthedocs.io](https://umap-learn.readthedocs.io).
#
#
# ## Coefficient of variation
# Using masspectrometry data, we can compute the coefficient of variation on the
# non-log transformed intensities. We do this for each group separately.
# First we undo the log transformation, which is something specific to this dataset.

# %%
data_exp = data.drop(columns=["group"]).apply(lambda x: np.exp2(x)).join(data["group"])
data_exp

# %%
res = ea.get_coefficient_variation(data=data_exp, group="group")
res

# %%
res.describe()

# %%
map_names = {"x": "mean_log2", "y": "coef_of_var", "group": "group"}
fig, ax = make_plot(res, **map_names)

# %% [markdown]
# ## Correlation analysis
# See [`acore.correlation_analysis`](acore.correlation_analysis) for more functions
# and details.
#
# The basic functionality is built into pandas, but you need to filter out columns
# which are not numeric for pearson correlation.
#
# Generally: Ordered categorical values can be used, assuming equal spacing between
# the categories. Otherwise, continous numeric values are required.

# %%
corr = data.drop(columns=["group"]).corr(method="pearson")
corr

# %% [markdown]
# Plot the correlation heatmap using seaborn

# %%
plt.rcParams["xtick.labelsize"], plt.rcParams["ytick.labelsize"] = 5, 5
fig, ax = plt.subplots(figsize=(7.1, 6))
heatmap = sns.heatmap(
    corr,
    cmap="vlag",
    center=0,
    square=True,
    linewidths=0.1,
    cbar_kws={"label": "Pearson r"},
    ax=ax,
)
ax.set(title="Correlation Heatmap")
fig.tight_layout()

# %%
# If you only want to keep the lower triangle of the correlation matrix to have
# unique values of interst, you can use the utility function:

# %%
lower_corr = ca.corr_lower_triangle(data.drop(columns=["group"]), method="pearson")
lower_corr

# %% [markdown]
# Plot the lower triangle correlations as a histrogram to see the distribution of
# correlation values

# %%
ax = lower_corr.stack().plot.hist(
    bins=50,
    grid=False,
    figsize=(6, 4),
    title="Distribution of Pearson correlation values",
    xlabel="Pearson r",
    ylabel="Frequency",
    xlim=(-1.02, 1.02),
)

# %% [markdown]
# or to find the strongest correlations, which you might want to filter further for
# uninteresting correlation between redundant features.

# %%
lower_corr_stack = lower_corr.stack()
idx_largerst_corr = lower_corr_stack.abs().sort_values(ascending=False).head(20).index
lower_corr_stack.loc[idx_largerst_corr]

# %% [markdown]
# This function can be used to compute multiple correlation methods at once
# and compare them, here for the first four features.
#
# It only works on numeric values.

# %%
corr = list()
for method in ["pearson", "spearman", "kendall"]:
    _corr = (
        ca.corr_lower_triangle(data.iloc[:, :4], method=method, numeric_only=True)
        .stack()
        .rename(method)
    )
    corr.append(_corr)
corr = pd.concat(corr, axis=1).sort_values(by="pearson", ascending=True)
corr.plot(
    style=".",
    ylim=(-1.05, 1.05),
    alpha=0.5,
    rot=45,
)

# %% [markdown]
# Filtering correlations based on p-values with multiple testing correction
# - the p-value depends on the number of samples
# - and the strenght of the correlation

# %%
res = ca.calculate_correlations(data.iloc[:, 0], data.iloc[:, 1], method="pearson")
print(res)

# %% [markdown]
# For the first four features, we would only keep one significant correlation
# after multiple testing correction with the Benjamini-Hochberg method.

# %%
correlation = ca.run_correlation(
    data.iloc[:, :4], alpha=0.05, group="group", method="pearson", correction="fdr_bh"
)
correlation

# %% [markdown]
# The efficient correlation calculation can be used to compute the correlation
# matrix and p-value matrix for larger datasets.

# %%
corr, p = ca.run_efficient_correlation(data.iloc[:, :4], method="spearman")
pd.DataFrame(p)

# %% [markdown]
# you can verify the results against [`scipy.stats.spearmanr`](scipy.stats.spearmanr)

# %%
r, p = scipy.stats.spearmanr(data.iloc[:, :4])
pd.DataFrame(p)

# %% [markdown]
# same for pearson correlation

# %%
r, p = ca.run_efficient_correlation(data.iloc[:, :3], method="pearson")
pd.DataFrame(p)

# %%
r_20, p_20 = scipy.stats.pearsonr(data.iloc[:, 0], data.iloc[:, 2])
r_20, p_20
assert r[2, 0] - r_20 < 1e-8
assert p[2, 0] - p_20 < 1e-8

# %% [markdown]
# To calculate p-values for the correlation matrix, you can use

# %%
res = ca.calculate_pvalue_correlation_sample_in_rows(
    data.iloc[:, :3].corr(method="pearson").values, n_obs=data.shape[0]
)
pd.DataFrame(res)

# %% [markdown]
# Done.
