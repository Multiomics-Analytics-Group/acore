# ---
# jupyter:
#   jupytext:
#     cell_data_filter: tags,-all
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Imputation data (MS-example)
#
# We will explore imputation of proteomics data using an Alzheimer dataset where the
# data was collected in four different sites.
# - k Nearest Neighbour imputation can also be used with other types of data
# - the replacement from the normal distribution on the sample level is typical to
#   normally distributed samples from mass spectrometer data (in the log2 space)
#
# Refers to the [`acore.imputation_analysis`](acore.imputation_analysis) module.
#
# Common shared parameters across: `imputation_KNN`, `imputation_normal_distribution` and
# `imputation_mixed_norm_KNN` function presented here:
#
# - `data`: `pd.DataFrame` with samples as rows and features as columns, which can 
#   can contain a `group` column.
# - `drop_cols`: optional iterable of column names excluded from imputation.


# %% tags=["hide-output"]
# %pip install acore

# %% tags=["hide-input"]
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import sklearn
import sklearn.impute
import sklearn.preprocessing
import vuecore.decomposition

import acore.decomposition
from acore.imputation_analysis import (
    imputation_KNN,
    imputation_mixed_norm_KNN,
    imputation_normal_distribution,
)


def plot_umap(X_scaled, y, meta_column=None, random_state=42) -> plt.Axes:
    """Fit and plot UMAP embedding with two components with colors defined by meta_column."""
    embedding = acore.decomposition.umap.run_umap(
        X_scaled, y, random_state=random_state
    )
    if meta_column is None:
        meta_column = y.name
    ax = embedding.plot.scatter("UMAP 1", "UMAP 2", c=meta_column, cmap="Paired")
    return ax


def standard_normalize(X: pd.DataFrame) -> pd.DataFrame:
    """Standard normalize data and keep indices of DataFrame."""
    X_scaled = (
        sklearn.preprocessing.StandardScaler()
        .set_output(transform="pandas")
        .fit_transform(X)
    )
    return X_scaled


def median_impute(X: pd.DataFrame) -> pd.DataFrame:
    X_imputed = (
        sklearn.impute.SimpleImputer(strategy="median")
        .set_output(transform="pandas")
        .fit_transform(X)
    )
    return X_imputed


def run_and_plot_pca(
    X_scaled,
    y,
    meta_column: Optional[str] = None,
    n_components: int = 4,
) -> tuple[pd.DataFrame, plt.Figure]:
    PCs, _ = acore.decomposition.pca.run_pca(X_scaled, n_components=n_components)
    PCs.columns = [s.replace("principal component", "PC") for s in PCs.columns]
    fig = vuecore.decomposition.pca_grid(
        PCs=PCs, meta_column=y, n_components=n_components, meta_col_name=meta_column
    )
    return PCs, fig


# %% [markdown]
#
# ## Set some parameters

# %% tags=["parameters"]
BASE = (
    "https://raw.githubusercontent.com/Multiomics-Analytics-Group/acore/"
    "main/example_data/alzheimer_proteomics/"
)
# data is already preprocessed: log2, filtered
fname: str = "alzheimer_example_omics_and_clinic.csv"  # combined omics and meta data
covariates: list[str] = ["age", "male"]
group: str = "collection_site"
subject_col: str = "Sample ID"
drop_cols: list[str] = ["AD"]
factor_and_covars: list[str] = [group, *covariates]
group_label: Optional[str] = "site"  # optional: rename target variable

# %%
# %% [markdown]
# ## Data loading
# Use combined dataset for ANCOVA analysis.

# %% tags=["hide-input"]
omics_and_meta = pd.read_csv(f"{BASE}/{fname}", index_col=subject_col).convert_dtypes()
omics_and_meta

# %% [markdown]
# Separate omics and the grouping variable

# %% tags=["hide-input"]
omics = omics_and_meta.drop(columns=[*factor_and_covars, *drop_cols])
na_counts = omics.isna().sum().sort_values(ascending=False)

fig, axes = plt.subplots(1, 2, figsize=(8, 4), constrained_layout=True)

na_counts.plot(
    ax=axes[0],
    rot=45,
    style=".",
    alpha=0.5,
    ylabel=f"Number of missing values of {omics.shape[0]} samples",
    title="Missing values (count)",
)

ax1 = axes[1]
ax2 = ax1.twinx()

(na_counts / omics.shape[0]).plot(
    ax=ax1,
    rot=45,
    style=".",
    alpha=0.5,
    color="C0",
    ylabel="Ratio of missing values",
    title="Missing values & Completeness (ratios)",
)
not_nan_counts = omics.notna().sum().sort_values(ascending=False)
(not_nan_counts / omics.shape[0]).plot(
    ax=ax2,
    rot=45,
    style=".",
    alpha=0.0,
    color="C1",
    ylabel="Ratio of non-missing values\n(completeness)",
)
ax1.tick_params(axis="y", labelcolor="C0")
ax2.tick_params(axis="y", labelcolor="C1")

# %% [markdown]
# Show samples and features with at least 3 missing values

# %% tags=["hide-input"]
omics_and_y = omics_and_meta.drop(columns=[*covariates, *drop_cols])
assert omics_and_meta.notna().any(axis=None), "Nothing to impute"
omics_and_y.loc[omics_and_y.isna().sum(axis=1) >= 3].loc[
    :, omics_and_y.isna().sum(axis=0) >= 3
]

# %% [markdown]
# ## KNN imputation
#
# > Can be generally applied
# - both by group and overall
# - returns the imputed data, per default only features that meet the criteria
#   based on the selected cutoff for the fraction of non-missing values for a
#   single feature (e.g. protein group).
# - setting `alone=False` will ensure that all features, imputed or not, are returned.
#   This can be useful for downstream analysis where you want to keep all features,
#   but only impute those that meet a minimal quality criteria.

# %% [markdown]
# ### overall

# %%
cutoff = 0.60
omics_and_y_imputed = imputation_KNN(
    data=omics_and_y,
    drop_cols=[],
    group=None,
    cutoff=cutoff,
    alone=False,
)
assert omics_and_y_imputed.isna().sum().sum() == 0
omics_and_y_imputed

# %% [markdown]
# As we have increase the threshold `cutoff` for the fraction of non-misisng
# values per feature, the more features will not be imputed and therefore have
# missing values.

# %%
cutoff = 0.90
omics_and_y_imputed = imputation_KNN(
    data=omics_and_y,
    drop_cols=[],
    group=None,
    cutoff=cutoff,
    alone=False,
)
n_still_missing = omics_and_y_imputed.isna().sum().sum()
print(f"Still missing features with cutoff of {cutoff}: {n_still_missing}")

# %% [markdown]
# ### Keep only imputed features
# Use the `alone=True` to only keep the imputed features. It is the default.

# %%
cutoff = 0.90
omics_and_y_imputed = imputation_KNN(
    data=omics_and_y,
    drop_cols=[],
    group=None,
    cutoff=cutoff,
    alone=True,
)
assert omics_and_y_imputed.isna().sum().sum() == 0
print("Shape of of input data: ", omics_and_y.shape)
print("Shape of imputed data: ", omics_and_y_imputed.shape)

# %% [markdown]
# ### By group
# Do the imputation separately for each group (e.g. target vs control) and
# then combine the results.
#
# Let's see the ratop of missing (left y-axis) and of non-missing (right y-axis) values
# per feature (e.g. protein group) for which no missing values for each group:

# %% tags=["hide-input"]
frac_na_by_group = (
    omics_and_y.groupby(
        group,
    )
    .apply(lambda x: x.isna().sum() / x.shape[0], include_groups=False)
    .T
).sort_values(by="Berlin")
frac_non_na_by_group = (
    omics_and_y.groupby(group)
    .apply(lambda x: x.notna().sum() / x.shape[0], include_groups=False)
    .T
).sort_values(by="Berlin")
fig, ax = plt.subplots(1, 1, figsize=(6, 4), constrained_layout=True)
ax2 = ax.twinx()

for g in frac_na_by_group.columns:
    frac_na_by_group[g].plot(
        ax=ax,
        rot=45,
        style=".",
        alpha=0.5,
        label=g,
        ylabel="Ratio of missing values",
        title="Missing values & Completeness (ratios)",
    )
    frac_non_na_by_group[g].plot(
        ax=ax2,
        rot=45,
        style=".",
        alpha=0.0,
        ylabel="Ratio of non-missing values\n(completeness)",
    )
ax.legend(loc="upper left")
ax.set(xlabel="Protein Groups (sorted by fraction of missing values in Berlin)")
ax.tick_params(axis="y", labelcolor="C0")
ax2.tick_params(axis="y", labelcolor="C1")

# %% [markdown]
# Clearly some protein groups only have missing values if combined from a certain
# collection site, and that the ratio can be different in each group. Therefore,
# imputation by KNN for a threshold of non-missing values per feature
# (e.g. protein group) per group can be a good option.

# %%
omics_and_y_imputed = imputation_KNN(
    data=omics_and_y,
    drop_cols=[],
    group=group,
    cutoff=0.65,  # selected to leave some missing values for demonstration
    alone=False,
)
omics_and_y_imputed.isna().sum().value_counts().sort_index()

# %% [markdown]
# If we look at the number of missing values still remaining by collection site,
# we see that Sweden has most of these missing values due to a higher fraction of
# missing values for these.

# %% tags=["hide-input"]
omics_and_y_imputed.groupby(group).apply(
    lambda x: x.isna().sum(), include_groups=False
).T.value_counts().sort_index()


# %% [markdown]
# As we increase the threshold `cutoff` for the fraction of non-misisng
# values per feature, the more features will not be imputed and therefore have
# missing values.

# %%
omics_and_y_imputed = imputation_KNN(
    data=omics_and_y,
    drop_cols=[],
    group=group,
    cutoff=0.90,
    alone=False,
)
n_still_missing = omics_and_y_imputed.isna().sum().sum()
print("Still missing features with cutoff of {cutoff}: {n_still_missing}")

# %% [markdown]
# ## Imputation from a shifted normal distribution per sample
# > Specific to massspectrometry based data in log2 space: normal distributed data,
# > with detection limit (if that applies it can be used)
# - based on mean and standard deviation missing values are replaced by drawing
#   random values from a shifted normal distribution
# - assumption is that missing values are due to falling below the detection limit
#   which can be revealed by the distribution of intensities
#
# Below you find a generated example highlighting the idea

# %% tags=["hide-input"]
mu = 25.0
stddev = 1.0

x = np.linspace(mu - 3, mu + 3, num=101)

y_normal = scipy.stats.norm.pdf(x, loc=mu, scale=stddev)

mu_shifted = mu - (1.8 * stddev)
stddev_shifted = 0.3 * stddev
print(f"Downshifted: {mu_shifted = }, {stddev_shifted = }")
y_impute = scipy.stats.norm.pdf(x, loc=mu - (1.8 * stddev), scale=0.3 * stddev)

fig, ax = plt.subplots(1, 1, figsize=(5, 4))

for i, (y, label) in enumerate(zip([y_normal, y_impute], ["original", "down shifted"])):
    ax.plot(x, y, color=f"C{i}", label=label)
    ax.fill_between(x, y, color=f"C{i}", alpha=0.5)
    ax.set_label(label)
ax.set_xlabel("log2 intensity distribution in sample")
ax.set_ylabel("density")
ax.legend()
fig.tight_layout()

# %% [markdown]
# This idea can be applied on a per sample basis using:

# %%
# does not account for groups as it is done on a per sample basis (along columns)
imputation_normal_distribution(
    data=omics_and_y,
    drop_cols=[group],
)

# %% [markdown]
# Note that using this type of imputation before differential regulation can lead to
# false positive and negativ results. If values are not due to assumed missing
# mechanism (Missing not-at-random due to low abundance), but are due to technical noise
# these values should not be replace.
#
# Therefore many use in proteomics a combined approach
# [Santos et al., 2020](https://www.nature.com/articles/s41587-021-01145-6):

# %% [markdown]
# ## Combining KNN based imputation and random imputation from a shifted random
# distribution
# - For features (e.g. protein groups) that are present across groups in high enough
#   frequency, use KNN-based imputation (which is deterministic)
# - for the remaining missing values, use based on the distribution of observed values
#   in a sample a shifted normal distribution to draw replacements (random,
#   but deterministic due to the set seed)
# See the methods section of
# [Santos et al., 2020](https://www.nature.com/articles/s41587-021-01145-6)
# for more details.

# %%
imputation_mixed_norm_KNN(data=omics_and_y, drop_cols=[], group=group, cutoff=0.9)

# %% [markdown]
# done.
