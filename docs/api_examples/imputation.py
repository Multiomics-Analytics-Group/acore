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
# - k Neareast Neighbour imputation can also be used with other types of data
# - the replacement from the normal distrubtion on the sample level is typical to
#   normally distributed samples from massspectrometer data (in the log2 space)
#
# Refers to the [`acore.imputation`](acore.imputation) module.

# %% tags=["hide-output"]
# %pip install acore

# %% tags=["hide-input"]
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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

# %%
omics = omics_and_meta.drop(columns=[*factor_and_covars, *drop_cols])
na_counts = omics.isna().sum().sort_values(ascending=False)
na_counts.plot(
    rot=45,
    style=".",
    alpha=0.5,
    ylabel=f"Number of missing values of {omics.shape[0]} samples",
)

# %%
(na_counts / omics.shape[0]).plot(
    rot=45,
    style=".",
    alpha=0.5,
    ylabel="Ratio of missing values",
)

# %%
not_nan_counts = omics.notna().sum().sort_values(ascending=False)
(not_nan_counts / omics.shape[0]).plot(
    rot=45,
    style=".",
    alpha=0.5,
    ylabel="Ratio of non-missing values\n(completeness)",
)


# %%
omics_and_y = omics_and_meta.drop(columns=[*covariates, *drop_cols])

# %%
omics_and_meta.notna().any(axis=None)

# %%
omics_and_y.loc[omics_and_y.isna().any(axis=1)].loc[:, omics_and_y.isna().any(axis=0)]

# %% [markdown]
# ## KNN imputation
#
# > Can be generally applied
# - both by group and overall

# %% [markdown]
# ### per group
# - does return the original shape, potentially with remaining
#   missing values

# %%
omics_and_y_imputed = imputation_KNN(
    data=omics_and_y,
    drop_cols=[],
    group=group,
    cutoff=0.65,  # selected to leave some missing values for demonstration
    alone=True,  # ! should return only columns without missing values
)
omics_and_y_imputed.isna().sum().value_counts()

# %% [markdown]
# As we have increase the threshold `cutoff` for the fraction of non-misisng
# values per feature, the more features will not be imputed and therefore have
# missing values.

# %%
view = omics_and_y_imputed.loc[:, omics_and_y_imputed.isna().any(axis=0)].loc[
    omics_and_y_imputed.isna().any(axis=1)
]
view

# %% [markdown]
# ### overall


# %%
omics_and_y_imputed = imputation_KNN(
    data=omics_and_y,
    drop_cols=[],
    group=None,
    cutoff=0.6,
    alone=False,
)
assert omics_and_y_imputed.isna().sum().sum() == 0
omics_and_y_imputed

# %%
omics_and_y_imputed = imputation_KNN(
    data=omics_and_y,
    drop_cols=[],
    group=None,
    cutoff=0.90,
    alone=False,
)
omics_and_y_imputed.isna().sum().sum()

# %% [markdown]
#

# %%
# ToDo: Plot per group missingness

# %%
omics_and_y_imputed = imputation_KNN(
    data=omics_and_y,
    drop_cols=[],
    group=group,
    cutoff=0.90,
    alone=False,
)
omics_and_y_imputed.isna().sum().sum()

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

# %%
# ToDo: add an simulated example

# %%
# ! does not account for groups
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
# Therefore many use in proteomics a combined approach (Santos ...)

# %% [markdown]
# ## Combining KNN based imputation and random imputation from a shifted random
# distribution
# - For features (e.g. protein groups) that are present across groups in high enough
#   frequency, use KNN-based imputation (which is deterministic)
# - for the remaining missing values, use based on the distribution of observed values
#   in a sample a shifted normal distribution to draw replacements (random,
#   but deterministic due to the set seed)

# %%
imputation_mixed_norm_KNN(data=omics_and_y, drop_cols=[], group=group, cutoff=0.9)

# %% [markdown]
# done.
