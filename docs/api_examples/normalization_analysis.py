# ---
# jupyter:
#   jupytext:
#     cell_data_filter: tags,-all
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Normalization of samples
#
# We will explore an Alzheimer dataset where the data was collected in four different sites.
# We will see that the sites have a an effect where the data is in principal component space
# and in UMAP space. We will then normalize the data and see how the effect on these plots.
#
# Refers to the [`acore.normalization`](acore.normalization) module.

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
import acore.normalization


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

# %% [markdown]
# ## Data loading
# Use combined dataset for ANCOVA analysis.

# %% tags=["hide-input"]
omics_and_meta = (
    pd.read_csv(f"{BASE}/{fname}", index_col=subject_col)
    .convert_dtypes()
    .dropna(subset=factor_and_covars)
)
omics_and_meta

# %% [markdown]
# Metadata here is of type integer. All floats are proteomics measurements.

# %% tags=["hide-input"]
omics_and_meta.dtypes.value_counts()

# %% tags=["hide-input"]
omics_and_meta[factor_and_covars]

# %%
omics = omics_and_meta.drop(columns=[*factor_and_covars, *drop_cols])
y = omics_and_meta[group].astype("category").rename(group_label)

# %% [markdown]
# For simplicity we normalize here all samples together, but normally you would need to
# apply the normalization from you training data to the test data. So see these examples
# here as a way to do it for your training data.

# %% [markdown]
# ### Fill missing values for preliminary plots

# %% [markdown]
# Impute using median to impute (before scaling, which can be changed).

# %% tags=["hide-input"]
omics_imputed = median_impute(omics)
assert omics_imputed.isna().sum().sum() == 0
omics_imputed.shape

# %% [markdown]
# Explained variance by first four principal components in data.

# %% tags=["hide-input"]
PCs, pca = acore.decomposition.pca.run_pca(omics_imputed, n_components=4)
ax = vuecore.decomposition.plot_explained_variance(pca)
ax.locator_params(axis="x", integer=True)

# %% [markdown]
# ## Normalization of samples in a dataset
# We will use the `acore.normalization` module to normalize the data.
#
# We will do it for each of the data on the omics dataset which is log transformed,
# but not yet imputed and normalized. Then we will reapply standard
# normalization before replotting the PCA and UMAP plots. The execption is combat as it
# need complete data.

# %% tags=["hide-input"]
omics


# %% [markdown]
# ## Median normalization
# Substracts a constant from all features of a sample. All samples will have the same
# global median.

# %%
# %%time
X = acore.normalization.normalize_data(omics, "median")
X

# %% tags=["hide-input"]
omics_imp = median_impute(X)
omics_imp_scaled = standard_normalize(omics_imp)
PCs, fig = run_and_plot_pca(omics_imp_scaled, y, y.name, n_components=4)
ax = plot_umap(omics_imp_scaled, y)

# %% [markdown]
# See change by substracting median normalized data from original data.

# %% tags=["hide-input"]
omics - X

# %% [markdown]
# ## Z-score normalization
# Normalize a sample by it's mean and standard deviation.

# %%
# %%time
X = acore.normalization.normalize_data(omics, "zscore")
X

# %% tags=["hide-input"]
omics_imp = median_impute(X)
omics_imp_scaled = standard_normalize(omics_imp)
PCs, fig = run_and_plot_pca(omics_imp_scaled, y, n_components=4)
ax = plot_umap(omics_imp_scaled, y)

# %% [markdown]
# See change by substracting z-score normalized data from original data.

# %% tags=["hide-input"]
omics_imp_scaled - X

# %% [markdown]
# ## Median Polish Normalization
# - normalize iteratively features and samples to have zero median.

# %%
# %%time
X = acore.normalization.normalize_data(omics, "median_polish")
X

# %% tags=["hide-input"]
omics_imp = median_impute(X)
omics_imp_scaled = standard_normalize(omics_imp)
PCs, fig = run_and_plot_pca(omics_imp_scaled, y, n_components=4)
ax = plot_umap(omics_imp_scaled, y)

# %% [markdown]
# See change by substracting median polish normalized data from original data.

# %% tags=["hide-input"]
omics_imp_scaled - X

# %% [markdown]
# ## Quantile normalization
# quantile normalize each feature column.

# %%
# %%time
X = acore.normalization.normalize_data(omics, "quantile")
X

# %% tags=["hide-input"]
omics_imp = median_impute(X)
omics_imp_scaled = standard_normalize(omics_imp)
PCs, fig = run_and_plot_pca(omics_imp_scaled, y, n_components=4)
ax = plot_umap(omics_imp_scaled, y)

# %%
omics - X


# %% [markdown]
# ## Linear normalization

# %%
# %%time
X = acore.normalization.normalize_data(omics, "linear")
X

# %% tags=["hide-input"]
omics_imp = median_impute(X)
omics_imp_scaled = standard_normalize(omics_imp)
PCs, fig = run_and_plot_pca(omics_imp_scaled, y, n_components=4)
ax = plot_umap(omics_imp_scaled, y)

# %% tags=["hide-input"]
omics - X

# %% [markdown]
# ## Summmary
# Besides the median polish normalization, the structure of the data is not changed
# too much by the normalization using this Alzheimer example. This notebook can be opened
# on colab and might be a good starting point for investigating the effect of normalization
# on your data - or to disect some approaches further.
