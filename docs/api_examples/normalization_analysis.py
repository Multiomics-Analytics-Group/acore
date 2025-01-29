# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: tags,-all
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Normalization of samples in a dataset example
#
# We will explore an Alzheimer dataset where the data was collected in four different sites.
# We will see that the sites have a an effect where the data is in principal component space
# and in UMAP space. We will then normalize the data and see how the effect on these plots.
#
# Refers to the `acore.normalization` module.

# %% tags=["hide-output"]
# %pip install acore

# %% tags=["hide-input"]
import itertools
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn
import sklearn.impute
import sklearn.preprocessing
import umap
from vuecore.decomposition import plot_explained_variance

import acore.normalization
import acore.sklearn
from acore.decomposition import pca as acore_pca  # ! to remove


def plot_umap(X_scaled, y, meta_column, random_state=42) -> plt.Axes:
    """Fit and plot UMAP embedding with two components with colors defined by meta_column."""
    reducer = umap.UMAP(random_state=random_state, n_jobs=1)
    embedding = reducer.fit_transform(X_scaled)
    embedding = pd.DataFrame(
        embedding, index=X_scaled.index, columns=["UMAP 1", "UMAP 2"]
    ).join(y.astype("category"))
    ax = embedding.plot.scatter("UMAP 1", "UMAP 2", c=meta_column, cmap="Paired")
    return ax


def standard_normalize(X: pd.DataFrame) -> pd.DataFrame:
    """Standard normalize data and keep indices of DataFrame."""
    scaler = sklearn.preprocessing.StandardScaler()
    X_scaled = acore.sklearn.transform_DataFrame(X, fct=scaler.fit_transform)
    return X_scaled


def median_impute(X: pd.DataFrame) -> pd.DataFrame:
    median_imputer = sklearn.impute.SimpleImputer(strategy="median")
    X_imputed = acore.sklearn.transform_DataFrame(X, median_imputer.fit_transform)
    return X_imputed


def run_and_plot_pca(
    X_scaled,
    y,
    meta_column,
    n_components=4,
) -> tuple[pd.DataFrame, plt.Figure]:
    PCs, _ = acore_pca.run_pca(X_scaled, n_components=n_components)
    PCs.columns = [s.replace("principal component", "PC") for s in PCs.columns]
    PCs = PCs.join(y.astype("category"))
    up_to = min(PCs.shape[-1], n_components)
    fig, axes = plt.subplots(up_to - 1, 2, figsize=(6, 8), layout="constrained")
    for k, (pos, ax) in enumerate(
        zip(itertools.combinations(range(up_to), 2), axes.flatten())
    ):
        i, j = pos
        plot_heatmap = bool(k % 2)
        PCs.plot.scatter(
            i, j, c=meta_column, cmap="Paired", ax=ax, colorbar=plot_heatmap
        )
    _ = PCs.pop(
        meta_column,
    )
    return PCs, fig


# %% [markdown]
#
# ## Set some parameters

# %% tags=["parameters"]
fname_metadata: str = (
    "https://raw.githubusercontent.com/RasmussenLab/"
    "njab/HEAD/docs/tutorial/data/alzheimer/meta.csv"  # clincial data
)
fname_omics: str = (
    "https://raw.githubusercontent.com/RasmussenLab/"
    "njab/HEAD/docs/tutorial/data/alzheimer/proteome.csv"  # omics data
)
METACOL: str = "_collection site"  # target column in fname_metadata dataset (binary)
METACOL_LABEL: Optional[str] = "site"  # optional: rename target variable
n_features_max: int = 5
freq_cutoff: float = 0.5  # Omics cutoff for sample completeness
VAL_IDS: str = ""  #
VAL_IDS_query: str = ""
weights: bool = True
FOLDER = "alzheimer"
model_name = "all"

# %% [markdown]
# ## Setup

# %% [markdown]
# ### Load proteomics (protein groups) data

# %%
if METACOL_LABEL is None:
    METACOL_LABEL = METACOL
metadata = (
    pd.read_csv(fname_metadata, usecols=["Sample ID", METACOL], index_col=0)
    .convert_dtypes()
    .rename(columns={METACOL: METACOL_LABEL})
)
omics = pd.read_csv(fname_omics, index_col=0)

# %% [markdown]
# Data shapes

# %%
omics.shape, metadata.shape

# %% [markdown]
# See how common omics features are and remove feature below choosen frequency cutoff

# %%
ax = omics.notna().sum().sort_values().plot(rot=90)

# %% tags=["hide-input"]
M_before = omics.shape[1]
omics = omics.dropna(thresh=int(len(omics) * freq_cutoff), axis=1)
M_after = omics.shape[1]
msg = (
    f"Removed {M_before-M_after} features with more "
    f"than {freq_cutoff*100}% missing values."
    f"\nRemaining features: {M_after} (of {M_before})"
)
print(msg)
# keep a map of all proteins in protein group, but only display first protein
# proteins are unique to protein groups
pg_map = {k: k.split(";")[0] for k in omics.columns}
omics = omics.rename(columns=pg_map)
# log2 transform raw intensity data:
omics = np.log2(omics + 1)
ax = (
    omics.notna()
    .sum()
    .sort_values()
    .plot(
        rot=90,
        ylabel="Number of samples",
        xlabel="Proteins (ranked by missing values)",
    )
)
omics

# %% [markdown]
# ### Sample metadata

# %%
metadata

# %% [markdown]
# Tabulate selected metadata and check for missing values

# %% tags=["hide-input"]
metadata[METACOL_LABEL].value_counts(dropna=False)

# %% tags=["hide-input"]
target_counts = metadata[METACOL_LABEL].value_counts()

if target_counts.sum() < len(metadata):
    print(
        "Target has missing values."
        f" Can only use {target_counts.sum()} of {len(metadata)} samples."
    )
    mask = metadata[METACOL_LABEL].notna()
    metadata, omics = metadata.loc[mask], omics.loc[mask]

if METACOL_LABEL is None:
    METACOL_LABEL = METACOL_LABEL
y = metadata[METACOL_LABEL].astype("category")

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
# ## Dimensionality reduction - unnormalized data
# on median imputed and standard normalized omics data.

# %% [markdown]
# ### Principal Components
# Plot first 4 PCs with categorical metadata as label annotating each sample.

# %% tags=["hide-input"]
omics_imp_scaled = standard_normalize(omics_imputed)

PCs, pca = acore_pca.run_pca(omics_imp_scaled, n_components=4)
ax = plot_explained_variance(pca)
ax.locator_params(axis="x", integer=True)
omics_imp_scaled.shape

# %% tags=["hide-input"]
pcs, fig = run_and_plot_pca(omics_imp_scaled, y, METACOL_LABEL)

# %% [markdown]
# ### UMAP
# of median imputed and normalized omics data:

# %% tags=["hide-input"]
ax = plot_umap(omics_imp_scaled, y, METACOL_LABEL)

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
# ## Combat normalization
# Correct for batch effects in the data using a robust regression approach normalizing
# mean and scale effetcs out for each feature by batch. Assumes normally distributed data.
#
# > ⚠️ Combat needs imputed data

# %%
# %%time
X = median_impute(omics)
X = acore.normalization.combat_batch_correction(
    X.join(y),
    batch_col="site",
)
X

# %% tags=["hide-input"]
omics_imp_scaled = standard_normalize(omics_imp_scaled)
PCs, fig = run_and_plot_pca(omics_imp_scaled, y, METACOL_LABEL, n_components=4)
ax = plot_umap(omics_imp_scaled, y, METACOL_LABEL)

# %% [markdown]
# See change by substracting combat normalized data from original data.

# %% tags=["hide-input"]
omics - X


# %% [markdown]
# ## Median normalization
# Substracts a constant from all features of a sample. All samples will have the same
# global median.

# %%
# %%time
X = acore.normalization.normalize_data(omics, "median")
X

# %% tags=["hide-input"]
omics_imp_scaled = median_impute(X)
omics_imp_scaled = standard_normalize(omics_imp_scaled)
PCs, fig = run_and_plot_pca(omics_imp_scaled, y, METACOL_LABEL, n_components=4)
ax = plot_umap(omics_imp_scaled, y, METACOL_LABEL)

# %% [markdown]
# See change by substracting median normalized data from original data.

# %% tags=["hide-input"]
omics - X

# %%
# %% [markdown]
# ## Z-score normalization
# Normalize a sample by it's mean and standard deviation.

# %%
# %%time
X = acore.normalization.normalize_data(omics, "zscore")
X

# %% tags=["hide-input"]
omics_imp_scaled = median_impute(X)
omics_imp_scaled = standard_normalize(omics_imp_scaled)
PCs, fig = run_and_plot_pca(omics_imp_scaled, y, METACOL_LABEL, n_components=4)
ax = plot_umap(omics_imp_scaled, y, METACOL_LABEL)

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
omics_imp_scaled = median_impute(X)
omics_imp_scaled = standard_normalize(omics_imp_scaled)
PCs, fig = run_and_plot_pca(omics_imp_scaled, y, METACOL_LABEL, n_components=4)
ax = plot_umap(omics_imp_scaled, y, METACOL_LABEL)

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
omics_imp_scaled = median_impute(X)
omics_imp_scaled = standard_normalize(omics_imp_scaled)
PCs, fig = run_and_plot_pca(omics_imp_scaled, y, METACOL_LABEL, n_components=4)
ax = plot_umap(omics_imp_scaled, y, METACOL_LABEL)

# %%
omics - X


# %% [markdown]
# ## Linear normalization

# %%
# %%time
X = acore.normalization.normalize_data(omics, "linear")
X

# %% tags=["hide-input"]
omics_imp_scaled = median_impute(X)
omics_imp_scaled = standard_normalize(omics_imp_scaled)
PCs, fig = run_and_plot_pca(omics_imp_scaled, y, METACOL_LABEL, n_components=4)
ax = plot_umap(omics_imp_scaled, y, METACOL_LABEL)

# %% tags=["hide-input"]
omics - X

# %% [markdown]
# ## Summmary
# Besides the median polish normalization, the structure of the data is not changed
# too much by the normalization using this Alzheimer example. This notebook can be opened
# on colab and might be a good starting point for investigating the effect of normalization
# on your data - or to disect some approaches further.
