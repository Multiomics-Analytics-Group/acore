# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: tags,-all
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
# # Batch correction of samples
#
# We will explore an Alzheimer dataset where the data was collected in four different sites.
# We will see that the sites have a an effect where the data is in principal component space
# and in UMAP space. We will then batch correct the data and see how the effect on these plots.
#
# Refers to the [`acore.batch_correction`](acore.batch_correction) module.


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

import acore.batch_correction
import acore.decomposition


def plot_umap(X_scaled, y, meta_column, random_state=42) -> plt.Axes:
    """Fit and plot UMAP embedding with two components with colors defined by meta_column."""
    embedding = acore.decomposition.umap.run_umap(
        X_scaled, y, random_state=random_state
    )
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
# ## Before batch correction
# Explore data in PCA and UMAP space before batch correction


# %% tags=["hide-input"]
omics_imp = median_impute(omics)
omics_imp_scaled = standard_normalize(omics_imp)
PCs, fig = run_and_plot_pca(omics_imp, y, METACOL_LABEL, n_components=4)
ax = plot_umap(omics_imp, y, METACOL_LABEL)


# %% [markdown]
# ## Combat batch correction
# Correct for batch effects in the data using a robust regression approach removing
# mean and scale effects out for each provided co-variate by batch.
# Assumes normally distributed data.
#
# > ⚠️ Combat needs imputed data

# %%
# %%time
X = median_impute(omics)
X = acore.batch_correction.combat_batch_correction(
    X.join(y),
    batch_col="site",
)
X

# %% [markdown]
# Plot PCA and UMAP after batch correction on standard normalized data

# %% tags=["hide-input"]
PCs, fig = run_and_plot_pca(standard_normalize(X), y, METACOL_LABEL, n_components=4)
ax = plot_umap(X, y, METACOL_LABEL)

# %% [markdown]
# See change by substracting combat corrected data from original data.
# - NAs in original data will remain NA below (no imputation done here)

# %% tags=["hide-input"]
omics - X
