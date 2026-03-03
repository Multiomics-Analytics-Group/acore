# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: tags,-all
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

import dsp_pandas
import matplotlib.pyplot as plt
import pandas as pd
import sklearn
import sklearn.impute
import sklearn.preprocessing
import vuecore.decomposition

import acore.batch_correction
import acore.decomposition


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


dsp_pandas.format.set_pandas_options(
    max_columns=9,
    max_colwidth=20,
)

# %% [markdown]
#
# ## Set some parameters

# %% tags=["parameters"]
BASE = (
    "https://raw.githubusercontent.com/Multiomics-Analytics-Group/acore/"
    "update_batch_corr_example/example_data/alzheimer_proteomics/"
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
# ## Before batch correction
# Explore data in PCA and UMAP space before batch correction

# %% tags=["hide-input"]
omics_imp = median_impute(omics)
omics_imp_scaled = standard_normalize(omics_imp)
PCs, fig = run_and_plot_pca(omics_imp, y, n_components=4)
ax = plot_umap(omics_imp, y)

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
    X.join(y.astype("category")),
    batch_col=y.name,
)
X

# %% [markdown]
# Plot PCA and UMAP after batch correction on standard normalized data

# %% tags=["hide-input"]
PCs, fig = run_and_plot_pca(standard_normalize(X), y, n_components=4)
ax = plot_umap(X, y)

# %% [markdown]
# See change by substracting combat corrected data from original data.
# - NAs in original data will remain NA below (no imputation done here)

# %% tags=["hide-input"]
omics - X

# %% [markdown]
# Done.
