# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.3
#   kernelspec:
#     display_name: acore-dev
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Imputation (Metabolomics example)

# %% [markdown]
# In this notebook, we will showcase the acore functions for imputing data which are specific for metabolomics data analysis. Namely, we will go through imputation with zeros and half-minimum imputation.
#
# For this, we will use a Diabetes example data set from this paper: [Barranco-Altirriba M et al., 2025](https://doi.org/10.3389/fendo.2025.1706886htt)
#
# This notebook refers to the [`acore.imputation_analysis`](acore.imputation_analysis) module.

# %% tags=["hide-output"]
# %pip install acore

import matplotlib.pyplot as plt
import numpy as np

# %% tags=["hide-input"]
import pandas as pd

from acore.imputation_analysis import (
    imputation_half_minimum,
    imputation_zeros,
)


def plot_feature_missingness(data):
    missing_features = data.isnull().mean() * 100
    missing_features_nonzero = missing_features[missing_features > 0]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Histogram
    ax = axes[0]
    ax.hist(
        missing_features_nonzero.values,
        bins=30,
        color="mediumvioletred",
        edgecolor="white",
        linewidth=0.5,
    )
    ax.set_xlabel("Missing (%)")
    ax.set_ylabel("Number of features")
    ax.set_title(
        f"Distribution of missingness\n({len(missing_features_nonzero)} features with any missing)"
    )
    ax.axvline(
        x=20, color="black", linestyle="--", linewidth=0.8, label="20% threshold"
    )
    ax.legend()

    # Dot plot, sorted by missingness
    ax = axes[1]
    sorted_missing = missing_features.sort_values(ascending=True).reset_index(drop=True)
    ax.scatter(
        sorted_missing.index,
        sorted_missing.values,
        s=4,
        color="mediumvioletred",
        alpha=0.6,
        linewidths=0,
    )
    ax.set_xlabel("Features (sorted by missingness)")
    ax.set_ylabel("Missing (%)")
    ax.set_title("Sorted missingness per feature")
    ax.axhline(
        y=20, color="black", linestyle="--", linewidth=0.8, label="20% threshold"
    )
    ax.legend()

    plt.tight_layout()
    plt.show()


def missingness_summary(df_before, df_after):
    total = df_before.size
    n_before = df_before.isnull().sum().sum()
    n_after = df_after.isnull().sum().sum()

    print(f"Total values      : {total:,}")
    print(f"Missing before    : {n_before:,}  ({100*n_before/total:.1f}%)")
    print(f"Missing after     : {n_after:,}  ({100*n_after/total:.1f}%)")
    print(
        f"Features affected : {(df_before.isnull().any()).sum()} / {df_before.shape[1]}"
    )
    print(
        f"Samples affected  : {(df_before.isnull().any(axis=1)).sum()} / {df_before.shape[0]}"
    )


def plot_intensity_distribution(data):
    values = data.values.flatten().astype(float)

    n_total = data.size
    n_missing = int(np.isnan(values).sum())
    pct_missing = n_missing / n_total * 100

    log_values = np.log10(values[(~np.isnan(values)) & (values > 0)])

    fig, (ax_nan, ax_hist) = plt.subplots(
        1, 2, figsize=(11, 5), gridspec_kw={"width_ratios": [1, 8]}, sharey=True
    )

    ax_nan.bar(0, n_missing, width=0.2, color="mediumvioletred", alpha=0.8)
    ax_nan.set_xlim(-0.5, 0.5)
    ax_nan.set_xticks([0])
    ax_nan.set_xticklabels([f"NaN\n({pct_missing:.1f}%)"], fontsize=9)
    ax_nan.set_ylabel("Count")

    ax_hist.hist(
        log_values, bins=100, color="cornflowerblue", edgecolor="none", alpha=0.8
    )
    ax_hist.set_xlabel("Intensity (log₁₀)")
    ax_hist.set_title("Intensity distribution (all features, all samples)")
    ax_hist.yaxis.set_visible(False)

    plt.tight_layout()
    plt.show()


# %% [markdown]
# ### Data Loading
#
# Load in your data and inspect the resulting dataframe. The example data set can be found in example_data/DidacMauricio_hilic.
#
# The data set has been filtered already, using the [`acore.filter_metabolomics`](acore.filter_metabolomics) module. That means that features with a lot of missingness have been filtered out already, meaning that the features that are remaining have limited missingness and the data set is ready for the imputation step.

# %%
data_path = (
    "https://raw.githubusercontent.com/Multiomics-Analytics-Group/acore/"
    "refs/heads/main/"
)
data_path = "../../example_data/DidacMauricio_hilic/DM_FIS2018_Hilic_pos_results2023_filtered.csv"
data_original = pd.read_csv(data_path, index_col=0)

# %%
data_original

# %% [markdown]
# In order to run our further analysis, including the filtering functions, we have to transform the data and remove metadata such as mass and retention time.

# %%
data = data_original.T
data = data.drop(
    ["Qidx", "SOIidx", "rtmed", "start", "end", "mass", "MaxInt", "formula", "anot"]
)

# %% tags=["hide-input"]
data

# %% [markdown]
# Check how much missingness there is in the data.

# %%
# Total missing count and percentage
print(
    f"Total count of missing cells: {data.isnull().sum().sum()}"
)  # total missing cells
print(
    f"Overall percentage of missingness: {data.isnull().mean().mean() * 100}\n"
)  # overall % missing

plot_feature_missingness(data)
plot_intensity_distribution(data)

# %% [markdown]
# As we can see, overall, 1.7% of the dataset is missing. There are some features that have a lot of missing values, whereas most have very few.
#
# Now that we have an overview, we can try two different methods of imputation.

# %% [markdown]
# ### Imputing with zeros
#
# In this method, which is commonly used in metabolomics and often automatically done by preprocessing softwared like MetaboIgniter, all missing values get filled in with zeros.
#
# Here, this method can be applied easily using the function [`imputation_zeros()`](acore.imputation_analysis.imputation_zeros).

# %% tags=["hide-input"]
help(imputation_zeros)

# %%
data_imputed_zeros = imputation_zeros(data=data)

# %%
# Total missing count and percentage
print(
    f"Total count of missing cells: {data_imputed_zeros.isnull().sum().sum()}"
)  # total missing cells
print(
    f"Overall percentage of missingness: {data_imputed_zeros.isnull().mean().mean() * 100}\n"
)  # overall % missing

print("SUMMARY of imputation changes:")
missingness_summary(data, data_imputed_zeros)
plot_intensity_distribution(data_imputed_zeros)

# %% [markdown]
# We have imputed all of our missing values with zeros.

# %% [markdown]
# ### Imputation with half minimum
#
# This method is also widely used across the metabolomics community. Here, missing values are imputed with half of the minimum value that has been recorded across the data set.
#
# This is done following the assumption of missing-not-at-random; that measurements may be missing not because they are truly absent in the biological sample but because they are for example below the limit of detection.
#
# Here, in acore, the function [`imputation_half_minimum()`](acore.imputation_analysis.imputation_half_minimum) is used for this.

# %% tags=["hide-input"]
help(imputation_half_minimum)

# %%
data_imputed_hm = imputation_half_minimum(data)

# %%
# Total missing count and percentage
print(
    f"Total count of missing cells: {data_imputed_hm.isnull().sum().sum()}"
)  # total missing cells
print(
    f"Overall percentage of missingness: {data_imputed_hm.isnull().mean().mean() * 100}\n"
)  # overall % missing

print("SUMMARY of imputation changes:")
missingness_summary(data, data_imputed_hm)
plot_intensity_distribution(data_imputed_hm)

# %% [markdown]
# Again, we can see that after imputation, no missing values are left in our data.
