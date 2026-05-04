# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: acore-dev
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Metabolomics drift correction
#

# %% [markdown]
# Within-batch correction of metabolomics data.
# This script shows how to correct for instrumental drift based on pooled QC data samples.

# %% [markdown]
# ## 1. LOESS smoothing-based drift correction

# %% [markdown]
# ### Load in example data

# %%
# #%pip install acore

# %%
import importlib
import os

import pandas as pd

import acore
from acore import drift_correction as dc

importlib.reload(acore)

# %% [markdown]
# ### Load in data

# %%
df = pd.read_excel(
    "../../example_data/aradopsis_seedling_lipids/kehelpannala_AnArabidopsisLipid_seedling.xlsx"
)
sample_order = pd.read_csv(
    "../../example_data/aradopsis_seedling_lipids/seedling_art_sample_order.csv"
)

# %%
df

# %% [markdown]
# We also have (artificial) data that has the order in which our samples were run. This information is crucial for the drift corrrection algorithm.

# %%
sample_order.sort_values("Sample ID")

# %% [markdown]
# ### Run drift correction with LOESS smoothing

# %% [markdown]
# We can now correct our data for experimental drift.
#
# With the acore LOESS drift-correction function, a LOESS (locally estimated regression) smoother is applied separately to the features in the data to model slow temporal trends, and the resulting smooth trend is used to correct the data.
#
# Before the function estimation and correction, the data can filtered, to remove features that have too many missing values in the QC samples.

# %%
# First, we can create a dictionary for our sample names, ordering them into groups, to make the upcoming function call easier.

data_groups = {
    "SD1": ["Sd1_1", "Sd1_2", "Sd1_3", "Sd1_4", "Sd1_5", "Sd1_6"],
    "SD2": [
        "Sd2_1",
        "Sd2_2",
        "Sd2_3",
        "Sd2_4",
        "Sd2_5",
        "Sd2_6",
        "Sd2_7",
        "Sd2_8",
        "Sd2_9",
    ],
    "SD3": [
        "Sd3-1",
        "Sd3-2",
        "Sd3-3",
        "Sd3-4",
        "Sd3-5",
        "Sd3-6",
        "Sd3-7",
        "Sd3-8",
        "Sd3-9",
    ],
    "QC": [
        "PBQC_Sd_1",
        "PBQC_Sd_2",
        "PBQC_Sd_3",
        "PBQC_Sd_4",
        "PBQC_Sd_5",
        "PBQC_Sd_6",
    ],
}

# %%
# Separate groups
sample_cols = data_groups["SD1"] + data_groups["SD2"] + data_groups["SD3"]
qc_cols = data_groups["QC"]

# Create a sub-df consisting only of the interesting columns, omitting metadata.
df_dc = df[sample_cols + qc_cols]

# %% [markdown]
# Now we can run the drift correction.

# %%
corrected_df, correction_info = dc.run_drift_correction(
    df,
    qc_cols,
    sample_cols,
    sample_order=sample_order,
    feature_name_col=None,
    filter_percent=0.5,
    print_logs=True,
)

### Explanation of the parameters chosen
# - feature_name_col = the name of the column containing feature names, if there is one.
#   This information is used for logging and showing outputs, it's not required for the functioning of the method.
#   Here, there is no feature name column available, so "None" is used.
# - filter_percent =  the minimum percentage of values that must be present for this feature to be retained.
#   If the percentage of non-missing is below this, the feature will be filtered out.
#   If this parameter is set to "None", no filtering will be done.
# - print_logs = whether there should be an output for the logs of the function. Like verbose.

# %% [markdown]
# Now we can inspect our results. First, the corrected output dataframe.

# %%
corrected_df

# %%
df

# %% [markdown]
# We can also look further into the correction_info object, to see the parameters that were chosen for each feature.
#
# For example, let's check the parameters used for the 200th feature.

# %%
correction_info[200]

# %% [markdown]
# ### Plot an example curve for one feature

# %% [markdown]
# We can also plot an example feature, to see how the values have changed and what the LOESS curve would look like for the data of this feature.

# %%
dc.loess_drift_correction.loess_example_curve(
    df=df,
    feature_idx=5,
    sample_cols=sample_cols,
    qc_cols=qc_cols,
    sample_order=sample_order,
)

# %% [markdown]
# In this plot, we can see all of the data points of this feature, ordered by the injection time. The red points are our QC samples, so we can see whether there was any instrumental drift over time. The LOESS curve is also calculated, with the smoothing value alpha chosen with leave-one-out cross validation, just like in the run_drift_correction function.
#
# Alternative values for the smoothing parameter alpha can be tested by adding the argument "alpha" and choosing a value, just like in the example below.

# %%
dc.loess_drift_correction.loess_example_curve(
    df=df,
    feature_idx=5,
    sample_cols=sample_cols,
    qc_cols=qc_cols,
    sample_order=sample_order,
    alpha=0.6,
)

# %% [markdown]
# ## 2. Common Principal Components Analysis-based drift correction

# %% [markdown]
# We can also use another method for correcting drift which is based on Common Principal Components Analysis (CPCA).
# This method is based on common principal components in defined groups of the data. It assumes that when calculating common principal components of QC samples, the drift contribution can be identified as the direction capturing maximum variance that simultaneously diagonalizes the covariance matrices of a set of classes.
#
# Therefore, the variability in the identified direction can be explained as caused by experimental drift and subtracted from all samples.
#
# Let's use different example data for demonstrating this method.

# %% [markdown]
# ### Load in data

# %%
# Load data
df = pd.read_csv(
    "../../example_data/DidacMauricio_hilic/DM_FIS2018_Hilic_pos_results2023_filled_imputed.csv"
)

# Define sample columns and qc columns
collist = list(df.columns.values)
sample_cols = []
qc_cols = []
for col in collist:
    if col.startswith("AAA"):
        sample_cols.append(col)
    elif col.startswith("QC"):
        qc_cols.append(col)

# %%
df

# %% [markdown]
# Seeing as this method is based on calculating principal components, as with PCA, there must not be any missing data.
#
# We will therefore first calculate missingness (NAs).
#
# We need to check for missing values in both the sample columns and the QC columns.

# %%
if dc.check_missingness(df, sample_cols + qc_cols):
    print(
        "There are missing values. Consider imputing first, or use LOESS drift correction instead."
    )
else:
    print("\nThere is no missingness. We can proceed with the CPCA drift correction.")

# %% [markdown]
# ### Visualise non-corrected data
#
# Now that we know we can proceed, let's visualise our data before drift correction with a PCA.

# %%
dc.pca_for_cpca_drift(
    df,
    sample_cols,  # list of col names, OR dict {group_name: [col names]}
    qc_cols,
    log_transform=True,
    title="PCA",
)

# %% [markdown]
# We see some clear indication of instrumental drift in the QC samples.

# %% [markdown]
# ### Run drift correction based on CPCA

# %%
df_corrected = dc.run_cpca_drift_correction(df, sample_cols, qc_cols, n_comps=1)

# %% [markdown]
# Let's plot the corrected data.

# %%
dc.pca_for_cpca_drift(
    df_corrected,
    sample_cols,  # list of col names, OR dict {group_name: [col names]}
    qc_cols,
    log_transform=True,
    title="PCA",
)

# %% [markdown]
# There is some change, but still a signficant amount of drift is clearly visible from the PCA plot.
#
# We can play around with the n_comps variable which decides the number of components.

# %%
df_corrected_2comps = dc.run_cpca_drift_correction(df, sample_cols, qc_cols, n_comps=2)

df_corrected_3comps = dc.run_cpca_drift_correction(df, sample_cols, qc_cols, n_comps=3)

df_corrected_4comps = dc.run_cpca_drift_correction(df, sample_cols, qc_cols, n_comps=4)

# %%
dc.pca_for_cpca_drift(
    df_corrected_2comps,
    sample_cols,  # list of col names, OR dict {group_name: [col names]}
    qc_cols,
    log_transform=True,
    title="PCA with 2 components",
)

dc.pca_for_cpca_drift(
    df_corrected_3comps,
    sample_cols,  # list of col names, OR dict {group_name: [col names]}
    qc_cols,
    log_transform=True,
    title="PCA with 3 components",
)

dc.pca_for_cpca_drift(
    df_corrected_4comps,
    sample_cols,  # list of col names, OR dict {group_name: [col names]}
    qc_cols,
    log_transform=True,
    title="PCA with 4 components",
)

# %% [markdown]
# The correction methods with three and four components are already looking better. Let's calculate the centroids of the QC principal components and the distance of the QC points to them, to objectively decide which number of n_comps is most favourable.

# %%
dc.cpca_centroid(df_corrected, sample_cols, qc_cols, log_transform=True)
dc.cpca_centroid(df_corrected_2comps, sample_cols, qc_cols, log_transform=True)
dc.cpca_centroid(df_corrected_3comps, sample_cols, qc_cols, log_transform=True)
dc.cpca_centroid(df_corrected_4comps, sample_cols, qc_cols, log_transform=True)

# %% [markdown]
# According to this, the CPCA method with three principal components is most favourable in this case.
#
# We can go ahead and continue our metabolomics data analysis with this data set.

# %%
df_corrected_3comps
