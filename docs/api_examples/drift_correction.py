# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.3
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
# ### Load in example data

# %%
# #%pip install acore

# %%
import acore
from acore import filter_metabolomics as fm

import pandas as pd
import os
import importlib

importlib.reload(acore)

# %%
from acore import drift_correction as dc

# %%
help(dc.loess_drift_correction)

# %%
help(dc.loess_drift_correction.run_drift_correction)

# %% [markdown]
# ## Load in data

# %%
df = pd.read_excel("../../example_data/aradopsis_seedling_lipids/kehelpannala_AnArabidopsisLipid_seedling.xlsx")
sample_order = pd.read_csv("../../example_data/aradopsis_seedling_lipids/seedling_art_sample_order.csv")

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
    "SD1": ['Sd1_1','Sd1_2', 'Sd1_3', 'Sd1_4', 'Sd1_5', 'Sd1_6'],
    "SD2": ['Sd2_1', 'Sd2_2', 'Sd2_3', 'Sd2_4', 'Sd2_5', 'Sd2_6', 'Sd2_7', 'Sd2_8', 'Sd2_9'],
    "SD3": ['Sd3-1', 'Sd3-2', 'Sd3-3', 'Sd3-4', 'Sd3-5', 'Sd3-6', 'Sd3-7', 'Sd3-8', 'Sd3-9'],
    "QC" : [ 'PBQC_Sd_1', 'PBQC_Sd_2', 'PBQC_Sd_3', 'PBQC_Sd_4', 'PBQC_Sd_5', 'PBQC_Sd_6']
}

# %%
# Separate groups
sample_cols = data_groups["SD1"] + data_groups["SD2"] + data_groups[ "SD3"]
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
    print_logs=True
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
