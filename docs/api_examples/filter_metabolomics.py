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
#     display_name: acore
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Metabolomics data filtering

# %%
# %pip install acore

# %%
import acore
from acore import filter_metabolomics as fm

import pandas as pd
import os
import importlib

importlib.reload(acore)

# %% [markdown]
# Load in your data. We will use an example data set from MetaboLights. It can be found in example_data/MTBLS733.

# %%
data_path = "../../example_data/MTBLS733/MetaboLights-MTBLS733-Nextflow4MS-DIAL.csv"
data = pd.read_csv(data_path)

# %% [markdown]
# Let's look more into the data.
#
# - The .dtypes function shows which columns are numeric (int64, float64) and which are categorical (object, bool).
#
# - The .describe function summarises the numeric columns.
#
# - With .columns, we can see which columns our data has.

# %%
data.dtypes

# %%
data.describe()

# %%
print(f"There are {data.shape[0]} rows and {data.shape[1]} columns in our data.")
print("Our data has the following columns:")
for colname in data.columns.values.tolist():
    print("\t",colname)


# %% [markdown]
# It looks like the m/z and RT columns contain categorical data, so we need to change that first before we can filter.

# %%
numeric_data = fm.make_numeric.convert_to_numeric(data, ["Average Mz", "Average Rt(min)"], print_na_summary=True)


# %%
numeric_data

# %%
numeric_data.dtypes

# %% [markdown]
# Now the Mz and RT columns are dtype float, so they are numeric. That means that now we can proceed with filtering.
#
# We want to filter based on RT, by removing all features that have a RT below a certain time, to exclude features that are at the dead volume. In our case, the cut-off will be 0.8 minutes.
#
# We also want to filter out features that have an m/z value below 600 an dhave m/z decimals between 0.3 and 0.9.
#
# We can do both of those things with the filtering function.

# %%
filtered_data, removed_features = fm.filter_mz_rt(
    numeric_data, 
    "Average Rt(min)", 
    "Average Mz", 
    mz_decimals=(0.3, 0.9), 
    mz_low=600,
    rt_dead_volume=0.8,
    save_removed = True)



# %% [markdown]
# Let's look at our filtered data.

# %%
print(f"There are {data.shape[0]} rows and {data.shape[1]} columns in our original data.")
print(f"There are {filtered_data.shape[0]} rows and {filtered_data.shape[1]} columns in our filtered data.")
print(f"{removed_features.shape[0]} features were removed from our data.")

# %% [markdown]
# We can also look into our removed features.

# %%
bdv = removed_features[removed_features["RemovalReason"] == "BelowDeadVolume"].shape
nbr = removed_features[removed_features["RemovalReason"] == "NotBiologicallyRelevant"].shape
print(f"{bdv[0]} features were removed because they were below the dead volume.")
print(f"{nbr[0]} features were removed because they were deemed not biologically relevant (m/z-based filtering).")

removed_features
