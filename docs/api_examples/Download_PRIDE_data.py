# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Download and Process PRIDE Data
# This notebook shows how `acore` can be used to download data from
# the Proteomics Identifications Database - PRIDE -
# ([ebi.ac.uk/pride/](https://www.ebi.ac.uk/pride/))
# and parse the data to be used in the analytics core.
# and quickly formated to start analyzing them with the functionality in the analytics core.
#
# > based on CKG recipe: [Download PRIDE Data](https://ckg.readthedocs.io/en/latest/notebooks/recipes/Download_PRIDE_data.html)

# %%
from pathlib import Path

import numpy as np
import pandas as pd

import acore.io

# %% [markdown]
# ## Parameters
# Specify the PRIDE identifier and file to be downloaded
# - and where to store intermediate files.

# %% tags=["parameters"]
pxd_id: str = "PXD008541"  # PRIDE identifier
fname = "SearchEngineResults_secretome.zip.rar"  # file to download
folder_downloads = Path("downloaded")  # folder to download the file
folder_unzipped = Path("unzipped")  # folder to uncompress the file

# %% [markdown]
# ## Specify the PRIDE identifier and file to be downloaded

# %% [markdown]
# ## Download data
#
# We can use functionality in `acore` to directly download data files from EBI's
# PRIDE database [ebi.ac.uk/pride/](https://www.ebi.ac.uk/pride/).
# For that you just need to specify the
# PRIDE identifier for the project (`PXD_...`) and the name of the file to download.
# In this case, the project identifier is `PXD008541` and the file we will use
# is `SearchEngineResults_secretome.zip.rar`,
# a RAR compressed file with the output files from MaxQuant.

# %%
ret = acore.io.download_PRIDE_data(pxd_id=pxd_id, file_name=fname, to=folder_downloads)
ret["acore_downloaded_file"]  # folder_downloads / fname

# %% d [markdown]
# ## Read Data In

# %% [markdown]
# ### Decompress File

# %%
# ! you need a system installation of a rar archive tool
acore.io.unrar(filepath=ret["acore_downloaded_file"], to=folder_unzipped)

# %% [markdown]
# The list of files within the compressed folder

# %%
list(folder_unzipped.iterdir())

# %% [markdown]
# We use the proteinGroups file that contains the proteomics data processed
# using MaxQuant software.

# %%
proteinGroups_file = folder_unzipped / "proteinGroups.txt"

# %%
index_cols = [
    "Majority protein IDs",
]
pgs = pd.read_csv(proteinGroups_file, index_col=index_cols, sep="\t")
pgs.sample(5)

# %% [markdown]
# Get ride of potential contaminants, reverse (decoys) and identified only by a
# modification site
# reference:
# [cox-labs.github.io/coxdocs/output_tables.html#protein-groups](https://cox-labs.github.io/coxdocs/output_tables.html#protein-groups)

# %%
filters = ["Reverse", "Only identified by site", "Contaminant"]
pgs[filters].describe()

# %%
mask = pgs[filters].isna().all(axis=1)
pgs = pgs.loc[mask]
pgs

# %% [markdown]
# Then we can filter the columns that contain the string `LFQ intensity`. The sample names
# are part of the column names (here: `LFQ intensity {sample_name}`)

# %%
stub_intensity = "LFQ intensity"
pgs = pgs.filter(like=stub_intensity)
pgs

# %% [markdown]
# No we can get rid of the common part `LFQ intensity` and keep only the sample names

# %%
pgs.columns = pgs.columns.str.replace(stub_intensity, "").str.strip()
pgs.columns.name = "sample"
pgs

# %% [markdown]
# The group could be defined in a sample metadata file, but here we just parse it from the
# sample names omitting the numbers at the end of the sample name.

# %%
pgs.columns.str.replace(r"\d", "", regex=True)


# %% [markdown]
# We add to the information as a MultiIndex of group and sample name to the columns
# (sample metadata)
# pgs.columns = pd.MultiIndex.from_arrays(
#     [pgs.columns.str.replace(r"\d", "", regex=True), pgs.columns],
#     names=["group", pgs.columns.name],
# )
# pgs

# %% [markdown]
# From here we can stack both levels, name the values intensity. If we reset the index we
# get the original CKG format.

# %%
pgs = pgs.stack([0, 1]).to_frame("intensity")
pgs

# %% [markdown]
# First we  `log2` transform the data. We first set the zeros to `np.nan` to avoid
# `-inf` values.

# %%
pgs = np.log2(pgs.replace(0.0, np.nan).dropna())
pgs

# %%
pgs.reset_index()

# %% [markdown]
# Data to be saved in the CKG format.
