# %% [markdown]
# # Process script for MTBLS13311 example dataset

import numpy as np
import pandas as pd

# %%
fname_meta = "s_MTBLS13411.txt"
fname_data = "m_MTBLS13411_LC-MS_alternating_hilic_metabolite_profiling_v2_maf.tsv"

# %%
metadata = pd.read_csv(fname_meta, index_col=0, sep="\t")
metadata

# %%
data = pd.read_csv(fname_data, sep="\t", index_col="metabolite_identification")
data

# %%
df = data[metadata.index].T
df = np.log2(df)
df

# %%
df.to_csv("MTBLS13411_processed_data.csv")

# %%
factor = "Factor Value[Strain type]"
group = "group"

metadata[group] = metadata.index.str.split("-").str[0]
metadata[[factor, group]]

# %%
metadata[[factor, group]].to_csv("MTBLS13411_meta_data.csv")

# %%
omics_and_meta = metadata[[factor, group]].join(df, how="inner")
omics_and_meta

# %%
omics_and_meta.to_csv("MTBLS13411_omics_and_meta.csv")

# %%
