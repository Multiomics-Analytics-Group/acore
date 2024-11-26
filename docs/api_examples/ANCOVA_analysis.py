# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: .venv
#     language: python
#     name: python3
# ---

# %% [markdown]
# # ANCOVA analysis

# %%
# include a PCA colored by groups as well as covariance factors 
# This is a new update

# %% [markdown]
# Import data.

# %%
import pandas as pd
import acore.differential_regulation as ad
from pathlib import Path
import numpy as np

folder_out = Path("data")

# %%
CLINIC: str = 'https://raw.githubusercontent.com/RasmussenLab/njab/HEAD/docs/tutorial/data/alzheimer/clinic_ml.csv'  # clincial data
OMICS: str = 'https://raw.githubusercontent.com/RasmussenLab/njab/HEAD/docs/tutorial/data/alzheimer/proteome.csv'  # omics data

# %%
clinic = pd.read_csv(CLINIC, index_col=0).convert_dtypes()
omics = pd.read_csv(OMICS, index_col=0)

# %%
clinic

# %%
omics

# %% [markdown]
# ### Filtering data

# %% [markdown]
# If data is already filtered and/or imputed, skip this step. 

# %%
# Filtering parameters
freq_cutoff = 0.7


# %%
M_before = omics.shape[1]
omics = omics.dropna(thresh=int(len(omics) * freq_cutoff), axis=1)
M_after = omics.shape[1]
msg = (
    f"Removed {M_before-M_after} features with more than {freq_cutoff*100}% missing values." # if theres 100 feat with >30% missing, how can there be 400 feat with >70% 
    f"\nRemaining features: {M_after} (of {M_before})")
print(msg)
# keep a map of all proteins in protein group, but only display first protein
# proteins are unique to protein groups
pg_map = {k: k.split(";")[0] for k in omics.columns}
omics = omics.rename(columns=pg_map)
# log2 transform raw intensity data:
omics = np.log2(omics + 1)
omics

# %%

# %% [markdown]
# Consider replacing with the filter from the acore package!

# %% [markdown]
# ### Preparing metadata

# %%
clinic['age'].info()

# %%
clinic_omics = omics.join(clinic)
clinic_omics

# %%
omics_group = clinic_omics.drop(columns = ['Kiel','Magdeburg','Sweden','male','age'])

# %%
omics_group

# %%

# %% [markdown]
# ### Checking missing data

# %%
data_completeness = omics_group.groupby("AD").count().divide(clinic['AD'].value_counts(), axis=0)
data_completeness

# %%
data_completeness.T.sort_values(0).plot(style='.',ylim=(0,1))

# %% [markdown]
# ### Running ANCOVA analysis

# %%
clinic_omics

# %%
clinic_omics.index.to_series().info() 

# %%
clinic_omics.dtypes.value_counts()

# %%
col='A0A024QZX5'
group="AD"
covariates=['male',]
clinic_omics[[group, col] + covariates]

# %%
ad.calculate_ancova(clinic_omics.astype('float'), column='A0A024QZX5', group="AD", covariates=['male',])

# %%
clinic_omics

# %%
clinic_omics = clinic_omics.astype(float) # this is no needed for run_ancova (the regex where groups are joined)
ad.calculate_ancova(clinic_omics, column='A0A024QZX5', group="AD", covariates=['male',])

# %%

ancova = ad.run_ancova(
                        clinic_omics.astype({'AD':str}),
                        # subject='Sample ID', # not used
                        drop_cols=['Kiel','Magdeburg','Sweden','age',
                                   ],
                        group='AD', # needs to be a string 
                        covariates=['male',]) # need to be floats?

# %%
anova = ad.run_anova(clinic_omics.reset_index(),
                        subject='Sample ID',
                        drop_cols=['Kiel','Magdeburg','Sweden','age','male'],
                        group='AD')

# %% [markdown]
# ### Running ANOVA analysis (optional)

# %%

# %%

# %% [markdown]
# ### Comparing ANOVA and ANCOVA results

# %%
