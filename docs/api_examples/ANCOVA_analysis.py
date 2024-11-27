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
#
# - [ ] include a PCA colored by groups as well as covariance factors

# %% tags=["hide-output"]
# %pip install acore

# %% tags=["hide-input"]
from pathlib import Path

import dsp_pandas
import numpy as np
import pandas as pd

import acore.differential_regulation as ad

dsp_pandas.format.set_pandas_options(
    max_columns=9,
    max_colwidth=20,
)

# %% tags=["parameters"]
BASE = (
    "https://raw.githubusercontent.com/RasmussenLab/njab/"
    "HEAD/docs/tutorial/data/alzheimer/"
)
CLINIC: str = "clinic_ml.csv"  # clincial data
OMICS: str = "proteome.csv"  # omics data
freq_cutoff: float = (
    0.7  # at least x percent of samples must have a value for a feature (here: protein group)
)

# %% [markdown]
# ## Load data.
# Clinical data:

# %% tags=["hide-input"]
clinic = pd.read_csv(f"{BASE}/{CLINIC}", index_col=0).convert_dtypes()
omics = pd.read_csv(f"{BASE}/{OMICS}", index_col=0)
clinic

# %% [markdown]
# Proteomics data:

# %% tags=["hide-input"]
omics

# %% [markdown]
# ### Filtering data

# %% [markdown]
# If data is already filtered and/or imputed, skip this step.

# %% [markdown]
# ## Filtering features

# %% tags=["hide-input"]
M_before = omics.shape[1]
omics = omics.dropna(thresh=int(len(omics) * freq_cutoff), axis=1)
M_after = omics.shape[1]
msg = (
    f"Removed {M_before-M_after} features "
    f"with more than {(1-freq_cutoff)*100:.2f}% missing values."
    f"\nRemaining features: {M_after} (of {M_before})"
)
print(msg)
# keep a map of all proteins in protein group, but only display first protein
# proteins are unique to protein groups
pg_map = {k: k.split(";")[0] for k in omics.columns}
omics = omics.rename(columns=pg_map)
# log2 transform raw intensity data:
omics = np.log2(omics + 1)
omics

# %% [markdown]
# Consider replacing with the filter from the acore package!

# %% [markdown]
# ## Preparing metadata
# add both relevant clinical information to the omics data

# %% tags=["hide-input"]
clinic[["age", "male", "AD"]].describe()

# %% tags=["hide-input"]
omics_and_clinic = omics.join(clinic[["age", "male", "AD"]])
omics_and_clinic


# %% [markdown]
# ### Checking missing data
# ... between two AD groups (after previous filtering)

# %% tags=["hide-input"]
data_completeness = (
    omics_and_clinic.groupby("AD").count().divide(clinic["AD"].value_counts(), axis=0)
)
data_completeness

# %% [markdown]
# Plot number of missing values per group, ordered by proportion of non-misisng values
# in non-Alzheimer disease group

# %% tags=["hide-input"]
data_completeness.T.sort_values(0).plot(style=".", ylim=(0, 1.05), alpha=0.5, rot=45)

# %% [markdown]
# Plot 20 protein groups with biggest difference in missing values between groups

# %% tags=["hide-input"]
idx_largerst_diff = (
    data_completeness.diff().dropna().T.squeeze().abs().nlargest(20).index
)
ax = (
    data_completeness.loc[:, idx_largerst_diff]
    .T.sort_values(0)
    .plot(
        style=".",
        ylim=(0, 1.05),
        alpha=0.5,
        rot=45,
    )
)
_ = ax.set_xticks(range(len(idx_largerst_diff)))
_ = ax.set_xticklabels(
    idx_largerst_diff,
    rotation=45,
    ha="right",
    fontsize=7,
)

# %% [markdown]
# ## Running ANCOVA analysis
# Use combined dataset for ANCOVA analysis.

# %% tags=["hide-input"]
omics_and_clinic

# %% [markdown]
# metadata here is of type integer. All floats are proteomics measurements.

# %% tags=["hide-input"]
omics_and_clinic.dtypes.value_counts()

# %% tags=["hide-input"]
group = "AD"
covariates = ["male", "age"]
omics_and_clinic[[group, *covariates]]


# %% [markdown]
# run ANCOVA analysis

# %%
# omics_and_clinic = omics_and_clinic.astype(float)
# # ? this is no needed for run_ancova (the regex where groups are joined)
ancova = (
    ad.run_ancova(
        omics_and_clinic.astype({"AD": str}),  # ! target needs to be of type str
        # subject='Sample ID', # not used
        drop_cols=[],
        group="AD",  # needs to be a string
        covariates=covariates,
    )
    .set_index("identifier")
    .sort_values(by="posthoc padj")
)  # need to be floats?
ancova_acore = ancova
ancova

# %% [markdown]
# The first columns contain group averages for each group for the specific
# protein group
# ancova.iloc[:, :6]

# %% [markdown]
# The others contain the test results (based on a linear model) for each protein group
# (on each row). Some information is duplicated.

# %% tags=["hide-input"]
regex_filter = "pval|padj|reject|post"
ancova.filter(regex=regex_filter)

# %% [markdown]
# The other information is about fold-changes and other information.

# %% tags=["hide-input"]
ancova.iloc[:, 6:].filter(regex=f"^(?!.*({regex_filter})).*$")

# %% [markdown]
# ## ANOVA analysis
# not controlling for covariates
# > To check: pvalues for proteins with missing mean values? some merging issue?

# %%
anova = (
    ad.run_anova(
        omics_and_clinic.reset_index(),
        subject="Sample ID",
        drop_cols=covariates,
        group="AD",
    )
    .set_index("identifier")
    .sort_values(by="padj")
)
anova

# %% [markdown]
# view averages per protein group

# %% tags=["hide-input"]
anova.iloc[:, 2:7]

# %% [markdown]
# ## Comparing ANOVA and ANCOVA results
# Cross tabulated results after FDR correction for both ANOVA and ANCOVA

# %% tags=["hide-input"]
pd.crosstab(
    anova.rejected.rename("rejected ANOVA"),
    ancova_acore.rejected.rename("rejected ANCOVA"),
)

# %% [markdown]
# The ANOVA and ANCOVA results are not identical. Control for relevant covariates
# as they can confound the results. Here we used age and biological sex.
