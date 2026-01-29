# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: all
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
# ---

# %% [markdown]
# # Differential regulation (ANCOVA for two groups)
#
# This API example shows the functionality in the
# [`acore.differential_regulation`](acore.differential_regulation) module:
#
# - Analysis of Covariance (ANCOVA) for two groups
#
# Then we can do the same for examples with three
# and more groups, where a omnibus analysis across groups
# is combined with posthoc anlysis for separate groups.
#
# The functions are the same for both cases. The `group1` and
# `group2` columns give the posthoc comparison.
#
# Using vuecore we add:
# - [ ] include a PCA colored by groups as well as covariance factors

# %% tags=["hide-output"]
# %pip install acore

# %% tags=["hide-input"]
import dsp_pandas
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
CLINIC_ML: str = "clinic_ml.csv"  # clinical data
OMICS: str = "proteome.csv"  # omics data
freq_cutoff: float = (
    0.7  # at least x percent of samples must have a value for a feature (here: protein group)
)
#
covariates: list[str] = ["age", "male"]
group: str = "AD"
subject_col: str = "Sample ID"
factor_and_covars: list[str] = [group, *covariates]

# BASE = (
#     "https://raw.githubusercontent.com/Multiomics-Analytics-Group/acore/"
#     "HEAD/example_data/MTBLS13311/"
#     ""
# )
# CLINIC_ML: str = "MTBLS13411_meta_data.csv"  # clinical data
# OMICS: str = "MTBLS13411_processed_data.csv"  # omics data
# covariates: list[str] = []
# group: str = "Factor Value[Strain type]"
# subject_col: str | int = 0
# factor_and_covars: list[str] = [group, *covariates]

# %% [markdown]
# # ANCOVA analysis for two groups
# Use combined dataset for ANCOVA analysis.

# %% tags=["hide-input"]
omics_and_clinic = pd.read_csv(
    "../../example_data/alzheimer_proteomics/alzheimer_example_omics_and_clinic.csv",
    index_col=subject_col,
)
omics_and_clinic

# %% [markdown]
# metadata here is of type integer. All floats are proteomics measurements.

# %% tags=["hide-input"]
omics_and_clinic.dtypes.value_counts()

# %% tags=["hide-input"]
omics_and_clinic[[group, *covariates]]

# %% [markdown]
# run ANCOVA analysis

# %%
# omics_and_clinic = omics_and_clinic.astype(float)
# # ? this is no needed for run_ancova (the regex where groups are joined)
ancova = (
    ad.run_ancova(
        omics_and_clinic.astype({group: str}),  # ! target needs to be of type str
        # subject=subject_col, # not used
        drop_cols=[],
        group=group,  # needs to be a string
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

# %%
ancova.iloc[:, :6]

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
