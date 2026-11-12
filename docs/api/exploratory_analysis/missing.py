# %% [markdown]
# # Missing values
#
# Using the Alzheimer's CSF proteomics dataset
# ([PXD016278](../../../example_data/alzheimer_proteomics/README.md)) as an example, we
# show how to find features with enough valid values, either by an absolute number of
# valid values or by a maximum fraction of missing values. Both can be applied across
# all samples or per group (here Alzheimer's disease vs. control).

# %%
# %pip install acore

# %%
import pandas as pd

import acore.exploratory_analysis as ea

# %% [markdown]
# ## Load proteomics example data
# The file contains the clinical data (including the group `AD`) and the
# protein group intensities for each sample.

# %%
data = pd.read_csv(
    "https://raw.githubusercontent.com/Multiomics-Analytics-Group/acore/"
    "refs/heads/main/"
    "example_data/alzheimer_proteomics/alzheimer_example_omics_and_clinic.csv",
    index_col=0,
)
clinic_cols = ["age", "male", "collection_site"]
group = "AD"  # Alzheimer's disease (1) vs. control (0)
data[group].value_counts()

# %% [markdown]
# The fraction of missing values per protein group is heterogeneous.

# %%
data.drop(columns=clinic_cols + [group]).isna().mean().describe()

# %% [markdown]
# ## Number of valid values
# Keep the features which have at least `min_valid` valid values across all samples.

# %%
keep = ea.missing.extract_number_missing(
    data, min_valid=170, drop_cols=clinic_cols, group=None
)
msg = (
    "{} features have at least {} valid values across"
    " all samples ( {} were dropped)."
)
print(
    msg.format(
        len(keep),
        170,
        data.shape[1] - len(clinic_cols) - 1 - len(keep),
    )
)

# %% [markdown]
# Per group, a feature is kept if at least one group has `min_valid` valid values.

# %%
keep = ea.missing.extract_number_missing(
    data, min_valid=80, drop_cols=clinic_cols, group=group
)
print(msg.format(len(keep), 170, data.shape[1] - len(clinic_cols) - 1 - len(keep)))


# %% [markdown]
# ## Percentage of missing values
# Keep the features with at most 10% missing values across all samples.

# %%
keep = ea.missing.extract_percentage_missing(
    data, missing_max=0.1, drop_cols=clinic_cols, group=None
)
msg = (
    "{} features have at most {:.1%} missing values across"
    " all samples ( {} were dropped)."
)
print(msg.format(len(keep), 0.1, data.shape[1] - len(clinic_cols) - 1 - len(keep)))


# %% [markdown]
# Per group, `how="all"` keeps a feature if it passes the threshold in at least one
# group, whereas `how="any"` requires that it passes in all groups.

# %%
for how in ["all", "any"]:
    keep = ea.missing.extract_percentage_missing(
        data, missing_max=0.1, drop_cols=clinic_cols, group=group, how=how
    )
    print(how, len(keep))

# %% [markdown]
# Done.
