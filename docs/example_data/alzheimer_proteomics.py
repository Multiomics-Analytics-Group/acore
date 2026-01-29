# %% [markdown]
# # Alzheimer's MS-based proteomics
#
# Proteome Profiling in Cerebrospinal Fluid Reveals Novel Biomarkers of Alzheimer's Disease
#
# - [PXD016278](https://www.ebi.ac.uk/pride/archive/projects/PXD016278)
# - [publication](https://www.embopress.org/doi/full/10.15252/msb.20199356)
# - [curated data version from omiclearn](https://github.com/MannLabs/OmicLearn/tree/master/omiclearn/data)
# - download `Alzheimer.xlsx` from repository and process using `data/prepare_alzheimer_excel.py`
# as provided by njab, see
# [njab](https://github.com/RasmussenLab/njab/tree/main/docs/tutorial/data)

# %%
import numpy as np
import pandas as pd

import acore
import acore.types

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
# # Data
# ## Clinical data:

# %% tags=["hide-input"]
clinic = pd.read_csv(f"{BASE}/{CLINIC_ML}", index_col=subject_col).convert_dtypes()
omics = pd.read_csv(f"{BASE}/{OMICS}", index_col=subject_col)
clinic

# %% [markdown]
# ## Proteomics data:

# %% tags=["hide-input"]
omics

# %% [markdown]
# # Filtering data

# %% [markdown]
# If data is already filtered and/or imputed, skip this step.

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
# Check if all values are numeric as this is required for differential analysis

# %%
acore.types.check_numeric_dataframe(omics)

# %% [markdown]
# Validate the schema of the omics DataFrame. Builds and then uses the schema on the
# same data frame (experimental)

# %%
acore.types.build_schema_all_floats(omics).validate(omics)

# %% [markdown]
# For easier inspection we just sample 100 protein groups. Remove this step in a
# real analysis.

# %% tags=["hide-input"]
omics = omics.sample(min(omics.shape[1], 100), axis=1, random_state=42)
omics

# %% [markdown]
# Consider replacing with the filter from the acore package!

# %% [markdown]
# ## Preparing metadata
# add both relevant clinical information to the omics data

# %% tags=["hide-input"]
clinic[factor_and_covars].describe()

# %% tags=["hide-input"]
omics_and_clinic = clinic[factor_and_covars].dropna().join(omics)
omics_and_clinic

# %% [markdown]
# Check that the added clinical metadata is numeric

# %%
acore.types.check_numeric_dataframe(omics_and_clinic)

# %% [markdown]
# ## Checking missing data
# ... between two AD groups (after previous filtering)

# %% tags=["hide-input"]
data_completeness = (
    omics_and_clinic.groupby(by=group)
    .count()
    .divide(clinic[group].value_counts(), axis=0)
)
data_completeness

# %% [markdown]
# Plot number of missing values per group, ordered by proportion of non-misisng values
# in non-Alzheimer disease group

# %% tags=["hide-input"]
sort_by = data_completeness.index[0]
ax = data_completeness.T.sort_values(sort_by).plot(
    style=".", ylim=(0, 1.05), alpha=0.5, rot=45
)

# %% [markdown]
# Plot 20 protein groups with biggest difference in missing values between groups

# %% tags=["hide-input"]
idx_largerst_diff = (
    data_completeness.diff().dropna().T.squeeze().abs().nlargest(20).index
)
ax = (
    data_completeness.loc[:, idx_largerst_diff]
    .T.sort_values(sort_by)
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
# # Save data for use in ANCOVA and ANOVA examples
#

# %%
omics_and_clinic.to_csv(
    "../../example_data/alzheimer_proteomics/alzheimer_example_omics_and_clinic.csv"
)
