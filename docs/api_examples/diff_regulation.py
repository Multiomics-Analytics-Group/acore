# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: .venv
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Differential regulation
#
# This API example shows the functionality in the [`acore.differential_regulation`](acore.differential_regulation) module:
#
# - ANOVA (for two groups)
# - ANCOVA (for two groups)
# - comparision between two results
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

# %% [markdown]
# # ANOVA analysis for two groups
# not controlling for covariates
# > To check: pvalues for proteins with missing mean values? some merging issue?

# %%
if not isinstance(subject_col, str):
    subject_col = omics_and_clinic.index.name or "index"
    omics_and_clinic.rename_axis(subject_col, axis=0, inplace=True)
anova = (
    ad.run_anova(
        omics_and_clinic.reset_index(),
        subject=subject_col,
        drop_cols=covariates,
        group=group,
    )
    .set_index("identifier")
    .sort_values(by="padj")
)
anova

# %%
# ToDo: Something is wrong with the mean and std dev calculations for each group
# FC is also calculated along means and std dev.
anova.describe().T

# %% [markdown]
# Set subject to None

# %%
anova = (
    ad.run_anova(
        omics_and_clinic,
        subject=None,
        drop_cols=covariates,
        group=group,
    )
    .set_index("identifier")
    .sort_values(by="padj")
)
anova

# %% [markdown]
# view averages per protein group

# %% tags=["hide-input"]
view = anova.iloc[:, 2:7]
viewed_cols = view.columns.to_list()
view

# %% [markdown]
# Test results

# %% tags=["hide-input"]
regex_filter = "pval|padj|reject|stat|FC"
view = anova.filter(regex=regex_filter)
viewed_cols.extend(view.columns)
view

# %% [markdown]
# Other information

# %% tags=["hide-input"]
anova.drop(columns=viewed_cols)

# %% [markdown]
# # Comparing ANOVA and ANCOVA results for two groups
# Cross tabulated results after FDR correction for both ANOVA and ANCOVA

# %% tags=["hide-input"]
pd.crosstab(
    anova.rejected.rename("rejected ANOVA"),
    ancova_acore.rejected.rename("rejected ANCOVA"),
)

# %% [markdown]
# The ANOVA and ANCOVA results are not identical. Control for relevant covariates
# as they can confound the results. Here we used age and biological sex.

# %% [markdown]
# # With three and more groups
# Acore make each combinatorial comparison between groups in the group column.

# %%
CLINIC: str = "meta.csv"  # clincial data
meta = (
    pd.read_csv(f"{BASE}/{CLINIC}", index_col=0)
    .convert_dtypes()
    .rename(
        {
            "_collection site": "site",
            "_age at CSF collection": "age",
            "_gender": "gender",
        },
        axis=1,
    )
)[["site", "age", "gender"]].astype(
    {
        "gender": "category",
        "site": "category",
    }
)
meta

# %% [markdown]
# Sample five protein groups (for easier inspection) and combine with metadata.

# %%
omics_and_clinic = (
    omics_and_clinic.drop(columns=["AD", "age", "male"])
    .sample(5, axis=1, random_state=42)
    .join(meta)
)
omics_and_clinic

# %%
anova = (
    ad.run_anova(
        omics_and_clinic,  # .reset_index(),
        subject="Sample ID",
        drop_cols=["age", "gender"],
        group="site",
    ).set_index(["identifier", "group1", "group2"])
    # .sort_values(by="padj")
)
anova

# %% [markdown]
# pairwise t-test results:

# %%
cols_pairwise_ttest = [
    # "group1",
    # "group2",
    "mean(group1)",
    "std(group1)",
    "mean(group2)",
    "std(group2)",
    "posthoc Paired",
    "posthoc Parametric",
    "posthoc T-Statistics",
    "posthoc dof",
    "posthoc tail",
    "posthoc pvalue",
    "posthoc BF10",
    "posthoc effsize",
    # "identifier",
    "log2FC",
    "FC",
    "efftype",
]
anova[cols_pairwise_ttest]

# %% [markdown]
# ANOVA results

# %%
anova.drop(columns=cols_pairwise_ttest)

# %% [markdown]
# Test results

# %% tags=["hide-input"]
regex_filter = "pval|padj|reject|stat|FC"
view = anova.filter(regex=regex_filter)
viewed_cols.extend(view.columns)
view

# %% [markdown]
# Done.
