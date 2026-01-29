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
# # Differential regulation (ANOVA for more than two groups)
#
# This API example shows the functionality in the [`acore.differential_regulation`](acore.differential_regulation) module:
#
# - Analysis of Variance (ANOVA) as omnibus test across more than two groups
# - posthoc t-tests between groups
#
# An omnibus analysis across groups
# is combined with posthoc analysis between each set of the separate groups.
#
# The function is the same as for the two groups case. The `group1` and
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
viewed_cols = ["mean(group1)", "mean(group2)", "std(group1)", "std(group2)", "log2FC"]
viewed_cols.extend(view.columns)
view
