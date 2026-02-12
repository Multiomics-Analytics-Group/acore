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
    "https://raw.githubusercontent.com/Multiomics-Analytics-Group/acore/"
    "updt_diff_reg_api_example/example_data/MTBLS13311/"
)
fname: str = "MTBLS13411_omics_and_meta.csv"  # combined omics and meta data
covariates: list[str] = []
group: str = "group"
drop_cols: list[str] = ["Factor Value[Strain type]"]
subject_col: str | int = 0
factor_and_covars: list[str] = [group, *covariates]

# %% [markdown]
# # ANOVA analysis for two groups
# Use combined dataset for ANOVA analysis.

# %% tags=["hide-input"]
omics_and_meta = (
    pd.read_csv(f"{BASE}/{fname}", index_col=subject_col)
    .convert_dtypes()
    .dropna(subset=factor_and_covars)
)
omics_and_meta

# %% [markdown]
# Drop unnecessary columns, if there are any specified in `drop_cols`.

# %% tags=["hide-input"]
if drop_cols:
    omics_and_meta.drop(columns=drop_cols, inplace=True)
omics_and_meta

# %% [markdown]
# Check data types of the columns. Metadata can be numeric, but also strings.

# %% tags=["hide-input"]
omics_and_meta.dtypes.value_counts()

# %% tags=["hide-input"]
omics_and_meta[factor_and_covars]


# %% [markdown]
# ## With four groups
# Acore make each combinatorial comparison between groups in the group column.


# %%
if isinstance(subject_col, int):
    subject_col = omics_and_meta.index.name
anova = (
    ad.run_anova(
        omics_and_meta,  # .reset_index(),
        subject=subject_col,
        drop_cols=[],
        group=group,
    ).set_index(["identifier", "group1", "group2"])
    # .sort_values(by="padj")
)
anova.head().T

# %% [markdown]
# ### pairwise t-test results:

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
view

# %%
