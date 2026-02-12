# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: all
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: .venv
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Differential regulation: T-test for two groups
#
# This API example shows the functionality in the [`acore.differential_regulation`](acore.differential_regulation) module:
#
# - ANOVA (for two groups) is equivalent to t-test between two groups
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
import vuecore.plots.basic.scatter

import acore.differential_regulation as ad

dsp_pandas.format.set_pandas_options(
    max_columns=9,
    max_colwidth=20,
)

# %% tags=["parameters"]
BASE = (
    "https://raw.githubusercontent.com/Multiomics-Analytics-Group/acore/"
    "updt_diff_reg_api_example/example_data/alzheimer_proteomics/"
)
# data is already preprocessed: log2, filtered
fname: str = "alzheimer_example_omics_and_clinic.csv"  # combined omics and meta data
covariates: list[str] = ["age", "male"]
group: str = "AD"
subject_col: str = "Sample ID"
drop_cols: list[str] = []
factor_and_covars: list[str] = [group, *covariates]

# %% [markdown]
# Alternatively you can use the metabolomics dataset MTBLS13311,
# which is also preprocessed and combined.

# %% tags=["hide-input"]
# BASE = (
#     "https://raw.githubusercontent.com/Multiomics-Analytics-Group/acore/"
#     "updt_diff_reg_api_example/example_data/MTBLS13311/"
# )
# fname: str = "MTBLS13411_omics_and_meta.csv"  # combined omics and meta data
# covariates: list[str] = []
# group: str = "Factor Value[Strain type]"
# drop_cols: list[str] = ["group"]
# subject_col: str | int = 0
# factor_and_covars: list[str] = [group, *covariates]

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
omics_and_meta[[group, *covariates]]


# %% [markdown]
# # ANOVA analysis for two groups
# - is not controlling for covariates

# %%
if not isinstance(subject_col, str):
    subject_col = omics_and_meta.index.name or "index"
    omics_and_meta.rename_axis(subject_col, axis=0, inplace=True)
anova = (
    ad.run_anova(
        omics_and_meta,
        subject=None,
        drop_cols=covariates,
        group=group,
    )
    .set_index("identifier")
    .sort_values(by="padj")
)
anova

# %%
# See summary statistics for numeric values in ANOVA results, such as for
# - pvalue (unadjusted)
# - padj (adjusted p-value, FDR corrected)
# - fold change (FC) of the mean of group1 vs group2 (group1/group2)
anova.describe().T

# %% [markdown]
# ## View Volcano plot for the ANOVA results
# - p-value (unadjusted) using the -log10 transformation on the y-axis
# - fold change (FC) in the log2 space on the x-axis (log2FC of 2 means a fold
#   change of 4, log2FC of -1 means a fold change of 0.5)
# - significant features are colored in red using the `reject` column,
#   which is True for significant features and False for non-significant features, and
#   is calculated based on the `padj` column and a significance threshold (default 0.05)

# %%
scatter_plot_adv = vuecore.plots.basic.scatter.create_scatter_plot(
    data=anova.reset_index(),
    x="log2FC",
    y="-log10 pvalue",
    color="rejected",
    title="Simple Volcano Plot",
    subtitle="Visualizing ANOVA results",
    labels={
        "log2FC": "Log2 Fold Change",
        "-log10 pvalue": "-log10(p-value)",
        "rejected": "FDR corrected Significant",
        "identifier": "Protein Identifier",
    },
    hover_data=["identifier"],
    # currently does not work:
    # color_discrete_map={False: "#2166AC", True: "#B2182B"},  # Blue  # Red
    color_discrete_sequence=["red", "blue"],
    opacity=1,
    marker_line_width=1,
    marker_line_color="darkgray",
    width=800,
    height=600,
)
scatter_plot_adv

# %% [markdown]
# ## Inspect ANOVA results
# - summary statistics for each group
# - test results
# - t-test statistic, FDR correction method, etc.

# %% [markdown]
# ## View averages per protein group

# %% tags=["hide-input"]
view = anova.iloc[:, 2:7]
viewed_cols = view.columns.to_list()
view

# %% [markdown]
# ### Test results

# %% tags=["hide-input"]
regex_filter = "pval|padj|reject|stat|FC"
view = anova.filter(regex=regex_filter)
viewed_cols.extend(view.columns)
view

# %% [markdown]
# ### Other information

# %% tags=["hide-input"]
anova.drop(columns=viewed_cols)

# %% [markdown]
# Done.
