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
# # Differential regulation: ANCOVA for two groups
#
# This API example shows the functionality in the
# [`acore.differential_regulation`](acore.differential_regulation) module:
#
# - Analysis of Covariance (ANCOVA) for two groups
# - Effect of covariates is taken into account in the linear model
#
# Then we can do the same for examples with three
# and more groups, where a omnibus analysis across groups
# is combined with posthoc anlysis for separate groups.
#
# The functions are the same for both cases. The `group1` and
# `group2` columns give the posthoc comparison.
#

# %% tags=["hide-output"]
# %pip install acore

# %% tags=["hide-input"]
import dsp_pandas
import pandas as pd
import sklearn.impute
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
# ## Data loading
# Use combined dataset for ANCOVA analysis.

# %% tags=["hide-input"]
omics_and_meta = (
    pd.read_csv(f"{BASE}/{fname}", index_col=subject_col)
    .convert_dtypes()
    .dropna(subset=factor_and_covars)
)
omics_and_meta

# %% [markdown]
# Metadata here is of type integer. All floats are proteomics measurements.

# %% tags=["hide-input"]
omics_and_meta.dtypes.value_counts()

# %% tags=["hide-input"]
omics_and_meta[[group, *covariates]]


# %% [markdown]
# ## Impute missing values
# Currently the `run_ancova` function does not handle missing values for measurments
# which are not given.
# > Request implementation when necessary.
# We will therefore use acore imputation to impute missing values for the proteomics
# measurements.


# %%
# import acore.imputation_analysis

# acore.imputation_analysis.imputation_mixed_norm_KNN(
#     omics_and_meta,
#     index_cols=[group, *covariates],
#     group=group,
# )
omics_and_meta = pd.DataFrame(
    sklearn.impute.KNNImputer(n_neighbors=3).fit_transform(omics_and_meta),
    index=omics_and_meta.index,
    columns=omics_and_meta.columns,
)
omics_and_meta

# %% [markdown]
# ## Run ANCOVA analysis

# %%
# omics_and_meta = omics_and_meta.astype(float)
# # ? this is no needed for run_ancova (the regex where groups are joined)
ancova = (
    ad.run_ancova(
        omics_and_meta.astype({group: str}),  # ! target needs to be of type str
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
# ## Inspect ANCOVA results
# - summary statistics for each group
# - test results
# - t-test statistic, FDR correction method, etc.

# %% [markdown]
# ### View averages per protein group
# The first columns contain group averages for each group for the specific
# protein group

# %%
ancova.iloc[:, :6]

# %% [markdown]
# ### Test results
# The others contain the test results (based on a linear model) for each protein group
# (on each row). Some information is duplicated.

# %% tags=["hide-input"]
regex_filter = "pval|padj|reject|post"
ancova.filter(regex=regex_filter)

# %% [markdown]
# ### Other information

# %% tags=["hide-input"]
ancova.iloc[:, 6:].filter(regex=f"^(?!.*({regex_filter})).*$")


# %% [markdown]
# ## See Volcano plot of ANCOVA results

# %%
scatter_plot_adv = vuecore.plots.basic.scatter.create_scatter_plot(
    data=ancova.reset_index(),
    x="log2FC",
    y="-log10 pvalue",
    color="rejected",
    title="Simple Volcano Plot",
    subtitle="Visualizing ANCOVA results",
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
