# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: all
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.3
# ---

# %% [markdown]
# # Differential regulation: ANOVA across more than two groups and posthoc t-tests
#
# This API example shows the functionality in the
# [`acore.differential_regulation`](acore.differential_regulation) module:
#
# - Analysis of Variance (ANOVA) as omnibus test across more than two groups
# - posthoc t-tests between each of the groups
#
# An omnibus analysis across groups
# is combined with posthoc analysis between each set of the separate groups.
#
# The function is the same as for the two groups case. The `group1` and
# `group2` columns give the posthoc comparison.
#

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
    "main/example_data/MTBLS13311/"
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
        subject=None,
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
    "posthoc pvalue adj",
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

# %% [markdown]
# ## Volcano plot of ANOVA results
# - volcano plot of ANOVA
# - more than one group present in the example.
# - double check what adj. pvalue `rejected` is based on


# %%
# ToDo: check why the pvalue is not used for log-tranformation of omnibus test
# anova["-log10 pvalue #2"] = -np.log10(anova["pvalue"]) # not current -log10 pvalue

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
    hover_data=["identifier", "group1", "group2"],
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

# %%
