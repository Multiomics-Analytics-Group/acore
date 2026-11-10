# %% [markdown]
# # Correlation analysis
#
# Using a metabolomics dataset ([MTBLS13411](../../../example_data/MTBLS13411.md)) as an
# example, we show how to compute and filter correlations between features.

# %%
# %pip install acore

# %%
import matplotlib.pyplot as plt
import pandas as pd
import scipy
import seaborn as sns

import acore.correlation_analysis as ca

# %% [markdown]
# ## Load metabolomics example data

# %%
data = (
    "https://raw.githubusercontent.com/Multiomics-Analytics-Group/acore/"
    "refs/heads/main/"
    "example_data/MTBLS13311/MTBLS13411_processed_data.csv"
)
data = pd.read_csv(data, index_col=0)
# specific to this data, we shorten some column names for better readability
data.columns = data.columns.str.split("(").str[-1].str.replace(")", "")
data

# %% [markdown]
# We add the group here based on the sample names. Alternatively you could merge it from
# the avilable metadata.

# %%
data["group"] = data.index.str.split("-").str[0]
data["group"].value_counts()

# %% [markdown]
# ## Correlation analysis
# See [`acore.correlation_analysis`](acore.correlation_analysis) for more functions
# and details.
#
# The basic functionality is built into pandas, but you need to filter out columns
# which are not numeric for pearson correlation.
#
# Generally: Ordered categorical values can be used, assuming equal spacing between
# the categories. Otherwise, continous numeric values are required.

# %%
corr = data.drop(columns=["group"]).corr(method="pearson")
corr

# %% [markdown]
# Plot the correlation heatmap using seaborn

# %%
plt.rcParams["xtick.labelsize"], plt.rcParams["ytick.labelsize"] = 5, 5
fig, ax = plt.subplots(figsize=(7.1, 6))
heatmap = sns.heatmap(
    corr,
    cmap="vlag",
    center=0,
    square=True,
    linewidths=0.1,
    cbar_kws={"label": "Pearson r"},
    ax=ax,
)
ax.set(title="Correlation Heatmap")
fig.tight_layout()

# %%
# If you only want to keep the lower triangle of the correlation matrix to have
# unique values of interst, you can use the utility function:
lower_corr = ca.corr_lower_triangle(data.drop(columns=["group"]), method="pearson")
lower_corr

# %% [markdown]
# Plot the lower triangle correlations as a histrogram to see the distribution of
# correlation values

# %%
ax = lower_corr.stack().plot.hist(
    bins=50,
    grid=False,
    figsize=(6, 4),
    title="Distribution of Pearson correlation values",
    xlabel="Pearson r",
    ylabel="Frequency",
    xlim=(-1.02, 1.02),
)

# %% [markdown]
# or to find the strongest correlations, which you might want to filter further for
# uninteresting correlation between redundant features.

# %%
lower_corr_stack = lower_corr.stack()
idx_largerst_corr = lower_corr_stack.abs().sort_values(ascending=False).head(20).index
lower_corr_stack.loc[idx_largerst_corr]

# %% [markdown]
# This function can be used to compute multiple correlation methods at once
# and compare them, here for the first four features.
#
# It only works on numeric values.

# %%
corr = list()
for method in ["pearson", "spearman", "kendall"]:
    _corr = (
        ca.corr_lower_triangle(data.iloc[:, :4], method=method, numeric_only=True)
        .stack()
        .rename(method)
    )
    corr.append(_corr)
corr = pd.concat(corr, axis=1).sort_values(by="pearson", ascending=True)
_ = corr.plot(
    style=".",
    ylim=(-1.05, 1.05),
    alpha=0.5,
    rot=45,
)

# %% [markdown]
# Filtering correlations based on p-values with multiple testing correction
# - the p-value depends on the number of samples
# - and the strenght of the correlation

# %%
res = ca.calculate_correlations(data.iloc[:, 0], data.iloc[:, 1], method="pearson")
print(res)

# %% [markdown]
# For the first four features, we would only keep one significant correlation
# after multiple testing correction with the Benjamini-Hochberg method.
#
# > Note that the p-value for the correlation between the same feature is set to one,
# > not zero as in scipy, to not consider these correlation in downstream analysis by
# > accident.

# %%
correlation = ca.run_correlation(
    data.iloc[:, :4], alpha=0.05, group="group", method="pearson", correction="fdr_bh"
)
correlation

# %% [markdown]
# The efficient correlation calculation can be used to compute the correlation
# matrix and p-value matrix for larger datasets.

# %%
corr, p = ca.run_efficient_correlation(data.iloc[:, :4], method="spearman")
pd.DataFrame(p)

# %% [markdown]
# you can verify the results against [`scipy.stats.spearmanr`](scipy.stats.spearmanr)

# %%
r, p = scipy.stats.spearmanr(data.iloc[:, :4])
pd.DataFrame(p)

# %% [markdown]
# same for pearson correlation

# %%
r, p = ca.run_efficient_correlation(data.iloc[:, :3], method="pearson")
pd.DataFrame(p)

# %%
r_20, p_20 = scipy.stats.pearsonr(data.iloc[:, 0], data.iloc[:, 2])
assert r[2, 0] - r_20 < 1e-8
assert p[2, 0] - p_20 < 1e-8
desc = (
    f"Correlation coefficient (r) in row three and the first column is {r_20:.3f} \n"
    f"and the corresponding p-value is {p_20:.3f}."
)
print(desc)

# %% [markdown]
# To calculate p-values for the correlation matrix, e.g. using pandas
# [`corr`](pandas.DataFrame.corr)
# you can use:

# %%
res = ca.calculate_pvalue_correlation(
    data.iloc[:, :3].corr(method="pearson").values, n_obs=data.shape[0]
)
pd.DataFrame(res)

# %% [markdown]
# Done.
