# %% [markdown]
# # Variation
#
# Using a metabolomics dataset ([MTBLS13411](../../../example_data/MTBLS13411.md)) as an
# example, we show how to use:
#
# - Coefficient of variation (CoV)
# - Histogram of the data
# - Summary statistics
#
# The focus of acore is not plotting, but getting the data in the right format for
# plotting. VueCore offers more plotting examples.

# %%
# %pip install acore

# %%
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import acore.exploratory_analysis as ea

# %% [markdown]
# Utility function for plotting


# %%
def make_plot(
    embeddings,
    x: str,
    y: str,
    annotation: Optional[dict[str, str]] = None,
    group: str = "group",
    **kwargs,
):
    """Utility function for static plot of dimensionality reductions."""
    fig, ax = plt.subplots()
    for i, (group, group_df) in enumerate(embeddings.groupby("group")):
        ax = group_df.rename(columns=map_names).plot.scatter(
            x=x,
            y=y,
            label=group,
            c=f"C{i}",
            ax=ax,
        )
    if annotation is not None:
        _ = ax.set(ylabel=annotation.y, xlabel=annotation.x)
    return fig, ax


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
# ## Coefficient of variation
# Using masspectrometry data, we can compute the coefficient of variation on the
# non-log transformed intensities. We do this for each group separately.
# First we undo the log transformation, which is something specific to this dataset.

# %%
data_exp = data.drop(columns=["group"]).apply(lambda x: np.exp2(x)).join(data["group"])
data_exp

# %%
res = ea.variation.get_coefficient_variation(data=data_exp, group="group")
res

# %%
res.describe()

# %%
map_names = {"x": "mean_log2", "y": "coef_of_var", "group": "group"}
fig, ax = make_plot(res, **map_names)


# %% [markdown]
# ## Histogram of the data
# We often want to plot the distribution of feature values.
# Sometimes you want to use custom bins, e.g. to align multiple histograms.
# Here we plot and compute a histogram frequencies with custom bins.

# %%
view = data.select_dtypes(include="number")
bins = np.arange(int(view.min(axis=None)), int(view.max(axis=None)) + 1, step=1)

for col in view.columns[:4]:
    ax = view[col].plot.hist(
        bins=bins,
        alpha=0.5,
        xlabel="Value",
        ylabel="Frequency",
    )
ax.title.set_text("Histogram with custom bins")
_ = ax.legend()

# %% [markdown]
# We could easily plot the histogram with custom bins, but if we want the count data
# of the histogram frequencies in a dataframe, `acore` has a utility function
# which can be used to compute the histogram frequencies for a feature.


# %%
hist_series = []
for col in view.columns[:4]:
    s = view[col]
    ret = ea.variation.get_histogram_series(s, bins)
    hist_series.append(ret.rename(col))

hist_df = pd.concat(hist_series, axis=1)
hist_df

# %% [markdown]
# ## Summary statistics
# You can use pandas to compute summary statistics for the data, specifying the
# percentiles you are interested in.

# %%
summary_stats = (
    data.select_dtypes(include="number").describe(percentiles=[0.5]).T
).rename(columns={"50%": "median"})
summary_stats

# %% [markdown]
# Done.
