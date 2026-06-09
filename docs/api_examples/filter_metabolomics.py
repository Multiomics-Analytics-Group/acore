# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.3
#   kernelspec:
#     display_name: acore-dev
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Metabolomics data filtering
#
# This example notebook shows how to use the acore filtering functions specific for metabolomics data.
#
# Scroll to read about the following methods:
# - 80%-rule
# - modified 80%-rule
# - CV-based filtering
# - Blanks-based filtering
#
# Apply whichever one fits your data best, or a combination of multiple. Adjust the thresholds to match your data and desired stringency.

# %% tags=["hide-output"]
# %pip install acore

# %% tags=["hide-input"]
import matplotlib.pyplot as plt
import pandas as pd

from acore import filter_metabolomics as fm


def plot_tic(
    data,
    blanks,
    threshold=0.5,
    figsize=(8, 4),
    color="#3266ad",
    ylim=None,
    title="TIC: Blanks intensities before filtering",
):
    data.loc[blanks].sum(axis=1).plot(kind="bar", figsize=figsize, color=color)

    plt.title(title)
    plt.xlabel("Sample")
    plt.ylabel("Total intensity")
    plt.xticks(rotation=0)
    if ylim is not None:
        plt.ylim(0, ylim)
    plt.tight_layout()
    plt.show()


# %% [markdown]
# #### Data preparation
# Load in your data. The example data set can be found in
# `example_data/DidacMauricio_hilic`.

# %%
data_path = (
    "https://raw.githubusercontent.com/Multiomics-Analytics-Group/acore/"
    "refs/heads/main/"
)
data_path = (
    "../../example_data/DidacMauricio_hilic/DM_FIS2018_Hilic_pos_results2023_filled.csv"
)
data_original = pd.read_csv(data_path)

# %% [markdown]
# We can now inspect our data:

# %% tags=["hide-input"]
data_original

# %% [markdown]
# In order to run our further analysis, including the filtering functions,
# we have to transform the data and remove metadata such as mass and retention time.

# %%
data = data_original.T
data = data.drop(
    ["Qidx", "SOIidx", "rtmed", "start", "end", "mass", "MaxInt", "formula", "anot"]
)

# %% [markdown]
# Let's see what our data looks like now:

# %% tags=["hide-input"]
data

# %% [markdown]
# The rows are our samples, with the row index being the sample names. The columns are our
# individual features. All metadata has been removed.

# %% [markdown]
# Let's also define a variable that contains all the names of samples that belong to each
# group. This will make it easier later-on to perform filtering based on different
# variables, controls or conditions.

# %% tags=["hide-input"]
blanks = ["Bf_1", "Bf_2", "Bi_1", "Bi_2"]
qcs = [idx for idx in data.index if idx.startswith("QC_")]
samples = [label for label in data.index if label not in blanks and label not in qcs]

print(
    "Checking... all group labels accounted for in our lists:"
    f" {(len(blanks)+len(qcs)+len(samples)) == len(data.index)}"
)

samples_a = [idx for idx in data.index if idx.startswith("AAA")]
samples_p = [idx for idx in data.index if idx.startswith("P")]

print(
    "Checking... all sample names accounted for in our lists: "
    f"{(len(samples_a)+len(samples_p)) == len(samples)}"
)

s_groups = pd.Series("other", index=data.index)
s_groups.loc[samples_a] = "samples_a"
s_groups.loc[samples_p] = "samples_p"
s_groups.loc[qcs] = "qcs"
s_groups.loc[blanks] = "blanks"

# %% [markdown]
# Now we are ready to filter our data.

# %% [markdown]
# ### Filtering by missingness: 80%-rule
#
# The 80%-rule filters out features with too much missingness from our data. More
# specifically, if for a feature, more than 20% of the data is missing across all sample
# columns, it will be removed, so features must have at least 80% of data present in order
# to be retained.
#
# Although it is called the 80%-rule, other thresholds can be used to make the filtering
# more lenient or more stringent.
#
# In acore, this method is implemented in the function filter_by_missingness. Let's first
# have a look at our function.


# %% tags=["hide-input"]
help(fm.filter_by_missingness)

# %% [markdown]
# Now that we know how to use the function, we can run it, using the default 80%.

# %%
# 80% rule classic with percent=80
data_missingness_80_classic = fm.filter_by_missingness(
    data, method="classic", samples=samples
)

# %% tags=["hide-input"]
print(
    f"Num. of features before filtering: {data.shape[1]}\n"
    f"Num. of features after filtering: {data_missingness_80_classic.shape[1]}"
)
print(f"Difference: {data.shape[1]-data_missingness_80_classic.shape[1]} rows removed.")

# %% [markdown]
# Let's see how our filtering changes when we apply a different threshold.

# %%
# 80% rule classic with percent=60
data_missingness_60_classic = fm.filter_by_missingness(
    data, percent=60, method="classic", samples=samples
)

# %% tags=["hide-input"]
print(
    f"Num. of features before filtering: {data.shape[1]}\n"
    f"Num. of features after filtering: {data_missingness_60_classic.shape[1]}"
)
print(f"Difference: {data.shape[1]-data_missingness_60_classic.shape[1]} rows removed.")

# %% [markdown]
# Now we can also use the modified 80%-rule. This method divides into sample groups and
# computes the missingness per group. A feature survives filtering if it meets the
# missingness requirements in at least one group. The idea behind this is making sure that
# if there is "the perfect biomarker" in our data, meaning that there is a feature which
# shows up very strongly in one experimental condition and not at all in another
# condition, it is not filtered out.

# %%
# 80% rule modified with percent=80
data_missingness_80_modified = fm.filter_by_missingness(
    data,
    percent=80,
    method="modified",
    groups={"samples_a": samples_a, "samples_p": samples_p},
)

# %% tags=["hide-input"]
print(
    "Num. of features after filtering with classic method: "
    f"{data_missingness_80_classic.shape[1]}\nNum. of features after filtering "
    f"with modified method: {data_missingness_80_modified.shape[1]}"
)
features_only_in_modified = data_missingness_80_modified.columns.difference(
    data_missingness_80_classic.columns
)
print(
    f"Difference: {len(features_only_in_modified)} more features retained"
    " by modified method than classic method."
)

# %% [markdown]
# Features retained by the modified 80%-rule but removed by the classic rule —
# shown in the original data (samples only, no metadata). It shows that the featues
# retained additionally have in one group missingness above 20%.
#
# > note that the group sizes are unqual in this example.

# %% tags=["hide-input"]
(
    data[features_only_in_modified]
    .notna()
    .groupby(s_groups)
    .mean()
    .T[["samples_a", "samples_p"]]
    .sort_values(by=["samples_a", "samples_p"], ascending=[False, True])
    .assign(
        global_average=lambda df: df.multiply(
            s_groups.value_counts().loc[["samples_a", "samples_p"]], axis="columns"
        )
        .sum(axis=1)
        .div(s_groups.value_counts().loc[["samples_a", "samples_p"]].sum())
    )
)

# %%
counts = s_groups.value_counts()
(
    data[features_only_in_modified]
    .notna()
    .groupby(s_groups)
    .mean()
    .multiply(counts, axis="index")  # weight each row by its group size
    .sum()
    .div(counts.sum())  # normalize to get overall weighted mean
)


# %% [markdown]
# ### Filtering by Coefficient of Variation (CV)

# %% [markdown]
# In this method, we are taking into account the quality control (QC) samples.
#
# The CV of the biological samples and the CV of the QC samples are calculated per
# feature, and if for a given feature the CV of the QC samples is larger than that of the
# biological samples, it is removed.
#
# In acore, this method is implemented in the function filter_cv.

# %% tags=["hide-input"]
help(fm.filter_cv)

# %%
data_cv = fm.filter_cv(data=data, samples=samples, qcs=qcs)

# %% tags=["hide-input"]
print(
    f"Num. of features before filtering: {data.shape[1]}\n"
    f"Num. of features after filtering: {data_cv.shape[1]}"
)
print(f"Difference: {data.shape[1]-data_cv.shape[1]} rows removed.")

# %% [markdown]
# ### Filtering with the blanks control: Removing background noise and carryover
#
# This method removes features that have too high intensities in the blanks control,
# measured by the ratio of the mean intensity in blanks to the mean intensity in
# biological samples. The default threshold is 0.5, meaning that a featuer gets removed if
# its mean intensity in the blanks is half as large as its mean intensity in samples.
#

# %% tags=["hide-input"]
help(fm.filter_blanks)

# %% [markdown]
# First, we can check whether there is signal in the blanks samples. For
# that, we can plot their total ion chromatograms (TICs).

# %%
plot_tic(data, blanks, ylim=150000000)

# %% [markdown]
# There is some signal. So we can run the blanks filtering, filtering out
# features with the default threshold.

# %%
data_blanks_05 = fm.filter_blanks(data, blanks=blanks, samples=samples)

# %% tags=["hide-input"]
print(
    f"Num. of features before filtering: {data.shape[1]}\n"
    f"Num. of features after filtering: {data_blanks_05.shape[1]}"
)
print(f"Difference: {data.shape[1]-data_blanks_05.shape[1]} rows removed.")

plot_tic(
    data_blanks_05,
    blanks,
    ylim=150000000,
    title="TIC: Blanks intensities after filtering with threshold=0.5",
)

# %% [markdown]
# Keeping the y axis at the same scale, we can see that the total intensity has changed.
# Now we can check how it looks if we use a different parameter than the default.

# %%
data_blanks_01 = fm.filter_blanks(data, blanks=blanks, samples=samples, threshold=0.1)

# %% tags=["hide-input"]
print(
    f"Num. of features before filtering: {data.shape[1]}\n"
    f"Num. of features after filtering: {data_blanks_01.shape[1]}"
)
print(f"Difference: {data.shape[1]-data_blanks_01.shape[1]} rows removed.")

plot_tic(
    data_blanks_01,
    blanks,
    ylim=150000000,
    title="TIC: Blanks intensities after filtering with threshold=0.1",
)

# %% [markdown]
# This threshold is more stringent so less features are retained which also
# means that there is less total intensity in the blanks samples.
