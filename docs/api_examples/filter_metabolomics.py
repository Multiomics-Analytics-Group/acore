# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
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


import matplotlib.pyplot as plt

# %%
import pandas as pd

from acore import filter_metabolomics as fm


# %% tags=["hide-input"]
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

# %% [markdown]
# Load in your data. The example data set can be found in example_data/DidacMauricio_hilic.

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
# In order to run our further analysis, including the filtering functions, we have to transform the data and remove metadata such as mass and retention time.

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
# The rows are our samples, with the row index being the sample names. The columns are our individual features. All metadata has been removed.

# %% [markdown]
# Let's also define a variable that contains all the names of samples that belong to each group. This will make it easier later-on to perform filtering based on different variables, controls or conditions.

# %%
blanks = ["Bf_1", "Bf_2", "Bi_1", "Bi_2"]
qcs = [
    "QC_00",
    "QC_01",
    "QC_02",
    "QC_03",
    "QC_04",
    "QC_05",
    "QC_06",
    "QC_07",
    "QC_08",
    "QC_09",
    "QC_10",
    "QC_11",
    "QC_12",
    "QC_13",
    "QC_14",
    "QC_15",
    "QC_16",
    "QC_17",
    "QC_18",
    "QC_19",
    "QC_20",
    "QC_21",
    "QC_22",
    "QC_23",
    "QC_24",
    "QC_25",
    "QC_26",
    "QC_27",
    "QC_28",
    "QC_29",
    "QC_30",
    "QC_31",
    "QC_32",
    "QC_33",
    "QC_34",
    "QC_35",
    "QC_36",
    "QC_37",
    "QC_38",
    "QC_39",
    "QC_40",
    "QC_41",
    "QC_42",
    "QC_43",
    "QC_44",
]
samples = [label for label in data.index if label not in blanks and label not in qcs]

print(
    f"Checking... all group labels accounted for in our lists: {(len(blanks)+len(qcs)+len(samples)) == len(data.index)}"
)

samples_a = [idx for idx in data.index if idx.startswith("AAA")]
samples_p = [idx for idx in data.index if idx.startswith("P")]

print(
    f"Checking... all sample names accounted for in our lists: {(len(samples_a)+len(samples_p)) == len(samples)}"
)


# %% [markdown]
# Now we are ready to filter our data.

# %% [markdown]
# ### Filtering by missingness: 80%-rule
#
# The 80%-rule filters out features with too much missingness from our data. More specifically, if for a feature, more than 20% of the data is missing across all sample columns, it will be removed, so features must have at least 80% of data present in order to be retained.
#
# Although it is called the 80%-rule, other thresholds can be used to make the filtering more lenient or more stringent.
#
# In acore, this method is implemented in the function filter_by_missingness. Let's first have a look at our function.
#
#

# %% tags=["hide-input"]
help(fm.filter_by_missingness)

# %% [markdown]
# Now that we know how to use the function, we can run it, using the default 80%.

# %%
# 80% rule classic with percent=80
data_missingness_80_classic = fm.filter_by_missingness(
    data, method="classic", samples=samples
)

# %%
print(
    f"Num. of features before filtering: {data.shape[1]}\nNum. of features after filtering: {data_missingness_80_classic.shape[1]}"
)
print(f"Difference: {data.shape[1]-data_missingness_80_classic.shape[1]} rows removed.")

# %% [markdown]
# Let's see how our filtering changes when we apply a different threshold.

# %%
# 80% rule classic with percent=60
data_missingness_60_classic = fm.filter_by_missingness(
    data, percent=60, method="classic", samples=samples
)

# %%
print(
    f"Num. of features before filtering: {data.shape[1]}\nNum. of features after filtering: {data_missingness_60_classic.shape[1]}"
)
print(f"Difference: {data.shape[1]-data_missingness_60_classic.shape[1]} rows removed.")

# %% [markdown]
# Now we can also use the modified 80%-rule. This method divides into sample groups and computes the missingness per group. A feature survives filtering if it meets the missingness requirements in at least one group.
# The idea behind this is making sure that if there is "the perfect biomarker" in our data, meaning that there is a feature which shows up very strongly in one experimental condition and not at all in another condition, it is not filtered out.

# %%
# 80% rule modified with percent=80
data_missingness_80_modified = fm.filter_by_missingness(
    data,
    percent=80,
    method="modified",
    groups={"samples_a": samples_a, "samples_p": samples_p},
)

# %%
print(
    f"Num. of features after filtering with classic method: {data_missingness_80_classic.shape[1]}\nNum. of features after filtering with modified method: {data_missingness_80_modified.shape[1]}"
)
print(
    f"Difference: {data_missingness_80_modified.shape[1]-data_missingness_80_classic.shape[1]} more rows retained than when using classic method."
)

# %% [markdown]
# ### Filtering by Coefficient of Variation (CV)

# %% [markdown]
# In this method, we are taking into account the quality control (QC) samples.
#
# The CV of the biological samples and the CV of the QC samples are calculated per feature, and if for a given feature the CV of the QC samples is larger than that of the biological samples, it is removed.
#
# In acore, this method is implemented in the function filter_cv.

# %% tags=["hide-input"]
help(fm.filter_cv)

# %%
data_cv = fm.filter_cv(data=data, samples=samples, qcs=qcs)

# %%
print(
    f"Num. of features before filtering: {data.shape[1]}\nNum. of features after filtering: {data_cv.shape[1]}"
)
print(f"Difference: {data.shape[1]-data_cv.shape[1]} rows removed.")

# %% [markdown]
# ### Filtering with the Blanks control: Removing background noise and carryover
#
# This method removes features that have too high intensities in the Blanks control, measured by the ratio of the mean intensity in Blanks to the mean intensity in biological samples. The default threshold is 0.5, meaning that a featuer gets removed if its mean intensity in the Blanks is half as large as its mean intensity in samples.
#

# %% tags=["hide-input"]
help(fm.filter_blanks)

# %% [markdown]
# First, we can check whether there is signal in the blanks samples. For that, we can plot their total ion chromatograms (TICs).

# %%
plot_tic(data, blanks, ylim=150000000)

# %% [markdown]
# There is some signal. So we can run the blanks filtering, filtering out features with the default threshold.

# %%
data_blanks_05 = fm.filter_blanks(data, blanks=blanks, samples=samples)

# %%
print(
    f"Num. of features before filtering: {data.shape[1]}\nNum. of features after filtering: {data_blanks_05.shape[1]}"
)
print(f"Difference: {data.shape[1]-data_blanks_05.shape[1]} rows removed.")

# %%
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

# %%
print(
    f"Num. of features before filtering: {data.shape[1]}\nNum. of features after filtering: {data_blanks_01.shape[1]}"
)
print(f"Difference: {data.shape[1]-data_blanks_01.shape[1]} rows removed.")

# %%
plot_tic(
    data_blanks_01,
    blanks,
    ylim=150000000,
    title="TIC: Blanks intensities after filtering with threshold=0.1",
)

# %% [markdown]
# This threshold is more stringent so less features are retained which also means that there is less total intensity in the blanks samples.
