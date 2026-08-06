# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.3
#   kernelspec:
#     display_name: .venv (3.11.0)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Permutation Tests
#
# In this notebook we will demonstrate how to use acore's permutation testing functions.
#
# ## Intro to permutation testing
#
# A hypothesis test, permutation tests compare an observed metric of your chosing (e.g.,
# t-statistic, mean difference) with metrics calculated when the dataset values are
# randomly shuffled permutations of the dataset.
#
# If we do 100 permutations of our data (although we should do a bunch more) and only 1 of
# those permutations falsely showed a larger effect size than the actual observed effect
# than it suggests there is a 1/100 chance (p value of 0.01) of the observed effect size
# having occurred by chance.
#
# ---
#
# ## The Dataset
#
# [Ju and colleagues (2018)](https://doi.org/10.1038/s41396-018-0277-8) collected
# microbiome samples from wastewater treatment plant (WWTP) influent
# ([MGYS00005056](https://www.ebi.ac.uk/metagenomics/studies/MGYS00005056/overview)) and
# effluent
# ([MGYS00005058](https://www.ebi.ac.uk/metagenomics/studies/MGYS00005058/overview)) from
# 12 WWTPs in Switzerland.
#
# In this demo we use a subset of the [GO
# term](https://geneontology.org/docs/go-annotations/) abundance data annotated by the
# [MGnify](https://www.ebi.ac.uk/metagenomics/) pipeline version 4.1. Specifically we look
# at [go term GO:0017001](https://www.ebi.ac.uk/QuickGO/term/GO:0017001) where it's
# expected that antibiotic catabolic processes to be higher in influent (INF) vs effluent
# (EFF) samples.
#
# The original datasets contained absolute abundance of selected GO terms for each sample,
# which we then transform to relative abundances and centred-log ratios. The [Data
# Preparation Details](#data-preparation-details) are in the dropdown below.
#
# The preprocessed and subsetted data for this example is provided in CSV,
# [`Ju2018_GO0017001_enf_inf_paired_demo_only.csv`](example_data/mgnify/Ju2018_GO0017001_enf_inf_paired_demo_only.csv).
# The data dictionary is below:
#
# | column            | description                                                                                                       | dtype |
# |-------------------|-------------------------------------------------------------------------------------------------------------------|-------|
# | sampling_location | The WWTP ID | str   |
# | eff_abundance     | The relative abundance of GO:0017001 for the EFF sample following preprocessing (i.e., CoDA and CLR)    | float |
# | inf_abundance     | The relative abundance of GO:0017001 for the INF sample following preprocessing (i.e., CoDA and CLR)    | float |
#
# ### Data Preparation Details
# <details>
# <summary style="color:green;"> More Info </summary>
#
# #### Downloading
# The analysed samples were downloaded via the [MGnify
# API](https://www.ebi.ac.uk/metagenomics/api/v2/). The influent (INF) and effluent
# (EFF) datasets have paired samples and we also needed to download the sample metadata
# (also available via MGnify API) to assign the correct pairing.
#
# #### Preprocessing of abundances
# - To account for technical variation due to sequencing technology limitations, we first
#   transform the abundance values so they are relative to the total reads for the sample
#   aka getting relative abundances.
# - The relative abundances are compositional data (CoDa) so we map them to unconstrained
#   vectors using centred log-ratio transformation
#   [`acore.transform.compositional.calc_clr`](acore.transform.compositional)
#   to not violate assumptions of any frequentist stats we do
#
# #### Preprocessing of the metadata
# - the sample metadata needed for this demo (sampling location) were available in their
#   "sample-desc"
# - the sample-desc for each sample in both INF and EFF were parsed and used for pairing
#   off
#
# #### Subsetted
# - We will only look at [go term GO:0017001](https://www.ebi.ac.uk/QuickGO/term/GO:0017001)
#
# #### Edit: A quick "fix" (not recommended)
# - Paired end sequences were treated as if separate observations/samples and annotated
#   using MGnify separately (e.g., analysed as if Read 1 is Sample 1, Read 2 is Sample 2
#   BUT this is incorrect because Read 1 and 2 are for the same sample and these reads
#   should have been combined prior to annotation)
# - Since this dataset is only for demonstrative purposes we crudely take a mean of the
#   already preprocessed read 1 and 2 abundances.
#
# <br>
# </details>
# <br>
#
# -----
#
# We will now proceed with reading in the prepared dataset as a `pandas.DataFrame`.

# %% [tags="hide-input"]
import pandas as pd

# reading in the data
df_data = pd.read_csv(
    # "https://raw.githubusercontent.com/Multiomics-Analytics-Group/acore/refs/heads/anglup-learning/"
    "example_data/mgnify/Ju2018_GO0017001_enf_inf_paired_demo_only.csv"
)
# sanity check
df_data

# %% [markdown]
# ---
#
# # Paired permutation test
#
# **Our null hypothesis:** The abundance of GO:0017001 is not different between INF and EFF samples.
#
# **Our alternative hypothesis:** The abundance of GO:0017001 is higher in INF than EFF samples.
#
#
# Since these are paired samples we will proceed with paired sample permutation test using
# [`acore.permutation_test.paired_permutation()`](acore.permutation_test.paired_permutation).
# Optional choice of random number generator for reproducibility.
#
# Here we will repeat the permutation test with 3 different metrics but this is for
# demonstrative purposes.

# %%
import numpy as np
from scipy.stats import ttest_rel

from acore.permutation_test import paired_permutation

# optional: use a random number generator for reproducibility
rng = np.random.default_rng(12345)

# for each metric, run the paired permutation test and print the results
for metric in [
    "t-statistic",
    ttest_rel,  # same as
    "mean",
    np.mean,
    "median",
    np.median,
]:

    # the paired permutation test
    result = paired_permutation(
        cond1=df_data["inf_abundance"],
        cond2=df_data["eff_abundance"],
        metric=metric,
        n_permutations=1000,
        rng=rng,
    )

    # verbosity
    print(result)

# %% [markdown]
# Based on the permutation tests by test statistic and mean difference, the probability of
# the observed metrics (t=4.770 and mean diff=0.535) occurring at random would be 0.001.
#
# ### The result
# Thus, we would reject the null hypothesis and accept the alternative.
#
# ---
#


# %% [markdown]
# ## The other tests
#


# %% [markdown]
#
# # Independent sample permutation test
#
# `paired_permutation` only makes sense when observations have a one-to-one
# correspondence, as our INF and EFF samples do. When two groups are unrelated
# (no pairing), we use
# [`acore.permutation_test.indep_permutation()`](acore.permutation_test.indep_permutation)
# instead. `indep_permutation` requires its inputs as numpy arrays.
#
# Now for demonstrative purposes we create observations by sampling from different (reject
# the null) and then the same distribution (fail to reject the null). We also use
# non-normal distributions to highlight the versality of permutation testing

# %% tags=["hide-input"]
# generating dummy continuous data
import matplotlib.pyplot as plt
import seaborn as sns

size = 100

first_dist = rng.uniform(low=1, high=10, size=size)
print(f"First distribution: mean = {np.mean(first_dist)}, std = {np.std(first_dist)}")

second_dist = rng.uniform(low=30, high=40, size=size)
print(
    f"Second distribution: mean = {np.mean(second_dist)}, std = {np.std(second_dist)}"
)

third_dist = rng.uniform(low=1, high=10, size=size)
print(f"Third distribution: mean = {np.mean(third_dist)}, std = {np.std(third_dist)}")

df_dist = pd.DataFrame([first_dist, second_dist, third_dist]).T
df_dist.columns = ["First", "Second", "Third"]

# plot the distributions
plt.figure(figsize=(10, 4))
sns.set_palette("colorblind")
sns.histplot(first_dist, label="first distribution", alpha=0.4, linestyle="solid")
sns.histplot(second_dist, label="second distribution", alpha=0.4, linestyle="dashed")
sns.histplot(third_dist, label="third distribution", alpha=0.4, linestyle="dotted")
plt.legend()
plt.show()

# %% ] [markdown]
# We first use the first and second distributions (different)

# %%
from scipy.stats import ks_2samp, ttest_ind

from acore.permutation_test import indep_permutation

metrics = ["t-statistic", ttest_ind, "mean", "median", ks_2samp]

for metric in metrics:
    result = indep_permutation(
        df_dist["First"],
        df_dist["Second"],
        metric=metric,
        n_permutations=1000,
        rng=rng,
    )
    # verbosity
    print(result)

# %% [markdown]
# and for example if testing 2 groups from similar distributions (first and third):

# %%
for metric in metrics:
    result = indep_permutation(
        df_dist["First"],
        df_dist["Third"],
        metric=metric,
        n_permutations=1000,
        rng=rng,
    )
    # verbosity
    print(result)

# %% [markdown]
# # Chi-squared permutation test
#
# [`acore.permutation_test.chi2_permutation()`](acore.permutation_test.chi2_permutation)
# tests whether the distribution of *categorical* observations differs between groups,
# using a permuted chi-squared statistic on a contingency table.
#
# Generate dummy categorical data:

# %% tags=["hide-input"]
first_group = rng.choice(["Happy", "Sad", "Neutral"], size=size, p=[0.5, 0.3, 0.2])
second_group = rng.choice(["Happy", "Sad", "Neutral"], size=size, p=[0.01, 0.19, 0.8])
third_group = rng.choice(["Happy", "Sad", "Neutral"], size=size, p=[0.5, 0.3, 0.2])

# display contingency table
df = pd.DataFrame(
    {
        "Group": (["First"] * size + ["Second"] * size + ["Third"] * size),
        "Mood": np.concatenate([first_group, second_group, third_group]),
    }
)
print(pd.crosstab(df["Group"], df["Mood"]))

# plot the distributions
plt.figure(figsize=(10, 4))
sns.set_palette("colorblind")
sns.histplot(first_group, label="first group", alpha=0.4, linestyle="solid")
sns.histplot(second_group, label="second group", alpha=0.4, linestyle="dashed")
sns.histplot(third_group, label="third group", alpha=0.4, linestyle="dotted")
plt.legend()
plt.show()

# %% [markdown]
# Run test

# %%
from acore.permutation_test import chi2_permutation

result = chi2_permutation(first_group, second_group, n_permutations=1000, rng=rng)
# verbosity
print("Results on dissimilar groups (1 and 2): ", result)

# now on similar groups (1 and 3)
result2 = chi2_permutation(first_group, third_group, n_permutations=1000, rng=rng)
print("Results on similar groups (1 and 3): ", result2)

# %% [markdown]
# Done.
