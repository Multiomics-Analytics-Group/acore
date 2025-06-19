# %% [markdown]
# # Example of running up/down regulation enrichment analysis on synthetic peptide data
# - load synthetic peptide data (can be manipulated to simulate different scenarios)
# - extend annotations to duplicate proteins with peptides
# %%
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import acore.enrichment_analysis as ea

folder = Path("data")
fname_reg_df = folder / "regulation_data.csv"
fname_annotations = folder / "annotations.csv"

# %% [markdown]
# The regulation data defines the total population of peptides:

# %%
reg_df = pd.read_csv(fname_reg_df)
reg_df

# %% [markdown]
# The pathways define the pathway specific set

# %%
annotations = pd.read_csv(fname_annotations)
annotations


# %% [markdown]
# Extend the data

# %%
annotations_extended = reg_df[["leading_protein", "identifier"]].join(
    annotations.set_index("identifier"), on="leading_protein"
)
annotations_extended

# %% [markdown]
# Run enrichment analysis

# %%
res = ea.run_up_down_regulation_enrichment(
    regulation_data=reg_df,
    annotation=annotations_extended,
    identifier="identifier",
    pval_col="padj",
    min_detected_in_set=1,  # ! default is 2, so more conservative
    lfc_cutoff=0.1,  # ! the default is 1
).reset_index(drop=True)
res

# %% [markdown]
# playing with the components of single test to understand
# fisher exact test and how it works:

# %%
from scipy.stats import fisher_exact

# %%
num_in_foreground = 4
num_in_background = 2
num_in_pathway = (
    num_in_background + num_in_foreground
)  # this is the number of peptides in the pathway
foreground_pop = 4  # constanst across pathways for up/down regulated)
background_pop = 19  # constanst across pathways for up/down regulated
print(f"{foreground_pop - num_in_foreground = }")
print(f"{background_pop - foreground_pop - num_in_background = }")
res = fisher_exact(
    [
        [num_in_foreground, foreground_pop - num_in_foreground],
        [num_in_background, background_pop - foreground_pop - num_in_background],
    ]
)
print(f"Fisher exact test p-value: {res.pvalue:.4f}")

# %% [markdown]
# We have a population of M peptides of which n are up-regulated (or down-regulated)
# defining the foreground.
# We find N peptides in a pathway of interest, of which x are in the foreground.
# - fixed for a certain foreground and background population (constant)
# - pathway size only inderectly taken into account through size of background population
# - foreground population and number of peptides in foreground are important

# %%
from scipy.stats import hypergeom

[M, n, N] = [background_pop, foreground_pop, num_in_background + num_in_foreground]
rv = hypergeom(M, n, N)
x = np.arange(0, n + 1)
pmf_dogs = rv.pmf(x)
print(f"pmf_dogs: {pmf_dogs}")
fig, ax = plt.subplots()
ax.plot(x, pmf_dogs, "bo")
ax.vlines(x, 0, pmf_dogs, lw=2)
ax.set_xlabel("# of peptides in foreground")
ax.set_ylabel("hypergeom PMF")
prb = hypergeom.cdf(num_in_foreground - 1, M, n, N)
print(
    f"Probability of finding {num_in_foreground} or more in foreground: {1 - prb:.4f}"
)
plt.show()

# %%
