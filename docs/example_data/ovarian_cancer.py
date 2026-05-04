# %% [markdown]
# (nb_ref_ovarian_data)=
# # Download from journal (ovarian cancer proteome)
# Download the ovarian cancer proteome data from the journal's website. It was
# provided as supplementary data. See the article here:
#
# > Fabian Coscia, Ernst Lengyel, Jaikumar Duraiswamy, Bradley Ashcroft, Michal Bassani-Sternberg, Michael Wierer, Alyssa Johnson, Kristen Wroblewski, Anthony Montag, S. Diane Yamada, Blanca López-Méndez, Jakob Nilsson, Andreas Mund, Matthias Mann, Marion Curtis,
# > Multi-level Proteomics Identifies CT45 as a Chemosensitivity Mediator and Immunotherapy Target in Ovarian Cancer,
# > Cell,
# > Volume 175, Issue 1,
# > 2018,
# >
# > https://doi.org/10.1016/j.cell.2018.08.065.

# %%
# %pip install acore openpyxl

# %%
from pathlib import Path

import pandas as pd

import acore.io

# %% [markdown]
# Specify the proteome file's url

# %%
furl_pg = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867418311668-mmc2.xlsx"

# %% [markdown]
# Load it using an acore function

# %%
local_filename = Path(furl_pg).name
acore.io.download_file(furl_pg, local_filename)

# %% [markdown]
# Open the excel file from the supplementary data given in the article.

# %%
data = pd.read_excel(
    local_filename,
    sheet_name="SupplementaryTable2_PatientProt",
)
data

# %% [markdown]
# We will use the first protein in the a protein group as identifier,
# which we verify to be unique.

# %%
data["first_prot"] = data["Majority protein Ids"].str.split(";").str[0]
data["first_prot"].nunique() == data["Majority protein Ids"].nunique()
data = data.set_index("first_prot")
assert data.index.is_unique
data

# %% [markdown]
# Filter intensity values for patients

# %%
pgs = data.filter(like="Patient")
pgs

# %% [markdown]
# There are two measurements for patient 11 in the data: 11 and 11B. In the methods of
# the paper it is stated:
#
# "We required a minimum peptide ratio count of 1 to report a quantitative readout
# and averaged the results from duplicate measurements of the same sample."
#
# We will do this manually for patient 11 measurements.

# %%
pgs.filter(like="Patient11").describe()

# %%
pgs = pgs.assign(Patient11=lambda df: df.filter(like="Patient11").mean(axis=1)).drop(
    ["Patient11B"], axis=1
)

# %% [markdown]
# Keep the other information of protein groups as additional annotations on protein groups.

# %%
meta_pgs = data.drop(pgs.columns, axis=1).drop("Patient11B", axis=1)
meta_pgs

# %% [markdown]
# View non-numeric columns of protein group metadata.

# %%
meta_pgs.describe(exclude="number")

# %% [markdown]
# Get rid of potential contaminants (marked with a +, so non missing have `NAN`).

# %%
mask = meta_pgs["Potential contaminant"].isna()
pgs = pgs.loc[mask].T
pgs

# %% [markdown]
# ## Patient metadata
# The patient metadata was only provided as a pdf file. We parsed it and
# saved it as a csv file. You can load it for our GitHub repository:

# %%
# ! ToDo: link file after it is available on the main branch

# %% [markdown]
# Done.
