# %% [metadata]
# # PXD010372
# %%
from pathlib import Path

import pandas as pd

import acore.io

# %%
fpath_meta = "meta.csv"  # in same folder as this script
furl_pg = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867418311668-mmc2.xlsx"

# %% [markdown]
# Where to save the dat

# %%
dir_proc = Path("processed")
dir_proc.mkdir(exist_ok=True)

# %%
local_filename = Path(furl_pg).name
acore.io.download_file(furl_pg, local_filename)

# %%
meta_patients = pd.read_csv(fpath_meta, index_col=0)
meta_patients

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

# %%
patient_int_map = (
    pgs.columns.to_frame()
    .squeeze()
    .str.extract(r"(\d+)")
    .squeeze()
    .astype(int)
    .to_frame("Patient")
    .reset_index()
    .set_index("Patient")
    .rename({"index": "Patient_id"}, axis=1)
)
patient_int_map

# %%
meta_patients = meta_patients.join(patient_int_map).set_index(
    "Patient_id",
)
fname = dir_proc / "meta_patients.csv"
meta_patients.to_csv(fname)
print("Saved to ", fname)
meta_patients

# %%
meta_pgs = data.drop(pgs.columns, axis=1).drop("Patient11B", axis=1)
fname = dir_proc / "meta_pgs.csv"
meta_pgs.to_csv(fname)
print("Saved to ", fname)
meta_pgs

# %%
meta_pgs.describe(exclude="number")

# %%
mask = meta_pgs["Potential contaminant"].isna()
mask.value_counts()

# %%
pgs = pgs.loc[mask].T
fname = dir_proc / "omics.csv"
pgs.to_csv(fname)
print("Saved to ", fname)
pgs

# %% [markdown]
# Need to decide how to handle the duplicated patient measurement for patient 11.
