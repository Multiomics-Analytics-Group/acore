# %% [metadata]
# # PXD010372
# %%
from pathlib import Path

import pandas as pd
import requests


def download_file(url, local_filename):
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        with open(local_filename, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)


# %%
fpath_meta = "meta.csv"  # in same folder as this script
furl_pg = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867418311668-mmc2.xlsx"

# %%
local_filename = Path(furl_pg).name
download_file(furl_pg, local_filename)

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
# 11 and 11B ?

# %%
pgs.filter(like="Patient11").describe()

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
meta_patients.join(patient_int_map)

# %%
meta_pgs = data.drop(pgs.columns, axis=1)
meta_pgs

# %%
meta_pgs.describe(exclude="number")

# %%
pgs.T

# %% [markdown]
# Need to decide how to handle the duplicated patient measurement for patient 11.
