# %%
from pathlib import Path

import acore.io as aio

id = "PXD008541"
fname = "SearchEngineResults_secretome.zip.rar"
folder_downloads = Path("downloaded")
folder_unzipped = Path("unzipped")

# %%
ret = aio.download_PRIDE_data(id, fname, to="downloaded")
ret

# %% [markdown]
# see ftp link for PRIDE project

# %%
ret["_links"]["datasetFtpUrl"]["href"]


# %%
# ! you need a system installation of a rar archive tool
aio.unrar(filepath=folder_downloads / fname, to=folder_unzipped)

# %%
list(folder_unzipped.iterdir())

# %%
