# %% [markdown]
# # Enrichment analysis
#
# - we need some groups of genes to compute clusters
# - we need functional annotations, i.e. a category summarizing a set of genes.
# -
# You can start with watching Lars Juhl Jensen's brief introduction to enrichment analysis
# on [youtube](https://www.youtube.com/watch?v=2NC1QOXmc5o).
#
# Use example data for ovarian cancer
# ([PXD010372](https://github.com/Multiomics-Analytics-Group/acore/tree/main/example_data/PXD010372))


# %% tags=["hide-output"]
# %pip install acore


# %%
from pathlib import Path

import pandas as pd

import acore
import acore.differential_regulation
import acore.enrichment_analysis

# %% [markdown]
# Parameters of this notebook

# %% tags=["parameters"]
base_path: str = (
    "https://raw.githubusercontent.com/Multiomics-Analytics-Group/acore/refs/heads/main/"
    "example_data/PXD010372/processed"
)
omics: str = f"{base_path}/omics.csv"
meta_pgs: str = f"{base_path}/meta_pgs.csv"
meta: str = f"{base_path}/meta_patients.csv"
N_to_sample: int = 1_000

# %% [markdown]
# # Load processed data

# %%
df_omics = pd.read_csv(omics, index_col=0)
df_meta_pgs = pd.read_csv(meta_pgs, index_col=0)
df_meta = pd.read_csv(meta, index_col=0)
df_omics

# %%
df_omics.notna().sum().sort_values(ascending=True).plot()

# %% [markdown]
# Keep only features with a certain amount of non-NaN values and select 100 of these
# for illustration. Add the ones which were differently regulated in the ANOVA using all
# the protein groups.

# %%
idx_always_included = ["Q5HYN5", "P39059", "O43432", "O43175"]
df_omics[idx_always_included]

# %%
df_omics = (
    df_omics
    # .dropna(axis=1)
    .drop(idx_always_included, axis=1)
    .dropna(thresh=18, axis=1)
    .sample(
        N_to_sample - len(idx_always_included),
        axis=1,
        random_state=42,
    )
    .join(df_omics[idx_always_included])
)
df_omics


# %%
df_meta


# %% [markdown]
# ## Compute up and downregulated genes
# These will be used to find enrichments in the set of both up and downregulated genes.

# %%
group = "Status"
covariates = ["PlatinumValue"]
diff_reg = acore.differential_regulation.run_anova(
    df_omics.join(df_meta[[group]]),
    drop_cols=[],
    subject=None,
    group=group,
)
diff_reg.describe(exclude=["float"])

# %%
diff_reg.query("rejected == True")

# %% [markdown]
# ## Find functional annotations, here pathways
#

# %%
from acore.io.uniprot import (
    check_id_mapping_results_ready,
    get_id_mapping_results_link,
    get_id_mapping_results_search,
    submit_id_mapping,
)


def fetch_annotations(ids: pd.Index | list) -> pd.DataFrame:
    """Fetch annotations for UniProt IDs. Combines several calls to the API of UniProt's
    knowledgebase (KB).

    Parameters
    ----------
    ids : pd.Index | list
        Iterable of UniProt IDs. Fetches annotations as speecified by the specified fields.
    fields : str, optional
        Fields to fetch, by default "accession,go_p,go_c. See for availble fields:
        https://www.uniprot.org/help/return_fields

    Returns
    -------
    pd.DataFrame
        DataFrame with annotations of the UniProt IDs.
    """
    job_id = submit_id_mapping(from_db="UniProtKB_AC-ID", to_db="UniProtKB", ids=ids)

    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        # add fields to the link to get more information
        # From and Entry (accession) are the same for UniProt IDs.
        results = get_id_mapping_results_search(
            link + "?fields=accession,go_p,go_c,go_f&format=tsv"
        )
    header = results.pop(0).split("\t")
    results = [line.split("\t") for line in results]
    df = pd.DataFrame(results, columns=header)
    return df


fname_annotations = "downloaded/annotations.csv"
fname = Path(fname_annotations)
try:
    annotations = pd.read_csv(fname, index_col=0)
    print(f"Loaded annotations from {fname}")
except FileNotFoundError:
    print(f"Fetching annotations for {df_omics.columns.size} UniProt IDs.")
    annotations = fetch_annotations(df_omics.columns)
    annotations = (
        annotations.set_index("Entry")
        .rename_axis("identifier")
        .drop("From", axis=1)
        .rename_axis("source", axis=1)
        .stack()
        .to_frame("annotation")
        .replace("", pd.NA)
        .dropna()
        .sort_values(["source", "annotation"])
        .reset_index()
    )
    fname.parent.mkdir(exist_ok=True, parents=True)
    annotations.to_csv(fname, index=True)

annotations

# %% [markdown]
# ## Enrichment analysis
#

# %%
ret = acore.enrichment_analysis.run_regulation_enrichment(
    regulation_data=diff_reg,
    annotation=annotations,
    correction_alpha=0.01,
)
ret


# %%
