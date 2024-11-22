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
import pandas as pd

import acore
import acore.differential_regulation
import acore.enrichment_analysis

# %% [markdown]
# Parameters of this notebook

# %% tags=["parameters"]
omics: str = (
    "https://raw.githubusercontent.com/Multiomics-Analytics-Group/acore/refs/heads/main/"
    "example_data/PXD010372/processed/omics.csv"
)
meta: str = (
    "https://raw.githubusercontent.com/Multiomics-Analytics-Group/acore/refs/heads/main/"
    "example_data/PXD010372/processed/meta_patients.csv"
)

# %% [markdown]
# # Load processed data

# %%
df_omics = pd.read_csv(omics, index_col=0)
df_meta = pd.read_csv(meta, index_col=0)
df_omics

# %%
df_omics.notna().sum().sort_values(ascending=True).plot()

# %% [markdown]
# Keep only features with a certain amount of non-NaN values and select 100 of these
# for illustration

# %%
df_omics = df_omics.dropna(axis=1)
# .sample(
#     500,
#     axis=1,
#     random_state=42,
# )

# %%
df_meta

# %% [markdown]
# ## Compute up and downregulated genes
# These will be used to find enrichments in the set of both up and downregulated genes.
#
# - only CT45 shows to be regulated

# %%
group = "Status"
covariates = ["PlatinumValue"]
diff_reg = acore.differential_regulation.run_anova(
    df_omics.join(df_meta[[group]]),
    # df_omics.join(df_meta[[group] + covariates]),
    # covariates=covariates,
    drop_cols=[],
    subject=None,
    group=group,
)
diff_reg.describe(exclude=["float"])

# %%
diff_reg.describe(exclude=[float])

# %% [markdown]
# ## Find functional annotations, here pathways
#

# %%
from dsp_phospho.ptms.uniprot_id_mapping import (
    check_id_mapping_results_ready,
    get_id_mapping_results_link,
    get_id_mapping_results_search,
    submit_id_mapping,
)

job_id = submit_id_mapping(
    from_db="UniProtKB_AC-ID", to_db="UniProtKB", ids=df_omics.columns
)

if check_id_mapping_results_ready(job_id):
    link = get_id_mapping_results_link(job_id)
    # add fields to the link to get more information
    # From and Entry (accession) are the same.
    results = get_id_mapping_results_search(
        link + "?fields=accession,go_p,go_c,go_f&format=tsv"
    )
header = results.pop(0).split("\t")
results = [line.split("\t") for line in results]
results = pd.DataFrame(results, columns=header)
results

# %%
annotation = (
    results.set_index("Entry")
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
annotation

# %% [markdown]
# ## Enrichment analysis
#

# %%
ret = acore.enrichment_analysis.run_regulation_enrichment(
    regulation_data=diff_reg, annotation=annotation
)
ret

# %%
