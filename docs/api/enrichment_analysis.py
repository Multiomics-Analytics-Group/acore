# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: tags,-all
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.3
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Enrichment analysis
#
# This API example shows the main functions of the
# [`acore.enrichment_analysis`](acore.enrichment_analysis) module:
#
# - overrepresentation analysis (ORA) of up- and downregulated features with
#   [`run_up_down_regulation_enrichment`](acore.enrichment_analysis.run_up_down_regulation_enrichment)
# - single sample GSEA with
#   [`run_ssgsea`](acore.enrichment_analysis.run_ssgsea)
# - GSEA on a ranked list with [`run_gsea`](acore.enrichment_analysis.run_gsea)
#
# Inputs are
#
# 1. a differential regulation table, here computed with
#    [`acore.differential_regulation.run_anova`](acore.differential_regulation.run_anova)
# 2. functional annotations in long format with the columns `annotation` and
#    `identifier`, here GO terms fetched from UniProt with
#    [`acore.io.uniprot`](acore.io.uniprot)
#
# A complete analysis with interpretation of the results is shown in
# [the ovarian cancer analysis](../analyses/ovarian_cancer_enrichment.ipynb).

# %% tags=["hide-output"]
# %pip install acore

# %%
import pandas as pd

import acore.differential_regulation
import acore.enrichment_analysis

# %% tags=["parameters"]
base_path: str = (
    "https://raw.githubusercontent.com/Multiomics-Analytics-Group/acore/refs/heads/main/"
    "example_data/PXD010372/processed"
)
features_to_sample: int = 100

# %% [markdown]
# ## Inputs
# Load the ovarian cancer proteomics data (samples as rows, protein groups as columns)
# and select the same 100 protein groups for which UniProt annotations are cached in the
# repository.

# %% tags=["hide-input"]
df_omics = pd.read_csv(f"{base_path}/omics.csv", index_col=0)
df_meta = pd.read_csv(f"{base_path}/meta_patients.csv", index_col=0)
idx_always_included = ["Q5HYN5", "P39059", "O43432", "O43175"]
df_omics = (
    df_omics.drop(idx_always_included, axis=1)
    .dropna(thresh=18, axis=1)
    .sample(features_to_sample - len(idx_always_included), axis=1, random_state=42)
    .join(df_omics[idx_always_included])
)
df_omics.shape

# %% [markdown]
# Functional annotations (fetch them with
# `acore.io.uniprot.fetch_annotations` and `process_annotations`):

# %%
annotations = pd.read_csv(
    f"{base_path}/annotations_{features_to_sample}.csv", index_col=0
)
annotations

# %% [markdown]
# Differential regulation between deceased and alive patients:

# %%
group = "Status"
diff_reg = acore.differential_regulation.run_anova(
    df_omics.join(df_meta[[group]]),
    drop_cols=[],
    subject=None,
    group=group,
)
diff_reg["rejected"] = diff_reg["rejected"].astype(bool)
diff_reg.head()

# %% [markdown]
# ## Overrepresentation analysis (ORA)
# Up- and downregulated features (adjusted p-value and absolute log2 fold change above
# `lfc_cutoff`) are tested separately against all features using Fisher's exact test.
# The output follows the
# [`EnrichmentAnalysisSchema`](acore.types.enrichment_analysis.EnrichmentAnalysisSchema).

# %%
ora = acore.enrichment_analysis.run_up_down_regulation_enrichment(
    regulation_data=diff_reg,
    annotation=annotations,
    pval_col="padj",
    min_detected_in_set=1,
    lfc_cutoff=0.5,
)
ora

# %% [markdown]
# Use the permutation-based empirical FDR instead of Benjamini-Hochberg by setting
# `correction="efdr"`. See the
# [comparison with mulea](nb_ref_efdr_mulea) for a discussion of both procedures.

# %%
ora_efdr = acore.enrichment_analysis.run_up_down_regulation_enrichment(
    regulation_data=diff_reg,
    annotation=annotations,
    pval_col="padj",
    min_detected_in_set=1,
    lfc_cutoff=0.5,
    correction="efdr",
    permutations=1000,
    seed=42,
)
ora_efdr[["terms", "direction", "pvalue", "padj", "rejected"]].head()

# %% [markdown]
# ## Single sample GSEA (ssGSEA)
# Computes an enrichment score per sample and annotation from the data matrix.

# %%
ssgsea = acore.enrichment_analysis.run_ssgsea(
    data=df_omics,
    annotation=annotations,
    min_size=2,
)
ssgsea

# %% [markdown]
# ## Gene set enrichment analysis (GSEA)
# Ranks all features by `ranking_col` (here the log2 fold change) and tests whether the
# features of an annotation accumulate at the top or bottom of the ranking.

# %%
gsea = acore.enrichment_analysis.run_gsea(
    regulation_data=diff_reg,
    annotation=annotations,
    ranking_col="log2FC",
    min_size=1,
    permutation_num=100,
    seed=42,
)
gsea[["Term", "NES", "NOM p-val", "FDR q-val", "rejected"]]
