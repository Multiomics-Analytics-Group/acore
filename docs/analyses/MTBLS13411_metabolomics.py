# %% [markdown]
# # Phage-resistant _E. coli_: UMAP and differential analysis of metabolites
#
# An exploratory and differential analysis of the energy metabolome of enterotoxigenic
# _Escherichia coli_ (ETEC) from MetaboLights
# [MTBLS13411](https://www.ebi.ac.uk/metabolights/editor/MTBLS13411/files). Two
# phage-resistant mutants (`WRP` and `rapZE227Stop`) are compared to the wild type
# (`WT`), with six biological replicates each and three pooled quality control (`QC`)
# samples. The authors report that metabolites of the glucosamine-6-phosphate (GlcN6P)
# pathway, such as GlcNAc-6P, GlcNAc-1P and UDP-GlcNAc, accumulate in the resistant
# mutants. See the [dataset description](../example_data/MTBLS13411.md) for details.
#
# The analysis has two parts:
#
# 1. **UMAP**: do the samples group by strain, and do the QC samples sit in the middle?
# 2. **Differential analysis**: which metabolites differ
#    - between phage-resistant strains and the wild type (two groups), and
#    - between the three strains (ANOVA with pairwise posthoc tests).
#
# The single functions are documented in the API guide under
# [decomposition](../api/exploratory_analysis/decomposition.ipynb) and [differential
# regulation](../api/diff_reg/differential_regulation_two_groups.ipynb).

# %% tags=["hide-output"] %pip install acore vuecore

# %%
from pathlib import Path

import dsp_pandas
import pandas as pd
import vuecore.interactive
from vuecore.dimensionality_reduction import get_embedding_plot_interactive

import acore.differential_regulation as ad
import acore.exploratory_analysis as ea
import acore.normalization as an

dsp_pandas.format.set_pandas_options(max_colwidth=60)

# %% [markdown] Parameters of this notebook

# %% tags=["parameters"] folder in the acore repository, searched from the working
# directory upwards
local_dir: str = "example_data/MTBLS13411"
# fallback outside the repository (e.g. Google Colab): download from the main branch
github_url: str = (
    "https://raw.githubusercontent.com/Multiomics-Analytics-Group/acore/"
    "refs/heads/main/example_data/MTBLS13411"
)
strain_type: str = "Factor Value[Strain type]"  # phage-resistant or Wild Type
strain: str = "group"  # strain derived from the sample names: WT, WRP, rapZE227Stop
qc_group: str = "QC"  # pooled quality control samples
n_neighbors: int = 5  # small, as there are only 21 samples
alpha: float = 0.05  # FDR threshold

# %% [markdown]
# ## Load data
# The combined file holds the metabolite intensities together with the two metadata
# columns: the strain type and the strain (`group`). The intensities are already on a
# log2 scale and have no missing values.
#
# Inside a clone of the acore repository the file is read from the repository, which is
# faster and works offline. Otherwise, e.g. on Google Colab, it is downloaded from
# GitHub.

# %%
local_paths = [
    parent / local_dir
    for parent in (Path.cwd(), *Path.cwd().parents)
    if (parent / local_dir).is_dir()
]
base_path = local_paths[0] if local_paths else github_url
print(f"Loading data from: {base_path}")
fname = f"{base_path}/MTBLS13411_omics_and_meta.csv"
omics_and_meta = pd.read_csv(fname, index_col=0)
meta_cols = [strain_type, strain]
omics = omics_and_meta.drop(columns=meta_cols)
omics_and_meta[meta_cols].value_counts(dropna=False).to_frame()

# %% [markdown] Summary of the metabolite intensities per sample. The distributions are
# similar across samples, so no further normalization is applied here.

# %% tags=["hide-input"]
omics.T.describe().T.join(omics_and_meta[strain])

# %% [markdown]
# ## UMAP
# The metabolites cover different intensity ranges, so each metabolite is standardized
# (zero mean, unit variance) before computing the embedding. Otherwise the most variable
# metabolites would dominate the distances between samples. All samples, including the QC
# samples, are embedded.
#
# With only 21 samples, we use a small number of neighbors (`n_neighbors`). Check the
# [umap-learn documentation](https://umap-learn.readthedocs.io) for the parameters.

# %%
omics_scaled = an.normalize_data(omics, method="zscore", normalize="features")
umap_result, umap_annotation = ea.decomposition.run_umap(
    omics_scaled.join(omics_and_meta[strain]),
    drop_cols=[],
    group=strain,
    n_neighbors=n_neighbors,
    min_dist=0.3,
    metric="euclidean",
    dropna=True,
)
umap_coords = umap_result["umap"]
umap_coords.index = omics.index
umap_coords

# %%
figure = get_embedding_plot_interactive(
    umap_coords,
    annotation=umap_annotation,
    group=strain,
    title="UMAP of MTBLS13411 metabolite profiles",
    width=700,
    height=600,
)
figure

# %% [markdown] The wild type and the `rapZE227Stop` replicates each form a tight
# cluster, separated from each other. Five of the six `WRP` replicates form a third
# cluster, while `WRP-1` lies with the `rapZE227Stop` samples. This sample could be
# checked further before drawing conclusions. The three QC samples lie close together,
# which indicates a stable measurement over the run. UMAP preserves local neighborhoods
# only: distances between clusters are not directly interpretable, and the exact
# coordinates can change between runs. The differential analysis below quantifies which
# metabolites drive the separation.

# %% [markdown]
# ## Differential analysis: phage-resistant versus wild type
# The QC samples have no strain type and are removed. We compare both resistant mutants
# together against the wild type with a two-group ANOVA (equivalent to a t-test) on each
# metabolite. P-values are FDR corrected (Benjamini-Hochberg) and `rejected` marks
# metabolites with an adjusted p-value (`padj`) below the threshold. As the data is log2
# transformed, `log2FC` is the difference of the group means (`group1` versus `group2`).

# %%
two_groups = omics_and_meta.dropna(subset=[strain_type]).drop(columns=[strain])
anova_two_groups = (
    ad.run_anova(
        two_groups,
        subject=None,
        drop_cols=[],
        group=strain_type,
        alpha=alpha,
    )
    .set_index("identifier")
    .sort_values("padj")
)
anova_two_groups[["group1", "group2"]].drop_duplicates()

# %% [markdown] Number of differentially abundant metabolites:

# %%
anova_two_groups["rejected"].value_counts()

# %% [markdown] The top hits by adjusted p-value. The GlcN6P pathway metabolites
# N-acetylglucosamine 1-phosphate, UDP-N-acetylglucosamine and N-acetylglucosamine
# 6-phosphate are among the most significant, and are higher in the phage-resistant
# strains, matching the findings of the study.

# %%
anova_two_groups.filter(regex="mean|std|log2FC|pvalue|padj|rejected").head(20)

# %% [markdown] The volcano plot shows the log2 fold change against the unadjusted
# p-value. Red metabolites are significant after FDR correction.

# %%
vuecore.interactive.get_differential_regulation_plot(
    data=anova_two_groups,
)

# %% [markdown]
# ## Differential analysis: three strains
# The two resistant mutants carry different mutations, so we also test the three strains
# separately. The omnibus ANOVA asks whether any strain differs (`padj`, `rejected`).
# Pairwise posthoc t-tests then show which strains differ and in which direction
# (`posthoc pvalue adj`, `posthoc rejected`). Each metabolite has one row per pairwise
# comparison.

# %%
three_groups = omics_and_meta.query(f"{strain} != @qc_group").drop(
    columns=[strain_type]
)
anova_three_groups = ad.run_anova(
    three_groups,
    subject=None,
    drop_cols=[],
    group=strain,
    alpha=alpha,
).set_index(["identifier", "group1", "group2"])
anova_three_groups.filter(
    regex="^padj|^rejected|posthoc pvalue adj|posthoc rejected|log2FC"
).head(9)

# %% [markdown] Number of metabolites with a significant omnibus test, and number of
# significant metabolites per pairwise comparison:

# %%
omnibus = anova_three_groups.reset_index().drop_duplicates("identifier")
print(f"Significant in omnibus ANOVA: {omnibus['rejected'].sum()} of {len(omnibus)}")
anova_three_groups.groupby("posthoc comparison")["posthoc rejected"].sum().to_frame()

# %% [markdown]
# ### Posthoc volcano plots
# One volcano plot per pairwise comparison. Metabolites are colored by the posthoc FDR
# decision and the direction of the change. Hover over a point to see the metabolite.

# %%
pairwise_results = anova_three_groups.reset_index()
results_by_identifier = pairwise_results.set_index("identifier")
figures = {}
for group1, group2 in pairwise_results[["group1", "group2"]].drop_duplicates().values:
    comparison = f"{group1}~{group2}"
    figures[comparison] = vuecore.interactive.get_differential_regulation_plot(
        data=results_by_identifier,
        comparison=comparison,
        title="Posthoc volcano plot",
        width=700,
        height=500,
    )
list(figures)

# %%
figures["WRP~WT"]

# %%
figures["WT~rapZE227Stop"]

# %%
figures["WRP~rapZE227Stop"]

# %% [markdown]
# ### Metabolites changed in both mutants
# Metabolites that are significantly different from the wild type in both resistant
# mutants are candidates for a shared resistance mechanism. Note that the comparison with
# `rapZE227Stop` has the wild type as `group1`, so its fold change has the opposite sign.

# %%
vs_wt = (
    pairwise_results.query("group1 == 'WT' or group2 == 'WT'")
    .assign(
        mutant=lambda df: df["group1"].where(df["group2"] == "WT", df["group2"]),
        log2FC_vs_WT=lambda df: df["log2FC"].where(df["group2"] == "WT", -df["log2FC"]),
    )
    .pivot(
        index="identifier",
        columns="mutant",
        values=["log2FC_vs_WT", "posthoc pvalue adj", "posthoc rejected"],
    )
)
shared = vs_wt["posthoc rejected"].all(axis=1)
vs_wt.loc[shared].sort_values(("log2FC_vs_WT", "WRP"), ascending=False)

# %% [markdown] Done.
