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
# %pip install acore vuecore


# %%
from pathlib import Path

import dsp_pandas
import pandas as pd

import acore
import acore.differential_regulation
import acore.enrichment_analysis
from acore.io.uniprot import fetch_annotations

dsp_pandas.format.set_pandas_options(max_colwidth=15)

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
features_to_sample: int = 100

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

# %% tags=["hide-input"]
df_omics = (
    df_omics
    # .dropna(axis=1)
    .drop(idx_always_included, axis=1)
    .dropna(thresh=18, axis=1)
    .sample(
        features_to_sample - len(idx_always_included),
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
diff_reg["rejected"] = diff_reg["rejected"].astype(bool)  # ! needs to be fixed in anova
diff_reg.query("rejected")

# %% [markdown]
# ## Find functional annotations, here pathways
#

# %%
fname_annotations = f"downloaded/annotations_{features_to_sample}.csv"
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
# See how many protein groups are associated with each annotation.

# %% tags=["hide-input"]
_ = (
    annotations.groupby("annotation")
    .size()
    .value_counts()
    .sort_index()
    .plot(kind="bar")
)

# %% [markdown]
# ## Enrichment analysis
#

# %%
ret = acore.enrichment_analysis.run_regulation_enrichment(
    regulation_data=diff_reg,
    annotation=annotations,
    min_detected_in_set=1,  # ! default is 2, so more conservative
    correction_alpha=0.01,
)
ret

# %% [markdown]
# ### For up- and downregulated genes separately

# %% tags=["hide-input"]
diff_reg.query("rejected")[
    [
        "identifier",
        "group1",
        "group2",
        "pvalue",
        "padj",
        "rejected",
        "log2FC",
        "FC",
    ]
].sort_values("log2FC")

# %% [markdown]
# - this additionally sets a fold change cutoff
# - and the fore and backgroud populations are changed due to the separation

# %%
ret = acore.enrichment_analysis.run_up_down_regulation_enrichment(
    regulation_data=diff_reg,
    annotation=annotations,
    min_detected_in_set=1,  # ! default is 2, so more conservative
    lfc_cutoff=1,
)
ret

# %% [markdown]
# here we see differences for the same set of differently regulated protein groups,
# which can be reset using lfc_cutoff=0.

# %%
ret = acore.enrichment_analysis.run_up_down_regulation_enrichment(
    regulation_data=diff_reg,
    annotation=annotations,
    min_detected_in_set=1,  # ! default is 2, so more conservative
    lfc_cutoff=0,
)
ret

# %% [markdown]
# ### Site specific enrichment analysis

# %% [markdown]
# The basic example uses a modified peptide sequence to
# demonstrate the enrichment analysis.
# - compare groups per amino acid modified (kinases targeting certain motifs?)

# %%
import re

regex = "(\\w+~.+)_\\w\\d+\\-\\w+"
identifier_ckg = "gnd~P00350_T10-WW"
match = re.search(regex, identifier_ckg)
match.group(1)

# %%
seq_mod = "_AAADQET(ph)DTDPEPQPVVGPDAADHRPTVM(ox)LLGGGALSR"
regex = "\(\w\w\)"
matches = re.findall(regex, seq_mod)
matches

# %%
# acore.enrichment_analysis.run_up_down_regulation_enrichment(

# %% [markdown]
# ## Single sample GSEA (ssGSEA)
# Run a gene set enrichment analysis (GSEA) for each sample,
# see [article](https://www.nature.com/articles/nature08460#Sec3) and
# the package [`gseapy`](https://gseapy.readthedocs.io/en/latest/run.html#gseapy.ssgsea)
# for more details.

# %%
enrichtments = acore.enrichment_analysis.run_ssgsea(
    data=df_omics,
    annotation=annotations,
    min_size=1,
)
enrichtments

# %%
enrichtments.iloc[0].to_dict()

# %%
enrichtments["NES"].plot.hist()

# %% [markdown]
# The normalised enrichment score (NES) can be used in a PCA plot to see if the samples
# cluster according to the enrichment of the gene sets.

# %%
nes = enrichtments.set_index("Term", append=True).unstack()["NES"].convert_dtypes()
nes

# %%
import acore.exploratory_analysis as ea

pca_result, pca_annotation = ea.run_pca(
    data=nes.join(df_meta[[group]]),
    drop_cols=[],
    annotation_cols=[],
    group=group,
    components=2,
    dropna=False,
)
resultDf, loadings, var_exp = pca_result
resultDf

# %% [markdown]
# The loadings show how the variables are correlated with the principal components.

# %%
loadings

# %% [markdown]
# We will plot both on the sample plot (samples on the first two principal components and
# loadings of variables). We use the
# [`vuecore` package](https://github.com/Multiomics-Analytics-Group/vuecore)
# for this, which is also developed by the Multiomics Analytics Group.

# %%
from plotly.offline import iplot
from vuecore import viz

args = {"factor": 1, "loadings": 10}
# #! pca_results has three items, but docstring requests only two -> double check
figure = viz.get_pca_plot(data=pca_result, identifier="PCA enrichment", args=args)
iplot(figure)

# %% [markdown]
# ## Compare two distributions - KS test
# The Kolmogorov-Smirnov test is a non-parametric test that compares two distributions.
# - we compare the distributions of the two differently upregulated protein groups
# This is not the best example for comparing distributions, but it shows how to use the
# KS test.

# %%
# plot two histograms of intensity values here
sel_pgs = ["O43175", "P39059"]
view = df_omics[sel_pgs].sub(df_omics[sel_pgs].mean())
view.plot.hist(bins=20, alpha=0.5)

# %% [markdown]
# Let us compare the two centered distributions using the KS test.

# %%
acore.enrichment_analysis.run_kolmogorov_smirnov(view[sel_pgs[0]], view[sel_pgs[1]])

# %% [markdown]
# The result suggests that the two distributions are from the same distribution.
