# %% [markdown]
# # Enrichment analysis
# requires
# - some cluster of proteins/genes (e.g. up- and downregulated proteins/genes)
# - functional annotations, i.e. a category summarizing a set of proteins/genes.
#
# You can start with watching Lars Juhl Jensen's brief introduction to enrichment analysis
# on [youtube](https://www.youtube.com/watch?v=2NC1QOXmc5o).
#
# Here we use as example data from an ovarian cancer dataset:
# [PXD010372](https://github.com/Multiomics-Analytics-Group/acore/tree/main/example_data/PXD010372)
#
# First make sure you have the required packages installed:


# %% tags=["hide-output"]
# %pip install acore vuecore 'plotly<6'

# %%
from pathlib import Path

import dsp_pandas
import pandas as pd

import acore
import acore.differential_regulation
import acore.enrichment_analysis

dsp_pandas.format.set_pandas_options(max_colwidth=60)

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
# from our repository. See details on obtaining the data under the example data section on
# [this page](nb_ref_ovarian_data)

# %%
df_omics = pd.read_csv(omics, index_col=0)
df_meta_pgs = pd.read_csv(meta_pgs, index_col=0)
df_meta = pd.read_csv(meta, index_col=0)
df_omics

# %%
ax = (
    df_omics.notna()
    .sum()
    .sort_values(ascending=True)
    .plot(xlabel="Protein groups", ylabel="Number of non-NaN values (samples)")
)

# %% [markdown]
# Keep only features with a certain amount of non-NaN values and select 100 of these
# for illustration. Add always four which were differently regulated in the ANOVA using all
# the protein groups.

# %%
idx_always_included = ["Q5HYN5", "P39059", "O43432", "O43175"]

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

# %% [markdown]
# And we have the following patient metadata, from which we will use the `Status` column as
# our dependent variable and the `PlatinumValue` as a covariate.

# %%
df_meta


# %% [markdown]
# # ANOVA: Compute up and downregulated genes
# These will be used to find enrichments in the set of both up and downregulated genes.

# %%
group = "Status"
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
# # Download functional annotations, here pathways, for the protein groups
# in our selection of the dataset.

# %%
from acore.io.uniprot import fetch_annotations, process_annotations

fname_annotations = f"downloaded/annotations_{features_to_sample}.csv"
fname = Path(fname_annotations)
try:
    annotations = pd.read_csv(fname, index_col=0)
    print(f"Loaded annotations from {fname}")
except FileNotFoundError:
    print(f"Fetching annotations for {df_omics.columns.size} UniProt IDs.")
    FIELDS = "go_p,go_c,go_f"
    annotations = fetch_annotations(df_omics.columns, fields=FIELDS)
    annotations = process_annotations(annotations, fields=FIELDS)
    # cache the annotations
    fname.parent.mkdir(exist_ok=True, parents=True)
    annotations.to_csv(fname, index=True)

annotations

# %% [markdown]
# See how many protein groups are associated with each annotation. We observe that most
# functional annotations are associated only to a single protein group in our dataset.

# %% tags=["hide-input"]
s_count_pg_per_annotation = (
    annotations.groupby("annotation").size().value_counts().sort_index()
)
_ = s_count_pg_per_annotation.plot(
    kind="bar",
    xlabel="Number of protein groups associated with annotation",
    ylabel="Number of annotations",
)
s_count_pg_per_annotation.to_frame("number of annotations").rename_axis(
    "N protein groups"
).T

# %%
annotations.groupby("annotation").size().value_counts(ascending=False)

# %% [markdown]
# # Enrichment analysis
# Is done separately for up- and downregulated genes as it's assumed that biological
# processes are regulated in one direction.

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
# Running the enrichment analysis for the up- and down regulated protein groups
# separately with the default settings of the function, i.e. a log2 fold change cutoff
# of 1 and at least 2 protein groups detected in the set of proteins
# defining the functional annotation.

# %%
ret = acore.enrichment_analysis.run_up_down_regulation_enrichment(
    regulation_data=diff_reg,
    annotation=annotations,
    pval_col="padj",
    min_detected_in_set=2,
    lfc_cutoff=1,
)
ret

# %% [markdown]
# we can decrease the cutoff for the log2 fold change to 0.5 and see that we retain
# more annotations.

# %%
ret = acore.enrichment_analysis.run_up_down_regulation_enrichment(
    regulation_data=diff_reg,
    annotation=annotations,
    pval_col="padj",
    min_detected_in_set=2,
    lfc_cutoff=0.5,  # ! the default is 1
)
ret

# %% [markdown]
# And even more if we do not restrict the analysis of finding at least two proteins
# of a functional set in our data set (i.e. we only need to find one match from the set).

# %%
ret = acore.enrichment_analysis.run_up_down_regulation_enrichment(
    regulation_data=diff_reg,
    annotation=annotations,
    pval_col="padj",
    min_detected_in_set=1,
    lfc_cutoff=0.5,  # ! the default is 1
)
ret

# %% [markdown]
# ## Site specific enrichment analysis

# %% [markdown]
# The basic example uses a modified peptide sequence to
# demonstrate the enrichment analysis.
# > TODO: The example on how to do that needs a PTM focused dataset.
# The details of how site specific enrichment analysis is done will depend on the
# dataset and the question at hand.
#
# If the identifiers contain PTMs this information is removed to match it to the annotation
# using a regular expression (in the function). For example:

# %%
import re

regex = "(\\w+~.+)_\\w\\d+\\-\\w+"
identifier_ckg = "gnd~P00350_T10-WW"
match = re.search(regex, identifier_ckg)
match.group(1)

# %%
# ToDo: Add example for site specific enrichment analysis

# %% [markdown]
# # Single sample GSEA (ssGSEA)
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
ax = enrichtments["NES"].plot.hist()

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
import plotly.graph_objects as go
from vuecore import viz

args = {"factor": 2, "loadings": 1}  # increase number of loadings or scaling factor
# #! pca_results has three items, but docstring requests only two -> double check
figure = viz.get_pca_plot(data=pca_result, identifier="PCA enrichment", args=args)
figure = go.Figure(data=figure["data"], layout=figure["layout"])
figure.show()

# %% [markdown]
# # Compare two distributions - KS test
#
# The Kolmogorov-Smirnov test is a non-parametric test that compares two distributions.
# - we compare the distributions of the two differently upregulated protein groups
# This is not the best example for comparing distributions, but it shows how to use the
# KS test.

# %%
# plot two histograms of intensity values here
sel_pgs = ["O43175", "P39059"]
view = df_omics[sel_pgs].sub(df_omics[sel_pgs].mean())
ax = view.plot.hist(bins=20, alpha=0.5)

# %% [markdown]
# Let us compare the two centered distributions using the KS test.

# %%
acore.enrichment_analysis.run_kolmogorov_smirnov(view[sel_pgs[0]], view[sel_pgs[1]])

# %% [markdown]
# The result suggests that the two distributions are from the same distribution.
