# Overview

Complete analyses of a dataset, chaining several acore modules and interpreting the
results. The individual functions are documented in the [API guide](../api/index.md). How
the datasets were obtained and prepared is described under [example
data](../example_data/overview.md).

## Ovarian Cancer Enrichment Analysis

[Enrichment analysis](ovarian_cancer_enrichment.ipynb) of the ovarian cancer proteome
from [Coscia et al. (2018)](https://doi.org/10.1016/j.cell.2018.08.065) (PRIDE
[PXD010372](https://www.ebi.ac.uk/pride/archive/projects/PXD010372)).

- data: [example_data/PXD010372](../example_data/PXD010372) (processed tables in
  [processed/](../example_data/PXD010372/processed))
- preparation: [ovarian cancer download notebook](../example_data/ovarian_cancer.ipynb)

## Phage-resistant _E. coli_ metabolomics

[UMAP and differential analysis](MTBLS13411_metabolomics.ipynb) of the metabolome of
phage-resistant enterotoxigenic _Escherichia coli_ (MetaboLights
[MTBLS13411](https://www.ebi.ac.uk/metabolights/editor/MTBLS13411/files)).

- data: [example_data/MTBLS13411](../example_data/MTBLS13411), with features and metadata
  combined in
  [MTBLS13411_omics_and_meta.csv](../example_data/MTBLS13411/MTBLS13411_omics_and_meta.csv)
- description and loading snippet: [MTBLS13411](../example_data/MTBLS13411.md)
