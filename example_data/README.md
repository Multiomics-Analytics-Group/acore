# Example datasets

The example datasets are not part of the Python package `acore`, but are used in the
documentation and examples. This is done by downloading the data from GitHub directly.
Folders with a README (linked below) describe the origin and processing of the data.

## Proteomics

- [PXD008541](PXD008541/README.md) - human primary brown and white fat cell secretome
  (used in the `Download_PRIDE_data` recipe)
- [PXD010372](PXD010372/README.md) - ovarian cancer proteome
  ([article](https://www.sciencedirect.com/science/article/pii/S0092867418311668),
  taken from supplementary tables); used in the enrichment examples. The MaxQuant output
  and sample annotations are also on PRIDE, see the
  [download notebook](../docs/example_data/ovarian_cancer.py)
- [PXD016278](PXD016278/README.md) - Alzheimer's disease cerebrospinal fluid proteome
  ([publication](https://www.embopress.org/doi/full/10.15252/msb.20199356)), curated
  version from OmicLearn
- [alzheimer_proteomics](alzheimer_proteomics/README.md) - processed Alzheimer's
  dataset, see the
  [processing notebook](../docs/example_data/alzheimer_proteomics.py)
- [syntethic_pep_enrichment](syntethic_pep_enrichment/README.md) - synthetic
  peptide-level differential analysis with protein-level pathway annotations (enrichment)

## Metabolomics and lipidomics

- [MTBLS13411](MTBLS13411/README.md) - metabolomics of phage resistance in
  enterotoxigenic *Escherichia coli* (used in the exploratory and differential
  regulation examples)
- [MTBLS733](MTBLS733/README.md) - untargeted metabolomics benchmark feature table
- [aradopsis_seedling_lipids](aradopsis_seedling_lipids/README.md) - untargeted
  Arabidopsis seedling lipids (three conditions, with QC samples)
- `DidacMauricio_hilic` (no README) - HILIC positive mode metabolomics, including
  injection order (e.g. for drift correction)

## Metagenomics

- `mgnify` (no README) - paired functional annotation data from
  [MGnify](https://www.ebi.ac.uk/metagenomics/), a resource for the analysis and
  archiving of metagenomic data
