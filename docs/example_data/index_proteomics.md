# Proteomics

## Download PRIDE data

Download results from PRIDE for project
[PXD008541](https://www.ebi.ac.uk/pride/archive/projects/PXD008541) on human primary
brown and white fat cell secretome.

## Ovarian cancer proteomics reanalysis

Download the ovarian cancer proteome data from the journal's website. It was provided as
supplementary data. See the article here:

> Fabian Coscia, Ernst Lengyel, Jaikumar Duraiswamy, Bradley Ashcroft, Michal
> Bassani-Sternberg, Michael Wierer, Alyssa Johnson, Kristen Wroblewski, Anthony Montag,
> S. Diane Yamada, Blanca López-Méndez, Jakob Nilsson, Andreas Mund, Matthias Mann,
> Marion Curtis, Multi-level Proteomics Identifies CT45 as a Chemosensitivity Mediator
> and Immunotherapy Target in Ovarian Cancer, Cell, Volume 175, Issue 1, 2018,
>
> [doi.org/10.1016/j.cell.2018.08.065](https://doi.org/10.1016/j.cell.2018.08.065)

## Alzheimer's proteomics reanalysis

```{include} ../../example_data/alzheimer_proteomics/README.md
:start-line: 5
```

You can use the dumped file with metadata and a sample of 100 proteins from the study
directly loading it from the acore repository using

```python
import pandas as pd

df = pd.read_csv(
    "https://raw.githubusercontent.com/Multiomics-Analytics-Group/acore/refs/heads/main/"
    "example_data/alzheimer_proteomics/alzheimer_example_omics_and_clinic.csv",
    index_col="Sample ID"
)

meta = df.iloc[:,:4]
omics= df.drop(meta.columns, axis=1)
```

```{toctree}
:maxdepth: 1
:hidden:

Download_PRIDE_data.ipynb
ovarian_cancer.ipynb
alzheimer_proteomics.ipynb
```
