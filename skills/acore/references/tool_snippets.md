# acore code snippets (one per analysis step)

Each section is a self-contained, runnable snippet for one analysis step, with the exact acore
call, its real defaults, and the input checks worth keeping. There is no acore MCP server: write
and run the code, so the user can read, edit, and re-run it.

Common preamble for every snippet:

```python
import pandas as pd

# Readers by extension. .parquet additionally needs `pyarrow` installed.
READERS = {
    ".csv": pd.read_csv,
    ".tsv": lambda p: pd.read_csv(p, sep="\t"),
    ".xlsx": pd.read_excel,
    ".xls": pd.read_excel,
    ".parquet": pd.read_parquet,
}


def load(path: str) -> pd.DataFrame:
    """Load a tabular omics file by extension."""
    import pathlib

    suffix = pathlib.Path(path).suffix.lower()
    if suffix not in READERS:
        raise ValueError(f"Unsupported extension {suffix!r}; expected {sorted(READERS)}")
    return READERS[suffix](path)
```

Write each intermediate result to disk (`df.to_csv(out, index=False)`) so later steps and the
`vuegen` skill can pick it up.

---

## 1. `normalize_data`

```python
from acore.normalization import normalize_data

METHODS = {"median", "median_polish", "median_zero", "quantile", "linear", "zscore"}

df = load("omics.csv")
method = "median"
assert method in METHODS, f"Invalid method {method!r}; allowed: {sorted(METHODS)}"

# normalize_data transforms every numeric column, including numeric metadata such as
# condition labels, so normalize only the feature columns and keep metadata unchanged.
metadata_cols = ["sample", "subject", "group"]
feature_cols = [c for c in df.columns if c not in metadata_cols]
normalized_features = normalize_data(df[feature_cols], method=method, normalize=None)
norm = pd.concat([df[metadata_cols], normalized_features], axis=1)
norm.to_csv("normalized.csv", index=False)
```

`normalize` is an axis, not a column subset: `None`/`"samples"` normalizes each row,
`"features"` each column (not every method honours it; `quantile`, `median_polish` ignore it).
Select feature columns on the DataFrame, as above. Normalize **before** differential testing.

## 2. `imputation_KNN`

```python
from acore.imputation_analysis import imputation_KNN

df = load("normalized.csv")
imputed = imputation_KNN(
    df,
    drop_cols=None,      # metadata columns to exclude, e.g. ["sample", "subject"]
    group=None,          # impute within group when set
    cutoff=0.6,          # min fraction of valid values: features missing in more than 40% are skipped
    alone=True,
    n_neighbors=3,
)
imputed.to_csv("imputed.csv", index=False)
```

Inspect missingness first (`extract_percentage_missing`). `cutoff` is the minimum fraction of
valid values, not a maximum missing fraction: with `cutoff=0.6` a feature missing in more than
40% of samples (per group, when `group` is set) is not imputed, and with `alone=True` it is then
dropped.

## 3. `run_anova`

```python
from acore.differential_regulation import run_anova

df = load("imputed.csv")
group = "group"
assert group in df.columns, f"Missing required group column: {group}"

res = run_anova(
    df,
    alpha=0.05,
    drop_cols=["sample", "subject"],   # REQUIRED: acore defaults to None (drops nothing)
    subject=None,                      # name the column for paired/repeated-measures designs;
                                       # a non-None subject that is not a column raises KeyError
    group=group,
    permutations=0,
    correction="fdr_bh",
    is_logged=True,                    # WRONG value silently corrupts log2FC
    non_par=False,
)
res.to_csv("anova.csv", index=False)
res.nsmallest(10, "padj")
```

## 4. `run_ttest` (two-group comparison)

```python
from acore.differential_regulation import run_ttest

res = run_ttest(
    load("imputed.csv"),
    condition1="control",
    condition2="treated",
    alpha=0.05,
    drop_cols=["sample"],              # REQUIRED: acore defaults to None (drops nothing).
                                       # Keep `subject` out: run_ttest drops it itself, and a
                                       # second drop raises KeyError
    subject=None,                      # unpaired; for paired designs pass the subject column
                                       # name and set paired=True
    group="group",
    paired=False,
    correction="fdr_bh",
    permutations=0,
    is_logged=True,
)
```

## 5. `apply_pvalue_correction`

```python
import numpy as np
from acore.multiple_testing import apply_pvalue_correction

res = load("anova.csv")
pvalues = np.asarray(res["pvalue"], dtype=float)   # None/missing -> np.nan

rejected, adjusted = apply_pvalue_correction(
    pvalues,
    alpha=0.05,
    method="fdr_bh",   # acore defaults to "bonferroni"; prefer fdr_bh
)
res["rejected"], res["padj"] = rejected, adjusted
# `rejected` is a float array with NaN where the p-value was missing, so count NaN-aware.
n_rejected = int(np.nansum(rejected))
print(f"{n_rejected} of {int(np.isfinite(pvalues).sum())} tested rejected at alpha=0.05")
```

Report `padj`, never raw `pvalue`.

## 6. `uniprot.fetch_annotations`

Performs a network request to the UniProt API.

```python
from acore.io import uniprot

ids = ["P01234", "Q99999"]     # must be a non-empty list
assert ids, "ids must be a non-empty list of UniProt IDs"

fields = "go_p,go_c,go_f"
annotations = uniprot.fetch_annotations(ids, fields=fields)
# Long-format (identifier, annotation, source) table the enrichment steps expect:
annotations = uniprot.process_annotations(annotations, fields=fields)
annotations.to_csv("annotations.csv", index=False)
```

Skip `process_annotations` only if you want the raw UniProt response.

## 7. `run_up_down_regulation_enrichment`

```python
from acore.enrichment_analysis import run_up_down_regulation_enrichment

enr = run_up_down_regulation_enrichment(
    regulation_data=load("anova.csv"),
    annotation=load("annotations.csv"),
    identifier="identifier",
    groups=("group1", "group2"),   # NAMES of columns in the regulation table (the
                                   # anova.csv columns identifying the compared groups),
                                   # not the condition values like "control"/"treated"
    annotation_col="annotation",
    pval_col="padj",
    group_col="group",
    log2fc_col="log2FC",
    method="fisher",
    min_detected_in_set=2,
    correction="fdr_bh",
    correction_alpha=0.05,
    lfc_cutoff=1.0,
)
enr.to_csv("enrichment.csv", index=False)
```

## 8. `run_enrichment` (over-representation)

```python
from acore.enrichment_analysis import run_enrichment

from acore.enrichment_analysis.annotate import annotate_features

annotations = load("annotations.csv")   # columns: identifier, source, annotation

# annotations.csv has no foreground/background membership; run_enrichment groups by
# `group_col`, so assign it from the selected feature sets first.
foreground_ids = ["P01234"]                       # e.g. significant proteins from anova.csv
background_ids = ["P01234", "Q99999", "P05067"]   # e.g. every protein that was tested
annotations["group"] = annotate_features(
    annotations["identifier"], foreground_ids, background_ids
)
annotations = annotations.dropna(subset=["group"])

res = run_enrichment(
    annotations,
    foreground_id="foreground",
    background_id="background",
    foreground_pop=len(foreground_ids),
    background_pop=len(background_ids),
    min_detected_in_set=2,
    annotation_col="annotation",
    group_col="group",
    identifier_col="identifier",
    method="fisher",
    correction="fdr_bh",
    correction_alpha=0.05,
)
```

## 9. `run_correlation`

```python
from acore.correlation_analysis import run_correlation

res = run_correlation(
    load("normalized.csv"),
    alpha=0.05,
    subject=None,          # set for repeated measures
    group="group",
    method="pearson",      # or "spearman"
    correction="fdr_bh",
    numeric_only=True,
    dropna=True,
)
res.to_csv("correlation.csv", index=False)
```

## 10. `combat_batch_correction`

```python
from acore.batch_correction import combat_batch_correction

df = load("normalized.csv")
batch_col = "batch"
assert isinstance(batch_col, str) and batch_col, "batch_col must be a non-empty string"
assert batch_col in df.columns, f"Missing required batch column: {batch_col}"
assert df[batch_col].nunique() > 1, "ComBat needs more than one batch"

corrected = combat_batch_correction(df, batch_col)
corrected.to_csv("batch_corrected.csv", index=False)
```

Check a PCA colored by batch before and after.

## 11. `calculate_coefficient_variation`

```python
from acore.exploratory_analysis.variation import calculate_coefficient_variation

df = load("normalized.csv")
drop_cols = ["sample", "subject", "group"]
df = df.drop(columns=[c for c in drop_cols if c in df.columns], errors="ignore")

numeric = df.select_dtypes("number")
assert not numeric.empty, "No numeric columns to compute coefficient of variation over"

cv = calculate_coefficient_variation(numeric).reset_index()
cv.columns = ["feature", "coef_of_var"]     # CV in %
cv.to_csv("cv.csv", index=False)
```

## 12. `get_enrichment_plots` (needs `vuecore` — see the `vuecore` skill)

```python
from vuecore.viz import get_enrichment_plots

enr = load("enrichment.csv")

# A CSV round-trip can turn the boolean `rejected` flag into strings, and
# get_enrichment_plots filters on it as a boolean mask.
if "rejected" in enr.columns and enr["rejected"].dtype == object:
    enr["rejected"] = enr["rejected"].astype(str).str.strip().str.lower() == "true"

figures = get_enrichment_plots(
    enr,
    identifier="enrichment",
    args={"title": "Enrichment Analysis", "width": 900, "height": 800},
)
for i, fig in enumerate(figures):
    fig.write_html(f"enrichment_{i}.html")   # or .write_image(...) for static output
```

---

## End-to-end example

acore's `example_data/` sits at the root of the **acore repository** and is not shipped inside the
installed package, so read it over HTTP (or point `load()` at the user's own file):

```python
import pandas as pd
from acore.normalization import normalize_data
from acore.differential_regulation import run_anova

URL = (
    "https://raw.githubusercontent.com/Multiomics-Analytics-Group/acore/main/"
    "example_data/alzheimer_proteomics/alzheimer_example_omics_and_clinic.csv"
)
df = pd.read_csv(URL)   # columns: Sample ID, AD, age, male, collection_site, <UniProt ids...>

# Normalize only the protein columns; normalize_data would also rewrite numeric
# metadata such as the `AD` labels.
metadata_cols = ["Sample ID", "AD", "age", "male", "collection_site"]
protein_cols = [c for c in df.columns if c not in metadata_cols]
normalized_proteins = normalize_data(df[protein_cols], method="median")
norm = pd.concat([df[metadata_cols], normalized_proteins], axis=1)
res = run_anova(
    norm,
    group="AD",                                             # this dataset's condition column
    drop_cols=["Sample ID", "age", "male", "collection_site"],
    subject=None,                                           # no subject column in this dataset
    correction="fdr_bh",
    is_logged=True,
)
res.nsmallest(10, "padj")
```

For the long tail acore offers beyond these steps (survival, networks, WGCNA, power analysis,
ssGSEA, ANCOVA/mixed ANOVA, decomposition), see [`acore_api.md`](acore_api.md).
