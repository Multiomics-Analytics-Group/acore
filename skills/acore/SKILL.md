---
name: acore
description: "Use for multi-omics preprocessing and statistical analysis with acore: normalization, imputation, differential expression (ANOVA/t-test), multiple-testing correction, correlation, PCA/UMAP, batch correction, and enrichment. Teaches the acore Python API code-first, with a runnable snippet per analysis step."
---

# acore skill

`acore` (Analytics Core) does preprocessing and statistical testing of multi-omics tabular
data. Write and run Python that calls it directly; there is no acore server or CLI.

Install it into the project's own environment and declare it in the project's manifest. Verify
with the submodule you need, e.g. `import acore.normalization`. The first import can take over
two minutes while numba/umap warm up; that is not a hang.

A ready snippet for every analysis step is in
[`references/tool_snippets.md`](references/tool_snippets.md). Start from the snippet for the step
you need, adapt the column names to the user's data, and run it.

## Critical gotcha: import submodules explicitly

`import acore` exposes almost nothing — the top-level `__init__` only sets pandas formatting.
Every capability lives in a submodule you must import by name:

```python
from acore.normalization import normalize_data
from acore.imputation_analysis import imputation_KNN
from acore.differential_regulation import run_anova, run_ttest
from acore.multiple_testing import apply_pvalue_correction
```

## DataFrame conventions (read before writing any call)

- **Shape**: wide numeric matrix (rows = samples, columns = features), plus metadata columns.
  Most `run_*` functions expect metadata columns named `sample`, `subject`, and `group`.
- **`drop_cols`**: non-numeric/metadata columns to exclude. Defaults to `None` (drop nothing) —
  **always pass it explicitly**, or metadata columns are treated as features. The `group` column
  is named separately and does not go in `drop_cols`.
- **`group`**: the column holding the condition/class label used for comparisons.
- **`correction`**: multiple-testing method; default `"fdr_bh"`. Keep it unless there is a reason.
- **`is_logged`**: whether values are already log-transformed. Differential functions assume
  `True`; setting it wrong silently corrupts log2FC.
- **Outputs** of most `run_*` functions are `pandera`-validated DataFrames (schemas in
  `acore.types`).

## The standard analysis pipeline

Load → normalize → impute → test → correct → enrich. Snippets are in
[`references/tool_snippets.md`](references/tool_snippets.md).

| Step | acore code | Snippet |
|---|---|---|
| Normalize | `normalization.normalize_data(df, method="median")` | §1 |
| Impute (KNN) | `imputation_analysis.imputation_KNN(df, ...)` | §2 |
| ANOVA | `differential_regulation.run_anova(df, group="group")` | §3 |
| t-test | `differential_regulation.run_ttest(df, c1, c2)` | §4 |
| p-value correction | `multiple_testing.apply_pvalue_correction(pvals, method="fdr_bh")` | §5 |
| UniProt annotations | `io.uniprot.fetch_annotations(ids, ...)` | §6 |
| Up/down enrichment | `enrichment_analysis.run_up_down_regulation_enrichment(...)` | §7 |
| Over-representation enrichment | `enrichment_analysis.run_enrichment(...)` | §8 |
| Correlation | `correlation_analysis.run_correlation(df, method="pearson")` | §9 |
| Batch correction | `batch_correction.combat_batch_correction(df, batch_col)` | §10 |
| Coefficient of variation | `exploratory_analysis.calculate_coefficient_variation(df)` | §11 |
| Enrichment plot | `vuecore.viz.get_enrichment_plots(...)` | §12 |
| PCA | `exploratory_analysis.run_pca(df, group="group")` | see `acore_api.md` |

Two defaults bite on data that isn't in acore's canonical shape: `drop_cols=None` drops nothing,
and `subject="subject"` makes `run_anova`/`run_ttest` group on a `subject` column — a `KeyError`
if the frame has none. Pass `subject=None` for unpaired data.

A fuller module→function map (survival, networks, WGCNA, power, ssGSEA, ANCOVA/mixed ANOVA,
decomposition) is in [`references/acore_api.md`](references/acore_api.md).

## Quick demo

Pulls the repository's example proteomics data, normalizes it, runs an ANOVA and prints the ten
most significant proteins. No user data is touched. `example_data/` is in the repository, not the
installed package, so read it over HTTP (pin the ref to the installed version's tag):

```python
import pandas as pd
from acore.normalization import normalize_data
from acore.differential_regulation import run_anova

URL = (
    "https://raw.githubusercontent.com/Multiomics-Analytics-Group/acore/main/"
    "example_data/alzheimer_proteomics/alzheimer_example_omics_and_clinic.csv"
)
df = pd.read_csv(URL)

# The condition column is `AD`; the rest of the metadata is clinical and there is no
# `subject` column — hence subject=None. Adapt all four to the user's data.
norm = normalize_data(df, method="median")
res = run_anova(
    norm,
    group="AD",
    drop_cols=["Sample ID", "age", "male", "collection_site"],
    subject=None,
    correction="fdr_bh",
    is_logged=True,
)
res.nsmallest(10, "padj")
```

## Reasoning guardrails (apply these before running)

1. **Normalize + log before differential testing.** `run_anova`/`run_ttest` assume roughly
   normalized, log-scaled intensities (`is_logged=True`). A wrong `is_logged` silently corrupts
   fold-changes.
2. **Handle missingness before imputation and before tests that dislike NaNs.** Inspect
   `extract_percentage_missing` / `calculate_coefficient_variation`; impute only features below
   the missingness `cutoff`.
3. **Keep FDR correction on.** Report `padj`, not raw `pvalue`.
4. **Warn on small n.** With few replicates per `group`, p-values are underpowered; effect sizes
   are the more trustworthy signal.
5. **Batch effects.** If a batch/run column exists, consider `combat_batch_correction` before
   differential testing, and check PCA colored by batch.

## Escape hatch

`acore` covers more than the snippets (decomposition, survival via `kaplan_meier_analysis`,
`network_analysis`, `wgcna_analysis`, `power_analysis`, ssGSEA). Look the function up in
[`references/acore_api.md`](references/acore_api.md) and write the call directly rather than
bending a snippet into a shape it does not fit.
