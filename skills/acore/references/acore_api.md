# acore API map

Compact module → function reference for the `acore` package (Analytics Core). Import each
capability from its submodule; `import acore` alone exposes nothing but `__version__`. Most
`run_*` functions return `pandera`-validated DataFrames (schemas in `acore.types`).

## Preprocessing

**`acore.normalization`**
- `normalize_data(data, method="median", normalize=None)` — method ∈ {`median`, `median_polish`,
  `median_zero`, `quantile`, `linear`, `zscore`}. `normalize` is an axis (`"samples"` rows or
  `"features"` columns), not a column subset; it normalizes every numeric column, so pass only the
  feature columns.
- `normalize_data_per_group(data, group, method="median", normalize=None)`
- strategies (`acore.normalization.strategies`): `median_normalization`, `median_zero_normalization`,
  `zscore_normalization`, `median_polish_normalization`, `quantile_normalization`, `linear_normalization`.

**`acore.imputation_analysis`**
- `imputation_KNN(data, drop_cols=None, group=None, cutoff=0.6, alone=True, n_neighbors=3)`
- `imputation_mixed_norm_KNN(...)`
- `imputation_normal_distribution(data, drop_cols=None, shift=1.8, nstd=0.3, random_state=112736)`

**`acore.batch_correction`**
- `combat_batch_correction(data, batch_col)` — ComBat via inmoose.

## Differential regulation / expression — `acore.differential_regulation`

- `run_anova(df, alpha=0.05, drop_cols=None, subject=None, group="group", permutations=0, correction="fdr_bh", is_logged=True, non_par=False)`
- `run_ttest(df, condition1, condition2, alpha=0.05, drop_cols=None, subject="subject", group="group", paired=False, correction="fdr_bh", permutations=0, is_logged=True, non_par=False)`
- `run_ancova(...)`, `run_repeated_measurements_anova(...)`, `run_mixed_anova(...)`,
  `run_two_way_anova(...)`, `run_diff_analysis(...)`
- underlying tests (`acore.differential_regulation.tests`): `calculate_ttest`, `calculate_anova`,
  `calculate_ancova`, `calculate_thsd`, `calculate_pairwise_ttest`,
  `calculate_repeated_measures_anova`, `eta_squared`.

## Multiple testing — `acore.multiple_testing`

- `apply_pvalue_correction(pvalues, alpha=0.05, method="bonferroni")` → `(rejected, adjusted)`.
- `apply_pvalue_fdrcorrection(pvalues, alpha=0.05, method="indep")`
- `apply_pvalue_twostage_fdrcorrection(pvalues, alpha=0.05, method="bh")`
- `apply_pvalue_permutation_fdrcorrection(...)`, `correct_pairwise_ttest(...)`

## Enrichment — `acore.enrichment_analysis`

- `run_enrichment(data, foreground_id, background_id, foreground_pop, background_pop, min_detected_in_set=2, annotation_col="annotation", group_col="group", identifier_col="identifier", method="fisher", correction="fdr_bh", correction_alpha=0.05)`
- `run_regulation_enrichment(...)`, `run_up_down_regulation_enrichment(...)`,
  `run_site_regulation_enrichment(...)`, `run_ssgsea(...)`
- annotations helper: `acore.enrichment_analysis.annotate`.

## Exploratory / decomposition

**`acore.exploratory_analysis`** — functions live in submodules, not the package namespace.
- `decomposition`: `run_pca(data, drop_cols=None, group="group", annotation_cols=None, components=2, dropna=True)`, `run_tsne(...)`, `run_umap(...)`
- `variation`: `calculate_coefficient_variation(df) -> pd.Series`
- `missing`: `extract_number_missing(data, min_valid, drop_cols=None, group="group")`, `extract_percentage_missing(data, missing_max, drop_cols=None, group="group", how="all")`

**`acore.decomposition`** — `pca.py`, `umap.py` (lower-level building blocks).

## Correlation — `acore.correlation_analysis`

- `run_correlation(df, alpha=0.05, subject=None, group="group", method="pearson", correction="fdr_bh", numeric_only=True, dropna=True)`
- `run_multi_correlation(...)`, `run_rm_correlation(df, alpha=0.05, subject="subject", correction="fdr_bh")`,
  `run_efficient_correlation(data, method="pearson")`, `calculate_correlations(x, y, method="pearson")`

## Other domains

- **`acore.kaplan_meier_analysis`** (survival, lifelines): `run_km`, `get_km_results`,
  `get_hazard_ratio_results`.
- **`acore.network_analysis`**: `get_network_communities`, `get_louvain_partitions`,
  `get_snf_clusters`, `run_snf`.
- **`acore.wgcna_analysis`** (PyWGCNA): weighted gene co-expression networks.
- **`acore.tda_analysis`** (kmapper): `run_mapper(...)`.
- **`acore.power_analysis`**: `power_analysis(...)`.
- **`acore.publications_analysis`**: `get_publications_abstracts(...)` (PubMed/Medline).
- **`acore.sklearn`**: `transform_DataFrame(X, fct)` — DataFrame-preserving sklearn helpers.
- **`acore.io`**: `uniprot` (`fetch_annotations`, `process_annotations`), `pride`
  (`download_PRIDE_data`), `ftp`, `http`, `uncompress`.

## Validated output schemas — `acore.types`

`differential_analysis` (e.g. `AnovaSchema`, `AnovaSchemaMultiGroup`), `enrichment_analysis`,
`exploratory_analysis`. Use these to know the exact result columns (`padj`, `log2FC`, `rejected`,
`terms`, etc.).
