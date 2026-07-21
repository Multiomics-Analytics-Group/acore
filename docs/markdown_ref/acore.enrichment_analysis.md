# acore.enrichment_analysis package

Enrichment Analysis Module. Contains different functions to perform enrichment
analysis.

Most things in this module are covered in [https://www.youtube.com/watch?v=2NC1QOXmc5o](https://www.youtube.com/watch?v=2NC1QOXmc5o)
by Lars Juhl Jensen.

### run_site_regulation_enrichment(regulation_data: DataFrame, annotation: DataFrame, identifier: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'identifier', groups: [list](https://docs.python.org/3/library/stdtypes.html#list)[[str](https://docs.python.org/3/library/stdtypes.html#str)] = ('group1', 'group2'), annotation_col: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'annotation', rejected_col: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'rejected', group_col: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'group', method: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'fisher', regex: [str](https://docs.python.org/3/library/stdtypes.html#str) = '(\\\\w+~.+)_\\\\w\\\\d+\\\\-\\\\w+', correction: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'fdr_bh', remove_duplicates: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[[EnrichmentAnalysisSchema](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema)]

This function runs a simple enrichment analysis for significantly
regulated protein sites in a dataset.

* **Parameters:**
  * **regulation_data** – pandas.DataFrame resulting from differential
    regulation analysis.
  * **annotation** – pandas.DataFrame with annotations for features
    (columns: ‘annotation’, ‘identifier’ (feature identifiers), and ‘source’).
  * **identifier** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – name of the column from annotation containing
    feature identifiers.
  * **groups** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – column names from regulation_data containing
    group identifiers.
  * **annotation_col** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – name of the column from annotation containing
    annotation terms.
  * **rejected_col** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – name of the column from regulation_data containing
    boolean for rejected null hypothesis.
  * **group_col** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column name for new column in annotation dataframe
    determining if feature belongs to foreground or background.
  * **method** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – method used to compute enrichment
    (only ‘fisher’ is supported currently).
  * **regex** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – how to extract the annotated identifier from the site identifier
* **Returns:**
  pandas.DataFrame with columns: ‘terms’, ‘identifiers’, ‘foreground’,
  ‘background’, foreground_pop, background_pop, ‘pvalue’, ‘padj’ and ‘rejected’.
* **Raises:**
  [**ValueError**](https://docs.python.org/3/library/exceptions.html#ValueError) – if regulation_data is None or empty.

Example:

```default
result = run_site_regulation_enrichment(regulation_data,
    annotation,
    identifier='identifier',
    groups=['group1', 'group2'],
    annotation_col='annotation',
    rejected_col='rejected',
    group_col='group',
    method='fisher',
    match="(\\w+~.+)_\\w\\d+\\-\\w+"
)
```

### run_up_down_regulation_enrichment(regulation_data: DataFrame, annotation: DataFrame, identifier: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'identifier', groups: [list](https://docs.python.org/3/library/stdtypes.html#list)[[str](https://docs.python.org/3/library/stdtypes.html#str)] = ('group1', 'group2'), annotation_col: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'annotation', pval_col: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'pval', group_col: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'group', log2fc_col: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'log2FC', method: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'fisher', min_detected_in_set: [int](https://docs.python.org/3/library/functions.html#int) = 2, correction: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'fdr_bh', correction_alpha: [float](https://docs.python.org/3/library/functions.html#float) = 0.05, lfc_cutoff: [float](https://docs.python.org/3/library/functions.html#float) = 1) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[[EnrichmentAnalysisSchema](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema)]

This function runs a simple enrichment analysis for significantly regulated proteins
distinguishing between up- and down-regulated.

* **Parameters:**
  * **regulation_data** ([*pandas.DataFrame*](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html#pandas.DataFrame)) – pandas.DataFrame resulting from differential regulation
    analysis (CKG’s regulation table).
  * **annotation** ([*pandas.DataFrame*](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html#pandas.DataFrame)) – pandas.DataFrame with annotations for features
    (columns: ‘annotation’, ‘identifier’ (feature identifiers), and ‘source’).
  * **identifier** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – name of the column from annotation containing feature identifiers.
  * **groups** ([*list*](https://docs.python.org/3/library/stdtypes.html#list) *[*[*str*](https://docs.python.org/3/library/stdtypes.html#str) *]*) – 

    column names from regulation_data containing group identifiers.
    See [pandas.DataFrame.groupby](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.groupby.html) for more information.
  * **annotation_col** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – name of the column from annotation containing annotation terms.
  * **rejected_col** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – name of the column from regulation_data containing boolean for
    rejected null hypothesis.
  * **group_col** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column name for new column in annotation dataframe determining
    if feature belongs to foreground or background.
  * **method** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – method used to compute enrichment
    (only ‘fisher’ is supported currently).
  * **correction** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – method to be used for multiple-testing correction
  * **alpha** ([*float*](https://docs.python.org/3/library/functions.html#float)) – adjusted p-value cutoff to define significance
  * **lfc_cutoff** ([*float*](https://docs.python.org/3/library/functions.html#float)) – log fold-change cutoff to define practical significance
* **Returns:**
  DataFrame adhering to EnrichmentAnalysisSchema
* **Return type:**
  DataFrame[[EnrichmentAnalysisSchema](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema)]

Example:

```default
result = run_up_down_regulation_enrichment(
    regulation_data,
    annotation,
    identifier='identifier',
    groups=['group1',
    'group2'],
    annotation_col='annotation',
    rejected_col='rejected',
    group_col='group',
    method='fisher',
    correction='fdr_bh',
    alpha=0.05,
    lfc_cutoff=1,
)
```

### run_fisher(group1: [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)], group2: [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)], alternative: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'two-sided') → [tuple](https://docs.python.org/3/library/stdtypes.html#tuple)[[float](https://docs.python.org/3/library/functions.html#float), [float](https://docs.python.org/3/library/functions.html#float)]

Run fisher’s exact test on two groups using [scipy.stats.fisher_exact](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html).

Example:

```default
# annotated   not-annotated
# group1      a               b
# group2      c               d


odds, pvalue = stats.fisher_exact(group1=[a, b],
                                  group2 =[c, d]
                )
```

### run_kolmogorov_smirnov(dist1: [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)], dist2: [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)], alternative: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'two-sided') → [tuple](https://docs.python.org/3/library/stdtypes.html#tuple)[[float](https://docs.python.org/3/library/functions.html#float), [float](https://docs.python.org/3/library/functions.html#float)]

Compute the Kolmogorov-Smirnov statistic on 2 samples.
See [scipy.stats.ks_2samp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ks_2samp.html)

* **Parameters:**
  * **dist1** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – sequence of 1-D ndarray (first distribution to compare)
    drawn from a continuous distribution
  * **dist2** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – sequence of 1-D ndarray (second distribution to compare)
    drawn from a continuous distribution
  * **alternative** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – defines the alternative hypothesis (default is ‘two-sided’):
    \* **‘two-sided’**
    \* **‘less’**
    \* **‘greater’**
* **Returns:**
  statistic float and KS statistic pvalue float Two-tailed p-value.

Example:

```default
result = run_kolmogorov_smirnov(dist1, dist2, alternative='two-sided')
```

## Subpackages

* [acore.enrichment_analysis.statistical_tests namespace](acore.enrichment_analysis.statistical_tests.md)
  * [Submodules](acore.enrichment_analysis.statistical_tests.md#submodules)
  * [acore.enrichment_analysis.statistical_tests.fisher module](acore.enrichment_analysis.statistical_tests.md#module-acore.enrichment_analysis.statistical_tests.fisher)
    * [`run_fisher()`](acore.enrichment_analysis.statistical_tests.md#acore.enrichment_analysis.statistical_tests.fisher.run_fisher)
  * [acore.enrichment_analysis.statistical_tests.kolmogorov_smirnov module](acore.enrichment_analysis.statistical_tests.md#module-acore.enrichment_analysis.statistical_tests.kolmogorov_smirnov)
    * [`run_kolmogorov_smirnov()`](acore.enrichment_analysis.statistical_tests.md#acore.enrichment_analysis.statistical_tests.kolmogorov_smirnov.run_kolmogorov_smirnov)

## Submodules

## acore.enrichment_analysis.annotate module

Put unique features into foreground, background or assign nan.

### annotate_features(features: Series, in_foreground: [set](https://docs.python.org/3/library/stdtypes.html#set)[[str](https://docs.python.org/3/library/stdtypes.html#str)] | [list](https://docs.python.org/3/library/stdtypes.html#list)[[str](https://docs.python.org/3/library/stdtypes.html#str)], in_background: [set](https://docs.python.org/3/library/stdtypes.html#set)[[str](https://docs.python.org/3/library/stdtypes.html#str)] | [list](https://docs.python.org/3/library/stdtypes.html#list)[[str](https://docs.python.org/3/library/stdtypes.html#str)]) → Series

Annotate features as foreground or background based on their presence in the
foreground and background lists.

* **Parameters:**
  * **features** – pandas.Series with features and their annotations.
  * **in_foreground** ([*set*](https://docs.python.org/3/library/stdtypes.html#set) *or* *list-like*) – list of features identifiers in the foreground.
  * **in_background** ([*set*](https://docs.python.org/3/library/stdtypes.html#set) *or* *list-like*) – list of features identifiers in the background.
* **Returns:**
  pandas.Series containing ‘foreground’ or ‘background’.
  missing values are preserved.

Example:

```default
result = _annotate_features(features, in_foreground, in_background)
```
