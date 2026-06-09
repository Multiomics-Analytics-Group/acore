# acore.differential_regulation package

Differential regulation module.

### run_anova(df: DataFrame, alpha: [float](https://docs.python.org/3/library/functions.html#float) = 0.05, drop_cols: [list](https://docs.python.org/3/library/stdtypes.html#list)[[str](https://docs.python.org/3/library/stdtypes.html#str)] = ['sample', 'subject'], subject: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'subject', group: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'group', permutations: [int](https://docs.python.org/3/library/functions.html#int) = 0, correction: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'fdr_bh', is_logged: [bool](https://docs.python.org/3/library/functions.html#bool) = True, non_par: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[[AnovaSchema](acore.types.md#acore.types.differential_analysis.AnovaSchema)] | [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[[AnovaSchemaMultiGroup](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup)]

Performs statistical test for each protein in a dataset.
Checks what type of data is the input (paired, unpaired or repeated measurements) and
performs posthoc tests for multiclass data (i.e., when there are more than two groups,
posthoc tests such as pairwise t-tests or Tukey’s HSD are used to determine which specific
groups differ after finding a significant overall effect).
Multiple hypothesis correction uses permutation-based
if permutations>0 and Benjamini/Hochberg if permutations=0.

* **Parameters:**
  * **df** (*pd.DataFrame*) – pandas dataframe with samples as rows and protein identifiers as columns
    (with additional columns ‘group’, ‘sample’ and ‘subject’).
  * **alpha** ([*float*](https://docs.python.org/3/library/functions.html#float)) – error rate for multiple hypothesis correction
  * **drop_cols** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – column labels to be dropped from the dataframe
  * **subject** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column with subject identifiers
  * **group** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column with group identifiers
  * **permutations** ([*int*](https://docs.python.org/3/library/functions.html#int)) – number of permutations used to estimate false discovery rates.
  * **correction** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – method of pvalue correction see apply_pvalue_correction for methods,
    use methods available in acore.multiple_testing
  * **is_logged** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – whether data is log-transformed
  * **non_par** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – if True, normality and variance equality assumptions are checked
    and non-parametric test Mann Whitney U test if not passed
* **Returns:**
  DataFrame adhering to AnovaSchema or AnovaSchemaMultiGroup.
* **Return type:**
  DataFrame[[AnovaSchema](acore.types.md#acore.types.differential_analysis.AnovaSchema)] | DataFrame[[AnovaSchemaMultiGroup](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup)]

Example:

```default
result = run_anova(df,
                   alpha=0.05,
                   drop_cols=["sample",'subject'],
                   subject='subject',
                   group='group',
                   permutations=50
        )
```

### run_ancova(df: DataFrame, covariates: [list](https://docs.python.org/3/library/stdtypes.html#list)[[str](https://docs.python.org/3/library/stdtypes.html#str)], alpha: [float](https://docs.python.org/3/library/functions.html#float) = 0.05, drop_cols: [list](https://docs.python.org/3/library/stdtypes.html#list)[[str](https://docs.python.org/3/library/stdtypes.html#str)] = ['sample', 'subject'], subject: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'subject', group: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'group', permutations: [int](https://docs.python.org/3/library/functions.html#int) = 0, correction: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'fdr_bh', is_logged: [bool](https://docs.python.org/3/library/functions.html#bool) = True, non_par: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[[AncovaSchema](acore.types.md#acore.types.differential_analysis.AncovaSchema)]

Performs statistical test for each protein in a dataset.
Checks what type of data is the input (paired, unpaired or repeated measurements)
and performs posthoc tests for multiclass data.
Multiple hypothesis correction uses permutation-based
if permutations>0 and Benjamini/Hochberg if permutations=0.

* **Parameters:**
  * **df** (*pd.DataFrame*) – Pandas DataFrame with samples as rows and protein identifiers and
    covariates as columns (with additional columns ‘group’, ‘sample’ and ‘subject’).
  * **covariates** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – list of covariates to include in the model (column in df)
  * **alpha** ([*float*](https://docs.python.org/3/library/functions.html#float)) – error rate for multiple hypothesis correction
  * **drop_cols** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – column labels to be dropped from the DataFrame
  * **subject** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column with subject identifiers
  * **group** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column with group identifiers
  * **permutations** ([*int*](https://docs.python.org/3/library/functions.html#int)) – number of permutations used to estimate false discovery rates.
  * **correction** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – method of pvalue correction see apply_pvalue_correction for methods,
    use methods available in acore.multiple_testing
  * **is_logged** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – whether data is log-transformed
  * **non_par** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – if True, normality and variance equality assumptions are checked
    and non-parametric test Mann Whitney U test if not passed
* **Returns:**
  DataFrame adhering to AncovaSchema
* **Return type:**
  DataFrame[[AncovaSchema](acore.types.md#acore.types.differential_analysis.AncovaSchema)]

Example:

```default
result = run_ancova(df,
                    covariates=['age'],
                    alpha=0.05,
                    drop_cols=["sample",'subject'],
                    subject='subject',
                    group='group',
                    permutations=50
        )
```

### run_diff_analysis(df: DataFrame, boolean_array: Series, event_names: [tuple](https://docs.python.org/3/library/stdtypes.html#tuple)[[str](https://docs.python.org/3/library/stdtypes.html#str), [str](https://docs.python.org/3/library/stdtypes.html#str)] = ('1', '0'), ttest_vars=('alternative', 'p-val', 'cohen-d')) → DataFrame

Differential analysis procedure between two groups. Calculaes
mean per group and t-test for each variable in vars between two groups.

### run_mixed_anova(df, alpha=0.05, drop_cols=['sample'], subject='subject', within='group', between='group2', correction='fdr_bh')

In statistics, a mixed-design analysis of variance model, also known as a split-plot
ANOVA, is used to test
for differences between two or more independent groups whilst subjecting participants
to repeated measures.
Thus, in a mixed-design ANOVA model, one factor (a fixed effects factor) is a
between-subjects variable and the other
(a random effects factor) is a within-subjects variable. Thus, overall, the model is a
type of mixed-effects model ([source](https://en.wikipedia.org/wiki/Mixed-design_analysis_of_variance))

* **Parameters:**
  * **df** (*pd.DataFrame*) – Pandas DataFrame with samples as rows and protein identifiers as columns
    (with additional columns ‘group’, ‘sample’ and ‘subject’).
  * **alpha** ([*float*](https://docs.python.org/3/library/functions.html#float)) – error rate for multiple hypothesis correction
  * **drop_cols** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – column labels to be dropped from the DataFrame
  * **subject** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column with subject identifiers
  * **within** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column with within factor identifiers
  * **between** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column with between factor identifiers
  * **correction** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – method of pvalue correction see apply_pvalue_correction for methods,
    use methods available in acore.multiple_testing
* **Returns:**
  Pandas DataFrame
* **Return type:**
  pd.DataFrame

Example:

```default
result = run_mixed_anova(df,
                         alpha=0.05,
                         drop_cols=['sample'],
                         subject='subject',
                         within='group',
                         between='group2',
        )
```

### run_repeated_measurements_anova(df, alpha=0.05, drop_cols=['sample'], subject='subject', within='group', permutations=50, correction='fdr_bh', is_logged=True) → DataFrame

Performs repeated measurements anova and pairwise posthoc tests for each protein in dataframe.

* **Parameters:**
  * **df** (*pd.DataFrame*) – Pandas DataFrame with samples as rows and protein identifiers as columns
    (with additional columns ‘group’, ‘sample’ and ‘subject’).
  * **alpha** ([*float*](https://docs.python.org/3/library/functions.html#float)) – error rate for multiple hypothesis correction
  * **drop_cols** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – column labels to be dropped from the DataFrame
  * **subject** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column with subject identifiers
  * **within** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column with within factor identifiers
  * **permutations** ([*int*](https://docs.python.org/3/library/functions.html#int)) – number of permutations used to estimate false discovery rates
  * **correction** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – method of pvalue correction see apply_pvalue_correction for methods,
    use methods available in acore.multiple_testing
  * **is_logged** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – whether data is log-transformed
* **Returns:**
  Pandas DataFrame

Example:

```default
result = run_repeated_measurements_anova(df,
                                         alpha=0.05,
                                         drop_cols=['sample'],
                                         subject='subject',
                                         within='group',
                                         permutations=50
        )
```

### run_ttest(df, condition1, condition2, alpha=0.05, drop_cols=['sample'], subject='subject', group='group', paired=False, correction='fdr_bh', permutations=0, is_logged=True, non_par=False)

Runs t-test (paired/unpaired) for each protein in dataset and performs
permutation-based (if permutations>0) or Benjamini/Hochberg (if permutations=0)
multiple hypothesis correction.

* **Parameters:**
  * **df** (*pd.DataFrame*) – Pandas DataFrame with samples as rows and protein identifiers as columns
    (with additional columns ‘group’, ‘sample’ and ‘subject’).
  * **condition1** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – first of two conditions of the independent variable
  * **condition2** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – second of two conditions of the independent variable
  * **alpha** ([*float*](https://docs.python.org/3/library/functions.html#float)) – error rate for multiple hypothesis correction
  * **drop_cols** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – column labels to be dropped from the DataFrame
  * **subject** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column with subject identifiers
  * **group** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column with group identifiers (independent variable)
  * **paired** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – paired or unpaired samples
  * **correction** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – method of pvalue correction see apply_pvalue_correction for methods
  * **permutations** ([*int*](https://docs.python.org/3/library/functions.html#int)) – number of permutations used to estimate false discovery rates.
  * **is_logged** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – data is log-transformed
  * **non_par** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – if True, normality and variance equality assumptions are checked
    and non-parametric test Mann Whitney U test if not passed
* **Returns:**
  Pandas DataFrame with columns ‘identifier’, ‘group1’, ‘group2’,
  ‘mean(group1)’, ‘mean(group2)’, ‘std(group1)’, ‘std(group2)’, ‘Log2FC’, ‘FC’,
  ‘rejected’, ‘T-statistics’, ‘p-value’, ‘correction’, ‘-log10 p-value’, and ‘method’.

Example:

```default
result = run_ttest(df,
                   condition1='group1',
                   condition2='group2',
                   alpha = 0.05,
                   drop_cols=['sample'],
                   subject='subject',
                   group='group',
                   paired=False,
                   correction='fdr_bh',
                   permutations=50
        )
```

### run_two_way_anova(df, drop_cols=['sample'], subject='subject', group=['group', 'secondary_group'])

Run a 2-way ANOVA when data[‘secondary_group’] is not empty

* **Parameters:**
  * **df** (*pd.DataFrame*) – processed pandas DataFrame with samples as rows,
    and proteins and groups as columns.
  * **drop_cols** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – column names to drop from DataFrame
  * **subject** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column name containing subject identifiers.
  * **group** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – column names corresponding to independent variable groups
* **Returns:**
  Two DataFrames, anova results and residuals.

Example:

```default
result = run_two_way_anova(data,
                           drop_cols=['sample'],
                           subject='subject',
                           group=['group', 'secondary_group']
        )
```

## Submodules

## acore.differential_regulation.format_test_table module

Formatting related functions for test tables.

## acore.differential_regulation.tests module

All the tests for differential regulation. Functions used in the user facing
function starting with run_.

### calc_means_between_groups(df: DataFrame, boolean_array: Series, event_names: [tuple](https://docs.python.org/3/library/stdtypes.html#tuple)[[str](https://docs.python.org/3/library/stdtypes.html#str), [str](https://docs.python.org/3/library/stdtypes.html#str)] = ('1', '0')) → DataFrame

Mean comparison between groups

### calc_ttest(df: DataFrame, boolean_array: Series, variables: [list](https://docs.python.org/3/library/stdtypes.html#list)[[str](https://docs.python.org/3/library/stdtypes.html#str)]) → DataFrame

Calculate t-test for each variable in variables between two groups defined
by boolean array.

### calculate_ttest(df, condition1, condition2, paired=False, is_logged=True, non_par=False, tail='two-sided', correction='auto', r=0.707)

Calculates the t-test for the means of independent samples belonging to two different
groups using [scipy.stats.ttest_ind](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html.).

* **Parameters:**
  * **df** – pandas dataframe with groups and subjects as rows and protein identifier
    as column.
  * **condition1** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – identifier of first group.
  * **condition2** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – ientifier of second group.
  * **is_logged** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – data is logged transformed
  * **non_par** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – if True, normality and variance equality assumptions are checked
    and non-parametric test Mann Whitney U test if not passed
* **Returns:**
  Tuple with t-statistics, two-tailed p-value, mean of first group,
  mean of second group and logfc.

Example:

```default
result = calculate_ttest(df, 'group1', 'group2')
```

### calculate_thsd(df, column, group='group', alpha=0.05, is_logged=True)

Pairwise Tukey-HSD posthoc test using [pingouin.pairwise_tukey](https://pingouin-stats.org/build/html/generated/pingouin.pairwise_tukey.html).

* **Parameters:**
  * **df** – pandas dataframe with group and protein identifier as columns
  * **column** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column containing the protein identifier
  * **group** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column label containing the between factor
  * **alpha** ([*float*](https://docs.python.org/3/library/functions.html#float)) – significance level
* **Returns:**
  Pandas dataframe.

Example:

```default
result = calculate_thsd(df, column='HBG2~P69892', group='group', alpha=0.05)
```

### calculate_pairwise_ttest(df, column, subject='subject', group='group', correction='none', is_logged=True)

Performs pairwise t-test using pingouin, as a posthoc test,
and calculates fold-changes using [pingouin.pairwise_ttests](https://pingouin-stats.org/build/html/generated/pingouin.pairwise_ttests.html.).

* **Parameters:**
  * **df** – pandas dataframe with subject and group as rows and protein identifier as column.
  * **column** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column label containing the dependant variable
  * **subject** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column label containing subject identifiers
  * **group** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column label containing the between factor
  * **correction** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – method used for testing and adjustment of p-values.
* **Returns:**
  Pandas dataframe with means, standard deviations, test-statistics,
  degrees of freedom and effect size columns.

Example:

```default
result = calculate_pairwise_ttest(df,
                                  'protein a',
                                  subject='subject',
                                  group='group',
                                  correction='none'
        )
```

### complement_posthoc(posthoc, identifier, is_logged)

Calculates fold-changes after posthoc test.

* **Parameters:**
  * **posthoc** – pandas dataframe from posthoc test. Should have at least columns
    ‘mean(group1)’ and ‘mean(group2)’.
  * **identifier** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – feature identifier.
* **Returns:**
  Pandas dataframe with additional columns ‘identifier’, ‘log2FC’ and ‘FC’.

### calculate_anova(df, column, group='group')

Calculates one-way ANOVA using pingouin.

* **Parameters:**
  * **df** – pandas dataframe with group as rows and protein identifier as column
  * **column** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – name of the column in df to run ANOVA on
  * **group** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column with group identifiers
* **Returns:**
  Tuple with t-statistics and p-value.

### calculate_ancova(data, column, group='group', covariates=[])

Calculates one-way ANCOVA using pingouin.

* **Parameters:**
  * **df** – pandas dataframe with group as rows and protein identifier as column
  * **column** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – name of the column in df to run ANOVA on
  * **group** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column with group identifiers
  * **covariates** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – list of covariates (columns in df)
* **Returns:**
  Tuple with column, F-statistics and p-value.

### calculate_repeated_measures_anova(df, column, subject='subject', within='group')

One-way and two-way repeated measures ANOVA using pingouin stats.

* **Parameters:**
  * **df** – pandas dataframe with samples as rows and protein identifier as column.
    Data must be in long-format for two-way repeated measures.
  * **column** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column label containing the dependant variable
  * **subject** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column label containing subject identifiers
  * **within** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column label containing the within factor
* **Returns:**
  Tuple with protein identifier, t-statistics and p-value.

Example:

```default
result = calculate_repeated_measures_anova(df,
                                          'protein a',
                                          subject='subject',
                                          within='group'
        )
```

### calculate_mixed_anova(df, column, subject='subject', within='group', between='group2')

One-way and two-way repeated measures ANOVA using pingouin stats.

* **Parameters:**
  * **df** – pandas dataframe with samples as rows and protein identifier as column.
    Data must be in long-format for two-way repeated measures.
  * **column** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column label containing the dependant variable
  * **subject** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column label containing subject identifiers
  * **within** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column label containing the within factor
  * **within** – column label containing the between factor
* **Returns:**
  Tuple with protein identifier, t-statistics and p-value.

Example:

```default
result = calculate_mixed_anova(df,
                               'protein a',
                               subject='subject',
                               within='group',
                               between='group2'
        )
```

### pairwise_ttest_with_covariates(df, column, group, covariates, is_logged)

Pairwise t-test with covariates using statsmodels.

### format_anova_table(df, aov_results, pairwise_results, pairwise_cols, group, permutations, alpha, correction)

Performs p-value correction (permutation-based and FDR) and converts pandas dataframe
into final format.

* **Parameters:**
  * **df** – pandas dataframe with samples as rows and protein identifiers as columns
    (with additional columns ‘group’, ‘sample’ and ‘subject’).
  * **aov_results** ([*list*](https://docs.python.org/3/library/stdtypes.html#list) *[*[*tuple*](https://docs.python.org/3/library/stdtypes.html#tuple) *]*) – list of tuples with anova results (one tuple per feature).
  * **pairwise_results** ([*list*](https://docs.python.org/3/library/stdtypes.html#list) *[**dataframes* *]*) – list of pandas dataframes with
    posthoc tests results
  * **group** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column with group identifiers
  * **alpha** ([*float*](https://docs.python.org/3/library/functions.html#float)) – error rate for multiple hypothesis correction
  * **permutations** ([*int*](https://docs.python.org/3/library/functions.html#int)) – number of permutations used to estimate false discovery rates
* **Returns:**
  Pandas dataframe

### calculate_pvalue_from_tstats(tstat, dfn)

Calculate two-tailed p-values from T- or F-statistics.

tstat: T/F distribution
dfn: degrees of freedrom *n* (values) per protein (keys),

> i.e. number of obervations - number of groups (dict)

### eta_squared(aov)

Calculates the effect size using Eta-squared.

* **Parameters:**
  **aov** – pandas dataframe with anova results from statsmodels.
* **Returns:**
  Pandas dataframe with additional Eta-squared column.

### omega_squared(aov)

Calculates the effect size using Omega-squared.

* **Parameters:**
  **aov** – pandas dataframe with anova results from statsmodels.
* **Returns:**
  Pandas dataframe with additional Omega-squared column.
