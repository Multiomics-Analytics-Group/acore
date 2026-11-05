# acore.differential_regulation package

Differential regulation module.

### run_ancova(df: DataFrame, covariates: [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)], alpha: [float](https://docs.python.org/3/builtins/functions.html#float) = 0.05, drop_cols: [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None, subject: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'subject', group: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'group', permutations: [int](https://docs.python.org/3/builtins/functions.html#int) = 0, correction: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'fdr_bh', is_logged: [bool](https://docs.python.org/3/builtins/functions.html#bool) = True, non_par: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[[AncovaSchema](acore.types.md#acore.types.differential_analysis.AncovaSchema)]

Compares groups for each protein with a one-way ANCOVA, adjusting for covariates.

Checks what type of data is the input (paired, unpaired or repeated measurements)
and performs posthoc tests for multiclass data.
Multiple hypothesis correction uses permutation-based
if permutations>0 and Benjamini/Hochberg if permutations=0.

* **Parameters:**
  * **df** (*pd.DataFrame*) – Pandas DataFrame with samples as rows and protein identifiers and
    covariates as columns (with additional columns ‘group’, ‘sample’ and ‘subject’).
  * **covariates** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – list of covariates to include in the model (column in df).
    Always a list of column names, also for one covariate:
    `covariates=['age']`, not `'age'`. The parameter is
    called `covariates`, not `covariate`. A covariate is a
    column of the table, so do not also list it in `drop_cols`.
  * **alpha** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – error rate for multiple hypothesis correction
  * **drop_cols** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – column labels to be dropped from the DataFrame. Pass `None` or `[]`
    to drop no columns.
    List every column that is
    not a measurement, except `group` and the covariates: a text
    column such as ‘sample’, ‘subject’ or another label column
    left in the table fails with `endog has evaluated to an array
    with multiple columns`, and a numeric one that is not a
    covariate would be tested like a protein. For example
    `drop_cols=['sample', 'subject']`.
  * **subject** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – not used by this function. A subject column is removed through
    `drop_cols` like any other identifier.
  * **group** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – name of the column that holds the group labels, as written in the
    table header. The default is ‘group’; pass another name when the
    labels are in a differently named column, for example ‘strain’ or
    ‘arm’.
  * **permutations** ([*int*](https://docs.python.org/3/builtins/functions.html#int)) – number of permutations used to estimate false discovery rates.
  * **correction** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – method of pvalue correction see apply_pvalue_correction for methods,
    use methods available in acore.multiple_testing
  * **is_logged** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool)) – whether data is log-transformed
  * **non_par** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool)) – if True, normality and variance equality assumptions are checked
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
                    permutations=0
        )

# columns: sample, subject, group, arm, age and the proteins; compare the arms
result = run_ancova(df,
                    covariates=['age'],
                    group='arm',
                    drop_cols=['sample', 'subject', 'group']
        )
```

### run_anova(df: DataFrame, alpha: [float](https://docs.python.org/3/builtins/functions.html#float) = 0.05, drop_cols: [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None, subject: [str](https://docs.python.org/3/builtins/stdtypes.html#str) | [None](https://docs.python.org/3/builtins/constants.html#None) = None, group: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'group', permutations: [int](https://docs.python.org/3/builtins/functions.html#int) = 0, correction: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'fdr_bh', is_logged: [bool](https://docs.python.org/3/builtins/functions.html#bool) = True, non_par: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[[AnovaSchema](acore.types.md#acore.types.differential_analysis.AnovaSchema)] | [DataFrame](https://pandera.readthedocs.io/en/stable/reference/generated/pandera.typing.DataFrame.html#pandera.typing.DataFrame)[[AnovaSchemaMultiGroup](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup)]

Compares two or more groups for each protein (t-test for two, one-way ANOVA for more).

Checks what type of data is the input (paired, unpaired or repeated measurements) and
performs posthoc tests for multiclass data (i.e., when there are more than two groups,
posthoc tests such as pairwise t-tests or Tukey’s HSD are used to determine which specific
groups differ after finding a significant overall effect).
Multiple hypothesis correction uses permutation-based
if permutations>0 and Benjamini/Hochberg if permutations=0.

* **Parameters:**
  * **df** (*pd.DataFrame*) – pandas dataframe with samples as rows and protein identifiers as columns
    (with additional columns ‘group’, ‘sample’ and ‘subject’).
  * **alpha** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – error rate for multiple hypothesis correction
  * **drop_cols** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – column labels to be dropped from the dataframe. Pass `None` or `[]`
    to drop no columns.
    List every column that is
    not a measurement, except `group`: a text column such as
    ‘sample’, ‘subject’ or ‘arm’ left in the table fails with
    `ufunc 'divide' not supported`, and a numeric one such as
    ‘age’ would be tested like a protein. There is one column per
    protein, so no ‘protein’ or ‘protein_id’ column exists to drop.
    For example `drop_cols=['sample', 'subject']`.
  * **subject** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column with subject identifiers
  * **group** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – name of the column that holds the group labels, as written in the
    table header. The default is ‘group’; pass another name when the
    labels are in a differently named column, for example ‘strain’ or
    ‘arm’.
  * **permutations** ([*int*](https://docs.python.org/3/builtins/functions.html#int)) – number of permutations used to estimate false discovery rates.
  * **correction** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – method of pvalue correction see apply_pvalue_correction for methods,
    use methods available in acore.multiple_testing
  * **is_logged** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool)) – whether data is log-transformed
  * **non_par** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool)) – if True, normality and variance equality assumptions are checked
    and non-parametric test Mann Whitney U test if not passed
* **Returns:**
  DataFrame adhering to AnovaSchema or AnovaSchemaMultiGroup. For more than
  two groups, `posthoc comparison` labels each pair as `group1~~group2`.
* **Return type:**
  DataFrame[[AnovaSchema](acore.types.md#acore.types.differential_analysis.AnovaSchema)] | DataFrame[[AnovaSchemaMultiGroup](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup)]

Example:

```default
result = run_anova(df,
                   alpha=0.05,
                   drop_cols=["sample",'subject'],
                   subject='subject',
                   group='group',
                   permutations=0
        )

# columns: sample, subject, group, arm, age and the proteins; compare the arms
result = run_anova(df,
                   group='arm',
                   drop_cols=['sample', 'subject', 'group', 'age']
        )
```

### run_diff_analysis(df: DataFrame, boolean_array: Series, event_names: [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[str](https://docs.python.org/3/builtins/stdtypes.html#str), [str](https://docs.python.org/3/builtins/stdtypes.html#str)] = ('1', '0'), ttest_vars=('alternative', 'p-val', 'cohen-d')) → DataFrame

Differential analysis procedure between two groups. Calculaes
mean per group and t-test for each variable in vars between two groups.

### run_mixed_anova(df: DataFrame, alpha: [float](https://docs.python.org/3/builtins/functions.html#float) = 0.05, drop_cols: [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None, subject: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'subject', within: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'group', between: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'group2', correction: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'fdr_bh') → DataFrame

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
  * **alpha** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – error rate for multiple hypothesis correction
  * **drop_cols** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – column labels to be dropped from the DataFrame. Pass `None` or `[]`
    to drop no columns.
  * **subject** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column with subject identifiers
  * **within** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column with within factor identifiers
  * **between** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column with between factor identifiers
  * **correction** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – method of pvalue correction see apply_pvalue_correction for methods,
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

### run_repeated_measurements_anova(df: DataFrame, alpha: [float](https://docs.python.org/3/builtins/functions.html#float) = 0.05, drop_cols: [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None, subject: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'subject', within: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'group', permutations: [int](https://docs.python.org/3/builtins/functions.html#int) = 50, correction: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'fdr_bh', is_logged: [bool](https://docs.python.org/3/builtins/functions.html#bool) = True) → DataFrame

Performs repeated measurements anova and pairwise posthoc tests for each protein in dataframe.

* **Parameters:**
  * **df** (*pd.DataFrame*) – Pandas DataFrame with samples as rows and protein identifiers as columns
    (with additional columns ‘group’, ‘sample’ and ‘subject’).
  * **alpha** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – error rate for multiple hypothesis correction
  * **drop_cols** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – column labels to be dropped from the DataFrame. Pass `None` or `[]`
    to drop no columns. The `subject` column must not be included here as
    it is required for the analysis.
    List every other column
    that is not a measurement, except `within`: a text column
    such as ‘sample’ or ‘arm’ left in the table fails with
    `DV must be numeric`, and a numeric one such as ‘age’ would
    be tested like a protein. For example
    `drop_cols=['sample', 'arm', 'age']`.
  * **subject** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column with subject identifiers, as written in the table header
    (for example ‘subject’, not ‘subject_id’ unless the table has it).
    Every subject has several rows, one per level of `within`.
  * **within** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – name of the column that holds the repeated conditions, such as
    the visits or time points T0, T1, T2. It must be a column of the
    table.
  * **permutations** ([*int*](https://docs.python.org/3/builtins/functions.html#int)) – number of permutations used to estimate false discovery
    rates. Pass `0` to use Benjamini/Hochberg correction
    instead. Note: the default of 50 fails in this
    function with `KeyError: ['F', 'p-unc'] not in index`.
  * **correction** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – method of pvalue correction see apply_pvalue_correction for methods,
    use methods available in acore.multiple_testing
  * **is_logged** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool)) – whether data is log-transformed
* **Returns:**
  Pandas DataFrame

Example:

```default
# columns: sample, subject, group (T0, T1, T2), arm, age and the proteins
result = run_repeated_measurements_anova(df,
                                         alpha=0.05,
                                         drop_cols=['sample', 'arm', 'age'],
                                         subject='subject',
                                         within='group',
                                         permutations=0
        )
```

### run_ttest(df: DataFrame, condition1: [str](https://docs.python.org/3/builtins/stdtypes.html#str), condition2: [str](https://docs.python.org/3/builtins/stdtypes.html#str), alpha: [float](https://docs.python.org/3/builtins/functions.html#float) = 0.05, drop_cols: [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None, subject: [str](https://docs.python.org/3/builtins/stdtypes.html#str) | [None](https://docs.python.org/3/builtins/constants.html#None) = 'subject', group: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'group', paired: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False, correction: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'fdr_bh', permutations: [int](https://docs.python.org/3/builtins/functions.html#int) = 0, is_logged: [bool](https://docs.python.org/3/builtins/functions.html#bool) = True, non_par: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → DataFrame

Runs a t-test between exactly two groups; needs both condition1 and condition2.

Runs t-test (paired/unpaired) for each protein in dataset and performs
permutation-based (if permutations>0) or Benjamini/Hochberg (if permutations=0)
multiple hypothesis correction.

This compares exactly two groups of one factor, so every call needs both
`condition1` and `condition2`, two values from the `group` column. To compare
each of several groups with a reference, call it once per group, for example
`run_ttest(df, condition1='WT', condition2='zwf', group='strain', subject=None,
drop_cols=['sample'])` for a table without a `subject` column, and again with
`condition2='sucA'`. Use
`run_anova` to compare three or more groups at once, and `run_two_way_anova` for
a design with two factors.

* **Parameters:**
  * **df** (*pd.DataFrame*) – Pandas DataFrame with samples as rows and protein identifiers as columns
    (with additional columns ‘group’, ‘sample’ and ‘subject’).
  * **condition1** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – required; one value of the `group` column, copied exactly as
    it is written in the table (for example ‘WT’, not ‘wild type’),
    for the first of the two groups. It must be a value that
    occurs in that column.
  * **condition2** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – required; another value of the `group` column, copied exactly
    as written, for the second of the two groups. To compare several
    groups with one reference, call this function once per pair.
  * **alpha** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – error rate for multiple hypothesis correction
  * **drop_cols** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – column labels to be dropped from the DataFrame. Pass `None` or `[]`
    to drop no columns. The `subject` column should not be included here
    as it is handled separately based on the `paired` parameter.
    List every other column
    that is not a measurement, except `group`: a text column such
    as ‘sample’ or ‘treatment’ left in the table fails with
    `ufunc 'divide' not supported`, and a numeric one such as
    ‘age’ would be tested like a protein. There is one column per
    protein, so no ‘protein’ or ‘protein_id’ column exists to drop.
    For example `drop_cols=['sample']`.
  * **subject** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column with subject identifiers. The default is ‘subject’; the
    column must exist in the table then, and it is dropped for an
    unpaired test, so never also list it in `drop_cols`. Pass
    `subject=None` when the table has no such column, as in a
    table with only ‘sample’, ‘strain’ and the measurements. Both
    mistakes fail with `['subject'] not found in axis`. A paired
    test needs it.
  * **group** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – name of the column that holds the group labels (independent
    variable), as written in the table header. The default is ‘group’;
    pass another name when the labels are in a differently named
    column, for example ‘strain’ or ‘genotype’. `condition1` and
    `condition2` are values in this
    column.
  * **paired** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool)) – paired or unpaired samples
  * **correction** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – method of pvalue correction see apply_pvalue_correction for methods
  * **permutations** ([*int*](https://docs.python.org/3/builtins/functions.html#int)) – number of permutations used to estimate false discovery rates.
  * **is_logged** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool)) – data is log-transformed
  * **non_par** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool)) – if True, normality and variance equality assumptions are checked
    and non-parametric test Mann Whitney U test if not passed
* **Returns:**
  Pandas DataFrame with columns ‘identifier’, ‘group1’, ‘group2’,
  ‘mean(group1)’, ‘mean(group2)’, ‘std(group1)’, ‘std(group2)’, ‘Log2FC’, ‘FC’,
  ‘rejected’, ‘T-statistics’, ‘p-value’, ‘correction’, ‘-log10 p-value’, and ‘method’.

Example:

```default
# columns: sample, subject, group ('sepsis' or 'healthy') and the proteins
result = run_ttest(df,
                   condition1='sepsis',
                   condition2='healthy',
                   alpha = 0.05,
                   drop_cols=['sample'],
                   subject='subject',
                   group='group',
                   paired=False,
                   correction='fdr_bh',
                   permutations=0
        )

# columns: sample, strain ('WT', 'zwf', ...) and the metabolites; no subject
result = run_ttest(df,
                   condition1='WT',
                   condition2='zwf',
                   group='strain',
                   subject=None,
                   drop_cols=['sample']
        )
```

### run_two_way_anova(df: DataFrame, drop_cols: [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None, subject: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'subject', group: [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[DataFrame, [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[str](https://docs.python.org/3/builtins/stdtypes.html#str), Series]]

Run a 2-way ANOVA when data[‘secondary_group’] is not empty

* **Parameters:**
  * **df** (*pd.DataFrame*) – processed pandas DataFrame with samples as rows,
    and proteins and groups as columns.
  * **drop_cols** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – column names to drop from DataFrame. Pass `None` or `[]`
    to drop no columns.
    List every column that
    is not a measurement, except `subject` and the two `group`
    columns: a text column such as ‘sample’ left in the table
    fails, and a numeric one such as ‘age’ would be tested like a
    protein. For example `drop_cols=['sample']`.
  * **subject** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column name containing subject identifiers. It must be a column
    of the table, because the function uses it as an index.
  * **group** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – names of the two columns that hold the labels of the two factors
    (independent variables), for example [‘genotype’, ‘treatment’].
    Defaults to [‘group’, ‘secondary_group’] if None, so name them when
    the table calls these columns anything else. Both must be columns
    of the table.
* **Returns:**
  Two DataFrames, anova results and residuals.

Example:

```default
# columns: sample, subject, genotype, treatment and the metabolites
result = run_two_way_anova(data,
                           drop_cols=['sample'],
                           subject='subject',
                           group=['genotype', 'treatment']
        )
```

## Submodules

## acore.differential_regulation.format_test_table module

Formatting related functions for test tables.

## acore.differential_regulation.tests module

All the tests for differential regulation. Functions used in the user facing
function starting with run_.

### calc_means_between_groups(df: DataFrame, boolean_array: Series, event_names: [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[str](https://docs.python.org/3/builtins/stdtypes.html#str), [str](https://docs.python.org/3/builtins/stdtypes.html#str)] = ('1', '0')) → DataFrame

Mean comparison between groups

### calc_ttest(df: DataFrame, boolean_array: Series, variables: [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)]) → DataFrame

Calculate t-test for each variable in variables between two groups defined
by boolean array.

### calculate_ttest(df: DataFrame, condition1: [str](https://docs.python.org/3/builtins/stdtypes.html#str), condition2: [str](https://docs.python.org/3/builtins/stdtypes.html#str), paired: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False, is_logged: [bool](https://docs.python.org/3/builtins/functions.html#bool) = True, non_par: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False, tail: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'two-sided', correction: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'auto', r: [float](https://docs.python.org/3/builtins/functions.html#float) = 0.707) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[float](https://docs.python.org/3/builtins/functions.html#float), [float](https://docs.python.org/3/builtins/functions.html#float), [float](https://docs.python.org/3/builtins/functions.html#float), [float](https://docs.python.org/3/builtins/functions.html#float), [float](https://docs.python.org/3/builtins/functions.html#float), [float](https://docs.python.org/3/builtins/functions.html#float), [float](https://docs.python.org/3/builtins/functions.html#float), [str](https://docs.python.org/3/builtins/stdtypes.html#str)]

Calculates the t-test for the means of independent samples belonging to two different
groups using [scipy.stats.ttest_ind](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html.).

* **Parameters:**
  * **df** – pandas dataframe with groups and subjects as rows and protein identifier
    as column.
  * **condition1** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – identifier of first group.
  * **condition2** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – ientifier of second group.
  * **is_logged** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool)) – data is logged transformed
  * **non_par** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool)) – if True, normality and variance equality assumptions are checked
    and non-parametric test Mann Whitney U test if not passed
* **Returns:**
  Tuple with t-statistics, two-tailed p-value, mean of first group,
  mean of second group and logfc.

Example:

```default
result = calculate_ttest(df, 'group1', 'group2')
```

### calculate_thsd(df, column, group='group', alpha=0.05, is_logged=True)

Pairwise Tukey-HSD posthoc test using [pingouin.pairwise_tukey](https://pingouin-stats.org/generated/pingouin.pairwise_tukey.html).

* **Parameters:**
  * **df** – pandas dataframe with group and protein identifier as columns
  * **column** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column containing the protein identifier
  * **group** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column label containing the between factor
  * **alpha** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – significance level
* **Returns:**
  Pandas dataframe.

Example:

```default
result = calculate_thsd(df, column='HBG2~P69892', group='group', alpha=0.05)
```

### calculate_pairwise_ttest(df: DataFrame, column: [str](https://docs.python.org/3/builtins/stdtypes.html#str), subject: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'subject', group: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'group', correction: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'none', is_logged: [bool](https://docs.python.org/3/builtins/functions.html#bool) = True) → DataFrame

Performs pairwise t-test using pingouin, as a posthoc test,
and calculates fold-changes using [pingouin.pairwise_tests](https://pingouin-stats.org/generated/pingouin.pairwise_tests.html.).

* **Parameters:**
  * **df** – pandas dataframe with subject and group as rows and protein identifier as column.
  * **column** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column label containing the dependant variable
  * **subject** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column label containing subject identifiers
  * **group** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column label containing the between factor
  * **correction** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – method used for testing and adjustment of p-values.
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
  * **identifier** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – feature identifier.
* **Returns:**
  Pandas dataframe with additional columns ‘identifier’, ‘log2FC’ and ‘FC’.

### calculate_anova(df: DataFrame, column: [str](https://docs.python.org/3/builtins/stdtypes.html#str), group: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'group')

Calculates one-way ANOVA using pingouin.

* **Parameters:**
  * **df** – pandas dataframe with group as rows and protein identifier as column
  * **column** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – name of the column in df to run ANOVA on
  * **group** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column with group identifiers
* **Returns:**
  Tuple with t-statistics and p-value.

### calculate_ancova(data: DataFrame, column: [str](https://docs.python.org/3/builtins/stdtypes.html#str), group: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'group', covariates: [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None)

Calculates one-way ANCOVA using pingouin.

* **Parameters:**
  * **data** – pandas dataframe with group as rows and protein identifier as column.
  * **column** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – name of the one measurement column in data to run ANOVA on, for
    example ‘CRP’. It must be a different column from `group` and
    from every covariate: `column='arm', group='arm'` fails with
    `categorical data cannot be >1-dimensional`.
  * **group** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column with group identifiers, as written in the table header
  * **covariates** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – list of covariates (columns in data), always a list of
    column names, for example `['age']`
* **Returns:**
  Tuple with column, F-statistics and p-value.

Example:

```default
# columns: arm, age, CRP, IL6 and more measurements
result = calculate_ancova(data, column='CRP', group='arm', covariates=['age'])
```

### calculate_repeated_measures_anova(df: DataFrame, column: [str](https://docs.python.org/3/builtins/stdtypes.html#str), subject: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'subject', within: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'group')

One-way and two-way repeated measures ANOVA using pingouin stats.

* **Parameters:**
  * **df** – pandas dataframe with samples as rows and protein identifier as column.
    Data must be in long-format for two-way repeated measures.
  * **column** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column label containing the dependant variable
  * **subject** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column label containing subject identifiers
  * **within** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column label containing the within factor
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

### calculate_mixed_anova(df: DataFrame, column: [str](https://docs.python.org/3/builtins/stdtypes.html#str), subject: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'subject', within: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'group', between: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'group2')

One-way and two-way repeated measures ANOVA using pingouin stats.

* **Parameters:**
  * **df** – pandas dataframe with samples as rows and protein identifier as column.
    Data must be in long-format for two-way repeated measures.
  * **column** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column label containing the dependant variable
  * **subject** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column label containing subject identifiers
  * **within** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column label containing the within factor
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

### pairwise_ttest_with_covariates(df: DataFrame, column: [str](https://docs.python.org/3/builtins/stdtypes.html#str), group: [str](https://docs.python.org/3/builtins/stdtypes.html#str), covariates: [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)], is_logged: [bool](https://docs.python.org/3/builtins/functions.html#bool))

Pairwise t-test with covariates using statsmodels.

### format_anova_table(df: DataFrame, aov_results: [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)], pairwise_results, pairwise_cols, group: [str](https://docs.python.org/3/builtins/stdtypes.html#str), permutations: [int](https://docs.python.org/3/builtins/functions.html#int), alpha: [float](https://docs.python.org/3/builtins/functions.html#float), correction: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → DataFrame

Performs p-value correction (permutation-based and FDR) and converts pandas dataframe
into final format.

* **Parameters:**
  * **df** – pandas dataframe with samples as rows and protein identifiers as columns
    (with additional columns ‘group’, ‘sample’ and ‘subject’).
  * **aov_results** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list) *[*[*tuple*](https://docs.python.org/3/builtins/stdtypes.html#tuple) *]*) – list of tuples with anova results (one tuple per feature).
  * **pairwise_results** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list) *[**dataframes* *]*) – list of pandas dataframes with
    posthoc tests results
  * **group** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column with group identifiers
  * **alpha** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – error rate for multiple hypothesis correction
  * **permutations** ([*int*](https://docs.python.org/3/builtins/functions.html#int)) – number of permutations used to estimate false discovery rates
* **Returns:**
  Pandas dataframe

### calculate_pvalue_from_tstats(tstat, dfn)

Calculate two-tailed p-values from T- or F-statistics.

tstat: T/F distribution
dfn: degrees of freedrom *n* (values) per protein (keys),

> i.e. number of obervations - number of groups (dict)

### eta_squared(aov: DataFrame) → DataFrame

Calculates the effect size using Eta-squared.

* **Parameters:**
  **aov** – pandas dataframe with anova results from statsmodels.
* **Returns:**
  Pandas dataframe with additional Eta-squared column.

### omega_squared(aov: DataFrame) → DataFrame

Calculates the effect size using Omega-squared.

* **Parameters:**
  **aov** – pandas dataframe with anova results from statsmodels.
* **Returns:**
  Pandas dataframe with additional Omega-squared column.
