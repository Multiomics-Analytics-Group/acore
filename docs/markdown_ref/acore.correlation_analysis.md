# acore.correlation_analysis package

### *class* CorrelationCoefficient(coefficient, pvalue)

Bases: [`NamedTuple`](https://docs.python.org/3/library/typing.html#typing.NamedTuple)

#### coefficient *: [float](https://docs.python.org/3/builtins/functions.html#float)*

Alias for field number 0

#### pvalue *: [float](https://docs.python.org/3/builtins/functions.html#float)*

Alias for field number 1

#### count(value,)

Return number of occurrences of value.

#### index(value, start=0, stop=9223372036854775807,)

Return first index of value.

Raises ValueError if the value is not present.

### corr_lower_triangle(df: DataFrame, \*\*kwargs) → DataFrame

Compute the correlation matrix, returning only unique values (lower triangle).
Passes kwargs to [pandas.DataFrame.corr](pandas.DataFrame.corr) method.

### calculate_correlations(x: [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray), y: [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray), method: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'pearson') → [CorrelationCoefficient](#acore.correlation_analysis.CorrelationCoefficient)

Calculates a Spearman (nonparametric)
or a Pearson (parametric) correlation coefficient
and p-value to test for non-correlation.

* **Parameters:**
  * **x** ([*numpy.ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – array 1
  * **y** ([*numpy.ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – array 2
  * **method** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – chooses which kind of correlation method to run
* **Returns:**
  Tuple with two floats, correlation coefficient and two-tailed p-value.

Example:

```default
result = calculate_correlations(x, y, method='pearson')
```

### run_correlation(df: DataFrame, alpha: [float](https://docs.python.org/3/builtins/functions.html#float) = 0.05, subject: [str](https://docs.python.org/3/builtins/stdtypes.html#str) | [None](https://docs.python.org/3/builtins/constants.html#None) = None, group: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'group', method: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'pearson', correction: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'fdr_bh', numeric_only: [bool](https://docs.python.org/3/builtins/functions.html#bool) = True, dropna: [bool](https://docs.python.org/3/builtins/functions.html#bool) = True) → DataFrame

Correlates every pair of numeric columns and returns an edge list of the results.

This function calculates pairwise correlations for columns in dataframe,
and returns it in the shape of a edge list with ‘weight’ as correlation score,
and the ajusted p-values. It takes no `x` or `y`: to read off one pair, pick the
row of the result whose ‘node1’ and ‘node2’ are the two columns.

* **Parameters:**
  * **df** – pandas dataframe with samples as rows and features as columns.
  * **subject** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – name of column containing subject identifiers. Keep it `None`
    (the default) unless you want a repeated measures correlation, used when every
    subject appears in more than one row. With `None`, text columns such as
    ‘sample’ or ‘subject’ are skipped by `numeric_only`; drop a numeric
    non-measurement column from `df` before the call. Note: with fewer than 200
    columns and more than two subjects, the repeated measures route fails with
    `KeyError: None of [Index([1.0, 0.0, ...])] are in the [columns]` (the 1.0
    and 0.0 are the ‘rejected’ flags of the result).
  * **group** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – name of column containing group identifiers. It is only used
    together with `subject`, to check whether the samples are paired.
  * **method** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – method to use for correlation calculation: ‘pearson’ or ‘spearman’
    only. Any other name, such as ‘kendall’, fails with `cannot access local
    variable 'p'`.
  * **alpha** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – error rate. Values velow alpha are considered significant.
  * **correction** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – type of correction see apply_pvalue_correction for methods
  * **numeric_only** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool)) – if True, only numeric columns are considered
    for correlation calculation.
  * **dropna** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool)) – if True, columns with NaN values are dropped before
    correlation calculation.
* **Returns:**
  Pandas dataframe with columns: ‘node1’, ‘node2’, ‘weight’,
  ‘padj’ and ‘rejected’.

Example:

```default
result = run_correlation(df, alpha=0.05, method='pearson', correction='fdr_bh')
```

### run_multi_correlation(df_dict, alpha=0.05, subject='subject', on=None, group='group', method='pearson', correction='fdr_bh')

This function merges all input dataframes and calculates pairwise correlations for all columns.

* **Parameters:**
  * **df_dict** ([*dict*](https://docs.python.org/3/builtins/stdtypes.html#dict)) – dictionary of pandas dataframes with samples as rows and features as columns.
  * **subject** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – name of the column containing subject identifiers.
  * **group** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – name of the column containing group identifiers.
  * **on** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – column names to join dataframes on (must be found in all dataframes),
    defaults to [‘subject’, ‘biological_sample’] for None.
  * **method** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – method to use for correlation calculation (‘pearson’, ‘spearman’).
  * **alpha** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – error rate. Values velow alpha are considered significant.
  * **correction** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – type of correction see apply_pvalue_correction for methods
* **Returns:**
  Pandas dataframe with columns: ‘node1’, ‘node2’, ‘weight’, ‘padj’ and ‘rejected’.

Example:

```default
result = run_multi_correlation(df_dict, alpha=0.05, subject='subject',
    on=['subject', 'biological_sample'],
    group='group', method='pearson', correction='fdr_bh')
```

### calculate_rm_correlation(df: DataFrame, x: [str](https://docs.python.org/3/builtins/stdtypes.html#str), y: [str](https://docs.python.org/3/builtins/stdtypes.html#str), subject: [str](https://docs.python.org/3/builtins/stdtypes.html#str))

Computes correlation and p-values between two columns a and b in df with
repeated measures (rm).

* **Parameters:**
  * **df** – pandas dataframe with subjects as rows and two features and columns.
  * **x** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – feature a name.
  * **y** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – feature b name.
  * **subject** – column name containing the covariate variable.
* **Returns:**
  Tuple with values for: feature a, feature b, correlation, p-value and degrees of freedom.

Example:

```default
result = calculate_rm_correlation(df, x='feature a', y='feature b', subject='subject')
```

### run_rm_correlation(df, alpha=0.05, subject='subject', correction='fdr_bh')

Computes pairwise repeated measurements correlations for all columns in dataframe,
and returns results as an edge list with ‘weight’ as correlation score,
p-values, degrees of freedom and ajusted p-values.

* **Parameters:**
  * **df** – pandas dataframe with samples as rows and features as columns.
  * **subject** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – name of column containing subject identifiers.
  * **alpha** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – error rate. Values velow alpha are considered significant.
  * **correction** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – type of correction type see apply_pvalue_correction for methods
* **Returns:**
  Pandas dataframe with columns: ‘node1’, ‘node2’,
  ‘weight’, ‘pvalue’, ‘dof’, ‘padj’ and ‘rejected’.

Example:

```default
result = run_rm_correlation(df, alpha=0.05, subject='subject', correction='fdr_bh')
```

### calculate_pvalue_correlation_old(r: DataFrame, n_obs: [int](https://docs.python.org/3/builtins/functions.html#int)) → DataFrame

Calculate p-values for Pearson correlation using a all values from the
correlation matrix.

* **Parameters:**
  * **r** (*pd.DataFrame*) – Correlation matrix.
  * **n_obs** ([*int*](https://docs.python.org/3/builtins/functions.html#int)) – Number of observations used to calculate the correlation matrix (assumes no
    missing values).
* **Returns:**
  p-value matrix assuming fixed number of observations.
* **Return type:**
  pd.DataFrame

### calculate_pvalue_correlation(r, n_obs)

Calculate p-values for Pearson correlation using a all values from the
correlation matrix.

Tested against Pearson correlation using a all values from the correlation
matrix, see
[https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.pearsonr.html](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.pearsonr.html)

* **Parameters:**
  * **r** (*pd.DataFrame*) – Correlation matrix.
  * **n_obs** ([*int*](https://docs.python.org/3/builtins/functions.html#int)) – Number of observations used to calculate the correlation matrix (assumes no
    missing values).
* **Returns:**
  p-value matrix assuming fixed number of observations.
* **Return type:**
  pd.DataFrame

### run_efficient_correlation(data: DataFrame, method: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'pearson') → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray), [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)]

Calculates pairwise correlations and returns lower triangle of the matrix with
correlation values and p-values. For pearson correlation, p-values are calculated
assuming a fixed number of observations.

* **Parameters:**
  * **data** – pandas dataframe with samples as index and features as columns (numeric data only).
  * **method** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – method to use for correlation calculation (‘pearson’, ‘spearman’).
* **Returns:**
  Two numpy arrays: correlation and p-values.

Example:

```default
result = run_efficient_correlation(data, method='pearson')
```
