# acore.correlation_analysis package

### *class* CorrelationCoefficient(coefficient, pvalue)

Bases: [`NamedTuple`](https://docs.python.org/3/library/typing.html#typing.NamedTuple)

#### coefficient *: [float](https://docs.python.org/3/library/functions.html#float)*

Alias for field number 0

#### pvalue *: [float](https://docs.python.org/3/library/functions.html#float)*

Alias for field number 1

#### count(value,)

Return number of occurrences of value.

#### index(value, start=0, stop=9223372036854775807,)

Return first index of value.

Raises ValueError if the value is not present.

### corr_lower_triangle(df: DataFrame, \*\*kwargs) → DataFrame

Compute the correlation matrix, returning only unique values (lower triangle).
Passes kwargs to [pandas.DataFrame.corr](pandas.DataFrame.corr) method.

### calculate_correlations(x, y, method='pearson')

Calculates a Spearman (nonparametric)
or a Pearson (parametric) correlation coefficient
and p-value to test for non-correlation.

* **Parameters:**
  * **x** ([*numpy.ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – array 1
  * **y** ([*numpy.ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – array 2
  * **method** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – chooses which kind of correlation method to run
* **Returns:**
  Tuple with two floats, correlation coefficient and two-tailed p-value.

Example:

```default
result = calculate_correlations(x, y, method='pearson')
```

### run_correlation(df, alpha=0.05, subject=None, group='group', method='pearson', correction='fdr_bh', numeric_only=True, dropna=True)

This function calculates pairwise correlations for columns in dataframe,
and returns it in the shape of a edge list with ‘weight’ as correlation score,
and the ajusted p-values.

* **Parameters:**
  * **df** – pandas dataframe with samples as rows and features as columns.
  * **subject** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – name of column containing subject identifiers.
  * **group** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – name of column containing group identifiers.
  * **method** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – method to use for correlation calculation (‘pearson’, ‘spearman’).
  * **alpha** ([*float*](https://docs.python.org/3/library/functions.html#float)) – error rate. Values velow alpha are considered significant.
  * **correction** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – type of correction see apply_pvalue_correction for methods
  * **numeric_only** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – if True, only numeric columns are considered
    for correlation calculation.
  * **dropna** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – if True, columns with NaN values are dropped before
    correlation calculation.
* **Returns:**
  Pandas dataframe with columns: ‘node1’, ‘node2’, ‘weight’,
  ‘padj’ and ‘rejected’.

Example:

```default
result = run_correlation(df, alpha=0.05, subject='subject', group='group',
            method='pearson', correction='fdr_bh')
```

### run_multi_correlation(df_dict, alpha=0.05, subject='subject', on=['subject', 'biological_sample'], group='group', method='pearson', correction='fdr_bh')

This function merges all input dataframes and calculates pairwise correlations for all columns.

* **Parameters:**
  * **df_dict** ([*dict*](https://docs.python.org/3/library/stdtypes.html#dict)) – dictionary of pandas dataframes with samples as rows and features as columns.
  * **subject** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – name of the column containing subject identifiers.
  * **group** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – name of the column containing group identifiers.
  * **on** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – column names to join dataframes on (must be found in all dataframes).
  * **method** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – method to use for correlation calculation (‘pearson’, ‘spearman’).
  * **alpha** ([*float*](https://docs.python.org/3/library/functions.html#float)) – error rate. Values velow alpha are considered significant.
  * **correction** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – type of correction see apply_pvalue_correction for methods
* **Returns:**
  Pandas dataframe with columns: ‘node1’, ‘node2’, ‘weight’, ‘padj’ and ‘rejected’.

Example:

```default
result = run_multi_correlation(df_dict, alpha=0.05, subject='subject',
    on=['subject', 'biological_sample'],
    group='group', method='pearson', correction='fdr_bh')
```

### calculate_rm_correlation(df, x, y, subject)

Computes correlation and p-values between two columns a and b in df with
repeated measures (rm).

* **Parameters:**
  * **df** – pandas dataframe with subjects as rows and two features and columns.
  * **x** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – feature a name.
  * **y** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – feature b name.
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
  * **subject** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – name of column containing subject identifiers.
  * **alpha** ([*float*](https://docs.python.org/3/library/functions.html#float)) – error rate. Values velow alpha are considered significant.
  * **correction** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – type of correction type see apply_pvalue_correction for methods
* **Returns:**
  Pandas dataframe with columns: ‘node1’, ‘node2’,
  ‘weight’, ‘pvalue’, ‘dof’, ‘padj’ and ‘rejected’.

Example:

```default
result = run_rm_correlation(df, alpha=0.05, subject='subject', correction='fdr_bh')
```

### calculate_pvalue_correlation_old(r: DataFrame, n_obs: [int](https://docs.python.org/3/library/functions.html#int)) → DataFrame

Calculate p-values for Pearson correlation using a all values from the
correlation matrix.

* **Parameters:**
  * **r** (*pd.DataFrame*) – Correlation matrix.
  * **n_obs** ([*int*](https://docs.python.org/3/library/functions.html#int)) – Number of observations used to calculate the correlation matrix (assumes no
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
  * **n_obs** ([*int*](https://docs.python.org/3/library/functions.html#int)) – Number of observations used to calculate the correlation matrix (assumes no
    missing values).
* **Returns:**
  p-value matrix assuming fixed number of observations.
* **Return type:**
  pd.DataFrame

### run_efficient_correlation(data, method='pearson')

Calculates pairwise correlations and returns lower triangle of the matrix with
correlation values and p-values. For pearson correlation, p-values are calculated
assuming a fixed number of observations.

* **Parameters:**
  * **data** – pandas dataframe with samples as index and features as columns (numeric data only).
  * **method** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – method to use for correlation calculation (‘pearson’, ‘spearman’).
* **Returns:**
  Two numpy arrays: correlation and p-values.

Example:

```default
result = run_efficient_correlation(data, method='pearson')
```
