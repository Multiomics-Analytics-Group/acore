# acore.normalization package

Module normalization for analysis. This module provides functions to normalize data for analysis.

The higher-level convience functions are normalize_data and normalize_data_by_group.
Combat normalization was added using the inmoose package.

The actual normalization functions are in strategies.py.

### normalize_data(data: DataFrame, method: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'median', normalize: [str](https://docs.python.org/3/library/stdtypes.html#str) = None)

This function normalizes the data using the selected method. Normalizes only nummeric
data, but keeps the non-numeric columns in the output DataFrame.

* **Parameters:**
  * **data** – DataFrame with the data to be normalized (samples x features)
  * **method** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – normalization method to choose among: median (default),
    median_polish, median_zero, quantile, linear, zscore
  * **normalize** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – whether the normalization should be done by ‘features’ (columns)
    or ‘samples’ (rows) (default None)
* **Returns:**
  pandas.DataFrame.

Example:

```default
result = normalize_data(data, method='median_polish')
```

### normalize_data_per_group(data: DataFrame, group: [str](https://docs.python.org/3/library/stdtypes.html#str) | [int](https://docs.python.org/3/library/functions.html#int) | [list](https://docs.python.org/3/library/stdtypes.html#list)[[str](https://docs.python.org/3/library/stdtypes.html#str) | [int](https://docs.python.org/3/library/functions.html#int)], method: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'median', normalize: [str](https://docs.python.org/3/library/stdtypes.html#str) = None) → DataFrame

This function normalizes the data by group using the selected method

* **Parameters:**
  * **data** – DataFrame with the data to be normalized (samples x features)
  * **group_col** – Column containing the groups, passed to pandas.DataFrame.groupby
  * **method** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – normalization method to choose among: median_polish, median,
    quantile, linear
  * **normalize** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – whether the normalization should be done by ‘features’ (columns) or ‘samples’ (rows) (default None)
* **Returns:**
  pandas.DataFrame.

Example:

```default
result = normalize_data_per_group(data, group='group' method='median')
```

## Submodules

## acore.normalization.strategies module

Strategies for data normalization used by normalize_data.

### median_zero_normalization(data, normalize='samples')

This function normalizes each sample by using its median.

* **Parameters:**
  * **data**
  * **normalize** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – whether the normalization should be done by ‘features’ (columns)
    or ‘samples’ (rows)
* **Returns:**
  pandas.DataFrame.

Example:

```default
data = pd.DataFrame({'a': [2,5,4,3,3], 'b':[4,4,6,5,3], 'c':[4,14,8,8,9]})
result = median_zero_normalization(data, normalize='samples')
result
        a         b         c
    0 -1.333333  0.666667  0.666667
    1 -2.666667 -3.666667  6.333333
    2 -2.000000  0.000000  2.000000
    3 -2.333333 -0.333333  2.666667
    4 -2.000000 -2.000000  4.000000
```

### median_normalization(data, normalize='samples')

This function normalizes each sample by using its median.

* **Parameters:**
  * **data**
  * **normalize** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – whether the normalization should be done by ‘features’ (columns)
    or ‘samples’ (rows)
* **Returns:**
  pandas.DataFrame.

Example:

```default
data = pd.DataFrame({'a': [2,5,4,3,3], 'b':[4,4,6,5,3], 'c':[4,14,8,8,9]})
result = median_normalization(data, normalize='samples')
result
        a         b         c
    0 -1.333333  0.666667  0.666667
    1 -2.666667 -3.666667  6.333333
    2 -2.000000  0.000000  2.000000
    3 -2.333333 -0.333333  2.666667
    4 -2.000000 -2.000000  4.000000
```

### zscore_normalization(data, normalize='samples')

This function normalizes each sample by using its mean and standard deviation
(mean=0, std=1).

* **Parameters:**
  * **data**
  * **normalize** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – whether the normalization should be done by ‘features’ (columns)
    or ‘samples’ (rows)
* **Returns:**
  pandas.DataFrame.

Example:

```default
data = pd.DataFrame({'a': [2,5,4,3,3], 'b':[4,4,6,5,3], 'c':[4,14,8,8,9]})
result = zscore_normalization(data, normalize='samples')
result
          a         b         c
        0 -1.154701  0.577350  0.577350
        1 -0.484182 -0.665750  1.149932
        2 -1.000000  0.000000  1.000000
        3 -0.927173 -0.132453  1.059626
        4 -0.577350 -0.577350  1.154701
```

### median_polish_normalization(data, max_iter=250)

This function iteratively normalizes each sample and each feature to its
median until medians converge.

* **Parameters:**
  * **data**
  * **max_iter** ([*int*](https://docs.python.org/3/library/functions.html#int)) – number of maximum iterations to prevent infinite loop.
* **Returns:**
  pandas.DataFrame.

Example:

```default
data = pd.DataFrame({'a': [2,5,4,3,3], 'b':[4,4,6,5,3], 'c':[4,14,8,8,9]})
result = median_polish_normalization(data, max_iter = 10)
result
        a    b     c
    0  2.0  4.0   7.0
    1  5.0  7.0  10.0
    2  4.0  6.0   9.0
    3  3.0  5.0   8.0
    4  3.0  5.0   8.0
```

### quantile_normalization(data) → DataFrame

Applies quantile normalization to each column in pandas.DataFrame.

* **Parameters:**
  **data** – pandas.DataFrame with features as columns and samples as rows.
* **Returns:**
  pandas.DataFrame

Example:

```default
data = pd.DataFrame({'a': [2,5,4,3,3], 'b':[4,4,6,5,3], 'c':[4,14,8,8,9]})
result = quantile_normalization(data)
result
        a    b    c
    0  3.2  4.6  4.6
    1  4.6  3.2  8.6
    2  3.2  4.6  8.6
    3  3.2  4.6  8.6
    4  3.2  3.2  8.6
```

### linear_normalization(data, method='l1', normalize='samples') → DataFrame

This function scales input data to a unit norm. For more information visit:
[https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.normalize.html](https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.normalize.html)

* **Parameters:**
  * **data** – pandas.DataFrame with samples as rows and features as columns.
  * **method** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – norm to use to normalize each non-zero sample or non-zero feature
    (depends on axis).
  * **normalize** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – axis used to normalize the data along. If ‘samples’,
    independently normalize each sample, if ‘features’ normalize each feature.
* **Returns:**
  pandas.DataFrame

Example:

```default
data = pd.DataFrame({'a': [2,5,4,3,3], 'b':[4,4,6,5,3], 'c':[4,14,8,8,9]})
result = linear_normalization(data, method = "l1", normalize = 'samples')
result
        a         b         c
    0  0.117647  0.181818  0.093023
    1  0.294118  0.181818  0.325581
    2  0.235294  0.272727  0.186047
    3  0.176471  0.227273  0.186047
    4  0.176471  0.136364  0.209302
```
