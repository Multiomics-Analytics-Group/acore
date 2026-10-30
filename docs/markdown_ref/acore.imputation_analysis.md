# acore.imputation_analysis package

### imputation_KNN(data: DataFrame, drop_cols: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None, group: [str](https://docs.python.org/3/builtins/stdtypes.html#str) | [None](https://docs.python.org/3/builtins/constants.html#None) = None, cutoff=0.6, alone=True, n_neighbors=3)

K-Nearest Neighbors imputation for pandas dataframes with missing data. For more
information visit [fancyimpute](https://github.com/iskandr/fancyimpute/blob/HEAD/fancyimpute/knn.py).

* **Parameters:**
  * **data** – pandas dataframe with samples as rows and protein identifiers as
    columns (with additional columns ‘group’, ‘sample’ and ‘subject’).
  * **group** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column label containing group identifiers, restricted to be one
    single column for now. If present, it is returned as the first
    column.
  * **drop_cols** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – column labels to be dropped. Final dataframe should only
    have gene/protein/etc identifiers as columns.
  * **cutoff** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – minimum fraction of valid values required to impute
    a each column.
  * **alone** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool)) – if True removes all columns with any missing values after initial
    imputation.
* **Returns:**
  Pandas dataframe with samples as rows and protein identifiers as columns.

Example:

```default
result = imputation_KNN(data,
            drop_cols=['group', 'sample', 'subject'],
            group='group', cutoff=0.6, alone=True
)
```

### imputation_half_minimum(data: DataFrame, drop_cols: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None, on_cols: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None, on_rows: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None, group: [str](https://docs.python.org/3/builtins/stdtypes.html#str) | [None](https://docs.python.org/3/builtins/constants.html#None) = None)

Replace missing values with half the per-column minimum of observed values.

* **Parameters:**
  * **data** – DataFrame with samples as rows and features as columns.
  * **drop_cols** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – columns to permanently drop before imputation. If a column
    appears in both `on_cols` and `drop_cols` it will be dropped
    and a warning is emitted.
  * **on_cols** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – columns to impute. If None, all numeric columns are used.
    Non-numeric columns in `on_cols` will raise a TypeError.
  * **on_rows** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – row index labels to restrict imputation to. If None, all rows are
    imputed. When provided, the per-column minimum is also computed
    from only those rows, so each subset gets its own half-minimum
    (e.g. blanks are imputed with half the blank-minimum).
  * **group** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column label containing group identifiers. If set, the minimum is
    computed within each group (and within `on_rows`, if given).
    Columns without observed values in a group are left as NaN for
    that group. Rows with a missing group identifier are not imputed.
    The group column itself is never imputed and is returned as the
    first column.
* **Returns:**
  DataFrame with missing values replaced by half the per-column minimum.

Example:

```default
result = imputation_half_minimum(data, on_cols=['featureA', 'featureB'])
result = imputation_half_minimum(data, on_rows=['blank1', 'blank2'])
result = imputation_half_minimum(data, group='group')
```

### imputation_mixed_norm_KNN(data: DataFrame, drop_cols: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None, shift: [float](https://docs.python.org/3/builtins/functions.html#float) = 1.8, nstd: [float](https://docs.python.org/3/builtins/functions.html#float) = 0.3, group: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'group', cutoff: [float](https://docs.python.org/3/builtins/functions.html#float) = 0.6, random_state: [int](https://docs.python.org/3/builtins/functions.html#int) = 112736, n_neighbors: [int](https://docs.python.org/3/builtins/functions.html#int) = 3)

Missing values are replaced in two steps:

1. using k-Nearest Neighbors we impute protein columns with a higher ratio of
   missing/valid values than the defined cutoff,
2. the remaining missing values are replaced by random numbers that are drawn
   from a normal distribution.

* **Parameters:**
  * **data** – pandas dataframe with samples as rows and protein identifiers as
    columns (with additional columns ‘group’, ‘sample’ and ‘subject’).
  * **group** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column label containing group identifiers. If present, it is
    returned as the first column.
  * **drop_cols** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – list of column labels to be set as dataframe index.
  * **shift** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – specifies the amount by which the distribution used for the
    random numbers is shifted downwards. This is in units of the
    standard deviation of the valid data.
  * **nstd** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – defines the width of the Gaussian distribution relative to the
    standard deviation of measured values. A value of 0.5 would mean
    that the width of the distribution used for drawing random
    numbers is half of the standard deviation of the data.
  * **cutoff** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – minimum ratio of missing/valid values required to
    impute in each column.
* **Returns:**
  Pandas dataframe with samples as rows and protein identifiers as columns.

Example:

```default
result = imputation_mixed_norm_KNN(data,
            drop_cols=['group', 'sample', 'subject'],
            shift = 1.8, nstd = 0.3, group='group', cutoff=0.6
)
```

### imputation_normal_distribution(data: DataFrame, drop_cols: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None, shift: [float](https://docs.python.org/3/builtins/functions.html#float) = 1.8, nstd: [float](https://docs.python.org/3/builtins/functions.html#float) = 0.3, random_state: [int](https://docs.python.org/3/builtins/functions.html#int) = 112736)

Missing values will be replaced by random numbers that are drawn from a normal
distribution. The imputation is done for each sample (across all proteins)
separately.
For more information visit [replacemissingfromgaussian](https://cox-labs.github.io/coxdocs/replacemissingfromgaussian.html) in coxdocs from MaxQuant.
The basic assumptions is that given normally distributed missing values in a sample
are low abundant and are therefore replace with a downshifted minimum value.
There is no control on if the drawn replacement values are below the absolute
observed minimum of that protein at all, which can lead to false positives or
negatives in differential expression analysis that considers imputed values.

* **Parameters:**
  * **data** – pandas dataframe with samples as rows and protein identifiers as
    columns (with additional columns ‘group’, ‘sample’ and ‘subject’).
  * **drop_cols** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – list of column labels to be dropped from the imputation.
  * **shift** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – specifies the amount by which the distribution used for the
    random numbers is shifted downwards. This is in units of the
    standard deviation of the valid data.
  * **nstd** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – defines the width of the Gaussian distribution relative to the
    standard deviation of measured values. A value of 0.5 would mean
    that the width of the distribution used for drawing random
    numbers is half of the standard deviation of the data.
* **Returns:**
  Pandas dataframe with samples as rows and protein identifiers as columns.

Example:

```default
result = imputation_normal_distribution(data,
            drop_cols=['group', 'sample', 'subject'],
            shift = 1.8, nstd = 0.3
)
```

### imputation_zeros(data: DataFrame, drop_cols: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None, on_cols: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None, on_rows: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None)

Replace missing values with zeros.

* **Parameters:**
  * **data** – DataFrame with samples as rows and features as columns.
  * **drop_cols** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – columns to permanently drop before imputation. If a column
    appears in both “on_cols” and “drop_cols” it will be dropped
    and a warning is emitted.
  * **on_cols** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – columns to fill with zeros. If None, all numeric columns are filled.
    Non-numeric columns in “on_cols” will raise a TypeError.
  * **on_rows** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – row index labels to restrict imputation to. If None, all rows are
    imputed. Useful for imputing only a subset of samples (e.g. QCs,
    blanks, controls) while leaving others untouched.
* **Returns:**
  DataFrame with missing values in the target columns replaced by zero.

### Example

result = imputation_zeros(data, on_cols=[‘featureA’, ‘featureB’])
result = imputation_zeros(data, on_rows=[‘QC1’, ‘QC2’, ‘blank1’])

### median_impute(X: DataFrame, drop_cols: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None, group: [str](https://docs.python.org/3/builtins/stdtypes.html#str) | [None](https://docs.python.org/3/builtins/constants.html#None) = None, keep_empty_features: [bool](https://docs.python.org/3/builtins/functions.html#bool) = True) → DataFrame

Replace missing values with the median of each column.

* **Parameters:**
  * **X** – Dataframe with samples as rows and features as columns.
  * **drop_cols** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – column labels to be dropped before imputation, e.g.
    non-numeric metadata columns. All remaining columns
    must be numeric.
  * **group** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column label containing group identifiers. If set, the median
    is computed within each group. Columns that are entirely missing
    within a group are left as NaN for that group. Rows with a
    missing group identifier are not imputed. The group column is
    excluded from imputation and returned as the first column.
  * **keep_empty_features** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool)) – If True, retain columns that are entirely missing
    (in all rows) and fill them with zero. Otherwise,
    drop them.
* **Returns:**
  Imputed dataframe with the original row index and all columns except
  drop_cols (and entirely missing columns if keep_empty_features
  is False).
* **Raises:**
  [**ValueError**](https://docs.python.org/3/builtins/exceptions.html#ValueError) – if columns remaining after dropping drop_cols are
  non-numeric, or group is not a column (or is in drop_cols).

## Submodules

## acore.imputation_analysis.common module

### select_imputation_targets(data: DataFrame, on_cols: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None), on_rows: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None), drop_cols: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None), group: [str](https://docs.python.org/3/builtins/stdtypes.html#str) | [None](https://docs.python.org/3/builtins/constants.html#None) = None) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[DataFrame, [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)], [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None)]

Copy data, drop drop_cols and resolve which columns and rows to impute.

* **Returns:**
  tuple of the copied DataFrame (without drop_cols), the numeric columns
  to impute and the row labels found in the index (None if on_rows
  is None).
* **Raises:**
  * [**TypeError**](https://docs.python.org/3/builtins/exceptions.html#TypeError) – if on_cols contains non-numeric columns.
  * [**ValueError**](https://docs.python.org/3/builtins/exceptions.html#ValueError) – if group is not a column of data after dropping drop_cols.

## acore.imputation_analysis.strategies module

Imputation strategies for missing values in tabular data.

The methods here depend a bit on the omics type. If the function has a group
parameter, you can impute per group.

### median_impute(X: DataFrame, drop_cols: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None, group: [str](https://docs.python.org/3/builtins/stdtypes.html#str) | [None](https://docs.python.org/3/builtins/constants.html#None) = None, keep_empty_features: [bool](https://docs.python.org/3/builtins/functions.html#bool) = True) → DataFrame

Replace missing values with the median of each column.

* **Parameters:**
  * **X** – Dataframe with samples as rows and features as columns.
  * **drop_cols** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – column labels to be dropped before imputation, e.g.
    non-numeric metadata columns. All remaining columns
    must be numeric.
  * **group** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column label containing group identifiers. If set, the median
    is computed within each group. Columns that are entirely missing
    within a group are left as NaN for that group. Rows with a
    missing group identifier are not imputed. The group column is
    excluded from imputation and returned as the first column.
  * **keep_empty_features** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool)) – If True, retain columns that are entirely missing
    (in all rows) and fill them with zero. Otherwise,
    drop them.
* **Returns:**
  Imputed dataframe with the original row index and all columns except
  drop_cols (and entirely missing columns if keep_empty_features
  is False).
* **Raises:**
  [**ValueError**](https://docs.python.org/3/builtins/exceptions.html#ValueError) – if columns remaining after dropping drop_cols are
  non-numeric, or group is not a column (or is in drop_cols).

### imputation_KNN(data: DataFrame, drop_cols: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None, group: [str](https://docs.python.org/3/builtins/stdtypes.html#str) | [None](https://docs.python.org/3/builtins/constants.html#None) = None, cutoff=0.6, alone=True, n_neighbors=3)

K-Nearest Neighbors imputation for pandas dataframes with missing data. For more
information visit [fancyimpute](https://github.com/iskandr/fancyimpute/blob/HEAD/fancyimpute/knn.py).

* **Parameters:**
  * **data** – pandas dataframe with samples as rows and protein identifiers as
    columns (with additional columns ‘group’, ‘sample’ and ‘subject’).
  * **group** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column label containing group identifiers, restricted to be one
    single column for now. If present, it is returned as the first
    column.
  * **drop_cols** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – column labels to be dropped. Final dataframe should only
    have gene/protein/etc identifiers as columns.
  * **cutoff** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – minimum fraction of valid values required to impute
    a each column.
  * **alone** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool)) – if True removes all columns with any missing values after initial
    imputation.
* **Returns:**
  Pandas dataframe with samples as rows and protein identifiers as columns.

Example:

```default
result = imputation_KNN(data,
            drop_cols=['group', 'sample', 'subject'],
            group='group', cutoff=0.6, alone=True
)
```

### imputation_mixed_norm_KNN(data: DataFrame, drop_cols: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None, shift: [float](https://docs.python.org/3/builtins/functions.html#float) = 1.8, nstd: [float](https://docs.python.org/3/builtins/functions.html#float) = 0.3, group: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = 'group', cutoff: [float](https://docs.python.org/3/builtins/functions.html#float) = 0.6, random_state: [int](https://docs.python.org/3/builtins/functions.html#int) = 112736, n_neighbors: [int](https://docs.python.org/3/builtins/functions.html#int) = 3)

Missing values are replaced in two steps:

1. using k-Nearest Neighbors we impute protein columns with a higher ratio of
   missing/valid values than the defined cutoff,
2. the remaining missing values are replaced by random numbers that are drawn
   from a normal distribution.

* **Parameters:**
  * **data** – pandas dataframe with samples as rows and protein identifiers as
    columns (with additional columns ‘group’, ‘sample’ and ‘subject’).
  * **group** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column label containing group identifiers. If present, it is
    returned as the first column.
  * **drop_cols** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – list of column labels to be set as dataframe index.
  * **shift** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – specifies the amount by which the distribution used for the
    random numbers is shifted downwards. This is in units of the
    standard deviation of the valid data.
  * **nstd** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – defines the width of the Gaussian distribution relative to the
    standard deviation of measured values. A value of 0.5 would mean
    that the width of the distribution used for drawing random
    numbers is half of the standard deviation of the data.
  * **cutoff** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – minimum ratio of missing/valid values required to
    impute in each column.
* **Returns:**
  Pandas dataframe with samples as rows and protein identifiers as columns.

Example:

```default
result = imputation_mixed_norm_KNN(data,
            drop_cols=['group', 'sample', 'subject'],
            shift = 1.8, nstd = 0.3, group='group', cutoff=0.6
)
```

### imputation_normal_distribution(data: DataFrame, drop_cols: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None, shift: [float](https://docs.python.org/3/builtins/functions.html#float) = 1.8, nstd: [float](https://docs.python.org/3/builtins/functions.html#float) = 0.3, random_state: [int](https://docs.python.org/3/builtins/functions.html#int) = 112736)

Missing values will be replaced by random numbers that are drawn from a normal
distribution. The imputation is done for each sample (across all proteins)
separately.
For more information visit [replacemissingfromgaussian](https://cox-labs.github.io/coxdocs/replacemissingfromgaussian.html) in coxdocs from MaxQuant.
The basic assumptions is that given normally distributed missing values in a sample
are low abundant and are therefore replace with a downshifted minimum value.
There is no control on if the drawn replacement values are below the absolute
observed minimum of that protein at all, which can lead to false positives or
negatives in differential expression analysis that considers imputed values.

* **Parameters:**
  * **data** – pandas dataframe with samples as rows and protein identifiers as
    columns (with additional columns ‘group’, ‘sample’ and ‘subject’).
  * **drop_cols** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – list of column labels to be dropped from the imputation.
  * **shift** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – specifies the amount by which the distribution used for the
    random numbers is shifted downwards. This is in units of the
    standard deviation of the valid data.
  * **nstd** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – defines the width of the Gaussian distribution relative to the
    standard deviation of measured values. A value of 0.5 would mean
    that the width of the distribution used for drawing random
    numbers is half of the standard deviation of the data.
* **Returns:**
  Pandas dataframe with samples as rows and protein identifiers as columns.

Example:

```default
result = imputation_normal_distribution(data,
            drop_cols=['group', 'sample', 'subject'],
            shift = 1.8, nstd = 0.3
)
```

### imputation_zeros(data: DataFrame, drop_cols: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None, on_cols: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None, on_rows: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None)

Replace missing values with zeros.

* **Parameters:**
  * **data** – DataFrame with samples as rows and features as columns.
  * **drop_cols** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – columns to permanently drop before imputation. If a column
    appears in both “on_cols” and “drop_cols” it will be dropped
    and a warning is emitted.
  * **on_cols** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – columns to fill with zeros. If None, all numeric columns are filled.
    Non-numeric columns in “on_cols” will raise a TypeError.
  * **on_rows** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – row index labels to restrict imputation to. If None, all rows are
    imputed. Useful for imputing only a subset of samples (e.g. QCs,
    blanks, controls) while leaving others untouched.
* **Returns:**
  DataFrame with missing values in the target columns replaced by zero.

### Example

result = imputation_zeros(data, on_cols=[‘featureA’, ‘featureB’])
result = imputation_zeros(data, on_rows=[‘QC1’, ‘QC2’, ‘blank1’])

### imputation_half_minimum(data: DataFrame, drop_cols: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None, on_cols: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None, on_rows: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None, group: [str](https://docs.python.org/3/builtins/stdtypes.html#str) | [None](https://docs.python.org/3/builtins/constants.html#None) = None)

Replace missing values with half the per-column minimum of observed values.

* **Parameters:**
  * **data** – DataFrame with samples as rows and features as columns.
  * **drop_cols** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – columns to permanently drop before imputation. If a column
    appears in both `on_cols` and `drop_cols` it will be dropped
    and a warning is emitted.
  * **on_cols** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – columns to impute. If None, all numeric columns are used.
    Non-numeric columns in `on_cols` will raise a TypeError.
  * **on_rows** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – row index labels to restrict imputation to. If None, all rows are
    imputed. When provided, the per-column minimum is also computed
    from only those rows, so each subset gets its own half-minimum
    (e.g. blanks are imputed with half the blank-minimum).
  * **group** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – column label containing group identifiers. If set, the minimum is
    computed within each group (and within `on_rows`, if given).
    Columns without observed values in a group are left as NaN for
    that group. Rows with a missing group identifier are not imputed.
    The group column itself is never imputed and is returned as the
    first column.
* **Returns:**
  DataFrame with missing values replaced by half the per-column minimum.

Example:

```default
result = imputation_half_minimum(data, on_cols=['featureA', 'featureB'])
result = imputation_half_minimum(data, on_rows=['blank1', 'blank2'])
result = imputation_half_minimum(data, group='group')
```
