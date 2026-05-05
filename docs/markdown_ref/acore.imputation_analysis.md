# acore.imputation_analysis package

### imputation_KNN(data: DataFrame, drop_cols: [Iterable](https://docs.python.org/3/library/typing.html#typing.Iterable)[[str](https://docs.python.org/3/library/stdtypes.html#str)] | [None](https://docs.python.org/3/library/constants.html#None) = None, group: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None, cutoff=0.6, alone=True, n_neighbors=3)

K-Nearest Neighbors imputation for pandas dataframes with missing data. For more
information visit [fancyimpute](https://github.com/iskandr/fancyimpute/blob/HEAD/fancyimpute/knn.py).

* **Parameters:**
  * **data** – pandas dataframe with samples as rows and protein identifiers as
    columns (with additional columns ‘group’, ‘sample’ and ‘subject’).
  * **group** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column label containing group identifiers, restricted to be one
    single column for now.
  * **drop_cols** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – column labels to be dropped. Final dataframe should only
    have gene/protein/etc identifiers as columns.
  * **cutoff** ([*float*](https://docs.python.org/3/library/functions.html#float)) – minimum fraction of valid values required to impute
    a each column.
  * **alone** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – if True removes all columns with any missing values after initial
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

### imputation_mixed_norm_KNN(data: DataFrame, drop_cols: [Iterable](https://docs.python.org/3/library/typing.html#typing.Iterable)[[str](https://docs.python.org/3/library/stdtypes.html#str)] | [None](https://docs.python.org/3/library/constants.html#None) = None, shift: [float](https://docs.python.org/3/library/functions.html#float) = 1.8, nstd: [float](https://docs.python.org/3/library/functions.html#float) = 0.3, group: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'group', cutoff: [float](https://docs.python.org/3/library/functions.html#float) = 0.6, random_state: [int](https://docs.python.org/3/library/functions.html#int) = 112736, n_neighbors: [int](https://docs.python.org/3/library/functions.html#int) = 3)

Missing values are replaced in two steps:

1. using k-Nearest Neighbors we impute protein columns with a higher ratio of
   missing/valid values than the defined cutoff,
2. the remaining missing values are replaced by random numbers that are drawn
   from a normal distribution.

* **Parameters:**
  * **data** – pandas dataframe with samples as rows and protein identifiers as
    columns (with additional columns ‘group’, ‘sample’ and ‘subject’).
  * **group** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column label containing group identifiers.
  * **drop_cols** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – list of column labels to be set as dataframe index.
  * **shift** ([*float*](https://docs.python.org/3/library/functions.html#float)) – specifies the amount by which the distribution used for the
    random numbers is shifted downwards. This is in units of the
    standard deviation of the valid data.
  * **nstd** ([*float*](https://docs.python.org/3/library/functions.html#float)) – defines the width of the Gaussian distribution relative to the
    standard deviation of measured values. A value of 0.5 would mean
    that the width of the distribution used for drawing random
    numbers is half of the standard deviation of the data.
  * **cutoff** ([*float*](https://docs.python.org/3/library/functions.html#float)) – minimum ratio of missing/valid values required to
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

### imputation_normal_distribution(data: DataFrame, drop_cols: [Iterable](https://docs.python.org/3/library/typing.html#typing.Iterable)[[str](https://docs.python.org/3/library/stdtypes.html#str)] | [None](https://docs.python.org/3/library/constants.html#None) = None, shift: [float](https://docs.python.org/3/library/functions.html#float) = 1.8, nstd: [float](https://docs.python.org/3/library/functions.html#float) = 0.3, random_state: [int](https://docs.python.org/3/library/functions.html#int) = 112736)

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
  * **drop_cols** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – list of column labels to be dropped from the imputation.
  * **shift** ([*float*](https://docs.python.org/3/library/functions.html#float)) – specifies the amount by which the distribution used for the
    random numbers is shifted downwards. This is in units of the
    standard deviation of the valid data.
  * **nstd** ([*float*](https://docs.python.org/3/library/functions.html#float)) – defines the width of the Gaussian distribution relative to the
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
