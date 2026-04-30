# acore.imputation_analysis package

### imputation_KNN(data, drop_cols=['group', 'sample', 'subject'], group='group', cutoff=0.6, alone=True)

k-Nearest Neighbors imputation for pandas dataframes with missing data. For more information visit [https://github.com/iskandr/fancyimpute/blob/master/fancyimpute/knn.py](https://github.com/iskandr/fancyimpute/blob/master/fancyimpute/knn.py).

* **Parameters:**
  * **data** – pandas dataframe with samples as rows and protein identifiers as columns (with additional columns ‘group’, ‘sample’ and ‘subject’).
  * **group** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column label containing group identifiers.
  * **drop_cols** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – column labels to be dropped. Final dataframe should only have gene/protein/etc identifiers as columns.
  * **cutoff** ([*float*](https://docs.python.org/3/library/functions.html#float)) – minimum ratio of missing/valid values required to impute in each column.
  * **alone** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – if True removes all columns with any missing values.
* **Returns:**
  Pandas dataframe with samples as rows and protein identifiers as columns.

Example:

```default
result = imputation_KNN(data, drop_cols=['group', 'sample', 'subject'], group='group', cutoff=0.6, alone=True)
```

### imputation_mixed_norm_KNN(data, index_cols=['group', 'sample', 'subject'], shift=1.8, nstd=0.3, group='group', cutoff=0.6)

Missing values are replaced in two steps: 1) using k-Nearest Neighbors we impute protein columns with a higher ratio of missing/valid values than the defined cutoff,     2) the remaining missing values are replaced by random numbers that are drawn from a normal distribution.

* **Parameters:**
  * **data** – pandas dataframe with samples as rows and protein identifiers as columns (with additional columns ‘group’, ‘sample’ and ‘subject’).
  * **group** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column label containing group identifiers.
  * **index_cols** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – list of column labels to be set as dataframe index.
  * **shift** ([*float*](https://docs.python.org/3/library/functions.html#float)) – specifies the amount by which the distribution used for the random numbers is shifted downwards. This is in units of the                         standard deviation of the valid data.
  * **nstd** ([*float*](https://docs.python.org/3/library/functions.html#float)) – defines the width of the Gaussian distribution relative to the standard deviation of measured values.                         A value of 0.5 would mean that the width of the distribution used for drawing random numbers is half of the standard deviation of the data.
  * **cutoff** ([*float*](https://docs.python.org/3/library/functions.html#float)) – minimum ratio of missing/valid values required to impute in each column.
* **Returns:**
  Pandas dataframe with samples as rows and protein identifiers as columns.

Example:

```default
result = imputation_mixed_norm_KNN(data, index_cols=['group', 'sample', 'subject'], shift = 1.8, nstd = 0.3, group='group', cutoff=0.6)
```

### imputation_normal_distribution(data, index_cols=['group', 'sample', 'subject'], shift=1.8, nstd=0.3)

Missing values will be replaced by random numbers that are drawn from a normal distribution. The imputation is done for each sample (across all proteins) separately.
For more information visit [http://www.coxdocs.org/doku.php?id=perseus:user:activities:matrixprocessing:imputation:replacemissingfromgaussian](http://www.coxdocs.org/doku.php?id=perseus:user:activities:matrixprocessing:imputation:replacemissingfromgaussian).

* **Parameters:**
  * **data** – pandas dataframe with samples as rows and protein identifiers as columns (with additional columns ‘group’, ‘sample’ and ‘subject’).
  * **index_cols** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – list of column labels to be set as dataframe index.
  * **shift** ([*float*](https://docs.python.org/3/library/functions.html#float)) – specifies the amount by which the distribution used for the random numbers is shifted downwards. This is in units of the standard deviation of the valid data.
  * **nstd** ([*float*](https://docs.python.org/3/library/functions.html#float)) – defines the width of the Gaussian distribution relative to the standard deviation of measured values. A value of 0.5 would mean that the width of the distribution used for drawing random numbers is half of the standard deviation of the data.
* **Returns:**
  Pandas dataframe with samples as rows and protein identifiers as columns.

Example:

```default
result = imputation_normal_distribution(data, index_cols=['group', 'sample', 'subject'], shift = 1.8, nstd = 0.3)
```
