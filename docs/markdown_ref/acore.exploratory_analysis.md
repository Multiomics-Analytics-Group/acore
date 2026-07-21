# acore.exploratory_analysis package

### get_histogram_series(s: Series, bins: [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) → Series

### calculate_coefficient_variation(df: DataFrame) → Series

Compute the coefficient of variation (CV) for each column in a DataFrame.

The coefficient of variation is defined as the ratio of the standard deviation
to the mean, expressed as a percentage. This function uses the biased standard
deviation (normalization by N) as implemented in scipy.stats.variation.

* **Parameters:**
  **df** ([*pandas.DataFrame*](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html#pandas.DataFrame)) – Input DataFrame containing numeric values (e.g., log2-transformed data).
  Each column will be processed independently.
* **Returns:**
  Series containing the coefficient of variation (in percent) for each column.
  The index corresponds to the columns of the input DataFrame.
* **Return type:**
  [pandas.Series](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.Series.html#pandas.Series)

#### SEE ALSO
`scipy.stats.variation`
: Function used to compute the coefficient of variation.

### Examples

```pycon
>>> import pandas as pd
>>> import numpy as np
>>> df = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})
>>> calculate_coefficient_variation(df)
A   40.825
B   16.330
Name: coef_of_var, dtype: float64
```

### calculate_coef_of_var_and_mean(df: DataFrame) → DataFrame

Calculate coefficient of variation and mean for each column in the dataframe,
the mean calculated on both log2 and linear scale.

* **Parameters:**
  **df** (*pd.DataFrame*) – The input dataframe containing the linear values (non-log transformed).
* **Returns:**
  A dataframe with columns ‘mean_log2’, ‘mean’ and ‘coef_of_var’ for each column
  in the input dataframe.
* **Return type:**
  pd.DataFrame

### get_coefficient_variation(data: DataFrame, drop_columns: [list](https://docs.python.org/3/library/stdtypes.html#list)[[str](https://docs.python.org/3/library/stdtypes.html#str)] | [None](https://docs.python.org/3/library/constants.html#None) = None, group: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'group')

Extracts the coefficients of variation in each group.

* **Parameters:**
  * **data** – pandas.DataFrame with samples as rows and protein identifiers as columns
    (with additional columns ‘group’, ‘sample’ and ‘subject’). The values
    should be the original intensities for massspectrometry-based
    measurements.
  * **drop_columns** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – column labels to be dropped from the dataframe
  * **group** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column label containing group identifiers.
* **Returns:**
  Pandas dataframe with columns ‘name’ (protein identifier),
  ‘x’ (coefficient of variation), ‘y’ (mean) and ‘group’.

Example:

```default
result = get_coefficient_variation(data, drop_columns=['sample', 'subject'], group='group')
```

### extract_number_missing(data, min_valid, drop_cols=['sample'], group='group')

Counts how many valid values exist in each column and filters column labels with more
valid values than the minimum threshold defined.

* **Parameters:**
  * **data** – pandas DataFrame with group as rows and protein identifier as column.
  * **group** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column label containing group identifiers.
    If None, number of valid values is counted across all samples,
    otherwise is counted per unique group identifier.
  * **min_valid** ([*int*](https://docs.python.org/3/library/functions.html#int)) – minimum number of valid values to be filtered.
  * **drop_columns** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – column labels to be dropped.
* **Returns:**
  List of column labels above the threshold.

Example:

```default
result = extract_number_missing(data, min_valid=3, drop_cols=['sample'], group='group')
```

### extract_percentage_missing(data, missing_max, drop_cols=['sample'], group='group', how='all')

Extracts ratio of missing/valid values in each column and filters column labels with
lower ratio than the minimum threshold defined.

* **Parameters:**
  * **data** – pandas dataframe with group as rows and protein identifier as column.
  * **group** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column label containing group identifiers.
    If None, ratio is calculated across all samples,
    otherwise is calculated per unique group identifier.
  * **missing_max** ([*float*](https://docs.python.org/3/library/functions.html#float)) – maximum ratio of missing/valid values to be filtered.
  * **how** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – define if labels with a higher percentage of missing values than the threshold
    in any group (‘any’) or in all groups (‘all’) should be filtered
* **Returns:**
  List of column labels below the threshold.

Example::
: result = extract_percentage_missing(data, missing_max=0.3,
  : drop_cols=[‘sample’], group=’group’)

### run_pca(data, drop_cols=['sample', 'subject'], group='group', annotation_cols=['sample'], components=2, dropna=True)

Performs principal component analysis and returns the values of each component for each sample
and each protein, and the loadings for each protein.

For information visit
[https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html).

* **Parameters:**
  * **data** – pandas dataframe with samples as rows and protein identifiers as columns
    (with additional columns ‘group’, ‘sample’ and ‘subject’).
  * **drop_cols** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – column labels to be dropped from the dataframe.
  * **group** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column label containing group identifiers.
  * **annotation_cols** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – list of columns to be added in the scatter plot annotation
  * **components** ([*int*](https://docs.python.org/3/library/functions.html#int)) – number of components to keep.
  * **dropna** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – if True removes all columns with any missing values.
* **Returns:**
  tuple: 1) three pandas dataframes: components, loadings and variance; 2)
  xaxis and yaxis titles with components loadings for plotly.

Example:

```default
result = run_pca(data, drop_cols=['sample', 'subject'], group='group',
                 components=2, dropna=True)
```

### run_tsne(data, drop_cols=['sample', 'subject'], group='group', annotation_cols=['sample'], components=2, perplexity=40, max_iter=1000, init='pca', dropna=True)

Performs t-distributed Stochastic Neighbor Embedding analysis.

For more information visit
[https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html).

* **Parameters:**
  * **data** – pandas dataframe with samples as rows and protein identifiers as columns
    (with additional columns ‘group’, ‘sample’ and ‘subject’).
  * **drop_cols** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – column labels to be dropped from the dataframe.
  * **group** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column label containing group identifiers.
  * **components** ([*int*](https://docs.python.org/3/library/functions.html#int)) – dimension of the embedded space.
  * **annotation_cols** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – list of columns to be added in the scatter plot annotation
  * **perplexity** ([*int*](https://docs.python.org/3/library/functions.html#int)) – related to the number of nearest neighbors that is used
    in other manifold learning algorithms.
    Consider selecting a value between 5 and 50.
  * **max_iter** ([*int*](https://docs.python.org/3/library/functions.html#int)) – maximum number of iterations for the optimization (at least 250).
  * **init** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – initialization of embedding (‘random’, ‘pca’ or
    numpy array of shape n_samples x n_components).
  * **dropna** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – if True removes all columns with any missing values.
* **Returns:**
  Two dictionaries:
  1) pandas dataframe with embedding vectors,
  2) xaxis and yaxis titles for plotly.

Example:

```default
result = run_tsne(data,
                  drop_cols=['sample', 'subject'],
                  group='group',
                  components=2,
                  perplexity=40,
                  max_iter=1000,
                  init='pca',
                  dropna=True
                )
```

### run_umap(data, drop_cols=['sample', 'subject'], group='group', annotation_cols=['sample'], n_neighbors=10, min_dist=0.3, metric='cosine', dropna=True)

Performs Uniform Manifold Approximation and Projection.

For more information vist [https://umap-learn.readthedocs.io](https://umap-learn.readthedocs.io).

* **Parameters:**
  * **data** – pandas dataframe with samples as rows and protein identifiers as columns
    (with additional columns ‘group’, ‘sample’ and ‘subject’).
  * **drop_cols** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – column labels to be dropped from the dataframe.
  * **group** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column label containing group identifiers.
  * **annotation_cols** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – list of columns to be added in the scatter plot annotation
  * **n_neighbors** ([*int*](https://docs.python.org/3/library/functions.html#int)) – number of neighboring points used
    in local approximations of manifold structure.
  * **min_dist** ([*float*](https://docs.python.org/3/library/functions.html#float)) – controls how tightly the embedding is allowed compress points together.
  * **metric** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – metric used to measure distance in the input space.
  * **dropna** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – if True removes all columns with any missing values.
* **Returns:**
  Two dictionaries:
  1) pandas dataframe with embedding of the training data in low-dimensional space,
  2) xaxis and yaxis titles for plotly.

Example:

```default
result = run_umap(data,
                  drop_cols=['sample', 'subject'],
                  group='group',
                  n_neighbors=10,
                  min_dist=0.3,
                  metric='cosine',
                  dropna=True
                )
```
