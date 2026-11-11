# acore.batch_correction package

### combat_batch_correction(data: DataFrame, batch_col: [str](https://docs.python.org/3/builtins/stdtypes.html#str), keep_batch_col: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → DataFrame

Corrects data for batch effects; the data must have no missing values.

This function corrects processed data for batch effects. For more information visit:
[https://github.com/epigenelabs/inmoose](https://github.com/epigenelabs/inmoose)

* **Parameters:**
  * **data** – pandas.DataFrame with samples as rows and protein identifiers as columns.
    Columns that are not numbers, such as ‘sample’ and ‘group’, are kept
    unchanged in the result.
    It must have no missing values (NaN): otherwise the call fails with
    `The count matrix contains NaN values`. Run an imputation first, for
    example `imputation_KNN`, and pass its result here, not the original
    table.
  * **batch_col** – column with the batch identifiers
  * **keep_batch_col** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool)) – keep the `batch_col` column in the result. Pass `True`
    to colour or group a later step, such as a PCA of the
    corrected data, by the batch it came from.
* **Returns:**
  pandas.DataFrame with samples as rows and protein identifiers as columns. The
  `batch_col` column is removed from the result unless `keep_batch_col`
  is `True`, so by default a later step cannot group by it; the other
  non-numeric columns are still there.

Example:

```default
result = combat_batch_correction(data, batch_col="batch")

# a table with missing values: impute first, then correct the imputed table
imputed = imputation_KNN(data, drop_cols=["sample"], group="group")
result = combat_batch_correction(imputed, batch_col="site")
```
