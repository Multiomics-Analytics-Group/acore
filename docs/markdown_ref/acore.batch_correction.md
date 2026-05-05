# acore.batch_correction package

### combat_batch_correction(data: DataFrame, batch_col: [str](https://docs.python.org/3/library/stdtypes.html#str)) → DataFrame

This function corrects processed data for batch effects. For more information visit:
[https://github.com/epigenelabs/inmoose](https://github.com/epigenelabs/inmoose)

* **Parameters:**
  * **data** – pandas.DataFrame with samples as rows and protein identifiers as columns.
  * **batch_col** – column with the batch identifiers
* **Returns:**
  pandas.DataFrame with samples as rows and protein identifiers as columns.

Example:

```default
result = combat_batch_correction(
            data,
            batch_col="batch",
            index_cols=["subject", "sample", "group"],
        )
```
