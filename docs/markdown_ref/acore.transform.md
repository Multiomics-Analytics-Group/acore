# acore.transform namespace

## Submodules

## acore.transform.compositional module

### calc_clr(x)

Calculate the centered log-ratio (CLR) transformation.

* **Parameters:**
  **x** (*array-like*) – Input data to transform.
* **Returns:**
  CLR transformed data.
* **Return type:**
  array-like

### coda_clr(df: DataFrame) → DataFrame

Apply CoDA and CLR transformations to the numeric columns of a DataFrame.
Non-numeric columns are retained without transformation.

* **Parameters:**
  **df** (*pd.DataFrame*) – Input DataFrame with numeric and non-numeric columns.
* **Returns:**
  DataFrame with transformed numeric columns and original non-numeric columns.
* **Return type:**
  pd.DataFrame
