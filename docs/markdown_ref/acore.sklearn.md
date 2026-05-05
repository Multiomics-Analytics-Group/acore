# acore.sklearn package

### transform_DataFrame(X: DataFrame, fct: callable) → DataFrame

Set index and columns of a DataFrame after applying a callable
which might only return a numpy array.

* **Parameters:**
  * **X** (*pd.DataFrame*) – Original DataFrame to be transformed
  * **fct** (*callable*) – Callable to be applied to every element in the DataFrame.
* **Returns:**
  Transformed DataFrame
* **Return type:**
  pd.DataFrame
