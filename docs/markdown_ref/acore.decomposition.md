# acore.decomposition package

## Submodules

## acore.decomposition.pca module

### run_pca(df_wide: DataFrame, n_components: [int](https://docs.python.org/3/library/functions.html#int) = 2) → [tuple](https://docs.python.org/3/library/stdtypes.html#tuple)[DataFrame, PCA]

Run PCA on DataFrame using 

```
:class:`sklearn.decomposition.PCA`_
```

.

* **Parameters:**
  * **df** (*pd.DataFrame*) – DataFrame in wide format to fit features on.
  * **n_components** ([*int*](https://docs.python.org/3/library/functions.html#int) *,* *optional*) – Number of Principal Components to fit, by default 2
* **Returns:**
  principal components of DataFrame with same indices as in original DataFrame,
  and fitted PCA model of sklearn
* **Return type:**
  Tuple[pd.DataFrame, sklearn.decomposition.PCA]

## acore.decomposition.umap module

Run UMAP on DataFrame and return result.

### run_umap(X_scaled, y, random_state=42) → DataFrame
