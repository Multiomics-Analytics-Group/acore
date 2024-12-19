import pandas as pd
import sklearn.decomposition


# ! also a version in exploratory analysis
def run_pca(
    df_wide: pd.DataFrame, n_components: int = 2
) -> tuple[pd.DataFrame, sklearn.decomposition.PCA]:
    """Run PCA on DataFrame and return result.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame in wide format to fit features on.
    n_components : int, optional
        Number of Principal Components to fit, by default 2

    Returns
    -------
    Tuple[pd.DataFrame, PCA]
        principal compoments of DataFrame with same indices as in original DataFrame,
        and fitted PCA model of sklearn
    """
    n_comp_max = None
    if n_components is not None:
        n_comp_max = min(df_wide.shape)
        n_comp_max = min(n_comp_max, n_components)
    pca = sklearn.decomposition.PCA(n_components=n_comp_max)
    pcs = pca.fit_transform(df_wide)
    cols = [
        f"principal component {i+1} ({var_explained*100:.2f} %)"
        for i, var_explained in enumerate(pca.explained_variance_ratio_)
    ]
    pcs = pd.DataFrame(pcs, index=df_wide.index, columns=cols)
    return pcs, pca
