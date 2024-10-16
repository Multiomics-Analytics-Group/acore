import pandas as pd


def transform_DataFrame(X: pd.DataFrame, fct: callable) -> pd.DataFrame:
    """Set index and columns of a DataFrame after applying a callable
    which might only return a numpy array.

    Parameters
    ----------
    X : pd.DataFrame
        Original DataFrame to be transformed
    fct : callable
        Callable to be applied to every element in the DataFrame.

    Returns
    -------
    pd.DataFrame
        Transformed DataFrame
    """
    ret = fct(X)
    ret = pd.DataFrame(ret, index=X.index, columns=X.columns)
    return ret
