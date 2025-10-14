# preprocessing
from sklearn.preprocessing import FunctionTransformer
import numpy as np
from sklearn.pipeline import make_pipeline
import pandas as pd

def calc_clr(x):
    """
    Calculate the centered log-ratio (CLR) transformation.
    
    Parameters
    ----------
    x : array-like
        Input data to transform.

    Returns
    -------
    array-like
        CLR transformed data.
    """
    # replace zeros with small 
    new_x = np.where(x > 0, x, 1e-10)
    return np.log(new_x) - np.mean(np.log(new_x), axis=0)


def coda_clr(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply CoDA and CLR transformations to the numeric columns of a DataFrame.
    Non-numeric columns are retained without transformation.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with numeric and non-numeric columns.

    Returns
    -------
    pd.DataFrame
        DataFrame with transformed numeric columns and original non-numeric columns.
    """
    # Define the CoDA and CLR transformers
    coda = FunctionTransformer(
        lambda x: x / x.sum(axis=0), 
        feature_names_out="one-to-one"
    ).set_output(transform="pandas")
    clr = FunctionTransformer(
        calc_clr, 
        feature_names_out="one-to-one"
    ).set_output(transform="pandas")
    # Create a pipeline with CoDA and CLR transformations
    pipe = make_pipeline(coda, clr).set_output(transform="pandas")
    # Apply the pipeline to numeric columns and concatenate with non-numeric columns
    transformed = pipe.fit_transform(df.select_dtypes(include='number'))
    return pd.concat([df.select_dtypes(include='object'), transformed], axis=1)