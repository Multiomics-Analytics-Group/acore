"""Module normalization for analysis. This module provides functions to normalize data for analysis.

The higher-level convience functions are `normalize_data` and `normalize_data_by_group`.
Combat normalization was added using the inmoose package.

The actual normalization functions are in strategies.py.
"""

from __future__ import annotations

import pandas as pd

# PyPI pycombat as described in the docstring was never used:
# https://pypi.org/project/pycombat/
# https://github.com/epigenelabs/inmoose is the update from combat.pycombat on PyPI
# from combat.pycombat import pycombat
from inmoose.pycombat import pycombat_norm

from .strategies import (
    linear_normalization,
    median_normalization,
    median_polish_normalization,
    median_zero_normalization,
    quantile_normalization,
    zscore_normalization,
)

__all__ = ["combat_batch_correction", "normalize_data", "normalize_data_per_group"]


def combat_batch_correction(
    data: pd.DataFrame,
    batch_col: str,
    # index_cols: list[str],
) -> pd.DataFrame:
    """
    This function corrects processed data for batch effects. For more information visit:
    https://github.com/epigenelabs/inmoose

    :param data: pandas.DataFrame with samples as rows and protein identifiers as columns.
    :param batch_col: column with the batch identifiers
    :return: pandas.DataFrame with samples as rows and protein identifiers as columns.
    Example::
        result = combat_batch_correction(
                    data,
                    batch_col="batch",
                    index_cols=["subject", "sample", "group"],
                )

    """
    # :param index_cols: list of columns that don't need to be corrected (i.e group)
    df_corrected = pd.DataFrame()
    # index_cols = [c for c in index_cols if c != batch_col]
    # data = data.set_index(index_cols)  # ? should this not be provided directly as data
    df = data.drop(batch_col, axis=1)
    df_numeric = df.select_dtypes("number")
    num_batches = len(data[batch_col].unique())
    if df_numeric.empty:
        raise ValueError("No numeric columns found in data.")
    if not num_batches > 1:
        raise ValueError("Only one batch found in data.")
    info_cols = df.columns.difference(df_numeric.columns)
    df_corrected = pd.DataFrame(
        pycombat_norm(df_numeric.T, data[batch_col]).T,
        index=df.index,
    )
    df_corrected = df_corrected.join(df[info_cols])
    # df_corrected = df_corrected  # .reset_index()  # ? would also not reset index here

    return df_corrected


def normalize_data_per_group(
    data: pd.DataFrame,
    group: str | int | list[str | int],
    method: str = "median",
    normalize: str = None,
) -> pd.DataFrame:
    """
    This function normalizes the data by group using the selected method

    :param data: DataFrame with the data to be normalized (samples x features)
    :param group_col: Column containing the groups, passed to pandas.DataFrame.groupby
    :param str method: normalization method to choose among: median_polish, median,
                        quantile, linear
    :param str normalize: whether the normalization should be done by 'features' (columns) or 'samples' (rows) (default None)
    :return: pandas.DataFrame.

    Example::

        result = normalize_data_per_group(data, group='group' method='median')
    """
    ndf = pd.DataFrame(columns=data.columns)
    for _, gdf in data.groupby(group):
        norm_group = normalize_data(gdf, method=method, normalize=normalize)
        ndf = ndf.append(norm_group)

    return ndf


def normalize_data(
    data: pd.DataFrame,
    method: str = "median",
    normalize: str = None,
):
    """
    This function normalizes the data using the selected method. Normalizes only nummeric
    data, but keeps the non-numeric columns in the output DataFrame.

    :param data: DataFrame with the data to be normalized (samples x features)
    :param str method: normalization method to choose among: median (default),
                       median_polish, median_zero, quantile, linear, zscore
    :param str normalize: whether the normalization should be done by 'features' (columns)
                          or 'samples' (rows) (default None)
    :return: pandas.DataFrame.

    Example::

        result = normalize_data(data, method='median_polish')
    """
    numeric_cols = data.select_dtypes(
        include=["int64", "float64"]
    )  # ! too restrictive?
    non_numeric_cols = data.select_dtypes(
        exclude=["int64", "float64"]
    )  # ! too restrictive?
    if numeric_cols.empty:
        raise ValueError("No numeric columns found in data.")

    if method == "median_polish":
        norm_data = median_polish_normalization(numeric_cols, max_iter=250)
    elif method == "median_zero":
        norm_data = median_zero_normalization(numeric_cols, normalize)
    elif method == "median":
        norm_data = median_normalization(numeric_cols, normalize)
    elif method == "quantile":
        norm_data = quantile_normalization(numeric_cols)
    elif method == "linear":
        norm_data = linear_normalization(numeric_cols, method="l1", normalize=normalize)
    elif method == "zscore":
        norm_data = zscore_normalization(numeric_cols, normalize)
    else:
        raise ValueError(
            "Invalid normalization method. Should be one of:"
            " 'median', 'median_polish', 'median_zero', 'quantile', 'linear', 'zscore'"
        )

    if non_numeric_cols is not None and not non_numeric_cols.empty:
        norm_data = norm_data.join(non_numeric_cols)

    return norm_data
