from __future__ import annotations

import pandas as pd

# PyPI pycombat as described in the docstring was never used:
# https://pypi.org/project/pycombat/
# https://github.com/epigenelabs/inmoose is the update from combat.pycombat on PyPI
# from combat.pycombat import pycombat
from inmoose.pycombat import pycombat_norm

__all__ = ["combat_batch_correction"]


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
