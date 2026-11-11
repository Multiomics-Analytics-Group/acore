from __future__ import annotations

import pandas as pd

from acore.utils import check_columns_in_data

# PyPI pycombat as described in the docstring was never used:
# https://pypi.org/project/pycombat/
# https://github.com/epigenelabs/inmoose is the update from combat.pycombat on PyPI
# from combat.pycombat import pycombat

__all__ = ["combat_batch_correction"]


def combat_batch_correction(
    data: pd.DataFrame,
    batch_col: str,
    # index_cols: list[str],
) -> pd.DataFrame:
    """
    Corrects data for batch effects; the data must have no missing values.

    This function corrects processed data for batch effects. For more information visit:
    https://github.com/epigenelabs/inmoose

    :param data: pandas.DataFrame with samples as rows and protein identifiers as columns.
                 Columns that are not numbers, such as 'sample' and 'group', are kept
                 unchanged in the result.
                 It must have no missing values (NaN): otherwise the call fails with
                 ``The count matrix contains NaN values``. Run an imputation first, for
                 example ``imputation_KNN``, and pass its result here, not the original
                 table.
    :param batch_col: column with the batch identifiers
    :return: pandas.DataFrame with samples as rows and protein identifiers as columns. The
             ``batch_col`` column is removed from the result, so a later step cannot
             group by it; the other non-numeric columns are still there.

    Example::

        result = combat_batch_correction(data, batch_col="batch")

        # a table with missing values: impute first, then correct the imputed table
        imputed = imputation_KNN(data, drop_cols=["sample"], group="group")
        result = combat_batch_correction(imputed, batch_col="site")

    """
    check_columns_in_data(data, batch_col=batch_col)
    from inmoose.pycombat import pycombat_norm  # slow import, so only when needed

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
