"""Module for filtering metabolomics feature table."""

import pandas as pd

from filter_biological_relevance import filter_biological_relevance
from make_numeric import convert_to_numeric


def filter_mz_rt(
    df: pd.DataFrame,
    rt_col: str = "Average Rt(min)",
    mz_col: str = "Average Mz",
    mz_decimals: tuple = None,
    mz_low: int = None,
    rt_dead_volume: float = None,
    save_removed: bool = True,
    print_na_summary: bool = True,
) -> tuple:
    """
    This function filters rows from a data frame based on retention time and m/z and checks data types.
    If specified by the user, it evaluates each row on whether there are NaN values and prints
    a statement for each.

    Data types:
        Tries to convert mz and RT columns to numeric.

    M/z filtering:
        Filters out all features that have m/z decimals in a given range that are below
        a certain m/z value.

    RT filtering:
        Filters out all features that have RT below a certain number (in minutes). Corresponds
        to dead volume.

    If save_removed==True, the removed features are saved to a separate DataFrame that also
    contains a new column, "RemovalReason", indicating the removal reason.



    :param df: Input DataFrame.
    :type df: pd.DataFrame
    :param rt_col: Column name for retention time.
    :type rt_col: str
    :param mz_col: Column name for m/z values.
    :type mz_col: str
    :param mz_decimals: Tuple specifying fractional range to filter for m/z (e.g., (0.3, 0.9)).
    :type mz_decimals: tuple, optional
    :param mz_low: Threshold for low m/z values.
    :type mz_low: int, optional
    :param rt_dead_volume: Minimum retention time threshold.
    :type rt_dead_volume: float, optional
    :param save_removed: Whether to return a DataFrame containing removed rows.
    :type save_removed: bool, default=True

    :return: Tuple containing the cleaned DataFrame and optionally the removed features DataFrame.
    :rtype: tuple
        - If `save_removed=True`: `(cleaned_df, removed_df)`
        - If `save_removed=False`: `(cleaned_df, None)`
    """

    cols_to_convert = [rt_col, mz_col]
    numeric_df = convert_to_numeric(
        df, cols_to_convert=cols_to_convert, print_na_summary=True
    )

    rt_dead_volume = 0.8  # in minutes
    mz_decimals = [0.3, 0.9]
    mz_low = 600

    filtered_df, removed_features = filter_biological_relevance(
        numeric_df,
        rt_col=rt_col,
        mz_col=mz_col,
        mz_decimals=mz_decimals,
        mz_low=mz_low,
        rt_dead_volume=rt_dead_volume,
        save_removed=True,
    )

    if rt_dead_volume is not None:
        print("Filtering based on RT completed.")
    if mz_decimals is not None and mz_low is not None:
        print("Filtering based on m/z was completed.")

    return filtered_df, removed_features
