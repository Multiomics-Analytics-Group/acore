"""Functions for filtering metabolomics feature table by biologically relevant features."""

import pandas as pd


def filter_biological_relevance(
    df: pd.DataFrame,
    rt_col: str = "Average Rt(min)",
    mz_col: str = "Average Mz",
    mz_decimals: tuple = None,
    mz_low: int = None,
    rt_dead_volume: float = None,
    save_removed: bool = True,
) -> tuple:
    """
    Cleans a DataFrame by filtering rows based on retention time and m/z.

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

    removed_list = []  # store removed rows, if save_removed=True

    if rt_dead_volume is not None:
        # Filter: rt >= 0.8
        mask_rt = df[rt_col] < rt_dead_volume
        if save_removed:
            removed_rt = df.loc[mask_rt].copy()
            removed_rt["RemovalReason"] = "BelowDeadVolume"
            removed_list.append(removed_rt)
        df = df.loc[~mask_rt]  # Remove from df

    if mz_decimals and mz_low is not None:
        # Filter: mz pattern
        low, high = mz_decimals
        mask_mz = (
            (df[mz_col] < mz_low) & (df[mz_col] % 1 >= low) & (df[mz_col] % 1 <= high)
        )
        if save_removed:
            removed_mz = df.loc[mask_mz].copy()
            removed_mz["RemovalReason"] = "NotBiologicallyRelevant"
            removed_list.append(removed_mz)
        df = df.loc[~mask_mz]

    # Combine all removed features
    removed_df = (
        pd.concat(
            removed_list, ignore_index=True
        )  # concat removed features from list to df
        if save_removed and removed_list
        else None
    )

    return (df, removed_df)
