import pandas as pd


def extract_number_missing(
    data: pd.DataFrame,
    min_valid: int,
    drop_cols: list[str] | None = None,
    group: str | None = "group",
) -> list[str]:
    """
    Counts how many valid values exist in each column and filters column labels with more
    valid values than the minimum threshold defined.

    :param data: pandas DataFrame with group as rows and protein identifier as column.
    :param str group: column label containing group identifiers.
                      If None, number of valid values is counted across all samples,
                      otherwise is counted per unique group identifier.
    :param int min_valid: minimum number of valid values to be filtered.
    :param list drop_columns: column labels to be dropped.
    :return: List of column labels above the threshold.

    Example::

        result = extract_number_missing(data, min_valid=3, drop_cols=['sample'], group='group')
    """
    if drop_cols is None:
        drop_cols = ["sample"]
    if group is None:
        groups = data.loc[:, data.notnull().sum(axis=0) >= min_valid]
    else:
        groups = data.copy()
        if len(set(drop_cols).intersection(groups.columns.tolist())) == len(drop_cols):
            groups = groups.drop(drop_cols, axis=1)
        groups = groups.set_index(group).notnull().groupby(level=0).sum()
        groups = groups[groups >= min_valid]

    groups = groups.dropna(how="all", axis=1)
    return groups.columns.unique().tolist()


def extract_percentage_missing(
    data: pd.DataFrame,
    missing_max: float,
    drop_cols: list[str] | None = None,
    group: str | None = "group",
    how: str = "all",
) -> list[str]:
    """
    Extracts ratio of missing/valid values in each column and filters column labels with
    lower ratio than the minimum threshold defined.

    :param data: pandas dataframe with group as rows and protein identifier as column.
    :param float missing_max: maximum ratio of missing/valid values to be filtered.
    :param list drop_cols: column labels to be dropped from the dataframe.
                           Default is ['sample'] if None.
    :param str group: column label containing group identifiers.
                      If None, ratio is calculated across all samples,
                      otherwise is calculated per unique group identifier.
    :param str how: define if labels with a higher percentage of missing values than the threshold
                    in any group ('any') or in all groups ('all') should be filtered
    :return: List of column labels below the threshold.

    Example::
        result = extract_percentage_missing(data, missing_max=0.3,
                                            drop_cols=['sample'], group='group')
    """
    if drop_cols is None:
        drop_cols = ["sample"]
    if group is None:
        groups = data.loc[:, data.isnull().mean() <= missing_max].columns
    else:
        groups = data.copy()
        groups = groups.drop(drop_cols, axis=1)
        groups = groups.set_index(group)
        groups = groups.isnull().groupby(level=0).mean()
        groups = groups[groups <= missing_max]
        if how == "all":
            groups = groups.dropna(how="all", axis=1).columns.unique().tolist()
        elif how == "any":
            groups = groups.dropna(how="any", axis=1).columns.unique().tolist()
        else:
            if how in groups.index:
                groups = groups.loc[how, :].dropna().index.unique().tolist()

    return groups
