import logging
from typing import Iterable, Optional

import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer

DROP_COLS_DEFAULT = ("group", "sample", "subject")

logger = logging.getLogger(__name__)


def imputation_KNN(
    data: pd.DataFrame,
    drop_cols: Optional[Iterable[str]] = None,
    group: Optional[str] = None,
    cutoff=0.6,
    alone=True,
    n_neighbors=3,
):
    """
    K-Nearest Neighbors imputation for pandas dataframes with missing data. For more
    information visit `fancyimpute`_.

    .. _fancyimpute: https://github.com/iskandr/fancyimpute/blob/HEAD/fancyimpute/knn.py

    :param data: pandas dataframe with samples as rows and protein identifiers as
                 columns (with additional columns 'group', 'sample' and 'subject').
    :param str group: column label containing group identifiers, restricted to be one
                      single column for now.
    :param list drop_cols: column labels to be dropped. Final dataframe should only
                           have gene/protein/etc identifiers as columns.
    :param float cutoff: minimum fraction of valid values required to impute
                         a each column.
    :param bool alone: if True removes all columns with any missing values after initial
                       imputation.
    :return: Pandas dataframe with samples as rows and protein identifiers as columns.

    Example::

        result = imputation_KNN(data,
                    drop_cols=['group', 'sample', 'subject'],
                    group='group', cutoff=0.6, alone=True
        )
    """
    df = data.copy()

    if drop_cols is None:
        drop_cols = []

    # Prefer explicit numeric selection over deprecated/private pandas helpers
    df_numeric = df.select_dtypes(include="number").copy()

    # Columns to pass through (keep original order) which are non-numeric or the group
    # variable
    passthrough_cols = [
        c for c in df.columns if c not in df_numeric.columns and c != group
    ]
    logger.debug(f"Non-numeric columns ignored for imputation: {passthrough_cols}")
    value_cols = [c for c in df_numeric.columns if c not in drop_cols]

    if group is not None and group in df.columns:
        df_numeric[group] = df[group]
        for g in df_numeric[group].dropna().unique():
            miss_df = df_numeric.loc[df_numeric[group] == g, value_cols]
            # fraction cutoff is applied per group
            mask = miss_df.notna().mean() >= cutoff
            miss_df = miss_df.loc[:, mask].dropna(axis="columns", how="all")

            if miss_df.isna().any(axis=None):
                X = miss_df.to_numpy(dtype=np.float64, copy=False)
                X_trans = KNNImputer(n_neighbors=n_neighbors).fit_transform(X)

                dfm = pd.DataFrame(
                    X_trans, index=miss_df.index, columns=miss_df.columns
                )
                df_numeric.update(dfm)

    else:
        # Fallback: no grouping column present -> impute across all rows
        miss_df = df_numeric.loc[:, value_cols]
        mask = miss_df.notna().mean() >= cutoff
        miss_df = miss_df.loc[:, mask].dropna(axis="columns", how="all")
        if miss_df.isna().any(axis=None):
            X = miss_df  # .to_numpy(dtype=np.float64, copy=False)
            X_trans = KNNImputer(n_neighbors=n_neighbors).fit_transform(X)
            df_numeric.update(
                pd.DataFrame(X_trans, index=miss_df.index, columns=miss_df.columns)
            )

    if alone:
        df_numeric = df_numeric.dropna(axis="columns")

    # Re-attach passthrough columns (including group/sample/subject, etc.)
    if passthrough_cols:
        df_numeric = df_numeric.join(df[passthrough_cols])

    return df_numeric


def imputation_mixed_norm_KNN(
    data: pd.DataFrame,
    drop_cols: Optional[Iterable[str]] = None,
    shift: float = 1.8,
    nstd: float = 0.3,
    group: str = "group",
    cutoff: float = 0.6,
    random_state: int = 112736,
    n_neighbors: int = 3,
):
    """
    Missing values are replaced in two steps:

    1) using k-Nearest Neighbors we impute protein columns with a higher ratio of
       missing/valid values than the defined cutoff,
    2) the remaining missing values are replaced by random numbers that are drawn
       from a normal distribution.

    :param data: pandas dataframe with samples as rows and protein identifiers as
                 columns (with additional columns 'group', 'sample' and 'subject').
    :param str group: column label containing group identifiers.
    :param list drop_cols: list of column labels to be set as dataframe index.
    :param float shift: specifies the amount by which the distribution used for the
                        random numbers is shifted downwards. This is in units of the
                        standard deviation of the valid data.
    :param float nstd: defines the width of the Gaussian distribution relative to the
                       standard deviation of measured values. A value of 0.5 would mean
                       that the width of the distribution used for drawing random
                       numbers is half of the standard deviation of the data.
    :param float cutoff: minimum ratio of missing/valid values required to
                         impute in each column.
    :return: Pandas dataframe with samples as rows and protein identifiers as columns.

    Example::

        result = imputation_mixed_norm_KNN(data,
                    drop_cols=['group', 'sample', 'subject'],
                    shift = 1.8, nstd = 0.3, group='group', cutoff=0.6
        )
    """
    if isinstance(drop_cols, tuple):
        drop_cols = list(drop_cols)

    df = imputation_KNN(
        data,
        drop_cols=drop_cols,
        group=group,
        cutoff=cutoff,
        alone=False,
        n_neighbors=n_neighbors,
    )
    # drop group column for normal data imputation
    if drop_cols is None and not drop_cols:
        if group is None:
            _drop_cols = []
        else:
            _drop_cols = [group]
    else:
        _drop_cols = drop_cols.copy()
        _drop_cols.append(group)
    df = imputation_normal_distribution(
        df,
        drop_cols=_drop_cols,
        shift=shift,
        nstd=nstd,
        random_state=random_state,
    )
    return df


# ? This one is not using the grouping?
def imputation_normal_distribution(
    data: pd.DataFrame,
    drop_cols: Optional[Iterable[str]] = None,
    shift: float = 1.8,
    nstd: float = 0.3,
    random_state: int = 112736,
):
    """
    Missing values will be replaced by random numbers that are drawn from a normal
    distribution. The imputation is done for each sample (across all proteins)
    separately.
    For more information visit `replacemissingfromgaussian`_ in coxdocs from MaxQuant.
    The basic assumptions is that given normally distributed missing values in a sample
    are low abundant and are therefore replace with a downshifted minimum value.
    There is no control on if the drawn replacement values are below the absolute
    observed minimum of that protein at all, which can lead to false positives or 
    negatives in differential expression analysis that considers imputed values.

    .. _replacemissingfromgaussian: https://cox-labs.github.io/coxdocs/\
replacemissingfromgaussian.html

    :param data: pandas dataframe with samples as rows and protein identifiers as
                 columns (with additional columns 'group', 'sample' and 'subject').
    :param list drop_cols: list of column labels to be dropped from the imputation.
    :param float shift: specifies the amount by which the distribution used for the
                        random numbers is shifted downwards. This is in units of the
                        standard deviation of the valid data.
    :param float nstd: defines the width of the Gaussian distribution relative to the
                       standard deviation of measured values. A value of 0.5 would mean
                       that the width of the distribution used for drawing random
                       numbers is half of the standard deviation of the data.
    :return: Pandas dataframe with samples as rows and protein identifiers as columns.

    Example::

        result = imputation_normal_distribution(data,
                    drop_cols=['group', 'sample', 'subject'],
                    shift = 1.8, nstd = 0.3
        )
    """
    if isinstance(drop_cols, tuple):
        drop_cols = list(drop_cols)

    rng = np.random.default_rng(random_state)

    df = data.copy()
    if drop_cols is not None and drop_cols:
        df = df.drop(columns=drop_cols)

    # ToDo:  Transposing is expensive
    # data_imputed = df.T.sort_index()

    # Only iterate columns that actually have missing values
    for r in df.index[df.isna().any(axis="columns")]:
        row = df.loc[r]
        missing_mask = row.isna()
        n_missing = int(missing_mask.sum())
        if n_missing == 0:
            continue

        std = row.std(skipna=True)
        mean = row.mean(skipna=True)
        sigma = std * nstd
        mu = mean - (std * shift)

        if (
            not (
                np.isfinite(std)
                and np.isfinite(mean)
                and np.isfinite(sigma)
                and np.isfinite(mu)
            )
            or sigma <= 0
        ):
            fill_values = np.full(n_missing, 0.0, dtype=float)
        else:
            fill_values = rng.normal(mu, sigma, size=n_missing)

        df.loc[r, missing_mask] = fill_values

    return df


def imputation_zeros(
    data: pd.DataFrame,
    on_cols: Optional[Iterable[str]] = None,
    on_rows: Optional[Iterable[str]] = None,
    drop_cols: Optional[Iterable[str]] = None,
):
    """
    Replace missing values with zeros.

    :param data: DataFrame with samples as rows and features as columns.
    :param list on_cols: columns to fill with zeros. If `None`, all numeric columns are filled.
                         Non-numeric columns in "on_cols" will raise a TypeError.
    :param list on_rows: row index labels to restrict imputation to. If `None`, all rows are
                         imputed. Useful for imputing only a subset of samples (e.g. QCs,
                         blanks, controls) while leaving others untouched.
    :param list drop_cols: columns to permanently drop before imputation. If a column
                           appears in both "on_cols" and "drop_cols" it will be dropped
                           and a warning is emitted.
    :return: DataFrame with missing values in the target columns replaced by zero.

    Example:

        result = imputation_zeros(data, on_cols=['featureA', 'featureB'])
        result = imputation_zeros(data, on_rows=['QC1', 'QC2', 'blank1'])
    """
    df = data.copy()

    if drop_cols is not None and drop_cols:
        if on_cols is not None:
            overlap = set(on_cols) & set(drop_cols)
            if overlap:
                logger.warning(
                    f"Columns in both on_cols and drop_cols will be dropped, not filled: {overlap}"
                )
        df = df.drop(columns=drop_cols)

    if on_cols is None:
        for col in df.columns:
            if df[col].dtype == object:
                df[col] = pd.to_numeric(df[col], errors="ignore")
        cols = df.select_dtypes(include="number").columns.tolist()
        non_numeric = df.select_dtypes(exclude="number").columns.tolist()
        if non_numeric:
            logger.warning(f"Non-numeric columns ignored for imputation: {non_numeric}")
    else:
        cols = [c for c in on_cols if c in df.columns]
        non_numeric = [c for c in cols if not pd.api.types.is_numeric_dtype(df[c])]
        if non_numeric:
            raise TypeError(f"Non-numeric columns passed to `on_cols`: {non_numeric}")

    if on_rows is not None:
        rows = [r for r in on_rows if r in df.index]
        df.loc[rows, cols] = df.loc[rows, cols].fillna(0)
    else:
        df[cols] = df[cols].fillna(0)

    return df


def imputation_half_minimum(
    data: pd.DataFrame,
    on_cols: Optional[Iterable[str]] = None,
    on_rows: Optional[Iterable[str]] = None,
    drop_cols: Optional[Iterable[str]] = None,
):
    """
    Replace missing values with half the per-column minimum of observed values.

    :param data: DataFrame with samples as rows and features as columns.
    :param list on_cols: columns to impute. If None, all numeric columns are used.
                         Non-numeric columns in ``on_cols`` will raise a TypeError.
    :param list on_rows: row index labels to restrict imputation to. If None, all rows are
                         imputed. When provided, the per-column minimum is also computed
                         from only those rows, so each subset gets its own half-minimum
                         (e.g. blanks are imputed with half the blank-minimum).
    :param list drop_cols: columns to permanently drop before imputation. If a column
                           appears in both ``on_cols`` and ``drop_cols`` it will be dropped
                           and a warning is emitted.
    :return: DataFrame with missing values replaced by half the per-column minimum.

    Example::

        result = imputation_half_minimum(data, on_cols=['featureA', 'featureB'])
        result = imputation_half_minimum(data, on_rows=['blank1', 'blank2'])
    """
    df = data.copy()

    if drop_cols is not None and drop_cols:
        if on_cols is not None:
            overlap = set(on_cols) & set(drop_cols)
            if overlap:
                logger.warning(
                    f"Columns in both on_cols and drop_cols will be dropped, not filled: {overlap}"
                )
        df = df.drop(columns=drop_cols)

    if on_cols is None:
        for col in df.columns:
            if df[col].dtype == object:
                df[col] = pd.to_numeric(df[col], errors="ignore")
        cols = df.select_dtypes(include="number").columns.tolist()
        non_numeric = df.select_dtypes(exclude="number").columns.tolist()
        if non_numeric:
            logger.warning(f"Non-numeric columns ignored for imputation: {non_numeric}")
    else:
        cols = [c for c in on_cols if c in df.columns]
        non_numeric = [c for c in cols if not pd.api.types.is_numeric_dtype(df[c])]
        if non_numeric:
            raise TypeError(f"Non-numeric columns passed to `on_cols`: {non_numeric}")

    if on_rows is not None:
        rows = [r for r in on_rows if r in df.index]
        if len(rows) != len(on_rows):
            logger.warning(
                f"Some rows in `on_rows` were not found in the DataFrame index and will be skipped: "
                f"{set(on_rows) - set(rows)}"
            )
        subset = df.loc[rows, cols]
        all_nan = [c for c in cols if subset[c].isna().all()]
        if all_nan:
            logger.warning(
                f"Columns with no observed values in the selected rows cannot be imputed and are left as NaN: {all_nan}"
            )
        half_col_mins = subset.min(axis=0, skipna=True) / 2
        df.loc[rows, cols] = subset.fillna(half_col_mins)
    else:
        all_nan = [c for c in cols if df[c].isna().all()]
        if all_nan:
            logger.warning(
                f"Columns with no observed values cannot be imputed and are left as NaN: {all_nan}"
            )
        half_col_mins = df[cols].min(axis=0, skipna=True) / 2
        df[cols] = df[cols].fillna(half_col_mins)

    return df
