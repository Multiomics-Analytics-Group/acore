"""Imputation strategies for missing values in tabular data.

The methods here depend a bit on the omics type. If the function has a group
parameter, you can impute per group.
"""

import logging
from collections.abc import Iterable

import numpy as np
import pandas as pd
import sklearn.impute
from sklearn.impute import KNNImputer

from . import common

logger = logging.getLogger(__name__)


def median_impute(
    X: pd.DataFrame,
    drop_cols: Iterable[str] | None = None,
    group: str | None = None,
    keep_empty_features: bool = True,
) -> pd.DataFrame:
    """Replace missing values with the median of each column.

    :param X: Dataframe with samples as rows and features as columns.
    :param list drop_cols: column labels to be dropped before imputation, e.g.
                           non-numeric metadata columns. All remaining columns
                           must be numeric.
    :param str group: column label containing group identifiers. If set, the median
                      is computed within each group. Columns that are entirely missing
                      within a group are left as NaN for that group. Rows with a
                      missing group identifier are not imputed. The group column is
                      excluded from imputation and returned as the first column.
    :param bool keep_empty_features: If True, retain columns that are entirely missing
                                     (in all rows) and fill them with zero. Otherwise,
                                     drop them.
    :return: Imputed dataframe with the original row index and all columns except
             `drop_cols` (and entirely missing columns if `keep_empty_features`
             is False).
    :raises ValueError: if columns remaining after dropping `drop_cols` are
                        non-numeric, or `group` is not a column (or is in `drop_cols`).
    """
    df = X.drop(columns=drop_cols) if drop_cols else X
    if group is not None:
        if group not in df.columns:
            raise ValueError(
                f"`group` column {group!r} not found in columns (or in `drop_cols`)."
            )
        labels = df[group]
        df = df.drop(columns=group)
    non_numeric = df.select_dtypes(exclude="number").columns.tolist()
    if non_numeric:
        raise ValueError(
            f"Non-numeric columns found, add them to `drop_cols`: {non_numeric}"
        )
    if df.shape[1] == 0:
        imputed = df.copy()
    elif group is None:
        imputed = (
            sklearn.impute.SimpleImputer(
                strategy="median", keep_empty_features=keep_empty_features
            )
            .set_output(transform="pandas")
            .fit_transform(df)
        )
    else:
        imputed = df.astype(float).fillna(df.groupby(labels).transform("median"))
        empty = df.columns[df.isna().all()]
        if keep_empty_features:
            imputed[empty] = 0.0
        else:
            imputed = imputed.drop(columns=empty)

    if group is not None:
        imputed.insert(0, group, labels)
    return imputed


def imputation_KNN(
    data: pd.DataFrame,
    drop_cols: Iterable[str] | None = None,
    group: str | None = None,
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
                      single column for now. If present, it is returned as the first
                      column.
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

    if group is not None and group in df_numeric.columns:
        df_numeric = df_numeric[[group, *df_numeric.columns.drop(group)]]

    return df_numeric


def imputation_mixed_norm_KNN(
    data: pd.DataFrame,
    drop_cols: Iterable[str] | None = None,
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
    _drop_cols = list(drop_cols) if drop_cols else []
    if group is not None and group not in _drop_cols:
        _drop_cols.append(group)
    df = imputation_normal_distribution(
        df,
        drop_cols=_drop_cols,
        shift=shift,
        nstd=nstd,
        random_state=random_state,
    )
    return df


# This one is not using the grouping as it is supposed to be per sample.
def imputation_normal_distribution(
    data: pd.DataFrame,
    drop_cols: Iterable[str] | None = None,
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
    drop_cols: Iterable[str] | None = None,
    on_cols: Iterable[str] | None = None,
    on_rows: Iterable[str] | None = None,
):
    """
    Replace missing values with zeros.

    :param data: DataFrame with samples as rows and features as columns.
    :param list drop_cols: columns to permanently drop before imputation. If a column
                           appears in both "on_cols" and "drop_cols" it will be dropped
                           and a warning is emitted.
    :param list on_cols: columns to fill with zeros. If `None`, all numeric columns are filled.
                         Non-numeric columns in "on_cols" will raise a TypeError.
    :param list on_rows: row index labels to restrict imputation to. If `None`, all rows are
                         imputed. Useful for imputing only a subset of samples (e.g. QCs,
                         blanks, controls) while leaving others untouched.
    :return: DataFrame with missing values in the target columns replaced by zero.

    Example:

        result = imputation_zeros(data, on_cols=['featureA', 'featureB'])
        result = imputation_zeros(data, on_rows=['QC1', 'QC2', 'blank1'])
    """
    df, cols, rows = common.select_imputation_targets(data, on_cols, on_rows, drop_cols)

    if on_rows is not None:
        df.loc[rows, cols] = df.loc[rows, cols].fillna(0)
    else:
        df[cols] = df[cols].fillna(0)

    return df


def imputation_half_minimum(
    data: pd.DataFrame,
    drop_cols: Iterable[str] | None = None,
    on_cols: Iterable[str] | None = None,
    on_rows: Iterable[str] | None = None,
    group: str | None = None,
):
    """
    Replace missing values with half the per-column minimum of observed values.

    :param data: DataFrame with samples as rows and features as columns.
    :param list drop_cols: columns to permanently drop before imputation. If a column
                           appears in both ``on_cols`` and ``drop_cols`` it will be dropped
                           and a warning is emitted.
    :param list on_cols: columns to impute. If None, all numeric columns are used.
                         Non-numeric columns in ``on_cols`` will raise a TypeError.
    :param list on_rows: row index labels to restrict imputation to. If None, all rows are
                         imputed. When provided, the per-column minimum is also computed
                         from only those rows, so each subset gets its own half-minimum
                         (e.g. blanks are imputed with half the blank-minimum).
    :param str group: column label containing group identifiers. If set, the minimum is
                      computed within each group (and within ``on_rows``, if given).
                      Columns without observed values in a group are left as NaN for
                      that group. Rows with a missing group identifier are not imputed.
                      The group column itself is never imputed and is returned as the
                      first column.
    :return: DataFrame with missing values replaced by half the per-column minimum.

    Example::

        result = imputation_half_minimum(data, on_cols=['featureA', 'featureB'])
        result = imputation_half_minimum(data, on_rows=['blank1', 'blank2'])
        result = imputation_half_minimum(data, group='group')
    """
    df, cols, rows = common.select_imputation_targets(
        data, on_cols, on_rows, drop_cols, group=group
    )

    target = df.loc[rows, cols] if on_rows is not None else df[cols]
    where = "in the selected rows " if on_rows is not None else ""
    if group is None:
        all_nan = [c for c in cols if target[c].isna().all()]
        if all_nan:
            logger.warning(
                f"Columns with no observed values {where}cannot"
                f" be imputed and are left as NaN: {all_nan}"
            )
        half_mins = target.min(axis=0, skipna=True) / 2
    else:
        labels = df.loc[target.index, group] if on_rows is not None else df[group]
        grouped = target.groupby(labels)
        observed = grouped.min()
        missing = observed.isna().stack()
        all_nan = missing[missing].index.tolist()
        if all_nan:
            logger.warning(
                f"(group, column) pairs with no observed values"
                f" {where}cannot be imputed and are left as NaN: {all_nan}"
            )
        half_mins = grouped.transform("min") / 2

    if on_rows is not None:
        df.loc[rows, cols] = target.fillna(half_mins)
    else:
        df[cols] = target.fillna(half_mins)

    if group is not None:
        df = df[[group, *df.columns.drop(group)]]
    return df
