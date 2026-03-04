import logging
from typing import Iterable, Optional

import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer

DROP_COLS_DEFAULT = ("group", "sample", "subject")

logger = logging.getLogger(__name__)


def imputation_KNN(
    data: pd.DataFrame,
    drop_cols: Iterable[str] = DROP_COLS_DEFAULT,
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
    :param float cutoff: minimum ratio of missing/valid values required to impute
                         in each column.
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

            miss_df = miss_df.loc[:, miss_df.notna().mean() >= cutoff].dropna(
                axis="columns", how="all"
            )

            if miss_df.isna().any(axis=None):
                X = miss_df.to_numpy(dtype=np.float64, copy=False)
                X_trans = KNNImputer(n_neighbors=n_neighbors).fit_transform(X)

                dfm = pd.DataFrame(
                    X_trans, index=miss_df.index, columns=miss_df.columns
                )
                df_numeric.update(dfm)

        return df_numeric
    else:
        # Fallback: no grouping column present -> impute across all rows
        miss_df = df_numeric.loc[:, value_cols]
        miss_df = miss_df.loc[:, miss_df.notna().mean() >= cutoff].dropna(
            axis="columns", how="all"
        )
        if miss_df.isna().to_numpy().any():
            X = miss_df.to_numpy(dtype=np.float64, copy=False)
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
    data,
    index_cols=["group", "sample", "subject"],
    shift=1.8,
    nstd=0.3,
    group="group",
    cutoff=0.6,
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
    :param list index_cols: list of column labels to be set as dataframe index.
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
                    index_cols=['group', 'sample', 'subject'],
                    shift = 1.8, nstd = 0.3, group='group', cutoff=0.6
        )
    """
    df = imputation_KNN(
        data, drop_cols=index_cols, group=group, cutoff=cutoff, alone=False
    )
    df = imputation_normal_distribution(
        df, index_cols=index_cols, shift=shift, nstd=nstd
    )

    return df


def imputation_normal_distribution(
    data, index_cols=["group", "sample", "subject"], shift=1.8, nstd=0.3
):
    """
    Missing values will be replaced by random numbers that are drawn from a normal
    distribution. The imputation is done for each sample (across all proteins)
    separately.
    For more information visit `replacemissingfromgaussian`_ in coxdocs from MaxQuant.

    .. _replacemissingfromgaussian: https://cox-labs.github.io/coxdocs/\
replacemissingfromgaussian.html

    :param data: pandas dataframe with samples as rows and protein identifiers as
                 columns (with additional columns 'group', 'sample' and 'subject').
    :param list index_cols: list of column labels to be set as dataframe index.
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
                    index_cols=['group', 'sample', 'subject'],
                    shift = 1.8, nstd = 0.3
        )
    """
    np.random.seed(112736)
    df = data.copy()
    if index_cols is not None:
        df = df.set_index(index_cols)

    data_imputed = df.T.sort_index()
    null_columns = data_imputed.isnull().any().index.tolist()
    for c in null_columns:
        missing = data_imputed[data_imputed[c].isnull()].index.tolist()
        std = data_imputed[c].std()
        mean = data_imputed[c].mean()
        sigma = std * nstd
        mu = mean - (std * shift)
        value = 0.0
        if (
            not math.isnan(std)
            and not math.isnan(mean)
            and not math.isnan(sigma)
            and not math.isnan(mu)
        ):
            value = np.random.normal(mu, sigma, size=len(missing))
        i = 0
        for m in missing:
            if not isinstance(value, np.ndarray):
                data_imputed.loc[m, c] = value
            else:
                data_imputed.loc[m, c] = value[i]
                i += 1

    return data_imputed.T
