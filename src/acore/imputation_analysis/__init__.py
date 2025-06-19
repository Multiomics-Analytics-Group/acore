import math

import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer


def imputation_KNN(
    data,
    drop_cols=["group", "sample", "subject"],
    group="group",
    cutoff=0.6,
    alone=True,
):
    """
    k-Nearest Neighbors imputation for pandas dataframes with missing data. For more information visit https://github.com/iskandr/fancyimpute/blob/master/fancyimpute/knn.py.

    :param data: pandas dataframe with samples as rows and protein identifiers as columns (with additional columns 'group', 'sample' and 'subject').
    :param str group: column label containing group identifiers.
    :param list drop_cols: column labels to be dropped. Final dataframe should only have gene/protein/etc identifiers as columns.
    :param float cutoff: minimum ratio of missing/valid values required to impute in each column.
    :param bool alone: if True removes all columns with any missing values.
    :return: Pandas dataframe with samples as rows and protein identifiers as columns.

    Example::

        result = imputation_KNN(data, drop_cols=['group', 'sample', 'subject'], group='group', cutoff=0.6, alone=True)
    """
    np.random.seed(112736)
    df = data.copy()
    cols = df.columns
    df = df._get_numeric_data()
    if group in data.columns:
        df[group] = data[group]
        cols = list(set(cols).difference(df.columns))
        value_cols = [c for c in df.columns if c not in drop_cols]
        for g in df[group].unique():
            missDf = df.loc[df[group] == g, value_cols]
            missDf = missDf.loc[:, missDf.notnull().mean() >= cutoff].dropna(
                axis=1, how="all"
            )
            if missDf.isnull().values.any():
                X = np.array(missDf.values, dtype=np.float64)
                X_trans = KNNImputer(n_neighbors=3).fit_transform(X)
                missingdata_df = missDf.columns.tolist()
                dfm = pd.DataFrame(
                    X_trans, index=list(missDf.index), columns=missingdata_df
                )
                df.update(dfm)
        if alone:
            df = df.dropna(axis=1)

        df = df.join(data[cols])

    return df


def imputation_mixed_norm_KNN(
    data,
    index_cols=["group", "sample", "subject"],
    shift=1.8,
    nstd=0.3,
    group="group",
    cutoff=0.6,
):
    """
    Missing values are replaced in two steps: 1) using k-Nearest Neighbors we impute protein columns with a higher ratio of missing/valid values than the defined cutoff, \
    2) the remaining missing values are replaced by random numbers that are drawn from a normal distribution.

    :param data: pandas dataframe with samples as rows and protein identifiers as columns (with additional columns 'group', 'sample' and 'subject').
    :param str group: column label containing group identifiers.
    :param list index_cols: list of column labels to be set as dataframe index.
    :param float shift: specifies the amount by which the distribution used for the random numbers is shifted downwards. This is in units of the \
                        standard deviation of the valid data.
    :param float nstd: defines the width of the Gaussian distribution relative to the standard deviation of measured values. \
                        A value of 0.5 would mean that the width of the distribution used for drawing random numbers is half of the standard deviation of the data.
    :param float cutoff: minimum ratio of missing/valid values required to impute in each column.
    :return: Pandas dataframe with samples as rows and protein identifiers as columns.

    Example::

        result = imputation_mixed_norm_KNN(data, index_cols=['group', 'sample', 'subject'], shift = 1.8, nstd = 0.3, group='group', cutoff=0.6)
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
    Missing values will be replaced by random numbers that are drawn from a normal distribution. The imputation is done for each sample (across all proteins) separately.
    For more information visit http://www.coxdocs.org/doku.php?id=perseus:user:activities:matrixprocessing:imputation:replacemissingfromgaussian.

    :param data: pandas dataframe with samples as rows and protein identifiers as columns (with additional columns 'group', 'sample' and 'subject').
    :param list index_cols: list of column labels to be set as dataframe index.
    :param float shift: specifies the amount by which the distribution used for the random numbers is shifted downwards. This is in units of the standard deviation of the valid data.
    :param float nstd: defines the width of the Gaussian distribution relative to the standard deviation of measured values. A value of 0.5 would mean that the width of the distribution used for drawing random numbers is half of the standard deviation of the data.
    :return: Pandas dataframe with samples as rows and protein identifiers as columns.

    Example::

        result = imputation_normal_distribution(data, index_cols=['group', 'sample', 'subject'], shift = 1.8, nstd = 0.3)
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
