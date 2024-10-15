"""Strategies for data normalization used by `normalize_data`."""

import pandas as pd
from sklearn import preprocessing


# ! update docstring: example is not correct
def median_zero_normalization(data, normalize="samples"):
    """
    This function normalizes each sample by using its median.

    :param data:
    :param str normalize: whether the normalization should be done by 'features' (columns)
                          or 'samples' (rows)
    :return: pandas.DataFrame.

    Example::
        data = pd.DataFrame({'a': [2,5,4,3,3], 'b':[4,4,6,5,3], 'c':[4,14,8,8,9]})
        result = median_normalization(data, normalize='samples')
        result
                a         b         c
            0 -1.333333  0.666667  0.666667
            1 -2.666667 -3.666667  6.333333
            2 -2.000000  0.000000  2.000000
            3 -2.333333 -0.333333  2.666667
            4 -2.000000 -2.000000  4.000000
    """
    if normalize is None or normalize == "samples":
        norm_data = data.sub(data.median(axis=1), axis=0)
    else:
        norm_data = data.sub(data.median(axis=0), axis=1)

    return norm_data


# ! Update docstring: example is not correct
def median_normalization(data, normalize="samples"):
    """
    This function normalizes each sample by using its median.

    :param data:
    :param str normalize: whether the normalization should be done by 'features' (columns)
                          or 'samples' (rows)
    :return: pandas.DataFrame.

    Example::
        data = pd.DataFrame({'a': [2,5,4,3,3], 'b':[4,4,6,5,3], 'c':[4,14,8,8,9]})
        result = median_normalization(data, normalize='samples')
        result
                a         b         c
            0 -1.333333  0.666667  0.666667
            1 -2.666667 -3.666667  6.333333
            2 -2.000000  0.000000  2.000000
            3 -2.333333 -0.333333  2.666667
            4 -2.000000 -2.000000  4.000000
    """
    if normalize is None or normalize == "samples":
        norm_data = data.sub(data.median(axis=1) - data.median(axis=1).median(), axis=0)
    else:
        norm_data = data.sub(data.median(axis=0) - data.median(axis=0).median(), axis=1)

    return norm_data


def zscore_normalization(data, normalize="samples"):
    """
    This function normalizes each sample by using its mean and standard deviation
    (mean=0, std=1).

    :param data:
    :param str normalize: whether the normalization should be done by 'features' (columns)
                          or 'samples' (rows)
    :return: pandas.DataFrame.

    Example::
        data = pd.DataFrame({'a': [2,5,4,3,3], 'b':[4,4,6,5,3], 'c':[4,14,8,8,9]})
        result = zscore_normalization(data, normalize='samples')
        result
                  a         b         c
                0 -1.154701  0.577350  0.577350
                1 -0.484182 -0.665750  1.149932
                2 -1.000000  0.000000  1.000000
                3 -0.927173 -0.132453  1.059626
                4 -0.577350 -0.577350  1.154701
    """
    if normalize is None or normalize == "samples":
        norm_data = data.sub(data.mean(axis=1), axis=0).div(data.std(axis=1), axis=0)

    else:
        norm_data = data.sub(data.mean(axis=0), axis=1).div(data.std(axis=0), axis=1)

    return norm_data


def median_polish_normalization(data, max_iter=250):
    """
    This function iteratively normalizes each sample and each feature to its
    median until medians converge.

    :param data:
    :param int max_iter: number of maximum iterations to prevent infinite loop.
    :return: pandas.DataFrame.

    Example::
        data = pd.DataFrame({'a': [2,5,4,3,3], 'b':[4,4,6,5,3], 'c':[4,14,8,8,9]})
        result = median_polish_normalization(data, max_iter = 10)
        result
                a    b     c
            0  2.0  4.0   7.0
            1  5.0  7.0  10.0
            2  4.0  6.0   9.0
            3  3.0  5.0   8.0
            4  3.0  5.0   8.0
    """
    mediandf = data.copy()
    for _ in range(max_iter):
        row_median = mediandf.median(axis=1)
        mediandf = mediandf.sub(row_median, axis=0)
        col_median = mediandf.median(axis=0)
        mediandf = mediandf.sub(col_median, axis=1)
        if (mediandf.median(axis=0) == 0).all() and (
            mediandf.median(axis=1) == 0
        ).all():
            break

    norm_data = data - mediandf

    return norm_data


def quantile_normalization(data) -> pd.DataFrame:
    """
    Applies quantile normalization to each column in pandas.DataFrame.

    :param data: pandas.DataFrame with features as columns and samples as rows.
    :return: pandas.DataFrame

    Example::
        data = pd.DataFrame({'a': [2,5,4,3,3], 'b':[4,4,6,5,3], 'c':[4,14,8,8,9]})
        result = quantile_normalization(data)
        result
                a    b    c
            0  3.2  4.6  4.6
            1  4.6  3.2  8.6
            2  3.2  4.6  8.6
            3  3.2  4.6  8.6
            4  3.2  3.2  8.6
    """
    rank_mean = (
        data.T.stack().groupby(data.T.rank(method="first").stack().astype(int)).mean()
    )
    normdf = data.T.rank(method="min").stack().astype(int).map(rank_mean).unstack().T

    return normdf


def linear_normalization(data, method="l1", normalize="samples") -> pd.DataFrame:
    """
    This function scales input data to a unit norm. For more information visit:
    https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.normalize.html

    :param data: pandas.DataFrame with samples as rows and features as columns.
    :param str method: norm to use to normalize each non-zero sample or non-zero feature
                        (depends on axis).
    :param str normalize: axis used to normalize the data along. If 'samples',
                independently normalize each sample, if 'features' normalize each feature.
    :return: pandas.DataFrame

    Example::
        data = pd.DataFrame({'a': [2,5,4,3,3], 'b':[4,4,6,5,3], 'c':[4,14,8,8,9]})
        result = linear_normalization(data, method = "l1", by = 'samples')
        result
                a         b         c
            0  0.117647  0.181818  0.093023
            1  0.294118  0.181818  0.325581
            2  0.235294  0.272727  0.186047
            3  0.176471  0.227273  0.186047
            4  0.176471  0.136364  0.209302
    """
    if normalize is None or normalize == "samples":
        normvalues = preprocessing.normalize(
            data.fillna(0).values, norm=method, axis=0, copy=True, return_norm=False
        )
    else:
        normvalues = preprocessing.normalize(
            data.fillna(0).values, norm=method, axis=1, copy=True, return_norm=False
        )

    normdf = pd.DataFrame(normvalues, index=data.index, columns=data.columns)

    return normdf
