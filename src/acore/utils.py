from collections.abc import Iterable

import numpy as np
import pandas as pd


def convertToEdgeList(data, cols):
    """
    This function converts a pandas dataframe to an edge list where index becomes the source nodes and columns the target nodes.

    :param data: pandas dataframe.
    :param list cols: names for dataframe columns.
    :return: Pandas dataframe with columns cols.
    """
    data.index.name = None
    edge_list = data.stack().reset_index()
    edge_list.columns = cols

    return edge_list


# TODO move to dsp_pandas
def check_is_paired(df, subject, group):
    """
    Check if samples are paired.

    :param df: pandas dataframe with samples as rows and protein identifiers as columns (with additional columns 'group', 'sample' and 'subject').
    :param str subject: column with subject identifiers
    :param str group: column with group identifiers
    :return: True if paired samples.
    :rtype: bool
    """
    is_pair = False
    if subject is not None:
        count_subject_groups = df.groupby(subject)[group].count()
        is_pair = (count_subject_groups > 1).all()

    return is_pair


def transform_into_wide_format(data, index, columns, values, extra=[]):
    """
    This function converts a Pandas DataFrame from long to wide format using
    pandas pivot_table() function.

    :param data: long-format Pandas DataFrame
    :param list index: columns that will be converted into the index
    :param str columns: column name whose unique values will become the new column names
    :param str values: column to aggregate
    :param list extra: additional columns to be kept as columns
    :return: Wide-format pandas DataFrame

    Example::

        result = transform_into_wide_format(df, index='index', columns='x', values='y', extra='group')

    """
    df = pd.DataFrame()
    extra_cols = None
    if data is not None:
        df = data.copy()
        if not df.empty:
            if index is None:
                df = df.reset_index()
                index = "index"
            cols = []
            append_to_list(cols, columns)
            append_to_list(cols, values)
            append_to_list(cols, index)
            if len(extra) > 0:
                append_to_list(extra, index)
                extra_cols = df[extra].set_index(index)
            df = df[cols]
            df = df.drop_duplicates()
            if isinstance(index, list):
                df = df.pivot_table(
                    index=index, columns=columns, values=values, aggfunc="first"
                )
            else:
                df = df.pivot(index=index, columns=columns, values=values)
            if extra_cols is not None:
                df = df.join(extra_cols)

            df = df.drop_duplicates()
            df = df.reset_index()

    return df


def transform_into_long_format(data, drop_columns, group, columns=["name", "y"]):
    """
    Converts a Pandas DataDrame from wide to long format using pd.melt()
    function.

    :param data: wide-format Pandas DataFrame
    :param list drop_columns: columns to be deleted
    :param group: column(s) to use as identifier variables
    :type group: str or list
    :param list columns: names to use for the 1)variable column, and for the 2)value column
    :return: Long-format Pandas DataFrame.

    Example::

        result = transform_into_long_format(df, drop_columns=['sample', 'subject'], group='group', columns=['name','y'])
    """
    long_data = pd.DataFrame()
    if data is not None:
        data = data.drop(drop_columns, axis=1)

        long_data = pd.melt(
            data, id_vars=group, var_name=columns[0], value_name=columns[1]
        )
        long_data = long_data.set_index(group)
        long_data.columns = columns

    return long_data


def remove_group(data):
    """
    Removes column with label 'group'.

    :param data: pandas dataframe with one column labelled 'group'
    :return: Pandas dataframe

    Example::

        result = remove_group(data)
    """
    data.drop(["group"], axis=1)
    return data


def calculate_fold_change(df, condition1, condition2):
    """
    Calculates fold-changes between two groups for all proteins in a dataframe.

    :param df: pandas dataframe with samples as rows and protein identifiers as columns.
    :param str condition1: identifier of first group.
    :param str condition2: identifier of second group.
    :return: Numpy array.

    Example::

        result = calculate_fold_change(data, 'group1', 'group2')
    """
    group1 = df[condition1]
    group2 = df[condition2]

    if isinstance(group1, np.float64):
        group1 = np.array(group1)
    else:
        group1 = group1.values
    if isinstance(group2, np.float64):
        group2 = np.array(group2)
    else:
        group2 = group2.values

    if np.isnan(group1).all() or np.isnan(group2).all():
        fold_change = np.nan
    else:
        fold_change = np.nanmedian(group1) - np.nanmedian(group2)

    return fold_change


def pooled_standard_deviation(sample1, sample2, ddof):
    """
    Calculates the pooled standard deviation.
    For more information visit https://www.hackdeploy.com/learn-what-is-statistical-power-with-python/.

    :param array sample1: numpy array with values for first group
    :param array sample2: numpy array with values for second group
    :param int ddof: degrees of freedom
    """
    # calculate the sample size
    n1, n2 = len(sample1), len(sample2)
    # calculate the variances
    var1, var2 = np.var(sample1, ddof=1), np.var(sample2, ddof=ddof)
    # calculate the pooled standard deviation
    numerator = ((n1 - 1) * var1) + ((n2 - 1) * var2)
    denominator = n1 + n2 - 2
    return np.sqrt(numerator / denominator)


def cohens_d(sample1, sample2, ddof):
    """
    Calculates Cohen's d effect size based on the distance between two means, measured in standard deviations.
    For more information visit https://www.hackdeploy.com/learn-what-is-statistical-power-with-python/.

    :param array sample1: numpy array with values for first group
    :param array sample2: numpy array with values for second group
    :param int ddof: degrees of freedom
    """
    u1, u2 = np.mean(sample1), np.mean(sample2)
    s_pooled = pooled_standard_deviation(sample1, sample2, ddof)

    return (u1 - u2) / s_pooled


def hedges_g(df, condition1, condition2, ddof=0):
    """
    Calculates Hedges’ g effect size (more accurate for sample sizes below 20 than Cohen's d).
    For more information visit https://docs.scipy.org/doc/numpy/reference/generated/numpy.nanstd.html.

    :param df: pandas dataframe with samples as rows and protein identifiers as columns.
    :param str condition1: identifier of first group.
    :param str condition2: identifier of second group.
    :param int ddof: means Delta Degrees of Freedom.
    :return: Numpy array.

    Example::

        result = hedges_g(data, 'group1', 'group2', ddof=0)
    """
    group1 = df[condition1]
    group2 = df[condition2]

    if isinstance(group1, np.float64):
        group1 = np.array(group1)
    else:
        group1 = group1.values
    if isinstance(group2, np.float64):
        group2 = np.array(group2)
    else:
        group2 = group2.values

    ng1 = group1.size
    ng2 = group2.size
    # dof = ng1 + ng2 - 2
    if np.isnan(group1).all() or np.isnan(group2).all():
        g = np.nan
    else:
        meang1 = np.nanmean(group1)
        meang2 = np.nanmean(group2)
        sdpooled = np.nanstd(np.concatenate([group1, group2]), ddof=ddof)

        # Correct bias small sample size
        if ng1 + ng2 < 50:
            g = (
                ((meang1 - meang2) / sdpooled)
                * ((ng1 + ng2 - 3) / (ng1 + ng2 - 2.25))
                * np.sqrt((ng1 + ng2 - 2) / (ng1 + ng2))
            )
        else:
            g = (meang1 - meang2) / sdpooled

    return g


def unit_vector(vector):
    """
    Returns the unit vector of the vector.
    :param tuple vector: vector
    :return tuple unit_vector: unit vector
    """
    return vector / np.linalg.norm(vector)


def flatten(t, my_list=[]):
    """
    Code from: https://gist.github.com/shaxbee/0ada767debf9eefbdb6e
    Acknowledgements: Zbigniew Mandziejewicz (shaxbee)
    Generator flattening the structure

    >>> list(flatten([2, [2, (4, 5, [7], [2, [6, 2, 6, [6], 4]], 6)]]))
    [2, 2, 4, 5, 7, 2, 6, 2, 6, 6, 4, 6]
    """

    for x in t:
        if not isinstance(x, Iterable) or isinstance(x, str):
            my_list.append(x)
        else:
            flatten(x, my_list)

    return my_list


def angle_between(v1, v2):
    """
    Returns the angle in radians between vectors 'v1' and 'v2'

    :param tuple v1: vector 1
    :param tuple v2: vector 2
    :return float angle: angle between two vectors in radians

    Example::
        angle = angle_between((1, 0, 0), (0, 1, 0))
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def append_to_list(mylist, myappend):
    if isinstance(myappend, list):
        mylist.extend(myappend)
    else:
        mylist.append(myappend)


def generator_to_dict(genvar):
    dictvar = {}
    for i, gen in enumerate(genvar):
        dictvar.update({n: i for n in gen})

    return dictvar
