import numpy as np
import pandas as pd
import scipy.stats


def get_histogram_series(s: pd.Series, bins: np.ndarray) -> pd.Series:
    hist_data, bin_edges = np.histogram(s, bins=bins)
    ret = pd.Series(
        hist_data,
        index=pd.MultiIndex.from_arrays(
            [bin_edges[:-1], bin_edges[1:], range(len(hist_data))],
            names=["bin_number", "bin_start", "bin_end"],
        ),
    )
    return ret


def calculate_coefficient_variation(df: pd.DataFrame) -> pd.Series:
    """
    Compute the coefficient of variation (CV) for each column in a DataFrame.

    The coefficient of variation is defined as the ratio of the standard deviation
    to the mean, expressed as a percentage. This function uses the biased standard
    deviation (normalization by N) as implemented in `scipy.stats.variation`.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame containing numeric values (e.g., log2-transformed data).
        Each column will be processed independently.

    Returns
    -------
    pandas.Series
        Series containing the coefficient of variation (in percent) for each column.
        The index corresponds to the columns of the input DataFrame.

    See Also
    --------
    scipy.stats.variation : Function used to compute the coefficient of variation.

    Examples
    --------
    >>> import pandas as pd
    >>> import numpy as np
    >>> df = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})
    >>> calculate_coefficient_variation(df)
    A   40.825
    B   16.330
    Name: coef_of_var, dtype: float64
    """
    cv = scipy.stats.variation(df, axis=0) * 100
    cv = pd.Series(cv, index=df.columns).rename("coef_of_var")
    return cv


def calculate_coef_of_var_and_mean(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Calculate coefficient of variation and mean for each column in the dataframe,
    the mean calculated on both log2 and linear scale.

    Parameters
    ----------
    df : pd.DataFrame
        The input dataframe containing the linear values (non-log transformed).

    Returns
    -------
    pd.DataFrame
        A dataframe with columns 'mean_log2', 'mean' and 'coef_of_var' for each column
        in the input dataframe.
    """
    cv = calculate_coefficient_variation(df)
    means = df.mean().rename("mean")
    means_logs = np.log2(df).mean().rename("mean_log2")

    return pd.concat(
        [
            means_logs,
            means,
            cv,
        ],
        axis=1,
    ).rename_axis("name")


def get_coefficient_variation(
    data: pd.DataFrame,
    drop_columns: list[str] | None = None,
    group: str = "group",
):
    """
    Extracts the coefficients of variation in each group.

    :param data: pandas.DataFrame with samples as rows and protein identifiers as columns
                 (with additional columns 'group', 'sample' and 'subject'). The values
                 should be the original intensities for massspectrometry-based
                 measurements.
    :param list drop_columns: column labels to be dropped from the dataframe
    :param str group: column label containing group identifiers.
    :return: Pandas dataframe with columns 'name' (protein identifier),
             'x' (coefficient of variation), 'y' (mean) and 'group'.

    Example::

        result = get_coefficient_variation(data, drop_columns=['sample', 'subject'], group='group')
    """
    formated_df = data
    if drop_columns is not None:
        formated_df = data.drop(drop_columns, axis=1)
    cvs = formated_df.groupby(group).apply(func=calculate_coef_of_var_and_mean)
    cvs = cvs.reset_index()[["name", "mean_log2", "mean", "coef_of_var", group]]
    cvs_df = cvs

    return cvs_df
