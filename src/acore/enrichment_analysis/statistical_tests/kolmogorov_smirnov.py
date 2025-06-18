from __future__ import annotations

from scipy import stats


def run_kolmogorov_smirnov(
    dist1: list[float], dist2: list[float], alternative: str = "two-sided"
) -> tuple[float, float]:
    """
    Compute the Kolmogorov-Smirnov statistic on 2 samples.
    See `scipy.stats.ks_2samp`_

    .. _scipy.stats.ks_2samp: https://docs.scipy.org/doc/scipy/reference/generated/\
scipy.stats.ks_2samp.html


    :param list dist1: sequence of 1-D ndarray (first distribution to compare)
        drawn from a continuous distribution
    :param list dist2: sequence of 1-D ndarray (second distribution to compare)
        drawn from a continuous distribution
    :param str alternative: defines the alternative hypothesis (default is ‘two-sided’):
        * **'two-sided'**
        * **'less'**
        * **'greater'**
    :return: statistic float and KS statistic pvalue float Two-tailed p-value.

    Example::

        result = run_kolmogorov_smirnov(dist1, dist2, alternative='two-sided')

    """

    result = stats.ks_2samp(dist1, dist2, alternative=alternative, mode="auto")

    return result
