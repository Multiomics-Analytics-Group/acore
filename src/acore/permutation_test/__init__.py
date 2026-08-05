"""Permutation test functions for statistical analysis.

Nonparametric approaches are attractive because they make no assumptions about the
underlying distribution of the data — the significance is derived empirically from the
data itself.

Core idea: Permutation tests are a non-parametric method for hypothesis testing that
involves randomly shuffling the data to create a null distribution of the test statistic:
We test whether an observed difference/effect between groups is larger than what you'd
expect by chance, without assuming any particular data distribution.

Concretely, the logic is:

Observe a statistic (e.g., difference of means, t-statistic, chi-squared) computed on the
real, correctly-labeled data — this is your reference value. Assume the null hypothesis is
true: that group labels (or pairings) are arbitrary and don't matter — i.e., there's no
real difference between conditions. Shuffle/permute the data under that null assumption
many times (e.g., 10,000 times) — reassigning which values belong to which group/condition
— and recompute the same statistic each time. This produces a distribution of statistic
values that could arise purely by chance. Compare the observed statistic to this
permuted/null distribution. The p-value is the proportion of permuted statistics that are
as extreme or more extreme than the observed one using all the data and/or correct labels.
 If the observed effect is unusually
large compared to the shuffled versions, that's evidence it's not just random noise.

paired_permutation: For paired samples. Since pairing must be preserved, it doesn't
shuffle group membership directly — instead it randomly flips the sign of each paired
difference (cond1 - cond2), simulating the null hypothesis that within each pair, which
value is "cond1" vs "cond2" is arbitrary.

chi2_permutation: For categorical data. It shuffles group membership (via _permute) and
recomputes the chi-squared statistic on the resulting contingency table each time.
indep_permutation: For independent samples. It pools and reshuffles values between group1
and group2 (via _permute) since under the null there's no real distinction between the
groups. In all three, the p-value is computed as the fraction of permuted statistics whose
absolute value is ≥ ≥ the observed absolute statistic — i.e., p_value = mean(permuted >=
observed).
"""

import warnings
from collections.abc import Callable

import numpy as np
from scipy.stats import (
    chi2_contingency,
    f_oneway,
    ttest_ind,
    ttest_rel,
    wilcoxon,
)

from acore.types.permutation_test import PermutationResult

from .internal_functions import _check_degeneracy, _contingency_table, _permute

warnings.simplefilter("always", UserWarning)

STATS_PAIRED: list[str] = [ttest_rel.__name__, wilcoxon.__name__]


def paired_permutation(
    cond1: np.ndarray,
    cond2: np.ndarray,
    metric: str | Callable = "t-statistic",
    n_permutations: int = 10000,
    rng: np.random.Generator = None,
    **kwargs,
) -> dict:
    """
    Perform a permutation test for paired samples.

    Parameters
    ----------
    cond1 : np.ndarray
        First condition (paired samples).
    cond2 : np.ndarray
        Second condition (paired samples).
    metric : str or callable, optional
        Metric to compute ('t-statistic', 'mean', 'median', or
        a custom function that takes `cond1-cond2` as input).
    n_permutations : int, optional
        Number of permutations to perform (default is 10000).
    rng : np.random.Generator, optional
        Random number generator (default None triggers np.random.default_rng(seed=12345)).
    **kwargs
        Additional arguments passed to the metric function.

    Returns
    -------
    dict
        Dictionary with keys:
        - 'metric': Metric function used.
        - 'observed': Observed metric value.
        - 'p_value': Permutation test p-value (np.nan if degenerate).
    """
    # Validate input
    if cond1.shape != cond2.shape:
        raise ValueError("Input arrays must have the same shape.")

    if rng is None:
        rng = np.random.default_rng(seed=12345)

    # paired differences
    diff = cond1 - cond2
    args = [diff]

    # compute observed metric
    if metric == "t-statistic":
        calculator = ttest_rel
        args = [cond1, cond2]
    elif metric == "wilcoxon":
        calculator = wilcoxon
        args = [cond1, cond2]
    elif metric == "mean":
        calculator = np.mean
    elif metric == "median":
        calculator = np.median
    elif callable(metric):
        calculator = metric
        if calculator.__name__ in STATS_PAIRED:
            args = [cond1, cond2]
    else:
        raise ValueError(
            "Invalid metric specified. Acceptable metrics are: "
            "'t-statistic', 'mean', 'median', or a custom function "
            "that takes `cond1-cond2` as input."
        )
    # ? t-test observed is the t-statistic incl. a single p-value.
    # ? should this be save or rather the metric itself.
    observed_metric = calculator(*args, **kwargs)
    if calculator.__name__ in STATS_PAIRED:
        observed_statistic = observed_metric.statistic
    else:
        observed_statistic = observed_metric
    abs_met = abs(observed_statistic)

    # Perform permutations
    permuted_f = []
    for _ in range(n_permutations):
        # randomly flip direction of differences (sign)
        permuted_diff = diff * rng.choice([-1, 1], size=diff.shape)
        new_cond1 = np.where((permuted_diff == diff), cond1, cond2)
        new_cond2 = np.where((permuted_diff == diff), cond2, cond1)
        # postcondition check
        if not (permuted_diff == (new_cond1 - new_cond2)).all():
            raise ArithmeticError(
                "Postcondition failed: Issue with permuted differences"
            )
        # prep args and compute permuted metric
        if calculator.__name__ in STATS_PAIRED:
            new_result = abs(calculator(*[new_cond1, new_cond2], **kwargs).statistic)
        else:
            new_result = abs(calculator(*[permuted_diff], **kwargs))

        permuted_f.append(new_result)

    if _check_degeneracy(diff):
        identical_warn = (
            "Degenerate conditions detected (the data are identical). "
            "Consider using a different statistical test. "
            "Results may be unreliable."
        )
        warnings.warn(identical_warn)

        val_result = PermutationResult.model_validate(
            {
                "metric": metric,
                "observed_statistic": float(observed_statistic),
                "p_value": np.nan,
            }
        )

    else:
        # Compute p-value
        p_value = np.mean(np.asarray(permuted_f) >= abs_met)

        val_result = PermutationResult.model_validate(
            {
                "metric": metric,
                "observed_statistic": float(observed_statistic),
                "p_value": float(p_value),
            }
        )

    return val_result.model_dump()


def chi2_permutation(
    *groups,
    n_permutations: int = 10000,
    rng: np.random.Generator = None,
) -> dict:
    """
    Perform a permutation test for categorical data using the chi-squared statistic.

    Parameters
    ----------
    *groups : array-like
        Arrays representing categorical groups.
    n_permutations : int, optional
        Number of permutations to perform (default is 10000).
    rng : np.random.Generator, optional
        Random number generator (default None triggers np.random.default_rng(seed=12345)).

    Returns
    -------
    dict
        Dictionary with keys:
        - 'observed_statistic': Observed chi-squared test result.
        - 'p_value': Permutation test p-value.
    """

    if rng is None:
        rng = np.random.default_rng(seed=12345)

    # generate contingency table
    cont_table = _contingency_table(*groups, to_np=True)

    # Compute observed chi-squared statistic
    observed_test = chi2_contingency(cont_table)
    observed_chi2 = observed_test.statistic

    # Perform permutations
    permuted_chi2 = []
    for _ in range(n_permutations):
        # shuffle
        perm_cont_table = _contingency_table(*_permute(*groups, rng=rng))

        # calculate permuted chi-squared statistic
        permuted_chi2.append(chi2_contingency(perm_cont_table).statistic)

    # Compute p-value
    p_value = np.mean(permuted_chi2 >= observed_chi2)

    val_result = PermutationResult.model_validate(
        {
            "metric": "Chi-squared",
            "observed_statistic": float(observed_chi2),
            "p_value": float(p_value),
        }
    )

    return val_result.model_dump(exclude_none=True)


def indep_permutation(
    group1: np.ndarray,
    group2: np.ndarray,
    metric: str | Callable = "t-statistic",
    n_permutations: int = 10000,
    rng: np.random.Generator = None,
    **kwargs,
) -> dict:
    """
    Perform a permutation test for independent samples.

    Parameters
    ----------
    group1 : np.ndarray
        First group of samples.
    group2 : np.ndarray
        Second group of samples.
    metric : str or callable, optional
        Metric to compute ('t-statistic', 'anova', 'mean', 'median',
        or a function that would take groups 1 and 2
        as positional arguments 1 and 2 such as those in scipy.stats).
    n_permutations : int, optional
        Number of permutations to perform (default is 10000).
    rng : np.random.Generator, optional
        Random number generator (default None triggers np.random.default_rng(seed=12345)).
    **kwargs
        Additional arguments passed to the metric function.

    Returns
    -------
    dict
        Dictionary with keys:
        - 'metric': Metric function used.
        - 'observed_statistic': Observed metric value.
        - 'p_value': Permutation test p-value.
    """

    if rng is None:
        rng = np.random.default_rng(seed=12345)

    stat = False
    # what metric to use
    if metric == "t-statistic":
        calculator = ttest_ind
        stat = True
    elif metric == "anova":
        calculator = f_oneway
        stat = True
    elif metric == "mean":
        calculator = np.mean
    elif metric == "median":
        calculator = np.median
    elif callable(metric):
        calculator = metric
        stat = True
    else:
        raise ValueError(
            "Invalid metric specified. Acceptable metrics are: "
            "'t-statistic', 'mean', 'median', or a custom function "
            "that takes each group as input."
        )

    # compute observed metric
    if stat:
        observed_statistic = calculator(group1, group2, **kwargs).statistic
        # observed_statistic = observed_statistic.statistic
    else:
        observed_statistic = calculator(group1) - calculator(group2)
    observed_statistic = abs(observed_statistic)

    # Perform permutations
    permuted_f = []
    for _ in range(n_permutations):
        new_group1, new_group2 = _permute(group1, group2, rng=rng)
        # prep args
        if stat:
            new_result = calculator(new_group1, new_group2, **kwargs).statistic
        else:
            new_result = calculator(new_group1) - calculator(new_group2)
        abs_met = abs(new_result)
        # compute permuted metric
        permuted_f.append(abs_met)

    # Compute p-value
    p_value = np.mean(permuted_f >= observed_statistic)

    val_result = PermutationResult.model_validate(
        {
            "metric": metric,
            "observed_statistic": float(observed_statistic),
            "p_value": float(p_value),
        }
    )

    return val_result.model_dump()
