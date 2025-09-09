import numpy as np
import pandas as pd
from scipy.stats import (
    chi2_contingency, 
    ttest_rel,
    ttest_ind,
    f_oneway,
)

from .internal_functions import _permute, _contingency_table


def paired_permutation(
    cond1: np.ndarray,
    cond2: np.ndarray,
    metric: str = "t-statistic",
    n_permutations: int = 10000,
    rng: np.random.Generator = np.random.default_rng(seed=12345),
    **kwargs
) -> dict:
    """Perform a permutation test for paired samples."""
    # Validate input
    if not isinstance(cond1, np.ndarray) or not isinstance(cond2, np.ndarray):
        raise TypeError("Input must be numpy arrays.")
    if cond1.shape != cond2.shape:
        raise ValueError("Input arrays must have the same shape.")

    # paired differences
    diff = cond1-cond2
    args = [diff]

    # compute observed metric
    if metric == "t-statistic":
        calculator = ttest_rel
        args = [cond1, cond2]
    elif metric == "mean":
        calculator = np.mean
    elif metric == "median":
        calculator = np.median
    elif callable(metric):
        calculator = metric
    else:
        raise ValueError(
            "Invalid metric specified. Acceptable metrics are: "
            "'t-statistic', 'mean', 'median', or a custom function"
            "that takes `cond1-cond2` as input."
        )

    observed_metric = calculator(*args, **kwargs)
    if metric == "t-statistic":
        abs_met = abs(observed_metric.statistic)
    else:
        abs_met = abs(observed_metric)

    # Perform permutations
    permuted_f = []
    for _ in range(n_permutations):
        # randomly flip direction of differences
        permuted_diff = diff * rng.choice([-1, 1], size=diff.shape)
        new_cond1 = np.where((permuted_diff==diff), cond1, cond2)
        new_cond2 = np.where((permuted_diff==diff), cond2, cond1)
        # postcondition check
        if not (permuted_diff == (new_cond1 - new_cond2)).all():
            raise ArithmeticError("Postcondition failed: Issue with permuted differences")
        # prep args
        if metric == "t-statistic":
            new_args = [new_cond1, new_cond2]
            new_result = abs(calculator(*new_args, **kwargs).statistic)
        else:
            new_args = [permuted_diff]
            new_result = abs(calculator(*new_args, **kwargs))
        # compute permuted metric
        permuted_f.append(new_result)

    # Compute p-value
    p_value = np.mean(permuted_f >= abs_met)

    return {
        "metric": calculator,
        "observed": observed_metric,
        "p_value": p_value
    }


def chi2_permutation(
    *groups,
    n_permutations: int = 10000,
    rng: np.random.Generator = np.random.default_rng(seed=12345)
) -> dict:
    """Perform a permutation test for chi-squared statistics."""

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

    return {
        "observed": observed_test,
        "p_value": p_value
    }


def indep_permutation(
    group1: np.ndarray,
    group2: np.ndarray,
    metric: str = "t-statistic",
    n_permutations: int = 10000,
    rng: np.random.Generator = np.random.default_rng(seed=12345),
    **kwargs
) -> dict:
    """
    TODO
    """

    ## PRECONDITIONS
    if not isinstance(group1, np.ndarray) or not isinstance(group2, np.ndarray):
        raise TypeError("Input must be numpy arrays.")

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
    else:
        raise ValueError(
            "Invalid metric specified. Acceptable metrics are: "
            "'t-statistic', 'mean', 'median', or a custom function"
            "that takes each group as input."
        )

    # compute observed metric
    if stat: 
        observed_metric = calculator(group1, group2, **kwargs)
    else:
        observed_metric = abs(calculator(group1) - calculator(group2))

    # Perform permutations
    permuted_f = []
    for _ in range(n_permutations):
        new_group1, new_group2 = _permute(group1, group2, rng=rng)
        # prep args
        if stat:
            new_result = abs(calculator(new_group1, new_group2, **kwargs).statistic)
            abs_met = abs(observed_metric.statistic)
        else:
            new_result = abs(calculator(new_group1) - calculator(new_group2))
            abs_met = abs(observed_metric)
        # compute permuted metric
        permuted_f.append(new_result)

    # Compute p-value
    p_value = np.mean(permuted_f >= abs_met)

    return {
        "metric": calculator,
        "observed": observed_metric,
        "p_value": p_value
    }




