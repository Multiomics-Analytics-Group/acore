from acore import permutation_test as pt
import pytest
import numpy as np
from scipy.stats import ttest_rel, ttest_ind, mannwhitneyu, wilcoxon


@pytest.fixture
def set_array():
    return np.arange(1, 10)


@pytest.fixture
def random_low():
    return np.random.negative_binomial(10, 0.3, 100)


@pytest.fixture
def random_high():
    return np.random.negative_binomial(20, 0.1, 100)


@pytest.fixture
def random_choice_equal():
    return np.random.choice(
        ["A", "B", "C", "D"], size=100, replace=True, p=[0.25, 0.25, 0.25, 0.25]
    )


@pytest.fixture
def random_choice_unequal():
    return np.random.choice(
        ["A", "B", "C", "D"], size=100, replace=True, p=[0.1, 0.25, 0.05, 0.6]
    )


def test_paired_permutation_same(set_array):
    cond1 = set_array
    cond2 = set_array
    result = pt.paired_permutation(cond1, cond2, metric=np.mean)
    assert result["metric"] == np.mean
    assert result["observed"] == 0.0
    assert result["p_value"] is np.nan


def test_paired_permutation_reject(random_low, random_high):
    cond1 = random_low
    cond2 = random_high
    for stat in ["t-statistic", "mean", "median"]:
        result = pt.paired_permutation(cond1, cond2, metric=stat)
        print(result)
        print(ttest_rel(cond1, cond2))
        assert result["p_value"] < 0.05


def test_paired_permutation_degen(set_array):
    for stat in ["t-statistic", "mean", "median"]:
        assert (
            pt.paired_permutation(set_array, set_array + 30, metric=stat)["p_value"]
            is np.nan
        )


def test_indep_permutation_reject(random_low, random_high):
    group1 = random_low
    group2 = random_high
    for stat in ["t-statistic", "anova", "median", "mean"]:
        result = pt.indep_permutation(group1, group2, metric=stat)
        assert result["p_value"] < 0.05


def test_indep_permutation_accept(random_low):
    for stat in ["t-statistic", "anova", "median", "mean"]:
        assert (
            pt.indep_permutation(random_low, random_low, metric=stat)["p_value"] > 0.05
        )


def test_chi2_permutation_accept(random_choice_equal):
    result = pt.chi2_permutation(random_choice_equal, random_choice_equal)
    assert result["p_value"] >= 0.05


def test_chi2_permutation_reject(random_choice_unequal, random_choice_equal):
    result = pt.chi2_permutation(random_choice_equal, random_choice_unequal)
    assert result["p_value"] < 0.05
