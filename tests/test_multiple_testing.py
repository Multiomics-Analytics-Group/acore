import numpy as np
import pytest

from acore.multiple_testing import apply_pvalue_correction

# benferoni correction
array = [
    0.1,
    0.01,
    0.2,
    0.4,
    np.nan,
    0.5,
]  # 0.04, 0.02, 0.0001]

test_res = [
    (
        "b",
        np.array([0.0, 1.0, 0.0, 0.0, np.nan, 0.0]),
        np.array([0.5, 0.05, 1.0, 1.0, np.nan, 1.0]),
    ),
    (
        "fdr_bh",
        np.array(
            [
                0.0,
                1.0,
                0.0,
                0.0,
                np.nan,
                0.0,
            ]
        ),
        np.array([0.25, 0.05, 0.33333333, 0.5, np.nan, 0.5]),
    ),
]


@pytest.mark.parametrize("method,exp_rejected,exp_pvalues", test_res)
def test_apply_pvalue_correction_alpha_5_percent(method, exp_rejected, exp_pvalues):
    act_rejected, act_pvalues = apply_pvalue_correction(
        array, alpha=0.05, method=method
    )
    np.testing.assert_array_equal(act_rejected, exp_rejected)
    np.testing.assert_array_almost_equal(act_pvalues, exp_pvalues)
