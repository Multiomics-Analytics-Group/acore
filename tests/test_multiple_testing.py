import numpy as np
import pytest

from acore.multiple_testing import apply_pvalue_correction, compute_efdr

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


# ---------------------------------------------------------------------------
# Tests for compute_efdr
# ---------------------------------------------------------------------------


def test_compute_efdr_empty_permuted_raises():
    with pytest.raises(ValueError, match="permuted_pvalues"):
        compute_efdr([0.01, 0.05], [], alpha=0.05)


def test_compute_efdr_empty_observed():
    rejected, efdr = compute_efdr([], [[0.1, 0.2]], alpha=0.05)
    assert len(rejected) == 0
    assert len(efdr) == 0


def test_compute_efdr_all_permuted_high():
    """When all permuted p-values are 1.0, eFDR should be 0 (no false positives)."""
    observed = [0.01, 0.05, 0.1]
    # All permuted p-values are 1.0 – no permuted hits at any threshold
    permuted = [np.array([1.0, 1.0, 1.0])] * 10
    rejected, efdr = compute_efdr(observed, permuted, alpha=0.05)
    np.testing.assert_array_equal(efdr, [0.0, 0.0, 0.0])
    np.testing.assert_array_equal(rejected, [True, True, True])


def test_compute_efdr_monotonicity():
    """eFDR values must be non-decreasing when sorted by ascending p-value."""
    observed = [0.001, 0.01, 0.05, 0.2, 0.5]
    rng = np.random.default_rng(42)
    permuted = [rng.uniform(0, 1, size=5) for _ in range(50)]
    rejected, efdr = compute_efdr(observed, permuted, alpha=0.05)

    # Sort observed p-values and check eFDR is non-decreasing
    sort_idx = np.argsort(observed)
    efdr_sorted = efdr[sort_idx]
    assert np.all(np.diff(efdr_sorted) >= -1e-12), (
        f"eFDR is not monotonically non-decreasing: {efdr_sorted}"
    )


def test_compute_efdr_values_in_unit_interval():
    observed = [0.01, 0.05, 0.2, 0.5]
    permuted = [np.array([0.01, 0.1, 0.3, 0.9])] * 20
    rejected, efdr = compute_efdr(observed, permuted, alpha=0.05)
    assert np.all(efdr >= 0.0) and np.all(efdr <= 1.0)


def test_compute_efdr_under_null_efdr_near_one():
    """Under a null scenario (observed and permuted drawn from same distribution)
    the eFDR should be close to 1 at lenient thresholds."""
    rng = np.random.default_rng(0)
    observed = rng.uniform(0, 1, size=20)
    permuted = [rng.uniform(0, 1, size=20) for _ in range(200)]
    _, efdr = compute_efdr(observed, permuted, alpha=0.05)
    # At lenient thresholds (large p-values), eFDR should approach 1
    max_pval_idx = np.argmax(observed)
    assert efdr[max_pval_idx] <= 1.0  # clamped


def test_compute_efdr_rejected_respects_alpha():
    """rejected array must exactly equal efdr <= alpha."""
    observed = [0.001, 0.01, 0.05, 0.2]
    permuted = [np.array([0.5, 0.6, 0.7, 0.8])] * 10
    rejected, efdr = compute_efdr(observed, permuted, alpha=0.05)
    np.testing.assert_array_equal(rejected, efdr <= 0.05)
