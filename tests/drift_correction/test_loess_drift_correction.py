import numpy as np
import pandas as pd
import pytest

from acore.drift_correction.loess_drift_correction import (
    filter_features_by_qc,
    qc_rlsc_loess,
    run_loess_drift_correction,
)

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

QC_ROWS = [f"qc{i}" for i in range(6)]
SAMPLE_ROWS = [f"s{i}" for i in range(4)]
ALL_ROWS = SAMPLE_ROWS + QC_ROWS


@pytest.fixture
def base_df():
    """Small intensity DataFrame with no missing values."""
    rng = np.random.default_rng(0)
    data = rng.uniform(10, 200, size=(10, 5))
    return pd.DataFrame(data, index=ALL_ROWS, columns=[f"f{i}" for i in range(5)])


@pytest.fixture
def sample_order():
    """Injection order mapping covering all rows in ALL_ROWS."""
    rows = ALL_ROWS
    return pd.DataFrame({"File Name": rows, "Sample ID": list(range(1, len(rows) + 1))})


# ---------------------------------------------------------------------------
# filter_features_by_qc
# ---------------------------------------------------------------------------


def test_filter_features_by_qc_keeps_complete_columns():
    df = pd.DataFrame(
        {"A": [1.0, 2.0, 3.0], "B": [4.0, 5.0, 6.0]},
        index=["q1", "q2", "q3"],
    )
    result = filter_features_by_qc(df, ["q1", "q2", "q3"], threshold=0.5)
    assert list(result.columns) == ["A", "B"]


def test_filter_features_by_qc_removes_sparse_column():
    df = pd.DataFrame(
        {"A": [1.0, np.nan, np.nan], "B": [4.0, 5.0, 6.0]},
        index=["q1", "q2", "q3"],
    )
    # threshold=0.5 requires ceil(3*0.5)=2 valid values; A has only 1
    result = filter_features_by_qc(df, ["q1", "q2", "q3"], threshold=0.5)
    assert "A" not in result.columns
    assert "B" in result.columns


def test_filter_features_by_qc_threshold_boundary():
    """A feature with exactly the minimum required valid values is kept."""
    df = pd.DataFrame(
        {"A": [1.0, np.nan, 3.0]},
        index=["q1", "q2", "q3"],
    )
    # ceil(3 * 0.5) = 2 required; A has 2 valid — should be kept
    result = filter_features_by_qc(df, ["q1", "q2", "q3"], threshold=0.5)
    assert "A" in result.columns


def test_filter_features_by_qc_invalid_threshold():
    df = pd.DataFrame({"A": [1.0, 2.0]}, index=["q1", "q2"])
    with pytest.raises(ValueError, match="Threshold"):
        filter_features_by_qc(df, ["q1", "q2"], threshold=1.5)


def test_filter_features_by_qc_non_qc_rows_ignored():
    """Only QC rows are counted; NaN in non-QC rows has no effect."""
    df = pd.DataFrame(
        {"A": [np.nan, 2.0, 3.0]},
        index=["sample", "q1", "q2"],
    )
    result = filter_features_by_qc(df, ["q1", "q2"], threshold=0.5)
    assert "A" in result.columns


# ---------------------------------------------------------------------------
# qc_rlsc_loess
# ---------------------------------------------------------------------------


def _linear_qc_data(n=8):
    """QC data with a known linear trend (easy to fit)."""
    x_qc = np.arange(1, n + 1, dtype=float)
    y_qc = 100.0 + 2.0 * x_qc
    x_all = np.arange(1, n + 3, dtype=float)  # a couple of extra samples
    return x_qc, y_qc, x_all


def test_qc_rlsc_loess_returns_correct_shapes():
    x_qc, y_qc, x_all = _linear_qc_data()
    drift_curve, best_alpha = qc_rlsc_loess(x_qc, y_qc, x_all)
    assert drift_curve.shape == x_all.shape
    assert isinstance(best_alpha, float)


def test_qc_rlsc_loess_no_negatives():
    x_qc, y_qc, x_all = _linear_qc_data()
    drift_curve, _ = qc_rlsc_loess(x_qc, y_qc, x_all)
    assert np.all(drift_curve >= 1e-6)


def test_qc_rlsc_loess_clamping_before_range():
    """Points before the first QC must equal the first LOESS value."""
    x_qc = np.array([5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
    y_qc = np.ones(6) * 100.0
    x_all = np.array([1.0, 2.0, 5.0, 10.0, 12.0])  # 1,2 before range; 12 after
    drift_curve, _ = qc_rlsc_loess(x_qc, y_qc, x_all, always_use_default=True)
    # Values before range should be clamped to the first fitted QC value
    assert drift_curve[0] == drift_curve[1]


def test_qc_rlsc_loess_clamping_after_range():
    """Points after the last QC must equal the last LOESS value."""
    x_qc = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    y_qc = np.ones(6) * 50.0
    x_all = np.array([1.0, 6.0, 8.0, 10.0])
    drift_curve, _ = qc_rlsc_loess(x_qc, y_qc, x_all, always_use_default=True)
    assert drift_curve[2] == drift_curve[3]


def test_qc_rlsc_loess_always_use_default_skips_loocv():
    """always_use_default=True should return the default alpha."""
    x_qc, y_qc, x_all = _linear_qc_data()
    _, best_alpha = qc_rlsc_loess(
        x_qc, y_qc, x_all, always_use_default=True, default=0.75
    )
    assert best_alpha == 0.75


def test_qc_rlsc_loess_alpha_in_candidates():
    """Selected alpha should come from the candidate list."""
    candidates = np.array([0.4, 0.6, 0.8, 1.0])
    x_qc, y_qc, x_all = _linear_qc_data()
    _, best_alpha = qc_rlsc_loess(x_qc, y_qc, x_all, alpha_candidates=candidates)
    assert best_alpha in candidates


# ---------------------------------------------------------------------------
# run_loess_drift_correction
# ---------------------------------------------------------------------------


def test_run_loess_output_shape(base_df, sample_order):
    result, info = run_loess_drift_correction(
        base_df, QC_ROWS, SAMPLE_ROWS, sample_order, always_use_default=True
    )
    assert result.shape == base_df.shape


def test_run_loess_index_columns_preserved(base_df, sample_order):
    result, _ = run_loess_drift_correction(
        base_df, QC_ROWS, SAMPLE_ROWS, sample_order, always_use_default=True
    )
    assert list(result.index) == list(base_df.index)
    assert list(result.columns) == list(base_df.columns)


def test_run_loess_values_changed(base_df, sample_order):
    result, _ = run_loess_drift_correction(
        base_df, QC_ROWS, SAMPLE_ROWS, sample_order, always_use_default=True
    )
    original = base_df.loc[ALL_ROWS].values
    corrected = result.loc[ALL_ROWS].values
    assert not np.allclose(original, corrected)


def test_run_loess_correction_info_keys(base_df, sample_order):
    _, info = run_loess_drift_correction(
        base_df, QC_ROWS, SAMPLE_ROWS, sample_order, always_use_default=True
    )
    for feature, entry in info.items():
        assert "status" in entry
        assert "alpha" in entry
        assert "rsd_qc" in entry


def test_run_loess_corrected_status(base_df, sample_order):
    _, info = run_loess_drift_correction(
        base_df, QC_ROWS, SAMPLE_ROWS, sample_order, always_use_default=True
    )
    for feature in base_df.columns:
        assert info[feature]["status"] == "corrected"


def test_run_loess_extra_rows_unchanged(base_df, sample_order):
    """Rows not listed in sample_rows or qc_rows pass through unmodified."""
    extra = pd.DataFrame([[999.0] * 5], index=["extra"], columns=base_df.columns)
    df_with_extra = pd.concat([base_df, extra])
    result, _ = run_loess_drift_correction(
        df_with_extra, QC_ROWS, SAMPLE_ROWS, sample_order, always_use_default=True
    )
    assert result.loc["extra"].tolist() == [999.0] * 5


def test_run_loess_skipped_when_too_few_qcs(sample_order):
    """Features with fewer valid QC values than qc_min_threshold are not corrected."""
    rng = np.random.default_rng(1)
    data = rng.uniform(10, 200, size=(10, 2))
    df = pd.DataFrame(data, index=ALL_ROWS, columns=["f0", "f1"])
    # Blank out most QC values so fewer than qc_min_threshold=4 remain
    df.loc[QC_ROWS[2:], "f0"] = np.nan

    _, info = run_loess_drift_correction(
        df,
        QC_ROWS,
        SAMPLE_ROWS,
        sample_order,
        qc_min_threshold=4,
        always_use_default=True,
    )
    assert info["f0"]["status"] == "skipped_due_to_few_qcs"
    assert info["f1"]["status"] == "corrected"


def test_run_loess_filter_percent_excludes_from_info(sample_order):
    """Features failing the QC completeness check are not corrected and absent from info."""
    rng = np.random.default_rng(2)
    data = rng.uniform(10, 200, size=(10, 3))
    df = pd.DataFrame(data, index=ALL_ROWS, columns=["f0", "f1", "f2"])
    # Make f0 have only 1 valid QC value (below 50% of 6 QCs)
    df.loc[QC_ROWS[1:], "f0"] = np.nan

    result, info = run_loess_drift_correction(
        df,
        QC_ROWS,
        SAMPLE_ROWS,
        sample_order,
        filter_percent=0.5,
        always_use_default=True,
    )
    # Output DataFrame always matches input shape; f0 retains original (uncorrected) values
    assert result.shape == df.shape
    # Filtered features are skipped entirely — not present in correction_info
    assert "f0" not in info
    assert "f1" in info and "f2" in info


def test_run_loess_invalid_sample_id_raises(base_df):
    """Non-numeric Sample ID values should raise a ValueError."""
    bad_order = pd.DataFrame(
        {"File Name": ALL_ROWS, "Sample ID": ["a"] * len(ALL_ROWS)}
    )
    with pytest.raises(ValueError, match="non-integer"):
        run_loess_drift_correction(
            base_df, QC_ROWS, SAMPLE_ROWS, bad_order, always_use_default=True
        )
