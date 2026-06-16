import numpy as np
import pandas as pd
import pytest

from acore.drift_correction.cpca_drift_correction import (
    check_missingness,
    cpca_centroid,
    run_cpca_drift_correction,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def base_df():
    """Small DataFrame (rows=samples, cols=features) with no missing values."""
    rng = np.random.default_rng(42)
    data = rng.uniform(1, 100, size=(8, 6))
    index = [f"s{i}" for i in range(4)] + [f"qc{i}" for i in range(4)]
    cols = [f"f{i}" for i in range(6)]
    return pd.DataFrame(data, index=index, columns=cols)


SAMPLE_ROWS = [f"s{i}" for i in range(4)]
QC_ROWS = [f"qc{i}" for i in range(4)]


# ---------------------------------------------------------------------------
# check_missingness
# ---------------------------------------------------------------------------


def test_check_missingness_no_na():
    df = pd.DataFrame({"a": [1.0, 2.0, 3.0], "b": [4.0, 5.0, 6.0]})
    assert check_missingness(df, [0, 1, 2]) is False


def test_check_missingness_with_na_in_checked_rows():
    df = pd.DataFrame({"a": [1.0, np.nan, 3.0], "b": [4.0, 5.0, 6.0]})
    assert check_missingness(df, [0, 1, 2]) is True


def test_check_missingness_na_only_outside_checked_rows():
    """NaN exists but not in the rows we check — should return False."""
    df = pd.DataFrame({"a": [1.0, np.nan, 3.0], "b": [4.0, 5.0, 6.0]})
    assert check_missingness(df, [0, 2]) is False


# ---------------------------------------------------------------------------
# run_cpca_drift_correction
# ---------------------------------------------------------------------------


def test_run_cpca_drift_correction_output_shape(base_df):
    result = run_cpca_drift_correction(base_df, SAMPLE_ROWS, QC_ROWS)
    assert result.shape == base_df.shape


def test_run_cpca_drift_correction_index_preserved(base_df):
    result = run_cpca_drift_correction(base_df, SAMPLE_ROWS, QC_ROWS)
    assert list(result.index) == list(base_df.index)
    assert list(result.columns) == list(base_df.columns)


def test_run_cpca_drift_correction_values_change(base_df):
    result = run_cpca_drift_correction(base_df, SAMPLE_ROWS, QC_ROWS)
    intensity_rows = SAMPLE_ROWS + QC_ROWS
    original = base_df.loc[intensity_rows].values
    corrected = result.loc[intensity_rows].values
    assert not np.allclose(original, corrected)


def test_run_cpca_drift_correction_n_comps(base_df):
    """n_comps=2 should still return a DataFrame with same shape."""
    result = run_cpca_drift_correction(base_df, SAMPLE_ROWS, QC_ROWS, n_comps=2)
    assert result.shape == base_df.shape


def test_run_cpca_drift_correction_nan_raises(base_df):
    df_nan = base_df.copy()
    df_nan.iloc[0, 0] = np.nan
    with pytest.raises(ValueError, match="NA values"):
        run_cpca_drift_correction(df_nan, SAMPLE_ROWS, QC_ROWS)


def test_run_cpca_drift_correction_extra_rows_unchanged(base_df):
    """Rows not in sample_rows or qc_rows must pass through untouched."""
    extra_row = pd.DataFrame([[999.0] * 6], index=["extra"], columns=base_df.columns)
    df_with_extra = pd.concat([base_df, extra_row])
    result = run_cpca_drift_correction(df_with_extra, SAMPLE_ROWS, QC_ROWS)
    assert result.loc["extra"].tolist() == [999.0] * 6


def test_run_cpca_drift_correction_n_comps_changes_values(base_df):
    """n_comps=1 and n_comps=2 must produce different corrected values."""
    intensity_rows = SAMPLE_ROWS + QC_ROWS
    result1 = run_cpca_drift_correction(base_df, SAMPLE_ROWS, QC_ROWS, n_comps=1)
    result2 = run_cpca_drift_correction(base_df, SAMPLE_ROWS, QC_ROWS, n_comps=2)
    assert not np.allclose(
        result1.loc[intensity_rows].values,
        result2.loc[intensity_rows].values,
    )


def test_run_cpca_drift_correction_integer_input(base_df):
    """Integer-typed DataFrame should be accepted and return float values."""
    df_int = base_df.astype(int)
    result = run_cpca_drift_correction(df_int, SAMPLE_ROWS, QC_ROWS)
    assert result.shape == df_int.shape
    assert result.dtypes.eq(float).all()


def test_run_cpca_drift_correction_returns_float_dtype(base_df):
    result = run_cpca_drift_correction(base_df, SAMPLE_ROWS, QC_ROWS)
    assert result.dtypes.eq(float).all()


# ---------------------------------------------------------------------------
# cpca_centroid
# ---------------------------------------------------------------------------


@pytest.fixture
def centroid_df():
    """Larger DataFrame with no missing values needed for PCA(n_components=2)."""
    rng = np.random.default_rng(0)
    data = rng.uniform(1, 100, size=(10, 8))
    index = [f"s{i}" for i in range(6)] + [f"qc{i}" for i in range(4)]
    cols = [f"f{i}" for i in range(8)]
    return pd.DataFrame(data, index=index, columns=cols)


CENTROID_SAMPLES = [f"s{i}" for i in range(6)]
CENTROID_QCS = [f"qc{i}" for i in range(4)]


def test_cpca_centroid_returns_nonneg_float(centroid_df):
    result = cpca_centroid(centroid_df, CENTROID_SAMPLES, CENTROID_QCS)
    assert isinstance(result, float)
    assert result >= 0.0


def test_cpca_centroid_list_vs_dict_sample_rows(centroid_df):
    """List and single-group dict should produce the same result."""
    result_list = cpca_centroid(centroid_df, CENTROID_SAMPLES, CENTROID_QCS)
    result_dict = cpca_centroid(
        centroid_df, {"Samples": CENTROID_SAMPLES}, CENTROID_QCS
    )
    assert pytest.approx(result_list) == result_dict


def test_cpca_centroid_no_log_transform(centroid_df):
    result = cpca_centroid(
        centroid_df, CENTROID_SAMPLES, CENTROID_QCS, log_transform=False
    )
    assert isinstance(result, float)
    assert result >= 0.0


def test_cpca_centroid_missing_rows_skipped(centroid_df):
    """Rows not in the DataFrame index are silently skipped."""
    sample_rows_with_ghost = CENTROID_SAMPLES + ["nonexistent_row"]
    result = cpca_centroid(centroid_df, sample_rows_with_ghost, CENTROID_QCS)
    expected = cpca_centroid(centroid_df, CENTROID_SAMPLES, CENTROID_QCS)
    assert pytest.approx(result) == expected


def test_cpca_centroid_multigroup_dict(centroid_df):
    """A dict with multiple sample groups must not raise."""
    half = len(CENTROID_SAMPLES) // 2
    groups = {
        "A": CENTROID_SAMPLES[:half],
        "B": CENTROID_SAMPLES[half:],
    }
    result = cpca_centroid(centroid_df, groups, CENTROID_QCS)
    assert isinstance(result, float)
    assert result >= 0.0


def test_cpca_centroid_nan_column_dropped(centroid_df):
    """A feature column that is all NaN should be silently dropped."""
    df_with_nan_col = centroid_df.copy()
    df_with_nan_col["nan_col"] = np.nan
    result = cpca_centroid(df_with_nan_col, CENTROID_SAMPLES, CENTROID_QCS)
    expected = cpca_centroid(centroid_df, CENTROID_SAMPLES, CENTROID_QCS)
    assert pytest.approx(result) == expected


def test_cpca_centroid_log_vs_no_log_differ(centroid_df):
    """log_transform=True and False should produce different centroid distances."""
    with_log = cpca_centroid(
        centroid_df, CENTROID_SAMPLES, CENTROID_QCS, log_transform=True
    )
    without_log = cpca_centroid(
        centroid_df, CENTROID_SAMPLES, CENTROID_QCS, log_transform=False
    )
    assert with_log != pytest.approx(without_log)
