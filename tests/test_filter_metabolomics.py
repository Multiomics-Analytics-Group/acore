import pandas as pd
import pytest

from acore.filter_metabolomics import filter_blanks, filter_by_missingness, filter_cv

# --- filter_by_missingness ---


def test_filter_by_missingness_classic_removes_sparse_features():
    data = pd.DataFrame(
        {
            "F1": [1.0, 1.0, 1.0, 1.0, 1.0],
            "F2": [float("nan"), float("nan"), float("nan"), 1.0, 1.0],
            "F3": [1.0, 1.0, float("nan"), float("nan"), float("nan")],
        },
        index=["S1", "S2", "S3", "S4", "S5"],
    )
    result = filter_by_missingness(
        data, percent=80, method="classic", samples=["S1", "S2", "S3", "S4", "S5"]
    )
    assert list(result.columns) == ["F1"]


def test_filter_by_missingness_classic_keeps_all_at_threshold():
    data = pd.DataFrame(
        {
            "F1": [1.0, 1.0, 1.0, 1.0, float("nan")],
        },
        index=["S1", "S2", "S3", "S4", "S5"],
    )
    # 4/5 = 80% non-NA → exactly at threshold → keep
    result = filter_by_missingness(
        data, percent=80, method="classic", samples=["S1", "S2", "S3", "S4", "S5"]
    )
    assert list(result.columns) == ["F1"]


def test_filter_by_missingness_classic_preserves_non_sample_rows():
    data = pd.DataFrame(
        {
            "F1": [1.0, 1.0, 5.0],
            "F2": [float("nan"), float("nan"), 5.0],
        },
        index=["S1", "S2", "QC1"],
    )
    result = filter_by_missingness(
        data, percent=80, method="classic", samples=["S1", "S2"]
    )
    assert list(result.index) == ["S1", "S2", "QC1"]
    assert list(result.columns) == ["F1"]


def test_filter_by_missingness_modified_keeps_condition_specific_features():
    # F1 only present in treatment, F2 only in control, F3 in neither
    data = pd.DataFrame(
        {
            "F1": [1.0, 1.0, 1.0, float("nan"), float("nan")],
            "F2": [float("nan"), float("nan"), float("nan"), 1.0, 1.0],
            "F3": [
                float("nan"),
                float("nan"),
                float("nan"),
                float("nan"),
                float("nan"),
            ],
        },
        index=["S1", "S2", "S3", "S4", "S5"],
    )
    groups = {"treatment": ["S1", "S2", "S3"], "control": ["S4", "S5"]}
    result = filter_by_missingness(data, percent=80, method="modified", groups=groups)
    assert set(result.columns) == {"F1", "F2"}


def test_filter_by_missingness_raises_without_samples_for_classic():
    data = pd.DataFrame({"F1": [1.0, 2.0]}, index=["S1", "S2"])
    with pytest.raises(ValueError, match="sample names must be provided"):
        filter_by_missingness(data, method="classic")


def test_filter_by_missingness_raises_without_groups_for_modified():
    data = pd.DataFrame({"F1": [1.0, 2.0]}, index=["S1", "S2"])
    with pytest.raises(ValueError, match="groups must be provided"):
        filter_by_missingness(data, method="modified")


def test_filter_by_missingness_raises_on_unknown_method():
    data = pd.DataFrame({"F1": [1.0, 2.0]}, index=["S1", "S2"])
    with pytest.raises(ValueError, match="method must be"):
        filter_by_missingness(data, method="unknown", samples=["S1", "S2"])


# --- filter_cv ---


def test_filter_cv_removes_features_with_low_biological_variability():
    # F1: sample CV > QC CV → keep; F2: sample CV < QC CV → remove; F3: sample CV > QC CV → keep
    data = pd.DataFrame(
        {
            "F1": [1.0, 2.0, 3.0, 1.5, 1.6, 1.4],
            "F2": [2.0, 2.0, 2.0, 1.8, 2.2, 1.9],
            "F3": [3.0, 6.0, 9.0, 4.0, 4.5, 5.0],
        },
        index=["S1", "S2", "S3", "QC1", "QC2", "QC3"],
    )
    result = filter_cv(data, samples=["S1", "S2", "S3"], qcs=["QC1", "QC2", "QC3"])
    assert set(result.columns) == {"F1", "F3"}


def test_filter_cv_keeps_all_when_samples_more_variable():
    data = pd.DataFrame(
        {
            "F1": [10.0, 20.0, 30.0, 15.0, 16.0, 14.0],
        },
        index=["S1", "S2", "S3", "QC1", "QC2", "QC3"],
    )
    result = filter_cv(data, samples=["S1", "S2", "S3"], qcs=["QC1", "QC2", "QC3"])
    assert list(result.columns) == ["F1"]


def test_filter_cv_raises_with_single_qc():
    data = pd.DataFrame({"F1": [1.0, 2.0]}, index=["S1", "QC1"])
    with pytest.raises(ValueError, match="more than 1 QC sample"):
        filter_cv(data, samples=["S1"], qcs=["QC1"])


def test_filter_cv_drops_features_with_undefined_qc_cv(caplog):
    # F1 QC mean = 0 → CV undefined → feature dropped
    data = pd.DataFrame(
        {
            "F1": [1.0, 2.0, 3.0, 0.0, 0.0],
            "F2": [1.0, 2.0, 3.0, 1.0, 2.0],
        },
        index=["S1", "S2", "S3", "QC1", "QC2"],
    )
    import logging

    with caplog.at_level(logging.WARNING, logger="acore.filter_metabolomics"):
        result = filter_cv(data, samples=["S1", "S2", "S3"], qcs=["QC1", "QC2"])
    assert "F1" not in result.columns
    assert "CV is undefined" in caplog.text


# --- filter_blanks ---


def test_filter_blanks_removes_blank_dominated_features():
    data = pd.DataFrame(
        {
            "F1": [100.0, 100.0, 60.0, 60.0],
            "F2": [100.0, 100.0, 5.0, 5.0],
            "F3": [100.0, 100.0, 1.0, 1.0],
        },
        index=["S1", "S2", "Blank1", "Blank2"],
    )
    result = filter_blanks(data, blanks=["Blank1", "Blank2"], samples=["S1", "S2"])
    # F1: ratio ~0.6 > 0.5 → removed; F2, F3: ratios < 0.5 → kept
    assert set(result.columns) == {"F2", "F3"}


def test_filter_blanks_custom_threshold():
    data = pd.DataFrame(
        {
            "F1": [100.0, 100.0, 40.0, 40.0],
            "F2": [100.0, 100.0, 10.0, 10.0],
        },
        index=["S1", "S2", "Blank1", "Blank2"],
    )
    # With threshold=0.3: F1 ratio ~0.4 > 0.3 → removed; F2 ratio ~0.1 ≤ 0.3 → kept
    result = filter_blanks(
        data, blanks=["Blank1", "Blank2"], samples=["S1", "S2"], threshold=0.3
    )
    assert list(result.columns) == ["F2"]


def test_filter_blanks_keeps_all_when_blanks_are_clean():
    data = pd.DataFrame(
        {
            "F1": [100.0, 100.0, 1.0, 1.0],
        },
        index=["S1", "S2", "Blank1", "Blank2"],
    )
    result = filter_blanks(data, blanks=["Blank1", "Blank2"], samples=["S1", "S2"])
    assert list(result.columns) == ["F1"]


def test_filter_blanks_preserves_all_rows():
    data = pd.DataFrame(
        {
            "F1": [100.0, 100.0, 1.0, 1.0],
            "F2": [100.0, 100.0, 80.0, 80.0],
        },
        index=["S1", "S2", "Blank1", "Blank2"],
    )
    result = filter_blanks(data, blanks=["Blank1", "Blank2"], samples=["S1", "S2"])
    assert list(result.index) == ["S1", "S2", "Blank1", "Blank2"]
