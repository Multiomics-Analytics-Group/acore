"""
# Test Data set was created from a sample by shuffling:

fraction_missing = proteins.notna().mean()

data = data[data.columns[fraction_missing > 0.4]]
N_FEAT = 200
N_FEAT_digits = len(str(N_FEAT))
data = data.sample(N_FEAT, axis=1)
data.columns = [f"P{i:0{N_FEAT_digits}d}" for i in range(N_FEAT)]
data.reset_index(inplace=True)
data.drop('index', axis=1, inplace=True)
data.apply(numpy.random.shuffle, axis=1)
data.to_csv('test_data.csv')
"""

import logging

import numpy as np
import pandas as pd
import pytest
from numpy import nan

from acore.imputation_analysis import (
    imputation_half_minimum,
    imputation_KNN,
    imputation_mixed_norm_KNN,
    imputation_normal_distribution,
    imputation_zeros,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def small_data():
    return pd.DataFrame(
        {"A": [1.0, 2.0, nan], "B": [4.0, nan, 6.0], "C": [7.0, 8.0, 9.0]},
        index=["s1", "s2", "s3"],
    )


@pytest.fixture
def example_data():
    """
    Fixture to load example data from a csv file for testing.
    """

    sampled_dict = {
        "P112": {
            330: 27.00485326118615,
            176: 25.65507252359171,
            267: 27.92473341489751,
            205: 25.943980808304044,
            259: nan,
            384: 26.841581622745515,
            285: 25.433867342789323,
            163: 28.02912669650785,
            20: nan,
            225: 26.843978923132443,
        },
        "P029": {
            330: 26.571683331700314,
            176: 24.9274600607433,
            267: 27.642030643303336,
            205: nan,
            259: nan,
            384: 28.393579976270082,
            285: 27.643889210065804,
            163: 27.85424216217362,
            20: nan,
            225: nan,
        },
        "P182": {
            330: 30.897064286120333,
            176: 29.955045721154608,
            267: 31.58221646471536,
            205: 30.979851081103423,
            259: 31.64521283878785,
            384: nan,
            285: nan,
            163: 32.21242631603555,
            20: 29.776493976945563,
            225: 30.909921527489317,
        },
        "P199": {
            330: 27.042495516317647,
            176: 25.940310908437592,
            267: 27.79847654746734,
            205: 25.572796661829862,
            259: nan,
            384: 28.178594116473366,
            285: 27.624497855753205,
            163: 28.205364168538438,
            20: 25.431159371751118,
            225: 26.058128507338505,
        },
        "P193": {
            330: 25.381559057128122,
            176: nan,
            267: nan,
            205: nan,
            259: 23.135677977278544,
            384: nan,
            285: nan,
            163: nan,
            20: 27.64814799012052,
            225: nan,
        },
        "P085": {
            330: 27.271685713222595,
            176: 25.64995759681752,
            267: 26.554572107056842,
            205: 27.496937259755054,
            259: 28.125192488032894,
            384: 29.310058517309297,
            285: 28.21062192804231,
            163: 28.962070186129182,
            20: 26.960248217828653,
            225: 27.553766847847818,
        },
        "P010": {
            330: nan,
            176: nan,
            267: 25.229456948890068,
            205: 26.941565094869404,
            259: nan,
            384: 28.182482784706128,
            285: nan,
            163: nan,
            20: nan,
            225: nan,
        },
        "P054": {
            330: 29.718556914654666,
            176: 28.15214037206457,
            267: 30.267963547476853,
            205: 29.026147745241193,
            259: 29.022607193017727,
            384: 29.85910951891364,
            285: 29.737376918025642,
            163: 30.355104722713655,
            20: 27.647736399121325,
            225: 29.17631294981031,
        },
    }

    df = pd.DataFrame.from_dict(sampled_dict)
    return df


# ---------------------------------------------------------------------------
# imputation_mixed_norm_KNN
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ["desired", "random_state"], [(2151.079993672961, 1), (2149.893269248854, 2)]
)
def test_imputation_normal_distribution(example_data, desired, random_state):
    actual = (
        imputation_normal_distribution(
            example_data, drop_cols=None, random_state=random_state
        )
        .sum()
        .sum()
    )
    np.testing.assert_almost_equal(actual, desired)


def test_imputation_normal_distribution_drop_cols(small_data):
    df = small_data.copy()
    df["meta"] = [1.0, 2.0, 3.0]
    result = imputation_normal_distribution(df, drop_cols=["meta"], random_state=1)
    assert "meta" not in result.columns
    assert not result[["A", "B", "C"]].isna().any().any()


def test_imputation_normal_distribution_drop_cols_tuple(small_data):
    df = small_data.copy()
    df["meta"] = [1.0, 2.0, 3.0]
    result_list = imputation_normal_distribution(df, drop_cols=["meta"], random_state=7)
    result_tuple = imputation_normal_distribution(
        df, drop_cols=("meta",), random_state=7
    )
    pd.testing.assert_frame_equal(result_list, result_tuple)


def test_imputation_normal_distribution_all_nan_row_fills_zero():
    df = pd.DataFrame(
        {"A": [1.0, nan], "B": [2.0, nan], "C": [3.0, nan]},
        index=["s1", "s2"],
    )
    result = imputation_normal_distribution(df, random_state=42)
    assert (result.loc["s2"] == 0.0).all()


def test_imputation_normal_distribution_single_valid_fills_zero():
    # std == 0  ->  sigma == 0  ->  condition `sigma <= 0` triggers zeros
    df = pd.DataFrame(
        {"A": [5.0, nan], "B": [nan, nan]},
        index=["s1", "s2"],
    )
    result = imputation_normal_distribution(df, random_state=42)
    assert result.loc["s1", "B"] == pytest.approx(0.0)


def test_imputation_normal_distribution_reproducible(example_data):
    r1 = imputation_normal_distribution(example_data, random_state=99)
    r2 = imputation_normal_distribution(example_data, random_state=99)
    pd.testing.assert_frame_equal(r1, r2)


# ---------------------------------------------------------------------------
# imputation_KNN
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ["desired", "cutoff"],
    [(2206.6954572254326, 0.3), (1684.3585255006, 0.6), (837.3284651386293, 0.9)],
)
def test_imputation_KNN(example_data, desired, cutoff):
    actual = imputation_KNN(example_data, cutoff=cutoff).sum().sum()
    np.testing.assert_almost_equal(actual, desired)


def test_imputation_KNN_with_group():
    df = pd.DataFrame(
        {
            "A": [1.0, 2.0, nan, 4.0, 5.0, nan],
            "B": [10.0, nan, 12.0, 13.0, nan, 15.0],
            "C": [100.0, 110.0, 120.0, 130.0, 140.0, 150.0],
            "group": ["X", "X", "X", "Y", "Y", "Y"],
        }
    )
    result = imputation_KNN(df, group="group", cutoff=0.3, alone=False, n_neighbors=2)
    assert not result[["A", "B", "C"]].isna().any().any()
    assert "group" in result.columns


def test_imputation_KNN_non_numeric_passthrough():
    df = pd.DataFrame(
        {
            "A": [1.0, 2.0, nan, 4.0],
            "B": [10.0, nan, 12.0, 13.0],
            "label": ["a", "b", "c", "d"],
        }
    )
    result = imputation_KNN(df, cutoff=0.3, alone=False)
    assert "label" in result.columns
    assert list(result["label"]) == ["a", "b", "c", "d"]
    assert not result[["A", "B"]].isna().any().any()


def test_imputation_KNN_alone_false_retains_nan_columns():
    # Column B has only 20% valid values — below cutoff=0.5, so it is not imputed.
    # alone=False must keep it; alone=True (default) must drop it.
    df = pd.DataFrame(
        {
            "A": [1.0, 2.0, nan, 4.0, 5.0],
            "B": [nan, nan, nan, nan, 6.0],
        }
    )
    result_drop = imputation_KNN(df, cutoff=0.5, alone=True)
    result_keep = imputation_KNN(df, cutoff=0.5, alone=False)

    assert "B" not in result_drop.columns
    assert "B" in result_keep.columns
    assert result_keep["B"].isna().any()


def test_imputation_KNN_no_missing_unchanged():
    df = pd.DataFrame({"A": [1.0, 2.0, 3.0], "B": [4.0, 5.0, 6.0]})
    result = imputation_KNN(df, cutoff=0.5)
    pd.testing.assert_frame_equal(result, df)


# ---------------------------------------------------------------------------
# imputation_mixed_norm_KNN
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ["desired", "random_state"], [(2185.916268157555, 1), (2184.6842864127066, 2)]
)
def test_imputation_mixed_norm_KNN(example_data, desired, random_state):
    # P10 and P193 will be imputed using normal distribution shifting
    actual = (
        imputation_mixed_norm_KNN(
            example_data,
            cutoff=0.6,
            group=None,
            random_state=random_state,
        )
        .sum()
        .sum()
    )
    np.testing.assert_almost_equal(actual, desired)


def test_imputation_mixed_norm_KNN_fully_imputed(example_data):
    result = imputation_mixed_norm_KNN(
        example_data, cutoff=0.6, group=None, random_state=42
    )
    assert result.isna().sum().sum() == 0


def test_imputation_mixed_norm_KNN_with_group():
    df = pd.DataFrame(
        {
            "A": [1.0, 2.0, nan, 4.0, 5.0, nan],
            "B": [10.0, nan, 12.0, 13.0, nan, 15.0],
            "C": [100.0, 110.0, 120.0, 130.0, 140.0, 150.0],
            "group": ["X", "X", "X", "Y", "Y", "Y"],
        }
    )
    result = imputation_mixed_norm_KNN(df, group="group", cutoff=0.3, random_state=42)
    assert result.isna().sum().sum() == 0
    assert "group" not in result.columns


# ---------------------------------------------------------------------------
# imputation_zeros
# ---------------------------------------------------------------------------


def test_imputation_zeros_fills_all_numeric():
    df = pd.DataFrame({"A": [1.0, nan], "B": [nan, 3.0]})
    result = imputation_zeros(df)
    assert result.loc[0, "A"] == pytest.approx(1.0)
    assert result.loc[1, "A"] == pytest.approx(0.0)
    assert result.loc[0, "B"] == pytest.approx(0.0)
    assert result.loc[1, "B"] == pytest.approx(3.0)


def test_imputation_zeros_skips_non_numeric():
    df = pd.DataFrame({"A": [1.0, nan], "label": ["x", "y"]})
    result = imputation_zeros(df)
    assert result.loc[1, "A"] == pytest.approx(0.0)
    assert list(result["label"]) == ["x", "y"]


def test_imputation_zeros_on_cols():
    df = pd.DataFrame({"A": [nan, 2.0], "B": [nan, 4.0], "C": [nan, nan]})
    result = imputation_zeros(df, on_cols=["A", "B"])
    assert result["A"].isna().sum() == 0
    assert result["B"].isna().sum() == 0
    assert result["C"].isna().all()  # not in on_cols — unchanged


def test_imputation_zeros_on_rows():
    df = pd.DataFrame(
        {"A": [nan, nan, nan], "B": [nan, nan, nan]},
        index=["r1", "r2", "r3"],
    )
    result = imputation_zeros(df, on_rows=["r1", "r2"])
    assert result.loc["r1", "A"] == pytest.approx(0.0)
    assert result.loc["r2", "B"] == pytest.approx(0.0)
    assert np.isnan(result.loc["r3", "A"])
    assert np.isnan(result.loc["r3", "B"])


def test_imputation_zeros_drop_cols():
    df = pd.DataFrame({"A": [nan, 2.0], "meta": [nan, 1.0]})
    result = imputation_zeros(df, drop_cols=["meta"])
    assert "meta" not in result.columns
    assert result.loc[0, "A"] == pytest.approx(0.0)


def test_imputation_zeros_on_cols_non_numeric_raises():
    df = pd.DataFrame({"A": [1.0, 2.0], "label": ["x", "y"]})
    with pytest.raises(TypeError, match="Non-numeric columns"):
        imputation_zeros(df, on_cols=["label"])


def test_imputation_zeros_overlap_drop_cols_warning(caplog):
    df = pd.DataFrame({"A": [nan, 2.0], "B": [nan, 4.0]})
    with caplog.at_level(logging.WARNING, logger="acore.imputation_analysis"):
        result = imputation_zeros(df, on_cols=["A", "B"], drop_cols=["A"])
    assert "A" not in result.columns
    assert result["B"].isna().sum() == 0
    assert any("A" in msg for msg in caplog.messages)


# ---------------------------------------------------------------------------
# imputation_half_minimum
# ---------------------------------------------------------------------------


def test_imputation_half_minimum_fills_with_half_min():
    df = pd.DataFrame({"A": [2.0, 4.0, nan], "B": [6.0, nan, 12.0]})
    result = imputation_half_minimum(df)
    assert result.loc[2, "A"] == pytest.approx(1.0)  # min(A)=2 -> 2/2=1
    assert result.loc[1, "B"] == pytest.approx(3.0)  # min(B)=6 -> 6/2=3
    assert result.loc[0, "A"] == pytest.approx(2.0)  # unchanged
    assert result.loc[0, "B"] == pytest.approx(6.0)  # unchanged


def test_imputation_half_minimum_on_cols():
    df = pd.DataFrame({"A": [2.0, nan], "B": [nan, 8.0], "C": [nan, nan]})
    result = imputation_half_minimum(df, on_cols=["A", "B"])
    assert result.loc[1, "A"] == pytest.approx(1.0)
    assert result.loc[0, "B"] == pytest.approx(4.0)
    assert result["C"].isna().all()  # not in on_cols — unchanged


def test_imputation_half_minimum_on_rows():
    df = pd.DataFrame(
        {"A": [10.0, nan, nan], "B": [nan, 20.0, nan]},
        index=["r1", "r2", "r3"],
    )
    result = imputation_half_minimum(df, on_rows=["r1", "r2"])
    # subset r1,r2: A min=10 -> half=5; B min=20 -> half=10
    assert result.loc["r1", "B"] == pytest.approx(10.0)
    assert result.loc["r2", "A"] == pytest.approx(5.0)
    assert np.isnan(result.loc["r3", "A"])
    assert np.isnan(result.loc["r3", "B"])


def test_imputation_half_minimum_drop_cols():
    df = pd.DataFrame({"A": [2.0, nan], "meta": [10.0, 20.0]})
    result = imputation_half_minimum(df, drop_cols=["meta"])
    assert "meta" not in result.columns
    assert result.loc[1, "A"] == pytest.approx(1.0)


def test_imputation_half_minimum_all_nan_col_stays_nan(caplog):
    df = pd.DataFrame({"A": [1.0, 2.0, nan], "B": [nan, nan, nan]})
    with caplog.at_level(logging.WARNING, logger="acore.imputation_analysis"):
        result = imputation_half_minimum(df)
    assert result["B"].isna().all()
    assert any("B" in msg for msg in caplog.messages)


def test_imputation_half_minimum_on_cols_non_numeric_raises():
    df = pd.DataFrame({"A": [1.0, 2.0], "label": ["x", "y"]})
    with pytest.raises(TypeError, match="Non-numeric columns"):
        imputation_half_minimum(df, on_cols=["label"])


def test_imputation_half_minimum_unknown_row_index_warning(caplog):
    # "ghost" is not in the index — a warning is emitted; valid rows are still imputed.
    # on_rows=["r1","r2","ghost"] -> rows=["r1","r2"]; subset A=[2.0,nan] -> half-min=1.0
    df = pd.DataFrame({"A": [2.0, nan]}, index=["r1", "r2"])
    with caplog.at_level(logging.WARNING, logger="acore.imputation_analysis"):
        result = imputation_half_minimum(df, on_rows=["r1", "r2", "ghost"])
    assert any("ghost" in msg for msg in caplog.messages)
    assert result.loc["r2", "A"] == pytest.approx(1.0)
