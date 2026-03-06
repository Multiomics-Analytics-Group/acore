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

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from numpy import nan

from acore.imputation_analysis import (
    imputation_KNN,
    imputation_mixed_norm_KNN,
    imputation_normal_distribution,
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


def test_imputation_normal_distribution(example_data):
    actual = imputation_normal_distribution(example_data, drop_cols=None).sum().sum()
    desired = 2151.1985850814845
    np.testing.assert_almost_equal(actual, desired)


@pytest.mark.parametrize(
    ["desired", "cutoff"],
    [(2206.6954572254326, 0.3), (1684.3585255006, 0.6), (837.3284651386293, 0.9)],
)
def test_test_imputation_KNN(example_data, desired, cutoff):
    actual = imputation_KNN(example_data, cutoff=cutoff).sum().sum()
    np.testing.assert_almost_equal(actual, desired)
