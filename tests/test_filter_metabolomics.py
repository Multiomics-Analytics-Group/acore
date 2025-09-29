import platform
import unittest

import pandas as pd

import acore.filter_metabolomics as fm

def test_filter_mz_rt():

    data = pd.DataFrame.from_dict(
        {
            0: {"FeatureID": "F1", "Average Rt(min)": 0.5, "Average Mz": 599.35},
            1: {"FeatureID": "F2", "Average Rt(min)": 1.2, "Average Mz": 601.8},
            2: {"FeatureID": "F3", "Average Rt(min)": 0.7, "Average Mz": 605.2},
            3: {"FeatureID": "F4", "Average Rt(min)": 2.0, "Average Mz": 599.5},
            4: {"FeatureID": "F5", "Average Rt(min)": 1.5, "Average Mz": 607.9},
        },
        orient="index",
    )
    
    expected_filtered = pd.DataFrame.from_dict(
        {
            1: {"FeatureID": "F2", "Average Rt(min)": 1.2, "Average Mz": 601.8},
            4: {"FeatureID": "F5", "Average Rt(min)": 1.5, "Average Mz": 607.9},
        },
        orient="index",
    )

    # Expected removed DataFrame
    expected_removed = pd.DataFrame.from_dict(
        {
            0: {"FeatureID": "F1", "Average Rt(min)": 0.5, "Average Mz": 599.35, "RemovalReason": "BelowDeadVolume"},
            1: {"FeatureID": "F3", "Average Rt(min)": 0.7, "Average Mz": 605.2, "RemovalReason": "BelowDeadVolume"},
            2: {"FeatureID": "F4", "Average Rt(min)": 2.0, "Average Mz": 599.5, "RemovalReason": "NotBiologicallyRelevant"},
        },
        orient="index",
    )

    actual_filtered, actual_removed = filter_mz_rt(
        df=data,
        rt_dead_volume=0.8,
        mz_decimals=(0.3, 0.9),
        mz_low=600,
        save_removed=True
    )

    # Convert to dict for easy assertion
    assert actual_filtered.to_dict(orient="index") == expected_filtered.to_dict(orient="index")
    assert actual_removed.to_dict(orient="index") == expected_removed.to_dict(orient="index")

