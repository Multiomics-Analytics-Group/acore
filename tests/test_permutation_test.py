from acore import permutation_test as pt
import pytest
import numpy as np

@pytest.fixture
def data1():
    return np.arange(1, 10)

@pytest.fixture
def data2():
    return np.arange(31, 40)

def test_paired_permutation(data1):
    cond1 = data1
    cond2 = data1
    expected = {
        "metric": np.mean,
        "observed": 0.0,
        "p_value": 1.0,
    }
    result = pt.paired_permutation(cond1, cond2, metric=np.mean)
    assert result["metric"] == expected["metric"]
    assert result["observed"] == expected["observed"]
    assert result["p_value"] == expected["p_value"]

def test_paired_permutation_reject(data1, data2):
    cond1 = data1
    cond2 = data2
    for stat in ["median"]:
        result = pt.paired_permutation(cond1, cond2, metric=stat)
        print(result)
        assert result["p_value"] < 0.05


