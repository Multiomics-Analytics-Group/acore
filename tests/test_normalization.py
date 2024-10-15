import pandas as pd

from acore import normalization as normalization


def test_combat_batch_correction():
    data = pd.DataFrame.from_dict(
        {
            0: {"a": 2, "b": 4, "c": 4, "batch": "A"},
            1: {"a": 5, "b": 4, "c": 14, "batch": "B"},
            2: {"a": 4, "b": 6, "c": 8, "batch": "A"},
            3: {"a": 3, "b": 5, "c": 8, "batch": "A"},
            4: {"a": 3, "b": 3, "c": 9, "batch": "B"},
        },
        orient="index",
    )
    expected = {
        0: {"a": 2.2879985660018765, "b": 3.682428536138729, "c": 5.293448325592079},
        1: {"a": 4.482588006427705, "b": 4.5907334328411995, "c": 11.70033324453306},
        2: {"a": 4.322327732200314, "b": 5.661304683114737, "c": 9.362703235564519},
        3: {"a": 3.3051631491010953, "b": 4.671866609626733, "c": 9.362703235564519},
        4: {"a": 2.537904319175784, "b": 3.589114860839033, "c": 6.885375054723476},
    }
    actual = normalization.combat_batch_correction(data, "batch").to_dict(
        orient="index"
    )
    assert actual == expected


def test_median_normalization_along_columns():
    data = pd.DataFrame(
        {
            0: {"a": 2, "b": 4, "c": 4},
            1: {"a": 5, "b": 4, "c": 14},
            2: {"a": 4, "b": 6, "c": 8},
            3: {"a": 3, "b": 5, "c": 8},
            4: {"a": 3, "b": 3, "c": 9},
        }
    ).T

    # alberto's version
    # expected_1 = {
    #     "a": [-1.333333, -2.666667, -2.000000, -2.333333, -2.000000],
    #     "b": [0.666667, -3.666667, 0.000000, -0.333333, -2.000000],
    #     "c": [0.666667, 6.333333, 2.000000, 2.666667, 4.000000],
    # }

    expected = {
        0: {"a": 3.0, "b": 5.0, "c": 5.0},
        1: {"a": 5.0, "b": 4.0, "c": 14.0},
        2: {"a": 3.0, "b": 5.0, "c": 7.0},
        3: {"a": 3.0, "b": 5.0, "c": 8.0},
        4: {"a": 5.0, "b": 5.0, "c": 11.0},
    }
    actual = normalization.median_normalization(data, normalize="samples").to_dict(
        orient="index"
    )

    assert actual == expected


def test_median_zero_normalization_along_columns():

    data = pd.DataFrame(
        {
            0: {"a": 2, "b": 4, "c": 4},
            1: {"a": 5, "b": 4, "c": 14},
            2: {"a": 4, "b": 6, "c": 8},
            3: {"a": 3, "b": 5, "c": 8},
            4: {"a": 3, "b": 3, "c": 9},
        }
    ).T
    expected = {
        0: {"a": -2.0, "b": 0.0, "c": 0.0},
        1: {"a": 0.0, "b": -1.0, "c": 9.0},
        2: {"a": -2.0, "b": 0.0, "c": 2.0},
        3: {"a": -2.0, "b": 0.0, "c": 3.0},
        4: {"a": 0.0, "b": 0.0, "c": 6.0},
    }
    actual = normalization.median_zero_normalization(data, normalize="samples").to_dict(
        orient="index"
    )

    assert actual == expected


def test_zscore_normalization_along_columns():
    data = pd.DataFrame(
        {
            0: {"a": 2, "b": 4, "c": 4},
            1: {"a": 5, "b": 4, "c": 14},
            2: {"a": 4, "b": 6, "c": 8},
            3: {"a": 3, "b": 5, "c": 8},
            4: {"a": 3, "b": 3, "c": 9},
        }
    ).T
    expected = {
        0: {"a": -1.1547005383792517, "b": 0.5773502691896256, "c": 0.5773502691896256},
        1: {
            "a": -0.48418202613504197,
            "b": -0.6657502859356828,
            "c": 1.1499323120707245,
        },
        2: {"a": -1.0, "b": 0.0, "c": 1.0},
        3: {
            "a": -0.9271726499455306,
            "b": -0.13245323570650427,
            "c": 1.0596258856520353,
        },
        4: {
            "a": -0.5773502691896258,
            "b": -0.5773502691896258,
            "c": 1.1547005383792517,
        },
    }
    actual = normalization.zscore_normalization(data, normalize="samples").to_dict(
        orient="index"
    )

    assert actual == expected


def test_median_polish_normalizaton():
    data = pd.DataFrame(
        {
            0: {"a": 2, "b": 4, "c": 4},
            1: {"a": 5, "b": 4, "c": 14},
            2: {"a": 4, "b": 6, "c": 8},
            3: {"a": 3, "b": 5, "c": 8},
            4: {"a": 3, "b": 3, "c": 9},
        }
    ).T
    expected = {
        0: {"a": 2.0, "b": 4.0, "c": 7.0},
        1: {"a": 5.0, "b": 7.0, "c": 10.0},
        2: {"a": 4.0, "b": 6.0, "c": 9.0},
        3: {"a": 3.0, "b": 5.0, "c": 8.0},
        4: {"a": 3.0, "b": 5.0, "c": 8.0},
    }
    actual = normalization.median_polish_normalization(data).to_dict(orient="index")

    assert actual == expected


def test_quantile_normalization_along_index():
    data = pd.DataFrame(
        {
            0: {"a": 2, "b": 4, "c": 4},
            1: {"a": 5, "b": 4, "c": 14},
            2: {"a": 4, "b": 6, "c": 8},
            3: {"a": 3, "b": 5, "c": 8},
            4: {"a": 3, "b": 3, "c": 9},
        }
    ).T
    expected = {
        0: {"a": 3.2, "b": 4.6, "c": 4.6},
        1: {"a": 4.6, "b": 3.2, "c": 8.6},
        2: {"a": 3.2, "b": 4.6, "c": 8.6},
        3: {"a": 3.2, "b": 4.6, "c": 8.6},
        4: {"a": 3.2, "b": 3.2, "c": 8.6},
    }
    actual = normalization.quantile_normalization(data).to_dict(orient="index")
    assert actual == expected


def test_linear_normalization_along_columns():
    data = pd.DataFrame(
        {
            0: {"a": 2, "b": 4, "c": 4},
            1: {"a": 5, "b": 4, "c": 14},
            2: {"a": 4, "b": 6, "c": 8},
            3: {"a": 3, "b": 5, "c": 8},
            4: {"a": 3, "b": 3, "c": 9},
        }
    ).T
    expected = {
        0: {
            "a": 0.11764705882352941,
            "b": 0.18181818181818182,
            "c": 0.09302325581395349,
        },
        1: {
            "a": 0.29411764705882354,
            "b": 0.18181818181818182,
            "c": 0.32558139534883723,
        },
        2: {
            "a": 0.23529411764705882,
            "b": 0.2727272727272727,
            "c": 0.18604651162790697,
        },
        3: {
            "a": 0.17647058823529413,
            "b": 0.22727272727272727,
            "c": 0.18604651162790697,
        },
        4: {
            "a": 0.17647058823529413,
            "b": 0.13636363636363635,
            "c": 0.20930232558139536,
        },
    }

    actual = normalization.linear_normalization(
        data, method="l1", normalize="samples"
    ).to_dict(orient="index")
    assert actual == expected
