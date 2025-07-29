import pytest


def test_check_numeric_dataframe():
    import pandas as pd

    from acore.types import check_numeric_dataframe

    # Create a valid numeric DataFrame
    df_valid = pd.DataFrame({"A": [1, 2, 3], "B": [4.5, 5.5, 6.5]})

    # Should not raise an error
    assert check_numeric_dataframe(df_valid) is df_valid

    # Create an invalid DataFrame with non-numeric data
    df_invalid = pd.DataFrame({"A": [1, 2, 3], "B": ["a", "b", "c"]})

    # Should raise a ValueError
    with pytest.raises(
        ValueError, match="DataFrame contains non-numeric columns: \['B'\]"
    ):
        check_numeric_dataframe(df_invalid)


def test_select_numeric_columns():
    import pandas as pd

    from acore.types import select_numeric_columns

    # Create a DataFrame with mixed types
    df_mixed = pd.DataFrame(
        {"A": [1, 2, 3], "B": ["a", "b", "c"], "C": [4.5, 5.5, 6.5]}
    )

    # Should return only numeric columns
    result = select_numeric_columns(df_mixed)
    expected = pd.DataFrame({"A": [1, 2, 3], "C": [4.5, 5.5, 6.5]})

    pd.testing.assert_frame_equal(result, expected)
