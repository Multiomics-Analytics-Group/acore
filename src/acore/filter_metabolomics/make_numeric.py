"""Cleans data frame to make selected columns numeric.
Deals with entries that display two values (e.g. 3.564_3.745) by computing
an average and replacing the range with it.

"""

import pandas as pd


def parse_average(value):
    """
    Converts a string or numeric value to a float.

    If the input value contains an underscore ('_'), it splits the string into
    two numbers and returns their average. If the input is already numeric, it
    returns it as a float. If conversion fails, returns None.

    Usage: parse_average(value)

    :param (int, float, or str) value: Input value to convert. Can be numeric or a string containing
                  a single number or two numbers separated by an underscore.

    :return (float or None): Float representation of the input, the average if two numbers are
             provided, or None if conversion is not possible.
    """

    if pd.isna(value):
        return None
    # If it's already numeric, return as float

    if isinstance(value, (int, float)):
        return float(value)
    # If it's a string with underscore, split and average

    if "_" in value:
        parts = value.split("_")
        try:
            nums = [float(p) for p in parts]
            return sum(nums) / len(nums)
        except ValueError:
            return None  # cannot convert, treat as NA

    # Otherwise, try to convert directly to float
    try:
        return float(value)

    except ValueError:  # If none of the above work
        return None


def convert_to_numeric(
    df: pd.DataFrame, cols_to_convert: list, print_na_summary: bool = False
) -> pd.DataFrame:
    """
    Converts specified columns of a DataFrame to numeric values using `parse_average`.

    For each column in `cols_to_convert`, this function applies `parse_average`, which:
    - Converts numeric strings to floats.
    - Averages values in strings containing two numbers separated by an underscore.
    - Converts invalid values to `NaN`.

    After conversion, it prints a summary for each column indicating how many `NaN`
    values remain.

    Usage: convert_to_numeric(df, cols_to_convert, print_na_summary=False)

    :param pd.DataFrame df: Input DataFrame to be modified in place.
    :param (list of str) cols_to_convert: List of column names to convert to numeric
    :return pd.DataFrame:Converted DataFrame.
    """

    # Convert all to numeric using function parse_average
    for col in cols_to_convert:
        df[col] = df[col].apply(parse_average)

        if print_na_summary:
            # Check for NAs that remain
            num_na = df[col].isna().sum()
            if num_na > 0:
                print(f"Column '{col}' has {num_na} NaN values after conversion")
            else:
                print(f"Column {col} has been converted successfully, no NaN values.")

    return df
