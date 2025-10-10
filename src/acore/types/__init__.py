"""Collect common types of pandas DataFrames used in the package.

Documentation of DataFrame Models API:
https://pandera.readthedocs.io/en/stable/dataframe_models.html
"""

import pandas as pd
import pandera.pandas as pa


# ? could be moved to dsp_pandas
def check_numeric_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Check if the DataFrame contains only numeric data.
    returns the DataFrame again if it is valid (allowing chaining).
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("Input must be a pandas DataFrame.")
    non_numeric_cols = df.select_dtypes(exclude="number").columns
    if not non_numeric_cols.empty:
        raise ValueError(
            f"DataFrame contains non-numeric columns: {non_numeric_cols.tolist()}"
        )
    return df


def select_numeric_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Select only numeric columns from the DataFrame."""
    if not isinstance(df, pd.DataFrame):
        raise TypeError("Input must be a pandas DataFrame.")
    ret = df.select_dtypes(include="number")

    return ret


# Schema: all columns must be numeric (int or float)
def build_schema_all_floats(df: pd.DataFrame) -> pa.DataFrameSchema:
    """Build a schema that checks if all columns are float, potentially
    containing NaN values."""
    columns = {col: pa.Column(float, nullable=True) for col in df.columns}
    schema_for_df = pa.DataFrameSchema(
        columns=columns,  # we do not know the column names
        # dtype=float,  # checks all columns have that dtype
        # checks=check_numeric_dataframe
    )
    return schema_for_df
