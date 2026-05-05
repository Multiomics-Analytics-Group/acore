# acore.filter_metabolomics package

Module for filtering metabolomics feature table.

### filter_mz_rt(df: DataFrame, rt_col: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'Average Rt(min)', mz_col: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'Average Mz', mz_decimals: [tuple](https://docs.python.org/3/library/stdtypes.html#tuple) = None, mz_low: [int](https://docs.python.org/3/library/functions.html#int) = None, rt_dead_volume: [float](https://docs.python.org/3/library/functions.html#float) = None, save_removed: [bool](https://docs.python.org/3/library/functions.html#bool) = True, print_na_summary: [bool](https://docs.python.org/3/library/functions.html#bool) = True) → [tuple](https://docs.python.org/3/library/stdtypes.html#tuple)

This function filters rows from a data frame based on retention time and m/z and checks data types.
If specified by the user, it evaluates each row on whether there are NaN values and prints
a statement for each.

Data types:
: Tries to convert mz and RT columns to numeric.

M/z filtering:
: Filters out all features that have m/z decimals in a given range that are below
  a certain m/z value.

RT filtering:
: Filters out all features that have RT below a certain number (in minutes). Corresponds
  to dead volume.

If save_removed==True, the removed features are saved to a separate DataFrame that also
contains a new column, “RemovalReason”, indicating the removal reason.

Usage: filter_mz_rt(df, rt_col, mz_col, mz_decimals. mz_low, rt_dead_volume, save_removed, print_na_summary)

* **Parameters:**
  * **df** (*pd.DataFrame*) – Input DataFrame.
  * **rt_col** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Column name for retention time.
  * **mz_col** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Column name for m/z values.
  * **mz_decimals** ([*tuple*](https://docs.python.org/3/library/stdtypes.html#tuple)) – optional - Tuple specifying fractional range to filter for m/z (e.g., (0.3, 0.9)).
  * **mz_low** ([*int*](https://docs.python.org/3/library/functions.html#int)) – optional - Threshold for low m/z values.
  * **rt_dead_volume** ([*float*](https://docs.python.org/3/library/functions.html#float)) – optional - Minimum retention time threshold.
  * **save_removed** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – default=True - Whether to return a DataFrame containing removed rows.
* **Return tuple:**
  Tuple containing the cleaned DataFrame and optionally the removed features DataFrame.
  - If save_removed=True: (cleaned_df, removed_df)
  - If save_removed=False: (cleaned_df, None)

## Submodules

## acore.filter_metabolomics.filter_data module

Functions for filtering metabolomics feature table by biologically relevant features.

### filter_biological_relevance(df: DataFrame, rt_col: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'Average Rt(min)', mz_col: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'Average Mz', mz_decimals: [tuple](https://docs.python.org/3/library/stdtypes.html#tuple) = None, mz_low: [int](https://docs.python.org/3/library/functions.html#int) = None, rt_dead_volume: [float](https://docs.python.org/3/library/functions.html#float) = None, save_removed: [bool](https://docs.python.org/3/library/functions.html#bool) = True) → [tuple](https://docs.python.org/3/library/stdtypes.html#tuple)

Cleans a DataFrame by filtering rows based on retention time and m/z.

M/z filtering:
: Filters out all features that have m/z decimals in a given range that are below
  a certain m/z value.

RT filtering:
: Filters out all features that have RT below a certain number (in minutes). Corresponds
  to dead volume.

If save_removed==True, the removed features are saved to a separate DataFrame that also
contains a new column, “RemovalReason”, indicating the removal reason.

Usage: filter_biological_relevance(df, rt_col, mz_col, mz_decimals, mz_low, rt_dead_volume, save_removed)

* **Parameters:**
  * **df** (*pd.DataFrame*) – Input DataFrame.
  * **rt_col** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Column name for retention time.
  * **mz_col** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Column name for m/z values.
  * **mz_decimals** ([*tuple*](https://docs.python.org/3/library/stdtypes.html#tuple)) – optional - Tuple specifying fractional range to filter for m/z (e.g., (0.3, 0.9)).
  * **mz_low** ([*int*](https://docs.python.org/3/library/functions.html#int)) – optional - Threshold for low m/z values.
  * **rt_dead_volume** ([*float*](https://docs.python.org/3/library/functions.html#float)) – optional - Minimum retention time threshold.
  * **save_removed** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – default=True - Whether to return a DataFrame containing removed rows.
* **Return tuple:**
  Tuple containing the cleaned DataFrame and optionally the removed features DataFrame.
  - If save_removed=True: (cleaned_df, removed_df)
  - If save_removed=False: (cleaned_df, None)

## acore.filter_metabolomics.make_numeric module

Cleans data frame to make selected columns numeric.
Deals with entries that display two values (e.g. 3.564_3.745) by computing
an average and replacing the range with it.

### parse_average(value)

Converts a string or numeric value to a float.

If the input value contains an underscore (‘_’), it splits the string into
two numbers and returns their average. If the input is already numeric, it
returns it as a float. If conversion fails, returns None.

Usage: parse_average(value)

* **Parameters:**
  **value** ( *(*[*int*](https://docs.python.org/3/library/functions.html#int) *,* [*float*](https://docs.python.org/3/library/functions.html#float) *, or* [*str*](https://docs.python.org/3/library/stdtypes.html#str) *)*) – Input value to convert. Can be numeric or a string containing
  a single number or two numbers separated by an underscore.
* **Return (float or None):**
  Float representation of the input, the average if two numbers are
  provided, or None if conversion is not possible.

### convert_to_numeric(df: DataFrame, cols_to_convert: [list](https://docs.python.org/3/library/stdtypes.html#list), print_na_summary: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → DataFrame

Converts specified columns of a DataFrame to numeric values using parse_average.

For each column in cols_to_convert, this function applies parse_average, which:
- Converts numeric strings to floats.
- Averages values in strings containing two numbers separated by an underscore.
- Converts invalid values to NaN.

After conversion, it prints a summary for each column indicating how many NaN
values remain.

Usage: convert_to_numeric(df, cols_to_convert, print_na_summary=False)

* **Parameters:**
  * **df** (*pd.DataFrame*) – Input DataFrame to be modified in place.
  * **cols_to_convert** ( *(*[*list*](https://docs.python.org/3/library/stdtypes.html#list) *of* [*str*](https://docs.python.org/3/library/stdtypes.html#str) *)*) – List of column names to convert to numeric

:return pd.DataFrame:Converted DataFrame.
